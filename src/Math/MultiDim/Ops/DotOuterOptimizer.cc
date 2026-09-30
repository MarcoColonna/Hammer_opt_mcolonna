///
/// @file  DotOuterOptimizer.cc
/// @brief Outer Tensor dot product contractions
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-

// This optimizer uses a graph-based approach to identify and group equivalent
// tensor contractions in outer product structures. The key steps are:
// 1. Build a bipartite graph where vertices represent tensor blocks and edges represent contractions
// 2. Group equivalent edges/subgraphs based on label signatures and contraction patterns
// 3. Identify identical subgraphs that share the same tensor data
// 4. Iteratively collapse identical subgraphs to minimize redundant computations
// The optimizer handles three types of contraction patterns:
//   - RANK1_EDGE: Contractions between two rank-1 tensors (vectors)
//   - STAR: Multiple rank-1 tensors contracting with a single higher-rank tensor
//   - GENERAL_EDGE: All other contraction patterns

#include <unordered_set>
#include <algorithm>

#include <boost/functional/hash.hpp>
#include <boost/iterator/filter_iterator.hpp>
#include <boost/container/small_vector.hpp>

using boost::container::small_vector;

#include "Hammer/Math/MultiDim/Ops/DotOuterOptimizer.hh"
#include "Hammer/Math/MultiDim/OuterContainer.hh"
#include "Hammer/Math/MultiDim/BlockIndexing.hh"

using namespace std;

namespace Hammer::MultiDimensional {

    using OTensor = OuterContainer;

    namespace Ops {

        // Vertex indexing uses a padding scheme to separate vertex types in a single ID space:
        // - LEFT vertices: [0, VertexPadding)
        // - RIGHT vertices: [VertexPadding, 2*VertexPadding)
        // - COLLAPSED vertices: [2*VertexPadding, 3*VertexPadding)
        // This allows quick type identification and ensures unique IDs across all vertex types.
        static constexpr DotOuterOptimizer::VertexId VertexPadding = 64;

        DotOuterOptimizer::VertexId DotOuterOptimizer::Vertex::getIndex(Type type, VertexId position) {
            switch (type) {
            case Type::LEFT:
                return position;
            case Type::RIGHT:
                return position + VertexPadding;
            case Type::COLLAPSED:
                return position + 2 * VertexPadding;
            case Type::UNUSED:
                throw Error("Cannot produce index of unused vertex");
            default:
                throw Error("Unknown SubGrGroup type");
            }
        }

        DotOuterOptimizer::VertexId DotOuterOptimizer::Vertex::getIndex() const {
            return Vertex::getIndex(_type, _outerPosition);
        }

        DotOuterOptimizer::Edge::Edge(VertexId from, VertexId to, IndexPairList localPositions, bool isHc)
            : _fromVertex{from}, _toVertex{to}, _localPositions{std::move(localPositions)}, _isHc{isHc},
              _contractionsUID{} {
            updateUID();
        }


        // The contraction UID is a hash of the local index positions being contracted.
        // Edges with the same UID perform equivalent contractions (same index pattern)
        // even if they connect different vertices.
        void DotOuterOptimizer::Edge::updateUID() {
            _contractionsUID = boost::hash<IndexPairList>()(_localPositions);
        }

        DotOuterOptimizer::DotOuterOptimizer(const OuterContainer& first,
                                             const OuterContainer& second, // NOLINT(bugprone-exception-escape)
                                             const IndexPairList& indices) noexcept
            : _first{first}, _second{second}, _indices{indices}, _vertexDisjointSet{} {
            collectInfo();
        }

        DotOuterOptimizer::VertexId DotOuterOptimizer::addVertex(Vertex::Type type, IndexType position,
                                                                 LabelSignature signature, uint8_t rank) {
            auto pos = static_cast<VertexId>(position);
            auto idx = Vertex::getIndex(type, pos);
            _vertices[idx] = Vertex{type, pos, rank, signature, false, SharedTensorData{}};
            if (type == Vertex::Type::LEFT || type == Vertex::Type::RIGHT) {
                _vertexDisjointSet.setParent(idx, idx);
            }
            return idx;
        }

        void DotOuterOptimizer::addEdge(VertexId leftId, VertexId rightId, IndexPairList localPositions, bool isHc) {
            _edges.emplace_back(leftId, rightId, localPositions, isHc);
        }

        void DotOuterOptimizer::assignDataToVertex(VertexId vertexId, SharedTensorData data, bool isHc) {
            _vertices[vertexId]._data = std::move(data);
            _vertices[vertexId]._isHc = isHc;
        }

        void DotOuterOptimizer::assignDataToVertices(const OuterContainer::EntryType& left,
                                                     const OuterContainer::EntryType& right) {
            removeCollapsedVertices();
            for (VertexId i = 0; i < static_cast<VertexId>(left.size()); ++i) {
                auto vertexId = Vertex::getIndex(Vertex::Type::LEFT, i);
                _vertices[vertexId]._data = left[i].first;
                _vertices[vertexId]._isHc = left[i].second;
                _vertexDisjointSet.setParent(vertexId, vertexId);
            }
            for (VertexId i = 0; i < static_cast<VertexId>(right.size()); ++i) {
                auto vertexId = Vertex::getIndex(Vertex::Type::RIGHT, i);
                _vertices[vertexId]._data = right[i].first;
                _vertices[vertexId]._isHc = right[i].second;
                _vertexDisjointSet.setParent(vertexId, vertexId);
            }
        }

        void DotOuterOptimizer::removeCollapsedVertices() {
            for (auto it = _vertices.begin(); it != _vertices.end();) {
                if (it->second._type == Vertex::Type::COLLAPSED) {
                    it = _vertices.erase(it);
                } else {
                    ++it;
                }
            }
            _nextVertexId = 0;
        }

        // This function groups edges into equivalence classes based on:
        // 1. The label signatures of the connected vertices (accounting for hermitian conjugation)
        // 2. The contraction pattern (which indices are being contracted)
        // Edges in the same equivalence class perform the same mathematical operation
        // and can potentially be batched together. The grouping process:
        // - Separates rank-1 edges (at least one endpoint is rank-1) from general edges
        // - Edges between two rank-1 tensors are treated separately
        // - For rank-1 edges, identifies "star" patterns where multiple rank-1 tensors
        //   contract with the same higher-rank tensor
        // - Groups equivalent stars together for batch processing
        void DotOuterOptimizer::fillEquivalentSubGrGroups() {
            using Key = tuple<LabelSignature, LabelSignature, UIDType>;
            map<Key, small_vector<EdgeId, 16>> edgeGroupsMap;

            for (EdgeId i = 0; i < static_cast<EdgeId>(_edges.size()); ++i) {
                const auto& e = _edges[i];
                const Vertex& left = _vertices.at(e._fromVertex);
                const Vertex& right = _vertices.at(e._toVertex);

                if (left._type == Vertex::Type::COLLAPSED || right._type == Vertex::Type::COLLAPSED) {
                    continue;
                }

                // Construct a canonical key for edge equivalence:
                // - Left vertex signature (always non-HC for canonicalization)
                // - Right vertex signature (flip HC if left was HC to maintain consistency)
                // - Contraction UID (pattern of indices being contracted)
                // This ensures edges with equivalent operations map to the same key.
                bool isHc = left._labelsUID.second;
                Key canonical = {{left._labelsUID.first, false},
                                 {right._labelsUID.first, isHc ? !right._labelsUID.second : right._labelsUID.second},
                                 e._contractionsUID};

                edgeGroupsMap[canonical].push_back(i);
                _edges[i]._isHc = isHc;
            }

            _equivalentSubGrGroups.clear();
            vector<SubGrGroup> genedges;
            vector<small_vector<EdgeId, 16>> genLR1edges;
            vector<small_vector<EdgeId, 16>> genRR1edges;
            for (auto& entry : edgeGroupsMap) {
                auto& group = entry.second;
                const Vertex& l = _vertices.at(_edges[group.front()]._fromVertex);
                const Vertex& r = _vertices.at(_edges[group.front()]._toVertex);
                if (l._rank == 1 && r._rank == 1) {
                    _equivalentSubGrGroups.emplace_back(SubGrGroup{SubGrGroup::Type::RANK1_EDGE, group});
                } else if (l._rank == 1) {
                    sort(group.begin(), group.end(),
                         [&](EdgeId a, EdgeId b) -> bool { return _edges[a]._toVertex < _edges[b]._toVertex; });
                    genLR1edges.push_back(std::move(group));
                } else if (r._rank == 1) {
                    sort(group.begin(), group.end(),
                         [&](EdgeId a, EdgeId b) -> bool { return _edges[a]._fromVertex < _edges[b]._fromVertex; });
                    genRR1edges.push_back(std::move(group));
                } else {
                    genedges.emplace_back(SubGrGroup{SubGrGroup::Type::GENERAL_EDGE, group});
                }
            }

            // prepare finding stars
            // Find "right star" patterns: multiple left rank-1 tensors contracting with
            // the same right higher-rank tensor. Group by the shared right vertex.
            small_vector<pair<VertexId, small_vector<std::pair<EdgeId, EquivId>, 8>>, 16> groupsLR1;
            for (size_t i = 0; i < genLR1edges.size(); ++i) {
                for (const auto& value : genLR1edges[i]) {
                    auto it = std::find_if(groupsLR1.begin(), groupsLR1.end(),
                                           [&](auto& g) { return g.first == _edges[value]._toVertex; });
                    if (it != groupsLR1.end()) {
                        it->second.emplace_back(value, i);
                    } else {
                        groupsLR1.emplace_back(_edges[value]._toVertex,
                                               small_vector<std::pair<EdgeId, EquivId>, 8>{{{value, i}}});
                    }
                }
            }
            sort(groupsLR1.begin(), groupsLR1.end(),
                 [&](const auto& a, const auto& b) -> bool { return a.second.size() > b.second.size(); });
            for (auto& elem : groupsLR1) {
                sort(elem.second.begin(), elem.second.end(),
                     [&](const pair<EdgeId, EquivId>& a, const pair<EdgeId, EquivId>& b) -> bool {
                         return a.second < b.second;
                     });
            }
            for (const auto& elem : groupsLR1) {
                _stars.emplace_back(Star{elem.second, false});
            }

            // Find "left star" patterns: multiple right rank-1 tensors contracting with
            // the same left higher-rank tensor. Group by the shared left vertex.
            small_vector<pair<VertexId, small_vector<std::pair<EdgeId, EquivId>, 8>>, 16> groupsRR1;
            for (size_t i = 0; i < genRR1edges.size(); ++i) {
                for (const auto& value : genRR1edges[i]) {
                    auto it = std::find_if(groupsRR1.begin(), groupsRR1.end(),
                                           [&](auto& g) { return g.first == _edges[value]._fromVertex; });
                    if (it != groupsRR1.end()) {
                        it->second.emplace_back(value, i);
                    } else {
                        groupsRR1.emplace_back(_edges[value]._fromVertex,
                                               small_vector<std::pair<EdgeId, EquivId>, 8>{{{value, i}}});
                    }
                }
            }
            sort(groupsRR1.begin(), groupsRR1.end(),
                 [&](const auto& a, const auto& b) { return a.second.size() > b.second.size(); });
            for (auto& elem : groupsRR1) {
                sort(elem.second.begin(), elem.second.end(),
                     [&](const pair<EdgeId, EquivId>& a, const pair<EdgeId, EquivId>& b) -> bool {
                         return a.second < b.second;
                     });
            }
            for (const auto& elem : groupsRR1) {
                _stars.emplace_back(Star{elem.second, true});
            }

            // now find equivalent stars
            // Group equivalent stars: stars are equivalent if they have the same
            // pattern of equivalence class IDs, even if they involve different tensors.
            // This allows batching stars with identical contraction structures.
            // right stars
            small_vector<small_vector<EquivId, 16>, 16> equivLR1;
            for (VertexId i = 0; i < static_cast<VertexId>(groupsLR1.size()); ++i) {
                bool placed = false;
                for (auto& group : equivLR1) {
                    const auto& otherItem = groupsLR1[group.front()];
                    if (equal(otherItem.second.begin(), otherItem.second.end(), groupsLR1[i].second.begin(),
                              groupsLR1[i].second.end(),
                              [](const pair<EdgeId, EquivId>& a, const pair<EdgeId, EquivId>& b) -> bool {
                                  return a.second == b.second;
                              })) {
                        group.push_back(i);
                        placed = true;
                        break;
                    }
                }
                if (!placed) {
                    equivLR1.emplace_back(small_vector<EquivId, 8>({i}));
                }
            }
            std::sort(equivLR1.begin(), equivLR1.end(),
                      [&](const small_vector<EquivId, 16>& a, const small_vector<EquivId, 16>& b) {
                          return (groupsLR1[a.front()].second.size() > groupsLR1[b.front()].second.size()) ||
                                 ((groupsLR1[a.front()].second.size() == groupsLR1[b.front()].second.size()) &&
                                  (a.size() > b.size()));
                      });


            // left stars
            // RR1 star i is stored at _stars[rrOffset + i] because LR1 stars occupy [0, rrOffset).
            const auto rrOffset = static_cast<EquivId>(groupsLR1.size());
            small_vector<small_vector<EquivId, 16>, 16> equivRR1;
            for (VertexId i = 0; i < static_cast<VertexId>(groupsRR1.size()); ++i) {
                bool placed = false;
                auto starId = static_cast<EquivId>(i + rrOffset);
                for (auto& group : equivRR1) {
                    const auto& otherItem = groupsRR1[group.front() - rrOffset];
                    if (equal(otherItem.second.begin(), otherItem.second.end(), groupsRR1[i].second.begin(),
                              groupsRR1[i].second.end(),
                              [](const pair<EdgeId, EquivId>& a, const pair<EdgeId, EquivId>& b) -> bool {
                                  return a.second == b.second;
                              })) {
                        group.push_back(starId);
                        placed = true;
                        break;
                    }
                }
                if (!placed) {
                    equivRR1.emplace_back(small_vector<EquivId, 8>({starId}));
                }
            }
            std::sort(equivRR1.begin(), equivRR1.end(),
                      [&](const small_vector<EquivId, 16>& a, const small_vector<EquivId, 16>& b) {
                          return (groupsRR1[a.front() - rrOffset].second.size() >
                                  groupsRR1[b.front() - rrOffset].second.size()) ||
                                 ((groupsRR1[a.front() - rrOffset].second.size() ==
                                   groupsRR1[b.front() - rrOffset].second.size()) &&
                                  (a.size() > b.size()));
                      });

            // add stars to equivalent subgroups
            for (const auto& elem : equivLR1) {
                _equivalentSubGrGroups.emplace_back(SubGrGroup{SubGrGroup::Type::STAR, elem});
            }
            for (const auto& elem : equivRR1) {
                _equivalentSubGrGroups.emplace_back(SubGrGroup{SubGrGroup::Type::STAR, elem});
            }

            // now add general edges
            for (const auto& group : genedges) {
                _equivalentSubGrGroups.emplace_back(group);
            }
        }


        void DotOuterOptimizer::createIdenticalSubGrGroups() {
            _identicalSubGrGroups.clear();

            for (const auto& equivGroup : _equivalentSubGrGroups) {
                small_vector<SubGrGroup, 32> localGroups;
                switch (equivGroup.type) {
                case SubGrGroup::Type::RANK1_EDGE:
                case SubGrGroup::Type::GENERAL_EDGE:
                    for (auto id : equivGroup.ids) {
                        bool placed = false;
                        const Edge& e = _edges[id];
                        for (auto& group : localGroups) {
                            const Edge& rep = _edges[group.ids.front()];
                            if (areEdgesIdentical(e, rep)) {
                                group.ids.push_back(id);
                                placed = true;
                                break;
                            }
                        }
                        if (!placed) {
                            localGroups.emplace_back(SubGrGroup{equivGroup.type, small_vector<uint8_t, 16>({id})});
                        }
                    }
                    break;
                case SubGrGroup::Type::STAR:
                    for (auto id : equivGroup.ids) {
                        bool placed = false;
                        const Star& s = _stars[id];
                        for (auto& group : localGroups) {
                            const Star& rep = _stars[group.ids.front()];
                            if (areStarsIdentical(s, rep)) {
                                group.ids.push_back(id);
                                placed = true;
                                break;
                            }
                        }
                        if (!placed) {
                            localGroups.emplace_back(SubGrGroup{equivGroup.type, small_vector<uint8_t, 16>({id})});
                        }
                    }
                    break;
                default:
                    throw Error("Unknown SubGrGroup type");
                }

                std::sort(localGroups.begin(), localGroups.end(),
                          [](const SubGrGroup& a, const SubGrGroup& b) { return a.ids.size() > b.ids.size(); });

                for (auto& group : localGroups) {
                    _identicalSubGrGroups.push_back(std::move(group));
                }
            }
        }

        // Two edges are "identical" (not just equivalent) if they:
        // 1. Connect tensors with the same data pointers
        // 2. Have consistent hermitian conjugation states
        // Only identical edges can be collapsed together in a single operation.
        bool DotOuterOptimizer::areEdgesIdentical(const Edge& a, const Edge& b) const {
            const Vertex& l1 = _vertices.at(a._fromVertex);
            const Vertex& r1 = _vertices.at(a._toVertex);
            const Vertex& l2 = _vertices.at(b._fromVertex);
            const Vertex& r2 = _vertices.at(b._toVertex);

            return (l1._data == l2._data && r1._data == r2._data &&
                    ((l1._isHc == l2._isHc && r1._isHc == r2._isHc) || (l1._isHc != l2._isHc && r1._isHc != r2._isHc)));
        }

        // Stars are identical if they have the same structure (size, center position)
        // and all corresponding edges are pairwise identical.
        bool DotOuterOptimizer::areStarsIdentical(const Star& a, const Star& b) const {
            if (a._edges.size() != b._edges.size() || a._isCenterLeft != b._isCenterLeft) {
                return false;
            }

            for (size_t i = 0; i < a._edges.size(); ++i) {
                if (a._edges[i].second != b._edges[i].second) {
                    return false;
                }
                if (!areEdgesIdentical(_edges[a._edges[i].first], _edges[b._edges[i].first])) {
                    return false;
                }
            }
            return true;
        }


        std::tuple<OuterContainer::EntryType, OuterContainer::EntryType, IndexPairList>
        DotOuterOptimizer::getNextSubGr() const {
            if (_identicalSubGrGroups.empty()) {
                return {{}, {}, {}};
            }

            const SubGrGroup& group = _identicalSubGrGroups.front();
            switch (group.type) {
            case SubGrGroup::Type::RANK1_EDGE: {
                OuterContainer::EntryType left;
                OuterContainer::EntryType right;
                left.emplace_back(_vertices.at(_edges[group.ids.front()]._fromVertex)._data,
                                  _vertices.at(_edges[group.ids.front()]._fromVertex)._isHc);
                right.emplace_back(_vertices.at(_edges[group.ids.front()]._toVertex)._data,
                                   _vertices.at(_edges[group.ids.front()]._toVertex)._isHc);
                return {left, right, {}};
            }
            case SubGrGroup::Type::STAR: {
                OuterContainer::EntryType left;
                OuterContainer::EntryType right;
                IndexPairList indices;
                auto cand = _stars[group.ids.front()];
                if (cand._isCenterLeft) {
                    left.emplace_back(_vertices.at(_edges[cand._edges.front().first]._fromVertex)._data,
                                      _vertices.at(_edges[cand._edges.front().first]._fromVertex)._isHc);
                    IndexType i = 0;
                    for (const auto& edge : cand._edges) {
                        right.emplace_back(_vertices.at(_edges[edge.first]._toVertex)._data,
                                           _vertices.at(_edges[edge.first]._toVertex)._isHc);
                        indices.emplace_back(_edges[edge.first]._localPositions.front().first, i);
                        ++i;
                    }
                } else {
                    right.emplace_back(_vertices.at(_edges[cand._edges.front().first]._toVertex)._data,
                                       _vertices.at(_edges[cand._edges.front().first]._toVertex)._isHc);
                    IndexType i = 0;
                    for (const auto& edge : cand._edges) {
                        left.emplace_back(_vertices.at(_edges[edge.first]._fromVertex)._data,
                                          _vertices.at(_edges[edge.first]._fromVertex)._isHc);
                        indices.emplace_back(i, _edges[edge.first]._localPositions.front().second);
                        ++i;
                    }
                }
                return {std::move(left), std::move(right), std::move(indices)};
            }
            case SubGrGroup::Type::GENERAL_EDGE: {
                OuterContainer::EntryType left;
                OuterContainer::EntryType right;
                left.emplace_back(_vertices.at(_edges[group.ids.front()]._fromVertex)._data,
                                  _vertices.at(_edges[group.ids.front()]._fromVertex)._isHc);
                right.emplace_back(_vertices.at(_edges[group.ids.front()]._toVertex)._data,
                                   _vertices.at(_edges[group.ids.front()]._toVertex)._isHc);
                return {std::move(left), std::move(right), _edges[group.ids.front()]._localPositions};
            }
            default:
                throw Error("Unknown SubGrGroup type");
            }
        }

        const DotOuterOptimizer::Vertex* DotOuterOptimizer::getVertex(VertexId vertexId) const {
            auto it = _vertices.find(vertexId);
            if (it != _vertices.end()) {
                return &(it->second);
            }
            return nullptr;
        }


        void DotOuterOptimizer::collapseSubGrGroup(const SharedTensorData& data, bool isHc,
                                                   const IndexPairList& leftIndexMap,
                                                   const IndexPairList& rightIndexMap) {
            if (_identicalSubGrGroups.empty()) {
                return;
            }

            VertexId startIdx = _nextVertexId;
            switch (_identicalSubGrGroups.front().type) {
            case SubGrGroup::Type::RANK1_EDGE:
            case SubGrGroup::Type::GENERAL_EDGE: {
                auto& group = _identicalSubGrGroups.front();
                bool baseHc = _vertices.at(_edges[group.ids.front()]._fromVertex)._isHc;
                for (auto eid : group.ids) {
                    Edge& e = _edges[eid];
                    VertexId u = _vertexDisjointSet.find(e._fromVertex);
                    VertexId v = _vertexDisjointSet.find(e._toVertex);

                    // Track hermitian conjugation through collapse: XOR of whether this vertex
                    // is a conjugate variant within the group and whether the provided data
                    // is itself the conjugate of the base-case result (flipped hash match).
                    const Vertex& vl = _vertices.at(u);
                    bool new_bool = (vl._isHc != baseHc) != isHc;

                    auto newIdx = Vertex::getIndex(Vertex::Type::COLLAPSED, _nextVertexId);
                    ++_nextVertexId;

                    LabelSignature signature = data->getLabelSignature();
                    if (isHc) {
                        signature.second = !signature.second; // flip the sign if isHc is true
                    }
                    Vertex collapsed_vertex{
                        Vertex::Type::COLLAPSED, newIdx, static_cast<uint8_t>(data->rank()), signature, new_bool, data,
                        vl._connectedID};
                    _vertices[newIdx] = collapsed_vertex;

                    _vertexDisjointSet.setParent(newIdx, newIdx);
                    _vertexDisjointSet.setParent(u, newIdx);
                    _vertexDisjointSet.setParent(v, newIdx);
                }
                break;
            }
            case SubGrGroup::Type::STAR: {
                auto& group = _identicalSubGrGroups.front();
                // check center for bool
                auto& s1 = _stars[group.ids.front()];
                auto& e1 = _edges[s1._edges.front().first];
                auto& v1 = _vertices.at((s1._isCenterLeft) ? e1._fromVertex : e1._toVertex);
                bool baseHc = v1._isHc;
                for (auto sid : group.ids) {
                    const Star& star = _stars[sid];
                    auto& center = _vertices.at((star._isCenterLeft) ? _edges[star._edges.front().first]._fromVertex
                                                                     : _edges[star._edges.front().first]._toVertex);
                    bool new_bool = (center._isHc != baseHc) != isHc;
                    auto newIdx = Vertex::getIndex(Vertex::Type::COLLAPSED, _nextVertexId);
                    ++_nextVertexId;
                    LabelSignature signature = data->getLabelSignature();
                    if (isHc) {
                        signature.second = !signature.second; // flip the sign if isHc is true
                    }
                    Vertex collapsed_vertex{
                        Vertex::Type::COLLAPSED, newIdx, static_cast<uint8_t>(data->rank()), signature, new_bool, data,
                        center._connectedID};
                    _vertices[newIdx] = collapsed_vertex;
                    _vertexDisjointSet.setParent(newIdx, newIdx);

                    for (const auto& edge : star._edges) {
                        VertexId u = _vertexDisjointSet.find(_edges[edge.first]._fromVertex);
                        VertexId v = _vertexDisjointSet.find(_edges[edge.first]._toVertex);

                        _vertexDisjointSet.setParent(u, newIdx);
                        _vertexDisjointSet.setParent(v, newIdx);
                    }
                }
                break;
            }
            default:
                throw Error("Unknown SubGrGroup type");
            }
            _identicalSubGrGroups.pop_front();

            // After collapsing a subgraph, remaining edges must update their index positions
            // to account for the reduced dimensionality of the collapsed vertices.
            // leftIndexMap and rightIndexMap provide the old->new index mapping.
            const VertexId startId = Vertex::getIndex(Vertex::Type::COLLAPSED, startIdx);
            const VertexId endId = Vertex::getIndex(Vertex::Type::COLLAPSED, _nextVertexId);
            auto findLeftIndex = [&](IndexType idx) {
                return find_if(leftIndexMap.begin(), leftIndexMap.end(),
                               [idx](const std::pair<IndexType, IndexType>& p) { return p.first == idx; })
                    ->second;
            };
            auto findRightIndex = [&](IndexType idx) {
                return find_if(rightIndexMap.begin(), rightIndexMap.end(),
                               [idx](const std::pair<IndexType, IndexType>& p) { return p.first == idx; })
                    ->second;
            };
            for (auto& group : _identicalSubGrGroups) {
                switch (group.type) {
                case SubGrGroup::Type::RANK1_EDGE:
                    break;
                case SubGrGroup::Type::GENERAL_EDGE: {
                    for (EdgeId eid : group.ids) {
                        Edge& e = _edges[eid];
                        e._fromVertex = _vertexDisjointSet.find(e._fromVertex);
                        e._toVertex = _vertexDisjointSet.find(e._toVertex);
                        bool shouldDoLeft = startId <= e._fromVertex && e._fromVertex < endId && !leftIndexMap.empty();
                        bool shouldDoRight = startId <= e._toVertex && e._toVertex < endId && !rightIndexMap.empty();
                        for (auto& elem : e._localPositions) {
                            if (shouldDoLeft) {
                                elem.first = findLeftIndex(elem.first);
                            }
                            if (shouldDoRight) {
                                elem.second = findRightIndex(elem.second);
                            }
                        }
                        e.updateUID();

#ifndef NDEBUG
                        auto tf = _vertices.at(e._fromVertex)._type;
                        auto tt = _vertices.at(e._toVertex)._type;
                        if ((tf != Vertex::Type::LEFT || tt != Vertex::Type::RIGHT) && tf != Vertex::Type::COLLAPSED &&
                            tt != Vertex::Type::COLLAPSED) {
                            throw Error("Warning: invalid edge types after collapse");
                        }
#endif // NDEBUG
                    }
                    break;
                }
                case SubGrGroup::Type::STAR: {
                    for (StarId sid : group.ids) {
                        Star& star = _stars[sid];
                        if ((leftIndexMap.empty() && star._isCenterLeft) ||
                            (rightIndexMap.empty() && !star._isCenterLeft)) {
                            continue;
                        }
                        for (auto& edge : star._edges) {
                            Edge& e = _edges[edge.first];
                            // only need to do the center vertex
                            if (star._isCenterLeft) {
                                e._fromVertex = _vertexDisjointSet.find(e._fromVertex);
                            } else {
                                e._toVertex = _vertexDisjointSet.find(e._toVertex);
                            }
                            bool shouldDoLeft = startId <= e._fromVertex && e._fromVertex < endId &&
                                                !leftIndexMap.empty() && star._isCenterLeft;
                            bool shouldDoRight = startId <= e._toVertex && e._toVertex < endId &&
                                                 !rightIndexMap.empty() && !star._isCenterLeft;
                            for (auto& elem : e._localPositions) {
                                if (shouldDoLeft) {
                                    elem.first = findLeftIndex(elem.first);
                                } else if (shouldDoRight) {
                                    elem.second = findRightIndex(elem.second);
                                }
                            }
                            e.updateUID();
#ifndef NDEBUG
                            auto tf = _vertices.at(e._fromVertex)._type;
                            auto tt = _vertices.at(e._toVertex)._type;
                            if ((tf != Vertex::Type::LEFT || tt != Vertex::Type::RIGHT) &&
                                tf != Vertex::Type::COLLAPSED && tt != Vertex::Type::COLLAPSED) {
                                throw Error("Warning: invalid edge types after collapse");
                            }
#endif // NDEBUG
                        }
                    }
                    break;
                }
                default:
                    throw Error("Unknown SubGrGroup type");
                }
            }
        }


        // Collect results after optimization:
        // 1. Find unique collapsed vertices (multiple may point to same result via disjoint set)
        // 2. Separate rank-0 results (scalars) - multiply into the weight
        // 3. Sort remaining tensors by connected component ID for consistent ordering
        // 4. Include any unused vertices that weren't involved in contractions
        pair<OuterContainer::EntryType, OuterContainer::ElementType> DotOuterOptimizer::retrieveNewData() const {
            /// TODO: review
            OuterContainer::EntryType results;
            results.reserve(_maxTensorCount);
            complex<double> weight = 1.0;

            std::unordered_set<VertexId> collapsedIdsSet;
            for (const auto& [id, vertex] : _vertices) {
                if (vertex._type == Vertex::Type::COLLAPSED) {
                    collapsedIdsSet.insert(_vertexDisjointSet.find(id));
                }
                if (vertex._type == Vertex::Type::UNUSED) {
                    results.emplace_back(vertex._data, vertex._isHc);
                }
            }
            std::vector<VertexId> collapsedIds;
            collapsedIds.reserve(collapsedIdsSet.size());
            for (const auto& id : collapsedIdsSet) {
                const auto& elem = _vertices.at(id);
                if (elem._data->rank() == 0) {
                    weight *= (elem._isHc ? conj(elem._data->element()) : elem._data->element());
                } else {
                    collapsedIds.push_back(id);
                }
            }
            std::sort(collapsedIds.begin(), collapsedIds.end(), [&](VertexId left, VertexId right) -> bool {
                return _vertices.at(left)._connectedID < _vertices.at(right)._connectedID;
            });
            for (auto id : collapsedIds) {
                const auto& vertex = _vertices.at(id);
                results.emplace_back(vertex._data, vertex._isHc);
            }
            return {results, weight};
        }

        // Use a separate disjoint set to identify connected components in the graph.
        // This determines how many independent tensor products will remain after
        // all contractions are performed, which is needed for result allocation.
        void DotOuterOptimizer::estimateTensorCount() {
            DisjointSet<VertexId> dsu;
            for (const auto& e : _edges) {
                dsu.unite(e._fromVertex, e._toVertex);
            }
            std::unordered_set<VertexId> roots;
            for (const auto& [id, vertex] : _vertices) {
                auto root = dsu.find(id);
                roots.insert(root);
                const_cast<Vertex&>(vertex)._connectedID = static_cast<ConnectedId>(root);
            }
            _maxTensorCount = roots.size();
        }

        // Build the initial bipartite graph from the input tensors:
        // 1. Create vertices for each block in the left and right outer containers
        // 2. Create edges for each index pair being contracted
        // 3. Mark vertices as UNUSED if they don't participate in any contractions
        // 4. Estimate the final tensor count for result allocation
        void DotOuterOptimizer::collectInfo() {
            vector<pair<VertexId, bool>> usedVertices;
            usedVertices.reserve(_first.getIndexing().numSubIndexing() + _second.getIndexing().numSubIndexing());
            for (VertexId i = 0; i < static_cast<VertexId>(_first.getIndexing().numSubIndexing()); ++i) {
                auto idx = addVertex(Vertex::Type::LEFT, i, _first.getIndexing().getBlockLabelSignature(i),
                                     static_cast<uint8_t>(_first.getIndexing().getSubIndexing(i).rank()));
                usedVertices.emplace_back(idx, false);
            }
            for (VertexId i = 0; i < static_cast<VertexId>(_second.getIndexing().numSubIndexing()); ++i) {
                auto idx = addVertex(Vertex::Type::RIGHT, i, _second.getIndexing().getBlockLabelSignature(i),
                                     static_cast<uint8_t>(_second.getIndexing().getSubIndexing(i).rank()));
                usedVertices.emplace_back(idx, false);
            }
            map<VertexIdPair, IndexPairList> protoEdges;
            for (const auto& elem : _indices) {
                IndexPair leftIndex = _first.getIndexing().getElementIndex(elem.first);
                IndexPair rightIndex = _second.getIndexing().getElementIndex(elem.second);
                protoEdges[{static_cast<VertexId>(leftIndex.first), static_cast<VertexId>(rightIndex.first)}]
                    .emplace_back(leftIndex.second, rightIndex.second);
            }
            for (const auto& [key, localPositions] : protoEdges) {
                VertexId leftId = Vertex::getIndex(Vertex::Type::LEFT, key.first);
                VertexId rightId = Vertex::getIndex(Vertex::Type::RIGHT, key.second);
                auto it = lower_bound(usedVertices.begin(), usedVertices.end(), leftId,
                                      [](const auto& a, auto b) -> bool { return a.first < b; });
                if (it != usedVertices.end()) {
                    it->second = true;
                }
                it = lower_bound(usedVertices.begin(), usedVertices.end(), rightId,
                                 [](const auto& a, auto b) -> bool { return a.first < b; });
                if (it != usedVertices.end()) {
                    it->second = true;
                }
                addEdge(leftId, rightId, localPositions);
            }
            for (auto& elem : usedVertices) {
                if (!elem.second) {
                    const_cast<Vertex*>(getVertex(elem.first))->_type = Vertex::Type::UNUSED;
                }
            }
            estimateTensorCount();
        }


    } // namespace Ops


} // namespace Hammer::MultiDimensional
