///
/// @file  DotOuterOptimizer.hh
/// @brief Outer Tensor dot product contractions
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_MATH_MULTIDIM_OPS_DOTOUTEROPTIMIZER
#define HAMMER_MATH_MULTIDIM_OPS_DOTOUTEROPTIMIZER

#ifndef DOTOUTEROPTIMIZER_FRIENDS
#define DOTOUTEROPTIMIZER_FRIENDS typedef bool DummyFriend
#endif

#include <cstdint>
#include <vector>
#include <deque>
#include <map>
#include <tuple>
#include <utility>
#include <boost/container/small_vector.hpp>

#include "Hammer/Tools/Utils.hh"
#include "Hammer/Math/MultiDimensional.fhh"
#include "Hammer/Math/MultiDim/OuterContainer.hh"

namespace Hammer::MultiDimensional::Ops {

    class DotOuterOptimizer final {

        DOTOUTEROPTIMIZER_FRIENDS;

    public:

        explicit DotOuterOptimizer(const OuterContainer& first, const OuterContainer& second,
                                   const IndexPairList& indices) noexcept;

        void collectInfo();


        using EdgeId = uint8_t;
        using VertexId = uint8_t;
        using StarId = uint8_t;
        using ConnectedId = int8_t;
        using UIDType = size_t;
        using EquivId = uint8_t;
        using VertexIdPair = std::pair<VertexId, VertexId>;

    private:

        const OuterContainer& _first;
        const OuterContainer& _second;
        const IndexPairList& _indices;

        struct Vertex {
            enum class Type : std::uint8_t { LEFT, RIGHT, COLLAPSED, UNUSED };

            Vertex() = default;

            Type _type;
            VertexId _outerPosition;
            uint8_t _rank;
            LabelSignature _labelsUID;
            bool _isHc;
            SharedTensorData _data;
            ConnectedId _connectedID = -1;

            [[nodiscard]] VertexId getIndex() const;
            static VertexId getIndex(Type type, VertexId position);
        };

        struct Edge {
            Edge(VertexId from, VertexId to, IndexPairList localPositions, bool isHc = false);
            Edge() = default;
            Edge(const Edge&) = default;
            Edge& operator=(const Edge&) = default;
            Edge(Edge&&) = default;
            Edge& operator=(Edge&&) = default;
            ~Edge() = default;

            void updateUID();

            VertexId _fromVertex;
            VertexId _toVertex;
            IndexPairList _localPositions;
            bool _isHc;
            UIDType _contractionsUID;
        };


        struct Star {
            boost::container::small_vector<std::pair<EdgeId, EquivId>, 8> _edges;
            bool _isCenterLeft;
        };


        std::map<VertexId, Vertex> _vertices;
        std::vector<Edge> _edges;

        std::map<VertexId, ConnectedId> _connectedIds;

        std::vector<Star> _stars;

        struct SubGrGroup {
            enum class Type : std::uint8_t { RANK1_EDGE, STAR, GENERAL_EDGE };
            Type type;
            boost::container::small_vector<uint8_t, 16> ids;
        };


        std::vector<SubGrGroup> _equivalentSubGrGroups;
        std::deque<SubGrGroup> _identicalSubGrGroups;

        mutable DisjointSet<VertexId> _vertexDisjointSet;

        VertexId _nextVertexId = 0;
        size_t _maxTensorCount = 0;


        VertexId addVertex(Vertex::Type type, IndexType position, LabelSignature signature, uint8_t rank);
        void addEdge(VertexId leftId, VertexId rightId, IndexPairList localPositions, bool isHc = false);

        void assignDataToVertex(VertexId vertexId, SharedTensorData data, bool isHc = false);
        void removeCollapsedVertices();

        const Vertex* getVertex(VertexId vertexId) const;

        void estimateTensorCount();

        bool areEdgesIdentical(const Edge& a, const Edge& b) const;

        bool areStarsIdentical(const Star& a, const Star& b) const;

    public:

        void fillEquivalentSubGrGroups();
        void assignDataToVertices(const OuterContainer::EntryType& left, const OuterContainer::EntryType& right);
        void createIdenticalSubGrGroups();

        std::tuple<OuterContainer::EntryType, OuterContainer::EntryType, IndexPairList> getNextSubGr() const;

        void collapseSubGrGroup(const SharedTensorData& data, bool isHc, const IndexPairList& leftIndexMap,
                                const IndexPairList& rightIndexMap);

        std::pair<OuterContainer::EntryType, OuterContainer::ElementType> retrieveNewData() const;


    private:
    };


    // public:

    //     static constexpr IndexType NewChunkOffset = 100;


    //     OuterContainer::EntryType performDot(const OuterContainer::EntryType& left,
    //                                          const OuterContainer::EntryType& right);

    //     std::vector<std::tuple<IndexPair, uint, uint>> getSinglesBoth(IndexType iA, IndexType iB) const;
    //     std::vector<std::tuple<IndexType, IndexPairList, uint, uint, IndexType>>
    //     getSinglesLeft(IndexType iA, IndexType iB) const;
    //     std::vector<std::tuple<IndexType, IndexPairList, uint, uint, IndexType>>
    //     getSinglesRight(IndexType iA, IndexType iB) const;

    //     void findRepetitions();

    // protected:

    //     IndexType calcShift(IndexType idx, const IndexList& contractedIndices, bool isLeft) const;


    //     ContractionGraph _graph;


    //     DotGroupList _chunks;

    //     // Edges _edges;
    //     // EdgeCounts _multiplicities;

    //     std::map<IndexPair, IndexPairList> _contractions;
    //     std::map<IndexType, IndexList> _contractionsWithSinglesLeft;
    //     std::map<IndexType, IndexList> _contractionsWithSinglesRight;
    //     IndexPairList _contractionsWithSinglesBoth;
    //     IndexPairList _contractionsWithSinglesNone;
    //     std::map<IndexType, std::map<IndexType, std::vector<std::pair<IndexType, bool>>>> _findsA;
    //     std::map<IndexType, std::map<IndexType, std::vector<std::pair<IndexType, bool>>>> _findsB;
    //};

} // namespace Hammer::MultiDimensional::Ops


#endif
