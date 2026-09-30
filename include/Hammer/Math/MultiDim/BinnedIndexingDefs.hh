#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
///
/// @file  BinnedIndexingDefs.hh
/// @brief Binned tensor (histogram) indexer template method definitions
///
#pragma clang diagnostic pop

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-

#include "Hammer/Exceptions.hh"
#include "Hammer/Math/Utils.hh"

namespace Hammer::MultiDimensional {

    template <class BasicIndexing>
    BinnedIndexing<BasicIndexing>::BinnedIndexing() : BasicIndexing{}, _hasUnderOverFlow{false} {
    }

    template <class BasicIndexing>
    BinnedIndexing<BasicIndexing>::BinnedIndexing(const IndexList& dimensions, const BinEdgeList& edges,
                                                  bool hasUnderOverFlow)
        : BasicIndexing{dimensions}, _edges{edges}, _hasUnderOverFlow{hasUnderOverFlow} {
        if (!isValid()) {
            throw RangeError("Length of list of bin edges is incompatible with length of dimensions: " +
                             std::to_string(edges.size()) + " vs " + std::to_string(dimensions.size()) +
                             " -> all your base are belong to us.");
        }
    }

    template <class BasicIndexing>
    BinnedIndexing<BasicIndexing>::BinnedIndexing(const IndexList& dimensions, const BinRangeList& ranges,
                                                  bool hasUnderOverFlow)
        : BasicIndexing{BinnedIndexing<BasicIndexing>::processRanges(dimensions, hasUnderOverFlow)},
          _hasUnderOverFlow{hasUnderOverFlow} {
        fillEdges(ranges);
    }

    template <class BasicIndexing>
    BinnedIndexing<BasicIndexing>::BinnedIndexing(const BinEdgeList& edges, bool hasUnderOverFlow)
        : BasicIndexing{BinnedIndexing<BasicIndexing>::processEdges(edges, hasUnderOverFlow)}, _edges{edges},
          _hasUnderOverFlow{hasUnderOverFlow} {
    }

    template <class BasicIndexing>
    const BinEdgeList& BinnedIndexing<BasicIndexing>::edges() const {
        return _edges;
    }

    template <class BasicIndexing>
    const BinValue& BinnedIndexing<BasicIndexing>::edge(IndexType pos) const {
        ASSERT(pos < _edges.size());
        return _edges[pos];
    }

    template <class BasicIndexing>
    IndexList BinnedIndexing<BasicIndexing>::valueToPos(const BinValue& point) const {
        auto itv = point.begin();
        auto itb = _edges.begin();
        IndexList result;
        int delta = _hasUnderOverFlow ? 0 : -1;
        for (; itb != _edges.end(); ++itv, ++itb) {
            auto itp = lower_bound(itb->begin(), itb->end(), *itv);
            result.push_back(static_cast<IndexType>(distance(itb->begin(), itp) + delta));
        }
        return result;
    }

    template <class BasicIndexing>
    bool BinnedIndexing<BasicIndexing>::isValid(const BinValue& point) const {
        return (_hasUnderOverFlow || this->checkValidIndices(valueToPos(point)));
    }

    template <class BasicIndexing>
    BinRangeList BinnedIndexing<BasicIndexing>::binEdges(IndexList position) const {
        ASSERT(this->checkValidIndices(position));
        BinRangeList result;
        result.reserve(position.size());
        auto itp = position.begin();
        auto itd = this->dims().begin();
        auto ite = _edges.begin();
        int delta = _hasUnderOverFlow ? -1 : 0;
        for (; itp != position.end(); ++itp, ++ite, ++itd) {
            if (_hasUnderOverFlow && *itp == 0) {
                result.emplace_back(-std::numeric_limits<double>::max(), (*ite)[0]);
            } else if (_hasUnderOverFlow && *itp == ((*itd) - 1)) {
                result.emplace_back(ite->back(), std::numeric_limits<double>::max());
            } else {
                result.emplace_back((*ite)[*itp + delta], (*ite)[*itp + delta + 1]);
            }
        }
        return result;
    }

    template <class BasicIndexing>
    BinRange BinnedIndexing<BasicIndexing>::binEdge(IndexType position, IndexType coord) const {
        ASSERT(coord < this->rank());
        ASSERT(position < this->dims(coord));
        const auto& binnings = _edges[coord];
        int delta = _hasUnderOverFlow ? -1 : 0;
        if (_hasUnderOverFlow && position == 0) {
            return {-std::numeric_limits<double>::max(), binnings[0]};
        }
        if (_hasUnderOverFlow && position == (this->dims(coord) - 1)) {
            return {binnings.back(), std::numeric_limits<double>::max()};
        }
        return {binnings[position + delta], binnings[position + delta + 1]};
    }


    template <class BasicIndexing>
    bool BinnedIndexing<BasicIndexing>::isSameBinShape(const BinEdgeList& otherEdges, const IndexList& otherIndices,
                                                       bool otherUnderOverFlow) const {
        bool result = (_hasUnderOverFlow == otherUnderOverFlow);
        if (!result) {
            return false;
        }
        result = this->isSameShape(otherIndices);
        if (!result) {
            return false;
        }
        result &= equal(_edges.begin(), _edges.end(), otherEdges.begin());
        return result;
    }

    template <class BasicIndexing>
    bool BinnedIndexing<BasicIndexing>::isSameBinShape(const BinRangeList& otherRanges, const IndexList& otherIndices,
                                                       bool otherUnderOverFlow) const {
        bool result = (_hasUnderOverFlow == otherUnderOverFlow);
        if (!result) {
            return false;
        }
        result = this->isSameShape(otherIndices);
        if (!result) {
            return false;
        }
        auto ite = _edges.begin();
        auto itr = otherRanges.begin();
        for (; ite != _edges.end() && itr != otherRanges.end(); ++ite, ++itr) {
            result &= isZero(ite->front() - itr->first) && isZero(ite->back() - itr->second);
            if (!result) {
                break;
            }
        }
        return result;
    }

    template <class BasicIndexing>
    template <typename S>
    bool BinnedIndexing<BasicIndexing>::isSameBinShape(const BinnedIndexing<S>& other) const {
        bool result = (_hasUnderOverFlow == other._hasUnderOverFlow);
        if (!result) {
            return false;
        }
        result = this->isSameShape(other);
        if (!result) {
            return false;
        }
        result &= equal(_edges.begin(), _edges.end(), other._edges.begin());
        return result;
    }

    template <class BasicIndexing>
    bool BinnedIndexing<BasicIndexing>::isValid() const {
        bool valid = (_edges.size() == this->rank());
        for (IndexType i = 0; i < static_cast<IndexType>(_edges.size()); ++i) {
            valid &= ((_hasUnderOverFlow && (_edges[i].size() + 1u == dim(i))) ||
                      (!_hasUnderOverFlow && (_edges[i].size() == dim(i) + 1u)));
        }
        return valid;
    }

    template <class BasicIndexing>
    void BinnedIndexing<BasicIndexing>::fillEdges(const BinRangeList& ranges) {
        int deltaIdx = _hasUnderOverFlow ? -1 : 1;
        if (!ranges.empty()) {
            ASSERT(ranges.size() == this->rank());
            auto itr = ranges.begin();
            auto itn = this->dims().begin();
            for (; itn != this->dims().end(); ++itr, ++itn) {
                std::vector<double> tmpvec(static_cast<size_t>(*itn + deltaIdx));
                const double delta = (itr->second - itr->first) / (1. * (*itn + deltaIdx - 1));
                generate(tmpvec.begin(), tmpvec.end(), [&, n = itr->first]() mutable {
                    double tmp = n;
                    n += delta;
                    return tmp;
                });
                _edges.push_back(tmpvec);
            }
        } else {
            for (auto itn = this->dims().begin(); itn != this->dims().end(); ++itn) {
                std::vector<double> tmpvec(static_cast<size_t>(*itn + deltaIdx));
                generate(tmpvec.begin(), tmpvec.end(), [n = 0]() mutable { return n++; });
                _edges.push_back(tmpvec);
            }
        }
    }

    template <class BasicIndexing>
    IndexList BinnedIndexing<BasicIndexing>::processEdges(const BinEdgeList& edges, bool hasUnderOverFlow) {
        IndexList result;
        result.reserve(edges.size());
        int delta = hasUnderOverFlow ? 1 : -1;
        for (const auto& elem : edges) {
            ASSERT(static_cast<int>(elem.size()) + delta > 0);
            result.push_back(static_cast<IndexType>(static_cast<int>(elem.size()) + delta));
        }
        return result;
    }

    template <class BasicIndexing>
    IndexList BinnedIndexing<BasicIndexing>::processRanges(const IndexList& dimensions, bool hasUnderOverFlow) {
        IndexList result;
        result.reserve(dimensions.size());
        int delta = hasUnderOverFlow ? 2 : 0;
        for (const auto& elem : dimensions) {
            ASSERT(elem + delta > 0);
            result.push_back(static_cast<IndexType>(elem + delta));
        }
        return result;
    }


    template <class BasicIndexing>
    bool BinnedIndexing<BasicIndexing>::hasUnderOverFlow() const {
        return _hasUnderOverFlow;
    }

    template <class BasicIndexing>
    typename BasicIndexing::StrideMap
    BinnedIndexing<BasicIndexing>::getBinStrides(const UniqueIndexList& positions) const {
        std::map<IndexType, long> innerStrides;
        for (const auto& elem : positions) {
            innerStrides.insert({elem, 0});
        }
        return this->buildStrideMap(innerStrides);
    }

} // namespace Hammer::MultiDimensional
