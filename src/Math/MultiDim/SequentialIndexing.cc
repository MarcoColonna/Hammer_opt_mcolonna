///
/// @file  SequentialIndexing.cc
/// @brief Non-sparse tensor indexer
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <numeric>
#include <functional>
#include <cstdlib>
#include <tuple>
#include <algorithm>

#include "Hammer/Math/MultiDim/SequentialIndexing.hh"
#include "Hammer/Exceptions.hh"

using namespace std;

namespace Hammer::MultiDimensional {

    SequentialIndexing::SequentialIndexing(IndexList dimensions) : _dimensions{std::move(dimensions)}, _maxIndex{} {
        calcPadding();
    }

    size_t SequentialIndexing::rank() const {
        return _dimensions.size();
    }

    IndexType SequentialIndexing::dim(IndexType index) const {
        if (index < _dimensions.size()) {
            return _dimensions[index];
        }
        throw Hammer::RangeError("Index out of range for multidimensional object of rank " +
                                 to_string(_dimensions.size()) + ": all your base are belong to us.");
    }

    const IndexList& SequentialIndexing::dims() const {
        return _dimensions;
    }

    PositionType SequentialIndexing::numValues() const {
        return _maxIndex + 1;
    }

    bool SequentialIndexing::checkValidIndices(const IndexList& indices) const {
        bool valid = (indices.size() == _dimensions.size());
        if (!valid) {
            return false;
        }
        valid = inner_product(indices.begin(), indices.end(), _dimensions.begin(), valid, logical_and<>(), less<>());
        return valid;
    }

    bool SequentialIndexing::checkValidIndices(IndexList::const_iterator first, IndexList::const_iterator last) const {
        bool valid = (distance(first, last) == static_cast<ptrdiff_t>(_dimensions.size()));
        if (!valid) {
            return false;
        }
        valid = inner_product(first, last, _dimensions.begin(), valid, logical_and<>(), less<>());
        return valid;
    }

    PositionType SequentialIndexing::indicesToPos(const IndexList& indices) const {
        if (rank() == 1) {
            return indices[0];
        }
        return indicesToPos(indices.begin(), indices.end());
    }

    PositionType SequentialIndexing::indicesToPos(IndexList::const_iterator first,
                                                  IndexList::const_iterator last) const {
        if (rank() == 1) {
            return *first;
        }
        PositionType index = 0ul;
        function<size_t(PositionType, IndexType)> opCombine = [](PositionType a, IndexType b) -> PositionType {
            return a * b;
        };
        index = inner_product(first, last, _strides.begin(), index, plus<>(), opCombine);
        return index;
    }

    void SequentialIndexing::posToIndices(PositionType position, IndexList& result) const {
        if (rank() == 1) {
            result[0] = static_cast<IndexType>(position);
        }
        pair<IndexType, PositionType> current{'\0', position};
        auto its = _strides.begin();
        auto itd = _dimensions.begin();
        auto itr = result.begin();
        for (; its != _strides.end(); ++its, ++itd, ++itr) {
            current = {static_cast<IndexType>(current.second / (*its)), current.second % (*its)};
            if (static_cast<PositionType>(current.first) < *itd) {
                *itr = current.first;
            } else {
                throw RangeError("Index #" + to_string(distance(itd, _dimensions.begin())) +
                                 " for multidimensional object "
                                 " of rank " +
                                 to_string(_dimensions.size()) + " out of range (" + to_string(current.first) + " vs " +
                                 to_string(*itd) + "): all your base are belong to us.");
            }
        }
    }

    IndexType SequentialIndexing::ithIndexInPos(PositionType position, IndexType indexPosition) const {
        if (indexPosition == '\0') {
            return static_cast<IndexType>(position / _strides[indexPosition]);
        }
        PositionType temp = position % _strides[indexPosition - 1];
        return static_cast<IndexType>(temp / _strides[indexPosition]);
    }

    PositionType SequentialIndexing::extendPosition(PositionType position, IndexType indexPosition,
                                                    IndexType indexValue) const {
        PositionType strideEntry = _strides[indexPosition];
        auto res = div(static_cast<long>(position), static_cast<long>(strideEntry));
        return (static_cast<PositionType>(res.quot) * _dimensions[indexPosition] + indexValue) * strideEntry +
               static_cast<PositionType>(res.rem);
    }


    PositionType SequentialIndexing::build2ndPosition(PositionType reducedPosition, PositionType innerPosition,
                                                      const PositionPairList& conversion) {
        PositionType out = 0ul;
        PositionType tmp = reducedPosition;
        for (const auto& elem : conversion) {
            auto res = div(static_cast<long>(tmp), static_cast<long>(elem.first));
            out += static_cast<PositionType>(res.quot) * elem.second;
            tmp = static_cast<PositionType>(res.rem);
        }
        return out + innerPosition;
    }

    PositionType SequentialIndexing::reducedNumValues(const IndexPairList& indices) const {
        const auto& dims = _dimensions;
        auto tmp = accumulate(indices.begin(), indices.end(), 1ul,
                              [dims](PositionType current, const IndexPair& elem) -> PositionType {
                                  return current * dims[elem.second];
                              });
        return numValues() / tmp;
    }

    PositionPair SequentialIndexing::splitPosition(PositionType position, const StrideMap& conversion) {
        PositionType out = 0ul;
        long in = 0ul;
        PositionType tmp = position;
        for (const auto& elem : conversion) {
            if (!get<2>(elem) && get<0>(elem) == 0) {
                continue;
            }
            auto res = div(static_cast<long>(tmp), static_cast<long>(get<0>(elem)));
            if (get<2>(elem)) {
                in += res.quot * get<1>(elem);
            } else {
                out += static_cast<PositionType>(res.quot * get<1>(elem));
            }
            tmp = static_cast<PositionType>(res.rem);
        }
        return make_pair(out, in);
    }

    SequentialIndexing::StrideMap SequentialIndexing::getInnerOuterStrides(const IndexPairList& positions,
                                                                           const PositionList& secondStrides,
                                                                           bool flipSecond) const {
        ASSERT(!positions.empty() &&
               positions.size() <= _dimensions.size()); // NOLINT(readability-simplify-boolean-expr)
        ASSERT(is_sorted(positions.begin(), positions.end(),
                         [](const pair<IndexType, IndexType>& a, const pair<IndexType, IndexType>& b) -> bool {
                             return a.second < b.second;
                         }));
        map<IndexType, long> innerStrides;
        for (const auto& elem : positions) {
            innerStrides.insert({elem.first, secondStrides[elem.second]});
            if (flipSecond) {
                innerStrides.insert({elem.second, -1 * static_cast<long>(secondStrides[elem.second])});
            }
        }
        return buildStrideMap(innerStrides);
    }

    SequentialIndexing::StrideMap SequentialIndexing::buildStrideMap(const map<IndexType, long>& innerMap) const {
        auto it = innerMap.end();
        PositionType outerStride = 1ul;
        StrideMap results;
        results.reserve(_dimensions.size());
        bool done = (it == innerMap.begin());
        if (!done) {
            --it;
        }
        for (auto i = static_cast<IndexType>(_dimensions.size()); i-- != 0;) {
            if (!done && i == it->first) {
                results.emplace_back(_strides[i], it->second, true);
                if (it == innerMap.begin()) {
                    done = true;
                } else {
                    --it;
                }
            } else {
                results.emplace_back(_strides[i], outerStride, false);
                outerStride *= _dimensions[i];
            }
        }
        reverse(results.begin(), results.end());
        for (size_t i = 0; i < results.size() - 1; ++i) {
            if (!get<2>(results[i + 1]) && !get<2>(results[i])) {
                results[i] = make_tuple(0, 0, false);
            }
        }
        return results;
    }

    PositionPairList SequentialIndexing::getOuterStrides2nd(const IndexPairList& positions) const {
        ASSERT(!positions.empty() &&
               positions.size() <= _dimensions.size()); // NOLINT(readability-simplify-boolean-expr)
        ASSERT(is_sorted(positions.begin(), positions.end(),
                         [](const IndexPair& a, const IndexPair& b) -> bool { return a.second < b.second; }));
        auto it = positions.end();
        PositionType outerStride = 1ul;
        PositionPairList results;
        results.reserve(_dimensions.size() - positions.size());
        bool done = (it == positions.begin());
        if (!done) {
            --it;
        }
        for (auto i = static_cast<IndexType>(_dimensions.size()); i-- != 0;) {
            if (!done && i == it->second) {
                if (it == positions.begin()) {
                    done = true;
                } else {
                    --it;
                }
            } else {
                results.emplace_back(outerStride, _strides[i]);
                outerStride *= _dimensions[i];
            }
        }
        reverse(results.begin(), results.end());
        return results;
    }


    PositionType SequentialIndexing::stride(IndexType index) const {
        return _strides[index];
    }

    const PositionList& SequentialIndexing::strides() const {
        return _strides;
    }

    void SequentialIndexing::calcPadding() {
        _strides.clear();
        _strides.push_back(1);
        partial_sum(_dimensions.rbegin(), _dimensions.rend(), back_inserter(_strides), multiplies<>());
        _maxIndex = _strides.back() - 1;
        _strides.pop_back();
        reverse(_strides.begin(), _strides.end());
    }


    bool SequentialIndexing::isSameShape(const IndexList& indices) const {
        if (rank() != indices.size()) {
            return false;
        }
        return equal(_dimensions.begin(), _dimensions.end(), indices.begin());
    }

} // namespace Hammer::MultiDimensional
