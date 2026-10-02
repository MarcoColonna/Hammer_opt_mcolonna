///
/// @file  AlignedIndexing.cc
/// @brief Sparse tensor indexer
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <numeric>
#include <functional>
#include <algorithm>

#include "Hammer/Math/MultiDim/AlignedIndexing.hh"
#include "Hammer/Exceptions.hh"
#include "Hammer/Math/Utils.hh"

using namespace std;

namespace Hammer::MultiDimensional {

    static IndexType nextPowerOf2(IndexType value) {
        if (value == 0) {
            return value;
        }
        if (value == 1) {
            return 0;
        }
        --value;
        value = static_cast<IndexType>(value | (value >> 1));
        value = static_cast<IndexType>(value | (value >> 2));
        value = static_cast<IndexType>(value | (value >> 4));
        value = static_cast<IndexType>(value | (value >> 8));
        ++value;
        return value;
    }

    AlignedIndexing::AlignedIndexing(IndexList dimensions)
        : _dimensions{std::move(dimensions)}, _maxIndex{0}, _maxAlignedIndex{0} {
        calc();
    }

    size_t AlignedIndexing::rank() const {
        return _dimensions.size();
    }

    IndexType AlignedIndexing::dim(IndexType index) const {
        if (index < _dimensions.size()) {
            return _dimensions[index];
        }
        throw Hammer::RangeError("Index out of range for multidimensional object of rank " +
                                 to_string(_dimensions.size()) + ": all your base are belong to us.");
    }

    const IndexList& AlignedIndexing::dims() const {
        return _dimensions;
    }

    PositionType AlignedIndexing::numValues() const {
        return _maxIndex + 1;
    }

    bool AlignedIndexing::checkValidIndices(const IndexList& indices) const {
        bool valid = (indices.size() == _dimensions.size());
        if (!valid) {
            return false;
        }
        valid = inner_product(indices.begin(), indices.end(), _dimensions.begin(), valid, logical_and<>(), less<>());
        return valid;
    }

    bool AlignedIndexing::checkValidIndices(IndexList::const_iterator first, IndexList::const_iterator last) const {
        bool valid = (distance(first, last) == static_cast<ptrdiff_t>(_dimensions.size()));
        if (!valid) {
            return false;
        }
        valid = inner_product(first, last, _dimensions.begin(), valid, logical_and<>(), less<>());
        return valid;
    }

    PositionType AlignedIndexing::indicesToPos(const IndexList& indices) const {
        return indicesToPos(indices.begin(), indices.end());
    }

    PositionType AlignedIndexing::indicesToPos(IndexList::const_iterator first, IndexList::const_iterator last) const {
        if (rank() == 1) {
            return *first;
        }
        PositionType index = 0ul;
        function<size_t(PositionType, IndexType)> opCombine = [](PositionType idx, IndexType pad) -> PositionType {
            return idx << pad;
        };
        index = inner_product(first, last, _alignPads.begin(), index, plus<>(), opCombine);
        return index;
    }

    void AlignedIndexing::posToIndices(PositionType alignedPosition, IndexList& result) const {
        ASSERT(alignedPosition <= _maxAlignedIndex);
        result.clear();
        auto itp = _alignPads.begin();
        auto itm = _alignMasks.begin();
        for (; itp != _alignPads.end(); ++itp, ++itm) {
            result.push_back(static_cast<IndexType>((alignedPosition >> (*itp)) & (*itm)));
        }
    }

    IndexType AlignedIndexing::ithIndexInPos(PositionType alignedPosition, IndexType indexPosition) const {
        ASSERT(indexPosition < _alignPads.size());
        return static_cast<IndexType>((alignedPosition >> _alignPads[indexPosition]) & _alignMasks[indexPosition]);
    }

    PositionType AlignedIndexing::splitPosition(PositionType alignedPosition,
                                                const IndexList& outerShiftsInnerPositions, const FlagList& isOuter,
                                                IndexList& innerList, FlagList& innerAdded,
                                                bool shouldCompare) const {
        fill(innerAdded.begin(), innerAdded.end(), false);
        PositionType newPos = 0ul;
        for (size_t i = 0; i < _dimensions.size(); ++i) {
            auto idx = static_cast<IndexType>((alignedPosition >> _alignPads[i]) & _alignMasks[i]);
            if (isOuter[i]) {
                newPos += static_cast<PositionType>(idx << outerShiftsInnerPositions[i]);
            } else {
                ASSERT(outerShiftsInnerPositions[i] < innerList.size());
                if (!shouldCompare) {
                    size_t innerPos = outerShiftsInnerPositions[i];
                    if (innerAdded[innerPos] && innerList[innerPos] != idx) {
                        return numeric_limits<size_t>::max();
                    }
                    innerAdded[innerPos] = true;
                    innerList[innerPos] = idx;
                } else if (innerList[outerShiftsInnerPositions[i]] != idx) {
                    return numeric_limits<size_t>::max();
                }
            }
        }
        return newPos;
    }

    std::tuple<IndexList, FlagList, PositionType> AlignedIndexing::processShifts(const IndexPairList& pairs,
                                                                                          IndexPairMember which) const {
        map<IndexType, IndexType> inners;
        for (IndexType i = 0; i < static_cast<IndexType>(pairs.size()); ++i) {
            switch (which) {
            case IndexPairMember::Left:
                inners.insert(make_pair(pairs[i].first, i));
                break;
            case IndexPairMember::Right:
                inners.insert(make_pair(pairs[i].second, i));
                break;
            case IndexPairMember::Both:
                inners.insert(make_pair(pairs[i].first, i));
                inners.insert(make_pair(pairs[i].second, i));
                break;
            }
        }
        IndexList newDims;
        newDims.reserve(_dimensions.size());
        for (IndexType i = 0; i < static_cast<IndexType>(_dimensions.size()); ++i) {
            newDims.push_back((inners.find(i) != inners.end()) ? '\0' : _dimensions[i]);
        }
        IndexList retIndices;
        FlagList retBools(_dimensions.size(), true);
        auto maxIdx = calcPadding(newDims, retIndices);
        for (auto elem : inners) {
            retIndices[elem.first] = elem.second;
            retBools[elem.first] = false;
        }
        return make_tuple(retIndices, retBools, maxIdx + 1);
    }

    PositionType AlignedIndexing::maxIndex(bool aligned) const {
        return aligned ? _maxAlignedIndex : _maxIndex;
    }

    PositionType AlignedIndexing::calcPadding(const IndexList& dimensions, IndexList& pads) {
        IndexList tmpPads;
        tmpPads.reserve(dimensions.size());
        transform(dimensions.rbegin(), dimensions.rend(), back_inserter(tmpPads),
                  [](IndexType val) -> IndexType { return minPadding(val); });
        pads.clear();
        pads.push_back(0);
        partial_sum(tmpPads.begin(), tmpPads.end(), back_inserter(pads));
        PositionType maxAlignedIndex = (1ul << pads.back()) - 1ul;
        pads.pop_back();
        reverse(pads.begin(), pads.end());
        return maxAlignedIndex;
    }

    void AlignedIndexing::calcMasks(const IndexList& dimensions, IndexList& masks) {
        masks.clear();
        transform(dimensions.rbegin(), dimensions.rend(), back_inserter(masks), [](IndexType val) -> IndexType {
            return (val == 1u) ? 0u : static_cast<IndexType>(nextPowerOf2(val) - 1u);
        });
        reverse(masks.begin(), masks.end());
    }

    void AlignedIndexing::calc() {
        if (_dimensions.empty()) {
            _maxIndex = 0;
            _maxAlignedIndex = 0;
            _alignMasks.clear();
            _alignPads.clear();
            _unalignedEntries.clear();
            return;
        }
        _maxIndex = accumulate(_dimensions.begin(), _dimensions.end(), 1ul, multiplies<>()) - 1ul;
        _maxAlignedIndex = calcPadding(_dimensions, _alignPads);
        calcMasks(_dimensions, _alignMasks);
        calcUnaligned(_dimensions, _alignPads, _unalignedEntries);
    }

    void AlignedIndexing::calcUnaligned(const IndexList& dimensions, const IndexList& pads,
                                        PosIndexPairList& unaligned) {
        unaligned.clear();
        auto itd = dimensions.rbegin();
        auto itp = pads.rbegin();
        ++itp;
        PositionType currentStride = 1;
        for (; itp != pads.rend(); ++itd, ++itp) {
            currentStride *= *itd;
            if (*itd != 1 && nextPowerOf2(*itd) != *itd) {
                unaligned.emplace_back(currentStride, *itp);
            }
        }
        reverse(unaligned.begin(), unaligned.end());
    }

    PositionType AlignedIndexing::alignedPosToPos(PositionType alignedPosition) const {
        PositionType result = alignedPosition;
        PositionType reminder = result;
        for (const auto& elem : _unalignedEntries) {
            PositionType tmp = reminder >> elem.second;
            reminder -= (tmp << elem.second);
            result += tmp * elem.first - (tmp << elem.second);
        }
        return result;
    }

    PositionType AlignedIndexing::posToAlignedPos(PositionType position) const {
        PositionType result = position;
        PositionType reminder = result;
        for (const auto& elem : _unalignedEntries) {
            PositionType tmp = reminder / elem.first;
            reminder -= (tmp * elem.first);
            result += (tmp << elem.second) - (tmp * elem.first);
        }
        return result;
    }

    PositionType AlignedIndexing::extendAlignedPosition(PositionType alignedPosition, IndexType indexPosition,
                                                        IndexType indexValue) const {
        if (indexPosition == 0) {
            return static_cast<PositionType>(indexValue << _alignPads[0]) + alignedPosition;
        }
        PositionType tmp1 = alignedPosition >> _alignPads[indexPosition];
        PositionType tmp2 = alignedPosition - (tmp1 << _alignPads[indexPosition]);
        return (tmp1 << _alignPads[indexPosition - 1]) +
               static_cast<PositionType>(indexValue << _alignPads[indexPosition]) + tmp2;
    }

    PositionType AlignedIndexing::extendPosition(PositionType position, PositionType stride, IndexType indexPosition,
                                                 IndexType indexValue) const {
        auto res = div(static_cast<long>(position), static_cast<long>(stride));
        return posToAlignedPos((static_cast<PositionType>(res.quot * _dimensions[indexPosition]) + indexValue) *
                                   stride +
                               static_cast<PositionType>(res.rem));
    }

    bool AlignedIndexing::isSameShape(const IndexList& indices) const {
        if (rank() != indices.size()) {
            return false;
        }
        return equal(_dimensions.begin(), _dimensions.end(), indices.begin());
    }

} // namespace Hammer::MultiDimensional
