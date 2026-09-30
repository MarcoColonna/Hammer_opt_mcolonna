///
/// @file  BlockIndexing.cc
/// @brief Outer product tensor indexer
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <numeric>

#include "Hammer/Math/MultiDim/BlockIndexing.hh"
#include "Hammer/Exceptions.hh"
#include "Hammer/Math/Utils.hh"

using namespace std;

namespace Hammer::MultiDimensional {

    BlockIndexing::BlockIndexing(const vector<IndexList>& dims, const vector<LabelsList>& labels) {
        auto itd = dims.begin();
        auto itl = labels.begin();
        for (; itd != dims.end(); ++itd, ++itl) {
            _subIndexing.emplace_back(*itd, *itl);
        }
        calc();
    }

    BlockIndexing::BlockIndexing(const LabeledIndexing<AlignedIndexing>& left,
                                 const LabeledIndexing<AlignedIndexing>& right)
        : _subIndexing{left, right} {
        calc();
    }

    size_t BlockIndexing::rank() const {
        return _globalIndexing.rank();
    }

    IndexType BlockIndexing::dim(IndexType index) const {
        try {
            return _globalIndexing.dim(index);
        } catch (Hammer::RangeError&) {
            throw Hammer::RangeError("Index out of range for multidimensional object of rank " + to_string(rank()) +
                                     ": all your base are belong to us.");
        }
    }

    IndexType BlockIndexing::dim(IndexLabel label) const {
        try {
            return _globalIndexing.dim(label);
        } catch (Hammer::IndexLabelError&) {
            throw Hammer::IndexLabelError("Label " + to_string(label) +
                                          " not found for multidimensional object of rank " + to_string(rank()) +
                                          ": all your base are belong to us.");
        }
    }

    IndexList BlockIndexing::dims() const {
        return _globalIndexing.dims();
    }

    LabelsList BlockIndexing::labels() const {
        return _globalIndexing.labels();
    }

    UniqueLabelsList BlockIndexing::spinIndices() const {
        return _globalIndexing.spinIndices();
    }

    PositionType BlockIndexing::numValues() const {
        return _globalIndexing.numValues();
    }

    bool BlockIndexing::checkValidIndices(const IndexList& indices) const {
        return _globalIndexing.checkValidIndices(indices);
    }

    bool BlockIndexing::checkValidIndices(IndexList::const_iterator first, IndexList::const_iterator last) const {
        return _globalIndexing.checkValidIndices(first, last);
    }

    bool BlockIndexing::checkValidIndices(const vector<IndexList>& splits) const {
        auto it = _subIndexing.begin();
        auto its = splits.begin();
        bool result = true;
        for (; it != _subIndexing.end() && result; ++it, ++its) {
            result &= it->checkValidIndices(*its);
        }
        return result;
    }

    vector<IndexList> BlockIndexing::splitIndices(const IndexList& indices) const {
        vector<IndexList> result;
        result.reserve(_subIndexing.size());
        auto ite = indices.begin();
        auto itb = ite;
        for (const auto& it : _subIndexing) {
            itb = ite;
            advance(ite, static_cast<ptrdiff_t>(it.rank()));
            result.emplace_back(itb, ite);
        }
        return result;
    }

    vector<IndexList::const_iterator> BlockIndexing::splitIndices(IndexList::const_iterator first,
                                                                  IndexList::const_iterator last) const {
        vector<IndexList::const_iterator> result;
        result.reserve(_subIndexing.size() + 1);
        auto itv = first;
        for (const auto& elem : _subIndexing) {
            result.push_back(itv);
            advance(itv, static_cast<ptrdiff_t>(elem.rank()));
        }
        result.push_back(last);
        return result;
    }

    IndexPairList BlockIndexing::getSameLabelPairs(const LabelsList& otherLabels, const UniqueLabelsList& indices,
                                                   bool sortedBySecond) const {
        return _globalIndexing.getSameLabelPairs(otherLabels, indices, sortedBySecond);
    }

    size_t BlockIndexing::maxSubRank() const {
        size_t result = 0;
        result = std::accumulate(_subIndexing.begin(), _subIndexing.end(), result,
                                 [](size_t current, const LabeledIndexing<AlignedIndexing>& elem) -> size_t {
                                     return max(current, elem.rank());
                                 });
        return result;
    }

    IndexPairList BlockIndexing::getOppositeLabelPairs(const UniqueLabelsList& indices) const {
        return _globalIndexing.getOppositeLabelPairs(indices);
    }

    IndexType BlockIndexing::labelIndex(IndexLabel label) const {
        return _globalIndexing.labelIndex(label);
    }

    bool BlockIndexing::isSameLabelShape(const LabelsList& otherLabels, const IndexList& otherIndices) const {
        return _globalIndexing.isSameLabelShape(otherLabels, otherIndices);
    }

    bool BlockIndexing::isSameLabelShape(const BlockIndexing& other, bool includeBlockShapes) const {
        if (includeBlockShapes) {
            if (_subIndexing.size() != other._subIndexing.size()) {
                return false;
            }
            auto it1 = _subIndexing.begin();
            auto it2 = other._subIndexing.begin();
            bool res = true;
            for (; it1 != _subIndexing.end() && res; ++it1, ++it2) {
                res &= it1->isSameLabelShape(*it2);
            }
            return res;
        }
        return isSameLabelShape(other.labels(), other.dims());
    }

    IndexPair BlockIndexing::getElementIndex(IndexType position) const {
        auto idx = static_cast<IndexType>(
            distance(_splitIndices.begin(), upper_bound(_splitIndices.begin(), _splitIndices.end(), position)) - 1);
        auto res = static_cast<IndexType>(position - _splitIndices[idx]);
        return make_pair(idx, res);
    }

    const LabeledIndexing<AlignedIndexing>& BlockIndexing::getSubIndexing(IndexType position) const {
        return _subIndexing[position];
    }

    size_t BlockIndexing::numSubIndexing() const {
        return _subIndexing.size();
    }

    void BlockIndexing::calc() {
        LabelsList tmpLabs;
        IndexList tmpIdxs;
        _splitIndices.clear();
        _splitIndices.reserve(_subIndexing.size() + 1);
        _splitIndices.push_back(0);
        IndexType baseIdx = 0;
        _splitPads.clear();
        _splitPads.reserve(_subIndexing.size() + 1);
        _splitPads.push_back(0);
        IndexList tmpPads;
        tmpPads.reserve(_subIndexing.size());
        for (auto& elem : _subIndexing) {
            tmpLabs.insert(tmpLabs.end(), elem.labels().begin(), elem.labels().end());
            tmpIdxs.insert(tmpIdxs.end(), elem.dims().begin(), elem.dims().end());
            baseIdx = static_cast<IndexType>(baseIdx + elem.rank());
            _splitIndices.push_back(baseIdx);
        }
        transform(_subIndexing.rbegin(), _subIndexing.rend(), back_inserter(tmpPads),
                  [](const LabeledIndexing<AlignedIndexing>& val) -> IndexType {
                      return static_cast<IndexType>(minPadding(val.maxIndex() + 1));
                  });
        partial_sum(tmpPads.begin(), tmpPads.end(), back_inserter(_splitPads));
        _splitPads.pop_back();
        reverse(_splitPads.begin(), _splitPads.end());
        _globalIndexing = LabeledIndexing<AlignedIndexing>{tmpIdxs, tmpLabs};
    }

    void BlockIndexing::flipLabels() {
        _globalIndexing.flipLabels();
        for (auto& elem : _subIndexing) {
            elem.flipLabels();
        }
    }

    vector<tuple<IndexList, vector<bool>, PositionType>> BlockIndexing::processShifts(const DotGroupList& chunks,
                                                                                      IndexPairMember which) const {
        vector<tuple<IndexList, vector<bool>, PositionType>> result;
        result.reserve(chunks.size() - 1);
        for (IndexType i = 1; i < static_cast<IndexType>(chunks.size()); ++i) { // first entry is the unassociated;
            auto tmp = _globalIndexing.processShifts(get<2>(chunks[i]), which);
            // now process the others: need to eliminate the paddings for the outer indices belonging to
            // non-contracted subtensors
            IndexList& shifts = get<0>(tmp);
            vector<bool>& outerchecks = get<1>(tmp);
            auto maxPad = static_cast<IndexType>(minPadding(
                get<2>(tmp))); // this is the sum of the pads of all the uncontracted indices (both spectator and non)
            IndexType currentPadSubtract = 0ul;
            IndexType previousPad = 0ul;
            bool recalc = false;
            const IndexList& used = (which == IndexPairMember::Left) ? get<0>(chunks[i]) : get<1>(chunks[i]);
            for (auto j = static_cast<IndexType>(shifts.size()); j-- > 0;) {
                if (!outerchecks[j]) {
                    continue;
                }
                auto subPos = getElementIndex(j);
                bool isUsed = (find(used.begin(), used.end(), subPos.first) != used.end());
                if (recalc) {
                    currentPadSubtract = static_cast<IndexType>(currentPadSubtract + shifts[j] - previousPad);
                    recalc = false;
                }
                if (!isUsed) {
                    previousPad = shifts[j];
                    recalc = true;
                    shifts[j] = maxPad;
                } else {
                    shifts[j] = static_cast<IndexType>(shifts[j] - currentPadSubtract);
                }
            }
            if (recalc) {
                currentPadSubtract = static_cast<IndexType>(currentPadSubtract + maxPad - previousPad);
            }
            get<2>(tmp) = ((get<2>(tmp)) /
                           (1 << currentPadSubtract)); // maxelement for full (i.e 2^(sum of pads of uncontracted
                                                       // indices)) divided by 2^(sum of pads of spectator indices)
            result.push_back(tmp);
        }
        return result;
    }

    PositionType BlockIndexing::splitPosition(const OuterElemIterator& currentPosition, const DotGroupType& chunk,
                                              const IndexList& outerShiftsInnerPositions, const vector<bool>& isOuter,
                                              IndexList& innerList, vector<bool>& innerAdded,
                                              bool shouldCompare) const {
        const IndexList& chunkIndices = shouldCompare ? get<1>(chunk) : get<0>(chunk);
        PositionType fullPos = buildFullPosition(currentPosition, chunkIndices);
        return _globalIndexing.splitPosition(fullPos, outerShiftsInnerPositions, isOuter, innerList, innerAdded,
                                             shouldCompare);
    }


    PositionType BlockIndexing::buildFullPosition(const OuterElemIterator& current,
                                                  const IndexList& chunkIndices) const {
        PositionType fullPos = 0u;
        for (IndexType i = 0; i < static_cast<IndexType>(chunkIndices.size()); ++i) {
            fullPos += (current.isAligned(i) ? current.position(i)
                                             : _subIndexing[chunkIndices[i]].posToAlignedPos(current.position(i)))
                       << _splitPads[chunkIndices[i]];
        }
        return fullPos;
    }

    LabelSignature BlockIndexing::getLabelSignature() const {
        return _globalIndexing.getLabelSignature();
    }

    LabelSignature BlockIndexing::getBlockLabelSignature(IndexType position) const {
        return _subIndexing[position].getLabelSignature();
    }

    OuterElemIterator::OuterElemIterator(const OuterElemIterator& other)
        : _entry{other._entry}, _containers{other._containers}, _dimensions{other._dimensions} {
    }

    OuterElemIterator::OuterElemIterator(const EntryType& entry)
        : _entry{entry}, _containers(entry.size()), _dimensions(entry.size()) {
        for (size_t i = 0; i < _entry.size(); ++i) {
            auto* pTry = dynamic_cast<ISingleContainer*>(entry[i].first.get());
            if (pTry != nullptr) {
                _containers[i] = pTry;
                _dimensions[i] = static_cast<IndexType>(pTry->dataSize() / pTry->entrySize());

            } else {
                throw Error("Invalid contained type");
            }
        }
        setInitialState();
    }

    OuterElemIterator OuterElemIterator::begin() const {
        OuterElemIterator it{*this};
        it.setInitialState();
        return it;
    }

    OuterElemIterator OuterElemIterator::end() const {
        OuterElemIterator it{*this};
        it.setFinalState();
        return it;
    }

    OuterElemIterator& OuterElemIterator::operator++() {
        incrementEntry(_it.size(), 1);
        return *this;
    }

    OuterElemIterator OuterElemIterator::operator++(int n) {
        OuterElemIterator it{*this};
        incrementEntry(_it.size(), n);
        return it;
    }

    IContainer::ElementType OuterElemIterator::operator*() {
        IContainer::ElementType result = 1.;
        for (size_t i = 0; i < _entry.size(); ++i) {
            result *= _entry[i].second ? conj(_it[i]->value()) : _it[i]->value();
        }
        return result;
    }

    PositionType OuterElemIterator::position(IndexType idx) const {
        return _it[idx]->position();
    }

    bool OuterElemIterator::isAligned(IndexType idx) const {
        return _it[idx]->isAligned();
    }

    bool OuterElemIterator::isSame(const OuterElemIterator& other) const {
        ASSERT(&_entry == &(other._entry));
        return _currentIdx == other._currentIdx;
    }

    void OuterElemIterator::setInitialState() {
        _it.resize(_containers.size());
        _currentIdx.assign(_containers.size(), 0);
        for (size_t i = 0; i < _containers.size(); ++i) {
            _it[i] = _containers[i]->firstNonZero();
        }
    }

    void OuterElemIterator::setFinalState() {
        _currentIdx = _dimensions;
    }

    void OuterElemIterator::incrementEntry(size_t position, int n) {
        if (position == _it.size()) {
            if (_currentIdx[0] >= _dimensions[0]) {
                throw RangeError("Invalid increment");
            }
            incrementEntry(_it.size() - 1, n);
            return;
        }
        if (position == 0) {
            if (n + static_cast<long>(_currentIdx[0]) >= static_cast<long>(_dimensions[0])) {
                _currentIdx = _dimensions;
                return;
            }
        }
        auto d1 = static_cast<long>(_currentIdx[position]);
        auto res = div(d1 + n, static_cast<long>(_dimensions[position]));
        _it[position]->next(static_cast<int>(res.rem - d1));
        _currentIdx[position] = static_cast<IndexType>(res.rem);
        if (res.quot != 0) {
            long newpos = static_cast<long>(position) - 1;
            if (newpos < 0) {
                _currentIdx = _dimensions;
                return;
            }
            incrementEntry(static_cast<size_t>(newpos), static_cast<int>(res.quot));
        }
    }

} // namespace Hammer::MultiDimensional
