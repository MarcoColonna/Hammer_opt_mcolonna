///
/// @file  BruteForceIterator.cc
/// @brief Generic tensor indexing iterator
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cstdlib>

#include "Hammer/Math/MultiDim/BruteForceIterator.hh"
#include "Hammer/Exceptions.hh"

using namespace std;

namespace Hammer::MultiDimensional {

    BruteForceIterator::BruteForceIterator() {
        setInitialState();
    }

    BruteForceIterator::BruteForceIterator(IndexList dimensions, IndexList fixed)
        : _dimensions{std::move(dimensions)}, _fixedMask{std::move(fixed)} {
        ASSERT(_fixedMask.empty() || equal(_dimensions.begin(), _dimensions.end(), _fixedMask.begin(), _fixedMask.end(),
                                           [](IndexType a, IndexType b) -> bool { return a >= b; }));
        setInitialState();
    }

    BruteForceIterator BruteForceIterator::begin() const {
        BruteForceIterator b{*this};
        b.setInitialState();
        return b;
    }

    BruteForceIterator BruteForceIterator::end() const {
        BruteForceIterator b{*this};
        b._state = _dimensions;
        return b;
    }

    BruteForceIterator& BruteForceIterator::operator++() {
        incrementEntry(_dimensions.size(), 1);
        return *this;
    }

    BruteForceIterator BruteForceIterator::operator++(int n) {
        BruteForceIterator b{*this};
        incrementEntry(_dimensions.size(), n);
        return b;
    }

    void BruteForceIterator::incrementEntry(size_t position, int n) {
        if (position == _dimensions.size()) {
            if (_dimensions == _state) {
                throw RangeError("Invalid increment");
            }
            long newpos = static_cast<long>(position) - 1;
            if (!_fixedMask.empty()) {
                auto itd = _dimensions.rbegin();
                auto itm = _fixedMask.rbegin();
                while (newpos >= 0 && *itd != *itm) {
                    ++itd;
                    ++itm;
                    --newpos;
                }
                if (newpos < 0) {
                    _state = _dimensions;
                    return;
                }
            }
            incrementEntry(static_cast<size_t>(newpos), n);
            return;
        }
        if (position == 0) {
            if (n + _state[0] >= _dimensions[0]) {
                _state = _dimensions;
                return;
            }
        }
        auto res = div(_state[position] + n, _dimensions[position]);
        _state[position] = static_cast<IndexType>(res.rem);
        if (res.quot != 0) {
            long newpos = static_cast<long>(position) - 1;
            if (!_fixedMask.empty()) {
                while (newpos >= 0 &&
                       _fixedMask[static_cast<size_t>(newpos)] != _dimensions[static_cast<size_t>(newpos)]) {
                    --newpos;
                }
            }
            if (newpos < 0) {
                _state = _dimensions;
                return;
            }
            incrementEntry(static_cast<size_t>(newpos), res.quot);
        }
    }

    IndexList BruteForceIterator::operator*() const {
        return _state;
    }

    void BruteForceIterator::setInitialState() {
        _state.clear();
        if (_fixedMask.empty()) {
            _state = IndexList(_dimensions.size(), 0);
        } else {
            _state.reserve(_dimensions.size());
            auto itd = _dimensions.begin();
            auto itm = _fixedMask.begin();
            for (; itd != _dimensions.end(); ++itd, ++itm) {
                if (*itd == *itm) {
                    _state.push_back(0);
                } else {
                    _state.push_back(*itm);
                }
            }
        }
    }

    bool BruteForceIterator::isSame(const BruteForceIterator& other) const {
        return _state == other._state;
    }


} // namespace Hammer::MultiDimensional
