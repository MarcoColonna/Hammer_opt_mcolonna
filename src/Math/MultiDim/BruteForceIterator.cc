///
/// @file  BruteForceIterator.cc
/// @brief Generic tensor indexing iterator
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-

#include <algorithm>
#include <numeric>
#include <utility>
#include <stdexcept>

#include "Hammer/Math/MultiDim/BruteForceIterator.hh"
#include "Hammer/Exceptions.hh"

// if compiler supports it, give hints on branch probabilities (for use in
// BruteForceIterator::operator++())
#if defined(__GNUC__)
#define LIKELY(x) __builtin_expect((x), 1)
#else
#define LIKELY(x) (x)
#endif

namespace Hammer::MultiDimensional {

    IndexList BruteForceIteratorRange::build_first(const IndexList& dimensions, const IndexList& fixed)
        {
            // build first state in sequence
            IndexList retVal;
            retVal.reserve(dimensions.size());
            if (fixed.empty()) {
                // if no fixed dimensions, first in sequence is just all zero 
                retVal.assign(dimensions.size(), 0);
            } else {
                // fixed dimensions
                //
                // start with some basic validation of input arguments
                if (dimensions.size() != fixed.size())
                    throw std::logic_error("arguments must have same size");
                if (!std::inner_product(
                            fixed.begin(), fixed.end(), dimensions.begin(), true,
                            [](bool a, bool b) { return a && b; },
                            [](const auto& a, const auto& b) { return a <= b; }))
                    throw std::logic_error(
                            "fixed elements must be <= dimensions elements");
                // replace the non-fixed dimensions with zero, and keep the fixed
                // dimensions as indicated by fixed
                auto it = dimensions.begin();
                std::replace_copy_if(
                        fixed.begin(), fixed.end(), std::back_inserter(retVal),
                        [&it](const auto& el) { return el == *it++; }, 0);
            }
            return retVal;
    }

    BruteForceIteratorRange::BruteForceIteratorRange(IndexList dimensions, IndexList fixed)
                : _first{build_first(dimensions, fixed)},
                  _last{fixed.empty() ? std::move(dimensions) : std::move(fixed)}
    {}

    BruteForceIterator BruteForceIteratorRange::begin() const
        {
            return {_first, *this};
        }

    BruteForceIterator BruteForceIteratorRange::end() const
        {
            return {_last, *this};
        }


    BruteForceIterator& BruteForceIterator::operator++() noexcept
        {
            auto i = _state.size();
            if (LIKELY(i)) {

                while (true) {
                    --i;
                    ++_state[i];
                    if (LIKELY(_state[i] < _parent._last[i])) {
                        break;
                    } else {
                        if (LIKELY(i)) {
                            std::copy(_parent._first.begin() + i, _parent._first.end(),
                                      _state.begin() + i);
                        } else {
                            std::copy(_parent._last.begin(), _parent._last.end(), _state.begin());
                            break;
                        }
                    }
                }
            }
            return *this;
        }

    BruteForceIterator BruteForceIterator::operator++(int /* unused */)
        {
            const auto retVal{*this};
            operator++();
            return retVal;

        }

    IndexList BruteForceIterator::operator*() const {
        return _state;
    }

    bool operator==(const BruteForceIterator& a, const BruteForceIterator& b) noexcept
    {
        return a._state == b._state;
    }

    bool operator<(const BruteForceIterator& a, const BruteForceIterator& b) noexcept
    {
        return std::lexicographical_compare(a._state.begin(), a._state.end(), b._state.begin(), b._state.end());
    }

} // namespace Hammer::MultiDimensional
