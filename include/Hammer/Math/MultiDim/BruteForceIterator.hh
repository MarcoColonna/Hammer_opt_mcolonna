///
/// @file  BruteForceIterator.hh
/// @brief Generic tensor indexing iterator
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_MATH_MULTIDIM_BRUTEFORCEITERATOR
#define HAMMER_MATH_MULTIDIM_BRUTEFORCEITERATOR

#include "Hammer/Math/MultiDimensional.fhh"

namespace Hammer::MultiDimensional {

    class BruteForceIterator {
    public:

        BruteForceIterator();
        explicit BruteForceIterator(IndexList dimensions, IndexList fixed = {});

        BruteForceIterator(const BruteForceIterator&) = default;
        BruteForceIterator(BruteForceIterator&&) = default;
        BruteForceIterator& operator=(const BruteForceIterator&) = default;
        BruteForceIterator& operator=(BruteForceIterator&&) = default;

        ~BruteForceIterator() = default;

        [[nodiscard]] BruteForceIterator begin() const;

        [[nodiscard]] BruteForceIterator end() const;

        BruteForceIterator& operator++();

        [[nodiscard]] BruteForceIterator operator++(int n);

        [[nodiscard]] IndexList operator*() const;

        [[nodiscard]] bool isSame(const BruteForceIterator& other) const;

    private:

        void incrementEntry(size_t position, int n);
        void setInitialState();

        IndexList _dimensions;
        IndexList _fixedMask; ///< the strides for each tensor index (necessary to convert coordinates to position
                              ///< in `_data`)
        IndexList _state;
    };


    inline bool operator==(const BruteForceIterator& lhs, const BruteForceIterator& rhs) {
        return lhs.isSame(rhs);
    }

    inline bool operator!=(const BruteForceIterator& lhs, const BruteForceIterator& rhs) {
        return !(lhs == rhs);
    }

} // namespace Hammer::MultiDimensional

namespace std {

    template <>
    struct iterator_traits<Hammer::MultiDimensional::BruteForceIterator> {
        using difference_type = ptrdiff_t;
        using value_type = Hammer::IndexList;
        using pointer = Hammer::IndexList*;
        using reference = Hammer::IndexList&;
        using iterator_category = forward_iterator_tag;
    };

} // namespace std

#endif
