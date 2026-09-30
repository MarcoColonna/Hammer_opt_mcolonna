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
         private:
             Hammer::IndexList _state;
             Hammer::IndexList _first;
             Hammer::IndexList _last;

             BruteForceIterator(const Hammer::IndexList& state,
		                const Hammer::IndexList& first,
				const Hammer::IndexList& last);
             static Hammer::IndexList build_first(
		     const Hammer::IndexList& dimensions,
		     const Hammer::IndexList& fixed);

         public:
             explicit BruteForceIterator(Hammer::IndexList dimensions,
		                         Hammer::IndexList fixed = {});

             BruteForceIterator begin() const;
             BruteForceIterator end() const;

             inline const Hammer::IndexList& operator*() const noexcept
	     { return _state; }

             BruteForceIterator& operator++() noexcept;
             BruteForceIterator operator++(int /* unused */);

             friend bool operator==(const BruteForceIterator& a,
                                    const BruteForceIterator& b) noexcept;
             friend inline bool operator!=(const BruteForceIterator& a,
                                           const BruteForceIterator& b) noexcept
             { return !(a == b); }

             friend bool operator<(const BruteForceIterator& a,
                                   const BruteForceIterator& b) noexcept;
             friend inline bool operator>(const BruteForceIterator& a,
                                          const BruteForceIterator& b) noexcept
             { return b < a; }
             friend inline bool operator<=(const BruteForceIterator& a,
                                           const BruteForceIterator& b) noexcept
             { return !(b < a); }
             friend inline bool operator>=(const BruteForceIterator& a,
                                           const BruteForceIterator& b) noexcept
             { return !(a < b); }
         };
    } // namespace MultiDimensional

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
