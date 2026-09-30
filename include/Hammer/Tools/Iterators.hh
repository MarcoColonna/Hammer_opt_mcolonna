///
/// @file  Iterators.hh
/// @brief Iterator to cycle over outer products of containers
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_TOOLS_ITERATORS
#define HAMMER_TOOLS_ITERATORS

#include <utility>
#include <iterator>
#include <type_traits>

namespace Hammer {

    template <typename OuterIterator, typename InnerIterator>
    class product_iterator {
    public:

        using outer_iterator = OuterIterator;
        using inner_iterator = InnerIterator;

        using outer_category = typename std::iterator_traits<outer_iterator>::iterator_category;
        using inner_category = typename std::iterator_traits<inner_iterator>::iterator_category;
        using common_category = typename std::common_type<outer_category, inner_category>::type;

        using iterator_category =
            typename std::conditional<std::is_same<common_category, std::random_access_iterator_tag>::value,
                                      std::bidirectional_iterator_tag, common_category>::type;

        using inner_value_type = typename std::iterator_traits<inner_iterator>::value_type;
        using outer_value_type = typename std::iterator_traits<outer_iterator>::value_type;
        using value_type = typename std::pair<outer_value_type, outer_value_type>;
        // using difference_type = typename std::iterator_traits<inner_iterator>::difference_type;
        using inner_reference = typename std::iterator_traits<inner_iterator>::reference;
        using outer_reference = typename std::iterator_traits<outer_iterator>::reference;
        using inner_pointer = typename std::iterator_traits<inner_iterator>::pointer;
        using outer_pointer = typename std::iterator_traits<outer_iterator>::pointer;
        using pointer = typename std::pair<outer_pointer, inner_pointer>;
        using reference = typename std::pair<outer_reference, inner_reference>;

        product_iterator() {
        }
        product_iterator(outer_iterator out_it, outer_iterator out_begin, outer_iterator out_end, inner_iterator in_it,
                         inner_iterator in_begin, inner_iterator in_end)
            : outer_it_(out_it), inner_it_(in_it), outer_begin_(out_begin), outer_end_(out_end), inner_begin_(in_begin),
              inner_end_(in_end) {
            if (outer_begin_ == outer_end_) {
                return;
            }
            if (inner_begin_ == inner_end_) {
                return;
            }

            if (outer_it_ == outer_end_) {
                return;
            }
            if (inner_it_ == outer_end_) {
                return;
            }
        }

        product_iterator(outer_iterator out_begin, outer_iterator out_end, inner_iterator in_begin,
                         inner_iterator in_end)
            : outer_it_(out_begin), inner_it_(in_begin), outer_begin_(out_begin), outer_end_(out_end),
              inner_begin_(in_begin), inner_end_(in_end) {
            if (outer_begin_ == outer_end_) {
                return;
            }
            if (inner_begin_ == inner_end_) {
                return;
            }

            if (outer_it_ == outer_end_) {
                return;
            }
            if (inner_it_ == outer_end_) {
                return;
            }
        }

        reference operator*() const {
            return make_pair(*outer_it_, *inner_it_);
        }
        pointer operator->() const {
            return make_pair(&*outer_it_, &*inner_it_);
        }

        product_iterator& operator++() {
            ++inner_it_;
            if (inner_it_ == inner_end_) {
                inner_it_ = inner_begin_;
                ++outer_it_;
            }

            return *this;
        }

        product_iterator operator++(int) {
            product_iterator it(*this);
            ++*this;
            return it;
        }

        product_iterator& operator--() {
            if (outer_it_ == outer_end_) {
                --outer_it_;
                inner_it_ = inner_end_;
            }
            if (inner_it_ == inner_begin_) {
                if (outer_it_ == outer_begin_) {
                    return *this;
                }
                inner_it_ = inner_end_;
                --outer_it_;
            }
            --inner_it_;

            return *this;
        }

        product_iterator operator--(int) {
            product_iterator it(*this);
            --*this;
            return it;
        }

        friend bool operator==(const product_iterator& a, const product_iterator& b) {
            return (a.outer_it_ == b.outer_it_) && (a.inner_it_ == b.inner_it_);
        }

        friend bool operator!=(const product_iterator& a, const product_iterator& b) {
            return !(a == b);
        }

    private:

        outer_iterator outer_it_;
        outer_iterator outer_begin_;
        outer_iterator outer_end_;
        inner_iterator inner_it_;
        inner_iterator inner_begin_;
        inner_iterator inner_end_;
    };

} // namespace Hammer

#endif
