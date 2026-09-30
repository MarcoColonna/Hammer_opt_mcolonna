///
/// @file  SequentialIndexing.hh
/// @brief Non-sparse tensor indexer
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_MATH_MULTIDIM_SEQUENTIALINDEXING
#define HAMMER_MATH_MULTIDIM_SEQUENTIALINDEXING

#ifndef SEQUENTIALINDEXING_FRIENDS
#define SEQUENTIALINDEXING_FRIENDS typedef bool DummyFriend
#endif

#include <map>

#include "Hammer/Math/MultiDimensional.fhh"

namespace Hammer::MultiDimensional {

    class SequentialIndexing {

        SEQUENTIALINDEXING_FRIENDS;

    public:

        SequentialIndexing() = default;

        explicit SequentialIndexing(IndexList dimensions);

        SequentialIndexing(const SequentialIndexing& /*unused*/) = default;
        SequentialIndexing(SequentialIndexing&& /*unused*/) noexcept = default;
        SequentialIndexing& operator=(const SequentialIndexing& /*unused*/) = default;
        SequentialIndexing& operator=(SequentialIndexing&& /*unused*/) noexcept = default;

        ~SequentialIndexing() noexcept = default;


        using StrideMap = std::vector<std::tuple<PositionType, long, bool>>;

        /// @brief  rank of the tensor
        /// @return the number of components
        [[nodiscard]] size_t rank() const;

        /// @brief  dimension of a specific component
        /// @param[in] index   the component index
        /// @return  the dimension
        [[nodiscard]] IndexType dim(IndexType index) const;

        /// @brief  get the dimensions of all the indices at once
        /// @return  the list of dimensions
        [[nodiscard]] const IndexList& dims() const;

        /// @brief  stride of a specific component, i.e. the product of the dimension of
        ///         the indices to right of the given one (row-major convention).
        ///         The last index has stride equal to 1.
        /// @param[in] index   the component index
        /// @return  the stride
        [[nodiscard]] PositionType stride(IndexType index) const;

        /// @brief  list of strides, i.e. the product of the dimension of
        ///         the indices to right of the given one (row-major convention).
        ///         The last index has stride equal to 1.
        /// @return  the strides
        [[nodiscard]] const PositionList& strides() const;

        /// @brief  the number of elements (product of all the dimensions)
        [[nodiscard]] PositionType numValues() const;

        /// @brief convert the indices into the position indicizing
        ///        a sparse tensor container organized as row-major
        /// @param[in] indices the element indices
        /// @return the absolute position
        [[nodiscard]] PositionType indicesToPos(const IndexList& indices) const;

        [[nodiscard]] PositionType indicesToPos(IndexList::const_iterator first, IndexList::const_iterator last) const;

        /// @brief  check that the indices are within range for each component
        /// @param[in] indices the element indices
        /// @return true if all indices are OK
        [[nodiscard]] bool checkValidIndices(const IndexList& indices) const;

        [[nodiscard]] bool checkValidIndices(IndexList::const_iterator first, IndexList::const_iterator last) const;

        /// @brief convert the absolute position (in row-major convention)
        ///        into the list of indices
        /// @param[in] position  the absolute position
        /// @param[out] result  the indices
        void posToIndices(PositionType position, IndexList& result) const;

        /// @brief extract the value of the i-th index from an absolute position
        /// @param[in] position  the absolute position
        /// @param[in] indexPosition  the index coordinate i
        /// @return the index value
        [[nodiscard]] IndexType ithIndexInPos(PositionType position, IndexType indexPosition) const;

        [[nodiscard]] static PositionType build2ndPosition(PositionType reducedPosition, PositionType innerPosition,
                                                           const PositionPairList& conversion);

        [[nodiscard]] static PositionPair splitPosition(PositionType position, const StrideMap& conversion);

        /// @brief extends an absolute position from a rank N-1 tensor to the corresponding
        ///        position for this rank N tensor given the missing index and its value.
        ///        Used in AddAt
        /// @param position the sub-tensor position
        /// @param indexPosition the missing index coordinate
        /// @param indexValue  the missing index value
        /// @return the corresponding absolute position
        [[nodiscard]] PositionType extendPosition(PositionType position, IndexType indexPosition,
                                                  IndexType indexValue) const;

        /// @brief Get the Inner Outer Strides object
        ///
        /// @param positions
        /// @param secondStrides
        /// @param flipSecond
        /// @return
        [[nodiscard]] StrideMap getInnerOuterStrides(const IndexPairList& positions, const PositionList& secondStrides,
                                                     bool flipSecond = false) const;
        // StrideMap getInnerOuterStrides(const IndexBoolPairList& positions) const;

        /// @brief Get the Outer Strides2nd object
        ///
        /// @param positions
        /// @return
        [[nodiscard]] PositionPairList getOuterStrides2nd(const IndexPairList& positions) const;

        [[nodiscard]] PositionType reducedNumValues(const IndexPairList& indices) const;


        template <typename BasicIndexing>
        [[nodiscard]] bool isSameShape(const BasicIndexing& other) const;

        [[nodiscard]] bool isSameShape(const IndexList& indices) const;

    protected:

        [[nodiscard]] StrideMap buildStrideMap(const std::map<IndexType, long>& innerMap) const;

    private:

        void calcPadding();

        IndexList _dimensions; ///< the dimensions of each tensor index
        PositionList _strides; ///< the strides for each tensor index (necessary to convert coordinates to position
                               ///< in `_data`)
        PositionType _maxIndex;
    };


    template <typename BasicIndexing>
    [[nodiscard]] bool SequentialIndexing::isSameShape(const BasicIndexing& other) const {
        return isSameShape(other.dims());
    }

} // namespace Hammer::MultiDimensional

#endif
