///
/// @file  AlignedIndexing.hh
/// @brief Sparse tensor indexer
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_MATH_MULTIDIM_ALIGNEDINDEXING
#define HAMMER_MATH_MULTIDIM_ALIGNEDINDEXING

#ifndef ALIGNEDINDEXING_FRIENDS
#define ALIGNEDINDEXING_FRIENDS typedef bool DummyFriend
#endif

#include <cstdint>
#include <map>

#include "Hammer/Math/MultiDimensional.fhh"

namespace Hammer {


    enum class IndexPairMember : std::uint8_t { Left, Right, Both };

    namespace MultiDimensional {

        class AlignedIndexing {

            ALIGNEDINDEXING_FRIENDS;

        public:

            AlignedIndexing() = default;

            explicit AlignedIndexing(IndexList dimensions);

            AlignedIndexing(const AlignedIndexing&) = default;
            AlignedIndexing(AlignedIndexing&&) = default;
            AlignedIndexing& operator=(const AlignedIndexing&) = default;
            AlignedIndexing& operator=(AlignedIndexing&&) = default;

            ~AlignedIndexing() = default;

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

            /// @brief  the number of elements (product of all the dimensions)
            [[nodiscard]] PositionType numValues() const;

            /// @brief convert the indices into the (aligned) position indicizing
            ///        a sparse tensor container organized as row-major
            /// @param[in] indices the element indices
            /// @return the absolute (aligned) position
            [[nodiscard]] PositionType indicesToPos(const IndexList& indices) const;

            [[nodiscard]] PositionType indicesToPos(IndexList::const_iterator first,
                                                    IndexList::const_iterator last) const;

            /// @brief convert the absolute aligned position to the absolute unaligned position
            ///        that can be used with a `SequentialIndexing` object
            /// @param[in] alignedPosition the absolute aligned position
            /// @return the absolute unaligned position
            [[nodiscard]] PositionType alignedPosToPos(PositionType alignedPosition) const;

            /// @brief convert the absolute unaligned position (e.g. coming from a `SequentialIndexing` object)
            ///        to the absolute aligned position
            /// @param[in] position the absolute unaligned position
            /// @return the absolute aligned position
            [[nodiscard]] PositionType posToAlignedPos(PositionType position) const;

            /// @brief  check that the indices are within range for each component
            /// @param[in] indices the element indices
            /// @return true if all indices are OK
            [[nodiscard]] bool checkValidIndices(const IndexList& indices) const;

            [[nodiscard]] bool checkValidIndices(IndexList::const_iterator first, IndexList::const_iterator last) const;

            /// @brief convert the absolute aligned position (in row-major convention)
            ///        into the list of indices
            /// @param[in] alignedPosition  the absolute aligned position
            /// @param[out] result  the indices
            void posToIndices(PositionType alignedPosition, IndexList& result) const;

            /// @brief extract the value of the i-th index from an absolute aligned position
            /// @param[in] alignedPosition  the absolute position
            /// @param[in] indexPosition  the index coordinate i
            /// @return the index value
            [[nodiscard]] IndexType ithIndexInPos(PositionType alignedPosition, IndexType indexPosition) const;

            /// @brief split an absolute aligned position by separating contracted and uncontracted indices
            ///        and returning the absolute position made of uncontracted indices and the coordinates
            ///        of the contracted ones properly reordered. This function is the workhorse of Dot and Trace.
            ///        The method accepts an optional coordinate list for the contracted indices: if that does not match
            ///        the coordinate list extracted from the position the function return 0 for the external position
            ///        and an empty list for the contracted indices.
            /// @param[in] alignedPosition the absolute aligned index
            /// @param[in] outerShiftsInnerPositions a list (of rank length) with a positional shift for each
            /// uncontracted index and the corresponding position for the contracted index
            /// @param[in] isOuter a list (of rank length) marking whether an index is not contracted
            /// @param[in] innerList
            /// @param[in] innerAdded if present, the list of values that the contracted coordinates should have
            /// @param[in] shouldCompare if present, the list of values that the contracted coordinates should have
            /// @return a pair (absolute aligned position made of the uncontracted indices, reduced index list for the
            /// contracted indices reordered according to the shifts)
            [[nodiscard]] PositionType splitPosition(PositionType alignedPosition,
                                                     const IndexList& outerShiftsInnerPositions,
                                                     const std::vector<bool>& isOuter, IndexList& innerList,
                                                     std::vector<bool>& innerAdded, bool shouldCompare = false) const;

            /// @brief build the input lists necessary for calling splitPosition based on which indices are being
            /// contracted
            ///        and whether the tensor indexed by this class is the on left or right of Dot or it is the argument
            ///        of Trace
            /// @param[in] pairs the list of contractions (pairs od index position being contracted
            /// @param[in] which whether the tensor is the first or second argument of Dot or it is the argument of
            /// Trace
            /// @return the lists to be fed to splitPosition  (outerShiftsInnerPositions and isOuter) and the maximum
            /// aligned index of the uncontracted indices
            [[nodiscard]] std::tuple<IndexList, std::vector<bool>, PositionType>
            processShifts(const IndexPairList& pairs, IndexPairMember which) const;

            /// @brief get the maximum allowed value for the absolute position
            /// @param[in] aligned whether the maximum absolute position is aligned or not
            /// @return the maximum allowed absolute position
            [[nodiscard]] PositionType maxIndex(bool aligned = true) const;

            /// @brief extends an unaligned absolute position from a rank N-1 tensor to the corresponding
            ///        aligned position for this rank N tensor given the missing index and its value.
            ///        Used in AddAt
            /// @param position the sub-tensor unaligned position
            /// @param stride the stride (in the sub-tensor) of the index preceding the missing index
            /// @param indexPosition the missing index coordinate
            /// @param indexValue  the missing index value
            /// @return the corresponding absolute aligned position
            [[nodiscard]] PositionType extendPosition(PositionType position, PositionType stride,
                                                      IndexType indexPosition, IndexType indexValue) const;

            /// @brief extends an aligned absolute position from a rank N-1 tensor to the corresponding
            ///        aligned position for this rank N tensor given the missing index and its value.
            ///        Used in AddAt
            /// @param alignedPosition the sub-tensor aligned position
            /// @param indexPosition the missing index coordinate
            /// @param indexValue  the missing index value
            /// @return the corresponding absolute aligned position
            [[nodiscard]] PositionType extendAlignedPosition(PositionType alignedPosition, IndexType indexPosition,
                                                             IndexType indexValue) const;

            // /// @brief
            // /// @param[in] index
            // /// @param[in] positions
            // /// @param[in] values
            // /// @param[in] result
            // /// @return
            // bool isIndexMatchingSecond(size_t index, const IndexPairList& positions, const IndexList& values,
            // IndexList& result) const;

            template <typename BasicIndexing>
            [[nodiscard]] bool isSameShape(const BasicIndexing& other) const;

            [[nodiscard]] bool isSameShape(const IndexList& indices) const;

        private:

            void calc();
            static PositionType calcPadding(const IndexList& dimensions, IndexList& pads);
            static void calcMasks(const IndexList& dimensions, IndexList& masks);
            static void calcUnaligned(const IndexList& dimensions, const IndexList& pads, PosIndexPairList& unaligned);

            IndexList _dimensions; ///< the dimensions of each tensor index
            IndexList _alignPads;  ///< the strides for each tensor index (necessary to convert coordinates to position
                                   ///< in `_data`)
            IndexList _alignMasks; ///< the strides for each tensor index (necessary to convert coordinates to position
                                   ///< in `_data`)
            PositionType _maxIndex;
            PositionType _maxAlignedIndex;
            PosIndexPairList _unalignedEntries;
        };

        template <typename BasicIndexing>
        [[nodiscard]] bool AlignedIndexing::isSameShape(const BasicIndexing& other) const {
            return isSameShape(other.dims());
        }

    } // namespace MultiDimensional

} // namespace Hammer

#endif
