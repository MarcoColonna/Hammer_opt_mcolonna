///
/// @file  BlockIndexing.hh
/// @brief Outer product tensor indexer
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_MATH_MULTIDIM_BLOCKINDEXING
#define HAMMER_MATH_MULTIDIM_BLOCKINDEXING

#include <vector>

#include "Hammer/Math/MultiDimensional.fhh"
#include "Hammer/Math/MultiDim/AlignedIndexing.hh"
#include "Hammer/Math/MultiDim/LabeledIndexing.hh"
#include "Hammer/Math/MultiDim/ISingleContainer.hh"

namespace Hammer::MultiDimensional {

    class OuterElemIterator {
    public:

        using EntryType = std::vector<std::pair<SharedTensorData, bool>>;

        explicit OuterElemIterator(const EntryType& entry);

        OuterElemIterator& operator=(const OuterElemIterator&) = delete;
        OuterElemIterator& operator=(OuterElemIterator&&) = delete;
        OuterElemIterator(OuterElemIterator&& other) = default;

    private:

        OuterElemIterator(const OuterElemIterator& other);

    public:

        [[nodiscard]] OuterElemIterator begin() const;

        [[nodiscard]] OuterElemIterator end() const;

        OuterElemIterator& operator++();

        [[nodiscard]] OuterElemIterator operator++(int n);

        [[nodiscard]] IContainer::ElementType operator*();
        [[nodiscard]] PositionType position(IndexType idx) const;
        [[nodiscard]] bool isAligned(IndexType idx) const;

        [[nodiscard]] bool isSame(const OuterElemIterator& other) const;


    private:

        void incrementEntry(size_t position, int n);
        void setInitialState();
        void setFinalState();

        using IteratorList = std::vector<ISingleContainer::NonZeroIt>;
        using ContainerList = std::vector<ISingleContainer*>;

        const EntryType& _entry;
        ContainerList _containers;
        IteratorList _it;
        IndexList _dimensions;
        IndexList _currentIdx;
    };

    inline bool operator==(const OuterElemIterator& lhs, const OuterElemIterator& rhs) {
        return lhs.isSame(rhs);
    }

    inline bool operator!=(const OuterElemIterator& lhs, const OuterElemIterator& rhs) {
        return !(lhs == rhs);
    }

    using DotGroupType = std::tuple<IndexList, IndexList, IndexPairList>;
    using DotGroupList = std::vector<DotGroupType>;

    class BlockIndexing {
    public:

        BlockIndexing() = default;
        explicit BlockIndexing(const std::vector<IndexList>& dims, const std::vector<LabelsList>& labels);
        explicit BlockIndexing(const LabeledIndexing<AlignedIndexing>& left,
                               const LabeledIndexing<AlignedIndexing>& right);

        BlockIndexing(const BlockIndexing&) = default;
        BlockIndexing(BlockIndexing&&) = default;
        BlockIndexing& operator=(const BlockIndexing&) = default;
        BlockIndexing& operator=(BlockIndexing&&) = default;

        ~BlockIndexing() = default;

        /// @brief  rank of the tensor
        /// @return the number of components
        [[nodiscard]] size_t rank() const;

        /// @brief  dimension of a specific component
        /// @param[in] index   the component index
        /// @return  the dimension
        [[nodiscard]] IndexType dim(IndexType index) const;

        /// @brief  dimension of a specific component by label.
        ///         If multiple components with the same label exists,
        ///         returns the dimension of first (they should all be the same
        ///         by consistency)
        /// @param[in] label   the component label
        /// @return  the dimension
        [[nodiscard]] IndexType dim(IndexLabel label) const;

        /// @brief   get the labels of all the indices at once
        /// @return  the list of labels
        [[nodiscard]] LabelsList labels() const;

        /// @brief  returns only the labels corresponding to spin indices
        /// @return   the set of labels of the spin indices (squashes repetitions if present)
        [[nodiscard]] UniqueLabelsList spinIndices() const;

        /// @brief  get the dimensions of all the indices at once
        /// @return  the list of dimensions
        [[nodiscard]] IndexList dims() const;

        /// @brief  the number of elements (product of all the dimensions)
        [[nodiscard]] PositionType numValues() const;

        [[nodiscard]] std::vector<IndexList> splitIndices(const IndexList& indices) const;

        [[nodiscard]] std::vector<IndexList::const_iterator> splitIndices(IndexList::const_iterator first,
                                                                          IndexList::const_iterator last) const;

        /// @brief  check that the indices are within range for each component
        /// @param[in] indices the element indices
        /// @return true if all indices are OK
        [[nodiscard]] bool checkValidIndices(const IndexList& indices) const;

        [[nodiscard]] bool checkValidIndices(IndexList::const_iterator first, IndexList::const_iterator last) const;

        [[nodiscard]] bool checkValidIndices(const std::vector<IndexList>& splits) const;

        /// @brief  returns the position of the indices in the two tensor (this and another) that can be
        ///         contracted together, given a set of allowed index labels
        /// @param[in] otherLabels the list of labels of the other tensor
        /// @param[in] indices  the list of labels of the allowed indices to be contracted
        /// @param[in] sortedBySecond whether the result should be sorted according to the second element of the
        /// pair
        /// @return pairs of indices corresponding to coordinates in this and other tensor that needs to be
        /// contracted together
        [[nodiscard]] IndexPairList getSameLabelPairs(const LabelsList& otherLabels, const UniqueLabelsList& indices,
                                                      bool sortedBySecond = true) const;

        /// @brief  returns the position of the indices that can be traced together, given a set of allowed index
        ///         labels
        /// @param[in] indices  the list of labels of the allowed indices to be traced
        /// @return pairs of coordinate indices that needs to be traced together
        [[nodiscard]] IndexPairList getOppositeLabelPairs(const UniqueLabelsList& indices) const;

        [[nodiscard]] IndexType labelIndex(IndexLabel label) const;

        [[nodiscard]] bool isSameLabelShape(const LabelsList& otherLabels, const IndexList& otherIndices) const;

        template <typename S>
        [[nodiscard]] bool isSameLabelShape(const LabeledIndexing<S>& other) const;

        [[nodiscard]] bool isSameLabelShape(const BlockIndexing& other, bool includeBlockShapes = false) const;

        void flipLabels();

        /// @brief returns the index of the block and the index within the block associated to a global index
        /// @param[in] position the global index
        /// @return the pair (block index, local index)
        [[nodiscard]] IndexPair getElementIndex(IndexType position) const;

        [[nodiscard]] const LabeledIndexing<AlignedIndexing>& getSubIndexing(IndexType position) const;

        /// @brief number of blocks
        /// @return the size of the block array
        [[nodiscard]] size_t numSubIndexing() const;

        [[nodiscard]] size_t maxSubRank() const;

        [[nodiscard]] std::vector<std::tuple<IndexList, FlagList, PositionType>>
        processShifts(const DotGroupList& chunks, IndexPairMember which) const;

        [[nodiscard]] PositionType splitPosition(const OuterElemIterator& currentPosition, const DotGroupType& chunk,
                                                 const IndexList& outerShiftsInnerPositions,
                                                 const FlagList& isOuter, IndexList& innerList,
                                                 FlagList& innerAdded, bool shouldCompare = false) const;

        [[nodiscard]] PositionType buildFullPosition(const OuterElemIterator& current,
                                                     const IndexList& chunkIndices) const;

        [[nodiscard]] LabelSignature getBlockLabelSignature(IndexType position) const;
        [[nodiscard]] LabelSignature getLabelSignature() const;

    private:

        void calc();

        std::vector<LabeledIndexing<AlignedIndexing>> _subIndexing;
        IndexList _splitIndices;
        PositionList _splitPads;
        LabeledIndexing<AlignedIndexing> _globalIndexing;
    };

    template <typename S>
    [[nodiscard]] bool BlockIndexing::isSameLabelShape(const LabeledIndexing<S>& other) const {
        return isSameLabelShape(other.labels(), other.dims());
    }

} // namespace Hammer::MultiDimensional

#endif
