///
/// @file  LabeledIndexing.hh
/// @brief Tensor labeled indexer
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_MATH_MULTIDIM_LABELEDINDEXING
#define HAMMER_MATH_MULTIDIM_LABELEDINDEXING

#include "Hammer/Math/MultiDimensional.fhh"
#include "Hammer/Tools/Utils.hh"

namespace Hammer::MultiDimensional {


    template <class BasicIndexing>
    class LabeledIndexing : public BasicIndexing {

    public:

        LabeledIndexing();

        LabeledIndexing(const IndexList& dimensions, const LabelsList& labels);

        LabeledIndexing(const LabeledIndexing&) = default;
        LabeledIndexing(LabeledIndexing&&) = default;
        LabeledIndexing& operator=(const LabeledIndexing&) = default;
        LabeledIndexing& operator=(LabeledIndexing&&) = default;

        ~LabeledIndexing() = default;

        /// @brief  dimension of a specific component by label.
        ///         If multiple components with the same label exists,
        ///         returns the dimension of first (they should all be the same
        ///         by consistency)
        /// @param[in] label   the component label
        /// @return  the dimension
        [[nodiscard]] IndexType dim(IndexLabel label) const;

        using BasicIndexing::dim;

        /// @brief   get the labels of all the indices at once
        /// @return  the list of labels
        [[nodiscard]] const LabelsList& labels() const;

        /// @brief  returns only the labels corresponding to spin indices
        /// @return   the set of labels of the spin indices (squashes repetitions if present)
        [[nodiscard]] UniqueLabelsList spinIndices() const;

    protected:

        /// @brief  returns the set of indices of this tensor that can be contracted
        ///         using dot with those of another given tensor
        /// @param[in] otherLabels the other tensor
        /// @return the set of labels
        [[nodiscard]] UniqueLabelsList sameIndices(const LabelsList& otherLabels) const;

        /// @brief returns the pairs of indices that can be traced over. Valid pairs
        ///        are composed by opposite sign labels which correspond to (reference)
        ///        spin indices
        /// @return the set of label pairs
        [[nodiscard]] LabelPairsSet oppositeIndices() const;

    public:

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

        [[nodiscard]] MaybeBool isSameLabelShapeNeedHc(const LabelsList& otherLabels,
                                                       const IndexList& otherIndices) const;

        template <typename S>
        [[nodiscard]] MaybeBool isSameLabelShapeNeedHc(const LabeledIndexing<S>& other) const;

        [[nodiscard]] bool canAddAt(const LabelsList& otherLabels, const IndexList& otherIndices, IndexLabel coord,
                                    IndexType position) const;

        void flipLabels();

        [[nodiscard]] LabelSignature getLabelSignature() const;

    private:

        LabelsList _labels; ///< the labels of each tensor index
    };

    LabelsList flipListOfLabels(LabelsList labels);
} // namespace Hammer::MultiDimensional

#include "Hammer/Math/MultiDim/LabeledIndexingDefs.hh"

#endif
