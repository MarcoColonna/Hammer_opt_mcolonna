///
/// @file  Dot.hh
/// @brief Tensor dot product algorithm
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_MATH_MULTIDIM_OPS_DOT
#define HAMMER_MATH_MULTIDIM_OPS_DOT

#include "Hammer/Math/MultiDimensional.fhh"
#include "Hammer/Math/MultiDim/BlockIndexing.hh"
#include "Hammer/Math/MultiDim/SequentialIndexing.hh"
#include "Hammer/Math/MultiDim/LabeledIndexing.hh"

namespace Hammer::MultiDimensional {

    class IContainer;
    class VectorContainer;
    class SparseContainer;
    class OuterContainer;

    namespace Ops {

        class Dot final {
        public:

            explicit Dot(IndexPairList indices, std::pair<bool, bool> shouldHC = {false, false}) noexcept;

            IContainer* operator()(VectorContainer& a, const VectorContainer& b);
            IContainer* operator()(SparseContainer& a, const SparseContainer& b);

            IContainer* operator()(SparseContainer& a, const VectorContainer& b);

            IContainer* operator()(OuterContainer& a, const OuterContainer& b);
            IContainer* operator()(OuterContainer& a, const SparseContainer& b);
            IContainer* operator()(SparseContainer& a, const OuterContainer& b);
            IContainer* operator()(VectorContainer& a, const OuterContainer& b);
            IContainer* operator()(OuterContainer& a, const IContainer& b);

            IContainer* operator()(IContainer& a, const IContainer& b);

            static IContainer* error(IContainer& /*unused*/, const IContainer& /*unused*/);

        private:

            [[nodiscard]] std::pair<IndexList, LabelsList> getNewIndexLabels(const IContainer& first,
                                                                             const IContainer& second) const;
            [[nodiscard]] static std::pair<IndexList, LabelsList>
            getNewIndexLabels(const LabeledIndexing<AlignedIndexing>& lhs, const BlockIndexing& rhs,
                              const DotGroupType& chunk);
            [[nodiscard]] static std::pair<IndexList, LabelsList>
            getNewIndexLabels(const LabeledIndexing<SequentialIndexing>& lhs, const BlockIndexing& rhs,
                              const DotGroupType& chunk);
            [[nodiscard]] static std::pair<IndexList, LabelsList>
            getNewIndexLabels(const BlockIndexing& lhs, const LabeledIndexing<AlignedIndexing>& rhs,
                              const DotGroupType& chunk);
            [[nodiscard]] static std::pair<IndexList, LabelsList>
            getNewIndexLabels(const std::pair<SharedTensorData, bool>& a, const IndexPairList& contracted);


            [[nodiscard]] IndexList combineIndex(const IndexList& a, const IndexList& b) const;

            static SharedTensorData calcSharedDot(SharedTensorData origin, const IContainer& other,
                                                  const IndexPairList& indices,
                                                  std::pair<bool, bool> shouldHC = {false, false});

        public:

            [[nodiscard]] static IContainer*
            contractStar(const std::pair<SharedTensorData, bool>& tensorA,
                         const std::vector<std::pair<SharedTensorData, bool>>& tensorsB,
                         const std::vector<IndexPair>& contractions);


            [[nodiscard]] static IContainer* contractSingles(const std::pair<SharedTensorData, bool>& tensorA,
                                                             const std::pair<SharedTensorData, bool>& tensorB);

            [[nodiscard]] static IContainer* contractR2Boomerang(
                const std::pair<SharedTensorData, bool>& tensorA,
                const std::pair<std::pair<SharedTensorData, bool>, std::pair<SharedTensorData, bool>>& tensorsB,
                const std::pair<IndexPair, IndexPair>& contractions);

            // DotGroupList partitionContractions(const BlockIndexing& lhs, const BlockIndexing& rhs) const;
            [[nodiscard]] DotGroupList partitionContractions(const LabeledIndexing<AlignedIndexing>& lhs,
                                                             const BlockIndexing& rhs) const;
            [[nodiscard]] DotGroupList partitionContractions(const BlockIndexing& lhs,
                                                             const LabeledIndexing<AlignedIndexing>& rhs) const;

            IndexPairList _indices;
            std::pair<bool, bool> _hc;
            UniqueIndexList _idxLeft;
            UniqueIndexList _idxRight;
        };

    } // namespace Ops

} // namespace Hammer::MultiDimensional


#endif
