///
/// @file  BinnedIndexing.hh
/// @brief Binned tensor (histogram) indexer
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_MATH_MULTIDIM_BinnedIndexing
#define HAMMER_MATH_MULTIDIM_BinnedIndexing

#include "Hammer/Math/MultiDimensional.fhh"

namespace Hammer::MultiDimensional {


    using BinValue = std::vector<double>;
    using BinEdgeList = std::vector<std::vector<double>>;
    using BinRange = std::pair<double, double>;
    using BinRangeList = std::vector<BinRange>;

    template <class BasicIndexing>
    class BinnedIndexing : public BasicIndexing {

    public:

        BinnedIndexing();

        BinnedIndexing(const IndexList& dimensions, const BinEdgeList& edges, bool hasUnderOverFlow);
        BinnedIndexing(const IndexList& dimensions, const BinRangeList& ranges, bool hasUnderOverFlow);
        BinnedIndexing(const BinEdgeList& edges, bool hasUnderOverFlow);

        BinnedIndexing(const BinnedIndexing&) = default;
        BinnedIndexing(BinnedIndexing&&) = default;
        BinnedIndexing& operator=(const BinnedIndexing&) = default;
        BinnedIndexing& operator=(BinnedIndexing&&) = default;

        ~BinnedIndexing() noexcept = default;

        using BasicIndexing::dim;

        /// @brief   get the labels of all the indices at once
        /// @return  the list of labels
        [[nodiscard]] const BinEdgeList& edges() const;

        [[nodiscard]] const BinValue& edge(IndexType pos) const;

        [[nodiscard]] IndexList valueToPos(const BinValue& point) const;
        [[nodiscard]] bool isValid(const BinValue& point) const;

        [[nodiscard]] BinRangeList binEdges(IndexList position) const;
        [[nodiscard]] BinRange binEdge(IndexType position, IndexType coord) const;

        [[nodiscard]] bool isSameBinShape(const BinEdgeList& otherEdges, const IndexList& otherIndices,
                                          bool otherUnderOverFlow) const;
        [[nodiscard]] bool isSameBinShape(const BinRangeList& otherRanges, const IndexList& otherIndices,
                                          bool otherUnderOverFlow) const;

        template <typename S>
        [[nodiscard]] bool isSameBinShape(const BinnedIndexing<S>& other) const;

        [[nodiscard]] bool hasUnderOverFlow() const;

        [[nodiscard]] typename BasicIndexing::StrideMap getBinStrides(const UniqueIndexList& positions) const;

    private:

        [[nodiscard]] bool isValid() const;
        void fillEdges(const BinRangeList& ranges);

        static IndexList processEdges(const BinEdgeList& edges, bool hasUnderOverFlow);
        static IndexList processRanges(const IndexList& dimensions, bool hasUnderOverFlow);

        BinEdgeList _edges; ///< the labels of each tensor index
        bool _hasUnderOverFlow;
    };

} // namespace Hammer::MultiDimensional

#include "Hammer/Math/MultiDim/BinnedIndexingDefs.hh"

#endif
