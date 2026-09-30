///
/// @file  TestHistogram.cc
/// @brief Tests for Histogram
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "Hammer/Math/Histogram.hh"
#include "Hammer/Math/Tensor.hh"
#include "Hammer/Math/MultiDim/VectorContainer.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

#include "gtest/gtest.h"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    TEST(HistogramTest, Create) {

        Tensor defaultTensor{"Value", MD::makeEmptySparse({3, 3}, {WILSON_BCENU, WILSON_BCENU_HC})};
        IndexList dims = {10, 20, 5};
        MD::BinRangeList ranges = {{0, 10}, {0, 20}, {0, 5}};
        HistogramDefinition def{false, dims, ranges, false, false};
        unique_ptr<Histogram> histo = makeHistogram("TestHisto", def, defaultTensor);

        EXPECT_EQ(histo->rank(), 3);
        EXPECT_EQ(histo->name(), "TestHisto");
        auto dimsTest = histo->dims();
        for (IndexType idx = 0; idx < 3; ++idx) {
            EXPECT_EQ(dimsTest[idx], dims[idx]);
        }
        for (IndexType idx = 0; idx < 3; ++idx) {
            EXPECT_EQ(histo->dim(idx), dims[idx]);
        }
    }

    TEST(HistogramTest, FillExtract) {

        Tensor defaultTensor{"Value", MD::makeEmptySparse({3, 3}, {WILSON_BCENU, WILSON_BCENU_HC})};
        IndexList dims = {10, 20, 5};
        MD::BinRangeList ranges = {{0, 10}, {0, 20}, {0, 5}};
        HistogramDefinition def{false, dims, ranges, false, false};
        unique_ptr<Histogram> histo = makeHistogram("TestHisto", def, defaultTensor);

        // fill Bins
        vector<Tensor> inTensor;
        //        vector<Tensor> inSqTensor;
        for (IndexType idx = 0; idx < 2; ++idx) {
            Tensor tmp{"Value", MD::makeVector({3, 3}, {WILSON_BCENU, WILSON_BCENU_HC},
                                               {1., 2., 1. + 1.i, 2., 4., 2. + 2.i, 1. - 1.i, 2. - 2.i, 2.})};
            tmp *= static_cast<double>(idx + 1);
            inTensor.push_back(tmp);
            auto idx2 = static_cast<IndexType>(3 + idx);
            histo->fillBin(tmp, {5, idx, idx2});
        }

        // Add one more set of events
        for (IndexType idx = 0; idx < 2; ++idx) {
            Tensor tmp{"Value", MD::makeVector({3, 3}, {WILSON_BCENU, WILSON_BCENU_HC},
                                               {2., -4., 2. - 2.i, -4., 8., -4. + 4.i, 2. + 2.i, -4. - 4.i, 4.})};
            tmp *= static_cast<double>(idx + 1);
            inTensor[idx] += tmp;
            auto idx2 = static_cast<IndexType>(3 + idx);
            histo->fillBin(tmp, {5, idx, idx2});
        }

        // read Bins
        for (IndexType idx = 0; idx < 2; ++idx) {
            Tensor readTensor = histo->getWeight(5, idx, 3 + idx);
            for (IndexType idx1 = 0; idx1 < 2; ++idx1) {
                for (IndexType idx2 = 0; idx2 < 2; ++idx2) {
                    EXPECT_EQ(readTensor.element({idx1, idx2}).real(), inTensor[idx].element({idx1, idx2}).real());
                    EXPECT_EQ(readTensor.element({idx1, idx2}).imag(), inTensor[idx].element({idx1, idx2}).imag());
                }
                EXPECT_EQ(histo->getNumEvents(5, idx, 3 + idx), 2);
            }
            //            Tensor readSqTensor = histo->getSquaredWeight(5,idx,3+idx);
            //            for(size_t idx1 = 0; idx1 < 2; ++idx1){
            //                for(size_t idx2 = 0; idx2 < 1; ++idx2){
            //                    for(size_t idx3 = 0; idx3 < 2; ++idx3){
            //                        for(size_t idx4 = 0; idx4 < 1; ++idx4){
            //                            EXPECT_EQ(readSqTensor.element({idx1,idx2,idx3,idx4}).real(),inSqTensor[idx].element({idx1,idx2,idx3,idx4}).real());
            //                            EXPECT_EQ(readSqTensor.element({idx1,idx2,idx3,idx4}).imag(),inSqTensor[idx].element({idx1,idx2,idx3,idx4}).imag());
            //                        }
            //                    }
            //                }
            //            }
        }

        Tensor wc{"wc", MD::makeVector({3}, {WILSON_BCENU}, {2., 1.i, 1.})};
        // Contraction wc.inTensor[0].wcconjugate = 30
        auto vecWgts = histo->evalToList({outerSquare(wc)});
        for (size_t idx = 0; idx < 2; ++idx) {
            EXPECT_EQ(vecWgts[5l * 20 * 5 + idx * 5 + 3 + idx].sumWi, (idx + 1) * 30);
            EXPECT_EQ(vecWgts[5l * 20 * 5 + idx * 5 + 3 + idx].n, 2);
        }
    }

} // namespace Hammer
