///
/// @file  TestAmplD2starDstarDGamPi.cc
/// @brief Tests for AmplD2starDstarDGamPi
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "Hammer/Amplitudes/AmplD2starDstarDGamPi.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/FourMomentum.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Math/MultiDim/VectorContainer.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"
#include "gtest/gtest.h"
// May need other stuff here

using namespace std;
using namespace complex_literals;

namespace Hammer {

    namespace MD = MultiDimensional;

    TEST(AmplD2starDstarDGamPiTest, evalTauDec) {

        FourMomentum pDssmes{2461., 0, 0, 0};
        FourMomentum pDsmes{2047.62210483543, 243.33950461271, 280.828758699094, -120.736812915758};
        FourMomentum pPmes{413.377895164567, -243.33950461271, -280.828758699094, 120.736812915758};
        FourMomentum pDmes{1930.33702358413, 269.34371344472, 376.7215687779, -183.05968498};
        FourMomentum kGmes{117.2850812513, -26.004208832009, -95.892810078805, 62.32287206424};
        FourMomentum kNuTau{103.70179046789, 31.76995444103, 97.77786575997, -13.57277821547};
        FourMomentum pTau{4126.6677747495, -31.76995444103, -97.77786575997, 4124.0504461171};

        // Evaluate at PS point
        AmplD2starDstarDGamPi ampl;
        ampl.eval(Particle(pDssmes, 425),
                  {Particle(pDsmes, 413), Particle(pPmes, -211), Particle(pDmes, 411), Particle(kGmes, 22)},
                  {Particle(kNuTau, 16), Particle(pTau, -15)});

        // Compare to direct evaluation
        Tensor amplEval{
            "amplEval",
            MD::makeVector({1, 5, 2}, {FF_DSSD2STAR, SPIN_DSSD2STAR, SPIN_GAMMA},
                           {-119191.64928254 - 563351.38120782 * 1i, -1.05238874596507e6 - 199197.11816746 * 1i,
                            473714.53185869 + 824303.1371519 * 1i, -396495.1911441 + 130961.48701221 * 1i,
                            476711.74731808 + 378114.09544055 * 1i, -476711.74731808 + 378114.09544055 * 1i,
                            396495.1911441 + 130961.48701221 * 1i, -473714.53185869 + 824303.1371519 * 1i,
                            1.05238874596507e6 - 199197.11816746 * 1i, 119191.64928254 - 563351.38120782 * 1i})};

        const Tensor& t = ampl.getTensor();
        // Check TEs.
        for (IndexType idx1 = 0; idx1 < 1; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 5; ++idx2) {
                for (IndexType idx3 = 0; idx3 < 2; ++idx3) {
                    double comRe =
                        compareVals(t.element({idx1, idx2, idx3}).real(), amplEval.element({idx1, idx2, idx3}).real());
                    double comIm =
                        compareVals(t.element({idx1, idx2, idx3}).imag(), amplEval.element({idx1, idx2, idx3}).imag());
                    EXPECT_NEAR(comRe, 1., 1e-4);
                    EXPECT_NEAR(comIm, 1., 1e-4);
                }
            }
        }
    }

    TEST(AmplD2starDstarDGamPiTest, evalTauDecNoRef) {

        // Momentum decs

        FourMomentum pDssmes{2461., 0, 0, 0};
        FourMomentum pDsmes{2047.62210483543, 243.33950461271, 280.828758699094, -120.736812915758};
        FourMomentum pPmes{413.377895164567, -243.33950461271, -280.828758699094, 120.736812915758};
        FourMomentum pDmes{1930.33702358413, 269.34371344472, 376.7215687779, -183.05968498};
        FourMomentum kGmes{117.2850812513, -26.004208832009, -95.892810078805, 62.32287206424};
        //        FourMomentum kNuTau{1., 0.707106781186548, 0, 0.707106781186548};
        //        FourMomentum pTau{1., -0.70710678118655, 0, 0.707106781186548};


        // Evaluate at PS point
        AmplD2starDstarDGamPi ampl;
        ampl.eval(Particle(pDssmes, 425),
                  {Particle(pDsmes, 413), Particle(pPmes, -211), Particle(pDmes, 411), Particle(kGmes, 22)}, {});

        // Compare to direct evaluation
        Tensor amplEval{
            "amplEval",
            MD::makeVector({1, 5, 2}, {FF_DSSD2STAR, SPIN_DSSD2STAR, SPIN_GAMMA},
                           {-234701.5638754 + 525819.93484642 * 1i, 734315.25181669 + 779732.43838483 * 1i,
                            930344.71081755 - 195805.61449139 * 1i, 2028.0233552 + 417558.66031267 * 1i,
                            476711.74731808 + 378114.09544055 * 1i, -476711.74731808 + 378114.09544055 * 1i,
                            -2028.0233552 + 417558.66031267 * 1i, -930344.71081755 - 195805.61449139 * 1i,
                            -734315.25181669 + 779732.43838483 * 1i, 234701.5638754 + 525819.93484642 * 1i})};

        const Tensor& t = ampl.getTensor();
        // Check TEs.
        for (IndexType idx1 = 0; idx1 < 1; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 5; ++idx2) {
                for (IndexType idx3 = 0; idx3 < 2; ++idx3) {
                    double comRe =
                        compareVals(t.element({idx1, idx2, idx3}).real(), amplEval.element({idx1, idx2, idx3}).real());
                    double comIm =
                        compareVals(t.element({idx1, idx2, idx3}).imag(), amplEval.element({idx1, idx2, idx3}).imag());
                    EXPECT_NEAR(comRe, 1., 1e-4);
                    EXPECT_NEAR(comIm, 1., 1e-4);
                }
            }
        }
    }

} // namespace Hammer
