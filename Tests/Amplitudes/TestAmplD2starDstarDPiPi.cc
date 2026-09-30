///
/// @file  TestAmplD2starDstarDPiPi.cc
/// @brief Tests for AmplD2starDstarDPiPi
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "Hammer/Amplitudes/AmplD2starDstarDPiPi.hh"
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

    TEST(AmplD2starDstarDPiPiTest, evalTauDec) {

        FourMomentum pDssmes{2461., 0, 0, 0};
        FourMomentum pDsmes{2047.62210483543, 243.33950461271, 280.828758699094, -120.736812915758};
        FourMomentum pPmes{413.377895164567, -243.33950461271, -280.828758699094, 120.736812915758};
        FourMomentum pDmes{1907.48119048676, 238.04888621381, 293.46680201094, -132.172069479921};
        FourMomentum pPDmes{140.140914348677, 5.290618398901, -12.63804331185, 11.4352565641629};
        FourMomentum kNuTau{103.70179046789, 31.76995444103, 97.77786575997, -13.57277821547};
        FourMomentum pTau{4126.6677747495, -31.76995444103, -97.77786575997, 4124.0504461171};


        // Evaluate at PS point
        AmplD2starDstarDPiPi ampl;
        ampl.eval(Particle(pDssmes, 425),
                  {Particle(pDsmes, 413), Particle(pPmes, -211), Particle(pDmes, 411), Particle(pPDmes, 111)},
                  {Particle(kNuTau, 16), Particle(pTau, -15)});

        // Compare to direct evaluation
        Tensor amplEval{"amplEval",
                        MD::makeVector({1, 5}, {FF_DSSD2STAR, SPIN_DSSD2STAR},
                                       {282791.34347149 + 724693.04214983 * 1i, -538428.45361942 + 675778.9257094 * 1i,
                                        740400.25975355 * 1i, 538428.45361942 + 675778.9257094 * 1i,
                                        -282791.34347149 + 724693.04214983 * 1i})};

        const Tensor& t = ampl.getTensor();
        // Check TEs.
        for (IndexType idx1 = 0; idx1 < 1; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 5; ++idx2) {
                double comRe = compareVals(t.element({idx1, idx2}).real(), amplEval.element({idx1, idx2}).real());
                double comIm = compareVals(t.element({idx1, idx2}).imag(), amplEval.element({idx1, idx2}).imag());
                EXPECT_NEAR(comRe, 1., 1e-4);
                EXPECT_NEAR(comIm, 1., 1e-4);
            }
        }
    }

    TEST(AmplD2starDstarDPiPiTest, evalTauDecNoRef) {

        // Momentum decs

        FourMomentum pDssmes{2461., 0, 0, 0};
        FourMomentum pDsmes{2047.62210483543, 243.33950461271, 280.828758699094, -120.736812915758};
        FourMomentum pPmes{413.377895164567, -243.33950461271, -280.828758699094, 120.736812915758};
        FourMomentum pDmes{1907.48119048676, 238.04888621381, 293.46680201094, -132.172069479921};
        FourMomentum pPDmes{140.140914348677, 5.290618398901, -12.63804331185, 11.4352565641629};
        //        FourMomentum kNuTau{1., 0.707106781186548, 0, 0.707106781186548};
        //        FourMomentum pTau{1., -0.70710678118655, 0, 0.707106781186548};


        // Evaluate at PS point
        AmplD2starDstarDPiPi ampl;
        ampl.eval(Particle(pDssmes, 425),
                  {Particle(pDsmes, 413), Particle(pPmes, -211), Particle(pDmes, 411), Particle(pPDmes, 111)}, {});

        // Compare to direct evaluation
        Tensor amplEval{"amplEval",
                        MD::makeVector({1, 5}, {FF_DSSD2STAR, SPIN_DSSD2STAR},
                                       {197180.87988408 - 752509.56797301 * 1i, 476320.40844743 + 720903.06185812 * 1i,
                                        740400.25975355 * 1i, -476320.40844743 + 720903.06185812 * 1i,
                                        -197180.87988408 - 752509.56797301 * 1i})};

        const Tensor& t = ampl.getTensor();
        // Check TEs.
        for (IndexType idx1 = 0; idx1 < 1; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 5; ++idx2) {
                double comRe = compareVals(t.element({idx1, idx2}).real(), amplEval.element({idx1, idx2}).real());
                double comIm = compareVals(t.element({idx1, idx2}).imag(), amplEval.element({idx1, idx2}).imag());
                EXPECT_NEAR(comRe, 1., 1e-4);
                EXPECT_NEAR(comIm, 1., 1e-4);
            }
        }
    }

} // namespace Hammer
