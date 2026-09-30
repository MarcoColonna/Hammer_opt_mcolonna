///
/// @file  TestAmplBOmega3PiLepNu.cc
/// @brief Tests for AmplBOmega3PiLepNu
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "Hammer/Amplitudes/AmplBOmega3PiLepNu.hh"
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

    TEST(AmplBOmega3PiLepNuTest, eval) {

        // Momentum decs
        FourMomentum pBmes{17576.6777493606, 0, 0, 16764.8799788388};
        FourMomentum pOmegames{782., 0, 0, 0};
        FourMomentum kNuMu{195.8915304338, 31.769954441031, 97.777865759968, 166.7442672416};
        FourMomentum pMu{16598.7862189268, -31.769954441031, -97.777865759968, 16598.1357115972};
        FourMomentum pPpmes{200., 35.9504867839088, 62.2680696665635, 124.536139333127};
        FourMomentum pPmmes{290.108282131125, -226.41984524796, -109.038050665746, 44.312216667374};
        FourMomentum pPzmes{291.891717868875, 190.469358464051, 46.7699809991826, -168.848356000501};

        //Fix sign convention in T2 and T23 definitions, FF indices 6 and 7. Amounts to switching signs in FFvec for purposes of this test.
        const double signfix = -1.;
        Tensor FFvec{"FFvec", MD::makeVector({8}, {FF_BOMEGA}, {0.7, 1.1, 0.9, 0.6, -0.4, 1.1, signfix*0.3, signfix*0.8})};

        // Evaluate at FF point
        AmplBOmega3PiLepNu ampl;
        ampl.eval(Particle(pBmes, 521),
                  {Particle(pOmegames, 223), Particle(kNuMu, 14), Particle(pMu, -13), Particle(pPpmes, 211),
                   Particle(pPmmes, -211), Particle(pPzmes, 111)},
                  {});
        auto t = ampl.getTensor();
        t.dot(FFvec, {FF_BOMEGA});

        // Compare to direct evaluation
        Tensor amplEval{"amplEval", MD::makeVector({2, 2, 11}, {SPIN_NUMU, SPIN_MUP, WILSON_BUMUNU},
                                                   {-1.5693489064268e10 + 2.27052437614954e11 * 1i,
                                                    -7.2937587704671e11 * 1i,
                                                    0,
                                                    7.2937587704671e11 * 1i,
                                                    0,
                                                    -1.5693489064268e10 - 2.27052437614954e11 * 1i,
                                                    0,
                                                    -1.5693489064268e10 + 2.27052437614954e11 * 1i,
                                                    0,
                                                    3.62414974124153e12 + 1.49719211237819e12 * 1i,
                                                    0,
                                                    -1.252709095505e12 - 1.35296540191613e12 * 1i,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    -1.53615879247846e11 + 1.07116811165e8 * 1i,
                                                    0,
                                                    -1.252709095505e12 - 1.35296540191613e12 * 1i,
                                                    0,
                                                    2.30937055446595e12 + 2.18815622786056e12 * 1i,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    -1.252709095505e12 + 1.35296540191613e12 * 1i,
                                                    0,
                                                    -1.53615879247846e11 - 1.07116811165e8 * 1i,
                                                    0,
                                                    2.30937055446595e12 - 2.18815622786056e12 * 1i,
                                                    0,
                                                    0,
                                                    7.2937587704671e11 * 1i,
                                                    0,
                                                    -7.2937587704671e11 * 1i,
                                                    0,
                                                    1.5693489064268e10 + 2.27052437614954e11 * 1i,
                                                    0,
                                                    1.5693489064268e10 - 2.27052437614954e11 * 1i,
                                                    0,
                                                    -3.62414974124153e12 + 1.49719211237819e12 * 1i})};
        amplEval *= GFermi;


        // Check TEs. Care that nuTau spin (second index) is double-copied; third index is ref qu.no.
        for (IndexType idx1 = 0; idx1 < 2; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 2; ++idx2) {
                for (IndexType idx3 = 0; idx3 < 11; ++idx3) {
                    double comRe = compareVals(t.element({idx3, idx1, idx1, idx2}).real(),
                                               amplEval.element({idx1, idx2, idx3}).real());
                    double comIm = compareVals(t.element({idx3, idx1, idx1, idx2}).imag(),
                                               amplEval.element({idx1, idx2, idx3}).imag());
                    EXPECT_NEAR(comRe, 1., 1e-4);
                    EXPECT_NEAR(comIm, 1., 1e-4);
                }
            }
        }
    }

} // namespace Hammer
