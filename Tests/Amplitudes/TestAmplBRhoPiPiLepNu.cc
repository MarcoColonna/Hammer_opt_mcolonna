///
/// @file  TestAmplBRhoPiPiLepNu.cc
/// @brief Tests for AmplBRhoPiPiLepNu
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "Hammer/Amplitudes/AmplBRhoPiPiLepNu.hh"
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

    TEST(AmplBRhoPiPiLepNuTest, eval) {

        // Momentum decs
        FourMomentum pBmes{17838.5064935065, 0, 0, 17039.1875956242};
        FourMomentum pRhomes{770., 0, 0, 0};
        FourMomentum kNuMu{198.6188302152, 31.769954441031, 97.777865759968, 169.9399855193};
        FourMomentum pMu{16869.8876632913, -31.769954441031, -97.777865759968, 16869.2476101049};
        FourMomentum pPpmes{385.17987012987, 89.80622624339, 155.548946689577, 311.097893379154};
        FourMomentum pPmmes{384.82012987013, -89.80622624339, -155.548946689577, -311.097893379154};

        //Fix sign convention in T2 and T23 definitions, FF indices 6 and 7. Amounts to switching signs in FFvec for purposes of this test.
        const double signfix = -1.;
        Tensor FFvec{"FFvec", MD::makeVector({8}, {FF_BRHO}, {0.7, 1.1, 0.9, 0.6, -0.4, 1.1, signfix*0.3, signfix*0.8})};

        // Evaluate at FF point
        AmplBRhoPiPiLepNu ampl;
        ampl.eval(Particle(pBmes, 521),
                  {Particle(pRhomes, 113), Particle(kNuMu, 14), Particle(pMu, -13), Particle(pPpmes, 211),
                   Particle(pPmmes, -211)},
                  {});
        auto t = ampl.getTensor();
        t.dot(FFvec, {FF_BRHO});

        // Compare to direct evaluation
        Tensor amplEval{"amplEval", MD::makeVector({2, 2, 11}, {SPIN_NUMU, SPIN_MUP, WILSON_BUMUNU},
                                                   {-1.90482750530001e7 - 61686.192496 * 1i,
                                                    5.88493243662213e7,
                                                    0,
                                                    -5.88493243662214e7,
                                                    0,
                                                    1.90482750530001e7 - 61686.192496 * 1i,
                                                    0,
                                                    -1.90482750530001e7 - 61686.192496 * 1i,
                                                    0,
                                                    -2.25440467931482e7 + 1.42171986514662e7 * 1i,
                                                    0,
                                                    -9.84125958220611e6 - 4.9138175706794e6 * 1i,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    -1.67461015654366e7 - 614005.6399094 * 1i,
                                                    0,
                                                    -9.84125958220611e6 - 4.9138175706794e6 * 1i,
                                                    0,
                                                    4.40758556203657e7 + 9.0592753288692e6 * 1i,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    9.8412595822061e6 - 4.9138175706794e6 * 1i,
                                                    0,
                                                    1.67461015654366e7 - 614005.6399094 * 1i,
                                                    0,
                                                    -4.40758556203657e7 + 9.0592753288692e6 * 1i,
                                                    0,
                                                    0,
                                                    -5.88493243662214e7,
                                                    0,
                                                    5.88493243662213e7,
                                                    0,
                                                    -1.90482750530001e7 + 61686.192496 * 1i,
                                                    0,
                                                    1.90482750530001e7 + 61686.192496 * 1i,
                                                    0,
                                                    -2.25440467931482e7 - 1.42171986514662e7 * 1i})};
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
