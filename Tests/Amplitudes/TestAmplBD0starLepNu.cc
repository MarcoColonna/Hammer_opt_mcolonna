///
/// @file  TestAmplBD0starLepNu.cc
/// @brief Tests for AmplBD0starLepNu
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "Hammer/Amplitudes/AmplBD0starLepNu.hh"
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

    TEST(AmplBD0starLepNuTest, eval) {

        // Momentum decs
        FourMomentum pBmes{5280., 0., 0., 0.};
        FourMomentum pD0starmes{3046.25000000000, 0, 0, -1997.40808612061};
        FourMomentum kNuEle{130.46662689653, -39.108606047977, -67.738092688291, -104.42034507554};
        FourMomentum pPos{2103.28337310347, 39.108606047977, 67.738092688291, 2101.82843119615};

        Tensor FFvec{"FFvec", MD::makeVector({4}, {FF_BDSSD0STAR}, {1.2, 0.2, 1.5, 1.7})};

        // Evaluate at FF point
        AmplBD0starLepNu ampl;
        ampl.eval(Particle(pBmes, 511), {Particle(pD0starmes, -10411), Particle(kNuEle, 12), Particle(pPos, -11)}, {});
        auto t = ampl.getTensor();
        t.dot(FFvec, {FF_BDSSD0STAR});

        // Compare to direct evaluation
        Tensor amplEval{"amplEval", MD::makeVector({2, 2, 11}, {SPIN_NUE, SPIN_EP, WILSON_BCTAUNU},
                                                   {579.443579381093,
                                                    -5.91394207507701e6,
                                                    0,
                                                    5.91394207507701e6,
                                                    0,
                                                    -579.443579381093,
                                                    0,
                                                    579.443579381093,
                                                    0,
                                                    2.87450870290811e7,
                                                    0,
                                                    -1.97777881405426e6,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    1.97777881405426e6,
                                                    0,
                                                    -1.97777881405426e6,
                                                    0,
                                                    -2326.46778032354,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    1.97777881405426e6,
                                                    0,
                                                    -1.97777881405426e6,
                                                    0,
                                                    2326.46778032354,
                                                    0,
                                                    0,
                                                    5.91394207507701e6,
                                                    0,
                                                    -5.91394207507701e6,
                                                    0,
                                                    579.443579381093,
                                                    0,
                                                    -579.443579381093,
                                                    0,
                                                    2.87450870290811e7})};
        amplEval *= GFermi;


        // Check TEs. Care that nuTau spin (second index) is double-copied; third index is ref qu.no.
        for (IndexType idx1 = 0; idx1 < 2; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 2; ++idx2) {
                for (IndexType idx3 = 0; idx3 < 11; ++idx3) {
                    for (IndexType idxSpin = 0; idxSpin < 1; ++idxSpin) {
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
    }

} // namespace Hammer
