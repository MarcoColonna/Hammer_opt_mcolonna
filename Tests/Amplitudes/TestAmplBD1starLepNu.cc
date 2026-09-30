///
/// @file  TestAmplBD1starLepNu.cc
/// @brief Tests for AmplBD1starLepNu
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "Hammer/Amplitudes/AmplBD1starLepNu.hh"
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

    TEST(AmplBD1starLepNuTest, eval) {

        // Momentum decs
        FourMomentum pBmes{5280., 0., 0., 0.};
        FourMomentum pD1starmes{3103.09933712121, 0, 0, -1933.62263537695};
        FourMomentum kNuEle{133.54203053243, -39.108606047977, -67.738092688291, -108.23835572781};
        FourMomentum pPos{2043.35863234636, 39.108606047977, 67.738092688291, 2041.86099110476};

        Tensor FFvec{"FFvec", MD::makeVector({8}, {FF_BDSSD1STAR}, {0.8, 0.5, 0.1, -1.1, 1.3, -0.9, 0.8, 0.2})};

        // Evaluate at FF point
        AmplBD1starLepNu ampl;
        ampl.eval(Particle(pBmes, 511), {Particle(pD1starmes, -20413), Particle(kNuEle, 12), Particle(pPos, -11)}, {});
        auto t = ampl.getTensor();
        t.dot(FFvec, {FF_BDSSD1STAR});

        // Compare to direct evaluation
        Tensor amplEval{"amplEval", MD::makeVector({3, 2, 2, 11}, {SPIN_DSSD1STAR, SPIN_NUE, SPIN_EP, WILSON_BCTAUNU},
                                                   {439.459799918904,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    -153.30229909109,
                                                    0,
                                                    439.459799918904,
                                                    0,
                                                    280376.705975637,
                                                    0,
                                                    1.09273310980195e7,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    -3.81191858860603e6,
                                                    0,
                                                    1.09273310980195e7,
                                                    0,
                                                    1820.44967681354,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    23610.8661420879,
                                                    0,
                                                    -67683.4370536655,
                                                    0,
                                                    281.508326908167,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    153.30229909109,
                                                    0,
                                                    -439.459799918904,
                                                    0,
                                                    6.99980907364613e6,
                                                    -410.913359227148,
                                                    3.22670099196555e6,
                                                    0,
                                                    3.22670099196555e6,
                                                    0,
                                                    -410.913359227148,
                                                    0,
                                                    -410.913359227148,
                                                    0,
                                                    -9.55370833084922e6,
                                                    0,
                                                    -1.08548847933052e6,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    -1.08548847933052e6,
                                                    0,
                                                    -1.08548847933052e6,
                                                    0,
                                                    773.224119719767,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    -1.08548847933052e6,
                                                    0,
                                                    -1.08548847933052e6,
                                                    0,
                                                    773.224119719767,
                                                    0,
                                                    0,
                                                    -3.22670099196555e6,
                                                    0,
                                                    -3.22670099196555e6,
                                                    0,
                                                    410.913359227148,
                                                    0,
                                                    410.913359227148,
                                                    0,
                                                    9.55370833084922e6,
                                                    -153.30229909109,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    439.459799918904,
                                                    0,
                                                    -153.30229909109,
                                                    0,
                                                    -6.99980907364613e6,
                                                    0,
                                                    23610.8661420879,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    -67683.4370536655,
                                                    0,
                                                    23610.8661420879,
                                                    0,
                                                    281.508326908167,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    1.09273310980195e7,
                                                    0,
                                                    -3.81191858860603e6,
                                                    0,
                                                    1820.44967681354,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    -439.459799918904,
                                                    0,
                                                    153.30229909109,
                                                    0,
                                                    -280376.705975637})};
        amplEval *= GFermi;


        // Check TEs. Care that nuTau spin (second index) is double-copied; third index is ref qu.no.
        for (IndexType idx1 = 0; idx1 < 2; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 2; ++idx2) {
                for (IndexType idx3 = 0; idx3 < 11; ++idx3) {
                    for (IndexType idxSpin = 0; idxSpin < 3; ++idxSpin) {
                        double comRe = compareVals(t.element({idx3, idxSpin, idx1, idx1, idx2}).real(),
                                                   amplEval.element({idxSpin, idx1, idx2, idx3}).real());
                        double comIm = compareVals(t.element({idx3, idxSpin, idx1, idx1, idx2}).imag(),
                                                   amplEval.element({idxSpin, idx1, idx2, idx3}).imag());
                        EXPECT_NEAR(comRe, 1., 1e-4);
                        EXPECT_NEAR(comIm, 1., 1e-4);
                    }
                }
            }
        }
    }

} // namespace Hammer
