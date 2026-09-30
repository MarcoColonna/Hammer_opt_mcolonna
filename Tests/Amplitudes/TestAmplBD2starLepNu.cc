///
/// @file  TestAmplBD2starLepNu.cc
/// @brief Tests for AmplBD2starLepNu
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "Hammer/Amplitudes/AmplBD2starLepNu.hh"
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

    TEST(AmplBD2starLepNuTest, eval) {

        // Momentum decs
        FourMomentum pBmes{5280., 0., 0., 0.};
        FourMomentum pD2starmes{3118.83721590909, 0, 0, -1915.88741301246};
        FourMomentum kNuEle{134.43152707986, -39.108606047977, -67.738092688291, -109.33390693174};
        FourMomentum pPos{2026.73125701105, 39.108606047977, 67.738092688291, 2025.22131994420};

        Tensor FFvec{"FFvec", MD::makeVector({8}, {FF_BDSSD2STAR}, {0.8, -1.9, 0.1, 1.4, -1.2, 0.9, 0.3, 0.2})};

        // Evaluate at FF point
        AmplBD2starLepNu ampl;
        ampl.eval(Particle(pBmes, 511), {Particle(pD2starmes, -415), Particle(kNuEle, 12), Particle(pPos, -11)}, {});
        auto t = ampl.getTensor();
        t.dot(FFvec, {FF_BDSSD2STAR});

        // Compare to direct evaluation
        Tensor amplEval{"amplEval", MD::makeVector({5, 2, 2, 11}, {SPIN_DSSD2STAR, SPIN_NUE, SPIN_EP, WILSON_BCTAUNU},
                                                   {0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    -449.57260366903,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    153.199329654279,
                                                    0,
                                                    -449.57260366903,
                                                    0,
                                                    -454783.425660257,
                                                    0,
                                                    -1.11787897181875e7,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    3.80935821532678e6,
                                                    0,
                                                    -1.11787897181875e7,
                                                    0,
                                                    -2952.84994301669,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    -23595.0072958496,
                                                    0,
                                                    69240.9613509596,
                                                    0,
                                                    -494.214814249867,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    -153.199329654279,
                                                    0,
                                                    449.57260366903,
                                                    0,
                                                    -1.22888348600966e7,
                                                    -555.694712850594,
                                                    2.01812846362868e6,
                                                    0,
                                                    -2.01812846362868e6,
                                                    0,
                                                    555.694712850594,
                                                    0,
                                                    -555.694712850594,
                                                    0,
                                                    -2.87019721727755e7,
                                                    0,
                                                    1.72783279944396e6,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    -1.72783279944396e6,
                                                    0,
                                                    1.72783279944396e6,
                                                    0,
                                                    2322.97830318449,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    -1.72783279944396e6,
                                                    0,
                                                    1.72783279944396e6,
                                                    0,
                                                    -2322.97830318449,
                                                    0,
                                                    0,
                                                    -2.01812846362868e6,
                                                    0,
                                                    2.01812846362868e6,
                                                    0,
                                                    -555.694712850594,
                                                    0,
                                                    555.694712850594,
                                                    0,
                                                    -2.87019721727755e7,
                                                    -153.199329654279,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    449.57260366903,
                                                    0,
                                                    -153.199329654279,
                                                    0,
                                                    -1.22888348600966e7,
                                                    0,
                                                    23595.0072958496,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    -69240.9613509596,
                                                    0,
                                                    23595.0072958496,
                                                    0,
                                                    494.214814249867,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    1.11787897181875e7,
                                                    0,
                                                    -3.80935821532678e6,
                                                    0,
                                                    2952.84994301669,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    -449.57260366903,
                                                    0,
                                                    153.199329654279,
                                                    0,
                                                    -454783.425660257,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0})};
        amplEval *= GFermi;


        // Check TEs. Care that nuTau spin (second index) is double-copied; third index is ref qu.no.
        for (IndexType idx1 = 0; idx1 < 2; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 2; ++idx2) {
                for (IndexType idx3 = 0; idx3 < 11; ++idx3) {
                    for (IndexType idxSpin = 0; idxSpin < 5; ++idxSpin) {
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
