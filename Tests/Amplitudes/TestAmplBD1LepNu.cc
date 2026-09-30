///
/// @file  TestAmplBD1LepNu.cc
/// @brief Tests for AmplBD1LepNu
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "Hammer/Amplitudes/AmplBD1LepNu.hh"
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

    TEST(AmplBD1LepNuTest, eval) {

        // Momentum decs
        FourMomentum pBmes{5280., 0., 0., 0.};
        FourMomentum pD1mes{3100.34479166667, 0, 0, -1936.72321905187};
        FourMomentum kNuEle{133.38809812752, -39.108606047977, -67.738092688291, -108.04838015451};
        FourMomentum pPos{2046.26711020582, 39.108606047977, 67.738092688291, 2044.77159920637};

        Tensor FFvec{"FFvec", MD::makeVector({8}, {FF_BDSSD1}, {-0.8, -1.1, -0.4, -0.2, 1.3, 0.9, -0.8, 0.5})};

        // Evaluate at FF point
        AmplBD1LepNu ampl;
        ampl.eval(Particle(pBmes, 511), {Particle(pD1mes, -10413), Particle(kNuEle, 12), Particle(pPos, -11)}, {});
        auto t = ampl.getTensor();
        t.dot(FFvec, {FF_BDSSD1});

        // Compare to direct evaluation
        Tensor amplEval{"amplEval", MD::makeVector({3, 2, 2, 11}, {SPIN_DSSD1, SPIN_NUE, SPIN_EP, WILSON_BCTAUNU},
                                                   {-17.1599972377016,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    -611.607844419305,
                                                    0,
                                                    -17.1599972377016,
                                                    0,
                                                    -280163.90913228,
                                                    0,
                                                    -426689.702885378,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    -1.52078561437225e7,
                                                    0,
                                                    -426689.702885378,
                                                    0,
                                                    -1819.06801444123,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    94196.8322174657,
                                                    0,
                                                    2642.89837908581,
                                                    0,
                                                    -283.126864839335,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    611.607844419305,
                                                    0,
                                                    17.1599972377016,
                                                    0,
                                                    -7.0400546202737e6,
                                                    -110.356868206529,
                                                    3.23587736357672e6,
                                                    0,
                                                    3.23587736357672e6,
                                                    0,
                                                    -110.356868206529,
                                                    0,
                                                    -110.356868206529,
                                                    0,
                                                    650617.417520752,
                                                    0,
                                                    -4.80133933799938e6,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    -4.80133933799938e6,
                                                    0,
                                                    -4.80133933799938e6,
                                                    0,
                                                    -52.6573622006435,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    -4.80133933799938e6,
                                                    0,
                                                    -4.80133933799938e6,
                                                    0,
                                                    -52.6573622006435,
                                                    0,
                                                    0,
                                                    -3.23587736357672e6,
                                                    0,
                                                    -3.23587736357672e6,
                                                    0,
                                                    110.356868206529,
                                                    0,
                                                    110.356868206529,
                                                    0,
                                                    -650617.417520752,
                                                    -611.607844419305,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    -17.1599972377016,
                                                    0,
                                                    -611.607844419305,
                                                    0,
                                                    7.0400546202737e6,
                                                    0,
                                                    94196.8322174657,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    2642.89837908581,
                                                    0,
                                                    94196.8322174657,
                                                    0,
                                                    -283.126864839335,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    -426689.702885378,
                                                    0,
                                                    -1.52078561437225e7,
                                                    0,
                                                    -1819.06801444123,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    17.1599972377016,
                                                    0,
                                                    611.607844419305,
                                                    0,
                                                    280163.90913228})};
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
