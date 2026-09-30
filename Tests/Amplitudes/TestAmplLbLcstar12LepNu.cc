///
/// @file  TestAmplLbLcstar12LepNu.cc
/// @brief Tests for AmplLbLcstar12LepNu
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "Hammer/Amplitudes/AmplLbLcstar12LepNu.hh"
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

    TEST(AmplLbLcstar12LepNuTest, eval) {

        // Momentum decs
        FourMomentum pLbmes{5620., 0., 0., 0.};
        FourMomentum pLcmes{3053.24065836299, 0, 0, -1608.80499684731};
        FourMomentum kNuTau{103.003318909325, -16.4797557212763, -28.5437742055744, -97.5876749994806};
        FourMomentum pTau{2463.75602272769, 16.4797557212763, 28.5437742055744, 1706.39267184679};

        // Create FF
        Tensor FFvec{"FFvec", MD::makeVector({12}, {FF_LBLCSTAR12},
                                             {1.1, 0.8, 0.9, 1.3, 0.3, 1.1, 0.8, 0.2, 1.2, 0.8, 0.2, 0.1})};

        // Evaluate
        AmplLbLcstar12LepNu ampl;
        ampl.eval(Particle(pLbmes, -5122), {Particle(pLcmes, -14122), Particle(kNuTau, 16), Particle(pTau, -15)}, {});
        auto t = ampl.getTensor();
        t.dot(FFvec, {FF_LBLCSTAR12});

        // Compare to direct evaluation
        Tensor amplEval{"amplEval", MD::makeVector({2, 2, 2, 2, 11},
                                                   {SPIN_LB, SPIN_LCSTAR12, SPIN_NUTAU, SPIN_TAUP, WILSON_BCTAUNU},
                                                   {-1.23016282612139e7,
                                                    -1.15166839763221e7,
                                                    0,
                                                    5.0344932952568e6,
                                                    0,
                                                    3.16323537831836e6,
                                                    0,
                                                    -1.23016282612139e7,
                                                    0,
                                                    1.03157883951303e8,
                                                    0,
                                                    -2.20784785978667e6,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    5.16770023196572e6,
                                                    0,
                                                    -2.20784785978667e6,
                                                    0,
                                                    -1.45157056982426e7,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    5.16770023196572e6,
                                                    0,
                                                    -2.20784785978667e6,
                                                    0,
                                                    -1.05844978074403e7,
                                                    0,
                                                    0,
                                                    1.15166839763221e7,
                                                    0,
                                                    -5.0344932952568e6,
                                                    0,
                                                    -3.16323537831836e6,
                                                    0,
                                                    1.23016282612139e7,
                                                    0,
                                                    -7.52202076289644e7,
                                                    -1.74436985460329e6,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    -843468.046294422,
                                                    0,
                                                    -1.74436985460329e6,
                                                    0,
                                                    5.48621542588267e6,
                                                    0,
                                                    -2.49477398424194e7,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    -1.20631650041514e7,
                                                    0,
                                                    -2.49477398424194e7,
                                                    0,
                                                    6.19315459279503e7,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    74718.745309588,
                                                    0,
                                                    154525.269172237,
                                                    0,
                                                    -2.55232553777963e6,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    843468.046294422,
                                                    0,
                                                    1.74436985460329e6,
                                                    0,
                                                    -3.65030118708228e7,
                                                    -843468.046294422,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    -1.74436985460329e6,
                                                    0,
                                                    -843468.046294422,
                                                    0,
                                                    3.65030118708228e7,
                                                    0,
                                                    74718.745309588,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    154525.269172237,
                                                    0,
                                                    74718.745309588,
                                                    0,
                                                    -2.55232553777963e6,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    -2.49477398424194e7,
                                                    0,
                                                    -1.20631650041514e7,
                                                    0,
                                                    6.19315459279503e7,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    1.74436985460329e6,
                                                    0,
                                                    843468.046294422,
                                                    0,
                                                    -5.48621542588267e6,
                                                    -3.16323537831836e6,
                                                    -5.0344932952568e6,
                                                    0,
                                                    1.15166839763221e7,
                                                    0,
                                                    1.23016282612139e7,
                                                    0,
                                                    -3.16323537831836e6,
                                                    0,
                                                    -7.52202076289644e7,
                                                    0,
                                                    -5.16770023196572e6,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    2.20784785978667e6,
                                                    0,
                                                    -5.16770023196572e6,
                                                    0,
                                                    1.05844978074403e7,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    2.20784785978667e6,
                                                    0,
                                                    -5.16770023196572e6,
                                                    0,
                                                    1.45157056982426e7,
                                                    0,
                                                    0,
                                                    5.0344932952568e6,
                                                    0,
                                                    -1.15166839763221e7,
                                                    0,
                                                    -1.23016282612139e7,
                                                    0,
                                                    3.16323537831836e6,
                                                    0,
                                                    1.03157883951303e8})};
        amplEval *= GFermi;


        // Check TEs. Care that nuTau spin is double-copied; second entry is ref qu.no.
        for (IndexType idx1 = 0; idx1 < 2; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 2; ++idx2) {
                for (IndexType idx3 = 0; idx3 < 2; ++idx3) {
                    for (IndexType idx4 = 0; idx4 < 2; ++idx4) {
                        for (IndexType idx5 = 0; idx5 < 11; ++idx5) {
                            double comRe = compareVals(t.element({idx5, idx1, idx2, idx3, idx3, idx4}).real(),
                                                       amplEval.element({idx1, idx2, idx3, idx4, idx5}).real());
                            double comIm = compareVals(t.element({idx5, idx1, idx2, idx3, idx3, idx4}).imag(),
                                                       amplEval.element({idx1, idx2, idx3, idx4, idx5}).imag());
                            EXPECT_NEAR(comRe, 1., 1e-4);
                            EXPECT_NEAR(comIm, 1., 1e-4);
                        }
                    }
                }
            }
        }
    }


} // namespace Hammer
