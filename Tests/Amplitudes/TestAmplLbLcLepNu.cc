///
/// @file  TestAmplLbLcLepNu.cc
/// @brief Tests for AmplLbLcLepNu
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "Hammer/Amplitudes/AmplLbLcLepNu.hh"
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

    TEST(AmplLbLcLepNuTest, eval) {

        // Momentum decs
        FourMomentum pLbmes{5620., 0., 0., 0.};
        FourMomentum pLcmes{2514.6000000000004, 380.2511160960518, -878.708483534587, -425.0854616098454};
        FourMomentum pTau{2456.13415321897, -581.6942711749232, 1418.910848905593, 723.6624190525239};
        FourMomentum kNuTau{649.2658467810297, 201.44315507887137, -540.2023653710064, -298.5769574426786};

        // Create FF
        Tensor FFvec{"FFvec",
                     MD::makeVector({12}, {FF_LBLC}, {0.8, 0.9, 1.0, 0.3, 0.3, 1.1, 0.4, 0.4, 0.9, 0.2, 0.2, 0.2})};

        // Evaluate
        AmplLbLcLepNu ampl;
        ampl.eval(Particle(pLbmes, -5122), {Particle(pLcmes, -4122), Particle(kNuTau, 16), Particle(pTau, -15)}, {});
        auto t = ampl.getTensor();
        t.dot(FFvec, {FF_LBLC});

        // Compare to direct evaluation
        Tensor amplEval{"amplEval",
                        MD::makeVector({11, 2, 2, 2, 2}, {WILSON_BCTAUNU, SPIN_LB, SPIN_LC, SPIN_NUTAU, SPIN_TAUP},
                                       {2.3495692863459468e7,
                                        -237366.33607951,
                                        0.,
                                        0.,
                                        1.0314384888211903e6,
                                        6.353336530176189e7,
                                        0.,
                                        0.,
                                        689949.6763784066,
                                        -30319.360069244118,
                                        0.,
                                        0.,
                                        2.840619479781007e6,
                                        2.714033714720155e6,
                                        0.,
                                        0.,
                                        -2.4028485902123664e7,
                                        0.,
                                        0.,
                                        0.,
                                        0.,
                                        0.,
                                        0.,
                                        0.,
                                        0.,
                                        0.,
                                        0.,
                                        0.,
                                        -1.4556146651755374e7,
                                        0.,
                                        0.,
                                        0.,
                                        0.,
                                        0.,
                                        0.,
                                        2.4028485902123664e7,
                                        0.,
                                        0.,
                                        0.,
                                        0.,
                                        0.,
                                        0.,
                                        0.,
                                        0.,
                                        0.,
                                        0.,
                                        0.,
                                        1.4556146651755374e7,
                                        -1.4556146651755374e7,
                                        0.,
                                        0.,
                                        0.,
                                        0.,
                                        0.,
                                        0.,
                                        0.,
                                        0.,
                                        0.,
                                        0.,
                                        0.,
                                        -2.4028485902123664e7,
                                        0.,
                                        0.,
                                        0.,
                                        0.,
                                        0.,
                                        0.,
                                        1.4556146651755374e7,
                                        0.,
                                        0.,
                                        0.,
                                        0.,
                                        0.,
                                        0.,
                                        0.,
                                        0.,
                                        0.,
                                        0.,
                                        0.,
                                        2.4028485902123664e7,
                                        2.840619479781007e6,
                                        2.714033714720155e6,
                                        0.,
                                        0.,
                                        -689949.6763784066,
                                        -4.249872901221634e7,
                                        0.,
                                        0.,
                                        -1.0314384888211903e6,
                                        45325.84911989305,
                                        0.,
                                        0.,
                                        2.3495692863459468e7,
                                        -237366.33607951,
                                        0.,
                                        0.,
                                        0.,
                                        0.,
                                        2.714033714720155e6,
                                        -2.840619479781007e6,
                                        0.,
                                        0.,
                                        30319.360069244118,
                                        689949.6763784066,
                                        0.,
                                        0.,
                                        -6.353336530176189e7,
                                        1.0314384888211903e6,
                                        0.,
                                        0.,
                                        -237366.33607951,
                                        -2.3495692863459468e7,
                                        2.3495692863459468e7,
                                        -237366.33607951,
                                        0.,
                                        0.,
                                        1.0314384888211903e6,
                                        6.353336530176189e7,
                                        0.,
                                        0.,
                                        689949.6763784066,
                                        -30319.360069244118,
                                        0.,
                                        0.,
                                        2.840619479781007e6,
                                        2.714033714720155e6,
                                        0.,
                                        0.,
                                        0.,
                                        0.,
                                        -237366.33607951,
                                        -2.3495692863459468e7,
                                        0.,
                                        0.,
                                        -45325.84911989305,
                                        -1.0314384888211903e6,
                                        0.,
                                        0.,
                                        4.249872901221634e7,
                                        -689949.6763784066,
                                        0.,
                                        0.,
                                        2.714033714720155e6,
                                        -2.840619479781007e6,
                                        9.678106248088317e7,
                                        -3.144646257544414e6,
                                        0.,
                                        0.,
                                        2.810061196688898e6,
                                        6.394596748846743e7,
                                        0.,
                                        0.,
                                        7.20354079951186e6,
                                        -116946.57131918453,
                                        0.,
                                        0.,
                                        -7.660221787240809e7,
                                        2.4889877376542683e6,
                                        0.,
                                        0.,
                                        0.,
                                        0.,
                                        2.4889877376542683e6,
                                        7.660221787240809e7,
                                        0.,
                                        0.,
                                        116946.57131918453,
                                        7.20354079951186e6,
                                        0.,
                                        0.,
                                        -6.394596748846743e7,
                                        2.810061196688898e6,
                                        0.,
                                        0.,
                                        -3.144646257544414e6,
                                        -9.678106248088317e7})};
        amplEval *= GFermi;


        // Check TEs. Care that nuTau spin is double-copied; second entry is ref qu.no.
        for (IndexType idx1 = 0; idx1 < 3; ++idx1) {
            for (IndexType idx3 = 0; idx3 < 2; ++idx3) {
                for (IndexType idx4 = 0; idx4 < 2; ++idx4) {
                    for (IndexType idx5 = 0; idx5 < 2; ++idx5) {
                        for (IndexType idx6 = 0; idx6 < 2; ++idx6) {
                            //                cout << idx1 << " " << idx3 << " " << idx4 << " " << idx5 << " " <<
                            //                idx6 << '\n'; cout << ampl.element({idx1, idx3, idx4, idx5, idx5,
                            //                idx6}).real() << '\n'; cout << amplEval.element({idx1, idx3, idx4,
                            //                idx5, idx6}).real() << '\n' << '\n';
                            double comRe = compareVals(t.element({idx1, idx3, idx4, idx5, idx5, idx6}).real(),
                                                       amplEval.element({idx1, idx3, idx4, idx5, idx6}).real());
                            double comIm = compareVals(t.element({idx1, idx3, idx4, idx5, idx5, idx6}).imag(),
                                                       amplEval.element({idx1, idx3, idx4, idx5, idx6}).imag());
                            EXPECT_NEAR(comRe, 1., 1e-4);
                            EXPECT_NEAR(comIm, 1., 1e-4);
                        }
                    }
                }
            }
        }
    }


} // namespace Hammer
