///
/// @file  TestAmplBDstarDGamLepNu.cc
/// @brief Tests for AmplBDstarDGamLepNu
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "Hammer/Amplitudes/AmplBDstarDGamLepNu.hh"
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

    TEST(AmplBDstarDGamLepNuTest, eval) {

        // Momentum decs
        FourMomentum pBmes{5280., 0., 0., 0.};
        FourMomentum pDstarmes{2258.9451, -542.95362175877925474, -718.43449427445361377, 513.65086086722812801};
        FourMomentum kNuTau{763.45512706233784082, 59.818100429037400761, 266.26625745054465575, 713.01318784585764875};
        FourMomentum pTau{2257.5997729376621592, 483.13552132974185397, 452.16823682390895801, -1226.6640487130857768};

        FourMomentum pDmes{2081.854507679764963, -560.22997131734918369, -646.12961981787364724, 352.91937433120433085};
        FourMomentum kPhoton{177.09059232023503702, 17.276349558569928956, -72.304874456579966531,
                             160.73148653602379716};


        // Create FF and kin factors
        double Mb = pBmes.mass();
        double Mds = pDstarmes.mass();
        double Mds2 = pDstarmes.mass2();
        double Md2 = pDmes.mass2();
        // double Mt = pTau.mass();
        double Sqq = Mb * Mb + Mds2 - 2. * (pBmes * pDstarmes);
        double eMuOnSqrt2MGam = sqrt(48. * pi * Mds2 / pow(Mds2 - Md2, 3.));

        // Form Factors Set Ff with model, equivalent to h_A1. Use HQET LO for others.
        double Ff = 2.;
        // Axial and vector FFs
        double MbdsSqq = pow(Mb + Mds, 2.) - Sqq;
        double Fg = Ff / MbdsSqq;
        double Fp = -Ff / MbdsSqq;
        double Fm = Ff / MbdsSqq;
        // NP form factors
        double Fs = -2 * Mds * Ff / MbdsSqq;
        double Fpt = -(Mb + Mds) * Ff / MbdsSqq;
        double Fmt = (Mb - Mds) * Ff / MbdsSqq;
        double Fzt = 0 * Ff;

        Tensor FFvec{"FFvec", MD::makeVector({8}, {FF_BDSTAR}, {Fs, Ff, Fg, Fm, Fp, Fzt, Fmt, Fpt})};


        // Evaluate at FF point
        AmplBDstarDGamLepNu ampl;
        ampl.eval(Particle(pBmes, 511),
                  {Particle(pDstarmes, -423), Particle(kNuTau, 16), Particle(pTau, -15), Particle(pDmes, -421),
                   Particle(kPhoton, 22)},
                  {});
        auto t = ampl.getTensor();
        t.dot(FFvec, {FF_BDSTAR});

        // Compare to direct evaluation
        Tensor amplEval{"amplEval", MD::makeVector({2, 2, 2, 11}, {SPIN_GAMMA, SPIN_NUTAU, SPIN_TAUP, WILSON_BCTAUNU},
                                                   {(-349496831.573 + 84129129.8695 * 1i),
                                                    279769540.919 * 1i,
                                                    0. * 1i,
                                                    -279769540.919 * 1i,
                                                    0. * 1i,
                                                    (409222464.192 - 348242839.295 * 1i),
                                                    0. * 1i,
                                                    (-349496831.573 + 84129129.8695 * 1i),
                                                    0. * 1i,
                                                    (-3286609166.59 + 1195921780.81 * 1i),
                                                    0. * 1i,
                                                    (-259842228.58 + 276053670.924 * 1i),
                                                    0. * 1i,
                                                    0. * 1i,
                                                    0. * 1i,
                                                    0. * 1i,
                                                    (-24675028.3509 - 357916465.154 * 1i),
                                                    0. * 1i,
                                                    (-259842228.58 + 276053670.924 * 1i),
                                                    0. * 1i,
                                                    (575050284.426 + 349820840.157 * 1i),
                                                    0. * 1i,
                                                    0. * 1i,
                                                    0. * 1i,
                                                    0. * 1i,
                                                    0. * 1i,
                                                    0. * 1i,
                                                    0. * 1i,
                                                    (765121675.258 - 2510458282.6 * 1i),
                                                    0. * 1i,
                                                    (-372536928.693 + 2114434040.66 * 1i),
                                                    0. * 1i,
                                                    (380606049.944 - 4575844547.9 * 1i),
                                                    0. * 1i,
                                                    0. * 1i,
                                                    -279769540.919 * 1i,
                                                    0. * 1i,
                                                    279769540.919 * 1i,
                                                    0. * 1i,
                                                    (-409222464.192 + 348242839.295 * 1i),
                                                    0. * 1i,
                                                    (349496831.573 - 84129129.8695 * 1i),
                                                    0. * 1i,
                                                    (-2306382398.81 - 3138755247.53 * 1i),
                                                    (-409222464.192 - 348242839.295 * 1i),
                                                    -279769540.919 * 1i,
                                                    0. * 1i,
                                                    279769540.919 * 1i,
                                                    0. * 1i,
                                                    (349496831.573 + 84129129.8695 * 1i),
                                                    0. * 1i,
                                                    (-409222464.192 - 348242839.295 * 1i),
                                                    0. * 1i,
                                                    (-2306382398.81 + 3138755247.53 * 1i),
                                                    0. * 1i,
                                                    (-765121675.258 - 2510458282.6 * 1i),
                                                    0. * 1i,
                                                    0. * 1i,
                                                    0. * 1i,
                                                    0. * 1i,
                                                    (372536928.693 + 2114434040.66 * 1i),
                                                    0. * 1i,
                                                    (-765121675.258 - 2510458282.6 * 1i),
                                                    0. * 1i,
                                                    (-380606049.944 - 4575844547.9 * 1i),
                                                    0. * 1i,
                                                    0. * 1i,
                                                    0. * 1i,
                                                    0. * 1i,
                                                    0. * 1i,
                                                    0. * 1i,
                                                    0. * 1i,
                                                    (259842228.58 + 276053670.924 * 1i),
                                                    0. * 1i,
                                                    (24675028.3509 - 357916465.154 * 1i),
                                                    0. * 1i,
                                                    (-575050284.426 + 349820840.157 * 1i),
                                                    0. * 1i,
                                                    0. * 1i,
                                                    279769540.919 * 1i,
                                                    0. * 1i,
                                                    -279769540.919 * 1i,
                                                    0. * 1i,
                                                    (-349496831.573 - 84129129.8695 * 1i),
                                                    0. * 1i,
                                                    (409222464.192 + 348242839.295 * 1i),
                                                    0. * 1i,
                                                    (-3286609166.59 - 1195921780.81 * 1i)})};
        amplEval *= GFermi * eMuOnSqrt2MGam;


        // Check TEs. Care that nuTau spin (third index) is double-copied; fourth index is ref qu.no.
        for (IndexType idx1 = 0; idx1 < 2; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 2; ++idx2) {
                for (IndexType idx3 = 0; idx3 < 2; ++idx3) {
                    for (IndexType idx4 = 0; idx4 < 11; ++idx4) {
                        double comRe = compareVals(t.element({idx4, idx1, idx2, idx2, idx3}).real(),
                                                   amplEval.element({idx1, idx2, idx3, idx4}).real());
                        double comIm = compareVals(t.element({idx4, idx1, idx2, idx2, idx3}).imag(),
                                                   amplEval.element({idx1, idx2, idx3, idx4}).imag());
                        EXPECT_NEAR(comRe, 1., 1e-4);
                        EXPECT_NEAR(comIm, 1., 1e-4);
                    }
                }
            }
        }
    }


} // namespace Hammer
