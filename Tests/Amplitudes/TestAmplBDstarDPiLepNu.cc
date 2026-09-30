///
/// @file  TestAmplBDstarDPiLepNu.cc
/// @brief Tests for AmplBDstarDPiLepNu
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "Hammer/Amplitudes/AmplBDstarDPiLepNu.hh"
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

    TEST(AmplBDstarDPiLepNuTest, eval) {

        // Momentum decs
        FourMomentum pBmes{5280., 0., 0., 0.};
        FourMomentum pDstarmes{2258.9451, -542.95362175877925474, -718.43449427445361377, 513.65086086722812801};
        FourMomentum kNuTau{763.45512706233784082, 59.818100429037400761, 266.26625745054465575, 713.01318784585764875};
        FourMomentum pTau{2257.5997729376621592, 483.13552132974185397, 452.16823682390895801, -1226.6640487130857768};

        FourMomentum pDmes{2092.5602758805729893, -521.50730436566028339, -660.56418912084064663,
                           438.48442587823866432};
        FourMomentum pPion{166.38482411942701071, -21.446317393118971344, -57.870305153612967138,
                           75.166434988989463689};


        // Create FF and kin factors
        double Mb = pBmes.mass();
        double Mds = pDstarmes.mass();
        double Mds2 = pDstarmes.mass2();
        double Md2 = pDmes.mass2();
        // double Mt = pTau.mass();
        double Mp2 = pPion.mass2();
        double Sqq = Mb * Mb + Mds * Mds - 2. * (pBmes * pDstarmes);

        double Ed = (Mds2 + Md2 - Mp2) / (2. * Mds);
        double Pp = sqrt(Ed * Ed - Md2);
        double gPiOnSqrt2MGam = sqrt(3. * pi * Mds / pow(Pp, 3.));

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
        AmplBDstarDPiLepNu ampl;
        ampl.eval(Particle(pBmes, 511),
                  {Particle(pDstarmes, -423), Particle(kNuTau, 16), Particle(pTau, -15), Particle(pDmes, -421),
                   Particle(pPion, 111)},
                  {});
        auto t = ampl.getTensor();
        t.dot(FFvec, {FF_BDSTAR});

        // Compare to direct evaluation
        Tensor amplEval{"amplEval", MD::makeVector({2, 2, 11}, {SPIN_NUTAU, SPIN_TAUP, WILSON_BCTAUNU},
                                                   {(218526.018848 - 38094.3734443 * 1i),
                                                    (-41834.873825 + 0. * 1i),
                                                    0. * 1i,
                                                    (41834.873825 + 0. * 1i),
                                                    0. * 1i,
                                                    (-218526.018848 - 38094.3734443 * 1i),
                                                    0. * 1i,
                                                    (218526.018848 - 38094.3734443 * 1i),
                                                    0. * 1i,
                                                    (1994452.81665 + 625211.034433 * 1i),
                                                    0. * 1i,
                                                    (130542.449024 - 322278.778667 * 1i),
                                                    0. * 1i,
                                                    0. * 1i,
                                                    0. * 1i,
                                                    0. * 1i,
                                                    (69836.9021198 + 253350.864021 * 1i),
                                                    0. * 1i,
                                                    (130542.449024 - 322278.778667 * 1i),
                                                    0. * 1i,
                                                    (-432564.098293 - 609539.450478 * 1i),
                                                    0. * 1i,
                                                    0. * 1i,
                                                    0. * 1i,
                                                    0. * 1i,
                                                    0. * 1i,
                                                    0. * 1i,
                                                    0. * 1i,
                                                    (-130542.449024 - 322278.778667 * 1i),
                                                    0. * 1i,
                                                    (-69836.9021198 + 253350.864021 * 1i),
                                                    0. * 1i,
                                                    (432564.098293 - 609539.450478 * 1i),
                                                    0. * 1i,
                                                    0. * 1i,
                                                    (41834.873825 + 0. * 1i),
                                                    0. * 1i,
                                                    (-41834.873825 + 0. * 1i),
                                                    0. * 1i,
                                                    (218526.018848 + 38094.3734443 * 1i),
                                                    0. * 1i,
                                                    (-218526.018848 + 38094.3734443 * 1i),
                                                    0. * 1i,
                                                    (1994452.81665 - 625211.034433 * 1i)})};
        amplEval *= GFermi * gPiOnSqrt2MGam;


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

    TEST(AmplBDstarDPiLepNuTest, eval2) {

        // Momentum decs
        int sgn = -1;

        Particle pB({5.280000000000000, 0., 0., 0.}, 511 * sgn);
        Particle pDs({2.615097065151515, 1.040032785070634, 0.2531480829708513, 1.285644007581599}, -413 * sgn);
        Particle pTau({2.520191631332623, -1.140513773639112, -0.1839438715862363, -1.363462446451798}, -15 * sgn);
        Particle pNuTau({0.14472, 0.1004809885684785, -0.06920421138461502, 0.07781843887019962}, 16 * sgn);
        Particle pD({2.454091581515532, 0.9743983145367637, 0.2651629525963963, 1.228674498662554}, -411 * sgn);
        Particle pPiZ({0.1610054836359843, 0.06563447053387053, -0.01201486962554486, 0.05696950891904523}, 111);


        // Evaluate at FF point
        AmplBDstarDPiLepNu ampl;
        ampl.setSignatureIndex(1);
        ampl.eval(pB, {pDs, pNuTau, pTau, pD, pPiZ}, {});

        /// Check for nans

        size_t count = 0;
        const Tensor& t = ampl.getTensor();
        for (IndexType idx1 = 0; idx1 < 11; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 8; ++idx2) {
                for (IndexType idx3 = 0; idx3 < 2; ++idx3) {
                    for (IndexType idx4 = 0; idx4 < 2; ++idx4) {
                        EXPECT_FALSE(is_ieee754_nan(t.element({idx1, idx2, idx3, idx3, idx4}).real()));
                        EXPECT_FALSE(is_ieee754_nan(t.element({idx1, idx2, idx3, idx3, idx4}).imag()));
                        if (!(isZero(t.element({idx1, idx2, idx3, idx3, idx4}).real()) &&
                              isZero(t.element({idx1, idx2, idx3, idx3, idx4}).imag()))) {
                            ++count;
                        }
                    }
                }
            }
        }

        /// Check for non-zeros
        EXPECT_NE(count, 0);
    }

} // namespace Hammer
