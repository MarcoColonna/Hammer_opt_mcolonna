///
/// @file  TestAmplBDstarLepNu.cc
/// @brief Tests for AmplBDstarLepNu
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "Hammer/Amplitudes/AmplBDstarLepNu.hh"
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

    TEST(AmplBDstarLepNuTest, eval) {

        // Momentum decs
        FourMomentum pBmes{5280., 0., 0., 0.};
        FourMomentum pDstarmes{2196.3156772683462121, 467.82967757107907740, -747.52221490115041077,
                               -77.772007996018766560};
        FourMomentum kNuEle{1866.4064989883714089, -78.510370590731764429, 1843.3032585533182967,
                            -282.03268991812669586};
        FourMomentum pPos{1217.2778237432823790, -389.31930698034731297, -1095.7810436521678859, 359.80469791414546242};


        // Create FF and kin factors
        double Mb = pBmes.mass();
        double Mds = pDstarmes.mass();
        // double Mds2 = pDstarmes.mass2();
        double Sqq = Mb * Mb + Mds * Mds - 2. * (pBmes * pDstarmes);


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
        AmplBDstarLepNu ampl;
        ampl.eval(Particle(pBmes, 511), {Particle(pDstarmes, -413), Particle(kNuEle, 12), Particle(pPos, -11)}, {});
        auto t = ampl.getTensor();
        t.dot(FFvec, {FF_BDSTAR});

        // Compare to direct evaluation
        Tensor amplEval{"amplEval", MD::makeVector({3, 2, 2, 11}, {SPIN_DSTAR, SPIN_NUE, SPIN_EP, WILSON_BCTAUNU},
                                                   {-0.8411650095,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    0.5486596758,
                                                    0.,
                                                    -0.8411650095,
                                                    0.,
                                                    -9441.733224,
                                                    0.,
                                                    -1907.522720,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    1244.203914,
                                                    0.,
                                                    -1907.522720,
                                                    0.,
                                                    -0.6407614565,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    -8084.597221,
                                                    0.,
                                                    12394.71497,
                                                    0.,
                                                    -11.52412191,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    -0.5486596758,
                                                    0.,
                                                    0.8411650095,
                                                    0.,
                                                    -26133.42702,
                                                    -1.923924495,
                                                    -1758.375357,
                                                    0.,
                                                    1758.375357,
                                                    0.,
                                                    1.923924495,
                                                    0.,
                                                    -1.923924495,
                                                    0.,
                                                    -24504.90650,
                                                    0.,
                                                    -6288.859262,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    6288.859262,
                                                    0.,
                                                    -6288.859262,
                                                    0.,
                                                    -3.931017839,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    6288.859262,
                                                    0.,
                                                    -6288.859262,
                                                    0.,
                                                    3.931017839,
                                                    0.,
                                                    0.,
                                                    1758.375357,
                                                    0.,
                                                    -1758.375357,
                                                    0.,
                                                    -1.923924495,
                                                    0.,
                                                    1.923924495,
                                                    0.,
                                                    -24504.90650,
                                                    -0.5486596758,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    0.8411650095,
                                                    0.,
                                                    -0.5486596758,
                                                    0.,
                                                    -26133.42702,
                                                    0.,
                                                    8084.597221,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    -12394.71497,
                                                    0.,
                                                    8084.597221,
                                                    0.,
                                                    11.52412191,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    1907.522720,
                                                    0.,
                                                    -1244.203914,
                                                    0.,
                                                    0.6407614565,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    -0.8411650095,
                                                    0.,
                                                    0.5486596758,
                                                    0.,
                                                    -9441.733224})};
        amplEval *= GFermi;


        // Check TEs. Care that nuTau spin (second index) is double-copied; third index is ref qu.no.
        for (IndexType idx1 = 0; idx1 < 2; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 2; ++idx2) {
                for (IndexType idx3 = 0; idx3 < 11; ++idx3) {
                    for (IndexType idx4 = 0; idx4 < 3; ++idx4) {
                        double comRe = compareVals(t.element({idx3, idx4, idx1, idx1, idx2}).real(),
                                                   amplEval.element({idx4, idx1, idx2, idx3}).real());
                        double comIm = compareVals(t.element({idx3, idx4, idx1, idx1, idx2}).imag(),
                                                   amplEval.element({idx4, idx1, idx2, idx3}).imag());
                        EXPECT_NEAR(comRe, 1., 1e-4);
                        EXPECT_NEAR(comIm, 1., 1e-4);
                    }
                }
            }
        }
    }

} // namespace Hammer
