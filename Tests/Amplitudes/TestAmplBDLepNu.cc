///
/// @file  TestAmplBDLepNu.cc
/// @brief Tests for AmplBDLepNu
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "Hammer/Amplitudes/AmplBDLepNu.hh"
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

    TEST(AmplBDLepNuTest, eval) {

        // Momentum decs
        FourMomentum pBmes{5280., 0., 0., 0.};
        FourMomentum pDmes{2381.3371929065643, 0., 0., -1475.671313782346};
        FourMomentum pTau{2034.777874233763, 553.6990600712329, -92.91177555543321, 819.1030156888944};
        FourMomentum kNuTau{863.8849328596731, -553.6990600712329, 92.91177555543321, 656.5682980934517};

        // Create FF
        double Mb = pBmes.mass();
        double Md = pDmes.mass();
        // double Mt = pTau.mass();
        double Sqq = Mb * Mb + Md * Md - 2. * (pBmes * pDmes);
        double Fp = 1.;
        double MbdSqq = pow(Mb + Md, 2) - Sqq;
        double Fz = Fp * MbdSqq / pow(Mb + Md, 2);
        double Fs = Fp * MbdSqq / (Mb + Md);
        double Ft = Fp / (Mb + Md);


        Tensor FFvec{"FFvec", MD::makeVector({4}, {FF_BD}, {Fs, Fz, Fp, Ft})};

        // Evaluate at FF point
        AmplBDLepNu ampl;
        ampl.eval(Particle(pBmes, 511), {Particle(pDmes, -411), Particle(kNuTau, 16), Particle(pTau, -15)}, {});
        auto t = ampl.getTensor();
        t.dot(FFvec, {FF_BD});

        // Compare to direct evaluation
        Tensor amplEval{"amplEval", MD::makeVector({2, 2, 11}, {SPIN_NUTAU, SPIN_TAUP, WILSON_BCTAUNU},
                                                   {-19651074.1655,
                                                    15558188.4079,
                                                    0.0,
                                                    15558188.4079,
                                                    0.0,
                                                    -19651074.1655,
                                                    0.0,
                                                    -19651074.1655,
                                                    0.0,
                                                    8843408.26834,
                                                    0.0,
                                                    -14122127.6492,
                                                    0.0,
                                                    0.0,
                                                    0.0,
                                                    0.0,
                                                    -14122127.6492,
                                                    0.0,
                                                    -14122127.6492,
                                                    0.0,
                                                    14033234.6929,
                                                    0.0,
                                                    0.0,
                                                    0.0,
                                                    0.0,
                                                    0.0,
                                                    0.0,
                                                    0.0,
                                                    -14122127.6492,
                                                    0.0,
                                                    -14122127.6492,
                                                    0.0,
                                                    14033234.6929,
                                                    0.0,
                                                    0.0,
                                                    -15558188.4079,
                                                    0.0,
                                                    -15558188.4079,
                                                    0.0,
                                                    19651074.1655,
                                                    0.0,
                                                    19651074.1655,
                                                    0.0,
                                                    -8843408.26834})};
        amplEval *= GFermi;


        // Check TEs. Care that nuTau spin is double-copied; second entry is ref qu.no.
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
