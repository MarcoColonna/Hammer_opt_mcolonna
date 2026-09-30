///
/// @file  TestAmplTauEllNuNu.cc
/// @brief Tests for AmplTauEllNuNu
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "Hammer/Amplitudes/AmplBDLepNu.hh"
#include "Hammer/Amplitudes/AmplTauEllNuNu.hh"
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

    TEST(AmplTauEllNuNuTest, evalTauDec) {

        // Momentum decs
        FourMomentum pBmes{5280., 0., 0., 0.};
        FourMomentum pDmes{2381.3371929065643, 0., 0., -1475.671313782346};
        FourMomentum kNuTau{863.8849328596731, -553.6990600712329, 92.91177555543321, 656.5682980934517};

        FourMomentum pTau{2034.777874233763, 553.6990600712329, -92.91177555543321, 819.1030156888944};
        FourMomentum kNuBarTau{751.3592069889307, -309.29260061671044, -268.9081311284792, 629.7357875679096};
        FourMomentum kNuMuon{532.7499030827818, 178.97752649767912, 105.73956803205624, 490.5187539695908};
        FourMomentum kMuon{750.6687641620506, 684.0141341902644, 70.25678754098979, -301.15152584860607};

        // Evaluate at PS point
        AmplTauEllNuNu ampl;
        ampl.eval(Particle(pTau, -15), {Particle(kNuBarTau, -16), Particle(kNuMuon, 14), Particle(kMuon, -13)},
                  {Particle(pDmes, -411), Particle(kNuTau, 16)});

        // Compare to direct evaluation
        Tensor amplEval{"amplEval", MD::makeVector({2, 2}, {SPIN_NUTAU_REF, SPIN_TAUP},
                                                   {-5091457.7666213075 + 1873392.2733436273 * 1i,
                                                    4804390.213406682 - 4709550.054420433 * 1i,
                                                    -563734.9263256947 + 5395807.977676054 * 1i,
                                                    -2108373.2013852084 - 6388801.8815306835 * 1i})};
        amplEval *= GFermi;

        // Check TEs.
        for (IndexType idx1 = 0; idx1 < 2; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 2; ++idx2) {
                const Tensor& t = ampl.getTensor();
                EXPECT_NEAR(t.element({idx1, idx2}).real(), amplEval.element({idx1, idx2}).real(), 1e-3);
                EXPECT_NEAR(t.element({idx1, idx2}).imag(), amplEval.element({idx1, idx2}).imag(), 1e-3);
            }
        }
    }

    TEST(AmplTauEllNuNuTest, evalTauDecNoRef) {

        // Momentum decs
        FourMomentum pDmes{1., 0., 0., 0.};
        FourMomentum kNuTau{1., 0., 0., 1.};
        FourMomentum pTau{2034.777874233763, 553.6990600712329, -92.91177555543321, 819.1030156888944};
        // FourMomentum pBmes = pDmes + pTau + pDmes;

        FourMomentum kNuBarTau{751.3592069889307, -309.29260061671044, -268.9081311284792, 629.7357875679096};
        FourMomentum kNuMuon{532.7499030827818, 178.97752649767912, 105.73956803205624, 490.5187539695908};
        FourMomentum kMuon{750.6687641620506, 684.0141341902644, 70.25678754098979, -301.15152584860607};

        // Evaluate at PS point
        AmplTauEllNuNu ampl;
        ampl.eval(Particle(pTau, -15), {Particle(kNuBarTau, -16), Particle(kNuMuon, 14), Particle(kMuon, -13)}, {});

        // Compare to direct evaluation
        Tensor amplEval{"amplEval", MD::makeVector({2, 2}, {SPIN_NUTAU_REF, SPIN_TAUP},
                                                   {-2626867.5979955015 + 29144.017772737312 * 1i,
                                                    5876614.5917958422 - 5767034.4207547056 * 1i,
                                                    1878989.0864361296 + 1835941.92875344 * 1i,
                                                    -8233340.3715081867 - 73426.239208280575 * 1i})};
        amplEval *= GFermi;

        // Check TEs.
        for (IndexType idx1 = 0; idx1 < 2; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 2; ++idx2) {
                const Tensor& t = ampl.getTensor();
                EXPECT_NEAR(t.element({idx1, idx2}).real(), amplEval.element({idx1, idx2}).real(), 1e-3);
                EXPECT_NEAR(t.element({idx1, idx2}).imag(), amplEval.element({idx1, idx2}).imag(), 1e-3);
            }
        }
    }

    TEST(AmplTauEllNuNuTest, evalAndContract) {
        // Momentum decs
        FourMomentum pBmes{5280., 0., 0., 0.};
        FourMomentum pDmes{2381.3371929065643, 0., 0., -1475.671313782346};
        FourMomentum kNuTau{863.8849328596731, -553.6990600712329, 92.91177555543321, 656.5682980934517};

        FourMomentum pTau{2034.777874233763, 553.6990600712329, -92.91177555543321, 819.1030156888944};
        FourMomentum kNuBarTau{751.3592069889307, -309.29260061671044, -268.9081311284792, 629.7357875679096};
        FourMomentum kNuMuon{532.7499030827818, 178.97752649767912, 105.73956803205624, 490.5187539695908};
        FourMomentum kMuon{750.6687641620506, 684.0141341902644, 70.25678754098979, -301.15152584860607};

        // Evaluate at PS point
        AmplTauEllNuNu ampl;
        ampl.eval(Particle(pTau, -15), {Particle(kNuBarTau, -16), Particle(kNuMuon, 14), Particle(kMuon, -13)},
                  {Particle(pDmes, -411), Particle(kNuTau, 16)});

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
        AmplBDLepNu amplParent;
        amplParent.eval(Particle(pBmes, -511), {Particle(pDmes, -411), Particle(kNuTau, 16), Particle(pTau, -15)}, {});
        auto t = amplParent.getTensor();
        t.dot(FFvec, {FF_BD});
        t.dot(ampl.getTensor(), {SPIN_TAUP, SPIN_NUTAU_REF});

        // Compare to direct evaluation
        Tensor amplParentEval{"amplParentEval", MD::makeVector({2, 11}, {SPIN_NUTAU, WILSON_BCTAUNU},
                                                               {(3.2204402312e+13 + 2.96946965342e+13 * 1i),
                                                                (-7.92138592042e+13 + 2.91465899507e+13 * 1i),
                                                                0. * 1i,
                                                                (-7.92138592042e+13 + 2.91465899507e+13 * 1i),
                                                                0. * 1i,
                                                                (3.2204402312e+13 + 2.96946965342e+13 * 1i),
                                                                0. * 1i,
                                                                (3.2204402312e+13 + 2.96946965342e+13 * 1i),
                                                                0. * 1i,
                                                                (2.23952957099e+13 - 4.95230484918e+13 * 1i),
                                                                0. * 1i,
                                                                0. * 1i,
                                                                0. * 1i,
                                                                (3.28024675014e+13 + 9.93981833738e+13 * 1i),
                                                                0. * 1i,
                                                                (3.28024675014e+13 + 9.93981833738e+13 * 1i),
                                                                0. * 1i,
                                                                (-3.3470661559e+13 - 2.01747108634e+14 * 1i),
                                                                0. * 1i,
                                                                (-3.3470661559e+13 - 2.01747108634e+14 * 1i),
                                                                0. * 1i,
                                                                (1.07341804761e+13 + 1.32219423093e+14 * 1i)})};
        amplParentEval *= pow(GFermi, 2.);


        // Check TEs.
        for (IndexType idx1 = 0; idx1 < 2; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 11; ++idx2) {
                double comRe = compareVals(t.element({idx2, idx1}).real(), amplParentEval.element({idx1, idx2}).real());
                double comIm = compareVals(t.element({idx2, idx1}).imag(), amplParentEval.element({idx1, idx2}).imag());
                EXPECT_NEAR(comRe, 1., 1e-4);
                EXPECT_NEAR(comIm, 1., 1e-4);
            }
        }
    }
} // namespace Hammer
