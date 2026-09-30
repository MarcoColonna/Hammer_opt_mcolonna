///
/// @file  TestAmplBcJpsiEllEllLepNu.cc
/// @brief Tests for AmplBcJpsiEllEllLepNu
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "Hammer/Amplitudes/AmplBcJpsiEllEllLepNu.hh"
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

    TEST(AmplBcJpsiEllEllLepNuTest, eval) {

        // Momentum decs
        FourMomentum pBmes{6275., 0., 0., 0.};
        FourMomentum pJPsimes{3821.58095617530, -857.35310207508, 0, -2069.83348677229};
        FourMomentum kNuEle{131.0024243058, -0.8276183320, 98.867858366606, -85.9424039919};
        FourMomentum pPos{2322.4166195189, 858.1807204071, -98.867858366606, 2155.7758907642};
        FourMomentum pMup{942.9153383107, 559.8448405981, 668.75990646868, 342.6343199059};
        FourMomentum pMum{2878.6656178646, -1417.1979426731, -668.75990646868, -2412.4678066782};

        Tensor FFvec{"FFvec",
                     MD::makeVector({8}, {FF_BCJPSI}, {-0.7, 5700, 0.0002, 0.0003, -0.0001, -6.e-9, 0.6, -1.1})};

        // Evaluate at FF point
        AmplBcJpsiEllEllLepNu ampl;
        ampl.eval(Particle(pBmes, 541),
                  {Particle(pJPsimes, 443), Particle(kNuEle, 12), Particle(pPos, -11), Particle(pMum, 13),
                   Particle(pMup, -13)},
                  {});
        auto t = ampl.getTensor();
        t.dot(FFvec, {FF_BCJPSI});

        // Compare to direct evaluation
        Tensor amplEval{"amplEval", MD::makeVector({2, 2, 2, 11}, {SPIN_MUM_REF, SPIN_NUE, SPIN_EP, WILSON_BCENU},
                                                   {-1245.929911257367 - 159.207667848831 * 1i,
                                                    -1.38000730704452e7,
                                                    0,
                                                    1.38000730704452e7,
                                                    0,
                                                    -8902.72499437946 + 2027.36579250475 * 1i,
                                                    0,
                                                    -1245.929911257367 - 159.207667848831 * 1i,
                                                    0,
                                                    5.83457215014948e8 - 8.8514177631392e7 * 1i,
                                                    0,
                                                    -2.5156109189465e7 - 2.7060272532613e6 * 1i,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    1.44458427859354e7 - 398452.0633446 * 1i,
                                                    0,
                                                    -2.5156109189465e7 - 2.7060272532613e6 * 1i,
                                                    0,
                                                    -43459.0141200144 + 4140.3554665915 * 1i,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    -1.65105138668831e8 + 3.7729200205593e7 * 1i,
                                                    0,
                                                    -1.105609175646636e7 - 225531.21179332 * 1i,
                                                    0,
                                                    -23100.9199515739 + 8111.9560730723 * 1i,
                                                    0,
                                                    0,
                                                    1.38000730704452e7,
                                                    0,
                                                    -1.38000730704452e7,
                                                    0,
                                                    8902.72499437946 - 2027.36579250475 * 1i,
                                                    0,
                                                    1245.929911257366 + 159.207667848831 * 1i,
                                                    0,
                                                    1.46528823774224e8 - 8.08466932918e6 * 1i,
                                                    8902.72499437946 + 2027.36579250475 * 1i,
                                                    -1.38000730704452e7,
                                                    0,
                                                    1.38000730704452e7,
                                                    0,
                                                    1245.929911257366 - 159.207667848831 * 1i,
                                                    0,
                                                    8902.72499437946 + 2027.36579250475 * 1i,
                                                    0,
                                                    1.46528823774224e8 + 8.08466932918e6 * 1i,
                                                    0,
                                                    1.65105138668831e8 + 3.7729200205593e7 * 1i,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    1.105609175646636e7 - 225531.21179332 * 1i,
                                                    0,
                                                    1.65105138668831e8 + 3.7729200205593e7 * 1i,
                                                    0,
                                                    23100.9199515739 + 8111.9560730723 * 1i,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    2.5156109189465e7 - 2.7060272532613e6 * 1i,
                                                    0,
                                                    -1.44458427859354e7 - 398452.0633446 * 1i,
                                                    0,
                                                    43459.0141200144 + 4140.3554665915 * 1i,
                                                    0,
                                                    0,
                                                    1.38000730704452e7,
                                                    0,
                                                    -1.38000730704452e7,
                                                    0,
                                                    -1245.929911257367 + 159.207667848831 * 1i,
                                                    0,
                                                    -8902.72499437946 - 2027.36579250475 * 1i,
                                                    0,
                                                    5.83457215014949e8 + 8.8514177631392e7 * 1i})};
        amplEval *= GFermi;


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
