///
/// @file  TestFFTauto3PiRCT.cc
/// @brief Tests for FFTauto3PiRCT
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#include <filesystem>
#include <fstream>

#define FORMFACTORBASE_FRIENDS              \
    FRIEND_TEST(FFTauto3PiRCTTest, eval);   \
    FRIEND_TEST(FFTauto3PiRCTTest, clone);  \
    FRIEND_TEST(FFTauto3PiRCTTest, addRefs)

#include "Hammer/FormFactors/RCT/FFTauto3PiRCT.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/FourMomentum.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Math/MultiDim/VectorContainer.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Tools/SettingsConsumer.hh"

using namespace std;
using namespace complex_literals;

namespace Hammer {

    namespace MD = MultiDimensional;

    TEST(FFTauto3PiRCTTest, eval) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs
        FourMomentum pBmes{9.00070555033828, 2.2348555435636683, -2.8664110899642328, -6.318497334542691};
        FourMomentum pDmes{5.29977696322721, 1.1366307865620362, -2.1310211510226202, -4.331430960882029};
        FourMomentum kNuTau{1.706348474218978, 1.0982247570016317, -0.7353899389416136, -1.0792261747574472};

        FourMomentum pTau{1.9945801128920924, 0., 0., -0.9078401989032162};
        FourMomentum kNuBarTau{0.9078401989032162, 0., 0., -0.9078401989032162};

        FourMomentum pPiPlus1{0.43768637864473114, -0.10734791320010291, -0.3089042414424765, 0.25554444075653515};
        FourMomentum pPiPlus2{0.483312990808882, 0.041088807936159696, 0.34705949502947275, -0.3035323110591347};
        FourMomentum pPiMinus{0.16574054453526274, 0.0662591052639432, -0.03815525358699622, 0.04798787030259953};

        // Evaluate the FF class
        FFTauto3PiRCT ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();
        ff.eval(Particle(pTau, -15),
                {Particle(pPiPlus1 *= 1e3, 211), Particle(pPiPlus2 *= 1e3, 211), Particle(pPiMinus *= 1e3, -211),
                 Particle(kNuBarTau *= 1e3, -16)},
                {Particle(pDmes *= 1e3, -411), Particle(kNuTau *= 1e3, 16)});

        // Compare to direct evaluation
        // F1, F2, F4(q^2, s1, s2)
        Tensor ffEval{"ffEval", MD::makeVector({3}, {FF_TAU3PI},
                                               {-14.34933926 - 34.29602840 * 1i, -3.85342367 - 40.89249792 * 1i,
                                                -0.01846726652 + 0.00480450641 * 1i})};
        ffEval *= 1e-3; // GeV^-1 to MeV^-1

        auto& t = ff.getTensor();
        for (IndexType idx1 = 0; idx1 < 3; ++idx1) {
            cout << t.element({idx1}).real() << " " << t.element({idx1}).imag() << '\n';
            double comRe = compareVals(t.element({idx1}).real(), ffEval.element({idx1}).real());
            double comIm = compareVals(t.element({idx1}).imag(), ffEval.element({idx1}).imag());
            //            cout << idx1 << "  " << setprecision (15) << comRe << "  " << comIm << "\n";
            EXPECT_NEAR(comRe, 1., 1e-4);
            EXPECT_NEAR(comIm, 1., 1e-4);
        }
    }

    TEST(FFTauto3PiRCTTest, clone) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFTauto3PiRCT ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto cloned = ff.clone("testLabel");
        ASSERT_NE(cloned, nullptr);
        EXPECT_NE(cloned->group().find("testLabel"), string::npos);

        FourMomentum pBmes{9.00070555033828, 2.2348555435636683, -2.8664110899642328, -6.318497334542691};
        FourMomentum pDmes{5.29977696322721, 1.1366307865620362, -2.1310211510226202, -4.331430960882029};
        FourMomentum kNuTau{1.706348474218978, 1.0982247570016317, -0.7353899389416136, -1.0792261747574472};

        FourMomentum pTau{1.9945801128920924, 0., 0., -0.9078401989032162};
        FourMomentum kNuBarTau{0.9078401989032162, 0., 0., -0.9078401989032162};

        FourMomentum pPiPlus1{0.43768637864473114, -0.10734791320010291, -0.3089042414424765, 0.25554444075653515};
        FourMomentum pPiPlus2{0.483312990808882, 0.041088807936159696, 0.34705949502947275, -0.3035323110591347};
        FourMomentum pPiMinus{0.16574054453526274, 0.0662591052639432, -0.03815525358699622, 0.04798787030259953};

        cloned->setSettingsHandler(set);
        cloned->initSettings();
        cloned->calcUnits();
        cloned->eval(Particle(pTau, -15),
                     {Particle(pPiPlus1 *= 1e3, 211), Particle(pPiPlus2 *= 1e3, 211), Particle(pPiMinus *= 1e3, -211),
                      Particle(kNuBarTau *= 1e3, -16)},
                     {Particle(pDmes *= 1e3, -411), Particle(kNuTau *= 1e3, 16)});

        // Reset momenta for the original ff call (already scaled above for cloned)
        FourMomentum pBmes2{9.00070555033828, 2.2348555435636683, -2.8664110899642328, -6.318497334542691};
        FourMomentum pDmes2{5.29977696322721, 1.1366307865620362, -2.1310211510226202, -4.331430960882029};
        FourMomentum kNuTau2{1.706348474218978, 1.0982247570016317, -0.7353899389416136, -1.0792261747574472};

        FourMomentum pTau2{1.9945801128920924, 0., 0., -0.9078401989032162};
        FourMomentum kNuBarTau2{0.9078401989032162, 0., 0., -0.9078401989032162};

        FourMomentum pPiPlus1b{0.43768637864473114, -0.10734791320010291, -0.3089042414424765, 0.25554444075653515};
        FourMomentum pPiPlus2b{0.483312990808882, 0.041088807936159696, 0.34705949502947275, -0.3035323110591347};
        FourMomentum pPiMinusb{0.16574054453526274, 0.0662591052639432, -0.03815525358699622, 0.04798787030259953};

        ff.eval(Particle(pTau2, -15),
                {Particle(pPiPlus1b *= 1e3, 211), Particle(pPiPlus2b *= 1e3, 211), Particle(pPiMinusb *= 1e3, -211),
                 Particle(kNuBarTau2 *= 1e3, -16)},
                {Particle(pDmes2 *= 1e3, -411), Particle(kNuTau2 *= 1e3, 16)});

        auto& tOrig = ff.getTensor();
        auto& tCloned = cloned->getTensor();
        for (IndexType idx = 0; idx < 3; ++idx) {
            EXPECT_NEAR(tOrig.element({idx}).real(), tCloned.element({idx}).real(), 1e-12);
            EXPECT_NEAR(tOrig.element({idx}).imag(), tCloned.element({idx}).imag(), 1e-12);
        }
    }

    TEST(FFTauto3PiRCTTest, addRefs) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFTauto3PiRCT ff;
        ff.setSettingsHandler(set);
        ff.initSettings();

        EXPECT_FALSE(set.checkReference("Shekhovtsova:2012ra"));
        EXPECT_FALSE(set.checkReference("Nugent:2013hxa"));

        static_cast<SettingsConsumer*>(&ff)->addRefs();
        EXPECT_TRUE(set.checkReference("Shekhovtsova:2012ra"));
        EXPECT_TRUE(set.checkReference("Nugent:2013hxa"));

        // Second call must be idempotent.
        static_cast<SettingsConsumer*>(&ff)->addRefs();
        EXPECT_TRUE(set.checkReference("Shekhovtsova:2012ra"));
        EXPECT_TRUE(set.checkReference("Nugent:2013hxa"));

        const string tmpfile = testing::TempDir() + "hammer_rct_refs_test.bib";
        set.saveReferences(tmpfile);

        ifstream reffile(tmpfile);
        ASSERT_TRUE(reffile.is_open());
        string content((istreambuf_iterator<char>(reffile)), istreambuf_iterator<char>());
        reffile.close();
        std::filesystem::remove(tmpfile);

        EXPECT_NE(content.find("Shekhovtsova:2012ra"), string::npos);
        EXPECT_NE(content.find("Nugent:2013hxa"), string::npos);
    }


} // namespace Hammer
