///
/// @file  TestFFBLPRBase.cc
/// @brief Tests for FFBLPRBase
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"
#include "Hammer/FormFactors/BLPR/FFBLPRBase.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/Utils.hh"
// May need other stuff here

using namespace std;
using namespace complex_literals;

namespace Hammer {

    class TestableFFBLPRBase : public FFBLPRBase {

    public:

        TestableFFBLPRBase() = default;

        using FFBLPRBase::defineSettings;
        using FFBLPRBase::Xi;
        using FFBLPRBase::ZofW;

        std::function<double(double)>& getIW(size_t id) {
            return _mIWs[id];
        }
        std::function<double(double)>& getLi(size_t id) {
            return _mLi[id];
        }

        void evalAtPSPoint(const vector<double>& /*point*/, const vector<double>& /*masses*/) override {
        }
        unique_ptr<FormFactorBase> clone(const string& /*label*/) override {
            return unique_ptr<FormFactorBase>{};
        }
    };

    TEST(FFBLPRBase, Z) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        // Evaluate the FF class
        TestableFFBLPRBase ff;
        ff.setSettingsHandler(set);
        ff.defineSettings();
        ff.Xi(1.);
        EXPECT_NEAR(ff.ZofW(1.), (1 - 1.509 / sqrt2) / (1 + 1.509 / sqrt2), 0.001);
        EXPECT_NEAR(ff.ZofW(1.5), (sqrt(2.5) - 1.509) / (sqrt(2.5) + 1.509), 0.001);
    }

    TEST(FFBLPRBase, Xi) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        double w = 1.274123699;
        double xiVal = 0.69672215903902468;


        // Evaluate the FF class
        TestableFFBLPRBase ff;
        ff.setSettingsHandler(set);
        ff.defineSettings();

        auto xi = ff.Xi(w);
        EXPECT_NEAR(xi.first / xi.second, xiVal, 0.001);
    }

    TEST(FFBLPRBase, IWs) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        // Evaluate the FF class
        TestableFFBLPRBase ff;
        ff.setSettingsHandler(set);
        ff.defineSettings();
        // chi2
        EXPECT_DOUBLE_EQ(ff.getIW(1)(1.), -0.06);
        EXPECT_DOUBLE_EQ(ff.getIW(1)(1.5), -0.06 + 0.5 * (0.0));
        // chi3
        EXPECT_DOUBLE_EQ(ff.getIW(2)(1.), 0.);
        EXPECT_DOUBLE_EQ(ff.getIW(2)(1.5), 0. + 0.5 * (0.05));
        // eta
        EXPECT_DOUBLE_EQ(ff.getIW(0)(1.), 0.3);
        EXPECT_DOUBLE_EQ(ff.getIW(0)(1.5), 0.3 + 0.5 * (-0.05));
    }

    TEST(FFBLPRBase, Lis) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        // Evaluate the FF class
        TestableFFBLPRBase ff;
        ff.setSettingsHandler(set);
        ff.defineSettings();
        EXPECT_DOUBLE_EQ(ff.getLi(1)(1.), 0.);
        EXPECT_DOUBLE_EQ(ff.getLi(1)(1.5), 4. * (3. * ff.getIW(2)(1.5) - 0.5 * ff.getIW(1)(1.5)));
        EXPECT_DOUBLE_EQ(ff.getLi(2)(1.), 0.);
        EXPECT_DOUBLE_EQ(ff.getLi(2)(1.5), -4. * ff.getIW(2)(1.5));
        EXPECT_DOUBLE_EQ(ff.getLi(3)(1.), 4. * ff.getIW(1)(1.));
        EXPECT_DOUBLE_EQ(ff.getLi(3)(1.5), 4. * ff.getIW(1)(1.5));
        EXPECT_DOUBLE_EQ(ff.getLi(4)(1.), 2. * ff.getIW(0)(1.) - 1.);
        EXPECT_DOUBLE_EQ(ff.getLi(4)(1.5), 2. * ff.getIW(0)(1.5) - 1.);
        EXPECT_DOUBLE_EQ(ff.getLi(5)(1.), -1.);
        EXPECT_DOUBLE_EQ(ff.getLi(5)(1.5), -1.);
        EXPECT_DOUBLE_EQ(ff.getLi(6)(1.), -(ff.getIW(0)(1.) + 1.));
        EXPECT_DOUBLE_EQ(ff.getLi(6)(1.5), -2. * (ff.getIW(0)(1.5) + 1.) / 2.5);
    }

} // namespace Hammer
