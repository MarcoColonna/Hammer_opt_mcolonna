///
/// @file  TestFFBLPRXPBase.cc
/// @brief Tests for FFBLPRXPBase
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"
#include "Hammer/FormFactors/BLPRXP/FFBLPRXPBase.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/Utils.hh"
// May need other stuff here

using namespace std;
using namespace complex_literals;

namespace Hammer {

    class TestableFFBLPRXPBase : public FFBLPRXPBase {

    public:

        TestableFFBLPRXPBase() = default;

        using FFBLPRXPBase::calcConstants;
        using FFBLPRXPBase::defineSettings;
        using FFBLPRXPBase::Xi;
        using FFBLPRXPBase::ZofW;

        std::function<double(double)>& getIW(size_t id) {
            return _mIWs[id];
        }
        std::function<double(double)>& getLi1(size_t id) {
            return _mLi1[id];
        }
        std::function<double(double)>& getLi2(size_t id) {
            return _mLi2[id];
        }
        std::function<double(double)>& getMi(size_t id) {
            return _mMi[id];
        }

        std::array<double, 10> getProtectedVars() const {
            return {_a, _zBC, _la2OverlaB2, _la1OverlaB2, _eb, _ec, _upsilonc, _upsilonb, _cmagc, _cmagb};
        }

        void evalAtPSPoint(const vector<double>& /*point*/, const vector<double>& /*masses*/) override {
        }
        unique_ptr<FormFactorBase> clone(const string& /*label*/) override {
            return unique_ptr<FormFactorBase>{};
        }
    };

    TEST(FFBLPRXPBase, Z) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        // Evaluate the FF class
        TestableFFBLPRXPBase ff;
        ff.setSettingsHandler(set);
        ff.defineSettings();
        ff.calcConstants();
        EXPECT_NEAR(ff.ZofW(1.), (1 - 1.509 / sqrt2) / (1 + 1.509 / sqrt2), 0.001);
        EXPECT_NEAR(ff.ZofW(1.5), (sqrt(2.5) - 1.509) / (sqrt(2.5) + 1.509), 0.001);
    }

    TEST(FFBLPRXPBase, CalcConstants) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        // Evaluate the FF class
        TestableFFBLPRXPBase ff;
        ff.setSettingsHandler(set);
        ff.defineSettings();
        ff.calcConstants();
        double laB = (4.707 * 5.313 - (4.707 - 3.407) * 1.973) / 3.407 - (2. * 4.707 - 3.407) +
                     (-0.363) / (4. * 4.707 * (4.707 - 3.407));
        double la1 = (2. * 4.707 * (4.707 - 3.407)) / 3.407 * (5.313 - 1.973 - 3.407) +
                     (-0.363) * (2. * 4.707 - 3.407) / (2. * 4.707 * (4.707 - 3.407));
        double c1S = 2. * pow(0.27 / 3. * 0.796, 2.);
        auto res = ff.getProtectedVars();
        EXPECT_DOUBLE_EQ(res[0], 1.509 / sqrt2);
        EXPECT_DOUBLE_EQ(res[1], 1. - 3.407 / 4.707);
        EXPECT_DOUBLE_EQ(res[2], 0.12 / (laB * laB));
        EXPECT_DOUBLE_EQ(res[3], la1 / (laB * laB));
        EXPECT_DOUBLE_EQ(res[4], laB / 4.707 / 2.);
        EXPECT_DOUBLE_EQ(res[5], laB / (4.707 - 3.407) / 2.);
        EXPECT_DOUBLE_EQ(res[6] - res[7], c1S * (1. - 4.707 / (4.707 - 3.407)));
        EXPECT_DOUBLE_EQ(res[7] + c1S, (c1S * 4.707 / (3.407) * (5.313 - 1.973) - 2 * c1S * 4.707) / laB);
        EXPECT_DOUBLE_EQ(res[8], -3. / 2. * (0.5 * log(1. - 3.407 / 4.707) - 13. / 9.));
        EXPECT_DOUBLE_EQ(res[9], 3. / 2. * (0.5 * log(1. - 3.407 / 4.707) + 13. / 9.));
    }

    TEST(FFBLPRXPBase, Xi) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        // Evaluate the FF class
        TestableFFBLPRXPBase ff;
        ff.setSettingsHandler(set);
        ff.defineSettings();
        ff.calcConstants();
        double z1 = ff.ZofW(1.);
        double z1p25 = ff.ZofW(1.25);
        double a2 = 1.509 * 1.509 / 2.;
        double a4 = a2 * a2;
        double RhoStSq = 1.104;
        double cSt = 2.392;
        EXPECT_NEAR(ff.Xi(1.), 1., 0.001);
        EXPECT_NEAR(ff.Xi(1.25),
                    (1. - 8. * a2 * RhoStSq * z1p25 + 16 * (2 * a4 * cSt - 1.13 * a2) * z1p25 * z1p25) /
                        (1. - 8. * a2 * RhoStSq * z1 + 16 * (2 * a4 * cSt - RhoStSq * a2) * z1 * z1),
                    0.001);
    }

    TEST(FFBLPRXPBase, IWs) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        // Evaluate the FF class
        TestableFFBLPRXPBase ff;
        ff.setSettingsHandler(set);
        ff.defineSettings();
        set.changeSetting<double>("", "chi21", 0.82);
        set.changeSetting<double>("", "chi2p", -0.69);
        set.changeSetting<double>("", "chi3p", 0.3);
        set.changeSetting<double>("", "eta1", 0.34);
        set.changeSetting<double>("", "phi1p", 0.26);
        ff.calcConstants();
        auto cs = ff.getProtectedVars();
        // chi2
        EXPECT_DOUBLE_EQ(ff.getIW(0)(1.), 0.82);
        EXPECT_DOUBLE_EQ(ff.getIW(0)(1.5), 0.82 + 0.5 * (-0.69));
        // chi3
        EXPECT_DOUBLE_EQ(ff.getIW(1)(1.), 0.);
        EXPECT_DOUBLE_EQ(ff.getIW(1)(1.5), 0. + 0.5 * (0.3));
        // eta
        EXPECT_DOUBLE_EQ(ff.getIW(2)(1.), 0.34);
        EXPECT_DOUBLE_EQ(ff.getIW(2)(1.5), 0.34 + 0.5 * 0.);
        // beta1
        EXPECT_DOUBLE_EQ(ff.getIW(3)(1.), cs[3] / 4.);
        EXPECT_DOUBLE_EQ(ff.getIW(3)(1.5), cs[3] / 4.);
        // beta2
        EXPECT_DOUBLE_EQ(ff.getIW(4)(1.), 0.0);
        EXPECT_DOUBLE_EQ(ff.getIW(4)(1.5), 0.0);
        // beta3
        EXPECT_DOUBLE_EQ(ff.getIW(5)(1.), cs[2] / 8.);
        EXPECT_DOUBLE_EQ(ff.getIW(5)(1.5), cs[2] / 8.);
        // phi1
        EXPECT_DOUBLE_EQ(ff.getIW(6)(1.), (cs[3] / 3. - cs[2] / 2.) / 2.);
        EXPECT_DOUBLE_EQ(ff.getIW(6)(1.5), (cs[3] / 3. - cs[2] / 2.) / 2. + 0.5 * (0.26));
        // phiq
        EXPECT_DOUBLE_EQ(ff.getIW(7)(1.), 0.26);
        EXPECT_DOUBLE_EQ(ff.getIW(7)(1.5), 0.26);
    }

    TEST(FFBLPRXPBase, Li1s) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        // Evaluate the FF class
        TestableFFBLPRXPBase ff;
        ff.setSettingsHandler(set);
        ff.defineSettings();
        ff.calcConstants();
        // auto cs = ff.getProtectedVars();
        EXPECT_DOUBLE_EQ(ff.getLi1(1)(1.), 0.);
        EXPECT_DOUBLE_EQ(ff.getLi1(1)(1.5), 4. * (3. * ff.getIW(1)(1.5) - 0.5 * ff.getIW(0)(1.5)));
        EXPECT_DOUBLE_EQ(ff.getLi1(2)(1.), 0.);
        EXPECT_DOUBLE_EQ(ff.getLi1(2)(1.5), -4. * ff.getIW(1)(1.5));
        EXPECT_DOUBLE_EQ(ff.getLi1(3)(1.), 4. * ff.getIW(0)(1.));
        EXPECT_DOUBLE_EQ(ff.getLi1(3)(1.5), 4. * ff.getIW(0)(1.5));
        EXPECT_DOUBLE_EQ(ff.getLi1(4)(1.), 2. * ff.getIW(2)(1.) - 1.);
        EXPECT_DOUBLE_EQ(ff.getLi1(4)(1.5), 2. * ff.getIW(2)(1.5) - 1.);
        EXPECT_DOUBLE_EQ(ff.getLi1(5)(1.), -1.);
        EXPECT_DOUBLE_EQ(ff.getLi1(5)(1.5), -1.);
        EXPECT_DOUBLE_EQ(ff.getLi1(6)(1.), -(ff.getIW(2)(1.) + 1.));
        EXPECT_DOUBLE_EQ(ff.getLi1(6)(1.5), -2. * (ff.getIW(2)(1.) + 1.) / 2.5);
    }

    TEST(FFBLPRXPBase, Li2s) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        // Evaluate the FF class
        TestableFFBLPRXPBase ff;
        ff.setSettingsHandler(set);
        ff.defineSettings();
        ff.calcConstants();
        auto cs = ff.getProtectedVars();
        EXPECT_DOUBLE_EQ(ff.getLi2(1)(1.), 2. * ff.getIW(3)(1.) + 4. * (3. * ff.getIW(5)(1.)));
        EXPECT_DOUBLE_EQ(ff.getLi2(1)(1.5),
                         2. * ff.getIW(3)(1.5) + 4. * (3. * ff.getIW(5)(1.5) - 0.5 * ff.getIW(4)(1.5)));
        EXPECT_DOUBLE_EQ(ff.getLi2(2)(1.), 2. * ff.getIW(3)(1.) - 4. * (ff.getIW(5)(1.)));
        EXPECT_DOUBLE_EQ(ff.getLi2(2)(1.5), 2. * ff.getIW(3)(1.5) - 4. * ff.getIW(5)(1.5));
        EXPECT_DOUBLE_EQ(ff.getLi2(3)(1.), 4. * ff.getIW(4)(1.));
        EXPECT_DOUBLE_EQ(ff.getLi2(3)(1.5), 4. * ff.getIW(4)(1.5));
        EXPECT_DOUBLE_EQ(ff.getLi2(4)(1.), 3. * cs[2] + 2. * (2.) * ff.getIW(6)(1.));
        EXPECT_DOUBLE_EQ(ff.getLi2(4)(1.5), 3. * cs[2] + 2. * (2.5) * ff.getIW(6)(1.5));
        EXPECT_DOUBLE_EQ(ff.getLi2(5)(1.), cs[2] + 2. * (2.) * ff.getIW(6)(1.));
        EXPECT_DOUBLE_EQ(ff.getLi2(5)(1.5), cs[2] + 2. * (2.5) * ff.getIW(6)(1.5));
        EXPECT_DOUBLE_EQ(ff.getLi2(6)(1.), 4. * ff.getIW(6)(1.));
        EXPECT_DOUBLE_EQ(ff.getLi2(6)(1.5), 4. * ff.getIW(6)(1.5));
    }

    TEST(FFBLPRXPBase, Mis) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        // Evaluate the FF class
        TestableFFBLPRXPBase ff;
        ff.setSettingsHandler(set);
        ff.defineSettings();
        ff.calcConstants();
        auto cs = ff.getProtectedVars();
        EXPECT_DOUBLE_EQ(ff.getMi(8)(1.), cs[3] + 6 * cs[2] / (2.));
        EXPECT_DOUBLE_EQ(ff.getMi(8)(1.5), cs[3] + 6 * cs[2] / (2.5) - 2. * (0.5) * ff.getIW(6)(1.5) -
                                               2. * (2. * ff.getIW(2)(1.5) - 1.) * 0.5 / 2.5);
        EXPECT_DOUBLE_EQ(ff.getMi(9)(1.), 3. * cs[2] / 2. + 2. * ff.getIW(6)(1.));
        EXPECT_DOUBLE_EQ(ff.getMi(9)(1.5),
                         3. * cs[2] / 2.5 + 2. * ff.getIW(6)(1.5) - (2. * ff.getIW(2)(1.5) - 1.) * 0.5 / 2.5);
        EXPECT_DOUBLE_EQ(ff.getMi(10)(1.), cs[3] / 3. - cs[2] * (5) / (2 * 2) + 2 * 3 * ff.getIW(7)(1.) -
                                               (2. * ff.getIW(2)(1) - 1.) / 2.);
        EXPECT_DOUBLE_EQ(ff.getMi(10)(1.5), cs[3] / 3. - cs[2] * (5.5) / (2 * 2.5) + 2 * 3.5 * ff.getIW(7)(1.5) -
                                                (2. * ff.getIW(2)(1.5) - 1.) / 2.5);
    }

} // namespace Hammer
