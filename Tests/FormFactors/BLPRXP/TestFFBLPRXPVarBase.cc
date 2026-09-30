///
/// @file  TestFFBLPRXPVarBase.cc
/// @brief Tests for FFBLPRXPVarBase
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"
#include "Hammer/FormFactors/BLPRXP/FFBLPRXPVarBase.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/Utils.hh"
// May need other stuff here

using namespace std;
using namespace complex_literals;

namespace Hammer {

    class TestableFFBLPRXPVarBase : public FFBLPRXPVarBase {

    public:

        TestableFFBLPRXPVarBase() = default;

        using FFBLPRXPVarBase::calcConstants;
        using FFBLPRXPVarBase::defineSettings;
        using FFBLPRXPVarBase::Xi;
        using FFBLPRXPVarBase::XiLinearCoeffs;
        using FFBLPRXPVarBase::ZofW;


        void evalAtPSPoint(const vector<double>& /*point*/, const vector<double>& /*masses*/) override {
        }
        unique_ptr<FormFactorBase> clone(const string& /*label*/) override {
            return unique_ptr<FormFactorBase>{};
        }
    };

    TEST(FFBLPRXPVarBase, XiExpansion) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        // Evaluate the FF class
        TestableFFBLPRXPVarBase ff;
        ff.setSettingsHandler(set);
        ff.defineSettings();
        ff.calcConstants();
        double z1 = ff.ZofW(1.);
        double z1p25 = ff.ZofW(1.25);
        double a2 = 1.509 * 1.509 / 2.;
        double a4 = a2 * a2;
        double RhoStSq = 1.104;
        double cSt = 2.392;
        auto v1 = ff.XiLinearCoeffs(1.);
        auto v1p5 = ff.XiLinearCoeffs(1.25);
        for (size_t i = 2; i < 10; ++i) {
            EXPECT_NEAR(v1[i], 0., 1.e-5);
            EXPECT_NEAR(v1p5[i], 0., 1.e-5);
        }
        EXPECT_NEAR(v1[0], 0., 1.e-5);
        EXPECT_NEAR(v1[1], 0., 1.e-5);
        double n = (1. - 8 * a2 * RhoStSq * z1p25 + 16. * (2 * cSt * a4 - RhoStSq * a2) * z1p25 * z1p25);
        double d = (1. - 8 * a2 * RhoStSq * z1 + 16. * (2 * cSt * a4 - RhoStSq * a2) * z1 * z1);
        EXPECT_DOUBLE_EQ(ff.Xi(1.25), n / d);
        EXPECT_NEAR(v1p5[0], -8 * a2 * (z1p25 * (1 + 2 * z1p25) / n - z1 * (1 + 2 * z1) / d), 1.e-4);
        EXPECT_NEAR(v1p5[1], 32 * a4 * (z1p25 * z1p25 / n - z1 * z1 / d), 1.e-6);
    }


} // namespace Hammer
