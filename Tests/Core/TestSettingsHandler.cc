///
/// @file  TestSettingsHandler.cc
/// @brief Tests for SettingsHandler
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include <fstream>
#include <vector>

#include "gtest/gtest.h"

#include "Hammer/Tools/Pdg.hh"
#include "Hammer/SettingsHandler.hh"

using namespace std;

namespace Hammer {

    TEST(SettingsHandlerTest, addgetSetting) {

        SettingsHandler seth;
        seth.addSetting<bool>("WriterHandler", "WriteEvents", true);
        seth.addSetting<int>("WriterHandler", "MaxEvents", 100);
        seth.addSetting<double>("SimulationHandler:BFieldValue", 8.5);
        seth.addSetting<std::string>("Run:RunName", "Run01");
        EXPECT_EQ(*seth.getNamedSettingValue<std::string>("Run", "RunName"), "Run01");
        EXPECT_EQ(seth.getNamedSettingValue<bool>("Run", "RunName"), static_cast<bool*>(nullptr));
        EXPECT_TRUE(seth.getNamedSettingValue<bool>("WriterHandler:WriteEvents"));
    }

    TEST(SettingsHandlerTest, saveDefaultSettings) {
        SettingsHandler seth;
        seth.addSetting<bool>("WriterHandler", "WriteEvents", true);
        seth.addSetting<int>("WriterHandler", "MaxEvents", 100);
        seth.addSetting<double>("SimulationHandler:BFieldValue", 8.5);
        seth.addSetting<std::string>("Run:RunName", "Run01");
        seth.saveSettings("test.yaml", true);
        std::ifstream ifs("test.yaml");
        std::ifstream ifs2("test_defaults.yaml");
        std::string content((std::istreambuf_iterator<char>(ifs)), (std::istreambuf_iterator<char>()));
        std::string content2((std::istreambuf_iterator<char>(ifs2)), (std::istreambuf_iterator<char>()));
        EXPECT_EQ(content, content2);
    }

    TEST(SettingsHandlerTest, saveSettings) {
        SettingsHandler seth;
        seth.addSetting<bool>("WriterHandler", "WriteEvents", true);
        seth.addSetting<int>("WriterHandler", "MaxEvents", 100);
        seth.addSetting<double>("SimulationHandler:BFieldValue", 8.5);
        seth.addSetting<std::string>("Run:RunName", "Run01");
        seth.saveSettings("test.yaml", false);
        std::ifstream ifs("test.yaml");
        std::ifstream ifs2("test_settings.yaml");
        std::string content((std::istreambuf_iterator<char>(ifs)), (std::istreambuf_iterator<char>()));
        std::string content2((std::istreambuf_iterator<char>(ifs2)), (std::istreambuf_iterator<char>()));
        EXPECT_EQ(content, content2);
    }

    TEST(SettingsHandlerTest, loadSettingsEmpty) {
        SettingsHandler seth;
        seth.readSettings("test_settings.yaml");
        EXPECT_EQ(*seth.getNamedSettingValue<std::string>("Run", "RunName"), "Run01");
        EXPECT_EQ(*seth.getNamedSettingValue<double>("SimulationHandler:BFieldValue"), 8.5);
    }

    TEST(SettingsHandlerTest, loadSettingsPresent) {
        SettingsHandler seth;
        seth.addSetting<bool>("WriterHandler", "WriteEvents", false);
        seth.addSetting<int>("WriterHandler", "MaxEvents", 1000);
        seth.addSetting<double>("SimulationHandler:BFieldValue", 6.5);
        seth.addSetting<std::string>("Run:RunName", "Run02");
        seth.readSettings("test_settings.yaml");
        EXPECT_EQ(*seth.getNamedSettingValue<int>("WriterHandler", "MaxEvents"), 100);
        EXPECT_EQ(*seth.getNamedSettingValue<double>("SimulationHandler:BFieldValue"), 8.5);
    }

    TEST(SettingsHandlerTest, saveSettingMatrix) {
        SettingsHandler seth;
        vector<vector<double>> tmp = {
            {5., 6., 7., 8., 9.}, {15., 16., 17., 18., 19.}, {25., 26., 27., 28., 29.}, {35., 36., 37., 38., 39.}};
        seth.addSetting<vector<vector<double>>>("Parent", "Matrix", tmp);
        auto* tmpset = seth.getNamedSettingValue<vector<vector<double>>>("Parent", "Matrix", WTerm::COMMON);
        EXPECT_FALSE(tmpset == nullptr);
        if (tmpset != nullptr) {
            double val = (*tmpset)[1][1];
            EXPECT_NEAR(val, 16., 1.e-5);
        }

        seth.saveSettings("test.yaml");
        std::ifstream ifs("test.yaml");
        std::ifstream ifs2("test_matrix.yaml");
        std::string content((std::istreambuf_iterator<char>(ifs)), (std::istreambuf_iterator<char>()));
        std::string content2((std::istreambuf_iterator<char>(ifs2)), (std::istreambuf_iterator<char>()));
        EXPECT_EQ(content, content2);
    }

    TEST(SettingsHandlerTest, loadSettingMatrix) {
        SettingsHandler seth;
        vector<vector<double>> tmp = {
            {5., 6., 7., 8., 9.}, {55., 56., 57., 58., 59.}, {25., 26., 27., 28., 29.}, {35., 36., 37., 38., 39.}};
        seth.addSetting<vector<vector<double>>>("Parent", "Matrix", tmp);
        auto* tmpset = seth.getNamedSettingValue<vector<vector<double>>>("Parent", "Matrix", WTerm::COMMON);
        EXPECT_FALSE(tmpset == nullptr);
        if (tmpset != nullptr) {
            double val = (*tmpset)[1][1];
            EXPECT_NEAR(val, 56., 1.e-5);
        }
        seth.readSettings("test_matrix.yaml");
        auto* tmpset2 = seth.getNamedSettingValue<vector<vector<double>>>("Parent", "Matrix");
        EXPECT_FALSE(tmpset2 == nullptr);
        if (tmpset2 != nullptr) {
            double val = (*tmpset2)[1][1];
            EXPECT_NEAR(val, 16., 1.e-5);
        }
    }

    TEST(SettingsHandlerTest, setMassWidth) {
        SettingsHandler seth;
        auto& pdg = PID::instance();
        EXPECT_NEAR(pdg.getMass("K0L"), 0.497611, 1.e-4);
        seth.readSettings("test_masses.yaml");
        EXPECT_NEAR(pdg.getMass("K0L"), 7.8, 1.e-4);
        seth.changeSetting<double>("Masses", "K0L", 4.3);
        EXPECT_NEAR(pdg.getMass("K0L"), 4.3, 1.e-4);
    }

    TEST(SettingsHandlerTest, changeSetting) {
        SettingsHandler seth;
        seth.addSetting<bool>("WriterHandler", "WriteEvents", true);
        EXPECT_TRUE(*seth.getNamedSettingValue<bool>("WriterHandler:WriteEvents"));
        seth.changeSetting<bool>("WriterHandler", "WriteEvents", false);
        EXPECT_FALSE(*seth.getNamedSettingValue<bool>("WriterHandler:WriteEvents"));
        EXPECT_THROW(seth.changeSetting<double>("WriterHandler", "WriteEvents", 3.1), InitializationError);
    }


} // namespace Hammer
