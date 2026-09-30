///
/// @file  TestDictionaryManager.cc
/// @brief Tests for DictionaryManager
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include <fstream>
#include <vector>

#include "gtest/gtest.h"

#define DICTMANAGER_FRIENDS                                \
    FRIEND_TEST(DictionaryManagerTest, readDecays);        \
    FRIEND_TEST(DictionaryManagerTest, loadFFSchemes);     \
    FRIEND_TEST(DictionaryManagerTest, saveDecays);        \
    FRIEND_TEST(DictionaryManagerTest, saveAndTestChanges)

#define PROCESSDEFINITIONS_FRIENDS                         \
    FRIEND_TEST(DictionaryManagerTest, readDecays);        \
    FRIEND_TEST(DictionaryManagerTest, loadFFSchemes);     \
    FRIEND_TEST(DictionaryManagerTest, saveDecays);        \
    FRIEND_TEST(DictionaryManagerTest, saveAndTestChanges)

#define SCHEMEDEFINITIONS_FRIENDS                          \
    FRIEND_TEST(DictionaryManagerTest, readDecays);        \
    FRIEND_TEST(DictionaryManagerTest, loadFFSchemes);     \
    FRIEND_TEST(DictionaryManagerTest, saveDecays);        \
    FRIEND_TEST(DictionaryManagerTest, saveAndTestChanges)

#include "Hammer/Tools/Pdg.hh"
#include "Hammer/DictionaryManager.hh"
#include "Hammer/ProcessDefinitions.hh"
#include "Hammer/SchemeDefinitions.hh"

using namespace std;

namespace Hammer {

    TEST(DictionaryManagerTest, readDecays) {
        DictionaryManager seth;
        DictionaryManager evalseth;
        seth.readDecays("test_decsettings.yaml");
        std::vector<std::set<HashId>> va = seth.processDefs().includedDecaySignatures();
        std::vector<std::set<HashId>> vf = seth.processDefs().forbiddenDecaySignatures();
        std::vector<std::vector<std::string>> includednames = {{"D*Dpi", "TauPiNu"}, {"B+DTau+Nu"}};
        std::vector<std::set<HashId>> Evalincluded = ProcessDefinitions::decaySignatures(includednames);
        std::vector<std::vector<std::string>> forbiddennames = {{"B+DTau+Nu", "TauEllNuNu"}};
        std::vector<std::set<HashId>> Evalforbidden = ProcessDefinitions::decaySignatures(forbiddennames);
        for (size_t i = 0; i < va.size(); ++i) {
            EXPECT_EQ(va[i], Evalincluded[i]);
        }
        for (size_t i = 0; i < vf.size(); ++i) {
            EXPECT_EQ(vf[i], Evalforbidden[i]);
        }
    }

    TEST(DictionaryManagerTest, loadFFSchemes) {
        DictionaryManager seth;
        seth.readDecays("test_decsettings.yaml");
        std::vector<std::string> v = seth.schemeDefs().getFFSchemeNames();
        std::vector<std::string> Eval;
        Eval.emplace_back("MixedScheme");
        Eval.emplace_back("AllHQET");
        std::sort(Eval.begin(), Eval.end());
        for (size_t i = 0; i < v.size(); ++i) {
            EXPECT_EQ(v[i], Eval[i]);
        }
    }

    TEST(DictionaryManagerTest, saveDecays) {
        DictionaryManager seth;
        seth.readDecays("test_decsettings.yaml");
        seth.saveDecays("testdec.yaml");
        DictionaryManager seth2;
        seth2.readDecays("testdec.yaml");
        std::vector<std::set<HashId>> aldecs = seth2.processDefs().includedDecaySignatures();
        std::vector<std::set<HashId>> fdecs = seth2.processDefs().forbiddenDecaySignatures();
        for (size_t i = 0; i < aldecs.size(); ++i) {
            EXPECT_EQ(seth.processDefs().includedDecaySignatures(), seth2.processDefs().includedDecaySignatures());
        }
        for (size_t j = 0; j < fdecs.size(); ++j) {
            EXPECT_EQ(seth.processDefs().forbiddenDecaySignatures(), seth2.processDefs().forbiddenDecaySignatures());
        }
        std::vector<std::string> c1 = seth.schemeDefs().getFFSchemeNames();
        std::vector<std::string> c2 = seth2.schemeDefs().getFFSchemeNames();
        for (size_t k = 0; k < c2.size(); ++k) {
            EXPECT_EQ(c1[k], c2[k]);
            EXPECT_EQ(seth.schemeDefs().getScheme(c2[k]), seth2.schemeDefs().getScheme(c2[k]));
        }
    }

    TEST(DictionaryManagerTest, saveAndTestChanges) {
        DictionaryManager seth;
        seth.readDecays("test_decsettings.yaml");
        std::vector<std::string> va = {"B+D*Tau+Nu"};
        seth.processDefs().addIncludedDecay(va);
        std::vector<std::string> vf = {"B+DMu+Nu"};
        seth.processDefs().addForbiddenDecay(vf);
        seth.schemeDefs().removeFFScheme("AllHQET");
        std::map<std::string, std::string> testschemes;
        testschemes.insert(std::pair<std::string, std::string>("BD", "ISGW2"));
        testschemes.insert(std::pair<std::string, std::string>("BD*", "CLN"));
        seth.schemeDefs().addFFScheme("TestScheme", testschemes);
        std::map<std::string, std::string> testdenom;
        testdenom.insert(std::pair<std::string, std::string>("BD", "HQET"));
        testdenom.insert(std::pair<std::string, std::string>("BD*", "HQET"));
        seth.schemeDefs().setFFInputScheme(testdenom);
        seth.saveDecays("testdec2.yaml");
        DictionaryManager seth2;
        seth2.readDecays("testdec2.yaml");
        DictionaryManager seth3;
        seth3.readDecays("test_decsettings2.yaml");
        std::vector<std::set<HashId>> aldecs = seth2.processDefs().includedDecaySignatures();
        std::vector<std::set<HashId>> fdecs = seth2.processDefs().forbiddenDecaySignatures();
        std::vector<std::string> c1 = seth2.schemeDefs().getFFSchemeNames();
        std::vector<std::string> c2 = seth3.schemeDefs().getFFSchemeNames();
        std::map<HashId, std::map<std::string, std::vector<std::string>>> m1 = seth2.schemeDefs().getFFDuplicates();
        std::map<HashId, std::map<std::string, std::vector<std::string>>> m2 = seth3.schemeDefs().getFFDuplicates();
        EXPECT_EQ(seth2.processDefs().includedDecaySignatures(), seth3.processDefs().includedDecaySignatures());
        EXPECT_EQ(seth2.processDefs().forbiddenDecaySignatures(), seth3.processDefs().forbiddenDecaySignatures());
        for (size_t k = 0; k < c2.size(); ++k) {
            EXPECT_EQ(c1[k], c2[k]);
            EXPECT_EQ(seth2.schemeDefs().getScheme(c2[k]), seth3.schemeDefs().getScheme(c2[k]));
        }
        EXPECT_EQ(m1, m2);
    }


} // namespace Hammer
