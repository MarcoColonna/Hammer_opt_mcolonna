///
/// @file  TestPdg.cc
/// @brief Tests for Pdg
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#define PID_FRIENDS                                  \
    FRIEND_TEST(PIDTest, fundamentalIDExtraBits);    \
    FRIEND_TEST(PIDTest, isDiQuarkAllFalseBranches); \
    FRIEND_TEST(PIDTest, isPentaquarkAllBranches)

#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Tools/ParticleUtils.hh"
#include <boost/functional/hash.hpp>
#include <algorithm>

using namespace std;
namespace Hammer {

    TEST(PIDTest, getMass) {
        auto& pdg = PID::instance();
        EXPECT_NEAR(pdg.getMass(PID::ANTIMUON), 0.106, 0.001);
        EXPECT_NEAR(pdg.getMass("Mu-"), 0.106, 0.001);
    }

    TEST(PIDTest, toPdg) {
        auto& pdg = PID::instance();
        EXPECT_EQ(pdg.toPdgCode("Mu-"), 13);
        EXPECT_EQ(pdg.toPdgCode("Ell"), 0);
        auto codes = pdg.toPdgList("Ell");
        EXPECT_EQ(codes.size(), 4);
        sort(codes.begin(), codes.end());
        EXPECT_EQ(codes[0], -13);
        EXPECT_EQ(codes[3], 13);
        EXPECT_EQ(codes[1], -11);
        EXPECT_EQ(codes[2], 11);
    }

    TEST(PIDTest, getCharges) {
        auto& pdg = PID::instance();
        EXPECT_EQ(pdg.getThreeCharge(PID::MUON), -3);
        EXPECT_EQ(pdg.getThreeCharge({PID::MUON, PID::PIPLUS, PID::NEUTRON}), 0);
        EXPECT_EQ(pdg.getLeptonNumber(PID::MUON), 1);
        EXPECT_EQ(pdg.getLeptonNumber({PID::MUON, PID::PIPLUS, PID::NEUTRON}), 1);
        EXPECT_EQ(pdg.getLeptonNumber({PID::MUON, PID::NU_EBAR}), 0);
        EXPECT_EQ(pdg.getBaryonNumber(PID::MUON), 0);
        EXPECT_EQ(pdg.getBaryonNumber({PID::MUON, PID::PIPLUS, PID::NEUTRON}), 1);
        EXPECT_EQ(pdg.getBaryonNumber({PID::NEUTRON, PID::ANTIPROTON}), 0);
        EXPECT_EQ(get<0>(pdg.getLeptonFlavorNumber(PID::MUON)), 0);
        EXPECT_EQ(get<1>(pdg.getLeptonFlavorNumber(PID::MUON)), 1);
        EXPECT_EQ(get<2>(pdg.getLeptonFlavorNumber(PID::MUON)), 0);
        EXPECT_EQ(get<0>(pdg.getLeptonFlavorNumber({PID::MUON, PID::NU_EBAR})), -1);
        EXPECT_EQ(get<1>(pdg.getLeptonFlavorNumber({PID::MUON, PID::NU_EBAR})), 1);
        EXPECT_EQ(get<2>(pdg.getLeptonFlavorNumber({PID::MUON, PID::NU_EBAR})), 0);
    }

    TEST(PIDTest, getValidVertices) {
        auto& pdg = PID::instance();
        auto list1 = pdg.expandToValidVertices("B+D0barNutTau+");
        auto list2 = pdg.expandToValidVertices("BDTauNu");
        auto list3 = pdg.expandToValidVertices("TauEllNuNu");
        auto list4 = pdg.expandToValidVertices("TauEllNutNu");
        auto list5 = pdg.expandToValidVertexUIDs("BD*", true);
        auto list6 = pdg.expandToValidVertexUIDs("TauPiPiPi", true);
        EXPECT_EQ(list1.size(), 1);
        EXPECT_EQ(list2.size(), 8);
        EXPECT_EQ(list3.size(), 4);
        EXPECT_EQ(list4.size(), 2);
        EXPECT_EQ(list5.size(), 14);
        EXPECT_EQ(list6.size(), 4);

        struct less_than_hash {
            bool operator()(const pair<PdgId, vector<PdgId>>& vertex1, const pair<PdgId, vector<PdgId>>& vertex2) {
                return (processID(vertex1.first, combineDaughters(vertex1.second)) <
                        processID(vertex2.first, combineDaughters(vertex2.second)));
            }
        };

        sort(list1.begin(), list1.end(), less_than_hash());
        sort(list2.begin(), list2.end(), less_than_hash());
        sort(list3.begin(), list3.end(), less_than_hash());
        sort(list4.begin(), list4.end(), less_than_hash());
        sort(list5.begin(), list5.end());
        sort(list6.begin(), list6.end());

        vector<pair<PdgId, vector<PdgId>>> list1Exp = {{521, {-421, 16, -15}}}; //, {-521, {421,-16,15}}};
        sort(list1Exp.begin(), list1Exp.end(), less_than_hash());
        for (size_t idx = 0; idx < list1.size(); ++idx) {
            auto parent = get<0>(list1[idx]);
            EXPECT_EQ(parent, get<0>(list1Exp[idx]));
            auto daughters = get<1>(list1[idx]);
            for (size_t idx1 = 0; idx1 < daughters.size(); ++idx1) {
                EXPECT_EQ(daughters[idx1], get<1>(list1Exp[idx])[idx1]);
            }
        }

        vector<pair<PdgId, vector<PdgId>>> list2Exp = {
            {511, {411, -16, 15}},  {-511, {-411, 16, -15}}, {511, {-411, 16, -15}},  {-511, {411, -16, 15}},
            {-521, {421, -16, 15}}, {521, {-421, 16, -15}},  {-521, {-421, -16, 15}}, {521, {421, 16, -15}}};
        sort(list2Exp.begin(), list2Exp.end(), less_than_hash());
        for (size_t idx = 0; idx < list2.size(); ++idx) {
            auto parent = get<0>(list2[idx]);
            EXPECT_EQ(parent, get<0>(list2Exp[idx]));
            auto daughters = get<1>(list2[idx]);
            for (size_t idx1 = 0; idx1 < daughters.size(); ++idx1) {
                EXPECT_EQ(daughters[idx1], get<1>(list2Exp[idx])[idx1]);
            }
        }

        vector<pair<PdgId, vector<PdgId>>> list3Exp = {
            {15, {16, -14, 13}}, {-15, {-16, 14, -13}}, {15, {16, -12, 11}}, {-15, {-16, 12, -11}}};
        sort(list3Exp.begin(), list3Exp.end(), less_than_hash());
        for (size_t idx = 0; idx < list3.size(); ++idx) {
            auto parent = get<0>(list3[idx]);
            EXPECT_EQ(parent, get<0>(list3Exp[idx]));
            auto daughters = get<1>(list3[idx]);
            for (size_t idx1 = 0; idx1 < daughters.size(); ++idx1) {
                EXPECT_EQ(daughters[idx1], get<1>(list3Exp[idx])[idx1]);
            }
        }

        //        vector<pair<PdgId, vector<PdgId>>> list4Exp = {{15, {16,-14,13}}, {-15, {-16,14,-13}},{15,
        //        {16,-12,11}}, {-15, {-16,12,-11}}};
        vector<pair<PdgId, vector<PdgId>>> list4Exp = {{15, {16, -14, 13}}, {15, {16, -12, 11}}};
        sort(list4Exp.begin(), list4Exp.end(), less_than_hash());
        for (size_t idx = 0; idx < list4.size(); ++idx) {
            auto parent = get<0>(list4[idx]);
            EXPECT_EQ(parent, get<0>(list4Exp[idx]));
            auto daughters = get<1>(list4[idx]);
            for (size_t idx1 = 0; idx1 < daughters.size(); ++idx1) {
                EXPECT_EQ(daughters[idx1], get<1>(list4Exp[idx])[idx1]);
            }
        }

        vector<pair<PdgId, vector<PdgId>>> list5Exp = {
            {511, {423}},  {-511, {-423}}, {511, {-423}},  {-511, {423}}, {-521, {-413}}, {521, {413}},   {511, {-413}},
            {-511, {413}}, {511, {413}},   {-511, {-413}}, {-521, {423}}, {521, {-423}},  {-521, {-423}}, {521, {423}}};
        sort(list5Exp.begin(), list5Exp.end(), less_than_hash());
        for (size_t idx = 0; idx < list5Exp.size(); ++idx) {
            auto daughters = get<1>(list5Exp[idx]);
            EXPECT_EQ(processID(get<0>(list5Exp[idx]), daughters), list5[idx]);
        }

        vector<pair<PdgId, vector<PdgId>>> list6Exp = {
            {15, {211, -211, -211}}, {-15, {211, 211, -211}}, {15, {-211, 111, 111}}, {-15, {211, 111, 111}}};
        sort(list6Exp.begin(), list6Exp.end(), less_than_hash());
        for (size_t idx = 0; idx < list6Exp.size(); ++idx) {
            auto daughters = get<1>(list6Exp[idx]);
            EXPECT_EQ(processID(get<0>(list6Exp[idx]), daughters), list6[idx]);
        }
    }


    TEST(PIDTest, getMassUnknown) {
        auto& pdg = PID::instance();
        EXPECT_NEAR(pdg.getMass(PID::ELECTRON), 5.11e-4, 1.e-6);
        EXPECT_NEAR(pdg.getMass(PID::TAU), 1.77686, 1.e-4);
        EXPECT_NEAR(pdg.getMass(PID::PIPLUS), 0.13957061, 1.e-6);
        EXPECT_NEAR(pdg.getMass(PID::BZERO), 5.27963, 1.e-4);
        EXPECT_NEAR(pdg.getMass(PID::LAMBDAB), 5.620, 1.e-4);
        EXPECT_NEAR(pdg.getMass(PID::BCPLUS), 6.275, 1.e-4);
        EXPECT_NEAR(pdg.getMass(PID::JPSI), 3.09690, 1.e-4);
        EXPECT_NEAR(pdg.getMass(-PID::BZERO), 5.27963, 1.e-4);
        EXPECT_EQ(pdg.getMass(99999), -1.);
    }

    TEST(PIDTest, getWidthByIdAndName) {
        auto& pdg = PID::instance();
        EXPECT_GT(pdg.getWidth(PID::TAU), 0.);
        EXPECT_NEAR(pdg.getWidth(PID::RHO0), 0.147, 1.e-4);
        EXPECT_NEAR(pdg.getWidth("Rho0"), 0.147, 1.e-4);
        EXPECT_EQ(pdg.getWidth(99999), -1.);
        EXPECT_EQ(pdg.getWidth("NotAParticle"), -1.);
    }

    TEST(PIDTest, setMassAndWidth) {
        auto& pdg = PID::instance();
        double origMass = pdg.getMass(PID::OMEGA);
        pdg.setMass(PID::OMEGA, 0.800);
        EXPECT_NEAR(pdg.getMass(PID::OMEGA), 0.800, 1.e-5);
        pdg.setMass(PID::OMEGA, origMass);

        double origMassStr = pdg.getMass("Omega");
        pdg.setMass("Omega", 0.801);
        EXPECT_NEAR(pdg.getMass("Omega"), 0.801, 1.e-5);
        pdg.setMass("Omega", origMassStr);

        double origWidth = pdg.getWidth(PID::OMEGA);
        pdg.setWidth(PID::OMEGA, 0.010);
        EXPECT_NEAR(pdg.getWidth(PID::OMEGA), 0.010, 1.e-5);
        pdg.setWidth(PID::OMEGA, origWidth);

        pdg.setWidth("Omega", 0.011);
        EXPECT_NEAR(pdg.getWidth("Omega"), 0.011, 1.e-5);
        pdg.setWidth("Omega", origWidth);
    }

    TEST(PIDTest, getSpinMultiplicity) {
        auto& pdg = PID::instance();
        EXPECT_EQ(pdg.getSpinMultiplicity(PID::ELECTRON), 2u);
        EXPECT_EQ(pdg.getSpinMultiplicity(PID::TAU), 2u);
        EXPECT_EQ(pdg.getSpinMultiplicity(PID::PHOTON), 2u);
        EXPECT_EQ(pdg.getSpinMultiplicity(PID::PIPLUS), 1u);
        EXPECT_EQ(pdg.getSpinMultiplicity(PID::DSTAR), 3u);
        EXPECT_EQ(pdg.getSpinMultiplicity(PID::PROTON), 2u);
        EXPECT_EQ(pdg.getSpinMultiplicity("Mu-"), 2u);
        EXPECT_EQ(pdg.getSpinMultiplicities({PID::PIPLUS, PID::ELECTRON}), 2u);
        EXPECT_EQ(pdg.getSpinMultiplicities(vector<string>{"Pi+", "Mu-"}), 2u);
    }

    TEST(PIDTest, toPdgListUnknown) {
        auto& pdg = PID::instance();
        auto result = pdg.toPdgList("NotAParticle");
        EXPECT_EQ(result.size(), 0u);
        EXPECT_EQ(pdg.toPdgCode("NotAParticle"), 0);
        EXPECT_EQ(pdg.toPdgCode("Ell"), 0);
        auto bList = pdg.toPdgList("B");
        EXPECT_EQ(bList.size(), 4u);
        auto piList = pdg.toPdgList("Pi");
        EXPECT_EQ(piList.size(), 3u);
    }

    TEST(PIDTest, isLepton) {
        EXPECT_TRUE(PID::isLepton(PID::ELECTRON));
        EXPECT_TRUE(PID::isLepton(PID::POSITRON));
        EXPECT_TRUE(PID::isLepton(PID::MUON));
        EXPECT_TRUE(PID::isLepton(PID::TAU));
        EXPECT_TRUE(PID::isLepton(PID::NU_E));
        EXPECT_TRUE(PID::isLepton(PID::NU_MU));
        EXPECT_TRUE(PID::isLepton(PID::NU_TAU));
        EXPECT_FALSE(PID::isLepton(PID::PIPLUS));
        EXPECT_FALSE(PID::isLepton(PID::PROTON));
        EXPECT_FALSE(PID::isLepton(PID::D0));
        EXPECT_FALSE(PID::isLepton(0));
    }

    TEST(PIDTest, isNeutrino) {
        EXPECT_TRUE(PID::isNeutrino(PID::NU_E));
        EXPECT_TRUE(PID::isNeutrino(PID::NU_EBAR));
        EXPECT_TRUE(PID::isNeutrino(PID::NU_MU));
        EXPECT_TRUE(PID::isNeutrino(PID::NU_MUBAR));
        EXPECT_TRUE(PID::isNeutrino(PID::NU_TAU));
        EXPECT_TRUE(PID::isNeutrino(PID::NU_TAUBAR));
        EXPECT_FALSE(PID::isNeutrino(PID::ELECTRON));
        EXPECT_FALSE(PID::isNeutrino(PID::MUON));
        EXPECT_FALSE(PID::isNeutrino(PID::PIPLUS));
    }

    TEST(PIDTest, isMesonIsBaryon) {
        EXPECT_TRUE(PID::isMeson(PID::PIPLUS));
        EXPECT_TRUE(PID::isMeson(PID::D0));
        EXPECT_TRUE(PID::isMeson(PID::DSTAR));
        EXPECT_TRUE(PID::isMeson(PID::BZERO));
        EXPECT_TRUE(PID::isMeson(PID::K0L));
        EXPECT_TRUE(PID::isMeson(PID::K0S));
        EXPECT_FALSE(PID::isMeson(PID::PROTON));
        EXPECT_FALSE(PID::isMeson(PID::ELECTRON));

        EXPECT_TRUE(PID::isBaryon(PID::PROTON));
        EXPECT_TRUE(PID::isBaryon(PID::NEUTRON));
        EXPECT_TRUE(PID::isBaryon(PID::LAMBDACPLUS));
        EXPECT_TRUE(PID::isBaryon(PID::LAMBDAB));
        EXPECT_FALSE(PID::isBaryon(PID::PIPLUS));
        EXPECT_FALSE(PID::isBaryon(PID::ELECTRON));

        EXPECT_TRUE(PID::isHadron(PID::PIPLUS));
        EXPECT_TRUE(PID::isHadron(PID::PROTON));
        EXPECT_FALSE(PID::isHadron(PID::ELECTRON));
    }

    TEST(PIDTest, getThreeChargeExtended) {
        EXPECT_EQ(PID::getThreeCharge(PID::ELECTRON), -3);
        EXPECT_EQ(PID::getThreeCharge(PID::POSITRON), 3);
        EXPECT_EQ(PID::getThreeCharge(PID::NU_E), 0);
        EXPECT_EQ(PID::getThreeCharge(PID::PIPLUS), 3);
        EXPECT_EQ(PID::getThreeCharge(PID::PIMINUS), -3);
        EXPECT_EQ(PID::getThreeCharge(PID::PI0), 0);
        EXPECT_EQ(PID::getThreeCharge(PID::PROTON), 3);
        EXPECT_EQ(PID::getThreeCharge(PID::ANTIPROTON), -3);
        EXPECT_EQ(PID::getThreeCharge(PID::NEUTRON), 0);
        EXPECT_EQ(PID::getThreeCharge(0), 0);
    }

    TEST(PIDTest, getLeptonFlavorNumberExtended) {
        auto tau = PID::getLeptonFlavorNumber(PID::TAU);
        EXPECT_EQ(get<0>(tau), 0);
        EXPECT_EQ(get<1>(tau), 0);
        EXPECT_EQ(get<2>(tau), 1);

        auto antitau = PID::getLeptonFlavorNumber(PID::ANTITAU);
        EXPECT_EQ(get<2>(antitau), -1);

        auto electron = PID::getLeptonFlavorNumber(PID::ELECTRON);
        EXPECT_EQ(get<0>(electron), 1);
        EXPECT_EQ(get<1>(electron), 0);
        EXPECT_EQ(get<2>(electron), 0);

        auto pion = PID::getLeptonFlavorNumber(PID::PIPLUS);
        EXPECT_EQ(get<0>(pion), 0);
        EXPECT_EQ(get<1>(pion), 0);
        EXPECT_EQ(get<2>(pion), 0);
    }

    TEST(PIDTest, getPartialWidths) {
        auto& pdg = PID::instance();
        auto pws = pdg.getPartialWidths();
        EXPECT_GT(pws.size(), 0u);
        for (const auto& elem : pws) {
            EXPECT_GT(elem.second, 0.);
        }
    }

    TEST(PIDTest, expandToValidVerticesLbLcTauNu) {
        auto& pdg = PID::instance();
        auto list = pdg.expandToValidVertices("LbLcTauNu");
        EXPECT_GT(list.size(), 0u);
    }

    TEST(PIDTest, expandToValidVerticesEmpty) {
        auto& pdg = PID::instance();
        auto list = pdg.expandToValidVertices("NotAParticleName");
        EXPECT_EQ(list.size(), 0u);
    }

    TEST(PIDTest, expandToValidVertexUIDsNoHadOnly) {
        auto& pdg = PID::instance();
        auto uids = pdg.expandToValidVertexUIDs("BDTauNu", false);
        EXPECT_GT(uids.size(), 0u);
    }

    TEST(PIDTest, isHadronExtraParticles) {
        EXPECT_FALSE(PID::isHadron(PID::ELECTRON));
        EXPECT_TRUE(PID::isHadron(PID::PIPLUS));
        EXPECT_TRUE(PID::isHadron(PID::PROTON));
        EXPECT_FALSE(PID::isHadron(0));
    }

    TEST(PIDTest, isMesonEdgeCases) {
        EXPECT_TRUE(PID::isMeson(PID::PIMINUS));
        EXPECT_FALSE(PID::isMeson(PID::ELECTRON));
        EXPECT_FALSE(PID::isMeson(0));
        EXPECT_TRUE(PID::isMeson(PID::K0L));
        EXPECT_FALSE(PID::isMeson(-PID::K0L));
        EXPECT_TRUE(PID::isMeson(150));
        EXPECT_TRUE(PID::isMeson(350));
        EXPECT_TRUE(PID::isMeson(510));
        EXPECT_TRUE(PID::isMeson(530));
        EXPECT_TRUE(PID::isMeson(110));
        EXPECT_TRUE(PID::isMeson(990));
        EXPECT_TRUE(PID::isMeson(9990));
    }

    TEST(PIDTest, isBaryonSpecialCases) {
        EXPECT_TRUE(PID::isBaryon(2110));
        EXPECT_TRUE(PID::isBaryon(2210));
        EXPECT_FALSE(PID::isBaryon(PID::ELECTRON));
        EXPECT_FALSE(PID::isBaryon(0));
    }

    TEST(PIDTest, getThreeChargeDiQuark) {
        EXPECT_EQ(PID::getThreeCharge(5501), -2);
        EXPECT_EQ(PID::getThreeCharge(5503), -2);
    }

    TEST(PIDTest, getThreeChargeIon) {
        EXPECT_EQ(PID::getThreeCharge(10000001), 0);
    }

    TEST(PIDTest, setMassUnknownParticle) {
        auto& pdg = PID::instance();
        pdg.setMass(99998, 1.0);
        EXPECT_NEAR(pdg.getMass(99998), 1.0, 1e-10);
        pdg.setWidth(99998, 0.01);
        EXPECT_NEAR(pdg.getWidth(99998), 0.01, 1e-10);
    }

    TEST(PIDTest, addBRAndGetPartialWidths) {
        auto& pdg = PID::instance();
        auto pws = pdg.getPartialWidths();
        EXPECT_GT(pws.size(), 0u);
    }

    TEST(PIDTest, isLeptonExtraBits) {
        EXPECT_FALSE(PID::isLepton(10000011));
    }

    TEST(PIDTest, isNeutrinoExtraBits) {
        EXPECT_FALSE(PID::isNeutrino(10000012));
    }

    TEST(PIDTest, toPdgCodeWithMultipleResults) {
        auto& pdg = PID::instance();
        EXPECT_EQ(pdg.toPdgCode("Tau"), 0);
        EXPECT_EQ(pdg.toPdgCode("Tau-"), PID::TAU);
    }

    TEST(PIDTest, getSpinMultiplicityBaryon) {
        auto& pdg = PID::instance();
        EXPECT_EQ(pdg.getSpinMultiplicity(PID::LAMBDAB), 2u);
        EXPECT_EQ(pdg.getSpinMultiplicity(PID::LAMBDACPLUS), 2u);
    }

    TEST(PIDTest, getSpinMultiplicityByName) {
        auto& pdg = PID::instance();
        EXPECT_EQ(pdg.getSpinMultiplicity("Pi+"), 1u);
        EXPECT_EQ(pdg.getSpinMultiplicity("D*0"), 3u);
    }

    TEST(PIDTest, getSpinMultiplicitiesStrings) {
        auto& pdg = PID::instance();
        EXPECT_EQ(pdg.getSpinMultiplicities(vector<string>{"Pi+", "D*0"}), 3u);
    }

    TEST(PIDTest, getLeptonNumberVector) {
        EXPECT_EQ(PID::getLeptonNumber(PID::PROTON), 0);
        EXPECT_EQ(PID::getLeptonNumber(-PID::ELECTRON), -1);
    }

    TEST(PIDTest, getBaryonNumberNegative) {
        EXPECT_EQ(PID::getBaryonNumber(-PID::PROTON), -1);
        EXPECT_EQ(PID::getBaryonNumber(PID::PIPLUS), 0);
    }

    TEST(PIDTest, getLeptonFlavorNumberVector) {
        auto res = PID::getLeptonFlavorNumber(vector<PdgId>{PID::TAU, PID::ANTITAU});
        EXPECT_EQ(get<2>(res), 0);
    }

    TEST(PIDTest, getSpinMultiplicityFallback) {
        auto& pdg = PID::instance();
        // WBOSON (24): abspid>18, not meson or baryon → fallback returns 1
        EXPECT_EQ(pdg.getSpinMultiplicity(PID::WBOSON), 1u);
    }

    TEST(PIDTest, getThreeChargeStrangeMeson) {
        // K+ (321): q2=3 (strange) → strange meson branch → charge 3
        EXPECT_EQ(PID::getThreeCharge(PID::KPLUS), 3);
    }

    TEST(PIDTest, getThreeChargeUnknownParticle) {
        // -RHO0 (-113): not meson (q2==q3 and pid<0), not diquark, not baryon → unknown → 0
        EXPECT_EQ(PID::getThreeCharge(-PID::RHO0), 0);
    }

    TEST(PIDTest, getThreeChargeKL) {
        // K0L (130): digit(nj)==0 → return 0 immediately (covers line 162 branch)
        EXPECT_EQ(PID::getThreeCharge(PID::K0L), 0);
        // pid=100: fundamentalID returns 100 → sid=100, hits ch100 table branch (line 160-161)
        EXPECT_EQ(PID::getThreeCharge(100), 0);
    }

    TEST(PIDTest, fundamentalIDExtraBits) {
        // extraBits > 0 → returns 0
        EXPECT_EQ(PID::fundamentalID(10000001), 0);
    }

    TEST(PIDTest, isMesonExtraBitsAndFundamentalID) {
        EXPECT_FALSE(PID::isMeson(10000001)); // extraBits > 0
        EXPECT_FALSE(PID::isMeson(10011));    // fundamentalID = 11, in (0, 100]
    }

    TEST(PIDTest, isBaryonExtraBitsAndFundamentalID) {
        EXPECT_FALSE(PID::isBaryon(10000001)); // extraBits > 0
        EXPECT_FALSE(PID::isBaryon(10011));    // fundamentalID = 11, in (0, 100]
    }

    TEST(PIDTest, isDiQuarkAllFalseBranches) {
        EXPECT_FALSE(PID::isDiQuark(10000001)); // extraBits > 0
        EXPECT_FALSE(PID::isDiQuark(11));       // abspid <= 100
        EXPECT_FALSE(PID::isDiQuark(10011));    // fundamentalID = 11, in (0, 100]
    }

    TEST(PIDTest, isHadronExtraBits) {
        EXPECT_FALSE(PID::isHadron(10000001)); // extraBits > 0
    }

    TEST(PIDTest, isPentaquarkAllBranches) {
        // extraBits > 0
        EXPECT_FALSE(PID::isPentaquark(10000001));
        // digit(n) != 9 (ordinary PID)
        EXPECT_FALSE(PID::isPentaquark(PID::PIPLUS));
        // digit(nr) == 9
        EXPECT_FALSE(PID::isPentaquark(9900001));
        // digit(nj) == 9
        EXPECT_FALSE(PID::isPentaquark(9110009));
        // digit(nq1) == 0
        EXPECT_FALSE(PID::isPentaquark(9110001));
        // digit(nq2) == 0
        EXPECT_FALSE(PID::isPentaquark(9111001));
        // digit(nq3) == 0
        EXPECT_FALSE(PID::isPentaquark(9111101));
        // digit(nj) == 0
        EXPECT_FALSE(PID::isPentaquark(9111110));
        // digit(nq2) > digit(nq1): nq2=3, nq1=1
        EXPECT_FALSE(PID::isPentaquark(9211321));
        // digit(nq1) > digit(nl): nq1=3, nl=2
        EXPECT_FALSE(PID::isPentaquark(9423221));
        // digit(nl) > digit(nr): nl=4, nr=3
        EXPECT_FALSE(PID::isPentaquark(9343211));
        // all checks pass: n=9,nr=4,nl=4,nq1=3,nq2=2,nq3=2,nj=1
        EXPECT_TRUE(PID::isPentaquark(9443221));
    }

} // namespace Hammer
