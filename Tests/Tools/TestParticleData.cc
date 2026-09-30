///
/// @file  TestParticleData.cc
/// @brief Tests for ParticleData
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#define PARTICLEDATA_FRIENDS FRIEND_TEST(ParticleDataTest, CheckValidSignatureDirect)

#include "Hammer/Tools/Utils.hh"
#include "Hammer/Tools/ParticleData.hh"
#include "Hammer/Tools/Pdg.hh"

using namespace std;
namespace Hammer {

    class TestableParticleData : public ParticleData {
    public:

        TestableParticleData() = default;
        void eval(const Particle& /*parent*/, const ParticleList& /*daughters*/,
                  const ParticleList& /*references*/) override {
        }
        using ParticleData::_signatureIndex;
        using ParticleData::_signatures;
        using ParticleData::addProcessSignature;
    };

    TEST(ParticleDataTest, SetSignature) {
        TestableParticleData p1;
        p1.addProcessSignature(PID::MUON, {PID::NU_MU, PID::NU_EBAR, PID::ELECTRON});
        TestableParticleData p2;
        p2.addProcessSignature(PID::K0L, {PID::PI0, PID::PI0}, {PID::GAMMA, PID::GAMMA, PID::EMINUS, PID::EPLUS});
        EXPECT_EQ(p1._signatures[0].parent, PID::MUON);
        EXPECT_EQ(p1._signatures[0].daughters.size(), 3);
        EXPECT_EQ(p2._signatures[0].parent, PID::K0L);
        EXPECT_EQ(p2._signatures[0].daughters.size(), 2 + 4);
        EXPECT_EQ(p1._signatures[0].daughters[0], PID::NU_MU);
        EXPECT_EQ(p1._signatures[0].daughters[2], PID::ELECTRON);
        EXPECT_EQ(p2._signatures[0].daughters[0], PID::PI0);
        EXPECT_EQ(p2._signatures[0].daughters[1], PID::PI0);
        EXPECT_EQ(p2._signatures[0].daughters[2], PID::GAMMA);
        EXPECT_EQ(p2._signatures[0].daughters[5], PID::EPLUS);
    }

    TEST(ParticleDataTest, CalcId) {
        TestableParticleData p1;
        p1.addProcessSignature(PID::MUON, {PID::NU_MU, PID::NU_EBAR, PID::ELECTRON});
        HashId seed = 0ul;
        combine_hash(seed, static_cast<int>(PID::MUON));
        combine_hash(seed, static_cast<int>(PID::NU_MU));
        combine_hash(seed, static_cast<int>(PID::NU_EBAR));
        combine_hash(seed, static_cast<int>(PID::ELECTRON));
        EXPECT_EQ(p1._signatures[0].id, seed);
    }

    TEST(ParticleDataTest, HadronicId) {
        TestableParticleData p1;
        p1.addProcessSignature(PID::BZERO, {PID::DSTARMINUS, PID::NU_MU, PID::ANTIMUON}, {-PID::D0, PID::PIMINUS});
        HashId seed = 0ul;
        combine_hash(seed, static_cast<int>(PID::BZERO));
        combine_hash(seed, static_cast<int>(PID::DSTARMINUS));
        EXPECT_EQ(p1._signatures[0].hadronicId, seed);
    }

    TEST(ParticleDataTest, IdAndHadronicIdAccessors) {
        TestableParticleData p1;
        p1.addProcessSignature(PID::MUON, {PID::NU_MU, PID::NU_EBAR, PID::ELECTRON});
        EXPECT_EQ(p1.id(), p1._signatures[0].id);
        EXPECT_EQ(p1.hadronicId(), p1._signatures[0].hadronicId);
    }

    TEST(ParticleDataTest, MassesAccessor) {
        TestableParticleData p1;
        p1.addProcessSignature(PID::MUON, {PID::NU_MU, PID::NU_EBAR, PID::ELECTRON});
        const auto& m = p1.masses();
        EXPECT_EQ(m.size(), 4u);
        EXPECT_NEAR(m[0], PID::instance().getMass(PID::MUON), 1e-8);
    }

    TEST(ParticleDataTest, NumSignatures) {
        TestableParticleData p1;
        EXPECT_EQ(p1.numSignatures(), 0u);
        p1.addProcessSignature(PID::MUON, {PID::NU_MU, PID::NU_EBAR, PID::ELECTRON});
        EXPECT_EQ(p1.numSignatures(), 1u);
        p1.addProcessSignature(PID::TAU, {PID::NU_TAU, PID::NU_EBAR, PID::ELECTRON});
        EXPECT_EQ(p1.numSignatures(), 2u);
    }

    TEST(ParticleDataTest, SetSignatureIndex) {
        TestableParticleData p1;
        p1.addProcessSignature(PID::MUON, {PID::NU_MU, PID::NU_EBAR, PID::ELECTRON});
        p1.addProcessSignature(PID::TAU, {PID::NU_TAU, PID::NU_EBAR, PID::ELECTRON});
        EXPECT_TRUE(p1.setSignatureIndex(1));
        EXPECT_EQ(p1._signatureIndex, 1u);
        EXPECT_FALSE(p1.setSignatureIndex(100));
        EXPECT_EQ(p1._signatureIndex, 0u);
    }

    TEST(ParticleDataTest, MultipleSignaturesMasses) {
        TestableParticleData p1;
        p1.addProcessSignature(PID::TAU, {PID::NU_TAU, PID::NU_EBAR, PID::ELECTRON});
        const auto& m = p1.masses();
        EXPECT_EQ(m.size(), 4u);
        EXPECT_NEAR(m[0], PID::instance().getMass(PID::TAU), 1e-4);
        EXPECT_NEAR(m[3], PID::instance().getMass(PID::ELECTRON), 1e-7);
    }

    TEST(ParticleDataTest, AddSignatureWithGranddaughtersNegativeParent) {
        TestableParticleData p1;
        p1.addProcessSignature(-PID::BZERO, {-PID::DSTARMINUS, PID::NU_MUBAR, PID::MUON}, {PID::D0, PID::PIPLUS});
        EXPECT_EQ(p1._signatures[0].parent, -PID::BZERO);
        EXPECT_EQ(p1._signatures[0].daughters.size(), 5u);
    }

    TEST(ParticleDataTest, CheckValidSignatureDirect) {
        TestableParticleData p;
        // valid: |NU_MU|=14 > |NU_EBAR|=12 > |ELECTRON|=11 matches combineDaughters order
        EXPECT_TRUE(p.checkValidSignature(PID::MUON, {PID::NU_MU, PID::NU_EBAR, PID::ELECTRON}, {}));
        // invalid: wrong order
        EXPECT_FALSE(p.checkValidSignature(PID::MUON, {PID::ELECTRON, PID::NU_MU, PID::NU_EBAR}, {}));
        // valid with granddaughters and negative parent (covers flipSign branch)
        EXPECT_TRUE(
            p.checkValidSignature(-PID::BZERO, {-PID::DSTARMINUS, PID::NU_MUBAR, PID::MUON}, {PID::D0, PID::PIPLUS}));
        // valid with granddaughters and positive parent (covers id-pass-through branch)
        EXPECT_TRUE(
            p.checkValidSignature(PID::BZERO, {PID::DSTARMINUS, PID::NU_MU, PID::ANTIMUON}, {-PID::D0, PID::PIMINUS}));
    }

} // namespace Hammer
