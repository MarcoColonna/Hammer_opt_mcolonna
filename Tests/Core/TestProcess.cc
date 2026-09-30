///
/// @file  TestProcess.cc
/// @brief Tests for Process
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "Hammer/Tools/Utils.hh"
#include "Hammer/Process.hh"
#include "Hammer/Tools/ParticleUtils.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Tools/SettingsConsumer.hh"
#include "Hammer/Tools/HammerSerial.hh"
#include "gtest/gtest.h"

using namespace std;

namespace Hammer {


    TEST(ProcessTest, Construction) {
        Process proc;

        Particle pBmes{{5280., 0., 0., 0.}, 511};
        Particle pDmes{{2381.3371929065643, 0., 0., -1475.671313782346}, -411};
        Particle pNuTau{{863.8849328596731, -553.6990600712329, 92.91177555543321, 656.5682980934517}, 16};
        Particle pTau{{2034.777874233763, 553.6990600712329, -92.91177555543321, 819.1030156888944}, -15};
        Particle pNuBarTau{{751.3592069889307, -309.29260061671044, -268.9081311284792, 629.7357875679096}, -16};
        Particle pNuMuon{{532.7499030827818, 178.97752649767912, 105.73956803205624, 490.5187539695908}, 14};
        Particle pMuon{{750.6687641620506, 684.0141341902644, 70.25678754098979, -301.15152584860607}, -13};

        auto idxBmes = proc.addParticle(pBmes);
        auto idxDmes = proc.addParticle(pDmes);
        auto idxTau = proc.addParticle(pTau);
        auto idxNuTau = proc.addParticle(pNuTau);
        auto idxNuBarTau = proc.addParticle(pNuBarTau);
        auto idxNuMuon = proc.addParticle(pNuMuon);
        auto idxMuon = proc.addParticle(pMuon);

        // out of canonical order
        proc.addVertex(idxBmes, {idxTau, idxNuTau, idxDmes});
        proc.addVertex(idxTau, {idxNuMuon, idxNuBarTau, idxMuon});

        EXPECT_EQ(proc.getParentId(idxDmes), idxBmes);
        EXPECT_EQ(proc.getSiblingsIds(idxDmes).size(), 2);

        EXPECT_EQ(proc.getParticle(idxTau).pdgId(), pTau.pdgId());
        // parent
        EXPECT_EQ(proc.getParent(idxTau).pdgId(), pBmes.pdgId());
        // daughters
        EXPECT_EQ(proc.getDaughters(idxTau, true)[0].pdgId(), pNuBarTau.pdgId());
        EXPECT_EQ(proc.getDaughters(idxTau, true)[1].pdgId(), pNuMuon.pdgId());
        EXPECT_EQ(proc.getDaughters(idxTau, true)[2].pdgId(), pMuon.pdgId());
        EXPECT_EQ(proc.getDaughters(idxTau, false)[1].pdgId(), pNuBarTau.pdgId());
        EXPECT_EQ(proc.getDaughters(idxTau, false)[0].pdgId(), pNuMuon.pdgId());
        EXPECT_EQ(proc.getDaughters(idxTau, false)[2].pdgId(), pMuon.pdgId());
        // siblings
        EXPECT_EQ(proc.getSiblings(idxTau, true)[0].pdgId(), pDmes.pdgId());
        EXPECT_EQ(proc.getSiblings(idxTau, true)[1].pdgId(), pNuTau.pdgId());
        EXPECT_EQ(proc.getSiblings(idxTau, false)[1].pdgId(), pDmes.pdgId());
        EXPECT_EQ(proc.getSiblings(idxTau, false)[0].pdgId(), pNuTau.pdgId());
    }

    TEST(ProcessTest, Hash) {
        Process proc;

        Particle pBmes{{5280., 0., 0., 0.}, 511};
        Particle pDmes{{2381.3371929065643, 0., 0., -1475.671313782346}, -411};
        Particle pNuTau{{863.8849328596731, -553.6990600712329, 92.91177555543321, 656.5682980934517}, 16};
        Particle pTau{{2034.777874233763, 553.6990600712329, -92.91177555543321, 819.1030156888944}, -15};
        Particle pNuBarTau{{751.3592069889307, -309.29260061671044, -268.9081311284792, 629.7357875679096}, -16};
        Particle pNuMuon{{532.7499030827818, 178.97752649767912, 105.73956803205624, 490.5187539695908}, 14};
        Particle pMuon{{750.6687641620506, 684.0141341902644, 70.25678754098979, -301.15152584860607}, -13};

        auto idxBmes = proc.addParticle(pBmes);
        auto idxDmes = proc.addParticle(pDmes);
        auto idxTau = proc.addParticle(pTau);
        auto idxNuTau = proc.addParticle(pNuTau);
        auto idxNuBarTau = proc.addParticle(pNuBarTau);
        auto idxNuMuon = proc.addParticle(pNuMuon);
        auto idxMuon = proc.addParticle(pMuon);

        proc.addVertex(idxBmes, {idxDmes, idxNuTau, idxTau});
        proc.addVertex(idxTau, {idxNuBarTau, idxNuMuon, idxMuon});

        proc.initialize();

        // Explicit hashes
        auto tmpB = combineDaughters({-411, 16, -15}, {});
        HashId idB = processID(511, tmpB);
        auto tmpT = combineDaughters({-16, 14, -13}, {});
        HashId idT = processID(-15, tmpT);
        HashId seed = 0ul;
        combine_hash(seed, idB);
        combine_hash(seed, idT);

        EXPECT_EQ(proc.getId(), seed);

        auto TauV = proc.getParticlesByVertex(-15, {-16, -13, 14});
        EXPECT_EQ(proc.getParticle(idxTau).pdgId(), get<0>(TauV).pdgId());
        EXPECT_EQ(proc.getDaughters(idxTau, true)[0].pdgId(), (get<1>(TauV))[0].pdgId());
        EXPECT_EQ(proc.getDaughters(idxTau, true)[1].pdgId(), (get<1>(TauV))[1].pdgId());
        EXPECT_EQ(proc.getDaughters(idxTau, true)[2].pdgId(), (get<1>(TauV))[2].pdgId());

        auto vertexId = proc.getVertexId("TauMuNuNu");
        EXPECT_EQ(idT, vertexId);
    }

    TEST(ProcessTest, Iterators) {
        Process proc;

        Particle pBmes{{5280., 0., 0., 0.}, 511};
        Particle pDmes{{2381.3371929065643, 0., 0., -1475.671313782346}, -411};
        Particle pNuTau{{863.8849328596731, -553.6990600712329, 92.91177555543321, 656.5682980934517}, 16};
        Particle pTau{{2034.777874233763, 553.6990600712329, -92.91177555543321, 819.1030156888944}, -15};
        Particle pNuBarTau{{751.3592069889307, -309.29260061671044, -268.9081311284792, 629.7357875679096}, -16};
        Particle pNuMuon{{532.7499030827818, 178.97752649767912, 105.73956803205624, 490.5187539695908}, 14};
        Particle pMuon{{750.6687641620506, 684.0141341902644, 70.25678754098979, -301.15152584860607}, -13};

        auto idxBmes = proc.addParticle(pBmes);
        auto idxDmes = proc.addParticle(pDmes);
        auto idxTau = proc.addParticle(pTau);
        auto idxNuTau = proc.addParticle(pNuTau);
        auto idxNuBarTau = proc.addParticle(pNuBarTau);
        auto idxNuMuon = proc.addParticle(pNuMuon);
        auto idxMuon = proc.addParticle(pMuon);

        proc.addVertex(idxBmes, {idxDmes, idxNuTau, idxTau});
        proc.addVertex(idxTau, {idxNuBarTau, idxNuMuon, idxMuon});
    }


    TEST(ProcessTest, PHOTOSTree) {
        SettingsHandler sh;
        sh.addSetting<string>("Hammer", "Units", "GeV");
        Process procR;
        procR.setSettingsHandler(sh);

        // Test rad photon expander on tree
        // Process B^0bar to D^*+ tau^- nubar gamma, D^*+ to D^0 pi^+, tau^- to mu^- nu nu
        Particle pBmes{
            {5.2969140979437572e+00, 3.2824493742660787e-01, -2.3429865618090287e-01, 1.4385195640293111e-01}, -511};
        Particle pDstar{
            {2.3470774909812628e+00, 8.0970717778699397e-01, 7.6332563305857859e-01, 4.7945687459445407e-01}, 413};
        Particle pTau{
            {2.0795408878883048e+00, -8.0737090039644110e-01, -6.1305090764429815e-01, 3.7379380121775035e-01}, 15};
        Particle pNuTauBar{
            {8.7026750024654076e-01, 3.2592376941100554e-01, -3.8455135916505445e-01, -7.0940783099892080e-01}, -16};
        Particle pGamma{
            {2.8218818201336615e-05, -1.5109375535821102e-05, -2.2022429710929912e-05, 9.1115893907028467e-06}, 22};
        Particle pD{{2.1652224797353812e+00, 7.5505747015863178e-01, 6.8723913013095317e-01, 4.1006189279903876e-01},
                    421};
        Particle pPi{{1.8185501124588163e-01, 5.4649707628362069e-02, 7.6086502927625407e-02, 6.9394981795415422e-02},
                     211};
        Particle pMu{{1.1411538114089392e+00, -5.8718874409586341e-01, -4.5464791308157315e-01, 8.5998429292620537e-01},
                     13};
        Particle pNuTau{
            {3.7918042426691972e-01, 2.7540552999632262e-01, -1.2352239877392272e-01, -2.2950338820118057e-01}, 16};
        Particle pNuMuBar{
            {5.5920665221244559e-01, -4.9558768629690025e-01, -3.4880595788802433e-02, -2.5668710350727436e-01}, -14};

        auto idxB = procR.addParticle(pBmes);
        auto idxDstar = procR.addParticle(pDstar);
        auto idxTau = procR.addParticle(pTau);
        auto idxNuTauBar = procR.addParticle(pNuTauBar);
        auto idxGamma = procR.addParticle(pGamma);
        auto idxD = procR.addParticle(pD);
        auto idxPi = procR.addParticle(pPi);
        auto idxMu = procR.addParticle(pMu);
        auto idxNuTau = procR.addParticle(pNuTau);
        auto idxNuMuBar = procR.addParticle(pNuMuBar);

        procR.addVertex(idxB, {idxDstar, idxNuTauBar, idxTau, idxGamma});
        procR.addVertex(idxDstar, {idxD, idxPi});
        procR.addVertex(idxTau, {idxNuTau, idxNuMuBar, idxMu});

        procR.initialize();

        auto daughters = procR.getDaughtersIds(idxB);
        EXPECT_EQ(daughters.size(), 3);
        vector<ParticleIndex> daughtersExp = {idxDstar, idxNuTauBar, idxTau};
        for (size_t idx1 = 0; idx1 < daughters.size(); ++idx1) {
            EXPECT_EQ(daughters[idx1], daughtersExp[idx1]);
        }

        auto pTot = procR.getParticle(idxB).p() - procR.getParticle(idxDstar).p() - procR.getParticle(idxTau).p() -
                    procR.getParticle(idxNuTauBar).p();
        auto pTotD = procR.getParticle(idxDstar).p() - procR.getParticle(idxD).p() - procR.getParticle(idxPi).p();
        auto pTotT = procR.getParticle(idxTau).p() - procR.getParticle(idxNuTau).p() -
                     procR.getParticle(idxNuMuBar).p() - procR.getParticle(idxMu).p();

        EXPECT_LT(pTot.E(), 1e-8);
        EXPECT_LT(pTot.px(), 1e-8);
        EXPECT_LT(pTot.py(), 1e-8);
        EXPECT_LT(pTot.pz(), 1e-8);

        EXPECT_LT(pTotD.E(), 1e-8);
        EXPECT_LT(pTotD.px(), 1e-8);
        EXPECT_LT(pTotD.py(), 1e-8);
        EXPECT_LT(pTotD.pz(), 1e-8);

        EXPECT_LT(pTotT.E(), 1e-8);
        EXPECT_LT(pTotT.px(), 1e-8);
        EXPECT_LT(pTotT.py(), 1e-8);
        EXPECT_LT(pTotT.pz(), 1e-8);
    }

    TEST(ProcessTest, Serialization) {
        SettingsHandler sh;
        sh.addSetting<string>("Hammer", "Units", "MeV");
        Process proc;
        proc.setSettingsHandler(sh);

        Particle pBmes{{5280., 0., 0., 0.}, 511};
        Particle pDmes{{2381.3371929065643, 0., 0., -1475.671313782346}, -411};
        Particle pNuTau{{863.8849328596731, -553.6990600712329, 92.91177555543321, 656.5682980934517}, 16};
        Particle pTau{{2034.777874233763, 553.6990600712329, -92.91177555543321, 819.1030156888944}, -15};
        Particle pNuBarTau{{751.3592069889307, -309.29260061671044, -268.9081311284792, 629.7357875679096}, -16};
        Particle pNuMuon{{532.7499030827818, 178.97752649767912, 105.73956803205624, 490.5187539695908}, 14};
        Particle pMuon{{750.6687641620506, 684.0141341902644, 70.25678754098979, -301.15152584860607}, -13};

        auto idxBmes = proc.addParticle(pBmes);
        auto idxDmes = proc.addParticle(pDmes);
        auto idxTau = proc.addParticle(pTau);
        auto idxNuTau = proc.addParticle(pNuTau);
        auto idxNuBarTau = proc.addParticle(pNuBarTau);
        auto idxNuMuon = proc.addParticle(pNuMuon);
        auto idxMuon = proc.addParticle(pMuon);

        proc.addVertex(idxBmes, {idxDmes, idxNuTau, idxTau});
        proc.addVertex(idxTau, {idxNuBarTau, idxNuMuon, idxMuon});

        proc.initialize();

        auto builder = make_unique<flatbuffers::FlatBufferBuilder>();
        flatbuffers::Offset<Serial::FBProcIDs> msg;
        builder->Clear();
        proc.write(builder.get(), &msg);
        builder->Finish(msg);
        uint8_t* buf = builder->GetBufferPointer();
        // unsigned int size = builder->GetSize();
        const auto* ser1 = flatbuffers::GetRoot<Serial::FBProcIDs>(buf);
        EXPECT_NE(ser1, nullptr);
        if (ser1 != nullptr) {
            Process proc2;
            proc2.initialize();
            proc2.read(ser1);
        }
    }

} // namespace Hammer
