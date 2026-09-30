///
/// @file  TestHammer.cc
/// @brief Tests for Hammer
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "Hammer/Hammer.hh"
#include "Hammer/Tools/ParticleUtils.hh"
#include "Hammer/Process.hh"
#include "Hammer/Math/MultiDimensional.fhh"
#include "gtest/gtest.h"

using namespace std;

namespace Hammer {

    TEST(HammerTest, Assignment) {

        Hammer ham;

        // Create included processes
        ham.includeDecay(vector<string>{"BD*TauNu", "D*DPi"});
        ham.includeDecay("BDEllNu");
        ham.addFFScheme("Scheme1", {{"BD", "BGL"}, {"BD*", "BLPR"}});
        ham.addFFScheme("Scheme2", {{"BD", "BLPR"}, {"BD*", "ISGW2"}});
        ham.setFFInputScheme({{"BD", "ISGW2"}, {"BD*", "ISGW2"}});
        IndexList bins = {10};
        ham.addHistogram("EPiP", bins, false);
        ham.initRun();
        ham.saveOptionCard("../../OptsA.yml", false);

        // Create D* process. Need to put in real momenta for an amplitude calculation. Test conj signatures with sgn
        Process procDs;
        int sgn = -1;

        Particle pB({5.280000000000000, 0., 0., 0.}, 511 * sgn);
        Particle pDs({2.660051210321970, 1.389723696643482, 0.4037520573116902, 0.9702704595830642}, -413 * sgn);
        Particle pTau({2.434889134265837, -1.420021405904506, -0.3418057134095529, -0.7985387637550296}, -15 * sgn);
        Particle pNuTau({0.18506, 0.03029770926102430, -0.06194634390213728, -0.1717316958280348}, 16 * sgn);
        Particle pD({2.501949325562050, 1.340019413229273, 0.3772064686434118, 0.9103036196031724}, -411 * sgn);
        Particle pPiZ({0.1581018847599200, 0.04970428341420862, 0.02654558866827830, 0.05996683997989181}, 111);
        Particle pPiP({1.681083587524002, -1.196486326353599, -1.024596160237620, -0.5703895226227106}, 211 * sgn);
        Particle pNuTauBar({0.75381, -0.2235350795509063, 0.6827904468280672, -0.2281492411323186}, -16 * sgn);

        auto idxB = procDs.addParticle(pB);
        auto idxDs = procDs.addParticle(pDs);
        auto idxTau = procDs.addParticle(pTau);
        auto idxNuTau = procDs.addParticle(pNuTau);
        auto idxD = procDs.addParticle(pD);
        auto idxPiZ = procDs.addParticle(pPiZ);
        auto idxPiP = procDs.addParticle(pPiP);
        auto idxNuTauBar = procDs.addParticle(pNuTauBar);

        procDs.addVertex(idxB, {idxDs, idxTau, idxNuTau});
        procDs.addVertex(idxDs, {idxD, idxPiZ});
        procDs.addVertex(idxTau, {idxPiP, idxNuTauBar});

        // Create D process
        Process procD;
        Particle qB({5.280000000000000, 0., 0., 0.}, 511);
        Particle qD({2.336261569054783, 0., 0., -1.401769281673099}, -411);
        Particle qMu({0.7919605358187830, 0.2826874295413168, 0.04579716682155486, -0.7308672430645772}, -13);
        Particle qNuMu({2.15178, -0.2826874295413168, -0.04579716682155486, 2.132636524737676}, 14);

        auto idx2B = procD.addParticle(qB);
        auto idx2D = procD.addParticle(qD);
        auto idx2Mu = procD.addParticle(qMu);
        auto idx2NuMu = procD.addParticle(qNuMu);

        procD.addVertex(idx2B, {idx2D, idx2Mu, idx2NuMu});

        // add Processes
        ham.addProcess(procDs);
        ham.addProcess(procD);

        // find the bin
        auto EPiPbin = static_cast<IndexType>(floor(10. * (pPiP.p().E() - 1.12465) / (1.80204 - 1.12465)));
        ham.setEventHistogramBin("EPiP", {EPiPbin});

        // Execute evaluations
        ham.processEvent();

        // Outputs
        auto histo = ham.getHistogram("EPiP", "Scheme1");
        auto weights0 = ham.getWeights("Denominator");
        auto weights1 = ham.getWeights("Scheme1");
        auto weights2 = ham.getWeights("Scheme2");

        for (auto& elem : weights0) {
            cout << "Weight0: " << elem.first << " " << elem.second << '\n';
        }
        for (auto& elem : weights1) {
            cout << "Weight1: " << elem.first << " " << elem.second << '\n';
        }
        for (auto& elem : weights2) {
            cout << "Weight2: " << elem.first << " " << elem.second << '\n';
        }

        cout << "Histo EPiP: {"; // << histo.first << " {";
        for (size_t idx = 0; idx < histo.size() - 1; ++idx) {
            cout << histo[idx].sumWi << ", ";
        }
        cout << histo[histo.size() - 1].sumWi << "}" << '\n';

        ham.setWilsonCoefficients("BtoCTauNu", {{"SM", 2.0}});
        weights0 = ham.getWeights("Denominator");
        weights1 = ham.getWeights("Scheme1");
        weights2 = ham.getWeights("Scheme2");
        for (auto& elem : weights0) {
            cout << "Weight0: " << elem.first << " " << elem.second << '\n';
            EXPECT_NEAR(elem.second, 1., 1e-4);
        }
        for (auto& elem : weights1) {
            cout << "Weight1: " << elem.first << " " << elem.second << '\n';
        }
        for (auto& elem : weights2) {
            cout << "Weight2: " << elem.first << " " << elem.second << '\n';
        }

        // Denominator scheme uses the den WCs, as in getExternalVectors
        auto weights0pre = ham.getWeight("Denominator");
        ham.setWilsonCoefficients("BtoCTauNu", {{"SM", 2.0}}, WTerm::DENOMINATOR);
        auto weights0post = ham.getWeight("Denominator");
        EXPECT_NEAR(weights0post / weights0pre, 4., 1e-4);

        // getRate
        auto vertexIdDs = procDs.getVertexId("BD*TauNu");
        auto vertexIdD = procD.getVertexId("BDMuNu");

        ham.setWilsonCoefficients("BtoCTauNu", {{"SM", 1.0}});
        auto rateDs0 = ham.getRate("B0barD*+TauNu", "Scheme1");
        auto rateDs1 = ham.getRate(vertexIdDs, "Scheme1");
        auto rateD1 = ham.getRate(vertexIdD, "Scheme1");
        auto rateDsDen1 = ham.getDenominatorRate(vertexIdDs);
        ham.setWilsonCoefficients("BtoCTauNu", {{"SM", 2.0}});
        auto rateDs2 = ham.getRate(511 * sgn, {-413 * sgn, -15 * sgn, 16 * sgn}, "Scheme1");
        auto rateD2 = ham.getRate(511, {-411, -13, 14}, "Scheme1");
        auto rateDsDen2 = ham.getDenominatorRate(vertexIdDs);

        EXPECT_NEAR(rateDs0 / rateDs1, 1., 1e-4);
        EXPECT_NEAR(rateDs1 / rateDs2, 0.25, 1e-4);
        EXPECT_NEAR(rateD1 / rateD2, 1., 1e-4);
        EXPECT_NEAR(rateDsDen1 / rateDsDen2, 1., 1e-4);

        auto vertexIdDPi = procDs.getVertexId("D*DPi");
        EXPECT_NEAR(ham.getRate(vertexIdDPi, "Scheme1") / ham.getRate(vertexIdDPi, "Denominator"), 1., 1e-4);

        // Test Units
        ham.setUnits("MeV");
        ham.setWilsonCoefficients("BtoCTauNu", {{"SM", 1.0}});
        auto rateDs0MeV = ham.getRate("B0barD*+TauNu", "Scheme1");
        EXPECT_NEAR(rateDs0MeV / (1e3 * rateDs0), 1., 1e-4);
    }

    TEST(HammerTest, Units) {
        Hammer ham;

        ham.setOptions("ProcessCalc: {Rates: false}");

        // Create included processes
        ham.includeDecay("BD*MuNu");
        ham.addFFScheme("Scheme1", {{"BD*", "BGL"}});
        ham.setFFInputScheme({{"BD*", "CLN"}});
        ham.setUnits("MeV");

        ham.initRun();

        // Create D* process.
        Process procDs;

        Particle pB({216824.0, 10723.3, 12380.7, 216140.0}, -511);
        Particle pDs({86708.9, 4436.91, 4908.47, 86432.7}, 413);
        Particle pNuMu({62072.52, 3730.09, 5010.67, 61757.4}, -14);
        Particle pMu({68042.3, 2556.34, 2461.54, 67949.6}, 13);

        Particle pD({79758.8, 4105.85, 4494.7, 79504.2}, 421);
        Particle pPiD;
        pPiD.setMomentum(pDs.p() - pD.p());
        pPiD.setPdgId(211);

        auto idxB = procDs.addParticle(pB);
        auto idxDs = procDs.addParticle(pDs);
        auto idxNuE = procDs.addParticle(pNuMu);
        auto idxE = procDs.addParticle(pMu);

        auto idxD = procDs.addParticle(pD);
        auto idxPiD = procDs.addParticle(pPiD);


        procDs.addVertex(idxDs, {idxD, idxPiD});
        procDs.addVertex(idxB, {idxDs, idxNuE, idxE});

        // add Processes
        ham.addProcess(procDs);

        // Execute evaluations
        ham.processEvent();

        // Outputs
        auto weight0 = ham.getWeight("Denominator");
        auto weight1 = ham.getWeight("Scheme1");
        cout << "Trivial: " << weight0 << "\n" << "CLN->BGL: " << weight1 << '\n';

        EXPECT_NEAR(weight0, 1., 1e-10);
        EXPECT_NEAR(weight1 / 0.984506, 1., 1e-5);
    }

    TEST(HammerTest, Tau3Pi) {

        Hammer ham;

        ham.setOptions("ProcessCalc: {Rates: false}");

        // Create included processes
        ham.includeDecay(vector<string>{"BDTauNu"});
        ham.addFFScheme("Scheme", {{"BD", "BGL"}, {"BD*", "BLPR"}, {"TauPiPiPi", "RCT"}});
        ham.setFFInputScheme({{"BD", "ISGW2"}, {"BD*", "ISGW2"}, {"TauPiPiPi", "RCT"}});
        ham.initRun();

        // Create process
        Process proc;

        // Particle decs
        Particle pBmes({7.27993, 4.69979, 1.60853, 0.666025}, 511);
        Particle pDmes({5.42945, 4.71853, 1.60413, 1.07146}, -411);
        Particle kNuTau({0.0251286, -0.018733, 0.00440374, 0.0161594}, 16);
        Particle pTau({1.82536, 0.0, 0.0, -0.421598}, -15);
        Particle kNuBarTau({0.421598, 0.0, 0.0, -0.421598}, -16);
        Particle pPiPlus1({0.554558, -0.525713, 0.102531, 0.03643}, 211);
        Particle pPiPlus2({0.360269, 0.169702, -0.0979545, 0.268475}, 211);
        Particle pPiMinus({0.48893, 0.356011, -0.00457666, -0.304905}, -211);

        auto idxBmes = proc.addParticle(pBmes);
        auto idxDmes = proc.addParticle(pDmes);
        auto idxNuTau = proc.addParticle(kNuTau);
        auto idxTau = proc.addParticle(pTau);
        auto idxNuBarTau = proc.addParticle(kNuBarTau);
        auto idxPiPlus1 = proc.addParticle(pPiPlus1);
        auto idxPiPlus2 = proc.addParticle(pPiPlus2);
        auto idxPiMinus = proc.addParticle(pPiMinus);

        proc.addVertex(idxBmes, {idxDmes, idxNuTau, idxTau});
        proc.addVertex(idxTau, {idxNuBarTau, idxPiPlus1, idxPiPlus2, idxPiMinus});

        ham.addProcess(proc);

        ham.processEvent();

        const std::vector<HashId>& processes{proc.getId()};
        auto weight = ham.getWeight("Scheme", processes);
        cout << "Tau Weight: " << weight << '\n';

        ham.setWilsonCoefficients("BtoCTauNu", {{"SM", 2.0}});
        const std::vector<HashId>& processes2{proc.getId()};
        auto weight2 = ham.getWeight("Scheme", processes2);

        EXPECT_NEAR(weight2 / weight, 4., 1e-4);
    }

    TEST(HammerTest, FFEigenvectors) {

        Hammer ham;

        // Create included processes
        ham.includeDecay("BDEllNu");
        ham.addFFScheme("Scheme1", {{"BD", "BGL_1"}});
        ham.addFFScheme("Scheme2", {{"BD", "BGL_2"}});
        ham.addFFScheme("SchemeVar", {{"BD", "BGLVar"}});
        ham.setFFInputScheme({{"BD", "ISGW2"}});

        ham.initRun();

        // Change BGL instances independently
        ham.setOptions("BtoDBGL_1: {ChiT: 0.01, ChiL: 0.002}");
        ham.setOptions("BtoDBGL_2: {ChiT: 0.02, ChiL: 0.004}");
        ham.setOptions("BtoDBGLVar: {ChiT: 0.01, ChiL: 0.002}");
        // Check eigenvectors. Set central values to uniform value
        ham.setOptions("BtoDBGLVar: {a0: [1., 1., 1., 1.], ap: [1., 1., 1., 1.]}");
        //        ham.saveOptionCard("../../OptsD.yml", false);

        Process procD;
        Particle qB({5.280000000000000, 0., 0., 0.}, 511);
        Particle qD({2.336261569054783, 0., 0., -1.401769281673099}, -411);
        Particle qMu({0.7919605358187830, 0.2826874295413168, 0.04579716682155486, -0.7308672430645772}, -13);
        Particle qNuMu({2.15178, -0.2826874295413168, -0.04579716682155486, 2.132636524737676}, 14);

        auto idxB = procD.addParticle(qB);
        auto idxD = procD.addParticle(qD);
        auto idxMu = procD.addParticle(qMu);
        auto idxNuMu = procD.addParticle(qNuMu);

        procD.addVertex(idxB, {idxD, idxMu, idxNuMu});

        // add Processes
        ham.addProcess(procD);

        ham.processEvent();

        // Check central values match
        auto weight1 = ham.getWeight("Scheme1");
        auto weight2 = ham.getWeight("Scheme2");
        auto weightVar = ham.getWeight("SchemeVar");

        cout << "BGL_1 weight: " << weight1 << '\n';
        cout << "BGL_2 weight: " << weight2 << '\n';
        EXPECT_NEAR(weight2 / weight1, 2., 1e-4);

        cout << "BGLVar central weight: " << weightVar << '\n';
        // Now effectively further double parameter values
        ham.setFFEigenvectors("BtoD", "BGLVar",
                              {{"delta_ap0", 1.},
                               {"delta_ap1", 1.},
                               {"delta_ap2", 1.},
                               {"delta_ap3", 1.},
                               {"delta_a00", 1.},
                               {"delta_a01", 1.},
                               {"delta_a02", 1.},
                               {"delta_a03", 1.}});
        auto weightPost = ham.getWeight("SchemeVar");

        cout << "BGLVar a's = 1 weight: " << weightVar << '\n';
        cout << "BGLVar a's = 2 weight: " << weightPost << '\n';
        EXPECT_NEAR(weightPost / weightVar, 4., 1e-4);
    }

    TEST(HammerTest, FFOptsBs) {

        Hammer ham;

        ham.setOptions("ProcessCalc: {Rates: false}");

        // Create included processes
        ham.includeDecay("BsDsEllNu");
        ham.includeDecay("BDTauNu");
        ham.addFFScheme("Scheme1", {{"BsDs", "BGL"}, {"BD", "BGL"}});
        ham.addFFScheme("Scheme2", {{"BsDs", "BGL_1"}, {"BD", "BGL_1"}});
        ham.addFFScheme("SchemeVar", {{"BsDs", "BGLVar"}, {"BD", "BGLVar"}});
        ham.setFFInputScheme({{"BD", "ISGW2"}, {"BsDs", "ISGW2"}});

        ham.initRun();

        // Change BGL instances independently
        ham.setOptions("BstoDsBGL: {ChiT: 0.01, ChiL: 0.002}");
        ham.setOptions("BtoDBGL: {ChiT: 0.01, ChiL: 0.002}");
        ham.setOptions("BtoDBGL_1: {ChiT: 0.03, ChiL: 0.006}");
        ham.setOptions("BstoDsBGL_1: {ChiT: 0.02, ChiL: 0.004}");
        ham.setOptions("BstoDsBGLVar: {ChiT: 0.01, ChiL: 0.002}");
        ham.setOptions("BtoDBGLVar: {ChiT: 0.01, ChiL: 0.002}");
        // Check eigenvectors. Set central values to uniform value
        ham.setOptions("BtoDBGLVar: {a0: [1., 1., 1., 1.], ap: [1., 1., 1., 1.]}");
        ham.setOptions("BstoDsBGLVar: {a0: [1., 1., 1., 1.], ap: [1., 1., 1., 1.]}");

        Process procDs;
        Particle qBs({5.280000000000000, 0., 0., 0.}, 531);
        Particle qDs({2.336261569054783, 0., 0., -1.401769281673099}, -431);
        Particle qMu({0.7919605358187830, 0.2826874295413168, 0.04579716682155486, -0.7308672430645772}, -13);
        Particle qNuMu({2.15178, -0.2826874295413168, -0.04579716682155486, 2.132636524737676}, 14);

        auto idxBs = procDs.addParticle(qBs);
        auto idxDs = procDs.addParticle(qDs);
        auto idxMu = procDs.addParticle(qMu);
        auto idxNuMu = procDs.addParticle(qNuMu);

        procDs.addVertex(idxBs, {idxDs, idxMu, idxNuMu});

        Process procD;
        Particle pB({5.280000000000000, 0., 0., 0.}, 511);
        Particle pD({2.660051210321970, 1.389723696643482, 0.4037520573116902, 0.9702704595830642}, -411);
        Particle pTau({2.434889134265837, -1.420021405904506, -0.3418057134095529, -0.7985387637550296}, -15);
        Particle pNuTau({0.18506, 0.03029770926102430, -0.06194634390213728, -0.1717316958280348}, 16);

        auto idxB = procD.addParticle(pB);
        auto idxD = procD.addParticle(pD);
        auto idxTau = procD.addParticle(pTau);
        auto idxNuTau = procD.addParticle(pNuTau);

        procD.addVertex(idxB, {idxD, idxTau, idxNuTau});

        // add Processes
        auto idD = ham.addProcess(procD);
        auto idDs = ham.addProcess(procDs);

        ham.processEvent();

        // Check central values match
        auto weights1 = ham.getWeights("Scheme1");
        auto weights2 = ham.getWeights("Scheme2");
        auto weightsVarPre = ham.getWeights("SchemeVar");

        cout << "B x Bs BGL weights: " << weights1[idD] << " " << weights1[idDs] << '\n';
        cout << "B x Bs BGL_1 weights: " << weights2[idD] << " " << weights2[idDs] << '\n';
        EXPECT_NEAR(weights2[idD] / weights1[idD], 3., 1e-4);
        EXPECT_NEAR(weights2[idDs] / weights1[idDs], 2., 1e-4);

        ham.setFFEigenvectors("BstoDs", "BGLVar", {1., 1., 1., 1., 1., 1., 1., 1.});
        auto weightsVarPost = ham.getWeights("SchemeVar");
        EXPECT_NEAR(weightsVarPost[idD] / weightsVarPre[idD], 1., 1e-4);
        EXPECT_NEAR(weightsVarPost[idDs] / weightsVarPre[idDs], 4., 1e-4);
    }

    TEST(HammerTest, FFEigenvectorsDstar) {

        Hammer ham;

        ham.setOptions("ProcessCalc: {Rates: false}");

        // Create included processes
        ham.includeDecay("BD*TauNu");
        ham.addFFScheme("Scheme1", {{"BD*", "BGL_1"}});
        ham.addFFScheme("Scheme2", {{"BD*", "BGL_2"}});
        ham.addFFScheme("SchemeVar", {{"BD*", "BGLVar"}});
        ham.setFFInputScheme({{"BD*", "ISGW2"}});

        ham.initRun();

        // Change BGL instances independently
        ham.setOptions("BtoD*BGL_1: {Chim: 0.01, Chip: 0.002, ChimL: 0.01}");
        ham.setOptions("BtoD*BGL_2: {Chim: 0.02, Chip: 0.004, ChimL: 0.02}");
        ham.setOptions("BtoD*BGLVar: {Chim: 0.01, Chip: 0.002, ChimL: 0.01}");
        // Check eigenvectors. Set central values to uniform value
        ham.setOptions("BtoD*BGLVar: {avec: [1., 1., 1.], bvec: [1., 1., 1.], cvec: [1., 1.], dvec: [1., 1.]}");
        //        ham.saveOptionCard("../../OptsDs.yml", false);

        Process procDs;
        Particle pB({5.280000000000000, 0., 0., 0.}, 511);
        Particle pDs({2.660051210321970, 1.389723696643482, 0.4037520573116902, 0.9702704595830642}, -413);
        Particle pTau({2.434889134265837, -1.420021405904506, -0.3418057134095529, -0.7985387637550296}, -15);
        Particle pNuTau({0.18506, 0.03029770926102430, -0.06194634390213728, -0.1717316958280348}, 16);

        auto idx2B = procDs.addParticle(pB);
        auto idx2Ds = procDs.addParticle(pDs);
        auto idx2Tau = procDs.addParticle(pTau);
        auto idx2NuTau = procDs.addParticle(pNuTau);

        procDs.addVertex(idx2B, {idx2Ds, idx2Tau, idx2NuTau});

        // add Processes
        ham.addProcess(procDs);

        ham.processEvent();

        // auto processed = std::chrono::system_clock::now();

        // Check central values match
        auto weight1 = ham.getWeight("Scheme1");
        auto weight2 = ham.getWeight("Scheme2");
        auto weightVar = ham.getWeight("SchemeVar");

        // auto wgt = std::chrono::system_clock::now();

        cout << "BGL_1 weight: " << weight1 << '\n';
        cout << "BGL_2 weight: " << weight2 << '\n';
        EXPECT_NEAR(weight2 / weight1, 2., 1e-4);

        cout << "BGLVar central weight: " << weightVar << '\n';
        // Now effectively further double parameter values
        ham.setFFEigenvectors("BtoD*", "BGLVar", {1., 1., 1., 1., 1., 1., 1., 1., 1., 1.});
        //        ham.setFFEigenvectors("BtoD*", "BGLVar",
        //        {{"delta_a0",1.},{"delta_a1",1.},{"delta_a2",1.},{"delta_b0",1.},{"delta_b1",1.},{"delta_b2",1.},{"delta_c1",1.},{"delta_c2",1.},{"delta_d0",1.},{"delta_d1",1.}});
        auto weightPost = ham.getWeight("SchemeVar");

        cout << "BGLVar a's = 1 weight: " << weightVar << '\n';
        cout << "BGLVar a's = 2 weight: " << weightPost << '\n';
        EXPECT_NEAR(weightPost / weightVar, 4., 1e-4);
    }

    TEST(HammerTest, FFEigenvectorsRename) {

        Hammer ham;

        ham.setOptions("ProcessCalc: {Rates: false}");

        // Create included processes
        ham.includeDecay("BD*TauNu");
        ham.addFFScheme("SchemeVar_1", {{"BD*", "BGLVar_1"}});
        ham.addFFScheme("SchemeVar_2", {{"BD*", "BGLVar_2"}});
        ham.setFFInputScheme({{"BD*", "ISGW2"}});
        cout << "A1" << '\n';

        ham.initRun();

        cout << "A2" << '\n';

        ham.renameFFEigenvectors("BtoD*", "BGLVar_1", {"dir_1", "dir_2", "dir_3"});
        ham.renameFFEigenvectors("BtoD*", "BGLVar_2", {"", "", "", "dir1", "dir2", "dir3"});

        cout << "A3" << '\n';
        // Change BGL instances independently
        ham.setOptions("BtoD*BGLVar_1: {Chim: 0.01, Chip: 0.002, ChimL: 0.01}");
        ham.setOptions("BtoD*BGLVar_2: {Chim: 0.02, Chip: 0.004, ChimL: 0.02}");
        // Check eigenvectors. Set central values to uniform value
        ham.setOptions("BtoD*BGLVar_1: {avec: [1., 1., 1.], bvec: [1., 1., 1.], cvec: [1., 1.], dvec: [1., 1.]}");
        ham.setOptions("BtoD*BGLVar_2: {avec: [1., 1., 1.], bvec: [1., 1., 1.], cvec: [1., 1.], dvec: [1., 1.]}");
        //        ham.saveOptionCard("../../OptsDs.yml", false);

        cout << "A4" << '\n';

        Process procDs;
        Particle pB({5.280000000000000, 0., 0., 0.}, 511);
        Particle pDs({2.660051210321970, 1.389723696643482, 0.4037520573116902, 0.9702704595830642}, -413);
        Particle pTau({2.434889134265837, -1.420021405904506, -0.3418057134095529, -0.7985387637550296}, -15);
        Particle pNuTau({0.18506, 0.03029770926102430, -0.06194634390213728, -0.1717316958280348}, 16);

        auto idx2B = procDs.addParticle(pB);
        auto idx2Ds = procDs.addParticle(pDs);
        auto idx2Tau = procDs.addParticle(pTau);
        auto idx2NuTau = procDs.addParticle(pNuTau);

        procDs.addVertex(idx2B, {idx2Ds, idx2Tau, idx2NuTau});

        // add Processes
        ham.addProcess(procDs);

        ham.processEvent();

        // auto processed = std::chrono::system_clock::now();

        // Check central values match
        auto weightVar1 = ham.getWeight("SchemeVar_1");
        auto weightVar2 = ham.getWeight("SchemeVar_2");

        // auto wgt = std::chrono::system_clock::now();

        cout << "BGLVar_1 weight: " << weightVar1 << '\n';
        cout << "BGLVar_2 weight: " << weightVar2 << '\n';
        EXPECT_NEAR(weightVar2 / weightVar1, 2., 1e-4);

        // Now effectively further double parameter values
        ham.setFFEigenvectors("BtoD*", "BGLVar_1", {1., 1., 1., 1., 1., 1., 1., 1., 1., 1.});
        auto weightPost1a = ham.getWeight("SchemeVar_1");
        ham.setFFEigenvectors("BtoD*", "BGLVar_1",
                              {{"dir_1", 1.},
                               {"dir_2", 1.},
                               {"dir_3", 1.},
                               {"delta_b0", 1.},
                               {"delta_b1", 1.},
                               {"delta_b2", 1.},
                               {"delta_c1", 1.},
                               {"delta_c2", 1.},
                               {"delta_d0", 1.},
                               {"delta_d1", 1.}});
        auto weightPost1b = ham.getWeight("SchemeVar_1");
        ham.setFFEigenvectors("BtoD*", "BGLVar_2", {1., 1., 1., 1., 1., 1., 1., 1., 1., 1.});
        auto weightPost2a = ham.getWeight("SchemeVar_2");
        ham.setFFEigenvectors("BtoD*", "BGLVar_2",
                              {{"delta_a0", 1.},
                               {"delta_a1", 1.},
                               {"delta_a2", 1.},
                               {"dir1", 1.},
                               {"dir2", 1.},
                               {"dir3", 1.},
                               {"delta_c1", 1.},
                               {"delta_c2", 1.},
                               {"delta_d0", 1.},
                               {"delta_d1", 1.}});
        auto weightPost2b = ham.getWeight("SchemeVar_2");

        cout << "BGLVar_1 a's = 1 weight: " << weightVar1 << '\n';
        cout << "BGLVar_1 a's = 2 weight: " << weightPost1a << " or " << weightPost1b << '\n';
        cout << "BGLVar_2 a's = 1 weight: " << weightVar2 << '\n';
        cout << "BGLVar_1 a's = 2 weight: " << weightPost2a << " or " << weightPost2b << '\n';
        EXPECT_NEAR(weightPost1a / weightPost1b, 1., 1e-4);
        EXPECT_NEAR(weightPost2a / weightPost2b, 1., 1e-4);
        EXPECT_NEAR(weightPost1a / weightVar1, 4., 1e-4);
        EXPECT_NEAR(weightPost2a / weightVar2, 4., 1e-4);
    }

    TEST(HammerTest, WCSpecialization) {
        Hammer ham;
        ham.setOptions("ProcessCalc: {Rates: false}");

        // Create included processes
        ham.includeDecay("BD*TauNu");
        ham.includeDecay("BD*MuNu");
        ham.addFFScheme("Scheme1", {{"BD*", "BLPR"}});
        ham.setFFInputScheme({{"BD*", "ISGW2"}});

        ham.addHistogram("pEllq2", {5, 5}, false, {{0., 2.5}, {0, 2.5}});

        ham.setUnits("GeV");

        ham.createWCSpecialization("JM", "BtoCTauNu",
                                   {"S1a", "R2"}); // default origin is SM or whatever WCs are set to when called
        ham.createWCSpecialization("JM", "BtoCMuNu", {"S1a"});
        ham.createWCSpecialization("RH", "BtoCTauNu", {}); // full specialization
        ham.createWCSpecialization("JC", "BtoCTauNu", {"R2"});
        ham.createWCSpecialization("TS", "BtoCTauNu", {"S1a", "R2"});

        // ham.setWCSpecializationOrigin("JM", "BtoCTauNu", {{"SM", 1}});
        ham.setWCSpecializationOrigin("RH", "BtoCTauNu",
                                      {{"SM", 1},
                                       {"S_qRlL", 1},
                                       {"T_qLlL", 0.125},
                                       {"T_qLqL", 0.125}}); // add a last wrong WC to check warn handling
        // ham.setWCSpecializationOrigin("JC", "BtoCTauNu", {{"SM", 1}});
        // ham.setWCSpecializationOrigin("TS", "BtoCTauNu", {{"SM", 1}});

        ham.setWCSpecializationBasis("JM", "BtoCTauNu",
                                     {{{"S_qRlL", -1}, {"S_qLlL", 1}}, {{"S_qRlL", 1}, {"T_qLlL", 0.125}}});
        ham.setWCSpecializationBasis("JM", "BtoCMuNu", {{{"S_qRlL", -1}, {"S_qLlL", 1}}});
        ham.setWCSpecializationBasis("JC", "BtoCTauNu", {{{"S_qRlL", 1}, {"T_qLlL", 0.125}}});
        ham.setWCSpecializationBasis("TS", "BtoCTauNu",
                                     {{{"S_qRlL", -1}, {"S_qLlL", 1}}, {{"S_qRlL", 1}, {"T_qLlL", 0.125}}});

        ham.applyWCSpecializationInWeights("JM"); // can work for weights after initRun, but not for histograms
        ham.applyWCSpecializationInWeights("RH");
        ham.applyWCSpecializationInWeights("JC");
        ham.applyWCSpecializationInHistogram("pEllq2", "TS"); // not in weights, should cause GeneralWeights to be kept
        ham.setOptions("Histos: {CalcGeneralHistogram_pEllq2: true}");

        ham.initRun();

        vector<pair<string, PdgId>> leps{make_pair("Tau", 15), make_pair("Mu", 13)};
        map<string, vector<complex<double>>> wcs{{"Tau", {1, 0}}, {"Mu", {1}}};

        map<string, double> wgtsJM{};

        for (auto& lep : leps) {
            ham.initEvent();

            Particle pB({5.280000000000000, 0., 0., 0.}, 511);
            Particle pDs({2.660051210321970, 1.389723696643482, 0.4037520573116902, 0.9702704595830642}, -413);
            Particle pEll({2.434889134265837, -1.420021405904506, -0.3418057134095529, -0.7985387637550296},
                          static_cast<PdgId>(-lep.second));
            Particle pNuEll({0.18506, 0.03029770926102430, -0.06194634390213728, -0.1717316958280348},
                            static_cast<PdgId>((lep.second) + 1));

            Process proc;

            auto idxB = proc.addParticle(pB);
            auto idxDs = proc.addParticle(pDs);
            auto idxEll = proc.addParticle(pEll);
            auto idxNuEll = proc.addParticle(pNuEll);

            proc.addVertex(idxB, {idxDs, idxEll, idxNuEll});

            ham.addProcess(proc);

            if (lep.first == "Tau") {
                ham.fillEventHistogram("pEllq2", {1.1, 1.1});
            } else {
                ham.fillEventHistogram("pEllq2", {0.6, 0.6});
            }

            ham.processEvent();

            if (lep.first == "Tau") {
                ham.setWilsonCoefficients("BtoCTauNu@JC", {{"R2", 1}});

                auto wgtRH = ham.getWeight("Scheme1", "RH");
                auto wgtJC = ham.getWeight("Scheme1", "JC");

                cout << wgtRH << " " << wgtJC << '\n';

                EXPECT_NEAR(wgtRH / wgtJC, 1., 1e-5);
            }

            ham.setWilsonCoefficients("BtoC" + lep.first + "Nu@JM", wcs[lep.first]);
            wgtsJM[lep.first] = ham.getWeight("Scheme1", "JM");
        }

        cout << wgtsJM["Tau"] << " " << wgtsJM["Mu"] << '\n';

        EXPECT_NEAR(wgtsJM["Tau"] / wgtsJM["Mu"], 1., 1e-5);

        ham.setWilsonCoefficients("BtoCMuNu", {{"S_qRlL", -1}, {"S_qLlL", 1}});
        ham.setWilsonCoefficients("BtoCTauNu@TS", {{"S1a", 1}});
        auto histoTS = ham.getHistogram("pEllq2", "Scheme1", "TS");
        auto histo = ham.getHistogram("pEllq2", "Scheme1");

        cout << histoTS[2 * 5 + 2].sumWi << " " << histo[1 * 5 + 1].sumWi << '\n';

        EXPECT_NEAR(histoTS[2 * 5 + 2].sumWi / wgtsJM["Tau"], 1., 1e-5);
        EXPECT_NEAR(histoTS[2 * 5 + 2].sumWi / histoTS[1 * 5 + 1].sumWi, 1., 1e-5);
    }

    TEST(HammerTest, Lambda) {

        Hammer ham;
        ham.setOptions("ProcessCalc: {Rates: false}");

        // Create included processes
        ham.includeDecay(vector<string>{"LbLcTauNu"});
        ham.addFFScheme("Scheme", {{"LbLc", "BLRS"}});
        ham.setFFInputScheme({{"LbLc", "PCR"}});
        ham.initRun();

        // Create process
        Process proc;

        // Particle decs
        Particle pLbmes(
            {6.0006081136124010e+00, -1.7205721657067963e+00, -1.0745365696783442e+00, 5.5992862596039261e-01}, 5122);
        Particle pLcmes(
            {2.6445298620131372e+00, -6.4721473009510044e-01, -6.5003404398021070e-01, -9.6135721471983748e-01}, 4122);
        Particle pTau(
            {2.3633597514416573e+00, -1.8573281291102936e-01, -8.1753561754710102e-01, 1.3135900572798234e+00}, 15);
        Particle kNuTauBar(
            {9.9271850015760621e-01, -8.8762462270066600e-01, 3.9303309184896784e-01, 2.0769578340040648e-01}, -16);
        Particle kNuTau(
            {1.0366235521488159e+00, -8.2115457296673355e-01, -3.2977241912807359e-01, 5.3994787502856967e-01}, 16);
        Particle kE({2.1228839569637614e-01, 7.8988747482349030e-03, -1.8848165295026625e-01, -9.7356952018238269e-02},
                    11);
        Particle kNuEBar(
            {1.1144478035964649e+00, 6.2752288530746925e-01, -2.9928154546876118e-01, 8.7099913426949183e-01}, -12);

        auto idxLbmes = proc.addParticle(pLbmes);
        auto idxLcmes = proc.addParticle(pLcmes);
        auto idxTau = proc.addParticle(pTau);
        auto idxNuTauBar = proc.addParticle(kNuTauBar);
        auto idxNuTau = proc.addParticle(kNuTau);
        auto idxE = proc.addParticle(kE);
        auto idxNuEBar = proc.addParticle(kNuEBar);

        proc.addVertex(idxLbmes, {idxLcmes, idxNuTauBar, idxTau});
        proc.addVertex(idxTau, {idxNuTau, idxE, idxNuEBar});

        ham.addProcess(proc);

        ham.processEvent();

        const std::vector<HashId>& processes{proc.getId()};
        auto weight = ham.getWeight("Scheme", processes);
        cout << "Weight: " << weight << '\n';

        ham.setWilsonCoefficients("BtoCTauNu", {{"SM", 2.0}});
        const std::vector<HashId>& processes2{proc.getId()};
        auto weight2 = ham.getWeight("Scheme", processes2);

        EXPECT_NEAR(weight2 / weight, 4., 1e-4);
    }

    TEST(HammerTest, D1) {

        Hammer ham;
        ham.setOptions("ProcessCalc: {Rates: false}");

        // Create included processes
        ham.includeDecay(vector<string>{"BD**1TauNu"});
        ham.addFFScheme("Scheme", {{"BD**1", "BLR_1"}});
        ham.setFFInputScheme({{"BD**1", "LLSW"}});
        ham.setUnits("GeV");
        ham.initRun();

        // Create process
        Process proc;

        // Particle decs
        Particle pBmes(
            {6.0006081136124010e+00, -1.7205721657067963e+00, -1.0745365696783442e+00, 5.5992862596039261e-01}, -511);
        Particle pD1mes(
            {2.6445298620131372e+00, -6.4721473009510044e-01, -6.5003404398021070e-01, -9.6135721471983748e-01}, 10413);
        Particle pTau(
            {2.3633597514416573e+00, -1.8573281291102936e-01, -8.1753561754710102e-01, 1.3135900572798234e+00}, 15);
        Particle kNuTauBar(
            {9.9271850015760621e-01, -8.8762462270066600e-01, 3.9303309184896784e-01, 2.0769578340040648e-01}, -16);
        Particle kNuTau(
            {1.0366235521488159e+00, -8.2115457296673355e-01, -3.2977241912807359e-01, 5.3994787502856967e-01}, 16);
        Particle kE({2.1228839569637614e-01, 7.8988747482349030e-03, -1.8848165295026625e-01, -9.7356952018238269e-02},
                    11);
        Particle kNuEBar(
            {1.1144478035964649e+00, 6.2752288530746925e-01, -2.9928154546876118e-01, 8.7099913426949183e-01}, -12);

        auto idxBmes = proc.addParticle(pBmes);
        auto idxD1mes = proc.addParticle(pD1mes);
        auto idxTau = proc.addParticle(pTau);
        auto idxNuTauBar = proc.addParticle(kNuTauBar);
        auto idxNuTau = proc.addParticle(kNuTau);
        auto idxE = proc.addParticle(kE);
        auto idxNuEBar = proc.addParticle(kNuEBar);

        proc.addVertex(idxBmes, {idxD1mes, idxNuTauBar, idxTau});
        proc.addVertex(idxTau, {idxNuTau, idxE, idxNuEBar});

        ham.addProcess(proc);

        ham.processEvent();

        const std::vector<HashId>& processes{proc.getId()};
        auto weight = ham.getWeight("Scheme", processes);
        cout << "Weight: " << weight << '\n';

        ham.setWilsonCoefficients("BtoCTauNu", {{"SM", 2.0}});
        const std::vector<HashId>& processes2{proc.getId()};
        auto weight2 = ham.getWeight("Scheme", processes2);

        EXPECT_NEAR(weight2 / weight, 4., 1e-4);
    }

    TEST(HammerTest, PHOTOSandOrdering) {

        Hammer ham;
        ham.setOptions("ProcessCalc: {Rates: false}");

        // Create included processes
        ham.includeDecay(vector<string>{"BD*TauNu"});
        ham.addFFScheme("Scheme", {{"BD*", "BLPR"}});
        ham.setFFInputScheme({{"BD*", "ISGW2"}});
        ham.setUnits("GeV");
        ham.initRun();

        // Create process
        Process proc1;
        Process proc2;

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

        auto idx1B = proc1.addParticle(pBmes);
        auto idx1Dstar = proc1.addParticle(pDstar);
        auto idx1Tau = proc1.addParticle(pTau);
        auto idx1NuTauBar = proc1.addParticle(pNuTauBar);
        auto idx1Gamma = proc1.addParticle(pGamma);
        auto idx1D = proc1.addParticle(pD);
        auto idx1Pi = proc1.addParticle(pPi);
        auto idx1Mu = proc1.addParticle(pMu);
        auto idx1NuTau = proc1.addParticle(pNuTau);
        auto idx1NuMuBar = proc1.addParticle(pNuMuBar);

        proc1.addVertex(idx1B, {idx1Dstar, idx1NuTauBar, idx1Tau, idx1Gamma});
        proc1.addVertex(idx1Dstar, {idx1D, idx1Pi});
        proc1.addVertex(idx1Tau, {idx1NuTau, idx1NuMuBar, idx1Mu});

        // Now add in jumbled order
        int flip = -1;
        Particle pBmes2{
            {5.2969140979437572e+00, 3.2824493742660787e-01, -2.3429865618090287e-01, 1.4385195640293111e-01},
            -511 * flip};
        Particle pDstar2{
            {2.3470774909812628e+00, 8.0970717778699397e-01, 7.6332563305857859e-01, 4.7945687459445407e-01},
            413 * flip};
        Particle pTau2{
            {2.0795408878883048e+00, -8.0737090039644110e-01, -6.1305090764429815e-01, 3.7379380121775035e-01},
            15 * flip};
        Particle pNuTauBar2{
            {8.7026750024654076e-01, 3.2592376941100554e-01, -3.8455135916505445e-01, -7.0940783099892080e-01},
            -16 * flip};
        Particle pGamma2{
            {2.8218818201336615e-05, -1.5109375535821102e-05, -2.2022429710929912e-05, 9.1115893907028467e-06}, 22};
        Particle pD2{{2.1652224797353812e+00, 7.5505747015863178e-01, 6.8723913013095317e-01, 4.1006189279903876e-01},
                     421 * flip};
        Particle pPi2{{1.8185501124588163e-01, 5.4649707628362069e-02, 7.6086502927625407e-02, 6.9394981795415422e-02},
                      211 * flip};
        Particle pMu2{
            {1.1411538114089392e+00, -5.8718874409586341e-01, -4.5464791308157315e-01, 8.5998429292620537e-01},
            13 * flip};
        Particle pNuTau2{
            {3.7918042426691972e-01, 2.7540552999632262e-01, -1.2352239877392272e-01, -2.2950338820118057e-01},
            16 * flip};
        Particle pNuMuBar2{
            {5.5920665221244559e-01, -4.9558768629690025e-01, -3.4880595788802433e-02, -2.5668710350727436e-01},
            -14 * flip};

        auto idx2Gamma = proc2.addParticle(pGamma2);
        auto idx2Dstar = proc2.addParticle(pDstar2);
        auto idx2Tau = proc2.addParticle(pTau2);
        auto idx2NuTauBar = proc2.addParticle(pNuTauBar2);
        auto idx2B = proc2.addParticle(pBmes2);
        auto idx2D = proc2.addParticle(pD2);
        auto idx2Pi = proc2.addParticle(pPi2);
        auto idx2Mu = proc2.addParticle(pMu2);
        auto idx2NuTau = proc2.addParticle(pNuTau2);
        auto idx2NuMuBar = proc2.addParticle(pNuMuBar2);

        proc2.addVertex(idx2Dstar, {idx2Pi, idx2D});
        proc2.addVertex(idx2B, {idx2Dstar, idx2Tau, idx2NuTauBar, idx2Gamma});
        proc2.addVertex(idx2Tau, {idx2NuTau, idx2NuMuBar, idx2Mu});

        auto procId1 = ham.addProcess(proc1);
        auto procId2 = ham.addProcess(proc2);

        ham.processEvent();

        ham.setWilsonCoefficients("BtoCTauNu", {{"S_qLlL", 0.5}});
        auto weights = ham.getWeights("Scheme");
        for (auto& elem : weights) {
            cout << "Id: " << elem.first << "\t Weight: " << elem.second << '\n';
        }

        EXPECT_EQ(weights[procId1] / weights[procId2], 1.);
    }

} // namespace Hammer
