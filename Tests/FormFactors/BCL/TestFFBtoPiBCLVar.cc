///
/// @file  TestFFBtoPiBCLVar.cc
/// @brief Tests for FFBtoPiBCLVar
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#include <filesystem>
#include <fstream>

#define FORMFACTORBASE_FRIENDS                  \
    FRIEND_TEST(FFBtoPiBCLVarTest, eval);       \
    FRIEND_TEST(FFBtoPiBCLVarTest, evalq2cons); \
    FRIEND_TEST(FFBtoPiBCLVarTest, clone);      \
    FRIEND_TEST(FFBtoPiBCLVarTest, addRefs)

#include "Hammer/FormFactors/BCL/FFBtoPiBCLVar.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/FourMomentum.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Math/MultiDim/VectorContainer.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Tools/SettingsConsumer.hh"
// May need other stuff here

using namespace std;
using namespace complex_literals;

namespace Hammer {

    namespace MD = MultiDimensional;

    TEST(FFBtoPiBCLVarTest, eval) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs (q^2 = 19.8958 GeV^2)
        FourMomentum pBmes{28783.8884892086, 0, 0, 28295.4737822005};
        FourMomentum pPimes{139., 0, 0, 0};
        FourMomentum kNuMu{483.601218826, 143.20931144041, 440.75294019035, 138.1957947273};
        FourMomentum pMu{28161.2872703827, -143.20931144041, -440.75294019035, 28157.2779874733};

        // Evaluate the FF class
        FFBtoPiBCLVar ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();
        ff.eval(Particle(pBmes, 511), {Particle(pPimes, -211), Particle(kNuMu, 14), Particle(pMu, -13)}, {});


        // Compare to direct evaluation
        Tensor ffEval{"ffEval", MD::makeVector({4, 9}, {FF_BPI, FF_BPI_VAR},
                                               {0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.50165560937469,
                                                -0.007821777518211096,
                                                -0.0004051044675213387,
                                                0.008995864646937461,
                                                0.010860431923679447,
                                                -0.0013621374270545716,
                                                0.009913274320449764,
                                                -0.00006568438560200482,
                                                -0.0002053165614030617,
                                                1.3961877732509,
                                                -0.0010469251541312754,
                                                -0.016688970856905454,
                                                -0.005614453604649242,
                                                -0.004432619020235696,
                                                0.024959265673772685,
                                                0.030772704276851758,
                                                -0.0029438756729662126,
                                                0.0007114724440708648,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.})};

        auto& t = ff.getTensor();
        for (IndexType idx1 = 0; idx1 < 4; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 9; ++idx2) {
                EXPECT_NEAR(compareVals(t.element({idx1, idx2}).real(), ffEval.element({idx1, idx2}).real()), 1.,
                            1.e-4);
                EXPECT_DOUBLE_EQ(t.element({idx1, idx2}).imag(), ffEval.element({idx1, idx2}).imag());
            }
        }
    }

    TEST(FFBtoPiBCLVarTest, evalq2cons) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs (q^2 = 19.8958 GeV^2)
        FourMomentum pBmes{28783.8884892086, 0, 0, 28295.4737822005};
        FourMomentum pPimes{139., 0, 0, 0};
        FourMomentum kNuMu{483.601218826, 143.20931144041, 440.75294019035, 138.1957947273};
        FourMomentum pMu{28161.2872703827, -143.20931144041, -440.75294019035, 28157.2779874733};

        // Evaluate the FF class
        FFBtoPiBCLVar ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        vector<double> apvec = {0.415, -0.488, -0.31};
        vector<double> a0vec = {0.500, -1.424};

        // Row basis is FP a012, FZ a01 2206.07501 Table 57
        vector<vector<double>> eigmat{{0.00619374, 0.000447436, 0.00473847, 0.00951309, 0.00666974},
                                      {0.0437911, 0.00822887, -0.0283399, -0.00407376, 0.00197006},
                                      {-0.179832, 0.00395287, -0.00663112, -0.000535919, 0.00067543},
                                      {-0.000228852, -0.0113757, -0.0104179, 0.0167898, -0.00301823},
                                      {0.00646603, 0.0533757, 0.00260018, 0.00416633, -0.00105292}};
        set.changeSetting<vector<double>>("BtoPiBCLVar", "ap", apvec);
        set.changeSetting<vector<double>>("BtoPiBCLVar", "a0", a0vec);
        set.changeSetting<vector<vector<double>>>("BtoPiBCLVar", "eigmatrix", eigmat);
        set.changeSetting<bool>("BtoPiBCLVar", "q2cons", true);
        ff.calcUnits();
        ff.eval(Particle(pBmes, 511), {Particle(pPimes, -211), Particle(kNuMu, 14), Particle(pMu, -13)}, {});


        // Compare to direct evaluation
        Tensor ffEval{"ffEval", MD::makeVector({4, 6}, {FF_BPI, FF_BPI_VAR},
                                               {0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.49302579294957,
                                                -0.00019703026564474668,
                                                -0.011112722166471082,
                                                -0.010403200355241908,
                                                0.016807378693426692,
                                                -0.003020143665495722,
                                                1.382905939023,
                                                0.02146916381031594,
                                                0.0016360327658263404,
                                                0.015413482126354398,
                                                0.031818602644046665,
                                                0.02238822563411089,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.})};

        auto& t = ff.getTensor();
        for (IndexType idx1 = 0; idx1 < 4; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 6; ++idx2) {
                EXPECT_NEAR(compareVals(t.element({idx1, idx2}).real(), ffEval.element({idx1, idx2}).real()), 1.,
                            1.e-4);
                EXPECT_DOUBLE_EQ(t.element({idx1, idx2}).imag(), ffEval.element({idx1, idx2}).imag());
            }
        }
    }

    TEST(FFBtoPiBCLVarTest, clone) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoPiBCLVar ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto cloned = ff.clone("testLabel");
        ASSERT_NE(cloned, nullptr);
        EXPECT_NE(cloned->group().find("testLabel"), string::npos);

        FourMomentum pBmes{28783.8884892086, 0, 0, 28295.4737822005};
        FourMomentum pPimes{139., 0, 0, 0};
        FourMomentum kNuMu{483.601218826, 143.20931144041, 440.75294019035, 138.1957947273};
        FourMomentum pMu{28161.2872703827, -143.20931144041, -440.75294019035, 28157.2779874733};

        cloned->setSettingsHandler(set);
        cloned->initSettings();
        cloned->calcUnits();
        cloned->eval(Particle(pBmes, 511), {Particle(pPimes, -211), Particle(kNuMu, 14), Particle(pMu, -13)}, {});

        ff.eval(Particle(pBmes, 511), {Particle(pPimes, -211), Particle(kNuMu, 14), Particle(pMu, -13)}, {});

        auto& tOrig = ff.getTensor();
        auto& tCloned = cloned->getTensor();
        for (IndexType idx1 = 0; idx1 < 4; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 9; ++idx2) {
                EXPECT_NEAR(compareVals(tOrig.element({idx1, idx2}).real(), tCloned.element({idx1, idx2}).real()), 1.,
                            1.e-4);
                EXPECT_DOUBLE_EQ(tOrig.element({idx1, idx2}).imag(), tCloned.element({idx1, idx2}).imag());
            }
        }
    }

    TEST(FFBtoPiBCLVarTest, addRefs) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoPiBCLVar ff;
        ff.setSettingsHandler(set);
        ff.initSettings();

        EXPECT_FALSE(set.checkReference("Aoki:2019cca"));
        EXPECT_FALSE(set.checkReference("FermilabLattice:2015mwy"));

        static_cast<SettingsConsumer*>(&ff)->addRefs();
        EXPECT_TRUE(set.checkReference("Aoki:2019cca"));
        EXPECT_TRUE(set.checkReference("FermilabLattice:2015mwy"));

        // Second call must be idempotent.
        static_cast<SettingsConsumer*>(&ff)->addRefs();
        EXPECT_TRUE(set.checkReference("Aoki:2019cca"));
        EXPECT_TRUE(set.checkReference("FermilabLattice:2015mwy"));

        const string tmpfile = testing::TempDir() + "hammer_bcl_btopivar_refs_test.bib";
        set.saveReferences(tmpfile);

        ifstream reffile(tmpfile);
        ASSERT_TRUE(reffile.is_open());
        string content((istreambuf_iterator<char>(reffile)), istreambuf_iterator<char>());
        reffile.close();
        std::filesystem::remove(tmpfile);

        EXPECT_NE(content.find("Aoki:2019cca"), string::npos);
        EXPECT_NE(content.find("FermilabLattice:2015mwy"), string::npos);
    }

} // namespace Hammer
