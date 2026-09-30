///
/// @file  TestWCSpecialization.cc
/// @brief Tests for WCSpecialization
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include <fstream>
#include <sstream>
#include <vector>

#include "gtest/gtest.h"

#define WCSPECIALIZATION_FRIENDS friend class PartialSpecFixture


#include "Hammer/WCSpecialization.hh"
#include "Hammer/Math/MultiDim/VectorContainer.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Amplitudes/AmplBToQLepNuBase.hh"
#include "Hammer/Math/MultiDim/Operations.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;


    class TestableAmplBToQLepNuBase : public AmplBToQLepNuBase {

    public:

        TestableAmplBToQLepNuBase() = default;

        using AmplBToQLepNuBase::addProcessSignature;
        using AmplBToQLepNuBase::addTensor;
        using AmplBToQLepNuBase::defineSettings;
        using AmplBToQLepNuBase::preProcessWCValues;
        using AmplBToQLepNuBase::updateWilsonCoeffLabelPrefix;

        void eval(const Particle& /*parent*/, const ParticleList& /*daughters*/,
                  const ParticleList& /*references*/) override {
            Tensor& t = getTensor();
            t.clearData();
        }
    };

    TEST(WCSpecializationTest, Creation) {
        WCSpecialization spec;
        WCSpecialization spec2{"BtoCMuNu", WILSON_BCMUNU, "RH", {"DIM1", "DIM2"}};
        WCSpecialization spec3{"BtoCTauNu", WILSON_BCTAUNU, "JM", {}};
        EXPECT_EQ(spec.getBaseLabel(), NONE);
        EXPECT_EQ(spec2.getBaseLabel(), WILSON_BCMUNU);
        EXPECT_EQ(spec3.getBaseLabel(), WILSON_BCTAUNU);
        EXPECT_EQ(spec.getPrefixId().id, "");
        EXPECT_EQ(spec2.getPrefixId().id, "RH");
        EXPECT_EQ(spec3.getPrefixId().id, "JM");
        EXPECT_EQ(spec.getPrefixId().prefix, "");
        EXPECT_EQ(spec2.getPrefixId().prefix, "BtoCMuNu");
        EXPECT_EQ(spec3.getPrefixId().prefix, "BtoCTauNu");
        auto pad2 = generateSpecializedPad("RH");
        EXPECT_EQ(spec.getPad(), 0ul);
        EXPECT_EQ(spec2.getPad(), pad2);
        EXPECT_EQ(spec3.getPad(), 0ul);
        EXPECT_EQ(spec.getFullLabel(), 0);
        EXPECT_EQ(spec2.getFullLabel(), specializeLabel(WILSON_BCMUNU, pad2));
        EXPECT_EQ(spec3.getFullLabel(), NONE);
        EXPECT_FALSE(spec.isPartialSpecialization());
        EXPECT_TRUE(spec2.isPartialSpecialization());
        EXPECT_FALSE(spec3.isPartialSpecialization());
        EXPECT_EQ(spec.getCoordinates().size(), 0ul);
        EXPECT_EQ(spec2.getCoordinates().size(), 2ul);
        EXPECT_EQ(spec3.getCoordinates().size(), 0ul);
        EXPECT_EQ(spec2.getCoordinates()[0], "DIM1");
        EXPECT_EQ(spec2.getCoordinates()[1], "DIM2");
    }

    TEST(WCSpecializationTest, GetSetProjectionTensor) {
        SettingsHandler sh;
        TestableAmplBToQLepNuBase ampl;
        TestableAmplBToQLepNuBase ampl2;
        ampl.setSettingsHandler(sh);
        ampl2.setSettingsHandler(sh);

        vector<IndexType> dims{{11, 4, 2, 2, 2}};
        string name{"AmplBDLepNu"};
        ampl.addProcessSignature(PID::BPLUS, {-PID::D0, PID::NU_TAU, PID::ANTITAU});
        ampl.addTensor(
            Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCTAUNU, FF_BD, SPIN_NUTAU, SPIN_NUTAU_REF, SPIN_TAUP})});
        ampl.updateWilsonCoeffLabelPrefix();
        ampl.defineSettings();

        ampl2.addProcessSignature(PID::BPLUS, {-PID::D0, PID::NU_MU, PID::ANTIMUON});
        ampl2.addTensor(
            Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCMUNU, FF_BD, SPIN_NUMU, SPIN_NUMU_REF, SPIN_MUP})});
        ampl2.updateWilsonCoeffLabelPrefix();
        ampl2.defineSettings();

        // WCSpecialization spec;
        WCSpecialization spec2{"BtoCMuNu", WILSON_BCMUNU, "RH", {"DIM1", "DIM2"}};
        WCSpecialization spec3{"BtoCTauNu", WILSON_BCTAUNU, "JM", {}};

        spec2.setAmplitude(&ampl2);
        spec3.setAmplitude(&ampl);

        spec2.initialize();
        spec3.initialize();

        // spec.setSubspaceOrigin({{"SM", 1}, {"S_qRlL", 1}});
        spec2.setSubspaceOrigin({{"SM", 1}, {"S_qRlL", 2}});
        spec3.setSubspaceOrigin({{"SM", 1}, {"S_qRlL", 3}});

        // spec.setSubspaceBasis({{{"T_qLlL", 4}, {"S_qRlL", 1}}, {{"T_qLlL", 3}, {"S_qRlL", -7}}});
        spec3.setSubspaceBasis({{{"T_qLlL", 4}, {"S_qRlL", 1}}, {{"T_qLlL", 3}, {"S_qRlL", -7}}});
        spec2.setSubspaceBasis({{{"T_qLlL", 4}, {"S_qRlL", 1}}, {{"T_qLlL", 3}, {"S_qRlL", -7}}});

        // EXPECT_FALSE(spec.getProjectionTensor());
        EXPECT_EQ(spec2.getProjectionTensor()->rank(), 2);
        EXPECT_EQ(spec3.getProjectionTensor()->rank(), 1);

        Tensor compt2{"proj",
                      MD::makeVector({11, 3},
                                     {WILSON_BCMUNU, specializeLabel(WILSON_BCMUNU, generateSpecializedPad("RH"))},
                                     {1, 0, 0, 0, 0, 0, 0, 0, 0, -2, -1, 7,  0,  0, 0, 0, 0,
                                      0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0,  -4, -3, 0, 0, 0})};

        Tensor compt3{"proj", MD::makeVector({11}, {WILSON_BCTAUNU}, {1, 0, 0, -3, 0, 0, 0, 0, 0, 0, 0})};
        auto t2 = spec2.getProjectionTensor();
        auto t3 = spec3.getProjectionTensor();
        for (IndexType idx1 = 0; idx1 < 11; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 3; ++idx2) {
                double com2Re = compareVals(t2->element({idx1, idx2}).real(), compt2.element({idx1, idx2}).real());
                EXPECT_DOUBLE_EQ(com2Re, 1.);
            }
            double com3Re = compareVals(t3->element({idx1}).real(), compt3.element({idx1}).real());
            EXPECT_DOUBLE_EQ(com3Re, 1.);
        }
        spec2.setSubspaceVector("DIM2", {{"T_qLlL", 3}, {"S_qRlL", 7}});
        EXPECT_DOUBLE_EQ(spec2.getProjectionTensor()->element({9, 2}).real(), -3.);
        EXPECT_DOUBLE_EQ(spec2.getProjectionTensor()->element({3, 2}).real(), -7.);
    }


} // namespace Hammer
