///
/// @file  TestSpecializationDefinitions.cc
/// @brief Tests for SpecializationDefinitions
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include <fstream>
#include <sstream>
#include <vector>

#include "gtest/gtest.h"

#include "Hammer/ProvidersRepo.hh"
#include "Hammer/SpecializationDefinitions.hh"
#include "Hammer/Math/MultiDim/VectorContainer.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Amplitudes/AmplBToQLepNuBase.hh"
#include "Hammer/Math/MultiDim/Operations.hh"
#include "Hammer/Tools/HammerSerial.hh"

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

    class DummyProvidersRepo : public IWCFFErrProviders {

    public:

        DummyProvidersRepo(AmplitudeBase* amp1, AmplitudeBase* amp2) : _amp1{amp1}, _amp2{amp2} {
        }
        DummyProvidersRepo(const DummyProvidersRepo&) = delete;
        DummyProvidersRepo& operator=(const DummyProvidersRepo&) = delete;
        DummyProvidersRepo(DummyProvidersRepo&&) = default;
        DummyProvidersRepo& operator=(DummyProvidersRepo&&) = default;
        ~DummyProvidersRepo() override {
            _amp1 = nullptr;
            _amp2 = nullptr;
        }

        [[nodiscard]] AmplitudeBase* getWCProvider(const std::string& wcPrefix) const override {
            if (wcPrefix == "BtoCMuNu") {
                return _amp1;
            }
            if (wcPrefix == "BtoCTauNu") {
                return _amp2;
            }
            return nullptr;
        }
        [[nodiscard]] IndexLabel getWCLabel(const std::string& wcPrefix) const override {
            if (wcPrefix == "BtoCMuNu") {
                return _amp1->getWCInfo().second;
            }
            if (wcPrefix == "BtoCTauNu") {
                return _amp2->getWCInfo().second;
            }
            return NONE;
        }

        [[nodiscard]] std::map<IndexLabel, AmplitudeBase*> getAllWCProviders() const override {
            return {};
        }

        [[nodiscard]] FormFactorBase* getFFErrProvider(const FFPrefixGroup& /*process*/) const override {
            return nullptr;
        }
        [[nodiscard]] IndexLabel getFFErrLabel(const FFPrefixGroup& /*process*/) const override {
            return NONE;
        }
        [[nodiscard]] std::map<IndexLabel, SchemeDict<FormFactorBase*>> getAllFFErrProviders() const override {
            return {};
        }

        [[nodiscard]] std::set<std::string>
        schemeNamesFromPrefixAndGroup(const FFPrefixGroup& /*value*/) const override {
            return {};
        }

        [[nodiscard]] AmplitudeBase* getAmplitude(PdgId /*parent*/, const std::vector<PdgId>& /*daughters*/,
                                                  const std::vector<PdgId>& /*granddaughters*/) const override {
            return nullptr;
        }

    private:

        AmplitudeBase* _amp1;
        AmplitudeBase* _amp2;
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

        WCSpecialization spec;
        WCSpecialization spec2{"BtoCMuNu", WILSON_BCMUNU, "RH", {"DIM1", "DIM2"}};
        WCSpecialization spec3{"BtoCTauNu", WILSON_BCTAUNU, "JM", {}};

        spec2.setAmplitude(&ampl2);
        spec3.setAmplitude(&ampl);

        spec2.initialize();
        spec3.initialize();

        spec.setSubspaceBasis({{{"T_qLlL", 4}, {"S_qRlL", 1}}, {{"T_qLlL", 3}, {"S_qRlL", -7}}});
        spec3.setSubspaceBasis({{{"T_qLlL", 4}, {"S_qRlL", 1}}, {{"T_qLlL", 3}, {"S_qRlL", -7}}});
        spec2.setSubspaceBasis({{{"T_qLlL", 4}, {"S_qRlL", 1}}, {{"T_qLlL", 3}, {"S_qRlL", -7}}});

        spec.setSubspaceOrigin({{"SM", 1}, {"S_qRlL", 1}});
        spec2.setSubspaceOrigin({{"SM", 1}, {"S_qRlL", 2}});
        spec3.setSubspaceOrigin({{"SM", 1}, {"S_qRlL", 3}});

        EXPECT_FALSE(spec.getProjectionTensor());
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


    TEST(SpecializationDefinitionsTest, addgetSetting) {

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

        DummyProvidersRepo provs{&ampl2, &ampl};
        SpecializationDefinitions specs{&provs};

        specs.addSpecialization("BtoCMuNu", "RH", {"DIM1", "DIM2"});
        specs.addSpecialization("BtoCTauNu", "JM", {});
        specs.addSpecialization("BtoCTauNu", "RH", {});

        const auto& allRH = specs.getSpecializations("RH");
        const auto& allJM = specs.getSpecializations("JM");

        EXPECT_EQ(allRH.size(), 2);
        EXPECT_EQ(allJM.size(), 1);

        EXPECT_EQ(allRH.find(WILSON_BCMUNU)->second.getCoordinates().size(), 2);
        EXPECT_EQ(allRH.find(WILSON_BCTAUNU)->second.getCoordinates().size(), 0);
        EXPECT_TRUE(allJM.find(WILSON_BCMUNU) == allJM.end());
        EXPECT_EQ(allJM.find(WILSON_BCTAUNU)->second.getCoordinates().size(), 0);

        auto specsIds = specs.specializationIds();
        EXPECT_EQ(specsIds.size(), 2);
        EXPECT_TRUE(specsIds.find("JM") != specsIds.end());
        EXPECT_TRUE(specsIds.find("RH") != specsIds.end());
        EXPECT_EQ(specs.specializationIds(WILSON_BCENU).size(), 0);
        EXPECT_EQ(specs.specializationIds(WILSON_BCMUNU).size(), 1);
        EXPECT_EQ(specs.specializationIds(WILSON_BCTAUNU).size(), 2);

        specs.removeSpecialization("RH", "BtoCENU");
        EXPECT_EQ(specs.getSpecializations("RH").size(), 2);
        specs.removeSpecialization("RH", WILSON_BCMUNU);
        EXPECT_EQ(specs.getSpecializations("RH").size(), 1);
        specs.removeSpecializations("JM");
        EXPECT_THROW(specs.getSpecializations("JM"), Error);

        specs.clear();
        EXPECT_EQ(specs.specializationIds().size(), 0);
    }

    TEST(SpecializationDefinitionsTest, Yaml) {
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

        DummyProvidersRepo provs{&ampl2, &ampl};
        SpecializationDefinitions specs{&provs};

        specs.addSpecialization("BtoCMuNu", "RH", {"DIM1", "DIM2"});
        specs.addSpecialization("BtoCTauNu", "JM", {});
        specs.addSpecialization("BtoCTauNu", "RH", {});

        specs.getSpecialization("JM", "BtoCTauNu").setAmplitude(&ampl);
        specs.getSpecialization("RH", "BtoCTauNu").setAmplitude(&ampl);
        specs.getSpecialization("RH", "BtoCMuNu").setAmplitude(&ampl2);

        specs.getSpecialization("JM", "BtoCTauNu").initialize();
        specs.getSpecialization("RH", "BtoCTauNu").initialize();
        specs.getSpecialization("RH", "BtoCMuNu").initialize();

        specs.getSpecialization("JM", "BtoCTauNu").setSubspaceOrigin({{"SM", 1}, {"S_qRlL", 1}});
        specs.getSpecialization("RH", "BtoCTauNu").setSubspaceOrigin({{"SM", 1}, {"S_qRlL", 2}});
        specs.getSpecialization("RH", "BtoCMuNu").setSubspaceOrigin({{"SM", 1}, {"S_qRlL", 3}});

        specs.getSpecialization("RH", "BtoCMuNu")
            .setSubspaceBasis({{{"T_qLlL", 4}, {"S_qRlL", 1}}, {{"T_qLlL", 3}, {"S_qRlL", -7}}});

        sh.saveSettings("test.yaml", false);
        std::ifstream ifs("test.yaml");
        std::ifstream ifs2("test_spec_settings.yaml");
        std::string content((std::istreambuf_iterator<char>(ifs)), (std::istreambuf_iterator<char>()));
        std::string content2((std::istreambuf_iterator<char>(ifs2)), (std::istreambuf_iterator<char>()));
        EXPECT_EQ(content, content2);

        ofstream file;
        YAML::Emitter emitter(file);
        file.open("test2.yaml");
        if (file.is_open()) {
            emitter << YAML::BeginMap;
            emitter << YAML::Key << "WCSpecializations";
            emitter << YAML::Value << specs;
            emitter << YAML::EndMap;
            emitter << YAML::Newline;
            file.close();
        }

        std::ifstream ifs3("test2.yaml");
        std::ifstream ifs4("test_specializations.yaml");
        std::string content3((std::istreambuf_iterator<char>(ifs3)), (std::istreambuf_iterator<char>()));
        std::string content4((std::istreambuf_iterator<char>(ifs4)), (std::istreambuf_iterator<char>()));
        EXPECT_EQ(content3, content4);
    }

    TEST(SpecializationDefinitionsTest, ReadWrite) {
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

        DummyProvidersRepo provs{&ampl2, &ampl};
        SpecializationDefinitions specs{&provs};

        specs.addSpecialization("BtoCMuNu", "RH", {"DIM1", "DIM2"});
        specs.addSpecialization("BtoCTauNu", "JM", {});
        specs.addSpecialization("BtoCTauNu", "RH", {});

        specs.getSpecialization("JM", "BtoCTauNu").setAmplitude(&ampl);
        specs.getSpecialization("RH", "BtoCTauNu").setAmplitude(&ampl);
        specs.getSpecialization("RH", "BtoCMuNu").setAmplitude(&ampl2);

        specs.getSpecialization("JM", "BtoCTauNu").initialize();
        specs.getSpecialization("RH", "BtoCTauNu").initialize();
        specs.getSpecialization("RH", "BtoCMuNu").initialize();

        specs.getSpecialization("JM", "BtoCTauNu").setSubspaceOrigin({{"SM", 1}, {"S_qRlL", 1}});
        specs.getSpecialization("RH", "BtoCTauNu").setSubspaceOrigin({{"SM", 1}, {"S_qRlL", 3}});
        specs.getSpecialization("RH", "BtoCMuNu").setSubspaceOrigin({{"SM", 1}, {"S_qRlL", 2}});

        specs.getSpecialization("RH", "BtoCMuNu")
            .setSubspaceBasis({{{"T_qLlL", 4}, {"S_qRlL", 1}}, {{"T_qLlL", 3}, {"S_qRlL", -7}}});


        auto builder = make_unique<flatbuffers::FlatBufferBuilder>();
        flatbuffers::Offset<Serial::FBHeader> msg;
        builder->Clear();
        vector<flatbuffers::Offset<Serial::FBSpecialization>> spcs;
        specs.write(builder.get(), &spcs);
        auto serialspecs = builder->CreateVector(spcs);
        Serial::FBHeaderBuilder serialheader{*builder};
        serialheader.add_specs(serialspecs);
        auto headoffset = serialheader.Finish();
        builder->Finish(headoffset);
        uint8_t* buf = builder->GetBufferPointer();
        // unsigned int size = builder->GetSize();
        const auto* ser1 = flatbuffers::GetRoot<Serial::FBHeader>(buf);
        EXPECT_NE(ser1, nullptr);
        if (ser1 != nullptr) {
            SpecializationDefinitions specs2{&provs};
            specs2.read(ser1, false);

            const auto& allRH = specs2.getSpecializations("RH");
            const auto& allJM = specs2.getSpecializations("JM");

            EXPECT_EQ(allRH.size(), 2);
            EXPECT_EQ(allJM.size(), 1);

            EXPECT_EQ(allRH.find(WILSON_BCMUNU)->second.getCoordinates().size(), 2);
            EXPECT_EQ(allRH.find(WILSON_BCTAUNU)->second.getCoordinates().size(), 0);
            EXPECT_TRUE(allJM.find(WILSON_BCMUNU) == allJM.end());
            EXPECT_EQ(allJM.find(WILSON_BCTAUNU)->second.getCoordinates().size(), 0);

            Tensor compt2{"proj",
                          MD::makeVector({11, 3},
                                         {WILSON_BCMUNU, specializeLabel(WILSON_BCMUNU, generateSpecializedPad("RH"))},
                                         {1, 0, 0, 0, 0, 0, 0, 0, 0, -2, -1, 7,  0,  0, 0, 0, 0,
                                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0,  -4, -3, 0, 0, 0})};

            Tensor compt3{"proj", MD::makeVector({11}, {WILSON_BCTAUNU}, {1, 0, 0, -3, 0, 0, 0, 0, 0, 0, 0})};
            auto t2 = specs2.getSpecialization("RH", "BtoCMuNu").getProjectionTensor();
            auto t3 = specs2.getSpecialization("RH", "BtoCTauNu").getProjectionTensor();
            for (IndexType idx1 = 0; idx1 < 11; ++idx1) {
                for (IndexType idx2 = 0; idx2 < 3; ++idx2) {
                    double com2Re = compareVals(t2->element({idx1, idx2}).real(), compt2.element({idx1, idx2}).real());
                    EXPECT_DOUBLE_EQ(com2Re, 1.);
                }
                double com3Re = compareVals(t3->element({idx1}).real(), compt3.element({idx1}).real());
                EXPECT_DOUBLE_EQ(com3Re, 1.);
            }
        }
    }


} // namespace Hammer
