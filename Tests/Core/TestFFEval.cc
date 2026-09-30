///
/// @file  TestFFEval.cc
/// @brief Tests for FFEval
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "Hammer/FFEval.hh"
#include "Hammer/Math/Tensor.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Tools/Utils.hh"
#include "Hammer/DictionaryManager.hh"
#include "Hammer/PurePhaseSpaceDefs.hh"
#include "Hammer/ProcessDefinitions.hh"
#include "Hammer/SchemeDefinitions.hh"
#include "Hammer/ProvidersRepo.hh"
#include "Hammer/FormFactorBase.hh"
#include "Hammer/ExternalData.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/IContainer.hh"
#include "Hammer/Math/MultiDim/VectorContainer.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"
#include "gtest/gtest.h"

#include <boost/functional/hash.hpp>

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    class DummyPurePhaseSpaceDefs : public PurePhaseSpaceDefs {
    public:

        DummyPurePhaseSpaceDefs() = default;

        ~DummyPurePhaseSpaceDefs() override = default;
    };
    
    class DummySchemeDefinitions : public SchemeDefinitions {
    public:

        DummySchemeDefinitions() {
           // _schemes.insert({"Scheme1", {{10000,0}}});
        }

        ~DummySchemeDefinitions() override = default;

        [[nodiscard]] SchemeDict<FFIndex> getFFSchemesForProcess(HashId /*id*/) const override {
            return {{"Scheme1", 0}};
        }

//    private:
//
//        SchemeDict<std::map<HashId, FFIndex>> _schemes;
    };
    
    class DummyProcessDefinitions : public ProcessDefinitions {
    public:

        DummyProcessDefinitions() = default;

        ~DummyProcessDefinitions() override = default;
        
    };

    class DummyFFBtoDX : public FormFactorBase { 
    public:   
        
        DummyFFBtoDX() : FormFactorBase{} {  
            addTensor(Tensor{"eval", MD::makeEmptySparse({2,2}, {FF_BD, FF_BD_VAR})});
            setSignatureIndex(0);
        };
        
        void evalAtPSPoint(const vector<double>& point, const vector<double>& masses) override {
            Tensor& result = getTensor();
            result.element({0,0}) = 2.;
            result.element({1,0}) = 1.;
            
            result.element({0,1}) = 0.5;
            result.element({1,1}) = -0.5;
        }
        
        void defineSettings() override { 
            
        }
        
        std::unique_ptr<FormFactorBase> clone(const std::string& label) override {
            return nullptr;
        }
        
        
        void eval(const Particle& parent, const ParticleList& daughters, const ParticleList& references) override {

        }
        
    };
    
    class DummyProvidersRepo : public ProvidersRepo {  
    public:

        DummyProvidersRepo() : ProvidersRepo{nullptr} {
        }

        [[nodiscard]] std::vector<FormFactorBase*> getFormFactor(HashId /*processId*/) const override {
            return {_ff};
        }
        
    private:
        
       FormFactorBase* _ff{new DummyFFBtoDX()};
        
    };

    class DummyExternalData : public ExternalData {
    public:

        DummyExternalData() : ExternalData{nullptr, nullptr} {
            _ext = MD::SharedTensorData{
                MD::makeVector({2}, {FF_BD_VAR}, {1.,1.}).release()};
        }

        MD::SharedTensorData getFFEigenVectors(FormFactorBase* ff, const std::string& schemeName) const override {
            return _ext;
        }

    private:

        MD::SharedTensorData _ext;
    };


    class DummyDictionaryManager : public DictionaryManager {

    public:

        DummyDictionaryManager() = default;

        const ProvidersRepo& providers() const override {
            return _prov;
        }

        PurePhaseSpaceDefs& purePSDefs() override {
            return _purePS;
        }

        const PurePhaseSpaceDefs& purePSDefs() const override {
            return _purePS;
        }

        SchemeDefinitions& schemeDefs() override {
            return _schemes;
        }

        const SchemeDefinitions& schemeDefs() const override {
            return _schemes;
        }

        ProcessDefinitions& processDefs() override {
            return _procs;
        }

        const ProcessDefinitions& processDefs() const override {
            return _procs;
        }

        const ExternalData& externalData() const override {
            return _ext;
        }

        ExternalData& externalData() override {
            return _ext;
        }

        // virtual SchemeDict<Tensor>* getProcessRates(HashId) {
        //     return nullptr;
        // }

    private:

        DummyProvidersRepo _prov;
        DummyPurePhaseSpaceDefs _purePS;
        DummySchemeDefinitions _schemes;
        DummyProcessDefinitions _procs;
        DummyExternalData _ext;
    };

    class TestableFFEval : public FFEval {
    public:

        TestableFFEval(DictionaryManager* dict) : FFEval{dict} {
        }

        using FFEval::evalFormFactors;
    };

    TEST(FFEvalTest, WithVar) {
        SettingsHandler sh;
        sh.addSetting<string>("Hammer", "Units", "MeV");

        // Create FFEval; needs dummyDictionaryManager and histo pointer
        DummyDictionaryManager dummyHammer;
        dummyHammer.setSettingsHandler(sh);

        TestableFFEval ffeval(&dummyHammer);
        ffeval.setSettingsHandler(sh);

        auto vals = ffeval.evalFormFactors("Scheme1", 10000, {10.}, {1.});

        EXPECT_EQ(vals.size(), 2);
        
        EXPECT_EQ(vals[0], 2.5);
        EXPECT_EQ(vals[1], 0.5);
    }

} // namespace Hammer
