///
/// @file  TestEvent.cc
/// @brief Tests for Event
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "Hammer/Event.hh"
#include "Hammer/DictionaryManager.hh"
#include "Hammer/PurePhaseSpaceDefs.hh"
#include "Hammer/SchemeDefinitions.hh"
#include "Hammer/ProcessDefinitions.hh"
#include "Hammer/ProvidersRepo.hh"
#include "Hammer/ExternalData.hh"
#include "Hammer/AmplitudeBase.hh"
#include "Hammer/RateBase.hh"
#include "Hammer/Tools/ParticleUtils.hh"
#include "Hammer/Histos.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"
#include "Hammer/Math/MultiDim/VectorContainer.hh"
#include "Hammer/Math/MultiDim/ScalarContainer.hh"
#include "gtest/gtest.h"

#include <boost/functional/hash.hpp>

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    class DummyPurePhaseSpaceDefs : public PurePhaseSpaceDefs {
    public:

        DummyPurePhaseSpaceDefs() = default;

        ~DummyPurePhaseSpaceDefs() override = default;

        [[nodiscard]] NumDenPair<bool> isPurePhaseSpace(PdgId /*parent*/,
                                                        const std::vector<PdgId>& /*daughters*/) const override {
            return {false, false};
        }
    };

    class DummySchemeDefinitions : public SchemeDefinitions {
    public:

        DummySchemeDefinitions() {
            _schemes.insert({"Scheme 1", {}});
            _schemes.insert({"Scheme 2", {}});
        }

        ~DummySchemeDefinitions() override = default;


        [[nodiscard]] std::set<FFIndex> getFormFactorIndices(HashId /*processId*/) const override {
            return {};
        }

        [[nodiscard]] FFIndex getDenominatorFormFactor(HashId /*processId*/) const override {
            return 0;
        }

        [[nodiscard]] const SchemeDict<std::map<HashId, FFIndex>>& getSchemeDefs() const override {
            return _schemes;
        }

        [[nodiscard]] SchemeDict<FFIndex> getFFSchemesForProcess(HashId /*id*/) const override {
            return {};
        }

    private:

        SchemeDict<std::map<HashId, FFIndex>> _schemes;
    };

    class DummyProcessDefinitions : public ProcessDefinitions {
    public:

        DummyProcessDefinitions() = default;

        ~DummyProcessDefinitions() override = default;

        bool isProcessIncluded(std::set<HashId> /*subamplitudes*/, HashId /*processId*/) const override {
            return true;
        }

        void forbidProcess(HashId /*processId*/) const override {
        }
    };

    class DummyProvidersRepo : public ProvidersRepo {
    public:

        DummyProvidersRepo() : ProvidersRepo{nullptr} {
        }

        [[nodiscard]] AmplitudeBase* getAmplitude(PdgId /*parent*/, const std::vector<PdgId>& /*daughters*/,
                                                  const std::vector<PdgId>& /*granddaughters*/) const override {
            return nullptr;
        }

        [[nodiscard]] RateBase* getRate(PdgId /*parent*/, const std::vector<PdgId>& /*daughters*/,
                                        const std::vector<PdgId>& /*granddaughters*/) const override {
            return nullptr;
        }


        [[nodiscard]] std::vector<FormFactorBase*> getFormFactor(HashId /*processId*/) const override {
            return {};
        }
    };

    class DummyExternalData : public ExternalData {
    public:

        DummyExternalData() : ExternalData{nullptr, nullptr} {
            Tensor wc{"wc", MD::makeVector({2}, {WILSON_BCENU}, {2., 1.i})};
            _ext = outerSquare(wc);
        }

        MD::SharedTensorData getWilsonCoefficients(PdgId /*parent*/, const std::vector<PdgId>& /*daughters*/,
                                                   const std::vector<PdgId>& /*granddaughters*/,
                                                   WTerm /*what*/) const override {
            return {};
        }

        MD::SharedTensorData getWilsonCoefficients(const AmplitudeBase* /*amp*/, WTerm /*what*/) const override {
            return {};
        }

        const Tensor& getExternalVectors(const std::string& /*schemeName*/, LabelsList /*labels*/,
                                         const std::string& /*specializationId*/) const override {
            return _ext;
        }

    private:

        Tensor _ext;
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

    class TestableProcManager : public ProcManager {

    public:

        explicit TestableProcManager(const Process& p) : ProcManager(p) {
            // define some dummy content
            Tensor ampNumerator{"Numerator", MD::makeEmptySparse({2, 3, 4}, {SPIN_TAUM, SPIN_GAMMA, WILSON_BCENU})};
            Tensor ampDenominator{"Denominator", MD::makeEmptySparse({2, 3}, {SPIN_TAUM, SPIN_GAMMA})};
            Tensor ampSqNumerator{"Numerator2", MD::makeEmptySparse({4, 4}, {WILSON_BCENU, WILSON_BCENU_HC})};
            Tensor ampSqDenominator{"Denominator2", MD::makeEmptyScalar()};
            Tensor ampWeight1{"Weight1", MD::makeVector({2, 2}, {WILSON_BCENU, WILSON_BCENU_HC}, {0., 1., 1., 0.})};
            Tensor ampWeight2{"Weight2", MD::makeVector({2, 2}, {WILSON_BCENU, WILSON_BCENU_HC}, {1., 1.i, -1.i, -1.})};
            // add some content to the tensors
            ampNumerator.element({0, 2, 0}) = 2.i;
            ampNumerator.element({0, 2, 1}) = 1.;
            ampDenominator.element({0, 2}) = 5.;
            ampDenominator.element({1, 2}) = 2.;
            ampSqNumerator.element({1, 2}) = 1.;
            ampSqNumerator.element({3, 2}) = 7.;
            ampSqDenominator.element({}) = 3.;
            // now fill the content
            modifiableResults().processAmplitude(WTerm::NUMERATOR) = ampNumerator;
            modifiableResults().processAmplitude(WTerm::DENOMINATOR) = ampDenominator;
            modifiableResults().processAmplitude(WTerm::NUMERATOR) = ampSqNumerator;
            modifiableResults().processAmplitude(WTerm::DENOMINATOR) = ampSqDenominator;
            modifiableResults().setProcessWeight("Scheme 1", ampWeight1);
            modifiableResults().setProcessWeight("Scheme 2", ampWeight2);

            modifiableInputs().initialize();
        }

        void initialize(DictionaryManager* /*dictionaries*/) override {
        }

    private:

        ProcResults& modifiableResults() {
            return const_cast<ProcResults&>(results());
        }

        Process& modifiableInputs() {
            return const_cast<Process&>(inputs());
        }
    };

    class TestableEvent : public Event {
    public:

        TestableEvent(Histos* histograms, DictionaryManager* dict) : Event{histograms, dict} {
        }

        ProcessUID addTestProcess(Process& p) {
            p.initialize();
            return addProcManager(TestableProcManager{p});
        }
        using Event::getProcManager;
    };

    TEST(EventTest, Assignment) {
        SettingsHandler sh;
        sh.addSetting<string>("Hammer", "Units", "GeV");

        // Signal process B^0 to D^- tau^+ nu, tau^+ to mu^+ nubar nu
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

        // Tag process B^0bar to (D^*+->D^+pi^0) e^- nubar
        Process procTag;
        procTag.setSettingsHandler(sh);

        Particle pBmesTag{{5280.0000000000000000, 0, 0, 0}, -511};
        Particle pDstarTag{
            {2024.1761928356350379, -52.772533876922062420, 195.31260942369715907, 127.50413291686430435}, 413};
        Particle pEleTag{{1701.1980842908242307, 366.10419906664868130, -191.85093019764307684, -1650.2228930910627331},
                         11};
        Particle pNuBarTag{
            {1554.6257228735407314, -313.33166518972661888, -3.4616792260540822306, 1522.7187601741984288}, -12};
        Particle pPi0Tag{
            {137.69341998653709623, 5.8597010263796357857, -0.94931993861716608795, -26.443157964291685866}, 111};
        Particle pDplTag{{1886.4827728490979417, -58.632234903301698206, 196.26192936231432516, 153.94729088115599022},
                         411};

        auto idxBmesTag = procTag.addParticle(pBmesTag);
        auto idxDstarTag = procTag.addParticle(pDstarTag);
        auto idxEleTag = procTag.addParticle(pEleTag);
        auto idxNuBarTag = procTag.addParticle(pNuBarTag);
        auto idxPi0Tag = procTag.addParticle(pPi0Tag);
        auto idxDplTag = procTag.addParticle(pDplTag);

        procTag.addVertex(idxBmesTag, {idxDstarTag, idxNuBarTag, idxEleTag});
        procTag.addVertex(idxDstarTag, {idxDplTag, idxPi0Tag});

        // Create event; needs dummyDictionaryManager and histo pointer
        DummyDictionaryManager dummyHammer;
        dummyHammer.setSettingsHandler(sh);
        Histos histos;
        histos.setSettingsHandler(sh);

        TestableEvent evt(&histos, &dummyHammer);
        evt.setSettingsHandler(sh);

        auto procId = evt.addTestProcess(proc);
        auto procTagId = evt.addTestProcess(procTag);

        EXPECT_EQ(proc.getId(), procId);
        EXPECT_EQ(procTag.getId(), procTagId);

        HashId idB = processID(511, combineDaughters({-411, 16, -15}, {}));
        HashId idT = processID(-15, combineDaughters({-16, 14, -13}, {}));
        std::set<HashId> fullId;
        fullId.clear();
        fullId.insert(idB);
        fullId.insert(idT);
        HashId id = combineProcessIDs(fullId);

        HashId idBTag = processID(-511, combineDaughters({413, -12, 11}, {}));
        HashId idDTag = processID(413, combineDaughters({411, 111}, {}));
        std::set<HashId> fullIdTag;
        fullIdTag.clear();
        fullIdTag.insert(idBTag);
        fullIdTag.insert(idDTag);
        HashId idTag = combineProcessIDs(fullIdTag);


        // auto procTest = evt.getProcess(id);
        // auto procTagTest = evt.getProcess(idTag);


        // EXPECT_EQ(proc.getId(), procTest.getId());
        // EXPECT_EQ(procTag.getId(), procTagTest.getId());

        // Test weights
        sh.changeSetting<bool>("Process:Weights", false);
        evt.init();

        double evt1wgt = evt.getWeight("Scheme 1", id);
        double evt2wgt = evt.getWeight("Scheme 2", id);

        EXPECT_EQ(evt1wgt, 0.);
        EXPECT_EQ(evt2wgt, 7.);

        auto evt1wgts = evt.getWeights("Scheme 1");
        auto evt2wgts = evt.getWeights("Scheme 2");

        EXPECT_EQ(evt1wgts.find(id)->second, 0.);
        EXPECT_EQ(evt2wgts.find(idTag)->second, 7.);
    }

    TEST(EventTest, Access) {
    }

} // namespace Hammer
