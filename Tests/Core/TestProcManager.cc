///
/// @file  TestProcManager.cc
/// @brief Tests for ProcManager
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#define FORMFACTORBASE_FRIENDS                  \
    FRIEND_TEST(ProcManagerTest, calc);         \
    FRIEND_TEST(ProcManagerTest, Serialization)


#include "Hammer/ProcManager.hh"
#include "Hammer/DictionaryManager.hh"
#include "Hammer/PurePhaseSpaceDefs.hh"
#include "Hammer/SchemeDefinitions.hh"
#include "Hammer/ProcessDefinitions.hh"
#include "Hammer/SpecializationDefinitions.hh"
#include "Hammer/ProvidersRepo.hh"
#include "Hammer/ExternalData.hh"
#include "Hammer/ProcRates.hh"
#include "Hammer/Amplitudes/AmplBDLepNu.hh"
#include "Hammer/Amplitudes/AmplBDstarDPiLepNu.hh"
#include "Hammer/Amplitudes/AmplTauEllNuNu.hh"
#include "Hammer/Rates/RateBDLepNu.hh"
#include "Hammer/Rates/RateBDstarLepNu.hh"
#include "Hammer/FormFactors/BLPR/FFBtoDBLPR.hh"
#include "Hammer/FormFactors/ISGW2/FFBtoDISGW2.hh"
#include "Hammer/FormFactors/ISGW2/FFBtoDstarISGW2.hh"
#include "Hammer/FormFactors/BLPR/FFBtoDstarBLPR.hh"
#include "Hammer/Math/MultiDim/IContainer.hh"
#include "Hammer/Math/Integrator.hh"
#include "Hammer/Tools/ParticleUtils.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Tools/SettingsConsumer.hh"
#include "Hammer/Tools/HammerSerial.hh"
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
            auto tmpB = combineDaughters({-411}, {});
            HashId idHadB = processID(511, tmpB);
            auto tmpB2 = combineDaughters({-413}, {});
            HashId idHadB2 = processID(511, tmpB2);
            _ffSchemes = {{"Denominator", {{idHadB, 0}, {idHadB2, 0}}},
                          {"ISGW2", {{idHadB, 1}, {idHadB2, 1}}},
                          {"BLPR", {{idHadB, 2}, {idHadB2, 2}}}};
        }

        ~DummySchemeDefinitions() override = default;


        [[nodiscard]] std::set<FFIndex> getFormFactorIndices(HashId processId) const override {
            auto tmpB = combineDaughters({-411}, {});
            HashId idB = processID(511, tmpB);
            // hash (hadronic) Ids always set on positive parent PdgId! (getAmplitude checks for flip signs)
            auto tmpB2 = combineDaughters({-413}, {});
            HashId idB2 = processID(511, tmpB2);
            if (processId == idB) {
                return {0, 1, 2};
            }
            if (processId == idB2) {
                return {0, 1, 2};
            }
            return {};
        }

        [[nodiscard]] FFIndex getDenominatorFormFactor(HashId /*processId*/) const override {
            return 0;
        }

        [[nodiscard]] const SchemeDict<std::map<HashId, FFIndex>>& getSchemeDefs() const override {
            return _ffSchemes;
        }

        [[nodiscard]] SchemeDict<FFIndex> getFFSchemesForProcess(HashId id) const override {
            HashId idB = processID(511, {-411});
            // hash (hadronic) Ids always set on positive parent PdgId! (getAmplitude checks for flip signs)
            HashId idB2 = processID(511, {-413});
            if (id == idB) {
                return {{"Denominator", 0}, {"ISGW2", 1}, {"BLPR", 2}};
            }
            if (id == idB2) {
                return {{"Denominator", 0}, {"ISGW2", 1}, {"BLPR", 2}};
            }
            return {};
        }

    private:

        SchemeDict<std::map<HashId, FFIndex>> _ffSchemes;
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

        explicit DummyProvidersRepo(SettingsHandler& sh)
            : ProvidersRepo{nullptr}, _ampB{new AmplBDLepNu()}, _ampT{new AmplTauEllNuNu()}, _rateB{new RateBDLepNu()},
              _ffBLPR{new FFBtoDBLPR()}, _ffISGW2{new FFBtoDISGW2()}, _ampB2{new AmplBDstarDPiLepNu()},
              _rateB2{new RateBDstarLepNu()}, _ffISGW22{new FFBtoDstarISGW2()}, _ffBLPR2{new FFBtoDstarBLPR()} {
            ProvidersRepo::setSettingsHandler(sh);
            _ffBLPR->setSettingsHandler(sh);
            _ffBLPR->setSignatureIndex(1);
            _ffBLPR->init();
            _ffISGW2->setSettingsHandler(sh);
            _ffISGW2->setSignatureIndex(1);
            _ffISGW2->init();
            _rateB->setSettingsHandler(sh);
            _rateB->setSignatureIndex(1);
            _rateB->init();
            _ampB->setSettingsHandler(sh);
            _ampB->setSignatureIndex(1);
            _ampB->init();
            _ampT->setSettingsHandler(sh);
            _ampT->init();


            _ffISGW22->setSettingsHandler(sh);
            _ffISGW22->setSignatureIndex(1);
            _ffISGW22->init();
            _ffBLPR2->setSettingsHandler(sh);
            _ffBLPR2->setSignatureIndex(1);
            _ffBLPR2->init();
            _rateB2->setSettingsHandler(sh);
            _rateB2->setSignatureIndex(5);
            _rateB2->init();
            _ampB2->setSettingsHandler(sh);
            _ampB2->setSignatureIndex(7);
            _ampB2->init();
        }

        DummyProvidersRepo(const DummyProvidersRepo&) = delete;
        DummyProvidersRepo& operator=(const DummyProvidersRepo&) = delete;
        DummyProvidersRepo(DummyProvidersRepo&&) = delete;
        DummyProvidersRepo& operator=(DummyProvidersRepo&&) = delete;

        [[nodiscard]] AmplitudeBase* getAmplitude(PdgId parent, const std::vector<PdgId>& daughters,
                                                  const std::vector<PdgId>& granddaughters = {}) const override {
            if (parent == 511 && daughters.size() == 3 && granddaughters.empty()) {
                return _ampB.get();
            }
            // getAmplitude checks for flip signs, here just hard code it
            if (parent == -511 && daughters.size() == 3 && granddaughters.size() == 2) {
                return _ampB2.get();
            }
            if (parent == -15 && daughters.size() == 3 && granddaughters.empty()) {
                return _ampT.get();
            }
            return nullptr;
        }

        [[nodiscard]] RateBase* getRate(PdgId parent, const std::vector<PdgId>& /*daughters*/,
                                        const std::vector<PdgId>& /*granddaughters*/) const override {
            //            auto tmpB = combineDaughters({-411, 16, -15}, {});
            //            HashId idB = processID(511, tmpB);
            //            auto tmpB2 = combineDaughters({413, -12, 11}, {411, 111});
            //            HashId idB2 = processID(-511, tmpB2);
            if (parent == 511) {
                return _rateB.get();
            }
            if (parent == -511) {
                return _rateB2.get();
            }
            return nullptr;
        }

        [[nodiscard]] const double* getPartialWidth(PdgId /*parent*/, const std::vector<PdgId>& /*daughters*/,
                                                    const std::vector<PdgId>& /*granddaughters*/ = {}) const override {
            return _pw;
        }

        [[nodiscard]] std::vector<FormFactorBase*> getFormFactor(HashId processId) const override {
            auto tmpB = combineDaughters({-411}, {});
            HashId idB = processID(511, tmpB);
            auto tmpB2 = combineDaughters({-413}, {});
            HashId idB2 = processID(511, tmpB2);
            if (processId == idB) {
                return {_ffISGW2.get(), _ffISGW2.get(), _ffBLPR.get()};
            }
            if (processId == idB2) {
                return {_ffISGW22.get(), _ffISGW22.get(), _ffBLPR2.get()};
            }
            return {};
        }

        AmplitudeBase* getAmpB(bool isSecond) {
            return isSecond ? _ampB2.get() : _ampB.get();
        }

    private:

        unique_ptr<AmplitudeBase> _ampB;
        unique_ptr<AmplitudeBase> _ampT;
        unique_ptr<RateBase> _rateB;
        unique_ptr<FormFactorBase> _ffBLPR;
        unique_ptr<FormFactorBase> _ffISGW2;
        unique_ptr<AmplitudeBase> _ampB2;
        unique_ptr<RateBase> _rateB2;
        unique_ptr<FormFactorBase> _ffISGW22;
        unique_ptr<FormFactorBase> _ffBLPR2;
        double pw = 1.;
        const double* _pw = &pw;
    };

    class DummyExternalData : public ExternalData {
    public:

        explicit DummyExternalData(DummyProvidersRepo* repo) : ExternalData{nullptr, nullptr}, _repo{repo} {
            const AmplitudeBase* ampB = _repo->getAmpB(false);
            const AmplitudeBase* ampB2 = _repo->getAmpB(true);
            ampB->updateWCTensor(ampB->getWCVectorFromSettings(WTerm::NUMERATOR), _wcTensors);
            ampB2->updateWCTensor(ampB2->getWCVectorFromSettings(WTerm::NUMERATOR), _wcTensors2);
        }

        DummyExternalData(const DummyExternalData&) = delete;
        DummyExternalData& operator=(const DummyExternalData&) = delete;
        DummyExternalData(DummyExternalData&&) = delete;
        DummyExternalData& operator=(DummyExternalData&&) = delete;

        MD::SharedTensorData getWilsonCoefficients(PdgId parent, const std::vector<PdgId>& /*daughters*/,
                                                   const std::vector<PdgId>& /*granddaughters*/,
                                                   WTerm /*what*/) const override {
            if (parent == 511) {
                return _wcTensors;
            }
            if (parent == -511) {
                return _wcTensors2;
            }
            return nullptr;
        }

        MD::SharedTensorData getWilsonCoefficients(const AmplitudeBase* amp, WTerm /*what*/) const override {
            if (amp == _repo->getAmpB(false)) {
                return _wcTensors;
            }
            if (amp == _repo->getAmpB(true)) {
                return _wcTensors2;
            }
            return nullptr;
        }

    private:

        MD::SharedTensorData _wcTensors;
        MD::SharedTensorData _wcTensors2;
        DummyProvidersRepo* _repo;
    };

    class DummyDictionaryManager : public DictionaryManager {

    public:

        explicit DummyDictionaryManager(SettingsHandler& sh) : _prov{sh}, _ext{&_prov}, _specs{&_prov} {
        }

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

        const SpecializationDefinitions& specializationDefs() const override {
            return _specs;
        }

        SpecializationDefinitions& specializationDefs() override {
            return _specs;
        }


    private:

        DummyProvidersRepo _prov;
        DummyPurePhaseSpaceDefs _purePS;
        DummySchemeDefinitions _schemes;
        DummyProcessDefinitions _procs;
        DummyExternalData _ext;
        SpecializationDefinitions _specs;
    };

    TEST(ProcManagerTest, calc) {
        SettingsHandler sh;
        sh.addSetting<string>("Hammer", "Units", "MeV");
        sh.addSetting<bool>("Hammer", "CalcGeneralWeights", true);
        sh.addSetting<bool>("Hammer", "GeneralRates", true);

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

        DummyDictionaryManager dummyHammer{sh};
        dummyHammer.setSettingsHandler(sh);
        ProcManager procm{proc};
        procm.setSettingsHandler(sh);
        procm.initialize(&dummyHammer);
        procm.calc();
        Tensor procAmpl = procm.results().processAmplitude();

        // Direct evaluation of amplitude
        AmplBDLepNu amplB;
        amplB.eval(pBmes, {pDmes, pNuTau, pTau}, {});
        AmplTauEllNuNu amplTau;
        amplTau.eval(pTau, {pNuBarTau, pNuMuon, pMuon}, {pDmes, pNuTau});
        Tensor tB = amplB.getTensor();
        tB.dot(amplTau.getTensor());


        // Check TEs.
        for (IndexType idx1 = 0; idx1 < 11; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 4; ++idx2) {
                for (IndexType idx3 = 0; idx3 < 2; ++idx3) {
                    double comRe =
                        compareVals(procAmpl.element({idx1, idx2, idx3}).real(), tB.element({idx1, idx2, idx3}).real());
                    double comIm =
                        compareVals(procAmpl.element({idx1, idx2, idx3}).imag(), tB.element({idx1, idx2, idx3}).imag());
                    EXPECT_NEAR(comRe, 1., 1e-4);
                    EXPECT_NEAR(comIm, 1., 1e-4);
                }
            }
        }

        // Test FFs
        //  auto tmpB = combineDaughters({-411}, {});
        //  HashId idHadB = processID(511, tmpB);
        Tensor procFFISGW2 = procm.results().processFormFactors("ISGW2")[0];
        Tensor procFFBLPR = procm.results().processFormFactors("BLPR")[0];

        FFBtoDISGW2 ffISGW2;
        ffISGW2.setSettingsHandler(sh);
        ffISGW2.initSettings();
        ffISGW2.setSignatureIndex(1);
        ffISGW2.eval(pBmes, {pDmes, pNuTau, pTau}, {});
        Tensor tISGW2 = ffISGW2.getTensor();

        FFBtoDBLPR ffBLPR;
        ffBLPR.setSettingsHandler(sh);
        ffBLPR.initSettings();
        ffBLPR.setSignatureIndex(1);
        ffBLPR.eval(pBmes, {pDmes, pNuTau, pTau}, {});
        Tensor tBLPR = ffBLPR.getTensor();


        // Check TEs.
        for (IndexType idx1 = 0; idx1 < 4; ++idx1) {
            double comRe = compareVals(procFFISGW2.element({idx1}).real(), tISGW2.element({idx1}).real());
            double comIm = compareVals(procFFISGW2.element({idx1}).imag(), tISGW2.element({idx1}).imag());
            EXPECT_NEAR(comRe, 1., 1e-4);
            EXPECT_NEAR(comIm, 1., 1e-4);
        }

        // Check TEs.
        for (IndexType idx1 = 0; idx1 < 4; ++idx1) {
            double comRe = compareVals(procFFBLPR.element({idx1}).real(), tBLPR.element({idx1}).real());
            double comIm = compareVals(procFFBLPR.element({idx1}).imag(), tBLPR.element({idx1}).imag());
            EXPECT_NEAR(comRe, 1., 1e-4);
            EXPECT_NEAR(comIm, 1., 1e-4);
        }

        // Check B process rates
        auto vertexId = procm.inputs().getVertexId("BDTauNu");
        Tensor& procRateISGW2 = dummyHammer.rates().getProcessRates(vertexId)->find("")->second.find("ISGW2")->second;
        Tensor& procRateBLPR = dummyHammer.rates().getProcessRates(vertexId)->find("")->second.find("BLPR")->second;
        EXPECT_FALSE(dummyHammer.rates().getProcessRates(vertexId)->find("") ==
                     dummyHammer.rates().getProcessRates(vertexId)->end());
        //        Tensor& procRateISGW2 = dummyHammer._rates.find("ISGW2")->second;
        //        Tensor& procRateBLPR = dummyHammer._rates.find("BLPR")->second;

        RateBDLepNu rate;
        // set correct signature index, B^0 -> D^- tau nu
        rate.setSignatureIndex(1);
        rate.init();
        rate.calcTensor();
        const auto& points = rate.getEvaluationPoints();
        Tensor integISGW2 = ffISGW2.getFFPSIntegrand(points);
        Tensor integBLPR = ffBLPR.getFFPSIntegrand(points);
        // Integrate
        Tensor rateT = rate.getTensor();
        integISGW2.dot(rateT);
        integBLPR.dot(rateT);
        //        double denRateVal = integISGW2.element().real();
        //        integISGW2 *= 1./denRateVal;
        //        integBLPR *= 1./denRateVal;

        // Check TEs.
        for (IndexType idx1 = 0; idx1 < 11; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 11; ++idx2) {
                double comRe =
                    compareVals(integISGW2.element({idx1, idx2}).real(), procRateISGW2.element({idx1, idx2}).real());
                double comIm =
                    compareVals(integISGW2.element({idx1, idx2}).imag(), procRateISGW2.element({idx1, idx2}).imag());
                EXPECT_NEAR(comRe, 1., 1e-4);
                EXPECT_NEAR(comIm, 1., 1e-4);
            }
        }

        // Check TEs.
        for (IndexType idx1 = 0; idx1 < 11; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 11; ++idx2) {
                double comRe =
                    compareVals(integBLPR.element({idx1, idx2}).real(), procRateBLPR.element({idx1, idx2}).real());
                double comIm =
                    compareVals(integBLPR.element({idx1, idx2}).imag(), procRateBLPR.element({idx1, idx2}).imag());
                EXPECT_NEAR(comRe, 1., 1e-4);
                EXPECT_NEAR(comIm, 1., 1e-4);
            }
        }


        // Test Weights
        const Tensor& procWgtISGW2 = procm.results().processWeight("ISGW2");
        const Tensor& procWgtBLPR = procm.results().processWeight("BLPR");

        Tensor wgt = outerSquare(tB);
        wgt.spinSum();
        tISGW2.outerSquare();
        tBLPR.outerSquare();
        Tensor wgtISGW2 = dot(wgt, tISGW2);
        Tensor wgtBLPR = dot(wgt, tBLPR);

        double denVal = wgtISGW2.element({0, 0}).real();
        wgtISGW2 *= 1. / denVal;
        wgtBLPR *= 1. / denVal;


        // Check TEs.
        for (IndexType idx1 = 0; idx1 < 11; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 11; ++idx2) {
                double comRe =
                    compareVals(procWgtISGW2.element({idx1, idx2}).real(), wgtISGW2.element({idx1, idx2}).real());
                double comIm =
                    compareVals(procWgtISGW2.element({idx1, idx2}).imag(), wgtISGW2.element({idx1, idx2}).imag());
                EXPECT_NEAR(comRe, 1., 1e-4);
                EXPECT_NEAR(comIm, 1., 1e-4);
            }
        }


        // Check TEs.
        for (IndexType idx1 = 0; idx1 < 11; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 11; ++idx2) {
                double comRe =
                    compareVals(procWgtBLPR.element({idx1, idx2}).real(), wgtBLPR.element({idx1, idx2}).real());
                double comIm =
                    compareVals(procWgtBLPR.element({idx1, idx2}).imag(), wgtBLPR.element({idx1, idx2}).imag());
                EXPECT_NEAR(comRe, 1., 1e-4);
                EXPECT_NEAR(comIm, 1., 1e-4);
            }
        }
    }

    TEST(ProcManagerTest, Edge) {
        SettingsHandler sh;
        sh.addSetting<string>("Hammer", "Units", "MeV");
        sh.addSetting<bool>("Hammer", "CalcGeneralWeights", true);
        sh.addSetting<bool>("Hammer", "GeneralRates", true);

        Process proc;
        proc.setSettingsHandler(sh);

        // Process B^0bar to (D^*+->D^+pi^0) e^- nubar
        Particle pBmes{{5280.0000000000000000, 0., 0., 0.}, -511};
        Particle pDstar{{2024.1761928356350379, -52.772533876922062420, 195.31260942369715907, 127.50413291686430435},
                        413};
        Particle pEle{{1701.1980842908242307, 366.10419906664868130, -191.85093019764307684, -1650.2228930910627331},
                      11};
        Particle pNuBar{{1554.6257228735407314, -313.33166518972661888, -3.4616792260540822306, 1522.7187601741984288},
                        -12};
        Particle pPi0{{137.69341998653709623, 5.8597010263796357857, -0.94931993861716608795, -26.443157964291685866},
                      111};
        Particle pDpl{{1886.4827728490979417, -58.632234903301698206, 196.26192936231432516, 153.94729088115599022},
                      411};

        auto idxBmes = proc.addParticle(pBmes);
        auto idxDstar = proc.addParticle(pDstar);
        auto idxEle = proc.addParticle(pEle);
        auto idxNuBar = proc.addParticle(pNuBar);
        auto idxPi0 = proc.addParticle(pPi0);
        auto idxDpl = proc.addParticle(pDpl);

        proc.addVertex(idxBmes, {idxDstar, idxNuBar, idxEle});
        proc.addVertex(idxDstar, {idxDpl, idxPi0});

        DummyDictionaryManager dummyHammer{sh};
        dummyHammer.setSettingsHandler(sh);
        ProcManager procm{proc};
        procm.setSettingsHandler(sh);
        procm.initialize(&dummyHammer);
        procm.calc();

        // Direct evaluation of amplitude
        Tensor procAmpl = procm.results().processAmplitude();
        AmplBDstarDPiLepNu amplB;
        amplB.eval(pBmes, {pDstar, pNuBar, pEle, pDpl, pPi0}, {});

        // Check TEs.
        for (IndexType idx1 = 0; idx1 < 2; ++idx1) {
            for (IndexType idx3 = 0; idx3 < 2; ++idx3) {
                for (IndexType idx4 = 0; idx4 < 11; ++idx4) {
                    for (IndexType idx5 = 0; idx5 < 8; ++idx5) {
                        const Tensor& tB = amplB.getTensor();
                        double comRe = compareVals(procAmpl.element({idx4, idx5, idx1, idx1, idx3}).real(),
                                                   tB.element({idx4, idx5, idx1, idx1, idx3}).real());
                        double comIm = compareVals(procAmpl.element({idx4, idx5, idx1, idx1, idx3}).imag(),
                                                   tB.element({idx4, idx5, idx1, idx1, idx3}).imag());
                        EXPECT_NEAR(comRe, 1., 1e-4);
                        EXPECT_NEAR(comIm, 1., 1e-4);
                    }
                }
            }
        }
    }

    TEST(ProcManagerTest, PHOTOS) {
        SettingsHandler sh;
        sh.addSetting<string>("Hammer", "Units", "GeV");
        sh.addSetting<bool>("Hammer", "CalcGeneralWeights", true);
        sh.addSetting<bool>("Hammer", "GeneralRates", true);
        Process procR;
        procR.setSettingsHandler(sh);

        // Test rad photon expander
        // Process B^0 to D^- e^+ nu gamma
        Particle pBmes{
            {5.2907140549134599e+00, 1.4076777041327704e-01, -2.8889170398557235e-01, -1.2006770023620371e-01}, 511};
        Particle pD{{2.9445868706022504e+00, -1.2787797729688957e+00, -3.7348260744627737e-01, -1.8440128309217103e+00},
                    -411};
        Particle pEle{{1.6950811161670416e+00, 1.1160408238116024e+00, 3.6598183916406024e-01, 1.2222151620128943e+00},
                      -11};
        Particle pNu{{6.4827796281968653e-01, 3.0168368427732462e-01, -2.8198754968283246e-01, 4.9973422290490982e-01},
                     12};
        Particle pGamma{
            {2.7681053249831080e-03, 1.8230352932592250e-03, 5.9661397944993388e-04, 1.9957457676908650e-03}, 22};

        auto idxB = procR.addParticle(pBmes);
        auto idxD = procR.addParticle(pD);
        auto idxEle = procR.addParticle(pEle);
        auto idxGamma = procR.addParticle(pGamma);
        auto idxNu = procR.addParticle(pNu);

        procR.addVertex(idxB, {idxD, idxNu, idxEle, idxGamma});

        DummyDictionaryManager dummyHammer{sh};
        dummyHammer.setSettingsHandler(sh);
        ProcManager procmR{procR};
        procmR.setSettingsHandler(sh);
        procmR.initialize(&dummyHammer);
        procmR.calc();

        const Process& procRnew = procmR.inputs();
        const_cast<Process&>(procRnew).setSettingsHandler(sh);
        auto pTot = procRnew.getParticle(idxB).p() - procRnew.getParticle(idxD).p() - procRnew.getParticle(idxEle).p() -
                    procRnew.getParticle(idxNu).p();
        EXPECT_LT(pTot.E(), 1e-8);
        EXPECT_LT(pTot.px(), 1e-8);
        EXPECT_LT(pTot.py(), 1e-8);
        EXPECT_LT(pTot.pz(), 1e-8);


        Tensor procRAmpl = procmR.results().processAmplitude();
        AmplBDLepNu amplBR;
        // amplBR.eval(pBmes, {pD, pNu, pEle}, {});
        amplBR.eval(procR.getParticle(idxB),
                    {procR.getParticle(idxD), procR.getParticle(idxNu), procR.getParticle(idxEle)}, {});
        // Check TEs.
        for (IndexType idx1 = 0; idx1 < 2; ++idx1) {
            for (IndexType idx3 = 0; idx3 < 2; ++idx3) {
                for (IndexType idx4 = 0; idx4 < 11; ++idx4) {
                    for (IndexType idx5 = 0; idx5 < 4; ++idx5) {
                        const Tensor& tBR = amplBR.getTensor();
                        double comRe = compareVals(procRAmpl.element({idx4, idx5, idx1, idx1, idx3}).real(),
                                                   tBR.element({idx4, idx5, idx1, idx1, idx3}).real());
                        double comIm = compareVals(procRAmpl.element({idx4, idx5, idx1, idx1, idx3}).imag(),
                                                   tBR.element({idx4, idx5, idx1, idx1, idx3}).imag());
                        EXPECT_NEAR(comRe, 1., 1e-4);
                        EXPECT_NEAR(comIm, 1., 1e-4);
                    }
                }
            }
        }
    }

    TEST(ProcManagerTest, Serialization) {
        SettingsHandler sh;
        sh.addSetting<string>("Hammer", "Units", "MeV");
        sh.addSetting<bool>("Hammer", "CalcGeneralWeights", true);
        sh.addSetting<bool>("Hammer", "SaveGeneralWeights", true);
        sh.addSetting<bool>("Hammer", "GeneralRates", true);
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

        DummyDictionaryManager dummyHammer{sh};
        dummyHammer.setSettingsHandler(sh);
        ProcManager procm{proc};
        procm.setSettingsHandler(sh);
        procm.initialize(&dummyHammer);
        procm.calc();

        auto builder = make_unique<flatbuffers::FlatBufferBuilder>();
        flatbuffers::Offset<Serial::FBProcess> msg;
        builder->Clear();
        procm.write(builder.get(), &msg);
        builder->Finish(msg);
        uint8_t* buf = builder->GetBufferPointer();
        // unsigned int size = builder->GetSize();
        const auto* ser1 = flatbuffers::GetRoot<Serial::FBProcess>(buf);
        EXPECT_NE(ser1, nullptr);
        if (ser1 != nullptr) {
            ProcManager procm2;
            procm2.initialize(&dummyHammer);
            cout << "AA" << '\n';
            procm2.read(ser1, false);
            cout << "BB" << '\n';
            Tensor procAmpl = procm2.results().processAmplitude();

            // Direct evaluation of amplitude
            AmplBDLepNu amplB;
            amplB.eval(pBmes, {pDmes, pNuTau, pTau}, {});
            AmplTauEllNuNu amplTau;
            amplTau.eval(pTau, {pNuBarTau, pNuMuon, pMuon}, {pDmes, pNuTau});
            Tensor tB = amplB.getTensor();
            tB.dot(amplTau.getTensor());

            cout << "CC" << '\n';

            // Check TEs.
            for (IndexType idx1 = 0; idx1 < 2; ++idx1) {
                for (IndexType idx2 = 0; idx2 < 11; ++idx2) {
                    for (IndexType idx3 = 0; idx3 < 4; ++idx3) {
                        double comRe = compareVals(procAmpl.element({idx2, idx3, idx1}).real(),
                                                   tB.element({idx2, idx3, idx1}).real());
                        double comIm = compareVals(procAmpl.element({idx2, idx3, idx1}).imag(),
                                                   tB.element({idx2, idx3, idx1}).imag());
                        EXPECT_NEAR(comRe, 1., 1e-4);
                        EXPECT_NEAR(comIm, 1., 1e-4);
                    }
                }
            }
            cout << "DD" << '\n';

            // Test FFs
            auto tmpB = combineDaughters({-411}, {});
            // HashId idHadB = processID(511, tmpB);
            Tensor procFFISGW2 = procm2.results().processFormFactors("ISGW2")[0];
            Tensor procFFBLPR = procm2.results().processFormFactors("BLPR")[0];

            FFBtoDISGW2 ffISGW2;
            ffISGW2.setSettingsHandler(sh);
            ffISGW2.initSettings();
            ffISGW2.setSignatureIndex(1);
            ffISGW2.eval(pBmes, {pDmes, pNuTau, pTau}, {});
            Tensor tISGW2 = ffISGW2.getTensor();

            FFBtoDBLPR ffBLPR;
            ffBLPR.setSettingsHandler(sh);
            ffBLPR.initSettings();
            ffBLPR.setSignatureIndex(1);
            ffBLPR.eval(pBmes, {pDmes, pNuTau, pTau}, {});
            Tensor tBLPR = ffBLPR.getTensor();

            cout << "EE" << '\n';
            // Check TEs.
            for (IndexType idx1 = 0; idx1 < 4; ++idx1) {
                double comRe = compareVals(procFFISGW2.element({idx1}).real(), tISGW2.element({idx1}).real());
                double comIm = compareVals(procFFISGW2.element({idx1}).imag(), tISGW2.element({idx1}).imag());
                EXPECT_NEAR(comRe, 1., 1e-4);
                EXPECT_NEAR(comIm, 1., 1e-4);
            }

            // Check TEs.
            for (IndexType idx1 = 0; idx1 < 4; ++idx1) {
                double comRe = compareVals(procFFBLPR.element({idx1}).real(), tBLPR.element({idx1}).real());
                double comIm = compareVals(procFFBLPR.element({idx1}).imag(), tBLPR.element({idx1}).imag());
                EXPECT_NEAR(comRe, 1., 1e-4);
                EXPECT_NEAR(comIm, 1., 1e-4);
            }

            cout << "FF" << '\n';
            // Test Weights
            Tensor procWgtISGW2 = procm2.results().processWeight("ISGW2");
            Tensor procWgtBLPR = procm2.results().processWeight("BLPR");

            cout << "GG" << '\n';
            Tensor wgt = outerSquare(tB);
            wgt.spinSum();
            tISGW2.outerSquare();
            tBLPR.outerSquare();
            cout << "HH" << '\n';
            Tensor wgtISGW2 = dot(wgt, tISGW2);
            Tensor wgtBLPR = dot(wgt, tBLPR);

            cout << "II" << '\n';
            double denVal = wgtISGW2.element({0, 0}).real();
            wgtISGW2 *= 1. / denVal;
            wgtBLPR *= 1. / denVal;

            // Check TEs.
            for (IndexType idx1 = 0; idx1 < 11; ++idx1) {
                for (IndexType idx2 = 0; idx2 < 11; ++idx2) {
                    double comRe =
                        compareVals(procWgtISGW2.element({idx1, idx2}).real(), wgtISGW2.element({idx1, idx2}).real());
                    double comIm =
                        compareVals(procWgtISGW2.element({idx1, idx2}).imag(), wgtISGW2.element({idx1, idx2}).imag());
                    EXPECT_NEAR(comRe, 1., 1e-4);
                    EXPECT_NEAR(comIm, 1., 1e-4);
                }
            }

            // Check TEs.
            for (IndexType idx1 = 0; idx1 < 11; ++idx1) {
                for (IndexType idx2 = 0; idx2 < 11; ++idx2) {
                    double comRe =
                        compareVals(procWgtBLPR.element({idx1, idx2}).real(), wgtBLPR.element({idx1, idx2}).real());
                    double comIm =
                        compareVals(procWgtBLPR.element({idx1, idx2}).imag(), wgtBLPR.element({idx1, idx2}).imag());
                    EXPECT_NEAR(comRe, 1., 1e-4);
                    EXPECT_NEAR(comIm, 1., 1e-4);
                }
            }
        }
    }

} // namespace Hammer
