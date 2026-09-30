///
/// @file  TestHistos.cc
/// @brief Tests for Histos
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "Hammer/Histos.hh"
#include "Hammer/DictionaryManager.hh"
#include "Hammer/ExternalData.hh"
#include "Hammer/SchemeDefinitions.hh"
#include "Hammer/Math/Histogram.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Tools/HammerSerial.hh"
#include "Hammer/Math/MultiDim/VectorContainer.hh"
#include "gtest/gtest.h"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    class DummyExternalData : public ExternalData {
    public:

        DummyExternalData()
            : ExternalData{nullptr, nullptr},
              _wcTensors{"WCS", MD::makeVector({3, 3}, {WILSON_BCENU, WILSON_BCENU_HC},
                                               {1., -1.i, 0., 1.i, 1., 1. - 1.i, 0., 1. + 1.i, 0.})} {
        }

        const Tensor& getExternalVectors(const std::string& /*schemeName*/, LabelsList /*labels*/,
                                         const std::string& /*specializationId*/) const override {
            return _wcTensors;
        }

    private:

        Tensor _wcTensors;
    };

    class DummyDictionaryManager : public DictionaryManager {

    public:

        DummyDictionaryManager() = default;

        const ExternalData& externalData() const override {
            return _ext;
        }

        ExternalData& externalData() override {
            return _ext;
        }


    private:

        DummyExternalData _ext;
    };

    TEST(HistosTest, Create) {
        SettingsHandler sh;
        DummyDictionaryManager dm;
        dm.setSettingsHandler(sh);
        Histos histos{&dm};
        histos.setSettingsHandler(sh);

        histos.addHistogramDefinition("Obs123", {10, 20, 5});
        histos.addHistogramDefinition("Obs13", {10, 5}, false);
        histos.addHistogramDefinition("Obs13", IndexList{10});

        PositionList bins1 = {1, 15, 4};
        PositionList bins2 = {2, 4};

        auto bins1Test = histos.getBinIndices("Obs123", {0.7, 14.1, 3.4});
        auto bins2Test = histos.getBinIndices("Obs13", {2.8, 4.4});

        for (size_t idx = 0; idx < 2; ++idx) {
            EXPECT_EQ(bins1Test[idx], bins1[idx]);
        }
        for (size_t idx = 0; idx < 1; ++idx) {
            EXPECT_EQ(bins2Test[idx], bins2[idx]);
        }

        EXPECT_EQ(histos.size(), 2);
    }

    TEST(HistosTest, InitFill) {
        SettingsHandler sh;
        sh.addSetting<string>("Hammer", "Units", "GeV");
        DummyDictionaryManager dm;
        dm.setSettingsHandler(sh);
        dm.schemeDefs().addFFScheme("Mixed", {{"BD", "HQET"}, {"BD*", "CLN"}});

        Histos histos{&dm};
        histos.setSettingsHandler(sh);

        histos.addHistogramDefinition("Obs123", {10, 20, 5}, false);
        histos.addHistogramDefinition("Obs13", {10, 5}, false);

        histos.init();

        // fill Histograms; 2 processes, one event each
        Tensor binTensor1{"Value1", MD::makeVector({3, 3}, {WILSON_BCENU, WILSON_BCENU_HC},
                                                   {1., 2., 1. + 1.i, 2., 4., 2. + 2.i, 1. - 1.i, 2. - 2.i, 2.})};
        Tensor binTensor2{"Value2", MD::makeVector({3, 3}, {WILSON_BCENU, WILSON_BCENU_HC},
                                                   {2., -4., 2. - 2.i, -4., 8., -4. + 4.i, 2. + 2.i, -4. - 4.i, 4.})};

        histos.setEventId({10001, 20000});
        histos.fillHisto("Obs123", Hammer::Spec::none(), "Mixed", {5, 1, 3}, binTensor1);
        histos.fillHisto("Obs13", Hammer::Spec::none(), "Mixed", {5, 3}, binTensor1);

        histos.setEventId({10002, 20000});
        histos.fillHisto("Obs123", Hammer::Spec::none(), "Mixed", {1, 4, 3}, binTensor2);
        histos.fillHisto("Obs13", Hammer::Spec::none(), "Mixed", {1, 3}, binTensor2);

        // test Event map
        auto histoSetTest123 = histos.getHistograms("Obs123", "Mixed", Hammer::Spec::none());
        auto histoSetTest13 = histos.getHistograms("Obs13", "Mixed", Hammer::Spec::none());

        for (size_t idx = 0; idx < 1; ++idx) {
            EXPECT_TRUE(histoSetTest123.find({{10001 + idx, 20000}}) != histoSetTest123.end());
            EXPECT_TRUE(histoSetTest13.find({{10001 + idx, 20000}}) != histoSetTest13.end());
        }

        auto histo123_1 = histoSetTest123.find({{10001, 20000}});
        auto histo123_2 = histoSetTest123.find({{10002, 20000}});
        auto histo13_1 = histoSetTest13.find({{10001, 20000}});
        auto histo13_2 = histoSetTest13.find({{10002, 20000}});
        EXPECT_NEAR(histo123_1->second[20 * 5 * 5 + 5 * 1 + 3].sumWi, 13., 1e-10);
        EXPECT_NEAR(histo123_2->second[20 * 5 * 1 + 5 * 4 + 3].sumWi, 10., 1e-10);
        EXPECT_NEAR(histo13_1->second[5 * 5 + 3].sumWi, 13., 1e-10);
        EXPECT_NEAR(histo13_2->second[5 * 1 + 3].sumWi, 10., 1e-10);

        // retrieve Histogram summed over events
        auto histo123 = histos.getHistogram("Obs123", "Mixed", Hammer::Spec::none());
        auto histo13 = histos.getHistogram("Obs13", "Mixed", Hammer::Spec::none());

        // Tensor readTensor1 = histo123.getWeight(5,1,3);
        // Tensor readTensor2 = histo123.getWeight(1,4,3);
        EXPECT_NEAR(histo123[20 * 5 * 5 + 5 * 1 + 3].sumWi, 13., 1e-10);
        EXPECT_NEAR(histo123[20 * 5 * 1 + 5 * 4 + 3].sumWi, 10., 1e-10);

        // readTensor1 = histo13.getWeight(5,3);
        // readTensor2 = histo13.getWeight(1,3);
        EXPECT_NEAR(histo13[5 * 5 + 3].sumWi, 13., 1e-10);
        EXPECT_NEAR(histo13[5 * 1 + 3].sumWi, 10., 1e-10);
    }

} // namespace Hammer
