///
/// @file  TestAmplBToQLepNuBase.cc
/// @brief Tests for AmplBToQLepNuBase
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "Hammer/Amplitudes/AmplBToQLepNuBase.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Math/MultiDim/VectorContainer.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"
#include "gtest/gtest.h"
// May need other stuff here

using namespace std;
using namespace complex_literals;

namespace Hammer {

    namespace MD = MultiDimensional;


    class TestableAmplBToQLepNuBase : public AmplBToQLepNuBase {

    public:

        TestableAmplBToQLepNuBase() = default;

        using AmplBToQLepNuBase::addTensor;
        using AmplBToQLepNuBase::defineSettings;
        using AmplBToQLepNuBase::updateWilsonCoeffLabelPrefix;

        void eval(const Particle& /*parent*/, const ParticleList& /*daughters*/,
                  const ParticleList& /*references*/) override {
            Tensor& t = getTensor();
            t.clearData();
        }
    };


    TEST(AmplBToQLepNuBaseTest, Mapping) {
        SettingsHandler sh;
        TestableAmplBToQLepNuBase ampl;
        ampl.setSettingsHandler(sh);

        vector<IndexType> dims{{11, 8, 3, 2, 2, 2}};
        ampl.addTensor(Tensor{"AmplBDstarLepNu", MD::makeEmptySparse(dims, {WILSON_BCTAUNU, FF_BDSTAR, SPIN_DSTAR,
                                                                            SPIN_NUTAU, SPIN_NUTAU_REF, SPIN_TAUP})});
        ampl.updateWilsonCoeffLabelPrefix();
        ampl.defineSettings();


        map<string, complex<double>> dict{{{"SM", 1},
                                           {"S_qLlL", 1i},
                                           {"S_qRlL", 2i},
                                           {"V_qLlL", 3i},
                                           {"V_qRlL", 4i},
                                           {"T_qLlL", 5i},
                                           {"S_qLlR", 6i},
                                           {"S_qRlR", 7i},
                                           {"V_qLlR", 8i},
                                           {"V_qRlR", 9i},
                                           {"T_qRlR", 10i}}};
        auto ovec = ampl.getWCVectorFromDict(dict);
        ampl.preProcessWCValues(ovec);

        vector<complex<double>> ovec2{1, -(-1i), -(-6i), -(-2i), -(-7i), (-4i), (-9i), (-3i), (-8i), -(-5i), -(-10i)};
        // Check
        for (IndexType idx1 = 0; idx1 < 11; ++idx1) {
            double comRe = compareVals(ovec[idx1].real(), ovec2[idx1].real());
            double comIm = compareVals(ovec[idx1].imag(), ovec2[idx1].imag());
            EXPECT_DOUBLE_EQ(comRe, 1.);
            EXPECT_DOUBLE_EQ(comIm, 1.);
        }


        ampl.preProcessWCValues(ovec, true);
        auto dict2 = ampl.getDictFromWCVector(ovec);

        for (auto& elem : dict) {
            auto it = dict2.find(elem.first);
            double comRe = compareVals(elem.second.real(), it->second.real());
            double comIm = compareVals(elem.second.imag(), it->second.imag());
            EXPECT_DOUBLE_EQ(comRe, 1.);
            EXPECT_DOUBLE_EQ(comIm, 1.);
        }
    }

    TEST(AmplBToQLepNuBaseTest, proj) {
        SettingsHandler sh;
        TestableAmplBToQLepNuBase ampl;
        ampl.setSettingsHandler(sh);

        vector<IndexType> dims{{11, 8, 3, 2, 2, 2}};
        ampl.addTensor(Tensor{"AmplBDstarLepNu", MD::makeEmptySparse(dims, {WILSON_BCTAUNU, FF_BDSTAR, SPIN_DSTAR,
                                                                            SPIN_NUTAU, SPIN_NUTAU_REF, SPIN_TAUP})});
        ampl.updateWilsonCoeffLabelPrefix();
        ampl.defineSettings();

        vector<map<string, complex<double>>> subspace{{{{"SM", 1}}, {{"T_qLlL", 4}, {"S_qRlL", 1}, {"X", 4}}}};
        map<string, complex<double>> origin{{{"SM", 1}, {"S_qRlL", 3}}};
        auto ovec = ampl.getWCVectorFromDict(origin);
        ampl.preProcessWCValues(ovec);
        auto otensor = MD::SharedTensorData{MD::makeVector({11}, {WILSON_BCTAUNU}, ovec)};
        auto proj = MD::SharedTensorData{};

        ampl.createWCProjectionTensor(subspace, specializeLabel(WILSON_BCTAUNU), otensor, proj);

        // Compare to diract projector
        Tensor projDir{"proj", MD::makeVector({11, 3}, {WILSON_BCTAUNU, specializeLabel(WILSON_BCTAUNU)},
                                              {1, 1, 0, 0, 0, 0, 0, 0, 0, -3, 0, -1, 0,  0, 0, 0, 0,
                                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,  -4, 0, 0, 0})};

        // Check
        for (IndexType idx1 = 0; idx1 < 11; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 3; ++idx2) {
                double comRe = compareVals(proj->element({idx1, idx2}).real(), projDir.element({idx1, idx2}).real());
                EXPECT_NEAR(comRe, 1., 1e-10);
                // cout << proj->element({idx1,idx2}).real() << " ";
            }
        }
    }

} // namespace Hammer
