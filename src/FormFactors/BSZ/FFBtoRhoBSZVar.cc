///
/// @file  FFBtoRhoBSZVar.cc
/// @brief \f$ B \rightarrow \rho \f$ BSZ form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/FormFactors/BSZ/FFBtoRhoBSZVar.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFBtoRhoBSZVar::FFBtoRhoBSZVar() {
        // Create tensor rank and dimensions
        IndexList dims = {8, 9}; // size is _mFFErrNames + 1 to include central values in zeroth component
        setGroup("BSZVar");              // override BSZ base class FF group setting
        string name{"FFBtoRhoBSZVar"};

        setPrefix("BtoRho");
        _mFFErrLabel = FF_BRHO_VAR;
        addProcessSignature(PID::BPLUS, {PID::RHO0});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BRHO, FF_BRHO_VAR})});

        addProcessSignature(PID::BZERO, {PID::RHOMINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BRHO, FF_BRHO_VAR})});

        setSignatureIndex();
    }

    void FFBtoRhoBSZVar::defineSettings() {
        setPath(getFFErrPrefixGroup().get());
        setUnits("GeV");

        // First 8 principal component directions
        defineAndAddErrSettings(
            {"delta_e1", "delta_e2", "delta_e3", "delta_e4", "delta_e5", "delta_e6", "delta_e7", "delta_e8"});

        addSetting<double>("mb", 4.8);   // GeV
        addSetting<double>("mu", 0.005); // GeV

        // Fit central values from 1503.05534
        vector<double> Vvec = {0.327414, -0.859703, 1.80233};
        vector<double> A0vec = {0.356161, -0.833013, 1.33088};
        vector<double> A1vec = {0.261755, 0.393552, 0.163472};
        vector<double> A12vec = {0.296777, 0.758805, 0.464519};
        vector<double> T1vec = {0.271842, -0.742455, 1.45328};
        vector<double> T2vec = {0.271842, 0.47089, 0.57618};
        vector<double> T23vec = {0.747176, 1.89557, 2.92992};

        addSetting<vector<double>>("Vvec", Vvec);
        addSetting<vector<double>>("A0vec", A0vec);
        addSetting<vector<double>>("A1vec", A1vec);
        addSetting<vector<double>>("A12vec", A12vec);
        addSetting<vector<double>>("T1vec", T1vec);
        addSetting<vector<double>>("T2vec", T2vec);
        addSetting<vector<double>>("T23vec", T23vec);

        addSetting<double>("mR0m", 5.279); // GeV
        addSetting<double>("mR1m", 5.325); // GeV
        addSetting<double>("mR1p", 5.724); // GeV

        // Row basis is V a012,A0 a012,A1 a012,A12 a012,T1 a012,T2 a012,T23 a012
        vector<vector<double>> eigmat{
            {-0.00282109, -0.00288928, -0.0159727, -0.0061639, 0.0100764, 0.00113332, 0.00433743, -0.000146681},
            {0.0227892, -0.0713633, -0.136704, -0.082366, -0.0186232, 0.0153861, 0.0168551, 0.00328266},
            {0.656741, 0.579098, 0.363103, -0.138168, -0.0648946, -0.0223208, -0.0942601, -0.00764847},
            {0.00576933, -0.0115232, 0.00861947, 0.00558892, -0.00570281, 0.0241398, -0.00493472, -0.00244883},
            {0.123107, -0.137028, 0.0356093, -0.00973585, -0.0250626, 0.0486688, -0.0117817, 0.0328683},
            {0.558697, 0.721709, -0.44177, 0.254257, -0.0890678, 0.0274326, -0.00590784, 0.000497942},
            {0., -0.00244184, -0.0137097, -0.00585193, 0.00736384, -0.000511847, -0.00121882, -0.00319619},
            {0.015818, -0.0367692, -0.10235, -0.0677045, -0.00922126, 0.0223016, 0.00379414, -0.00555292},
            {0.129931, -0.0401178, -0.269996, -0.240288, -0.115366, 0.0120762, 0.0156546, -0.0526391},
            {0.00480737, -0.00960191, 0.00718231, 0.00465706, -0.00475197, 0.0201148, -0.00411193, -0.00204048},
            {0.0841703, -0.152548, 0.0260669, 0.0255198, -0.0433473, 0.0646377, -0.0114549, -0.00971878},
            {0.359664, -0.525657, 0.170986, 0.139286, -0.342898, -0.0118739, 0.00383481, 0.000588265},
            {-0.00528002, -0.00511728, -0.0122345, -0.00482501, 0.00776372, -0.00172393, -0.00466985, -0.000720388},
            {0.0200716, -0.0576696, -0.106913, -0.064177, -0.0158466, 0.00567579, -0.0101603, 0.00626081},
            {0.491998, 0.480117, 0.31518, -0.091334, -0.0389457, 0.0548155, 0.113458, 0.0076113},
            {-0.00528002, -0.00511728, -0.0122345, -0.00482501, 0.00776372, -0.00172393, -0.00466985, -0.000720388},
            {0.0128065, -0.0413877, -0.102376, -0.0630881, -0.0107815, 0.0153273, -0.00518128, 0.00396859},
            {0.238221, 0.00794078, -0.289818, -0.248235, -0.0997581, -0.0243646, -0.0078635, 0.0483295},
            {0.01066, -0.00897258, 0.0147762, -0.00393163, 0., 0.0466346, -0.0115295, 0.00412783},
            {0.309366, -0.238442, 0.00322217, -0.0161546, 0.0633549, 0.157124, -0.0387314, -0.00039289},
            {1.75136, -0.411024, -0.0593208, 0.0244927, 0.14957, -0.0462118, 0.0126765, -0.00410159}};

        addSetting<vector<vector<double>>>("eigmatrix", eigmat);

        setInitialized();
    }

    void FFBtoRhoBSZVar::evalAtPSPoint(const vector<double>& point, const vector<double>& masses) {
        Tensor& result = getTensor();
        result.clearData();

        if (!isInitialized()) {
            MSG_WARNING("Warning, Settings have not been defined!");
        }

        double Mb = 0.;
        double Mu = 0.;
        double unitres = 1.;
        tie(Mb, Mu, unitres) = getParentDaughterHadMasses(masses);

        const double Sqq = point[0];
        const double Mb2 = Mb * Mb;
        const double Mu2 = Mu * Mu;
        const double rU = Mu / Mb;

        double w = getW(Sqq, Mb, Mu);
        const double wmax = (Mb2 + Mu2) / (2. * Mb * Mu);
        // safety measure if w==1.0
        if (isZero(w - 1.0)) {
            w += 1e-6;
        }

        // Optimized expansion as in 1503.05534
        const size_t nmax = 4;
        const double z = (pow(rU, 0.25) * sqrt(w + 1) - sqrt(rU + 1)) / (pow(rU, 0.25) * sqrt(w + 1) + sqrt(rU + 1));
        const double z0 =
            (pow(rU, 0.25) * sqrt(wmax + 1) - sqrt(rU + 1)) / (pow(rU, 0.25) * sqrt(wmax + 1) + sqrt(rU + 1));
        vector<double> zpow{1.};
        for (size_t n = 1; n < nmax; ++n) {
            zpow.push_back(zpow[n - 1] * (z - z0));
        }

        // parameters
        const double mb = (*getSetting<double>("mb")) * unitres;
        const double mu = (*getSetting<double>("mu")) * unitres;

        const double mR0m = (*getSetting<double>("mR0m")) * unitres;
        const double mR1m = (*getSetting<double>("mR1m")) * unitres;
        const double mR1p = (*getSetting<double>("mR1p")) * unitres;

        // Poles
        const double P0m = (1 - Sqq / (mR0m * mR0m));
        const double P1m = (1 - Sqq / (mR1m * mR1m));
        const double P1p = (1 - Sqq / (mR1p * mR1p));

        // Parameters.
        vector<pair<double, vector<double>>> PoleFFvec{};
        PoleFFvec.emplace_back(P1m, (*getSetting<vector<double>>("Vvec")));
        PoleFFvec.emplace_back(P0m, (*getSetting<vector<double>>("A0vec")));
        PoleFFvec.emplace_back(P1p, (*getSetting<vector<double>>("A1vec")));
        PoleFFvec.emplace_back(P1p, (*getSetting<vector<double>>("A12vec")));
        PoleFFvec.emplace_back(P1m, (*getSetting<vector<double>>("T1vec")));
        PoleFFvec.emplace_back(P1p, (*getSetting<vector<double>>("T2vec")));
        PoleFFvec.emplace_back(P1p, (*getSetting<vector<double>>("T23vec")));

        // Principal components
        const vector<vector<double>>& eigmat = (*getSetting<vector<vector<double>>>("eigmatrix"));

        // Compute and set elements
        // V,A0,A1,A12,T1,T2,T23
        for (size_t i = 0; i < PoleFFvec.size(); ++i) {
            // Central values
            double temp = 0;
            auto vec = PoleFFvec[i].second;
            for (size_t n = 0; n < vec.size(); ++n) {
                temp += vec[n] * zpow[n];
            }
            temp /= PoleFFvec[i].first;

            const auto idx = static_cast<IndexType>(i + 1u); // 0th element is Ap
            result.element({idx, 0}) = temp;

            // Now compute remaining FF tensor entries
            // Logic: FF^x = a^x_n z^n, a^x_n = M_{nk} delta_k, delta_k = {1, delta_...}
            // n indexes rows for a0 to a2 of each FF
            for (size_t k = 0; k < 8; ++k) {
                double entry = 0;
                for (size_t n = 0; n < 3; ++n) {
                    entry += eigmat[3 * i + n][k] * zpow[n];
                }
                entry /= PoleFFvec[i].first;
                // if(!isZero(entry)) {
                const auto idxp1 = static_cast<IndexType>(k + 1u);
                result.element({idx, idxp1}) = entry;
                //}
            }
        }

        // Equation of motion for Ap pseudoscalar FF
        for (IndexType idx = 0; idx < 9; ++idx) {
            result.element({0, idx}) = -2 * Mu / (mb + mu) * result.element({2, idx});
        }
    }

    std::unique_ptr<FormFactorBase> FFBtoRhoBSZVar::clone(const std::string& label) {
        MAKE_CLONE(FFBtoRhoBSZVar, label);
    }

} // namespace Hammer
