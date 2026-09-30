///
/// @file  FFBtoOmegaBSZVar.cc
/// @brief \f$ B \rightarrow \omega \f$ BSZ form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/FormFactors/BSZ/FFBtoOmegaBSZVar.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFBtoOmegaBSZVar::FFBtoOmegaBSZVar() {
        // Create tensor rank and dimensions
        IndexList dims = {8, 9}; // size is _mFFErrNames + 1 to include central values in zeroth component
        setGroup("BSZVar");              // override BSZ base class FF group setting
        string name{"FFBtoOmegaBSZVar"};

        setPrefix("BtoOmega");
        _mFFErrLabel = FF_BOMEGA_VAR;
        addProcessSignature(PID::BPLUS, {PID::OMEGA});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BOMEGA, FF_BOMEGA_VAR})});

        setSignatureIndex();
    }

    void FFBtoOmegaBSZVar::defineSettings() {
        setPath(getFFErrPrefixGroup().get());
        setUnits("GeV");

        // First 8 principal component directions
        defineAndAddErrSettings(
            {"delta_e1", "delta_e2", "delta_e3", "delta_e4", "delta_e5", "delta_e6", "delta_e7", "delta_e8"});

        addSetting<double>("mb", 4.8);   // GeV
        addSetting<double>("mu", 0.005); // GeV

        // Fit central values from 1503.05534
        vector<double> Vvec = {0.304239, -0.832336, 1.71741};
        vector<double> A0vec = {0.327866, -0.826659, 1.42113};
        vector<double> A1vec = {0.242567, 0.342791, 0.0926204};
        vector<double> A12vec = {0.270411, 0.65589, 0.277567};
        vector<double> T1vec = {0.25137, -0.715314, 1.41463};
        vector<double> T2vec = {0.25137, 0.412936, 0.460141};
        vector<double> T23vec = {0.682582, 1.64618, 2.4695};

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
            {-0.00238092, -0.0170978, 0.019075, -0.00373672, 0.00647176, 0.00208374, 0.00486367, 0.0015144},
            {-0.0175667, -0.18057, 0.145736, -0.165148, -0.0159879, -0.0117861, 0.0105216, -0.0194018},
            {0.474534, 1.06846, -0.313637, -0.246112, -0.000300221, 0.00373613, -0.138488, -0.00572043},
            {0.0189962, -0.00883243, -0.0109694, -0.00211439, -0.00754161, -0.0214784, 0.00150921, 0.00754659},
            {0.249214, -0.126036, -0.048206, -0.0207358, -0.0345583, -0.0660209, -0.00904269, -0.0253452},
            {-0.246756, 0.750601, 0.924939, 0.234163, -0.133554, -0.0464815, -0.00957123, 0.000757252},
            {0., -0.0136197, 0.0147743, -0.00481301, 0.00505258, 0.0017803, 0.000216096, 0.00217686},
            {-0.0284596, -0.138634, 0.124505, -0.133789, 0., -0.0229801, 0.000840447, -0.0125997},
            {-0.0157134, -0.171096, 0.260853, -0.458041, -0.100725, 0.0221401, 0.00240573, 0.0848575},
            {0.0156673, -0.00728463, -0.00904707, -0.00174386, -0.00621993, -0.0177143, 0.00124476, 0.00622414},
            {0.199791, -0.115806, -0.0593745, -0.0100734, -0.0564935, -0.0726485, -0.00223445, 0.00856552},
            {0.792763, -0.206004, -0.265965, 0.124226, -0.456783, 0.0084103, 0.000279959, -0.00536494},
            {-0.0024714, -0.0160236, 0.0131722, -0.0032687, 0.00466392, 0.00237563, -0.00377785, 0.000247321},
            {-0.00904051, -0.141579, 0.111806, -0.127391, -0.0128419, -0.00946817, -0.00752146, -0.0218896},
            {0.338472, 0.884461, -0.25917, -0.134626, 0.00163485, -0.0324706, 0.173197, 0.000806209},
            {-0.0024714, -0.0160236, 0.0131722, -0.0032687, 0.00466392, 0.00237563, -0.00377785, 0.000247321},
            {-0.0276567, -0.13757, 0.12065, -0.124623, -0.00317173, -0.0176447, -0.005348, -0.0204866},
            {0.0449456, -0.125036, 0.298551, -0.452008, -0.06596, 0.0315177, 0.0210719, -0.0614232},
            {0.0346527, -0.0164055, -0.0196033, -0.0111825, 0.000812639, -0.0463261, -0.00195006, 0.0107199},
            {0.518158, -0.253387, -0.0140923, -0.0770334, 0.0682118, -0.196477, -0.0190175, 0.00737968},
            {2.15566, -0.130147, 0.32894, 0.0793229, 0.145875, 0.0573646, 0.00749703, 0.0044929}};

        addSetting<vector<vector<double>>>("eigmatrix", eigmat);

        setInitialized();
    }

    void FFBtoOmegaBSZVar::evalAtPSPoint(const vector<double>& point, const vector<double>& masses) {
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

        double w = (Mb2 + Mu2 - Sqq) / (2. * Mb * Mu);
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

    std::unique_ptr<FormFactorBase> FFBtoOmegaBSZVar::clone(const std::string& label) {
        MAKE_CLONE(FFBtoOmegaBSZVar, label);
    }

} // namespace Hammer
