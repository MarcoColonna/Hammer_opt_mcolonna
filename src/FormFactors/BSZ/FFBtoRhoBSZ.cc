///
/// @file  FFBtoRhoBSZ.cc
/// @brief \f$ B \rightarrow \rho \f$ BSZ form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/FormFactors/BSZ/FFBtoRhoBSZ.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFBtoRhoBSZ::FFBtoRhoBSZ() {
        // Create tensor rank and dimensions
        IndexList dims = {8};
        string name{"FFBtoRhoBSZ"};

        setPrefix("BtoRho");
        addProcessSignature(PID::BPLUS, {PID::RHO0});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BRHO})});

        addProcessSignature(PID::BZERO, {PID::RHOMINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BRHO})});

        setSignatureIndex();
    }

    void FFBtoRhoBSZ::defineSettings() {
        //_mFFErrNames = ;
        setPath(getFFErrPrefixGroup().get());
        setUnits("GeV");

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

        setInitialized();
    }

    void FFBtoRhoBSZ::evalAtPSPoint(const vector<double>& point, const vector<double>& masses) {
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

        // Compute and set elements
        // V,A0,A1,A12,T1,T2,T23
        for (size_t i = 0; i < PoleFFvec.size(); ++i) {
            double temp = 0;
            auto vec = PoleFFvec[i].second;
            for (size_t n = 0; n < vec.size(); ++n) {
                temp += vec[n] * zpow[n];
            }
            temp /= PoleFFvec[i].first;

            const auto idx = static_cast<IndexType>(i + 1u); // 0th element is Ap
            result.element({idx}) = temp;
        }
        // Equation of motion for Ap pseudoscalar FF
        result.element({0}) = -2 * Mu / (mb + mu) * result.element({2});
    }

    std::unique_ptr<FormFactorBase> FFBtoRhoBSZ::clone(const std::string& label) {
        MAKE_CLONE(FFBtoRhoBSZ, label);
    }

} // namespace Hammer
