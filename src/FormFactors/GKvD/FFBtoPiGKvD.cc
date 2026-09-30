///
/// @file  FFBtoPiGKvD.cc
/// @brief \f$ B \rightarrow \pi \f$ GKvD form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/FormFactors/GKvD/FFBtoPiGKvD.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFBtoPiGKvD::FFBtoPiGKvD() {
        // Create tensor rank and dimensions
        IndexList dims = {4};
        setGroup("GKvD"); // override BSZ base class FF group setting
        string name{"FFBtoPiGKvD"};

        setPrefix("BtoPi");
        addProcessSignature(PID::BPLUS, {PID::PI0});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BPI})});

        addProcessSignature(PID::BZERO, {PID::PIMINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BPI})});

        setSignatureIndex();
    }

    void FFBtoPiGKvD::defineSettings() {
        //_mFFErrNames = ;
        setPath(getFFErrPrefixGroup().get());
        setUnits("GeV");

        addSetting<double>("mb", 4.710); // GeV
        addSetting<double>("mu", 0.005); // GeV

        // Fit central values from 1811.00983
        vector<double> Fzvec = {-0.014468, -0.195297};
        vector<double> Fpvec = {0.19595, -0.776357, -0.196605};
        vector<double> Ftvec = {0.173117, -0.409037, -0.258036};

        addSetting<vector<double>>("Fzvec", Fzvec);
        addSetting<vector<double>>("Fpvec", Fpvec);
        addSetting<vector<double>>("Ftvec", Ftvec);

        addSetting<double>("mR0p", 5.540); // GeV
        addSetting<double>("mR1m", 5.325); // GeV
        // addSetting<double>("mR1p", 5.724); //GeV

        setInitialized();
    }

    void FFBtoPiGKvD::evalAtPSPoint(const vector<double>& point, const vector<double>& masses) {
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

        const double mR0p = (*getSetting<double>("mR0p")) * unitres;
        const double mR1m = (*getSetting<double>("mR1m")) * unitres;
        // const double mR1p = (*getSetting<double>("mR1p"))*unitres;

        // Poles
        const double P0p = (1 - Sqq / (mR0p * mR0p));
        const double P1m = (1 - Sqq / (mR1m * mR1m));
        // const double P1p = (1 - Sqq/(mR1p*mR1p));

        // Parameters. Enforce f0(q^2=0) = f+(q^2=0)
        vector<double> Fzvec = {(*getSetting<vector<double>>("Fpvec"))[0]};
        for (auto& az : (*getSetting<vector<double>>("Fzvec"))) {
            Fzvec.push_back(az);
        }
        vector<pair<double, vector<double>>> PoleFFvec{};
        PoleFFvec.emplace_back(P0p, Fzvec);
        PoleFFvec.emplace_back(P1m, (*getSetting<vector<double>>("Fpvec")));
        PoleFFvec.emplace_back(P1m * Mb * (1 + rU), (*getSetting<vector<double>>("Ftvec")));

        // Compute and set elements
        // Fz, Fp, Ft
        for (size_t i = 0; i < PoleFFvec.size(); ++i) {
            double temp = 0;
            auto vec = PoleFFvec[i].second;
            for (size_t n = 0; n < vec.size(); ++n) {
                temp += vec[n] * zpow[n];
            }
            temp /= PoleFFvec[i].first;

            const auto idx = static_cast<IndexType>(i + 1u); // 0th element is Fs
            result.element({idx}) = temp;
        }
        // Equation of motion for Fs scalar FF
        result.element({0}) = Mb2 * (1 - rU * rU) / (mb + mu) * result.element({1});
    }

    std::unique_ptr<FormFactorBase> FFBtoPiGKvD::clone(const std::string& label) {
        MAKE_CLONE(FFBtoPiGKvD, label);
    }

    void FFBtoPiGKvD::addRefs() const {
        if (!getSettingsHandler()->checkReference("Gubernari:2018wyi")) {
            string ref =
                "@article{Gubernari:2018wyi,\n"
                "    author = \"Gubernari, Nico and Kokulu, Ahmet and van Dyk, Danny\",\n"
                "    title = \"{$B \\to P$ and $B \\to V$ Form Factors from $B$-Meson Light-Cone Sum Rules beyond "
                "Leading Twist}\",\n"
                "    eprint = \"1811.00983\",\n"
                "    archivePrefix = \"arXiv\",\n"
                "    primaryClass = \"hep-ph\",\n"
                "    reportNumber = \"EOS-2018-02, TUM-HEP-1172/18\",\n"
                "    doi = \"10.1007/JHEP01(2019)150\",\n"
                "    journal = \"JHEP\",\n"
                "    volume = \"01\",\n"
                "    pages = \"150\",\n"
                "    year = \"2019\"\n"
                "}\n";
            getSettingsHandler()->addReference("Gubernari:2018wyi", ref);
        }
    }

} // namespace Hammer
