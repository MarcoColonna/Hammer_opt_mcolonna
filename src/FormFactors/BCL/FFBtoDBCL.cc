///
/// @file  FFBtoDBCL.cc
/// @brief \f$ B \rightarrow D \f$ BCL form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/FormFactors/BCL/FFBtoDBCL.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFBtoDBCL::FFBtoDBCL() {
        // Create tensor rank and dimensions
        IndexList dims = {4};
        string name{"FFBtoDBCL"};

        setPrefix("BtoD");
        addProcessSignature(PID::BPLUS, {-PID::D0});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BD})});

        addProcessSignature(PID::BZERO, {PID::DMINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BD})});

        setSignatureIndex();
    }

    void FFBtoDBCL::defineSettings() {
        //_mFFErrNames = ;
        setPath(getFFErrPrefixGroup().get());
        setUnits("GeV");

        // if(getFFErrPrefixGroup().prefix == "BtoD"){
        // Using 2506.15256 Table IV
        vector<double> apvec = {0.8959, -8.03, 49.3};
        vector<double> a0vec = {0.7813, -3.38};
        addSetting<vector<double>>("ap", apvec);
        addSetting<vector<double>>("a0", a0vec);

        double m1m = 6.329; // GeV (see BGL class poles)
        double m0p = 6.716; // GeV
        addSetting<double>("m1m", m1m);
        addSetting<double>("m0p", m0p);
        addSetting<bool>("WithFpPole", false); // set pole factors to unity (per 2506.15256 and FLAG 2411.04268)
        addSetting<bool>("WithF0Pole", false);

        double sqtp = 0.; // GeV t+ = (mB + mD)^2 branch point, when set to zero uses the event masses.
        addSetting<double>("tp", sqtp * sqtp);

        addSetting<bool>("q2cons", true); // impose f_+(q^2=0) = f_0(q^2=0)
        //}

        setInitialized();
    }

    void FFBtoDBCL::evalAtPSPoint(const vector<double>& point, const vector<double>& masses) {
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
        // const double Mu2 = Mu*Mu;
        // const double rU = Mu/Mb;
        // const double rU2 = rU*rU;
        // const double sqrU = sqrt(rU);

        //        double w = getW(Sqq, Mb, Mu);
        //        const double wmax = (Mb*Mb + Mu*Mu)/(2.*Mb*Mu);
        //        //safety measure if w==1.0
        //        if(isZero(w - 1.0)) w += 1e-6;

        // Poles
        const double& m1m = (*getSetting<double>("m1m")) * unitres;
        const bool& withFpPole = (*getSetting<bool>("WithFpPole"));
        const double P1m = withFpPole ? (1. - Sqq / (m1m * m1m)) : 1.;
        const double& m0p = (*getSetting<double>("m0p")) * unitres;
        const bool& withF0Pole = (*getSetting<bool>("WithF0Pole"));
        const double P0p = withF0Pole ? (1. - Sqq / (m0p * m0p)) : 1.;

        // Parameters
        const vector<double>& ap = (*getSetting<vector<double>>("ap"));
        const vector<double>& a0 = (*getSetting<vector<double>>("a0"));

        const bool& q2cons = (*getSetting<bool>("q2cons"));

        // Optimized expansion, z and z(q^2=0)
        const size_t Nz = ap.size();
        const size_t N0 = a0.size();
        const size_t Nmax = Nz > N0 ? Nz : N0;

        double tp = (*getSetting<double>("tp")) * unitres * unitres; //(Mb + Mu)*(Mb + Mu);
        if (isZero(tp)) {
            tp = (Mb + Mu) * (Mb + Mu);
        }
        const double tm = (Mb - Mu) * (Mb - Mu);
        const double t0 = tp * (1. - sqrt(1. - tm / tp));

        const double z = (sqrt(tp - Sqq) - sqrt(tp - t0)) / (sqrt(tp - Sqq) + sqrt(tp - t0));
        const double z0 = (sqrt(tp) - sqrt(tp - t0)) / (sqrt(tp) + sqrt(tp - t0));
        //        const double z = (pow(rU, 0.25)*sqrt(w+1) - sqrt(rU + 1))/(pow(rU, 0.25)*sqrt(w+1) + sqrt(rU+1));
        //        const double z0 = (pow(rU, 0.25)*sqrt(wmax+1) - sqrt(rU + 1))/(pow(rU, 0.25)*sqrt(wmax+1) +
        //        sqrt(rU+1));
        vector<double> zpow{1.};
        vector<double> z0pow{1.};
        for (size_t n = 1; n < Nmax + 1; ++n) {
            zpow.push_back(zpow[n - 1] * z);
            z0pow.push_back(z0pow[n - 1] * z0);
        }

        // N = nz expansion.
        double Fp = 0;
        double F0 = 0;

        auto Nv = static_cast<double>(Nz);
        for (size_t n = 0; n < Nz; ++n) {
            auto nv = static_cast<double>(n);
            Fp += ap[n] * (zpow[n] - pow(-1, nv - Nv) * (nv / Nv) * zpow[Nz]);
        }
        Fp /= P1m;

        for (size_t n = 0; n < N0; ++n) {
            F0 += a0[n] * zpow[n];
        }

        if (q2cons) {
            double Fpq2 = 0;
            double F0q2 = 0;
            for (size_t n = 0; n < Nz; ++n) {
                auto nv = static_cast<double>(n);
                Fpq2 += ap[n] * (z0pow[n] - pow(-1, nv - Nv) * (nv / Nv) * z0pow[Nz]);
            }

            for (size_t n = 0; n < N0; ++n) {
                F0q2 += a0[n] * z0pow[n];
            }

            F0 += (Fpq2 - F0q2) / z0pow[N0] * zpow[N0];
        }
        F0 /= P0p;

        // set elements
        // Fs (dim +1)
        // result.element({0}) = 0;
        // Fz (dim 0)
        result.element({1}) = F0;
        // Fp (dim 0)
        result.element({2}) = Fp;
        // Ft (dim -1)
        // result.element({3}) = 0;
    }

    unique_ptr<FormFactorBase> FFBtoDBCL::clone(const string& label) {
        MAKE_CLONE(FFBtoDBCL, label);
    }

    void FFBtoDBCL::addRefs() const {
        if (!getSettingsHandler()->checkReference("Belle-II:2025rna")) {
            string ref =
                "@article{Belle-II:2025rna,\n"
                "    author = \"Adachi, I. and others\",\n"
                "    collaboration = \"Belle-II\",\n"
                "    title = \"{Determination of |Vcb| using "
                "B{\textrightarrow}D{\\ensuremath{\\ell}}{\\ensuremath{\\nu}}l decays at Belle II}\",\n"
                "    eprint = \"2506.15256\",\n"
                "    archivePrefix = \"arXiv\",\n"
                "    primaryClass = \"hep-ex\",\n"
                "    reportNumber = \"Belle II Preprint 2025-004, KEK Preprint 2025-1\",\n"
                "    doi = \"10.1103/vs8k-259v\",\n"
                "    journal = \"Phys. Rev. D\",\n"
                "    volume = \"112\",\n"
                "    number = \"11\",\n"
                "    pages = \"112009\",\n"
                "    year = \"2025\"\n"
                "}\n";
            getSettingsHandler()->addReference("Belle-II:2025rna", ref);
        }
        if (!getSettingsHandler()->checkReference("FlavourLatticeAveragingGroupFLAG:2024oxs")) {
            string ref =
                "@article{FlavourLatticeAveragingGroupFLAG:2024oxs,\n"
                "    author = \"Aoki, Y. and others\",\n"
                "    collaboration = \"Flavour Lattice Averaging Group (FLAG)\",\n"
                "    title = \"{FLAG review 2024}\",\n"
                "    eprint = \"2411.04268\",\n"
                "    archivePrefix = \"arXiv\",\n"
                "    primaryClass = \"hep-lat\",\n"
                "    reportNumber = \"CERN-TH-2024-192, FERMILAB-PUB-24-0785-T\",\n"
                "    doi = \"10.1103/nfzp-p5dn\",\n"
                "    journal = \"Phys. Rev. D\",\n"
                "    volume = \"113\",\n"
                "    number = \"1\",\n"
                "    pages = \"014508\",\n"
                "    year = \"2026\"\n"
                "}\n";
            getSettingsHandler()->addReference("FlavourLatticeAveragingGroupFLAG:2024oxs", ref);
        }
    }

} // namespace Hammer
