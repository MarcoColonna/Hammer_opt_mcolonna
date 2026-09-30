///
/// @file  FFBstoKBCL.cc
/// @brief \f$ B_s \rightarrow K \f$ BCL form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/FormFactors/BCL/FFBstoKBCL.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFBstoKBCL::FFBstoKBCL() {
        // Create tensor rank and dimensions
        IndexList dims = {4};
        string name{"FFBstoKBCL"};

        setPrefix("BstoK");
        addProcessSignature(PID::BS, {PID::KMINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BSK})});

        setSignatureIndex();
    }

    void FFBstoKBCL::defineSettings() {
        //_mFFErrNames = ;
        setPath(getFFErrPrefixGroup().get());
        setUnits("GeV");

        // Using 2111.09849 Table 48

        vector<double> apvec = {0.374, -0.672, 0.07, 1.34};
        vector<double> a0vec = {0.2203, 0.089, 0.24};
        addSetting<vector<double>>("ap", apvec);
        addSetting<vector<double>>("a0", a0vec);

        double m1m = 5.325; // GeV
        double m0p = 5.68;  // GeV
        addSetting<double>("m1m", m1m);
        addSetting<double>("m0p", m0p);

        double sqtp = 5.280 + 0.13957; // GeV t+ = (mB + mPi)^2 branch point
        addSetting<double>("tp", sqtp * sqtp);

        addSetting<bool>("q2cons", true); // impose f_+(q^2=0) = f_0(q^2=0)

        setInitialized();
    }

    void FFBstoKBCL::evalAtPSPoint(const vector<double>& point, const vector<double>& masses) {
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

        // Pole
        const double& m1m = (*getSetting<double>("m1m")) * unitres;
        const double P1m = (1 - Sqq / (m1m * m1m));
        const double& m0p = (*getSetting<double>("m0p")) * unitres;
        const double P0p = (1 - Sqq / (m0p * m0p));

        // Parameters
        const vector<double>& ap = (*getSetting<vector<double>>("ap"));
        const vector<double>& a0 = (*getSetting<vector<double>>("a0"));

        const bool& q2cons = (*getSetting<bool>("q2cons"));

        // Optimized expansion, z and z(q^2=0)
        const size_t Nz = ap.size();
        const size_t N0 = a0.size();
        const size_t Nmax = Nz > N0 ? Nz : N0;

        const double tp = (*getSetting<double>("tp")) * unitres * unitres; //(Mb + Mu)*(Mb + Mu);
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

    unique_ptr<FormFactorBase> FFBstoKBCL::clone(const string& label) {
        MAKE_CLONE(FFBstoKBCL, label);
    }

    void FFBstoKBCL::addRefs() const {
        if (!getSettingsHandler()->checkReference("FlavourLatticeAveragingGroupFLAG:2021npn")) {
            string ref =
                "@article{FlavourLatticeAveragingGroupFLAG:2021npn,\n"
                "    author = \"Aoki, Y. and others\",\n"
                "    collaboration = \"Flavour Lattice Averaging Group (FLAG)\",\n"
                "    title = \"{FLAG Review 2021}\",\n"
                "    eprint = \"2111.09849\",\n"
                "    archivePrefix = \"arXiv\",\n"
                "    primaryClass = \"hep-lat\",\n"
                "    reportNumber = \"CERN-TH-2021-191, JLAB-THY-21-3528, FERMILAB-PUB-21-620-SCD-T\",\n"
                "    doi = \"10.1140/epjc/s10052-022-10536-1\",\n"
                "    journal = \"Eur. Phys. J. C\",\n"
                "    volume = \"82\",\n"
                "    number = \"10\",\n"
                "    pages = \"869\",\n"
                "    year = \"2022\"\n"
                "}\n";
            getSettingsHandler()->addReference("FlavourLatticeAveragingGroupFLAG:2021npn", ref);
        }
    }

} // namespace Hammer
