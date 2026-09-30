///
/// @file  FFBstoKBCLVar.cc
/// @brief \f$ B_s \rightarrow K \f$ BCL form factors with variations
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/FormFactors/BCL/FFBstoKBCLVar.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFBstoKBCLVar::FFBstoKBCLVar() {
        // Create tensor rank and dimensions
        IndexList dims = {4, 8}; // size is _mFFErrNames + 1 to include central values in zeroth component
        setGroup("BCLVar");              // override BCL base class FF group setting
        string name{"FFBstoKBCLVar"};

        setPrefix("BstoK");
        _mFFErrLabel = FF_BSK_VAR;

        addProcessSignature(PID::BS, {PID::KMINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BSK, FF_BSK_VAR})});

        setSignatureIndex();
    }

    void FFBstoKBCLVar::defineSettings() {
        setPath(getFFErrPrefixGroup().get());
        setUnits("GeV");

        // First 8 principal component directions
        defineAndAddErrSettings({"delta_e1", "delta_e2", "delta_e3", "delta_e4", "delta_e5", "delta_e6", "delta_e7"});

        // Fit central values using 2111.09849 Table 48

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

        // Row basis is FP a0123, FZ a012 2111.09849 Table 48
        vector<vector<double>> eigmat{
            {0.00219874, 0.00503931, 0.00545253, -0.0027279, -0.00802577, 0.00347875, -0.000291095},
            {-0.0200244, 0.0486776, 0.0211644, 0.0289072, 0.00640553, 0.000962853, -0.0000563584},
            {-0.281887, 0.0934452, -0.088924, 0.000503329, 0.000439422, 0.000193383, -0.0000115601},
            {-0.508664, -0.104292, 0.0279807, 0.0012435, 0.0000479501, 0.0000248239, -2.64461e-6},
            {-0.00237505, 0.00227081, 0.000895574, 0.00158555, -0.0039699, 0.00403484, 0.000294166},
            {-0.0368324, 0.0371708, 0.0108756, 0.0109032, -0.016378, -0.00230361, 0.0000509403},
            {-0.15814, 0.154108, 0.0633564, -0.0111585, 0.00201402, -0.0000221714, 0.000015919}};

        addSetting<vector<vector<double>>>("eigmatrix", eigmat);

        addSetting<bool>("q2cons", true); // impose f_+(q^2=0) = f_0(q^2=0)

        setInitialized();
    }

    void FFBstoKBCLVar::evalAtPSPoint(const vector<double>& point, const vector<double>& masses) {
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

        // Principal components
        const vector<vector<double>>& eigmat = (*getSetting<vector<vector<double>>>("eigmatrix"));

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

        // Fz (dim 0)
        result.element({1, 0}) = F0;
        // Fp (dim 0)
        result.element({2, 0}) = Fp;

        // Now compute remaining FF tensor entries
        // Logic: FF^x = a^x_n z^n, a^x_n = M_{nk} delta_k, delta_k = {1, delta_...}
        for (size_t k = 0; k < a0.size() + ap.size(); ++k) {
            double Fpentry = 0;
            double F0entry = 0;
            for (size_t n = 0; n < Nz; ++n) {
                auto nv = static_cast<double>(n);
                Fpentry += eigmat[n][k] * (zpow[n] - pow(-1, nv - Nv) * (nv / Nv) * zpow[Nz]);
            }
            Fpentry /= P1m;
            for (size_t n = 0; n < N0; ++n) {
                F0entry += eigmat[Nz + n][k] * zpow[n];
            }

            if (q2cons) {
                double Fpq2entry = 0;
                double F0q2entry = 0;
                for (size_t n = 0; n < Nz; ++n) {
                    auto nv = static_cast<double>(n);
                    Fpq2entry += eigmat[n][k] * (z0pow[n] - pow(-1, nv - Nv) * (nv / Nv) * z0pow[Nz]);
                }
                for (size_t n = 0; n < N0; ++n) {
                    F0q2entry += eigmat[Nz + n][k] * z0pow[n];
                }
                F0entry += (Fpq2entry - F0q2entry) / z0pow[N0] * zpow[N0];
            }
            F0entry /= P0p;

            const auto idxk = static_cast<IndexType>(k + 1u);
            // Fz (dim 0)
            result.element({1, idxk}) = F0entry;
            // Fp (dim 0)
            result.element({2, idxk}) = Fpentry;
        }
    }

    std::unique_ptr<FormFactorBase> FFBstoKBCLVar::clone(const std::string& label) {
        MAKE_CLONE(FFBstoKBCLVar, label);
    }

    void FFBstoKBCLVar::addRefs() const {
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
