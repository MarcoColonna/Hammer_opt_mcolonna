///
/// @file  FFBtoPiBCLVar.cc
/// @brief \f$ B \rightarrow \pi \f$ BCL form factors with variations
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/FormFactors/BCL/FFBtoPiBCLVar.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFBtoPiBCLVar::FFBtoPiBCLVar() {
        // Create tensor rank and dimensions
        IndexList dims = {4, 9}; // size is _mFFErrNames + 1 to include central values in zeroth component
        setGroup("BCLVar");              // override BCL base class FF group setting
        string name{"FFBtoPiBCLVar"};

        setPrefix("BtoPi");
        _mFFErrLabel = FF_BPI_VAR;
        addProcessSignature(PID::BPLUS, {PID::PI0});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BPI, FF_BPI_VAR})});

        addProcessSignature(PID::BZERO, {PID::PIMINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BPI, FF_BPI_VAR})});

        setSignatureIndex();
    }

    void FFBtoPiBCLVar::defineSettings() {
        // First 8 principal component directions
        setPath(getFFErrPrefixGroup().get());
        setUnits("GeV");

        defineAndAddErrSettings(
            {"delta_e1", "delta_e2", "delta_e3", "delta_e4", "delta_e5", "delta_e6", "delta_e7", "delta_e8"});

        // Fit central values 1503.07839 Table XIX
        vector<double> apvec = {0.419, -0.495, -0.43, 0.22};
        vector<double> a0vec = {0.510, -1.700, 1.53, 4.52};
        addSetting<vector<double>>("ap", apvec);
        addSetting<vector<double>>("a0", a0vec);

        double m1m = 5.325; // GeV
        addSetting<double>("m1m", m1m);

        // Row basis is FP a0123, FZ a0123 1503.07839 Table XIX
        vector<vector<double>> eigmat{
            {-0.000317026, -0.00474679, -0.00163407, -0.00129409, 0.00755911, 0.00916734, -0.000892999, 0.000211949},
            {0.000815168, -0.0477344, -0.0084102, -0.00571681, -0.0227319, 0.00276067, 0.00299688, 0.000064176},
            {0.0274469, 0.123528, 0.01989, -0.00766265, -0.020577, 0.000159396, -0.00311267, 9.60587e-6},
            {-0.0684033, 0.302285, 0.00383034, 0.000658676, 0.00511887, 0.000536951, 0.00175753, 8.90054e-6},
            {-0.00751714, -0.000351954, 0.00898363, 0.0111163, -0.00138701, 0.00990563, -0.0000664247, -0.000205036},
            {-0.0622473, -0.0109081, 0.00329588, -0.0518798, 0.00504893, 0.00154931, 0.000152337, -0.0000568998},
            {0.0898965, 0.0259681, -0.165351, -0.00139355, -0.00107535, 0.00039169, -0.000444994, -0.0000166918},
            {0.829619, 0.0172467, 0.0179113, -0.00332793, 0.00161083, 0.000203347, 0.000303651, -3.88438e-6}};

        addSetting<vector<vector<double>>>("eigmatrix", eigmat);

        addSetting<bool>("q2cons", false); // impose f_+(q^2=0) = f_0(q^2=0)

        setInitialized();
    }

    void FFBtoPiBCLVar::evalAtPSPoint(const vector<double>& point, const vector<double>& masses) {
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

        const double tp = (Mb + Mu) * (Mb + Mu);
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

            const auto idxk = static_cast<IndexType>(k + 1u);
            // Fz (dim 0)
            result.element({1, idxk}) = F0entry;
            // Fp (dim 0)
            result.element({2, idxk}) = Fpentry;
        }
    }

    std::unique_ptr<FormFactorBase> FFBtoPiBCLVar::clone(const std::string& label) {
        MAKE_CLONE(FFBtoPiBCLVar, label);
    }

    void FFBtoPiBCLVar::addRefs() const {
        if (!getSettingsHandler()->checkReference("Aoki:2019cca")) {
            string ref =
                "@article{Aoki:2019cca,\n"
                "   author = \"Aoki, S. and others\",\n"
                "   collaboration = \"Flavour Lattice Averaging Group\",\n"
                "   title = \"{FLAG Review 2019: Flavour Lattice Averaging Group (FLAG)}\",\n"
                "   eprint = \"1902.08191\",\n"
                "   archivePrefix = \"arXiv\",\n"
                "   primaryClass = \"hep-lat\",\n"
                "   reportNumber = \"FERMILAB-PUB-19-077-T\",\n"
                "   doi = \"10.1140/epjc/s10052-019-7354-7\",\n"
                "   journal = \"Eur. Phys. J. C\",\n"
                "   volume = \"80\",\n"
                "   number = \"2\",\n"
                "   pages = \"113\",\n"
                "   year = \"2020\"\n"
                "}\n";
            getSettingsHandler()->addReference("Aoki:2019cca", ref);
        }
        if (!getSettingsHandler()->checkReference("FermilabLattice:2015mwy")) {
            string ref =
                "@article{FermilabLattice:2015mwy,\n"
                "    author = \"Bailey, Jon A. and others\",\n"
                "    collaboration = \"Fermilab Lattice, MILC\",\n"
                "    title = \"{$|V_{ub}|$ from $B\\to\\pi\\ell\\nu$ decays and (2+1)-flavor lattice QCD}\",\n"
                "    eprint = \"1503.07839\",\n"
                "    archivePrefix = \"arXiv\",\n"
                "    primaryClass = \"hep-lat\",\n"
                "    reportNumber = \"FERMILAB-PUB-15-108-T\",\n"
                "    doi = \"10.1103/PhysRevD.92.014024\",\n"
                "    journal = \"Phys. Rev. D\",\n"
                "    volume = \"92\",\n"
                "    number = \"1\",\n"
                "    pages = \"014024\",\n"
                "    year = \"2015\"\n"
                "}\n";
            getSettingsHandler()->addReference("FermilabLattice:2015mwy", ref);
        }
    }

} // namespace Hammer
