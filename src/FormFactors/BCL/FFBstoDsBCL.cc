///
/// @file  FFBstoDsBCL.cc
/// @brief \f$ B_s \rightarrow D_s \f$ BCL form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/FormFactors/BCL/FFBstoDsBCL.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFBstoDsBCL::FFBstoDsBCL() {
        // Create tensor rank and dimensions
        IndexList dims = {4};
        string name{"FFBstoDsBCL"};

        setPrefix("BstoDs");
        addProcessSignature(PID::BS, {PID::DSMINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BSDS})});

        setSignatureIndex();
    }

    void FFBstoDsBCL::defineSettings() {
        //_mFFErrNames = ;
        setPath(getFFErrPrefixGroup().get());
        setUnits("GeV");

        // Using 1906.00701 App A

        vector<double> apvec = {0.66574, -3.23599, 0.07};
        vector<double> a0vec = {0.66574, -0.25944, -0.10636};
        addSetting<vector<double>>("ap", apvec);
        addSetting<vector<double>>("a0", a0vec);

        double mBc0 = 6.704; // GeV
        double mBcs = 6.329; // GeV
        addSetting<double>("mBc0", mBc0);
        addSetting<double>("mBcs", mBcs);

        double sqtp = 5.3669 + 1.9683; // GeV t+ = (mBs + mDs)^2 branch point
        addSetting<double>("tp", sqtp * sqtp);

        addSetting<bool>("q2cons", true); // impose f_+(q^2=0) = f_0(q^2=0)

        setInitialized();
    }

    void FFBstoDsBCL::evalAtPSPoint(const vector<double>& point, const vector<double>& masses) {
        Tensor& result = getTensor();
        result.clearData();

        if (!isInitialized()) {
            MSG_WARNING("Warning, Settings have not been defined!");
        }

        double Mb = 0.;
        double Mc = 0.;
        double unitres = 1.;
        tie(Mb, Mc, unitres) = getParentDaughterHadMasses(masses);

        const double Sqq = point[0];

        // Pole
        const double& mBcs = (*getSetting<double>("mBcs")) * unitres;
        const double P1m = (1 - Sqq / (mBcs * mBcs));
        const double& mBc0 = (*getSetting<double>("mBc0")) * unitres;
        const double P0p = (1 - Sqq / (mBc0 * mBc0));

        // Parameters
        const vector<double>& ap = (*getSetting<vector<double>>("ap"));
        const vector<double>& a0 = (*getSetting<vector<double>>("a0"));

        const bool& q2cons = (*getSetting<bool>("q2cons"));

        // Optimized expansion, z and z(q^2=0)
        const size_t Nz = ap.size();
        const size_t N0 = a0.size();
        const size_t Nmax = Nz > N0 ? Nz : N0;

        const double& tp = (*getSetting<double>("tp")) * unitres * unitres; //(Mb + Mc)*(Mb + Mc);
        //        const double tm = (Mb - Mc)*(Mb - Mc);
        const double t0 = 0.; // sets z(q^2=0) = 0

        const double z = (sqrt(tp - Sqq) - sqrt(tp - t0)) / (sqrt(tp - Sqq) + sqrt(tp - t0));
        const double z0 = (sqrt(tp) - sqrt(tp - t0)) / (sqrt(tp) + sqrt(tp - t0)); // z at q2=0

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

            F0 += (Fpq2 - F0q2) / z0pow[0] *
                  zpow[0]; // add additional constant term to a0[0] (different to usual BCL, that constrains a0[N0])
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

    unique_ptr<FormFactorBase> FFBstoDsBCL::clone(const string& label) {
        MAKE_CLONE(FFBstoDsBCL, label);
    }

    void FFBstoDsBCL::addRefs() const {
        if (!getSettingsHandler()->checkReference("McLean:2019qcx")) {
            string ref =
                "@article{McLean:2019qcx,\n"
                "    author = \"McLean, E. and Davies, C. T. H. and Koponen, J. and Lytle, A. T.\",\n"
                "    title = \"{$B_s\\to D_s \\ell\\nu$ Form Factors for the full $q^2$ range from Lattice QCD with "
                "non-perturbatively normalized currents}\",\n"
                "    eprint = \"1906.00701\",\n"
                "    archivePrefix = \"arXiv\",\n"
                "    primaryClass = \"hep-lat\",\n"
                "    doi = \"10.1103/PhysRevD.101.074513\",\n"
                "    journal = \"Phys. Rev. D\",\n"
                "    volume = \"101\",\n"
                "    number = \"7\",\n"
                "    pages = \"074513\",\n"
                "    year = \"2020\"\n"
                "}\n";
            getSettingsHandler()->addReference("McLean:2019qcx", ref);
        }
    }

} // namespace Hammer
