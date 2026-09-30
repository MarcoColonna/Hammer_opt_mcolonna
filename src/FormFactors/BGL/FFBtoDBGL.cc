///
/// @file  FFBtoDBGL.cc
/// @brief \f$ B \rightarrow D \f$ BGL form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/FormFactors/BGL/FFBtoDBGL.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFBtoDBGL::FFBtoDBGL() {
        // Create tensor rank and dimensions
        IndexList dims = {4};
        string name{"FFBtoDBGL"};

        setPrefix("BtoD");
        addProcessSignature(PID::BPLUS, {-PID::D0});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BD})});

        addProcessSignature(PID::BZERO, {PID::DMINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BD})});

        setPrefix("BstoDs");
        addProcessSignature(PID::BS, {PID::DSMINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BSDS})});

        setSignatureIndex();
    }

    void FFBtoDBGL::defineSettings() {
        //_mFFErrNames = ;
        setPath(getFFErrPrefixGroup().get());
        setUnits("GeV");

        // Using 1606.08030
        addSetting<double>("ChiT", 5.131e-4); // GeV^-2
        addSetting<double>("ChiL", 6.332e-3);
        // approximate pole passes
        addSetting<double>("mb", 4.9);       // GeV
        addSetting<double>("mc", 4.9 - 3.4); // GeV (using mb-mc = 3.4 GeV)

        vector<double> apvec = {0.01565, -0.0353, -0.043, 0.194};
        vector<double> a0vec = {0.07932, -0.214, 0.17, -0.958};
        addSetting<vector<double>>("ap", apvec);
        addSetting<vector<double>>("a0", a0vec);

        vector<double> BcStatesp = {6.329, 6.920, 7.020}; // GeV
        vector<double> BcStates0 = {6.716, 7.121};        // GeV
        addSetting<vector<double>>("BcStatesp", BcStatesp);
        addSetting<vector<double>>("BcStates0", BcStates0);

        setInitialized();
    }

    void FFBtoDBGL::evalAtPSPoint(const vector<double>& point, const vector<double>& masses) {
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
        const double Mb2 = Mb * Mb;
        const double rC = Mc / Mb;
        const double rC2 = rC * rC;
        const double sqrC = sqrt(rC);

        double w = getW(Sqq, Mb, Mc);
        // safety measure if w==1.0
        if (isZero(w - 1.0)) {
            w += 1e-6;
        }

        const size_t nmax = 4;
        const double z = (sqrt(w + 1) - sqrt2) / (sqrt(w + 1) + sqrt2);
        vector<double> zpow{1.};
        for (size_t n = 1; n < nmax; ++n) {
            zpow.push_back(zpow[n - 1] * z);
        }

        const vector<double>& BcStatesp = (*getSetting<vector<double>>("BcStatesp"));
        const vector<double>& BcStates0 = (*getSetting<vector<double>>("BcStates0"));

        const vector<double>& ap = (*getSetting<vector<double>>("ap"));
        const vector<double>& a0 = (*getSetting<vector<double>>("a0"));

        const double nc = 2.6;
        const double chiT = (*getSetting<double>("ChiT")) / (unitres * unitres);
        const double chiL = (*getSetting<double>("ChiL"));
        const double tp = (Mb + Mc) * (Mb + Mc) / Mb2;
        const double tm = (Mb - Mc) * (Mb - Mc) / Mb2;
        const double sqtptm = sqrt(tp - tm);

        const double mb = (*getSetting<double>("mb")) * unitres;
        const double mc = (*getSetting<double>("mc")) * unitres;

        double Pp = 1.;
        for (double state : BcStatesp) {
            double sqtpBc = sqrt(tp - pow(state * unitres / Mb, 2));
            Pp *= ((z - ((sqtpBc - sqtptm) / (sqtpBc + sqtptm))) / (1 - z * ((sqtpBc - sqtptm) / (sqtpBc + sqtptm))));
        }

        double P0 = 1.;
        for (double state : BcStates0) {
            double sqtpBc = sqrt(tp - pow(state * unitres / Mb, 2));
            P0 *= ((z - ((sqtpBc - sqtptm) / (sqtpBc + sqtptm))) / (1 - z * ((sqtpBc - sqtptm) / (sqtpBc + sqtptm))));
        }

        double phip = 32. * sqrt(nc / (6. * pi * chiT * Mb2)) * rC2 * pow(1 + z, 2) * pow(1 - z, 0.5) /
                      pow((1 + rC) * (1 - z) + 2 * sqrC * (1 + z), 5);

        double phi0 = (1 - rC2) * 8. * sqrt(nc / (8. * pi * chiL)) * rC * (1 + z) * pow(1 - z, 1.5) /
                      pow((1 + rC) * (1 - z) + 2 * sqrC * (1 + z), 4);

        double Fp = 0;
        for (size_t n = 0; n < ap.size(); ++n) {
            Fp += ap[n] * zpow[n];
        }
        Fp /= (Pp * phip);

        double F0 = 0;
        for (size_t n = 0; n < a0.size(); ++n) {
            F0 += a0[n] * zpow[n];
        }
        F0 /= (P0 * phi0);

        // Scalar from eqn of motion
        const double Fs = F0 * Mb2 * (1. - rC * rC) / (mb - mc);

        // set elements, mapping to amplitude FF basis
        // Fs (dim +1)
        result.element({0}) = Fs;
        // Fz (dim 0)
        result.element({1}) = F0;
        // Fp (dim 0)
        result.element({2}) = Fp;
        // Ft (dim -1)
        // result.element({3}) = 0;
    }

    unique_ptr<FormFactorBase> FFBtoDBGL::clone(const string& label) {
        MAKE_CLONE(FFBtoDBGL, label);
    }

} // namespace Hammer
