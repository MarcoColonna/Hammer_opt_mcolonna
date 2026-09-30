///
/// @file  FFBLPRBase.cc
/// @brief Hammer base class for BLPR form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/FormFactors/BLPRXP/FFBLPRXPVarBase.hh"

using namespace std;

namespace Hammer {

    FFBLPRXPVarBase::FFBLPRXPVarBase() {
        setGroup("BLPRXPVar");
    }

    void FFBLPRXPVarBase::defineSettings() {

        defineAndAddErrSettings({"delta_RhoSq", "delta_cSt", "delta_chi21", "delta_chi2p", "delta_chi3p", "delta_eta1",
                                 "delta_etap", "delta_phi1p", "delta_beta21", "delta_beta3p"});

        FFBLPRXPBase::defineSettings(); // common parameters for BLPRXP are defined in FFBLPRXPBase

        // correlation matrix and set zero error eigenvectors
        // Row basis is RhoSq, cStar, chi2(1), chi2(1)prime, chi3(1)prime, eta(1), eta(1)prime, phi_1(1)prime,
        // beta_2(1), beta_3(1)prime
        vector<vector<double>> sliwmat{
            {1., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 1., 0., 0., 0., 0., 0., 0., 0., 0.},
            {0., 0., 1., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 0., 1., 0., 0., 0., 0., 0., 0.},
            {0., 0., 0., 0., 1., 0., 0., 0., 0., 0.}, {0., 0., 0., 0., 0., 1., 0., 0., 0., 0.},
            {0., 0., 0., 0., 0., 0., 1., 0., 0., 0.}, {0., 0., 0., 0., 0., 0., 0., 1., 0., 0.},
            {0., 0., 0., 0., 0., 0., 0., 0., 1., 0.}, {0., 0., 0., 0., 0., 0., 0., 0., 0., 1.}};
        addSetting<vector<vector<double>>>("sliwmatrix", sliwmat);
    }

    array<double, 10> FFBLPRXPVarBase::XiLinearCoeffs(double w) const {
        // Variables for leading IW function (derived from G(1))
        const double a2 = _a * _a;
        const double a4 = a2 * a2;
        const double cSt = *getSetting<double>("cSt");
        const double RhoSq = *getSetting<double>("RhoStSq");
        const double zCon = ZofW(w);
        const double zCon1 = (1 - _a) / (1 + _a);
        const double num = (1. - 8 * a2 * RhoSq * zCon + 16. * (2 * cSt * a4 - RhoSq * a2) * zCon * zCon);
        const double den = (1. - 8 * a2 * RhoSq * zCon1 + 16. * (2 * cSt * a4 - RhoSq * a2) * zCon1 * zCon1);
        // Linearization of LO IW functions
        array<double, 10> result = {-(8 * a2 * zCon * (1 + 2. * zCon) / num - 8 * a2 * zCon1 * (1 + 2. * zCon1) / den),
                                    (32 * a4 * zCon * zCon / num - 32 * a4 * zCon1 * zCon1 / den),
                                    0.,
                                    0.,
                                    0.,
                                    0.,
                                    0.,
                                    0.,
                                    0.,
                                    0.};

        return result;
    }


} // namespace Hammer
