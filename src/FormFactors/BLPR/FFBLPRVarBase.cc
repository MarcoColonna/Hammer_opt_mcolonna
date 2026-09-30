///
/// @file  FFBLPRBase.cc
/// @brief Hammer base class for BLPR form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/FormFactors/BLPR/FFBLPRVarBase.hh"
#include "Hammer/Math/Constants.hh"

using namespace std;

namespace Hammer {

    FFBLPRVarBase::FFBLPRVarBase() {
        setGroup("BLPRVar");
    }

    void FFBLPRVarBase::defineSettings() {

        defineAndAddErrSettings(
            {"delta_RhoSq", "delta_chi21", "delta_chi2p", "delta_chi3p", "delta_eta1", "delta_etap", "delta_dV20"});

        FFBLPRBase::defineSettings(); // common parameters for BLPR are defined in FFBLPRBase

        // correlation matrix and set zero error eigenvectors
        // Row basis is RhoSq, chi2(1), chi2(1)prime, chi3(1)prime, eta(1), eta(1)prime, V20
        vector<vector<double>> sliwmat{{1., 0., 0., 0., 0., 0., 0.}, {0., 1., 0., 0., 0., 0., 0.},
                                       {0., 0., 1., 0., 0., 0., 0.}, {0., 0., 0., 1., 0., 0., 0.},
                                       {0., 0., 0., 0., 1., 0., 0.}, {0., 0., 0., 0., 0., 1., 0.},
                                       {0., 0., 0., 0., 0., 0., 1.}};
        addSetting<vector<vector<double>>>("sliwmatrix", sliwmat);
    }

    vector<double> FFBLPRVarBase::XiLinearCoeffs(double w) const {
        // Variables for leading IW function (derived from G(1))
        const double V21 = 57.0;
        // const double V20 = 7.5 + (*getSetting<double>("dV20"));
        const double RhoSq = *getSetting<double>("RhoSq");
        const double rD = 1867. / 5280.;
        const double eb = (*getSetting<double>("la")) / (*getSetting<double>("mb") * 2.);
        const double ec = (*getSetting<double>("la")) / (*getSetting<double>("mc") * 2.);
        const double as = (*getSetting<double>("as")) / pi;
        const double chi21 = (*getSetting<double>("chi21"));
        const double chi2p = (*getSetting<double>("chi2p"));
        const double chi3p = (*getSetting<double>("chi3p"));
        const double etap = (*getSetting<double>("etap"));
        const double a2 = _a * _a;
        const double a4 = a2 * a2;
        const double a6 = a4 * a2;
        const double zCon = ZofW(w);
        const double zConSq = zCon * zCon;
        auto [xiIW, xiIWN] = Xi(w);
        const double Xi = 64 * a4 * RhoSq - 16 * a2 - V21;
        // const double wm1Sq = (w-1)*(w-1);
        const double am1Sq = (_a - 1) * (_a - 1);
        const double ap1Sq = (_a + 1) * (_a + 1);

        // Linearization of LO IW functions
        const vector<double> xiIWL = {
            -8 * a2 * zCon +
                zConSq * (V21 +
                          (64 * a4 * _mCv1zp + (64 * a4 * _mCv3z) / (1 + rD) + (128 * a6 * _mCv3zp) / (1 + rD) +
                           (64 * a4 * _mCv2z * rD) / (1 + rD) + (128 * a6 * _mCv2zp * rD) / (1 + rD)) *
                              as -
                          256 * a4 * (eb + ec) * chi21 - 16 * (-64 * a4 + 64 * a6) * (eb + ec) * chi2p +
                          768 * a4 * (eb + ec) * chi3p - (128 * a4 * (-1 + rD) * (eb - ec) * etap) / (1 + rD)),
            -4 * zConSq * (eb + ec) * Xi,
            -16 * zConSq * (eb + ec) * (-(a2 * (-16 + V21)) + V21 + 64 * a6 * RhoSq - 32 * a4 * (1 + 2 * RhoSq)),
            12 * zConSq * (eb + ec) * Xi,
            0,
            (-2 * (-1 + rD) * zConSq * (eb - ec) * Xi) / (1 + rD),
            -zConSq};

        const vector<double> xiIWLN = {
            (8 * a4 + a2 * (-8 + V21) + V21 - 2 * _a * V21) / ap1Sq +
                ((64 * am1Sq * a4 * _mCv1zp) / ap1Sq + (64 * am1Sq * a4 * _mCv3z) / (ap1Sq * (1 + rD)) +
                 (128 * am1Sq * a6 * _mCv3zp) / (ap1Sq * (1 + rD)) +
                 (64 * am1Sq * a4 * _mCv2z * rD) / (ap1Sq * (1 + rD)) +
                 (128 * am1Sq * a6 * _mCv2zp * rD) / (ap1Sq * (1 + rD))) *
                    as -
                (256 * am1Sq * a4 * (eb + ec) * chi21) / ap1Sq -
                (16 * am1Sq * (-64 * a4 + 64 * a6) * (eb + ec) * chi2p) / ap1Sq +
                (768 * am1Sq * a4 * (eb + ec) * chi3p) / ap1Sq -
                (128 * am1Sq * a4 * (-1 + rD) * (eb - ec) * etap) / (ap1Sq * (1 + rD)),
            (-4 * am1Sq * (eb + ec) * Xi) / ap1Sq,
            (-16 * am1Sq * (eb + ec) * (-(a2 * (-16 + V21)) + V21 + 64 * a6 * RhoSq - 32 * a4 * (1 + 2 * RhoSq))) /
                ap1Sq,
            (12 * am1Sq * (eb + ec) * Xi) / ap1Sq,
            0,
            (-2 * am1Sq * (-1 + rD) * (eb - ec) * Xi) / (ap1Sq * (1 + rD)),
            (-1 + 2 * _a - a2) / ap1Sq};

        vector<double> result(7);
        for (size_t i = 0; i < 7; ++i) {
            result[i] = xiIWL[i] / xiIW - xiIWLN[i] / xiIWN;
        }
        return result;
    }


} // namespace Hammer
