///
/// @file  FFBtoD2starBLRVar.cc
/// @brief \f$ B \rightarrow D_2^* \f$ BLR form factors with variations
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/FormFactors/BLR/FFBtoD2starBLRVar.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFBtoD2starBLRVar::FFBtoD2starBLRVar() {
        // Create tensor rank and dimensions
        IndexList dims = {8, 8}; // size is _mFFErrNames + 1 to include central values in zeroth component
        setGroup("BLRVar");              // override BLPR base class FF group setting
        string name{"FFBtoD2starBLRVar"};

        setPrefix("BtoD**2*");
        _mFFErrLabel = FF_BDSSD2STAR_VAR;
        addProcessSignature(PID::BPLUS, {-PID::DSSD2STAR});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BDSSD2STAR, FF_BDSSD2STAR_VAR})});

        addProcessSignature(PID::BZERO, {PID::DSSD2STARMINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BDSSD2STAR, FF_BDSSD2STAR_VAR})});

        setPrefix("BstoDs**2*");
        _mFFErrLabel = FF_BSDSSDS2STAR_VAR;
        addProcessSignature(PID::BS, {PID::DSSDS2STARMINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BSDSSDS2STAR, FF_BSDSSDS2STAR_VAR})});

        setSignatureIndex();
    }

    void FFBtoD2starBLRVar::defineSettings() {
        setPath(getFFErrPrefixGroup().get());
        setUnits("GeV");

        defineAndAddErrSettings(
            {"delta_t1", "delta_tp", "delta_tau1", "delta_tau2", "delta_eta1", "delta_eta2", "delta_eta3"});

        // 1S scheme: mb1S = 4710, mbBar = 5313, lambda1 = -3 10^5 MeV^2, delta mb-mc = 3400, alpha_s = 26/100
        addSetting<double>("as", 0.26);
        addSetting<double>("mb", 4.710);
        addSetting<double>("mc", 4.710 - 3.400);
        addSetting<double>("t1", 0.7);
        addSetting<double>("tp", -1.6);
        addSetting<double>("tau1", -0.5);
        addSetting<double>("tau2", 2.9);
        addSetting<double>("eta1", 0.);
        addSetting<double>("eta2", 0.);
        addSetting<double>("eta3", 0.);
        addSetting<double>("laB", 0.4);
        addSetting<double>("laP", 0.8);

        // correlation matrix and set zero error eigenvectors
        // Row basis is t1, tp, tau1, tau2, eta1(1), eta2(1), eta3(1)
        vector<vector<double>> sliwmat{{1., 0., 0., 0., 0., 0., 0.}, {0., 1., 0., 0., 0., 0., 0.},
                                       {0., 0., 1., 0., 0., 0., 0.}, {0., 0., 0., 1., 0., 0., 0.},
                                       {0., 0., 0., 0., 1., 0., 0.}, {0., 0., 0., 0., 0., 1., 0.},
                                       {0., 0., 0., 0., 0., 0., 1.}};
        addSetting<vector<vector<double>>>("sliwmatrix", sliwmat);

        setInitialized();
    }

    void FFBtoD2starBLRVar::evalAtPSPoint(const vector<double>& point, const vector<double>& masses) {
        Tensor& result = getTensor();
        result.clearData();

        if (!isInitialized()) {
            MSG_WARNING("Warning, Settings have not been defined!");
        }

        double Mb = 0.;
        double Mc = 0.;
        double unitres = 1.;
        tie(Mb, Mc, unitres) = getParentDaughterHadMasses(masses);


        double Sqq = point[0];

        // const double sqSqq = sqrt(Sqq);
        double w = getW(Sqq, Mb, Mc);
        // safety measure if w==1.0
        if (isZero(w - 1.0)) {
            w += 1e-6;
        }

        // BLR parameters
        const double zBC = (*getSetting<double>("mc")) / (*getSetting<double>("mb"));
        const double eB = 1. / (*getSetting<double>("mb") * 2.);
        const double eC = 1. / (*getSetting<double>("mc") * 2.);
        const double as = (*getSetting<double>("as")) / pi;
        const double t1 = (*getSetting<double>("t1"));
        const double tp = (*getSetting<double>("tp"));
        const double tau1 = (*getSetting<double>("tau1"));
        const double tau2 = (*getSetting<double>("tau2"));
        const double eta1 = (*getSetting<double>("eta1"));
        const double eta2 = (*getSetting<double>("eta2"));
        const double eta3 = (*getSetting<double>("eta3"));
        const double laB = (*getSetting<double>("laB"));
        const double laP = (*getSetting<double>("laP"));

        const vector<vector<double>>& sliwmat = (*getSetting<vector<vector<double>>>("sliwmatrix"));

        // const double LambdaD32 = -laB + laP*w;
        const double Fb = laB + laP - tau2 - tau1 * (1 + 2 * w);
        const double LOIWtau = t1 + (w - 1) * t1 * tp;

        // QCD correction functions
        //  const double Cs  = CS(w, zBC);
        const double Cps = _asCorrections.CP(w, zBC);
        const double Cv1 = _asCorrections.CV1(w, zBC);
        // const double Cv2 = CV2(w, zBC);
        // const double Cv3 = CV3(w, zBC);
        const double Ca1 = _asCorrections.CA1(w, zBC);
        const double Ca2 = _asCorrections.CA2(w, zBC);
        const double Ca3 = _asCorrections.CA3(w, zBC);
        const double Ct1 = _asCorrections.CT1(w, zBC);
        const double Ct2 = _asCorrections.CT2(w, zBC);
        const double Ct3 = _asCorrections.CT3(w, zBC);

        // Form factors
        const double Kp =
            1 + eC * (-2 * eta1 - 2 * (w - 1) * eta2 + eta3 + (1 + 2 * w) * tau1 + tau2) + as * Cps + eB * Fb;
        const double Ka1 = -(eC * (-((1 + w) * (2 * eta1 - eta3)) + (w - 1) * (tau1 - tau2))) +
                           (-1 - w) * (1 + as * Ca1) - (w - 1) * eB * Fb;
        const double Ka2 = -2 * eC * (eta2 + tau1) + as * Ca2;
        const double Ka3 = 1 - eC * (2 * eta1 - 2 * eta2 - eta3 + tau1 + tau2) + as * (Ca1 + Ca3) + eB * Fb;
        const double Kv = -1 - eC * (-2 * eta1 + eta3 + tau1 - tau2) - as * Cv1 - eB * Fb;
        const double Kt1 = 1 - eC * (2 * eta1 - eta3) + as * (Ct1 + ((w - 1) * (Ct2 - Ct3)) / 2.);
        const double Kt2 = -(eC * (tau1 - tau2)) + ((1 + w) * as * (Ct2 + Ct3)) / 2. + eB * Fb;
        const double Kt3 = 2 * eC * (-eta2 + tau1) - as * Ct2;

        // defined as 1/LOIWtau * D[LOIWtau * f, var]
        const array<double, 7> KpDer = {
            Kp / t1, Kp / LOIWtau * t1 * (w - 1), -(eB - eC) * (1 + 2 * w), -eB + eC, -2 * eC, -2 * eC * (-1 + w), eC};
        const array<double, 7> Ka1Der = {Ka1 / t1,
                                         Ka1 / LOIWtau * t1 * (w - 1),
                                         (-1 + w) * (eB - eC + 2 * eB * w),
                                         (eB + eC) * (-1 + w),
                                         2 * eC * (1 + w),
                                         0,
                                         -eC * (1 + w)};
        const array<double, 7> Ka2Der = {Ka2 / t1, Ka2 / LOIWtau * t1 * (w - 1), -2 * eC, 0., 0., -2 * eC, 0.};
        const array<double, 7> Ka3Der = {
            Ka3 / t1, Ka3 / LOIWtau * t1 * (w - 1), -eC - eB * (1 + 2 * w), -eB - eC, -2 * eC, 2 * eC, eC};
        const array<double, 7> KvDer = {Kv / t1, Kv / LOIWtau * t1 * (w - 1), eB - eC + 2 * eB * w, eB + eC, 2 * eC, 0.,
                                        -eC};
        const array<double, 7> Kt1Der = {Kt1 / t1, Kt1 / LOIWtau * t1 * (w - 1), 0., 0., -2 * eC, 0., eC};
        const array<double, 7> Kt2Der = {
            Kt2 / t1, Kt2 / LOIWtau * t1 * (w - 1), -eC - eB * (1 + 2 * w), -eB + eC, 0., 0., 0.};
        const array<double, 7> Kt3Der = {Kt3 / t1, Kt3 / LOIWtau * t1 * (w - 1), 2 * eC, 0., 0., -2 * eC, 0.};

        // Set elements
        result.element({0, 0}) = Kp;
        result.element({1, 0}) = Ka1;
        result.element({2, 0}) = Ka2;
        result.element({3, 0}) = Ka3;
        result.element({4, 0}) = Kv;
        result.element({5, 0}) = Kt1;
        result.element({6, 0}) = Kt2;
        result.element({7, 0}) = Kt3;

        for (IndexType i1 = 1; i1 <= 7; ++i1) {
            for (IndexType i2 = 0; i2 < 7; ++i2) {
                result.element({0, i1}) += sliwmat[i2][i1 - 1] * KpDer[i2];
                result.element({1, i1}) += sliwmat[i2][i1 - 1] * Ka1Der[i2];
                result.element({2, i1}) += sliwmat[i2][i1 - 1] * Ka2Der[i2];
                result.element({3, i1}) += sliwmat[i2][i1 - 1] * Ka3Der[i2];
                result.element({4, i1}) += sliwmat[i2][i1 - 1] * KvDer[i2];
                result.element({5, i1}) += sliwmat[i2][i1 - 1] * Kt1Der[i2];
                result.element({6, i1}) += sliwmat[i2][i1 - 1] * Kt2Der[i2];
                result.element({7, i1}) += sliwmat[i2][i1 - 1] * Kt3Der[i2];
            }
        }

        result *= LOIWtau;
        result.toVector();
    }

    std::unique_ptr<FormFactorBase> FFBtoD2starBLRVar::clone(const std::string& label) {
        MAKE_CLONE(FFBtoD2starBLRVar, label);
    }

} // namespace Hammer
