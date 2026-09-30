///
/// @file  FFBtoD1BLR.cc
/// @brief \f$ B \rightarrow D_1 \f$ BLR form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/FormFactors/BLR/FFBtoD1BLR.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFBtoD1BLR::FFBtoD1BLR() {
        // Create tensor rank and dimensions
        IndexList dims = {8};
        string name{"FFBtoD1BLR"};

        setPrefix("BtoD**1");
        addProcessSignature(PID::BPLUS, {-PID::DSSD1});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BDSSD1})});

        addProcessSignature(PID::BZERO, {PID::DSSD1MINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BDSSD1})});

        setPrefix("BstoDs**1");
        addProcessSignature(PID::BS, {PID::DSSDS1MINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BSDSSDS1})});

        setSignatureIndex();
    }

    void FFBtoD1BLR::defineSettings() {
        //_mFFErrNames = ;
        setPath(getFFErrPrefixGroup().get());
        setUnits("GeV");

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

        setInitialized();
    }

    void FFBtoD1BLR::evalAtPSPoint(const vector<double>& point, const vector<double>& masses) {
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

        const double LambdaD32 = -laB + laP * w;
        const double Fb = laB + laP - tau2 - tau1 * (1 + 2 * w);
        const double LOIWtau = t1 + (w - 1) * t1 * tp;

        // QCD correction functions
        const double Cs = _asCorrections.CS(w, zBC);
        // const double Cps = CP(w, zBC);
        const double Cv1 = _asCorrections.CV1(w, zBC);
        const double Cv2 = _asCorrections.CV2(w, zBC);
        const double Cv3 = _asCorrections.CV3(w, zBC);
        const double Ca1 = _asCorrections.CA1(w, zBC);
        // const double Ca2 = CA2(w, zBC);
        // const double Ca3 = CA3(w, zBC);
        const double Ct1 = _asCorrections.CT1(w, zBC);
        const double Ct2 = _asCorrections.CT2(w, zBC);
        const double Ct3 = _asCorrections.CT3(w, zBC);

        // Form factors
        const double Fsc = (-(eC * (4 * LambdaD32 + 2 * (1 + w) * (6 * eta1 + 2 * (w - 1) * eta2 - eta3) -
                                    2 * (w - 1) * ((1 + 2 * w) * tau1 + tau2))) -
                            2 * (1 + w) * (1 + as * Cs) - 2 * (w - 1) * eB * Fb) /
                           sqrt(6.);
        const double Fv1 =
            (-(eC * (4 * (1 + w) * LambdaD32 - (-1 + w * w) * (2 * eta1 + 3 * eta3 + 3 * tau1 - 3 * tau2))) +
             (1 - w * w) * (1 + as * Cv1) - (-1 + w * w) * eB * Fb) /
            sqrt(6.);
        const double Fv2 = (-3 - eC * (10 * eta1 + 4 * (w - 1) * eta2 - 5 * eta3 + (-1 + 4 * w) * tau1 + 5 * tau2) -
                            as * (3 * Cv1 + 2 * (1 + w) * Cv2) - 3 * eB * Fb) /
                           sqrt(6.);
        const double Fv3 = (-2 + w +
                            eC * (4 * LambdaD32 - 2 * (6 + w) * eta1 - 4 * (w - 1) * eta2 - (-2 + 3 * w) * eta3 +
                                  (2 + w) * tau1 + (2 + 3 * w) * tau2) -
                            as * ((2 - w) * Cv1 + 2 * (1 + w) * Cv3) + (2 + w) * eB * Fb) /
                           sqrt(6.);
        const double Fa = (-(eC * (4 * LambdaD32 - (1 + w) * (2 * eta1 + 3 * eta3) - 3 * (w - 1) * (tau1 - tau2))) +
                           (-1 - w) * (1 + as * Ca1) - (w - 1) * eB * Fb) /
                          sqrt(6.);
        const double Ft1 = (-(eC * (4 * LambdaD32 + (1 + w) * (2 * eta1 + 3 * eta3) - 3 * (w - 1) * (tau1 - tau2))) +
                            (1 + w) * (1 + as * (Ct1 + (w - 1) * Ct2)) + (w - 1) * eB * Fb) /
                           sqrt(6.);
        const double Ft2 = (-(eC * (4 * LambdaD32 - (1 + w) * (2 * eta1 + 3 * eta3) - 3 * (w - 1) * (tau1 - tau2))) +
                            (-1 - w) * (1 + as * (Ct1 - (w - 1) * Ct3)) + (w - 1) * eB * Fb) /
                           sqrt(6.);
        const double Ft3 = (3 - eC * (-10 * eta1 - 4 * (w - 1) * eta2 + 5 * eta3 + (-1 + 4 * w) * tau1 + 5 * tau2) +
                            as * (3 * Ct1 - (2 - w) * Ct2 + 3 * Ct3) + 3 * eB * Fb) /
                           sqrt(6.);


        // Set elements
        result.element({0}) = Fsc;
        result.element({1}) = Fv1;
        result.element({2}) = Fv2;
        result.element({3}) = Fv3;
        result.element({4}) = Fa;
        result.element({5}) = Ft1;
        result.element({6}) = Ft2;
        result.element({7}) = Ft3;

        result *= LOIWtau;
    }

    std::unique_ptr<FormFactorBase> FFBtoD1BLR::clone(const std::string& label) {
        MAKE_CLONE(FFBtoD1BLR, label);
    }

} // namespace Hammer
