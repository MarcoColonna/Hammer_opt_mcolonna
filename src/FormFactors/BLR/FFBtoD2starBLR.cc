///
/// @file  FFBtoD2starBLR.cc
/// @brief \f$ B \rightarrow D_2^* \f$ BLR form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/FormFactors/BLR/FFBtoD2starBLR.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFBtoD2starBLR::FFBtoD2starBLR() {
        // Create tensor rank and dimensions
        IndexList dims = {8};
        string name{"FFBtoD2starBLR"};

        setPrefix("BtoD**2*");
        addProcessSignature(PID::BPLUS, {-PID::DSSD2STAR});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BDSSD2STAR})});

        addProcessSignature(PID::BZERO, {PID::DSSD2STARMINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BDSSD2STAR})});

        setPrefix("BstoDs**2*");
        addProcessSignature(PID::BS, {PID::DSSDS2STARMINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BSDSSDS2STAR})});

        setSignatureIndex();
    }

    void FFBtoD2starBLR::defineSettings() {
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

    void FFBtoD2starBLR::evalAtPSPoint(const vector<double>& point, const vector<double>& masses) {
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


        // Set elements
        result.element({0}) = Kp;
        result.element({1}) = Ka1;
        result.element({2}) = Ka2;
        result.element({3}) = Ka3;
        result.element({4}) = Kv;
        result.element({5}) = Kt1;
        result.element({6}) = Kt2;
        result.element({7}) = Kt3;

        result *= LOIWtau;
    }

    std::unique_ptr<FormFactorBase> FFBtoD2starBLR::clone(const std::string& label) {
        MAKE_CLONE(FFBtoD2starBLR, label);
    }

} // namespace Hammer
