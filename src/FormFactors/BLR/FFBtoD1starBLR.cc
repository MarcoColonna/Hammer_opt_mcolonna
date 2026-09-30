///
/// @file  FFBtoD1starBLR.cc
/// @brief \f$ B \rightarrow D_1^* \f$ BLR form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/FormFactors/BLR/FFBtoD1starBLR.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFBtoD1starBLR::FFBtoD1starBLR() {
        // Create tensor rank and dimensions
        IndexList dims = {8};
        string name{"FFBtoD1starBLR"};

        setPrefix("BtoD**1*");
        addProcessSignature(PID::BPLUS, {-PID::DSSD1STAR});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BDSSD1STAR})});

        addProcessSignature(PID::BZERO, {PID::DSSD1STARMINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BDSSD1STAR})});

        setPrefix("BstoDs**1*");
        addProcessSignature(PID::BS, {PID::DSSDS1STARMINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BSDSSDS1STAR})});

        setSignatureIndex();
    }

    void FFBtoD1starBLR::defineSettings() {
        //_mFFErrNames = ;
        setPath(getFFErrPrefixGroup().get());
        setUnits("GeV");

        // 1S scheme: mb1S = 4710, mbBar = 5313, lambda1 = -3 10^5 MeV^2, delta mb-mc = 3400, alpha_s = 26/100
        addSetting<double>("as", 0.26);
        addSetting<double>("mb", 4.710);
        addSetting<double>("mc", 4.710 - 3.400);
        addSetting<double>("zt1", 0.7);
        addSetting<double>("ztp", 0.2);
        addSetting<double>("zeta1", 0.6);
        addSetting<double>("chi1", 0.);
        addSetting<double>("chi2", 0.);
        addSetting<double>("laB", 0.4);
        addSetting<double>("laS", 0.76);

        setInitialized();
    }

    void FFBtoD1starBLR::evalAtPSPoint(const vector<double>& point, const vector<double>& masses) {
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
        const double zt1 = (*getSetting<double>("zt1"));
        const double ztp = (*getSetting<double>("ztp"));
        const double zeta1 = (*getSetting<double>("zeta1"));
        const double chi1 = (*getSetting<double>("chi1"));
        const double chi2 = (*getSetting<double>("chi2"));
        const double laB = (*getSetting<double>("laB"));
        const double laS = (*getSetting<double>("laS"));

        const double LambdaD12 = -laB + laS * w;
        const double Gb = (-(laB * (2 + w)) + laS * (1 + 2 * w)) / (1 + w) - 2 * (w - 1) * zeta1;
        const double LOIWzeta = zt1 + (w - 1) * zt1 * ztp;

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
        const double Gs =
            1 - eC * (LambdaD12 / (1 + w) - 2 * (w - 1) * zeta1 + 2 * chi1 - 2 * (1 + w) * chi2) + as * Cs - eB * Gb;
        const double Gv1 = eC * (LambdaD12 - 2 * (w - 1) * chi1) + (w - 1) * (1 + as * Cv1) - (1 + w) * eB * Gb;
        const double Gv2 = eC * (2 * zeta1 - 2 * chi2) - as * Cv2;
        const double Gv3 =
            -1 - eC * (LambdaD12 / (1 + w) + 2 * zeta1 - 2 * chi1 + 2 * chi2) - as * (Cv1 + Cv3) + eB * Gb;
        const double Ga = 1 + eC * (LambdaD12 / (1 + w) - 2 * chi1) + as * Ca1 - eB * Gb;
        const double Gt1 = -1 + eC * (LambdaD12 / (1 + w) + 2 * chi1) - as * (Ct1 + (w - 1) * Ct2) + eB * Gb;
        const double Gt2 = 1 + eC * (LambdaD12 / (1 + w) - 2 * chi1) + as * (Ct1 - (w - 1) * Ct3) + eB * Gb;
        const double Gt3 = eC * (2 * zeta1 + 2 * chi2) - as * Ct2;


        // Set elements
        result.element({0}) = Gs;
        result.element({1}) = Gv1;
        result.element({2}) = Gv2;
        result.element({3}) = Gv3;
        result.element({4}) = Ga;
        result.element({5}) = Gt1;
        result.element({6}) = Gt2;
        result.element({7}) = Gt3;

        result *= LOIWzeta;
    }

    std::unique_ptr<FormFactorBase> FFBtoD1starBLR::clone(const std::string& label) {
        MAKE_CLONE(FFBtoD1starBLR, label);
    }

} // namespace Hammer
