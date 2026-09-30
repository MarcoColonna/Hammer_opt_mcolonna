///
/// @file  FFBtoD2starLLSW.cc
/// @brief \f$ B \rightarrow D_2^* \f$ LLSW form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/FormFactors/LLSW/FFBtoD2starLLSW.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFBtoD2starLLSW::FFBtoD2starLLSW() {
        // Create tensor rank and dimensions
        IndexList dims = {8};
        string name{"FFBtoD2starLLSW"};

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

    void FFBtoD2starLLSW::defineSettings() {
        //_mFFErrNames = ;
        setPath(getFFErrPrefixGroup().get());
        setUnits("GeV");

        // Pole mass scheme as in EvtGen
        addSetting<double>("mb", 4.2);
        addSetting<double>("mc", 1.4);
        addSetting<double>("t1", 0.71);
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

    void FFBtoD2starLLSW::evalAtPSPoint(const vector<double>& point, const vector<double>& masses) {
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

        // LLSW parameters
        //  const double zBC = (*getSetting<double>("mc"))/(*getSetting<double>("mb"));
        const double eB = 1. / (*getSetting<double>("mb") * 2.);
        const double eC = 1. / (*getSetting<double>("mc") * 2.);
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
        const double LOIWtau = t1 + tp * (w - 1);

        // Form factors
        const double Ka1 = -1 - w - eC * (-((1 + w) * (2 * eta1 - eta3)) + (w - 1) * (tau1 - tau2)) - (w - 1) * eB * Fb;
        const double Ka2 = -2 * eC * (eta2 + tau1);
        const double Ka3 = 1 - eC * (2 * eta1 - 2 * eta2 - eta3 + tau1 + tau2) + eB * Fb;
        const double Kv = -1 - eC * (-2 * eta1 + eta3 + tau1 - tau2) - eB * Fb;


        // Set elements
        result.element({1}) = Ka1;
        result.element({2}) = Ka2;
        result.element({3}) = Ka3;
        result.element({4}) = Kv;

        result *= LOIWtau;
    }

    std::unique_ptr<FormFactorBase> FFBtoD2starLLSW::clone(const std::string& label) {
        MAKE_CLONE(FFBtoD2starLLSW, label);
    }

} // namespace Hammer
