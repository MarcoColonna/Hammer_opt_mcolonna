///
/// @file  FFBtoD1LLSW.cc
/// @brief \f$ B \rightarrow D_1 \f$ LLSW form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/FormFactors/LLSW/FFBtoD1LLSW.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFBtoD1LLSW::FFBtoD1LLSW() {
        // Create tensor rank and dimensions
        IndexList dims = {8};
        string name{"FFBtoD1LLSW"};

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

    void FFBtoD1LLSW::defineSettings() {
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

    void FFBtoD1LLSW::evalAtPSPoint(const vector<double>& point, const vector<double>& masses) {
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

        const double LambdaD32 = -laB + laP * w;
        const double Fb = laB + laP - tau2 - tau1 * (1 + 2 * w);
        const double LOIWtau = t1 + tp * (w - 1);

        // Form factors
        const double Fv1 =
            (1 - w * w - eC * (4 * (1 + w) * LambdaD32 - (-1 + w * w) * (2 * eta1 + 3 * eta3 + 3 * tau1 - 3 * tau2)) -
             (-1 + w * w) * eB * Fb) /
            sqrt(6.);
        const double Fv2 =
            (-3 - eC * (10 * eta1 + 4 * (w - 1) * eta2 - 5 * eta3 + (-1 + 4 * w) * tau1 + 5 * tau2) - 3 * eB * Fb) /
            sqrt(6.);
        const double Fv3 = (-2 + w +
                            eC * (4 * LambdaD32 - 2 * (6 + w) * eta1 - 4 * (w - 1) * eta2 - (-2 + 3 * w) * eta3 +
                                  (2 + w) * tau1 + (2 + 3 * w) * tau2) +
                            (2 + w) * eB * Fb) /
                           sqrt(6.);
        const double Fa =
            (-1 - w - eC * (4 * LambdaD32 - (1 + w) * (2 * eta1 + 3 * eta3) - 3 * (w - 1) * (tau1 - tau2)) -
             (w - 1) * eB * Fb) /
            sqrt(6.);


        // Set elements
        result.element({1}) = Fv1;
        result.element({2}) = Fv2;
        result.element({3}) = Fv3;
        result.element({4}) = Fa;

        result *= LOIWtau;
    }

    std::unique_ptr<FormFactorBase> FFBtoD1LLSW::clone(const std::string& label) {
        MAKE_CLONE(FFBtoD1LLSW, label);
    }

} // namespace Hammer
