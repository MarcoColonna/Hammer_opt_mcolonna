///
/// @file  FFBtoD1starLLSW.cc
/// @brief \f$ B \rightarrow D_1^* \f$ LLSW form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/FormFactors/LLSW/FFBtoD1starLLSW.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFBtoD1starLLSW::FFBtoD1starLLSW() {
        // Create tensor rank and dimensions
        IndexList dims = {8};
        string name{"FFBtoD1starLLSW"};

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

    void FFBtoD1starLLSW::defineSettings() {
        //_mFFErrNames = ;
        setPath(getFFErrPrefixGroup().get());
        setUnits("GeV");

        // Pole mass scheme as in EvtGen
        addSetting<double>("mb", 4.2);
        addSetting<double>("mc", 1.4);
        addSetting<double>("zt1", 0.68);
        addSetting<double>("ztp", -0.2);
        addSetting<double>("zeta1", 0.3);
        addSetting<double>("chi1", 0.);
        addSetting<double>("chi2", 0.);
        addSetting<double>("laB", 0.4);
        addSetting<double>("laS", 0.76);

        setInitialized();
    }

    void FFBtoD1starLLSW::evalAtPSPoint(const vector<double>& point, const vector<double>& masses) {
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
        const double zt1 = (*getSetting<double>("zt1"));
        const double ztp = (*getSetting<double>("ztp"));
        const double zeta1 = (*getSetting<double>("zeta1"));
        const double chi1 = (*getSetting<double>("chi1"));
        const double chi2 = (*getSetting<double>("chi2"));
        const double laB = (*getSetting<double>("laB"));
        const double laS = (*getSetting<double>("laS"));

        const double LambdaD12 = -laB + laS * w;
        const double Gb = (-(laB * (2 + w)) + laS * (1 + 2 * w)) / (1 + w) - 2 * (w - 1) * zeta1;
        const double LOIWzeta = zt1 + (w - 1) * ztp;

        // Form factors
        const double Gv1 = -1 + w + eC * (LambdaD12 - 2 * (w - 1) * chi1) - (1 + w) * eB * Gb;
        const double Gv2 = eC * (2 * zeta1 - 2 * chi2);
        const double Gv3 = -1 - eC * (LambdaD12 / (1 + w) + 2 * zeta1 - 2 * chi1 + 2 * chi2) + eB * Gb;
        const double Ga = 1 + eC * (LambdaD12 / (1 + w) - 2 * chi1) - eB * Gb;


        // Set elements
        result.element({1}) = Gv1;
        result.element({2}) = Gv2;
        result.element({3}) = Gv3;
        result.element({4}) = Ga;

        result *= LOIWzeta;
    }

    std::unique_ptr<FormFactorBase> FFBtoD1starLLSW::clone(const std::string& label) {
        MAKE_CLONE(FFBtoD1starLLSW, label);
    }

} // namespace Hammer
