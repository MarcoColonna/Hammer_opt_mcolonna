///
/// @file  FFBtoDCLN.cc
/// @brief \f$ B \rightarrow D \f$ CLN form factors
/// @brief Ported from EvtGen HQET3 (custom class; F Bernlochner)
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/FormFactors/CLN/FFBtoDCLN.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFBtoDCLN::FFBtoDCLN() {
        // Create tensor rank and dimensions
        IndexList dims = {4};
        string name{"FFBtoDCLN"};

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

    void FFBtoDCLN::defineSettings() {
        //_mFFErrNames = ;
        setPath(getFFErrPrefixGroup().get());
        setUnits("GeV");

        addSetting<double>("a", 1.0);       // zero recoil expansion
        addSetting<double>("RhoSq", 1.186); // taken from HFLAV
        addSetting<double>("G1", 1.082);    // this value was taken from a BCL-Lattice fit!
        addSetting<double>("Delta", 1.);    // large errors!
        addSetting<double>("as", 0.26);     // mbmc scale
        addSetting<double>("la", 0.48);     // CLN
        addSetting<double>("mcOnmb", 0.29); // CLN
        addSetting<bool>("isHQET2", false); // emulate EvtGen HQET2 class
        setInitialized();
    }

    void FFBtoDCLN::evalAtPSPoint(const vector<double>& point, const vector<double>& masses) {
        Tensor& result = getTensor();
        result.clearData();

        if (!isInitialized()) {
            MSG_WARNING("Warning, Settings have not been defined!");
        }

        double Mb = 0.;
        double Mc = 0.;
        double unitres = 1.;
        tie(Mb, Mc, unitres) = getParentDaughterHadMasses(masses);
        const double Mb2 = Mb * Mb;

        double Sqq = point[0];

        // const double sqMbMc = sqrt(Mb*Mc);
        double w = getW(Sqq, Mb, Mc);
        // safety measure if w==1.0
        if (isZero(w - 1.0)) {
            w += 1e-6;
        }
        const double a = *getSetting<double>("a");
        const double zCon = (sqrt(w + 1) - sqrt2 * a) / (sqrt(w + 1) + sqrt2 * a);

        const double RhoSq = *getSetting<double>("RhoSq");
        const double G1par = *getSetting<double>("G1");
        const double Delta = *getSetting<double>("Delta");
        const double as = *getSetting<double>("as") / pi;
        const double LambdaBar = *getSetting<double>("la");
        const double zR = *getSetting<double>("mcOnmb");
        const bool isHQET2 = *getSetting<bool>("isHQET2");

        const double rC = Mc / Mb;

        // MSbar masses
        const double mbMu = Mb - unitres * LambdaBar - Mb * as * (4. / 3. + log(zR));
        const double mcMu = Mc - unitres * LambdaBar - Mc * as * (4. / 3. - log(zR));

        // define the variables. Added normalization to Fp.
        const double V1wOnV1 =
            1. - 8. * RhoSq * zCon + (51. * RhoSq - 10.) * zCon * zCon - (252. * RhoSq - 84.) * zCon * zCon * zCon;
        const double Fp = ((1. + rC) / (2. * sqrt(rC)) * V1wOnV1) * G1par;
        const double S1wOnV1w =
            1 +
            Delta *
                (-0.019 + 0.041 * (w - 1.) -
                 0.015 * pow(w - 1.,
                             2.)); // From 1005.4306 as in EVtGen HQET3 (custom class; F Bernlochner); may need updating
        const double Fz = isHQET2 ? 0.0 : (sqrt(rC) / (1. + rC) * (1 + w) * S1wOnV1w * V1wOnV1) * G1par;
        // Scalar from eqn of motion
        const double Fs = Fz * Mb2 * (1. - rC * rC) / (mbMu - mcMu);
        // double Fz =
        // (pow(2.*sqrt(rC)/(1.+rC),2.)*(1+w)*1.0036*(1.-0.0068*(w-1)+0.0017*pow((w-1),2.)-0.0013*pow((w-1),3.)))/(2.);


        // set elements, mapping to amplitude FF basis
        // Fs
        result.element({0}) = Fs;
        // Fz
        result.element({1}) = Fz;
        // Fp
        result.element({2}) = Fp;
        // Ft
        // result.element({3}) = 0;
    }

    std::unique_ptr<FormFactorBase> FFBtoDCLN::clone(const std::string& label) {
        MAKE_CLONE(FFBtoDCLN, label);
    }

} // namespace Hammer
