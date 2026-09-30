///
/// @file  FFBtoDstarCLN.cc
/// @brief \f$ B \rightarrow D^* \f$ CLN form factors
/// @brief Ported from EvtGen HQET2
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/FormFactors/CLN/FFBtoDstarCLN.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFBtoDstarCLN::FFBtoDstarCLN() {
        // Create tensor rank and dimensions
        IndexList dims = {8};
        string name{"FFBtoDstarCLN"};

        setPrefix("BtoD*");
        addProcessSignature(PID::BPLUS, {-PID::DSTAR});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BDSTAR})});

        addProcessSignature(PID::BZERO, {PID::DSTARMINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BDSTAR})});

        setPrefix("BstoDs*");
        addProcessSignature(PID::BS, {PID::DSSTARMINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BSDSSTAR})});

        setSignatureIndex();
    }

    void FFBtoDstarCLN::defineSettings() {
        //_mFFErrNames = ;
        setPath(getFFErrPrefixGroup().get());
        setUnits("GeV");

        addSetting<double>("a", 1.0); // zero recoil expansion
        addSetting<double>("RhoSq", 1.207);
        addSetting<double>("F1", 0.908);
        addSetting<double>("R1", 1.401);
        addSetting<double>("R2", 0.854);
        addSetting<double>("R0", 1.15);
        addSetting<double>("as", 0.26);     // mbmc scale
        addSetting<double>("la", 0.48);     // CLN
        addSetting<double>("mcOnmb", 0.29); // CLN
        setInitialized();
    }

    void FFBtoDstarCLN::evalAtPSPoint(const vector<double>& point, const vector<double>& masses) {
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
        const double R0par = *getSetting<double>("R0");
        const double R1par = *getSetting<double>("R1");
        const double R2par = *getSetting<double>("R2");
        const double F1par = *getSetting<double>("F1");
        const double as = *getSetting<double>("as") / pi;
        const double LambdaBar = *getSetting<double>("la");
        const double zR = *getSetting<double>("mcOnmb");

        const double rC = Mc / Mb;
        const double rC2 = rC * rC;

        // MSbar masses
        const double mbMu = Mb - unitres * LambdaBar - Mb * as * (4. / 3. + log(zR));
        const double mcMu = Mc - unitres * LambdaBar - Mc * as * (4. / 3. - log(zR));

        // define the variables
        const double Ha1 = F1par * (1. - 8. * RhoSq * zCon + (53. * RhoSq - 15.) * pow(zCon, 2.) -
                                    (231. * RhoSq - 91.) * pow(zCon, 3.));
        const double R0 = R0par - 0.11 * (w - 1.) + 0.01 * pow(w - 1., 2.);
        const double R1 = R1par - 0.12 * (w - 1.) + 0.05 * pow(w - 1., 2.);
        const double R2 = R2par + 0.11 * (w - 1.) - 0.06 * pow(w - 1., 2.);

        const double Ff = (Ha1 * (Mb2 * pow(1 + rC, 2) - Sqq)) / (2. * Mb * sqrt(rC));
        const double Fg = (Ha1 * R1) / (2. * Mb * sqrt(rC));
        const double Fm = -(Ha1 * (R2 * (-1 + rC2) - 2 * rC * (-1 + R0 + R0 * rC - w))) /
                          (2. * Mb * sqrt(rC) * (1 + rC2 - 2 * rC * w));
        const double Fp = -(Ha1 * R2) / (2. * Mb * sqrt(rC));
        // Pseudoscalar from eqn of motion
        const double Fps = -(Ff + Mb2 * (Fp * (1 - rC2) + Fm * (1 + rC2 - 2 * rC * w))) / (mbMu + mcMu);

        // set elements, mapping to amplitude FF basis
        // Fps
        result.element({0}) = Fps;
        // Ff
        result.element({1}) = Ff;
        // Fg
        result.element({2}) = Fg;
        // Fm
        result.element({3}) = Fm;
        // Fp
        result.element({4}) = Fp;
        // Fzt
        // result.element({5}) = 0;
        // Fmt
        // result.element({6}) = 0;
        // Fpt
        // result.element({7}) = 0;
    }

    std::unique_ptr<FormFactorBase> FFBtoDstarCLN::clone(const std::string& label) {
        MAKE_CLONE(FFBtoDstarCLN, label);
    }

} // namespace Hammer
