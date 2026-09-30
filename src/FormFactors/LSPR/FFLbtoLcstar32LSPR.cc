///
/// @file  FFLbtoLcstar32LSPR.cc
/// @brief \f$ \Lambda_b \rightarrow L_c^*(2625) \f$ LSPR form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/FormFactors/LSPR/FFLbtoLcstar32LSPR.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFLbtoLcstar32LSPR::FFLbtoLcstar32LSPR() {
        // Create tensor rank and dimensions
        IndexList dims = {16};
        string name{"FFLbtoLcstar32LSPR"};

        setPrefix("LbtoLc*2625");
        addProcessSignature(-PID::LAMBDAB, {PID::LAMBDACSTAR32MINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_LBLCSTAR32})});

        setSignatureIndex();
    }

    void FFLbtoLcstar32LSPR::defineSettings() {
        //_mFFErrNames = ;
        setPath(getFFErrPrefixGroup().get());
        setUnits("GeV");

        FFLSPRBase::defineSettings(); // common parameters for LSPR are defined in FFLSPRBase

        setInitialized();
    }

    void FFLbtoLcstar32LSPR::eval(const Particle& parent, const ParticleList& daughters,
                                  const ParticleList& /*references*/) {
        // Momenta
        const FourMomentum& pLbmes = parent.momentum();
        const FourMomentum& pLcmes = daughters[0].momentum();
        // const FourMomentum& pTau = daughters[2].momentum();

        // kinematic objects
        const double Mb = pLbmes.mass();
        const double Mc = pLcmes.mass();
        // const double Mt = pTau.mass();
        const double Sqq = Mb * Mb + Mc * Mc - 2. * (pLbmes * pLcmes);

        evalAtPSPoint({Sqq}, {Mb, Mc});
    }

    void FFLbtoLcstar32LSPR::evalAtPSPoint(const vector<double>& point, const vector<double>& masses) {
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
        const double w2m1 = w * w - 1;

        // 1S scheme
        const double aS = (*getSetting<double>("as")) / pi;
        const double mb1S = unitres * (*getSetting<double>("mb1S"));
        const double dmbmc = unitres * (*getSetting<double>("dmbmc"));
        const double mc1S = mb1S - dmbmc;

        // Explicit expansion to O(eps)
        const double eps = 2. * pow(aS * pi, 2.) / 9.;
        const double zBC = mc1S / mb1S;
        //        const double mb = mb1S*(1 + 2.*pow(aS*pi,2.)/9.);
        //        const double mc = mb - dmbmc;
        const double eB = 1. / (2. * mb1S) * (1 - eps);
        const double eC = 1. / (2. * mc1S) * (1. - eps / zBC);

        //(Sub)leading IW function parameters
        const double s1 = (*getSetting<double>("s1"));
        const double sp1 = (*getSetting<double>("sp1"));
        const double s11 = unitres * (*getSetting<double>("s11"));
        const double s1p1 = unitres * (*getSetting<double>("s1p1"));
        const double pc1 = unitres * (*getSetting<double>("pc1"));
        const double pb1 = unitres * (*getSetting<double>("pb1"));

        const double Lbar = unitres * (*getSetting<double>("barLambda"));
        const double Lbarp = unitres * (*getSetting<double>("barLambdap"));

        // Leading and subleading IW functions
        const double sigmaIW = s1 + s1 * sp1 * (w - 1);
        const double sigma1 = s11 + s1p1 * (w - 1);
        const double pc = pc1;
        const double pb = pb1;
        ////Tensor current kernel
        const double lam = 0.;

        // QCD correction functions
        const double Cs = _asCorrections.CS(w, zBC);
        const double Cps = _asCorrections.CP(w, zBC);
        const double Cv1 = _asCorrections.CV1(w, zBC);
        const double Cv2 = _asCorrections.CV2(w, zBC);
        const double Cv3 = _asCorrections.CV3(w, zBC);
        const double Ca1 = _asCorrections.CA1(w, zBC);
        const double Ca2 = _asCorrections.CA2(w, zBC);
        const double Ca3 = _asCorrections.CA3(w, zBC);
        const double Ct1 = _asCorrections.CT1(w, zBC);
        const double Ct2 = _asCorrections.CT2(w, zBC);
        const double Ct3 = _asCorrections.CT3(w, zBC);

        // Create hatted FFs
        const double lS = 1 + aS * Cs + eB * (Lbar - Lbarp - pb + sigma1 * (-1 + w)) + eC * (pc + sigma1 * (-1 + w));
        const double lP = 1 + aS * Cps + eB * (Lbar + Lbarp + pb + sigma1 * (-1 - w)) + eC * (pc + sigma1 * (1 + w));
        const double lV1 = 1 + aS * Cv1 + eB * (Lbar + Lbarp + pb + sigma1 * (-1 - w)) + eC * (pc + sigma1 * (1 + w));
        const double lV2 = aS * Cv2 - 2 * eC * sigma1;
        const double lV3 = aS * Cv3 + eB * (-2 * Lbarp - 2 * pb + 2 * sigma1 * w);
        const double lV4 = eB * (-2 * Lbar + 2 * pb * (-1 + w) + 2 * Lbarp * w + sigma1 * (2 - 2 * w * w));
        const double lA1 = 1 + aS * Ca1 + eB * (Lbar - Lbarp - pb + sigma1 * (-1 + w)) + eC * (pc + sigma1 * (-1 + w));
        const double lA2 = aS * Ca2 - 2 * eC * sigma1;
        const double lA3 = aS * Ca3 + eB * (2 * Lbarp + 2 * pb - 2 * sigma1 * w);
        const double lA4 = eB * (-2 * pb * (1 + w) + 2 * (Lbar - Lbarp * w) + 2 * sigma1 * w2m1);
        const double lT1 =
            1 - lam + aS * Ct1 + eB * (Lbar - Lbarp - pb + sigma1 * (-1 + w)) + eC * (pc + sigma1 * (-1 + w));
        const double lT2 = aS * Ct2 - 2 * eC * sigma1;
        const double lT3 = lam + aS * Ct3 + eB * (2 * Lbarp + 2 * pb - 2 * sigma1 * w);
        const double lT4 = -(lam * (1 + w)) + eB * (-2 * pb * (1 + w) + 2 * (Lbar - Lbarp * w) + 2 * sigma1 * w2m1);
        // const double lT5 = lam;
        const double lT6 = lam + 4 * eB * pb;
        const double lT7 = 0;

        // set elements
        result.element({0}) = lS;
        result.element({1}) = lP;
        result.element({2}) = lV1;
        result.element({3}) = lV2;
        result.element({4}) = lV3;
        result.element({5}) = lV4;
        result.element({6}) = lA1;
        result.element({7}) = lA2;
        result.element({8}) = lA3;
        result.element({9}) = lA4;
        result.element({10}) = lT1;
        result.element({11}) = lT2;
        result.element({12}) = lT3;
        result.element({13}) = lT4;
        result.element({14}) = lT6;
        result.element({15}) = lT7;

        result *= sigmaIW;
    }

    std::unique_ptr<FormFactorBase> FFLbtoLcstar32LSPR::clone(const std::string& label) {
        MAKE_CLONE(FFLbtoLcstar32LSPR, label);
    }

} // namespace Hammer
