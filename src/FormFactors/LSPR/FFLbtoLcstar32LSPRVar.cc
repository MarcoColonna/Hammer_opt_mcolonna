///
/// @file  FFLbtoLcstar32LSPRVar.cc
/// @brief \f$ \Lambda_b \rightarrow L_c^*(2625) \f$ LSPRVar form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/FormFactors/LSPR/FFLbtoLcstar32LSPRVar.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFLbtoLcstar32LSPRVar::FFLbtoLcstar32LSPRVar() {
        // Create tensor rank and dimensions
        IndexList dims = {16, 7}; // size is _mFFErrNames + 1 to include central values in zeroth component
        string name{"FFLbtoLcstar32LSPRVar"};

        setPrefix("LbtoLc*2625");
        _mFFErrLabel = FF_LBLCSTAR32_VAR;
        addProcessSignature(-PID::LAMBDAB, {PID::LAMBDACSTAR32MINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_LBLCSTAR32, FF_LBLCSTAR32_VAR})});

        setSignatureIndex();
    }

    void FFLbtoLcstar32LSPRVar::defineSettings() {
        //_mFFErrNames = ;
        setPath(getFFErrPrefixGroup().get());
        setUnits("GeV");

        FFLSPRVarBase::defineSettings(); // common parameters for LSPRVar are defined in FFLSPRVarBase

        setInitialized();
    }

    void FFLbtoLcstar32LSPRVar::eval(const Particle& parent, const ParticleList& daughters,
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

    void FFLbtoLcstar32LSPRVar::evalAtPSPoint(const vector<double>& point, const vector<double>& masses) {
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
        const double wm1sq = (w - 1.) * (w - 1.);

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

        const vector<vector<double>>& sliwmat = (*getSetting<vector<vector<double>>>("sliwmatrix"));

        const double Lbar = unitres * (*getSetting<double>("barLambda"));
        const double Lbarp = unitres * (*getSetting<double>("barLambdap"));

        // Leading and subleading IW functions
        const double sigmaIW = s1 + s1 * sp1 * (w - 1);
        const double sigmaIWp = s1 * (w - 1) / sigmaIW;
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
        result.element({0, 0}) = lS;
        result.element({1, 0}) = lP;
        result.element({2, 0}) = lV1;
        result.element({3, 0}) = lV2;
        result.element({4, 0}) = lV3;
        result.element({5, 0}) = lV4;
        result.element({6, 0}) = lA1;
        result.element({7, 0}) = lA2;
        result.element({8, 0}) = lA3;
        result.element({9, 0}) = lA4;
        result.element({10, 0}) = lT1;
        result.element({11, 0}) = lT2;
        result.element({12, 0}) = lT3;
        result.element({13, 0}) = lT4;
        result.element({14, 0}) = lT6;
        result.element({15, 0}) = lT7;

        // Now create remaining FF tensor entries
        // n indexes rows ie FFs; idx indexes the columns that contract with delta

        // Linearization of hatted FFs
        vector<vector<double>> FFmat{
            {lS / s1, lS * sigmaIWp, (eB + eC) * (-1 + w), (eB + eC) * wm1sq, eC, -eB},
            {lP / s1, lP * sigmaIWp, -((eB - eC) * (1 + w)), -((eB - eC) * w2m1), eC, eB},
            {lV1 / s1, lV1 * sigmaIWp, -((eB - eC) * (1 + w)), -((eB - eC) * w2m1), eC, eB},
            {lV2 / s1, lV2 * sigmaIWp, -2 * eC, -2 * eC * (-1 + w), 0, 0},
            {lV3 / s1, lV3 * sigmaIWp, 2 * eB * w, 2 * eB * (-1 + w) * w, 0, -2 * eB},
            {lV4 / s1, lV4 * sigmaIWp, -2 * eB * w2m1, -2 * eB * (-1 + w) * w2m1, 0, 2 * eB * (-1 + w)},
            {lA1 / s1, lA1 * sigmaIWp, (eB + eC) * (-1 + w), (eB + eC) * wm1sq, eC, -eB},
            {lA2 / s1, lA2 * sigmaIWp, -2 * eC, -2 * eC * (-1 + w), 0, 0},
            {lA3 / s1, lA3 * sigmaIWp, -2 * eB * w, -2 * eB * (-1 + w) * w, 0, 2 * eB},
            {lA4 / s1, lA4 * sigmaIWp, 2 * eB * w2m1, 2 * eB * (-1 + w) * w2m1, 0, -2 * eB * (1 + w)},
            {lT1 / s1, lT1 * sigmaIWp, (eB + eC) * (-1 + w), (eB + eC) * wm1sq, eC, -eB},
            {lT2 / s1, lT2 * sigmaIWp, -2 * eC, -2 * eC * (-1 + w), 0, 0},
            {lT3 / s1, lT3 * sigmaIWp, -2 * eB * w, -2 * eB * (-1 + w) * w, 0, 2 * eB},
            {lT4 / s1, lT4 * sigmaIWp, 2 * eB * w2m1, 2 * eB * (-1 + w) * w2m1, 0, -2 * eB * (1 + w)},
            {lT6 / s1, lT6 * sigmaIWp, 0, 0, 0, 4 * eB},
            {lT7 / s1, lT7 * sigmaIWp, 0, 0, 0, 0}};

        for (size_t n = 0; n < 16; ++n) {
            auto ni = static_cast<IndexType>(n);
            for (size_t idx = 0; idx < 6; ++idx) {
                complex<double> entry = 0;
                auto idxi = static_cast<IndexType>(idx + 1u);
                for (size_t idx1 = 0; idx1 < 6; ++idx1) {
                    double unitsl =
                        (idx1 >= 2) ? unitres : 1.; // delta's for subleading IW functions to internal class units (GeV)
                    entry += sliwmat[idx1][idx] * FFmat[n][idx1] * unitsl;
                }
                result.element({ni, idxi}) = entry;
            }
        }

        result *= sigmaIW;
    }

    std::unique_ptr<FormFactorBase> FFLbtoLcstar32LSPRVar::clone(const std::string& label) {
        MAKE_CLONE(FFLbtoLcstar32LSPRVar, label);
    }

} // namespace Hammer
