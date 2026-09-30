///
/// @file  FFLbtoLcstar12LSPRVar.cc
/// @brief \f$ \Lambda_b \rightarrow L_c^*(2595) \f$ LSPRVar form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/FormFactors/LSPR/FFLbtoLcstar12LSPRVar.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFLbtoLcstar12LSPRVar::FFLbtoLcstar12LSPRVar() {
        // Create tensor rank and dimensions
        IndexList dims = {12, 7}; // size is _mFFErrNames + 1 to include central values in zeroth component
        string name{"FFLbtoLcstar12LSPRVar"};

        setPrefix("LbtoLc*2595");
        _mFFErrLabel = FF_LBLCSTAR12_VAR;
        addProcessSignature(-PID::LAMBDAB, {PID::LAMBDACSTAR12MINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_LBLCSTAR12, FF_LBLCSTAR12_VAR})});

        setSignatureIndex();
    }

    void FFLbtoLcstar12LSPRVar::defineSettings() {
        //_mFFErrNames = ;
        setPath(getFFErrPrefixGroup().get());
        setUnits("GeV");

        FFLSPRVarBase::defineSettings(); // common parameters for LSPRVar are defined in FFLSPRVarBase

        setInitialized();
    }

    void FFLbtoLcstar12LSPRVar::eval(const Particle& parent, const ParticleList& daughters,
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

    void FFLbtoLcstar12LSPRVar::evalAtPSPoint(const vector<double>& point, const vector<double>& masses) {
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
        const double dS =
            (1 + w) / sqrt3 + (aS * Cs * (1 + w)) / sqrt3 +
            eC * ((-2 * pc * (1 + w)) / sqrt3 + sqrt3 * (-Lbar + Lbarp * w) - (2 * sigma1 * w2m1) / sqrt3) +
            eB * ((2 * pb * (1 + w)) / sqrt3 + (Lbar * (-2 + w) + Lbarp * (-1 + 2 * w)) / sqrt3 -
                  (2 * sigma1 * w2m1) / sqrt3);
        const double dP =
            (-1 + w) / sqrt3 + (aS * Cps * (-1 + w)) / sqrt3 +
            eC * ((-2 * pc * (-1 + w)) / sqrt3 + sqrt3 * (-Lbar + Lbarp * w) - (2 * sigma1 * w2m1) / sqrt3) +
            eB * ((-2 * pb * (-1 + w)) / sqrt3 + (Lbar * (2 + w) - Lbarp * (1 + 2 * w)) / sqrt3 +
                  (2 * sigma1 * w2m1) / sqrt3);
        const double dV1 =
            (-1 + w) / sqrt3 + (aS * Cv1 * (-1 + w)) / sqrt3 + (eB * (-Lbarp + Lbar * w)) / sqrt3 +
            eC * ((-2 * pc * (-1 + w)) / sqrt3 + sqrt3 * (-Lbar + Lbarp * w) - (2 * sigma1 * w2m1) / sqrt3);
        const double dV2 = -2 / sqrt3 + (4 * eC * pc) / sqrt3 - (aS * (2 * Cv1 + Cv2 + Cv2 * w)) / sqrt3 +
                           eB * ((-2 * (Lbar + Lbarp)) / sqrt3 - (2 * pb) / sqrt3 + (2 * sigma1 * (1 + w)) / sqrt3);
        const double dV3 = -((aS * Cv3 * (1 + w)) / sqrt3) +
                           eB * ((2 * (Lbar + Lbarp)) / sqrt3 - (2 * pb) / sqrt3 - (2 * sigma1 * (1 + w)) / sqrt3);
        const double dA1 =
            (1 + w) / sqrt3 + (aS * Ca1 * (1 + w)) / sqrt3 + (eB * (-Lbarp + Lbar * w)) / sqrt3 +
            eC * ((-2 * pc * (1 + w)) / sqrt3 + sqrt3 * (-Lbar + Lbarp * w) - (2 * sigma1 * w2m1) / sqrt3);
        const double dA2 = -2 / sqrt3 + (4 * eC * pc) / sqrt3 +
                           eB * ((2 * (-Lbar + Lbarp)) / sqrt3 + (2 * pb) / sqrt3 - (2 * sigma1 * (-1 + w)) / sqrt3) +
                           (aS * (-2 * Ca1 + Ca2 - Ca2 * w)) / sqrt3;
        const double dA3 = eB * ((2 * (-Lbar + Lbarp)) / sqrt3 - (2 * pb) / sqrt3 - (2 * sigma1 * (-1 + w)) / sqrt3) -
                           (aS * Ca3 * (-1 + w)) / sqrt3;
        const double dT1 =
            (1 + w) / sqrt3 + (aS * Ct1 * (1 + w)) / sqrt3 +
            eC * ((-2 * pc * (1 + w)) / sqrt3 + sqrt3 * (-Lbar + Lbarp * w) - (2 * sigma1 * w2m1) / sqrt3) +
            eB * ((-2 * pb * (1 + w)) / sqrt3 + (Lbar * (2 + w) - Lbarp * (1 + 2 * w)) / sqrt3 +
                  (2 * sigma1 * w2m1) / sqrt3);
        const double dT2 = -2 / sqrt3 + (4 * eC * pc) / sqrt3 +
                           eB * ((2 * (-Lbar + Lbarp)) / sqrt3 + (2 * pb) / sqrt3 - (2 * sigma1 * (-1 + w)) / sqrt3) +
                           (aS * (-2 * Ct1 + Ct2 - Ct2 * w)) / sqrt3;
        const double dT3 = eB * ((2 * (-Lbar + Lbarp)) / sqrt3 + (2 * pb) / sqrt3 - (2 * sigma1 * (-1 + w)) / sqrt3) -
                           (aS * Ct3 * (-1 + w)) / sqrt3;
        const double dT4 =
            (-2 * aS * Ct3) / sqrt3 + eB * ((-4 * Lbarp) / sqrt3 - (4 * pb) / sqrt3 + (4 * sigma1 * w) / sqrt3);

        // set elements
        result.element({0, 0}) = dS;
        result.element({1, 0}) = dP;
        result.element({2, 0}) = dV1;
        result.element({3, 0}) = dV2;
        result.element({4, 0}) = dV3;
        result.element({5, 0}) = dA1;
        result.element({6, 0}) = dA2;
        result.element({7, 0}) = dA3;
        result.element({8, 0}) = dT1;
        result.element({9, 0}) = dT2;
        result.element({10, 0}) = dT3;
        result.element({11, 0}) = dT4;

        // Now create remaining FF tensor entries
        // n indexes rows ie FFs; idx indexes the columns that contract with delta

        // Linearization of hatted FFs
        vector<vector<double>> FFmat{
            {dS / s1, dS * sigmaIWp, (-2 * (eB + eC) * w2m1) / sqrt3, (-2 * (eB + eC) * wm1sq * (1 + w)) / sqrt3,
             (-2 * eC * (1 + w)) / sqrt3, (2 * eB * (1 + w)) / sqrt3},
            {dP / s1, dP * sigmaIWp, (2 * (eB - eC) * w2m1) / sqrt3, (2 * (eB - eC) * wm1sq * (1 + w)) / sqrt3,
             (-2 * eC * (-1 + w)) / sqrt3, (-2 * eB * (-1 + w)) / sqrt3},
            {dV1 / s1, dV1 * sigmaIWp, (-2 * eC * w2m1) / sqrt3, (-2 * eC * wm1sq * (1 + w)) / sqrt3,
             (-2 * eC * (-1 + w)) / sqrt3, 0},
            {dV2 / s1, dV2 * sigmaIWp, (2 * eB * (1 + w)) / sqrt3, (2 * eB * w2m1) / sqrt3, (4 * eC) / sqrt3,
             (-2 * eB) / sqrt3},
            {dV3 / s1, dV3 * sigmaIWp, (-2 * eB * (1 + w)) / sqrt3, (-2 * eB * w2m1) / sqrt3, 0, (-2 * eB) / sqrt3},
            {dA1 / s1, dA1 * sigmaIWp, (-2 * eC * w2m1) / sqrt3, (-2 * eC * wm1sq * (1 + w)) / sqrt3,
             (-2 * eC * (1 + w)) / sqrt3, 0},
            {dA2 / s1, dA2 * sigmaIWp, (-2 * eB * (-1 + w)) / sqrt3, (-2 * eB * wm1sq) / sqrt3, (4 * eC) / sqrt3,
             (2 * eB) / sqrt3},
            {dA3 / s1, dA3 * sigmaIWp, (-2 * eB * (-1 + w)) / sqrt3, (-2 * eB * wm1sq) / sqrt3, 0, (-2 * eB) / sqrt3},
            {dT1 / s1, dT1 * sigmaIWp, (2 * (eB - eC) * w2m1) / sqrt3, (2 * (eB - eC) * wm1sq * (1 + w)) / sqrt3,
             (-2 * eC * (1 + w)) / sqrt3, (-2 * eB * (1 + w)) / sqrt3},
            {dT2 / s1, dT2 * sigmaIWp, (-2 * eB * (-1 + w)) / sqrt3, (-2 * eB * wm1sq) / sqrt3, (4 * eC) / sqrt3,
             (2 * eB) / sqrt3},
            {dT3 / s1, dT3 * sigmaIWp, (-2 * eB * (-1 + w)) / sqrt3, (-2 * eB * wm1sq) / sqrt3, 0, (2 * eB) / sqrt3},
            {dT4 / s1, dT4 * sigmaIWp, (4 * eB * w) / sqrt3, (4 * eB * (-1 + w) * w) / sqrt3, 0, (-4 * eB) / sqrt3}};

        for (size_t n = 0; n < 12; ++n) {
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

    std::unique_ptr<FormFactorBase> FFLbtoLcstar12LSPRVar::clone(const std::string& label) {
        MAKE_CLONE(FFLbtoLcstar12LSPRVar, label);
    }

} // namespace Hammer
