///
/// @file  FFLbtoLcBLRS.cc
/// @brief \f$ \Lambda_b \rightarrow Lambda_c \f$ BLRS form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/FormFactors/BLRS/FFLbtoLcBLRS.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFLbtoLcBLRS::FFLbtoLcBLRS() {
        // Create tensor rank and dimensions
        IndexList dims = {12};
        string name{"FFLbtoLcBLRS"};

        setPrefix("LbtoLc");
        addProcessSignature(-PID::LAMBDAB, {PID::LAMBDACMINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_LBLC})});

        setSignatureIndex();
    }

    void FFLbtoLcBLRS::defineSettings() {
        //_mFFErrNames = ;
        setPath(getFFErrPrefixGroup().get());
        setUnits("GeV");

        FFBLRSBase::defineSettings(); // common parameters for BLRS are defined in FFBLRSBase

        setInitialized();
    }

    void FFLbtoLcBLRS::eval(const Particle& parent, const ParticleList& daughters, const ParticleList& /*references*/) {
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

    void FFLbtoLcBLRS::evalAtPSPoint(const vector<double>& point, const vector<double>& masses) {
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
        // const double wSq = w*w;

        calcConstants(unitres);
        fillKisMis(unitres);
        // const double mc2 = pow(_mc1S,2.);

        //(Sub)leading IW function parameters
        const double z1 = (*getSetting<double>("z1"));
        const double z2 = (*getSetting<double>("z2"));

        // Leading IW function
        const double zetaIW = 1. + z1 * (w - 1) + 0.5 * z2 * pow(w - 1, 2.);

        // QCD correction functions
        const double Cs = _asCorrections.CS(w, _zBC);
        const double Cps = _asCorrections.CP(w, _zBC);
        const double Cv1 = _asCorrections.CV1(w, _zBC);
        const double Cv2 = _asCorrections.CV2(w, _zBC);
        const double Cv3 = _asCorrections.CV3(w, _zBC);
        const double Ca1 = _asCorrections.CA1(w, _zBC);
        const double Ca2 = _asCorrections.CA2(w, _zBC);
        const double Ca3 = _asCorrections.CA3(w, _zBC);
        const double Ct1 = _asCorrections.CT1(w, _zBC);
        const double Ct2 = _asCorrections.CT2(w, _zBC);
        const double Ct3 = _asCorrections.CT3(w, _zBC);
        // Derivatives (approx)
        const double Csp = _asCorrections.derCS(w, _zBC);
        const double Cpsp = _asCorrections.derCP(w, _zBC);
        const double Cv1p = _asCorrections.derCV1(w, _zBC);
        const double Cv2p = _asCorrections.derCV2(w, _zBC);
        const double Cv3p = _asCorrections.derCV3(w, _zBC);
        const double Ca1p = _asCorrections.derCA1(w, _zBC);
        const double Ca2p = _asCorrections.derCA2(w, _zBC);
        const double Ca3p = _asCorrections.derCA3(w, _zBC);
        const double Ct1p = _asCorrections.derCT1(w, _zBC);
        const double Ct2p = _asCorrections.derCT2(w, _zBC);
        const double Ct3p = _asCorrections.derCT3(w, _zBC);

        // Create hatted h FFs
        // LO + NLO
        double hs = 1. + _aS * Cs + (_eB + _eC) * _mKi1[1](w) - (_eB + _eC) * _mKi1[2](w) * (-1. + w);
        double hp = 1. + _aS * Cps + (_eB + _eC) * _mKi1[1](w) - (_eB + _eC) * _mKi1[2](w) * (1. + w);
        double f1 = 1. + _aS * Cv1 + (_eB + _eC) * _mKi1[1](w) - (_eB + _eC) * _mKi1[2](w) * (1. + w);
        double f2 = _aS * Cv2 + 2. * _eC * _mKi1[2](w);
        double f3 = _aS * Cv3 + 2. * _eB * _mKi1[2](w);
        double g1 = 1. + _aS * Ca1 + (_eB + _eC) * _mKi1[1](w) - (_eB + _eC) * _mKi1[2](w) * (-1. + w);
        double g2 = _aS * Ca2 + 2. * _eC * _mKi1[2](w);
        double g3 = _aS * Ca3 - 2. * _eB * _mKi1[2](w);
        double h1 = 1. + _aS * Ct1 + (_eB + _eC) * _mKi1[1](w) - (_eB + _eC) * _mKi1[2](w) * (-1. + w);
        double h2 = _aS * Ct2 + 2. * _eC * _mKi1[2](w);
        double h3 = _aS * Ct3 - 2. * _eB * _mKi1[2](w);
        double h4 = 0.;

        // as x 1/m
        hs += _aS * (Cs * ((_eB0 + _eC0) * _mKi1[1](w) - (_eB0 + _eC0) * _mKi1[2](w) * (w - 1)) +
                     2. * (_eB0 + _eC0) * (w - 1) * Csp);
        hp += _aS * (Cps * ((_eB0 + _eC0) * _mKi1[1](w) - (_eB0 + _eC0) * _mKi1[2](w) * (w + 1)) +
                     2. * (_eB0 + _eC0) * (w - 1) * Cpsp);
        f1 += _aS * (Cv1 * ((_eB0 + _eC0) * _mKi1[1](w) - (_eB0 + _eC0) * _mKi1[2](w) * (w + 1)) +
                     2. * (_eB0 + _eC0) * (w - 1) * Cv1p);
        f2 += _aS *
              (2. * Cv1 * _eC0 * _mKi1[2](w) - (2. * Cv3 * _eC0) / (w + 1) +
               Cv2 * ((_eB0 + _eC0) * _mKi1[1](w) - (_eB0 + _eC0) * _mKi1[2](w) * (w - 1) + (2. * _eB0 * w) / (w + 1)) +
               2. * (_eB0 + _eC0) * (w - 1) * Cv2p);
        f3 += _aS *
              (2. * Cv1 * _eB0 * _mKi1[2](w) - (2. * Cv2 * _eB0) / (w + 1) +
               Cv3 * ((_eB0 + _eC0) * _mKi1[1](w) - (_eB0 + _eC0) * _mKi1[2](w) * (w - 1) + (2. * _eC0 * w) / (w + 1)) +
               2. * (_eB0 + _eC0) * (w - 1) * Cv3p);
        g1 += _aS * (Ca1 * ((_eB0 + _eC0) * _mKi1[1](w) - (_eB0 + _eC0) * _mKi1[2](w) * (w - 1)) +
                     2. * (_eB0 + _eC0) * (w - 1) * Ca1p);
        g2 += _aS *
              (2. * Ca1 * _eC0 * _mKi1[2](w) - (2. * Ca3 * _eC0) / (w + 1) +
               Ca2 * ((_eB0 + _eC0) * _mKi1[1](w) + (2. * _eB0 * w) / (w + 1) - (_eB0 + _eC0) * _mKi1[2](w) * (w + 1)) +
               2. * (_eB0 + _eC0) * (w - 1) * Ca2p);
        g3 += _aS *
              (-2. * Ca1 * _eB0 * _mKi1[2](w) - (2. * Ca2 * _eB0) / (w + 1) +
               Ca3 * ((_eB0 + _eC0) * _mKi1[1](w) + (2. * _eC0 * w) / (w + 1) - (_eB0 + _eC0) * _mKi1[2](w) * (w + 1)) +
               2. * (_eB0 + _eC0) * (w - 1) * Ca3p);
        h1 += _aS * (Ct1 * ((_eB0 + _eC0) * _mKi1[1](w) - (_eB0 + _eC0) * _mKi1[2](w) * (w - 1)) +
                     2. * (_eB0 + _eC0) * (w - 1) * Ct1p);
        h2 += _aS *
              (2. * Ct1 * _eC0 * _mKi1[2](w) - (2. * Ct3 * _eC0) / (w + 1) +
               Ct2 * ((_eB0 + _eC0) * _mKi1[1](w) + (2. * _eB0 * w) / (w + 1) - (_eB0 + _eC0) * _mKi1[2](w) * (w + 1)) +
               2. * (_eB0 + _eC0) * (w - 1) * Ct2p);
        h3 += _aS *
              (-2. * Ct1 * _eB0 * _mKi1[2](w) - (2. * Ct2 * _eB0) / (w + 1) +
               Ct3 * ((_eB0 + _eC0) * _mKi1[1](w) + (2. * _eC0 * w) / (w + 1) - (_eB0 + _eC0) * _mKi1[2](w) * (w + 1)) +
               2. * (_eB0 + _eC0) * (w - 1) * Ct3p);
        h4 += _aS * (2. * Ct2 * _eB0 * _mKi1[2](w) - 2. * Ct3 * _eC0 * _mKi1[2](w));

        // NNLO 1/mc^2 only
        const double eC02 = _eC0 * _eC0;
        hs += eC02 * _mKi2[1](w) - eC02 * _mKi2[2](w) * (-1. + w);
        hp += eC02 * _mKi2[1](w) - eC02 * _mKi2[2](w) * (1. + w);
        f1 += eC02 * _mKi2[1](w) - eC02 * _mKi2[2](w) * (1. + w);
        f2 += 2. * eC02 * _mKi2[2](w);
        f3 += 0.;
        g1 += eC02 * _mKi2[1](w) - eC02 * _mKi2[2](w) * (-1. + w);
        g2 += 2. * eC02 * _mKi2[2](w);
        g3 += 0.;
        h1 += eC02 * _mKi2[1](w) - eC02 * _mKi2[2](w) * (-1. + w);
        h2 += 2. * eC02 * _mKi2[2](w);
        h3 += 0.;
        h4 += 0.;

        // set elements
        result.element({0}) = hs;
        result.element({1}) = hp;
        result.element({2}) = f1;
        result.element({3}) = f2;
        result.element({4}) = f3;
        result.element({5}) = g1;
        result.element({6}) = g2;
        result.element({7}) = g3;
        result.element({8}) = h1;
        result.element({9}) = h2;
        result.element({10}) = h3;
        result.element({11}) = h4;

        result *= zetaIW;
    }

    std::unique_ptr<FormFactorBase> FFLbtoLcBLRS::clone(const std::string& label) {
        MAKE_CLONE(FFLbtoLcBLRS, label);
    }

} // namespace Hammer
