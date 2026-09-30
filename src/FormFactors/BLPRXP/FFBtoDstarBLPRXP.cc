///
/// @file  FFBtoDstarBLPRXP.cc
/// @brief \f$ B \rightarrow D^* \f$ BLPRXP form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/FormFactors/BLPRXP/FFBtoDstarBLPRXP.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFBtoDstarBLPRXP::FFBtoDstarBLPRXP() {
        // Create tensor rank and dimensions
        IndexList dims = {8};
        string name{"FFBtoDstarBLPRXP"};

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

    void FFBtoDstarBLPRXP::defineSettings() {
        //_mFFErrNames = ;
        setPath(getFFErrPrefixGroup().get());
        setUnits("GeV");

        FFBLPRXPBase::defineSettings(); // common parameters for BLPRXP are defined in FFBLPRXPBase

        addSetting<bool>("WithHA1As2", true);

        setInitialized();
    }

    void FFBtoDstarBLPRXP::evalAtPSPoint(const vector<double>& point, const vector<double>& masses) {
        Tensor& result = getTensor();
        result.clearData();

        if (!isInitialized()) {
            MSG_WARNING("Warning, Settings have not been defined!");
        }

        double Mb = 0.;
        double Mc = 0.;
        double unitres = 1.;
        tie(Mb, Mc, unitres) = getParentDaughterHadMasses(masses);

        const double Mb3 = Mb * Mb * Mb;

        double Sqq = point[0];

        // const double sqSqq = sqrt(Sqq);
        const double sqMbMc = sqrt(Mb * Mc);
        double w = getW(Sqq, Mb, Mc);

        // safety measure if w==1.0
        if (isZero(w - 1.0)) {
            w += 1e-6;
        }

        // recalc BLPRXP derived expansion parameters if settings changed
        calcConstants();
        const double as = (*getSetting<double>("as")) / pi;

        // LO IW function
        const double h = Xi(w);

        // QCD correction functions
        //  const double Cs  = CS(w, zBC);
        const double Cps = _asCorrections.CP(w, _zBC);
        const double Cv1 = _asCorrections.CV1(w, _zBC);
        const double Ca1 = _asCorrections.CA1(w, _zBC);
        const double Ca2 = _asCorrections.CA2(w, _zBC);
        const double Ca3 = _asCorrections.CA3(w, _zBC);
        const double Ct1 = _asCorrections.CT1(w, _zBC);
        const double Ct2 = _asCorrections.CT2(w, _zBC);
        const double Ct3 = _asCorrections.CT3(w, _zBC);

        // Create hatted h FFs
        // LO
        double Hps = 1.;
        double Hv = 1.;
        double Ha1 = 1.;
        double Ha2 = 0.;
        double Ha3 = 1.;
        double Ht1 = 1.;
        double Ht2 = 0.;
        double Ht3 = 0.;

        // as
        Hps += as * Cps;
        Hv += as * Cv1;
        Ha1 += as * Ca1;
        Ha2 += as * Ca2;
        Ha3 += as * (Ca1 + Ca3);
        Ht1 += as * (Ct1 + 0.5 * (w - 1) * (Ct2 - Ct3));
        Ht2 += as * 0.5 * (w + 1) * (Ct2 + Ct3);
        Ht3 += as * Ct2;

        // 1/m
        Hps += _ec * (_mLi1[2](w) + _mLi1[3](w) * (w - 1) + _mLi1[5](w) - _mLi1[6](w) * (w + 1)) +
               _eb * (_mLi1[1](w) - _mLi1[4](w));
        Hv += _ec * (_mLi1[2](w) - _mLi1[5](w)) + _eb * (_mLi1[1](w) - _mLi1[4](w));
        Ha1 += _ec * (_mLi1[2](w) - _mLi1[5](w) * (w - 1) / (w + 1)) +
               _eb * (_mLi1[1](w) - _mLi1[4](w) * (w - 1) / (w + 1));
        Ha2 += _ec * (_mLi1[3](w) + _mLi1[6](w));
        Ha3 += _ec * (_mLi1[2](w) - _mLi1[3](w) - _mLi1[5](w) + _mLi1[6](w)) + _eb * (_mLi1[1](w) - _mLi1[4](w));
        Ht1 += _ec * (_mLi1[2](w)) + _eb * (_mLi1[1](w));
        Ht2 += _ec * (_mLi1[5](w)) - _eb * (_mLi1[4](w));
        Ht3 += _ec * (_mLi1[6](w) - _mLi1[3](w));


        // Upsilon expansion
        const double corrb = _eb * _upsilonb;
        const double corrc = _ec * _upsilonc;
        Hps += (corrc + corrb);
        Hv += (corrc + corrb);
        Ha1 += (corrc + corrb) * (w - 1) / (w + 1);
        Ha2 += -2. * corrc / (w + 1.);
        Ha3 += (corrc * (w - 1.) / (w + 1.) + corrb);
        Ht1 += 0.;
        Ht2 += -(corrc - corrb);
        Ht3 += -2. * corrc / (w + 1.);

        // 1/m^2
        const double ec2 = _ec * _ec;
        Hps += ec2 * (_mLi2[2](w) + _mLi2[3](w) * (w - 1) + _mLi2[5](w) - _mLi2[6](w) * (w + 1));
        Hv += ec2 * (_mLi2[2](w) - _mLi2[5](w));
        Ha1 += ec2 * (_mLi2[2](w) - _mLi2[5](w) * (w - 1) / (w + 1));
        Ha2 += ec2 * (_mLi2[3](w) + _mLi2[6](w));
        Ha3 += ec2 * (_mLi2[2](w) - _mLi2[3](w) - _mLi2[5](w) + _mLi2[6](w));
        Ht1 += ec2 * (_mLi2[2](w));
        Ht2 += ec2 * (_mLi2[5](w));
        Ht3 += ec2 * (_mLi2[6](w) - _mLi2[3](w));

        if (*getSetting<bool>("With1OverMb2")) {
            const double eb2 = _eb * _eb;
            Hps += eb2 * (_mLi2[1](w) - _mLi2[4](w));
            Hv += eb2 * (_mLi2[1](w) - _mLi2[4](w));
            Ha1 += eb2 * (_mLi2[1](w) - _mLi2[4](w) * (w - 1) / (w + 1));
            Ha2 += 0.;
            Ha3 += eb2 * (_mLi2[1](w) - _mLi2[4](w));
            Ht1 += eb2 * (_mLi2[1](w));
            Ht2 += -eb2 * (_mLi2[4](w));
            Ht3 += 0.;
        }
        if (*getSetting<bool>("With1OverMbMc")) {
            const double eceb = _eb * _ec;
            Hv += eceb * ((_mMi[2](w) + _mMi[9](w)) - (_mMi[16](w) + _mMi[18](w)));
            Ha1 += eceb * ((_mMi[2](w) + _mMi[9](w)) - (w - 1.) / (w + 1.) * (_mMi[16](w) + _mMi[18](w)));
            Ha2 += eceb * ((_mMi[3](w) - _mMi[10](w)) + (_mMi[17](w) - _mMi[19](w)));
            Ha3 += eceb * ((_mMi[2](w) + _mMi[9](w)) - (_mMi[3](w) - _mMi[10](w)) - (_mMi[16](w) + _mMi[18](w)) +
                           (_mMi[17](w) - _mMi[19](w)));
            Hps += eceb * (_mMi[2](w) - _mMi[9](w) + (w - 1.) * (_mMi[3](w) + _mMi[10](w)) + _mMi[16](w) - _mMi[18](w) -
                           (1. + w) * (_mMi[17](w) + _mMi[19](w)));
            Ht1 += eceb * (_mMi[2](w) - _mMi[9](w));
            Ht2 += eceb * (_mMi[16](w) - _mMi[18](w));
            Ht3 += -eceb * ((_mMi[3](w) + _mMi[10](w)) - (_mMi[17](w) + _mMi[19](w)));
        }

        // as/m
        if (*getSetting<bool>("WithAsOverM")) {
            const double asec = as * _ec;
            const double aseb = as * _eb;
            const double Cv2 = _asCorrections.CV2(w, _zBC);
            const double Cv3 = _asCorrections.CV3(w, _zBC);
            const double dCps = _asCorrections.derCP(w, _zBC);
            const double dCv1 = _asCorrections.derCV1(w, _zBC);
            const double dCa1 = _asCorrections.derCA1(w, _zBC);
            const double dCa2 = _asCorrections.derCA2(w, _zBC);
            const double dCa3 = _asCorrections.derCA3(w, _zBC);
            const double dCt1 = _asCorrections.derCT1(w, _zBC);
            const double dCt2 = _asCorrections.derCT2(w, _zBC);
            const double dCt3 = _asCorrections.derCT3(w, _zBC);
            Hv += asec * (_cmagc * _mLi1[2](w) + (_mLi1[2](w) - _mLi1[5](w)) * Cv1 - (_mLi1[4](w) - _mLi1[5](w)) * Cv3 +
                          2 * (-1 + w) * dCv1);
            Ha1 += asec * (_cmagc * _mLi1[2](w) + (_mLi1[2](w) - (_mLi1[5](w) * (-1 + w)) / (1 + w)) * Ca1 +
                           ((_mLi1[4](w) - _mLi1[5](w)) * (-1 + w) * Ca3) / (1 + w) + 2 * (-1 + w) * dCa1);
            Ha2 += asec * (_cmagc * _mLi1[3](w) + (_mLi1[3](w) + _mLi1[6](w)) * Ca1 +
                           (_mLi1[2](w) + _mLi1[5](w) + _mLi1[3](w) * (-1 + w) - _mLi1[6](w) * (1 + w)) * Ca2 -
                           ((_mLi1[4](w) - 3 * _mLi1[5](w)) * Ca3) / (1 + w) + 2 * (-1 + w) * dCa2);
            Ha3 += asec * (_cmagc * (-_mLi1[3](w) + _mLi1[2](w)) +
                           (_mLi1[2](w) - _mLi1[3](w) - _mLi1[5](w) + _mLi1[6](w)) * Ca1 +
                           ((_mLi1[4](w) - 3 * _mLi1[5](w)) * w * Ca3) / (1 + w) +
                           (_mLi1[2](w) + _mLi1[5](w) + _mLi1[3](w) * (-1 + w) - _mLi1[6](w) * (1 + w)) * Ca3 +
                           2 * (-1 + w) * (dCa1 + dCa3));
            Hps += asec * ((_cmagc * (_mLi1[3](w) * (-1 + w) + _mLi1[2](w))) +
                           (_mLi1[2](w) + _mLi1[5](w) + _mLi1[3](w) * (-1 + w) - _mLi1[6](w) * (1 + w)) * Cps +
                           2 * (-1 + w) * dCps);
            Ht1 += asec * (_cmagc * _mLi1[2](w) + _mLi1[2](w) * (Ct1 + ((-1 + w) * (Ct2 - Ct3)) / 2.) -
                           (_mLi1[5](w) * (-1 + w) * (Ct2 - Ct3)) / 2. + _mLi1[5](w) * (-1 + w) * Ct3 +
                           2 * (-1 + w) * (dCt1 + ((-1 + w) * (dCt2 - dCt3)) / 2.));
            Ht2 += asec * ((_mLi1[4](w) - _mLi1[5](w) * w) * Ct3 + (_mLi1[2](w) * (1 + w) * (Ct2 + Ct3)) / 2. +
                           _mLi1[5](w) * (Ct1 - ((1 + w) * (Ct2 + Ct3)) / 2.) + (-1 + pow(w, 2)) * (dCt2 + dCt3));
            Ht3 += asec *
                   (-(_cmagc * _mLi1[3](w)) - (_mLi1[3](w) - _mLi1[6](w)) * Ct1 + (_mLi1[2](w) - _mLi1[5](w)) * Ct2 -
                    ((_mLi1[4](w) - 3 * _mLi1[5](w)) * Ct3) / (1 + w) + 2 * (-1 + w) * dCt2);
            Hv += aseb * (_cmagb * _mLi1[1](w) + (_mLi1[1](w) - _mLi1[4](w)) * Cv1 - (_mLi1[4](w) - _mLi1[5](w)) * Cv2 +
                          2 * (-1 + w) * dCv1);
            Ha1 += aseb * (_cmagb * _mLi1[1](w) + (_mLi1[1](w) - (_mLi1[4](w) * (-1 + w)) / (1 + w)) * Ca1 +
                           ((_mLi1[4](w) - _mLi1[5](w)) * (-1 + w) * Ca2) / (1 + w) + 2 * (-1 + w) * dCa1);
            Ha2 += aseb * ((_mLi1[1](w) - _mLi1[4](w) - 2 * _mLi1[5](w)) * Ca2 -
                           ((_mLi1[4](w) - 3 * _mLi1[5](w)) * Ca2) / (1 + w) + 2 * (-1 + w) * dCa2);
            Ha3 += aseb * (_cmagb * _mLi1[1](w) + 2 * _mLi1[5](w) * Ca2 +
                           ((_mLi1[4](w) - 3 * _mLi1[5](w)) * w * Ca2) / (1 + w) +
                           (_mLi1[1](w) - _mLi1[4](w)) * (Ca1 + Ca3) + 2 * (-1 + w) * (dCa1 + dCa3));
            Hps += aseb * (_cmagb * _mLi1[1](w) + (_mLi1[1](w) - _mLi1[4](w)) * Cps + 2 * (-1 + w) * dCps);
            Ht1 += aseb *
                   (_cmagb * _mLi1[1](w) - _mLi1[5](w) * (-1 + w) * Ct2 +
                    _mLi1[1](w) * (Ct1 + ((-1 + w) * (Ct2 - Ct3)) / 2.) - (_mLi1[4](w) * (-1 + w) * (Ct2 - Ct3)) / 2. +
                    2 * (-1 + w) * (dCt1 + ((-1 + w) * (dCt2 - dCt3)) / 2.));
            Ht2 += aseb * ((_mLi1[4](w) - _mLi1[5](w) * w) * Ct2 + (_mLi1[1](w) * (1 + w) * (Ct2 + Ct3)) / 2. -
                           _mLi1[4](w) * (Ct1 + ((1 + w) * (Ct2 + Ct3)) / 2.) + (-1 + pow(w, 2)) * (dCt2 + dCt3));
            Ht3 += aseb * ((_mLi1[1](w) - _mLi1[4](w) - 2 * _mLi1[5](w)) * Ct2 -
                           ((_mLi1[4](w) - 3 * _mLi1[5](w)) * Ct2) / (1 + w) + 2 * (-1 + w) * dCt2);
        }


        if (*getSetting<bool>("WithHA1As2")) {
            Ha1 += -0.944 * (4. / 3.) * as * as;
        }

        // set elements, mapping to amplitude FF basis
        // Fs
        result.element({0}) = -((Hps * Mc) / sqMbMc);
        // Ff
        result.element({1}) = (Ha1 * ((Mb + Mc) * (Mb + Mc) - Sqq)) / (2. * sqMbMc);
        // Fg
        result.element({2}) = Hv / (2. * sqMbMc);
        // Fm
        result.element({3}) = -(Ha2 * sqrt(Mc / Mb3)) / 2. + Ha3 / (2. * sqMbMc);
        // Fp
        result.element({4}) = -(Ha2 * sqrt(Mc / Mb3)) / 2. - Ha3 / (2. * sqMbMc);
        // Fzt
        result.element({5}) = (Ht3 * Mc) / (2. * pow(Mb * Mc, 1.5));
        // Fmt
        result.element({6}) = (Ht1 * (Mb - Mc)) / (2. * sqMbMc) - (Ht2 * (Mb + Mc)) / (2. * sqMbMc);
        // Fpt
        result.element({7}) = (Ht2 * (Mb - Mc)) / (2. * sqMbMc) - (Ht1 * (Mb + Mc)) / (2. * sqMbMc);

        result *= h;
    }

    std::unique_ptr<FormFactorBase> FFBtoDstarBLPRXP::clone(const std::string& label) {
        MAKE_CLONE(FFBtoDstarBLPRXP, label);
    }

} // namespace Hammer
