///
/// @file  FFBtoDstarBLPRXPVar.cc
/// @brief \f$ B \rightarrow D^* \f$ BLPRXPVar form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/FormFactors/BLPRXP/FFBtoDstarBLPRXPVar.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFBtoDstarBLPRXPVar::FFBtoDstarBLPRXPVar() {
        // Create tensor rank and dimensions
        IndexList dims = {8, 11}; // size is _mFFErrNames + 1 to include central values in zeroth component
        string name{"FFBtoDstarBLPRXPVar"};

        setPrefix("BtoD*");
        _mFFErrLabel = FF_BDSTAR_VAR;
        addProcessSignature(PID::BPLUS, {-PID::DSTAR});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BDSTAR, FF_BDSTAR_VAR})});

        addProcessSignature(PID::BZERO, {PID::DSTARMINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BDSTAR, FF_BDSTAR_VAR})});

        setPrefix("BstoDs*");
        _mFFErrLabel = FF_BSDSSTAR_VAR;
        addProcessSignature(PID::BS, {PID::DSSTARMINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BSDSSTAR, FF_BSDSSTAR_VAR})});

        setSignatureIndex();
    }

    void FFBtoDstarBLPRXPVar::defineSettings() {
        //_mFFErrNames = ;
        setPath(getFFErrPrefixGroup().get());
        setUnits("GeV");

        FFBLPRXPVarBase::defineSettings(); // common parameters for BLPRXPVar are defined in FFBLPRXPVarBase

        addSetting<bool>("WithHA1As2", true);

        setInitialized();
    }

    void FFBtoDstarBLPRXPVar::evalAtPSPoint(const vector<double>& point, const vector<double>& masses) {
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
        const double Mb3 = Mb * Mb * Mb;

        double Sqq = point[0];

        // const double sqSqq = sqrt(Sqq);
        const double sqMbMc = sqrt(Mb * Mc);
        double w = getW(Sqq, Mb, Mc);

        // safety measure if w==1.0
        if (isZero(w - 1.0)) {
            w += 1e-6;
        }

        // recalc BLPRXPVar derived expansion parameters if settings changed
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

        const bool b1Mb2 = *getSetting<bool>("With1OverMb2");
        const double eb2 = b1Mb2 ? _eb * _eb : 0.;
        ;
        if (b1Mb2) {
            Hps += eb2 * (_mLi2[1](w) - _mLi2[4](w));
            Hv += eb2 * (_mLi2[1](w) - _mLi2[4](w));
            Ha1 += eb2 * (_mLi2[1](w) - _mLi2[4](w) * (w - 1) / (w + 1));
            Ha2 += 0.;
            Ha3 += eb2 * (_mLi2[1](w) - _mLi2[4](w));
            Ht1 += eb2 * (_mLi2[1](w));
            Ht2 += -eb2 * (_mLi2[4](w));
            Ht3 += 0.;
        }
        const bool b1MbMc = *getSetting<bool>("With1OverMbMc");
        const double eceb = b1MbMc ? _eb * _ec : 0.;
        if (b1MbMc) {
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
        const bool bAsM = *getSetting<bool>("WithAsOverM");
        const double asec = as * _ec;
        const double aseb = as * _eb;
        const double Cv2 = bAsM ? _asCorrections.CV2(w, _zBC) : 0.;
        const double Cv3 = bAsM ? _asCorrections.CV3(w, _zBC) : 0.;
        if (bAsM) {
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
        result.element({0, 0}) = -((Hps * Mc) / sqMbMc);
        // Ff
        result.element({1, 0}) = (Ha1 * ((Mb + Mc) * (Mb + Mc) - Sqq)) / (2. * sqMbMc);
        // Fg
        result.element({2, 0}) = Hv / (2. * sqMbMc);
        // Fm
        result.element({3, 0}) = -(Ha2 * sqrt(Mc / Mb3)) / 2. + Ha3 / (2. * sqMbMc);
        // Fp
        result.element({4, 0}) = -(Ha2 * sqrt(Mc / Mb3)) / 2. - Ha3 / (2. * sqMbMc);
        // Fzt
        result.element({5, 0}) = (Ht3 * Mc) / (2. * pow(Mb * Mc, 1.5));
        // Fmt
        result.element({6, 0}) = (Ht1 * (Mb - Mc)) / (2. * sqMbMc) - (Ht2 * (Mb + Mc)) / (2. * sqMbMc);
        // Fpt
        result.element({7, 0}) = (Ht2 * (Mb - Mc)) / (2. * sqMbMc) - (Ht1 * (Mb + Mc)) / (2. * sqMbMc);

        // Now create remaining FF tensor entries
        // n indexes rows ie FFs; idx indexes the columns that contract with delta.
        const double rC = Mc / Mb;
        const double wSq = w * w;
        const double sqrC = sqrt(rC);
        const double kHp = -sqrC;
        const double kHa1 = Mb * (w + 1) * sqrC;
        const double kHv = 1. / (2. * Mb * sqrC);
        const double kHa2 = -sqrC / (2. * Mb);
        const double kHa3 = 1. / (2. * Mb * sqrC);
        const double kHt3 = 1. / (2. * Mb2 * sqrC);
        const double kHt16 = (1 - rC) / (2. * sqrC);
        const double kHt26 = -(1 + rC) / (2. * sqrC);

        const double wm1Sq = (w - 1) * (w - 1);
        const double wm1Owp1 = (w - 1) / (w + 1);
        const double wm1SqOwp1 = wm1Sq / (w + 1);
        const double wSqm1 = wSq - 1;

        const vector<vector<double>>& sliwmat = (*getSetting<vector<vector<double>>>("sliwmatrix"));

        // Linearization of hatted FFs
        array<array<double, 10>, 8> FFmat{};

        // derivatives of Hp, Ha1, Hv, Ha2, Ha3, Ht3, Ht1, Ht2
        array<array<double, 10>, 8> DerMat1 = {
            {{0, 0, 4 * (-_eb + _ec) * (-1 + w), -4 * (_eb - _ec) * wm1Sq, 4 * (3 * _eb - _ec) * (-1 + w),
              -2 * _eb + 2 * _ec, -2 * (_eb - _ec) * (-1 + w), -2 * (eb2 + ec2 - eceb) * wSqm1,
              4 * (-eb2 + ec2) * (-1 + w), 4 * (3 * eb2 - ec2) * (-1 + w)},
             {0, 0, -4 * _eb * (-1 + w), -4 * _eb * wm1Sq, 4 * (3 * _eb - _ec) * (-1 + w), -2 * (_eb + eceb) * wm1Owp1,
              -2 * (_eb + eceb) * wm1SqOwp1, 2 * (-1 + w) * (eb2 + ec2 + eceb - (eb2 + ec2) * w), -4 * eb2 * (-1 + w),
              4 * (3 * eb2 - ec2) * (-1 + w)},
             {0, 0, -4 * _eb * (-1 + w), -4 * _eb * wm1Sq, 4 * (3 * _eb - _ec) * (-1 + w),
              -2 * (_eb + eceb) + (4 * eceb) / (1 + w), (-1 + w) * (-2 * (_eb + eceb) + (4 * eceb) / (1 + w)),
              -2 * (-1 + w) * (eb2 + ec2 - eceb + (eb2 + ec2) * w), -4 * eb2 * (-1 + w),
              4 * (3 * eb2 - ec2) * (-1 + w)},
             {0, 0, 4 * _ec, 4 * _ec * (-1 + w), 0, (-2 * (_ec - eceb)) / (1 + w), -2 * (_ec - eceb) * wm1Owp1,
              4 * ec2 * (-1 + w) - 2 * eceb * (2 + w), 4 * ec2, 0},
             {0, 0, -4 * (_ec + _eb * (-1 + w)), -4 * (_ec + _eb * (-1 + w)) * (-1 + w), 4 * (3 * _eb - _ec) * (-1 + w),
              (-2 * (_eb + _ec + (_eb + eceb) * w)) / (1 + w), -2 * (_eb + _ec + (_eb + eceb) * w) * wm1Owp1,
              2 * (eb2 - ec2 + eceb + 2 * (ec2 + eceb) * w - (eb2 + ec2) * wSq), -4 * (ec2 + eb2 * (-1 + w)),
              4 * (3 * eb2 - ec2) * (-1 + w)},
             {0, 0, -4 * _ec, -4 * _ec * (-1 + w), 0, (-2 * (_ec - eceb)) / (1 + w), -2 * (_ec - eceb) * wm1Owp1,
              4 * ec2 * (-1 + w) - 2 * eceb * (2 + w), -4 * ec2, 0},
             {0, 0, -4 * _eb * (-1 + w), -4 * _eb * wm1Sq, 4 * (3 * _eb - _ec) * (-1 + w), 2 * eceb * wm1Owp1,
              2 * eceb * wm1SqOwp1, -2 * eceb * (-1 + w), -4 * eb2 * (-1 + w), 4 * (3 * eb2 - ec2) * (-1 + w)},
             {0, 0, 0, 0, 0, -2 * _eb, -2 * _eb * (-1 + w), -2 * (eb2 - ec2) * wSqm1, 0, 0}}};

        array<array<double, 10>, 8> DerMat2{};
        if (bAsM) {
            DerMat2 = {
                {{0, 0, -4 * ((_cmagb + Cps) * _eb - (_cmagc + Cps) * _ec) * (-1 + w),
                  -4 * ((_cmagb + Cps) * _eb - (_cmagc + Cps) * _ec) * wm1Sq,
                  4 * (3 * (_cmagb + Cps) * _eb - (_cmagc + Cps) * _ec) * (-1 + w), 2 * Cps * (-_eb + _ec),
                  -2 * Cps * (_eb - _ec) * (-1 + w), 0, 0, 0},
                 {0, 0, -4 * (Ca1 + _cmagb) * _eb * (-1 + w), -4 * (Ca1 + _cmagb) * _eb * wm1Sq,
                  4 * (3 * (Ca1 + _cmagb) * _eb - (Ca1 + _cmagc) * _ec) * (-1 + w),
                  (-2 * (Ca1 * _eb - Ca2 * _eb - Ca3 * _ec) * (-1 + w)) / (1 + w),
                  (-2 * (Ca1 * _eb - Ca2 * _eb - Ca3 * _ec) * wm1Sq) / (1 + w), 0, 0, 0},
                 {0, 0, -4 * (_cmagb + Cv1) * _eb * (-1 + w), -4 * (_cmagb + Cv1) * _eb * wm1Sq,
                  4 * (3 * (_cmagb + Cv1) * _eb - (_cmagc + Cv1) * _ec) * (-1 + w),
                  -2 * ((Cv1 + Cv2) * _eb + Cv3 * _ec), -2 * ((Cv1 + Cv2) * _eb + Cv3 * _ec) * (-1 + w), 0, 0, 0},
                 {0, 0, 4 * (Ca1 + _cmagc) * _ec - 4 * Ca2 * (_eb - _ec) * (-1 + w),
                  -4 * (-((Ca1 + _cmagc) * _ec) + Ca2 * (_eb - _ec) * (-1 + w)) * (-1 + w),
                  4 * Ca2 * (3 * _eb - _ec) * (-1 + w),
                  (-2 * ((Ca1 + Ca3) * _ec - Ca2 * _ec * (1 + w) + Ca2 * _eb * (2 + w))) / (1 + w),
                  (-2 * (-1 + w) * ((Ca1 + Ca3) * _ec - Ca2 * _ec * (1 + w) + Ca2 * _eb * (2 + w))) / (1 + w), 0, 0, 0},
                 {0, 0,
                  -4 * (_cmagc * _ec + Ca1 * (_ec + _eb * (-1 + w)) + _cmagb * _eb * (-1 + w) +
                        Ca3 * (_eb - _ec) * (-1 + w)),
                  -4 *
                      (_cmagc * _ec + Ca1 * (_ec + _eb * (-1 + w)) + _cmagb * _eb * (-1 + w) +
                       Ca3 * (_eb - _ec) * (-1 + w)) *
                      (-1 + w),
                  4 * (3 * (Ca1 + Ca3 + _cmagb) * _eb - (Ca1 + Ca3 + _cmagc) * _ec) * (-1 + w),
                  (-2 * (-(Ca2 * _eb * w) + Ca1 * (_eb + _ec + _eb * w) + Ca3 * (_eb - _ec + _eb * w - 2 * _ec * w))) /
                      (1 + w),
                  (-2 * (-1 + w) *
                   (-(Ca2 * _eb * w) + Ca1 * (_eb + _ec + _eb * w) + Ca3 * (_eb - _ec + _eb * w - 2 * _ec * w))) /
                      (1 + w),
                  0, 0, 0},
                 {0, 0, -4 * ((_cmagc + Ct1) * _ec + Ct2 * _eb * (-1 + w)),
                  -4 * ((_cmagc + Ct1) * _ec + Ct2 * _eb * (-1 + w)) * (-1 + w), 4 * Ct2 * (3 * _eb - _ec) * (-1 + w),
                  (-2 * ((Ct1 + Ct3) * _ec + Ct2 * _eb * (2 + w))) / (1 + w),
                  (-2 * (-1 + w) * ((Ct1 + Ct3) * _ec + Ct2 * _eb * (2 + w))) / (1 + w), 0, 0, 0},
                 {0, 0, -2 * _eb * (2 * _cmagb + 2 * Ct1 + (Ct2 - Ct3) * (-1 + w)) * (-1 + w),
                  -2 * _eb * (2 * _cmagb + 2 * Ct1 + (Ct2 - Ct3) * (-1 + w)) * wm1Sq,
                  2 * (6 * _cmagb * _eb - 2 * _cmagc * _ec + (3 * _eb - _ec) * (2 * Ct1 + (Ct2 - Ct3) * (-1 + w))) *
                      (-1 + w),
                  -((Ct2 - Ct3) * _eb * (-1 + w)), -((Ct2 - Ct3) * _eb * wm1Sq), 0, 0, 0},
                 {0, 0, -2 * (Ct2 + Ct3) * _eb * wSqm1, -2 * (Ct2 + Ct3) * _eb * (1 + w) * wm1Sq,
                  2 * (Ct2 + Ct3) * (3 * _eb - _ec) * wSqm1,
                  -2 * Ct1 * _eb + Ct2 * _eb - Ct3 * _eb + 2 * Ct3 * _ec - (Ct2 + Ct3) * _eb * w,
                  -((-1 + w) * (2 * Ct1 * _eb + Ct2 * _eb * (-1 + w) + Ct3 * (_eb - 2 * _ec + _eb * w))), 0, 0, 0}}};
        }

        for (size_t i = 0; i < 10; ++i) {
            FFmat[0][i] = kHp * (DerMat1[0][i] + as * DerMat2[0][i]);
            FFmat[1][i] = kHa1 * (DerMat1[1][i] + as * DerMat2[1][i]);
            FFmat[2][i] = kHv * (DerMat1[2][i] + as * DerMat2[2][i]);
            FFmat[3][i] = kHa2 * (DerMat1[3][i] + as * DerMat2[3][i]) + kHa3 * (DerMat1[4][i] + as * DerMat2[4][i]);
            FFmat[4][i] = kHa2 * (DerMat1[3][i] + as * DerMat2[3][i]) - kHa3 * (DerMat1[4][i] + as * DerMat2[4][i]);
            FFmat[5][i] = kHt3 * (DerMat1[5][i] + as * DerMat2[5][i]);
            FFmat[6][i] = kHt16 * (DerMat1[6][i] + as * DerMat2[6][i]) + kHt26 * (DerMat1[7][i] + as * DerMat2[7][i]);
            FFmat[7][i] = kHt26 * (DerMat1[6][i] + as * DerMat2[6][i]) + kHt16 * (DerMat1[7][i] + as * DerMat2[7][i]);
        }

        // Linearization of LO IW functions
        const array<double, 10> xiIWL = XiLinearCoeffs(w);

        for (size_t n = 0; n < 8; ++n) {
            auto ni = static_cast<IndexType>(n);
            for (size_t idx = 0; idx < 10; ++idx) {
                complex<double> entry = 0;
                auto idxi = static_cast<IndexType>(idx + 1u);
                for (size_t idx1 = 0; idx1 < 10; ++idx1) {
                    entry += sliwmat[idx1][idx] * (FFmat[n][idx1] + result.element({ni, 0}) * xiIWL[idx1]);
                }
                result.element({ni, idxi}) = entry;
            }
        }

        result *= h;
        result.toVector();
    }

    std::unique_ptr<FormFactorBase> FFBtoDstarBLPRXPVar::clone(const std::string& label) {
        MAKE_CLONE(FFBtoDstarBLPRXPVar, label);
    }

} // namespace Hammer
