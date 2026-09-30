///
/// @file  FFBtoDBLPRXPVar.cc
/// @brief \f$ B \rightarrow D \f$ BLPRXPVar form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>
#include <array>

#include "Hammer/FormFactors/BLPRXP/FFBtoDBLPRXPVar.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFBtoDBLPRXPVar::FFBtoDBLPRXPVar() {
        // Create tensor rank and dimensions
        IndexList dims = {4, 11}; // size is _mFFErrNames + 1 to include central values in zeroth component
        string name{"FFBtoDBLPRXPVar"};

        setPrefix("BtoD");
        _mFFErrLabel = FF_BD_VAR;
        addProcessSignature(PID::BPLUS, {-PID::D0});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BD, FF_BD_VAR})});

        addProcessSignature(PID::BZERO, {PID::DMINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BD, FF_BD_VAR})});

        setPrefix("BstoDs");
        _mFFErrLabel = FF_BSDS_VAR;
        addProcessSignature(PID::BS, {PID::DSMINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BSDS, FF_BSDS_VAR})});

        setSignatureIndex();
    }

    void FFBtoDBLPRXPVar::defineSettings() {
        //_mFFErrNames = ;
        setPath(getFFErrPrefixGroup().get());
        setUnits("GeV");

        FFBLPRXPVarBase::defineSettings(); // common parameters for BLPRXPVar are defined in FFBLPRXPVarBase

        setInitialized();
    }

    void FFBtoDBLPRXPVar::evalAtPSPoint(const vector<double>& point, const vector<double>& masses) {
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
        const double sqMbMc = sqrt(Mb * Mc);
        double w = getW(Sqq, Mb, Mc);

        // recalc BLPRXPVar derived expansion parameters if settings changed
        calcConstants();
        const double as = (*getSetting<double>("as")) / pi;

        // LO IW function
        const double h = Xi(w);

        // QCD correction functions
        const double Cs = _asCorrections.CS(w, _zBC);
        // const double Cps = CP(w, zBC);
        const double Cv1 = _asCorrections.CV1(w, _zBC);
        const double Cv2 = _asCorrections.CV2(w, _zBC);
        const double Cv3 = _asCorrections.CV3(w, _zBC);
        // const double Ca1 = CA1(w, zBC);
        // const double Ca2 = CA2(w, zBC);
        // const double Ca3 = CA3(w, zBC);
        const double Ct1 = _asCorrections.CT1(w, _zBC);
        const double Ct2 = _asCorrections.CT2(w, _zBC);
        const double Ct3 = _asCorrections.CT3(w, _zBC);


        // Create hatted h FFs
        // LO
        double Hs = 1.;
        double Hp = 1.;
        double Hm = 0.;
        double Ht = 1.;

        // as
        Hs += as * Cs;
        Hp += as * (Cv1 + 0.5 * (w + 1) * (Cv2 + Cv3));
        Hm += as * 0.5 * (w + 1) * (Cv2 - Cv3);
        Ht += as * (Ct1 - Ct2 + Ct3);

        // 1/m
        Hs += (_ec + _eb) * (_mLi1[1](w) - _mLi1[4](w) * (w - 1) / (w + 1));
        Hp += (_ec + _eb) * _mLi1[1](w);
        Hm += (_ec - _eb) * _mLi1[4](w);
        Ht += (_ec + _eb) * (_mLi1[1](w) - _mLi1[4](w));

        // Upsilon expansion
        const double corrb = _eb * _upsilonb;
        const double corrc = _ec * _upsilonc;
        Hs += (corrc + corrb) * (w - 1) / (w + 1);
        Hp += 0.;
        Hm += -(corrc - corrb);
        Ht += (corrc + corrb);

        // 1/m^2
        const double ec2 = _ec * _ec;
        Hs += ec2 * (_mLi2[1](w) - _mLi2[4](w) * (w - 1) / (w + 1));
        Hp += ec2 * _mLi2[1](w);
        Hm += ec2 * _mLi2[4](w);
        Ht += ec2 * (_mLi2[1](w) - _mLi2[4](w));

        const bool b1Mb2 = *getSetting<bool>("With1OverMb2");
        const double eb2 = b1Mb2 ? _eb * _eb : 0.;
        ;
        if (b1Mb2) {
            Hs += +eb2 * (_mLi2[1](w) - _mLi2[4](w) * (w - 1) / (w + 1));
            Hp += +eb2 * _mLi2[1](w);
            Hm += -eb2 * _mLi2[4](w);
            Ht += +eb2 * (_mLi2[1](w) - _mLi2[4](w));
        }
        const bool b1MbMc = *getSetting<bool>("With1OverMbMc");
        const double eceb = b1MbMc ? _eb * _ec : 0.;
        if (b1MbMc) {
            Hs += eceb * (_mMi[1](w) + _mMi[8](w) - 2. * (w - 1.) / (w + 1.) * _mMi[15](w));
            Hp += eceb * (_mMi[1](w) - _mMi[8](w));
            Hm += 0.;
            Ht += eceb * (_mMi[1](w) + _mMi[8](w) - 2. * _mMi[15](w));
        }

        // as/m
        const bool bAsM = *getSetting<bool>("WithAsOverM");
        const double asec = as * _ec;
        const double aseb = as * _eb;
        if (bAsM) {
            const double dCs = _asCorrections.derCS(w, _zBC);
            const double dCv1 = _asCorrections.derCV1(w, _zBC);
            const double dCv2 = _asCorrections.derCV2(w, _zBC);
            const double dCv3 = _asCorrections.derCV3(w, _zBC);
            const double dCt1 = _asCorrections.derCT1(w, _zBC);
            const double dCt2 = _asCorrections.derCT2(w, _zBC);
            const double dCt3 = _asCorrections.derCT3(w, _zBC);
            Hs += asec *
                  (_cmagc * _mLi1[1](w) + (_mLi1[1](w) - (_mLi1[4](w) * (-1 + w)) / (1 + w)) * Cs + 2 * (-1 + w) * dCs);
            Hp += asec * (_cmagc * _mLi1[1](w) + _mLi1[1](w) * Cv1 - _mLi1[5](w) * (-1 + w) * Cv3 +
                          ((1 + w) * (_mLi1[1](w) - (_mLi1[4](w) * (-1 + w)) / (1 + w)) * (Cv2 + Cv3)) / 2. +
                          2 * (-1 + w) * (dCv1 + ((1 + w) * (dCv2 + dCv3)) / 2.));
            Hm += asec * (_mLi1[4](w) * Cv1 + _mLi1[5](w) * (1 + w) * Cv3 +
                          ((1 + w) * ((_mLi1[1](w) - (_mLi1[4](w) * (-1 + w)) / (1 + w)) * (Cv2 - Cv3) +
                                      2 * (-1 + w) * (dCv2 - dCv3))) /
                              2.);
            Ht += asec *
                  (_cmagc * _mLi1[1](w) - _mLi1[4](w) * Ct1 - 2 * _mLi1[5](w) * Ct3 + _mLi1[1](w) * (Ct1 - Ct2 + Ct3) +
                   _mLi1[4](w) * (Ct2 + Ct3) + 2 * (-1 + w) * (dCt1 - dCt2 + dCt3));
            Hs += aseb *
                  (_cmagb * _mLi1[1](w) + (_mLi1[1](w) - (_mLi1[4](w) * (-1 + w)) / (1 + w)) * Cs + 2 * (-1 + w) * dCs);
            Hp += aseb * (_cmagb * _mLi1[1](w) + _mLi1[1](w) * Cv1 - _mLi1[5](w) * (-1 + w) * Cv2 +
                          ((1 + w) * (_mLi1[1](w) - (_mLi1[4](w) * (-1 + w)) / (1 + w)) * (Cv2 + Cv3)) / 2. +
                          2 * (-1 + w) * (dCv1 + ((1 + w) * (dCv2 + dCv3)) / 2.));
            Hm += aseb * (-(_mLi1[4](w) * Cv1) - _mLi1[5](w) * (1 + w) * Cv2 +
                          ((1 + w) * ((_mLi1[1](w) - (_mLi1[4](w) * (-1 + w)) / (1 + w)) * (Cv2 - Cv3) +
                                      2 * (-1 + w) * (dCv2 - dCv3))) /
                              2.);
            Ht += aseb *
                  (_cmagb * _mLi1[1](w) - _mLi1[4](w) * Ct1 + 2 * _mLi1[5](w) * Ct2 + _mLi1[1](w) * (Ct1 - Ct2 + Ct3) -
                   _mLi1[4](w) * (Ct2 + Ct3) + 2 * (-1 + w) * (dCt1 - dCt2 + dCt3));
        }

        // set elements, mapping to amplitude FF basis
        // Fs
        result.element({0, 0}) = (Hs * ((Mb + Mc) * (Mb + Mc) - Sqq)) / (2. * sqMbMc);
        // Fz
        result.element({1, 0}) = (Hp * ((Mb + Mc) * (Mb + Mc) - Sqq)) / (2. * sqMbMc * (Mb + Mc)) +
                                 (Hm * (-(Mb - Mc) * (Mb - Mc) + Sqq)) / (2. * (Mb - Mc) * sqMbMc);
        // Fp
        result.element({2, 0}) = (Hm * (-Mb + Mc)) / (2. * sqMbMc) + (Hp * (Mb + Mc)) / (2. * sqMbMc);
        // Ft
        result.element({3, 0}) = Ht / (2. * sqMbMc);

        // Now create remaining FF tensor entries
        // n indexes rows ie FFs; idx indexes the columns that contract with delta.
        const double rC = Mc / Mb;
        const double sqrC = sqrt(rC);
        const double kHs = Mb * (w + 1) * sqrC;
        const double kHp1 = (w + 1) * sqrC / (1 + rC);
        const double kHm1 = -(-1 + w) * sqrC / (1 - rC);
        const double kHp2 = (1 + rC) / (2. * sqrC);
        const double kHm2 = -(1 - rC) / (2. * sqrC);
        const double kHt = 1. / (2. * Mb * sqrC);

        const double wSq = w * w;
        const double wm1Sq = (w - 1) * (w - 1);
        const double wm1Owp1 = (w - 1) / (w + 1);
        const double wm1SqOwp1 = wm1Sq / (w + 1);
        const double wSqm1 = wSq - 1;

        const vector<vector<double>>& sliwmat = (*getSetting<vector<vector<double>>>("sliwmatrix"));

        // Linearization of hatted FFs
        array<array<double, 10>, 4> FFmat;

        // derivatives of Hs, Hp, Hm, Ht
        array<array<double, 10>, 4> DerMat1 = {
            {{0, 0, -4 * (_eb + _ec) * (-1 + w), -4 * (_eb + _ec) * wm1Sq, 12 * (_eb + _ec) * (-1 + w),
              -2 * (_eb + _ec + 2 * eceb) * wm1Owp1, -2 * (_eb + _ec + 2 * eceb) * wm1SqOwp1,
              -2 * (eb2 + ec2 + eceb) * wm1Sq, -4 * (eb2 + ec2) * (-1 + w), 12 * (eb2 + ec2) * (-1 + w)},
             {0, 0, -4 * (_eb + _ec) * (-1 + w), -4 * (_eb + _ec) * wm1Sq, 12 * (_eb + _ec) * (-1 + w),
              4 * eceb * wm1Owp1, 4 * eceb * wm1SqOwp1, 2 * eceb * wm1Sq, -4 * (eb2 + ec2) * (-1 + w),
              12 * (eb2 + ec2) * (-1 + w)},
             {0, 0, 0, 0, 0, -2 * _eb + 2 * _ec, -2 * (_eb - _ec) * (-1 + w), -2 * (eb2 - ec2) * wSqm1, 0, 0},
             {0, 0, -4 * (_eb + _ec) * (-1 + w), -4 * (_eb + _ec) * wm1Sq, 12 * (_eb + _ec) * (-1 + w),
              -2 * (_eb + _ec + 2 * eceb) + (8 * eceb) / (1 + w),
              (-1 + w) * (-2 * (_eb + _ec + 2 * eceb) + (8 * eceb) / (1 + w)),
              -2 * (-1 + w) * (eb2 + ec2 - eceb + (eb2 + ec2 + eceb) * w), -4 * (eb2 + ec2) * (-1 + w),
              12 * (eb2 + ec2) * (-1 + w)}}};

        array<array<double, 10>, 4> DerMat2{};
        if (bAsM) {
            DerMat2 = {
                {{0, 0, -4 * ((_cmagb + Cs) * _eb + (_cmagc + Cs) * _ec) * (-1 + w),
                  -4 * ((_cmagb + Cs) * _eb + (_cmagc + Cs) * _ec) * wm1Sq,
                  12 * ((_cmagb + Cs) * _eb + (_cmagc + Cs) * _ec) * (-1 + w),
                  (-2 * Cs * (_eb + _ec) * (-1 + w)) / (1 + w), (-2 * Cs * (_eb + _ec) * wm1Sq) / (1 + w), 0, 0, 0},
                 {0, 0,
                  -2 * (-1 + w) *
                      (2 * _cmagb * _eb + 2 * _cmagc * _ec + (_eb + _ec) * (2 * Cv1 + (Cv2 + Cv3) * (1 + w))),
                  -2 * (2 * _cmagb * _eb + 2 * _cmagc * _ec + (_eb + _ec) * (2 * Cv1 + (Cv2 + Cv3) * (1 + w))) * wm1Sq,
                  6 * (-1 + w) *
                      (2 * _cmagb * _eb + 2 * _cmagc * _ec + (_eb + _ec) * (2 * Cv1 + (Cv2 + Cv3) * (1 + w))),
                  -((Cv2 + Cv3) * (_eb + _ec) * (-1 + w)), -((Cv2 + Cv3) * (_eb + _ec) * wm1Sq), 0, 0, 0},
                 {0, 0, -2 * (Cv2 - Cv3) * (_eb + _ec) * wSqm1, -2 * (Cv2 - Cv3) * (_eb + _ec) * (1 + w) * wm1Sq,
                  6 * (Cv2 - Cv3) * (_eb + _ec) * wSqm1, 2 * Cv1 * (-_eb + _ec) - (Cv2 - Cv3) * (_eb + _ec) * (-1 + w),
                  -((2 * Cv1 * (_eb - _ec) + (Cv2 - Cv3) * (_eb + _ec) * (-1 + w)) * (-1 + w)), 0, 0, 0},
                 {0, 0, -4 * (_cmagb * _eb + _cmagc * _ec + (Ct1 - Ct2 + Ct3) * (_eb + _ec)) * (-1 + w),
                  -4 * (_cmagb * _eb + _cmagc * _ec + (Ct1 - Ct2 + Ct3) * (_eb + _ec)) * wm1Sq,
                  12 * (_cmagb * _eb + _cmagc * _ec + (Ct1 - Ct2 + Ct3) * (_eb + _ec)) * (-1 + w),
                  -2 * ((Ct2 + Ct3) * (_eb - _ec) + Ct1 * (_eb + _ec)),
                  -2 * ((Ct2 + Ct3) * (_eb - _ec) + Ct1 * (_eb + _ec)) * (-1 + w), 0, 0, 0}}};
        }

        for (size_t i = 0; i < 10; ++i) {
            FFmat[0][i] = kHs * (DerMat1[0][i] + as * DerMat2[0][i]);
            FFmat[1][i] = kHp1 * (DerMat1[1][i] + as * DerMat2[1][i]) + kHm1 * (DerMat1[2][i] + as * DerMat2[2][i]);
            FFmat[2][i] = kHp2 * (DerMat1[1][i] + as * DerMat2[1][i]) + kHm2 * (DerMat1[2][i] + as * DerMat2[2][i]);
            FFmat[3][i] = kHt * (DerMat1[3][i] + as * DerMat2[3][i]);
        }

        // Linearization of LO IW functions
        const array<double, 10> xiIWL = XiLinearCoeffs(w);

        for (size_t n = 0; n < 4; ++n) {
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

    std::unique_ptr<FormFactorBase> FFBtoDBLPRXPVar::clone(const std::string& label) {
        MAKE_CLONE(FFBtoDBLPRXPVar, label);
    }

} // namespace Hammer
