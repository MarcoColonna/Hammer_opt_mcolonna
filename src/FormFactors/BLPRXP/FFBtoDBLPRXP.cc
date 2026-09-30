///
/// @file  FFBtoDBLPRXP.cc
/// @brief \f$ B \rightarrow D \f$ BLPRXP form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/FormFactors/BLPRXP/FFBtoDBLPRXP.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFBtoDBLPRXP::FFBtoDBLPRXP() {
        // Create tensor rank and dimensions
        IndexList dims = {4};
        string name{"FFBtoDBLPRXP"};

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

    void FFBtoDBLPRXP::defineSettings() {
        //_mFFErrNames = ;
        setPath(getFFErrPrefixGroup().get());
        setUnits("GeV");

        FFBLPRXPBase::defineSettings(); // common parameters for BLPRXP are defined in FFBLPRXPBase

        setInitialized();
    }

    void FFBtoDBLPRXP::evalAtPSPoint(const vector<double>& point, const vector<double>& masses) {
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

        // recalc BLPRXP derived expansion parameters if settings changed
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

        if (*getSetting<bool>("With1OverMb2")) {
            const double eb2 = _eb * _eb;
            Hs += +eb2 * (_mLi2[1](w) - _mLi2[4](w) * (w - 1) / (w + 1));
            Hp += +eb2 * _mLi2[1](w);
            Hm += -eb2 * _mLi2[4](w);
            Ht += +eb2 * (_mLi2[1](w) - _mLi2[4](w));
        }
        if (*getSetting<bool>("With1OverMbMc")) {
            const double eceb = _eb * _ec;
            Hs += eceb * (_mMi[1](w) + _mMi[8](w) - 2. * (w - 1.) / (w + 1.) * _mMi[15](w));
            Hp += eceb * (_mMi[1](w) - _mMi[8](w));
            Hm += 0.;
            Ht += eceb * (_mMi[1](w) + _mMi[8](w) - 2. * _mMi[15](w));
        }

        // as/m
        if (*getSetting<bool>("WithAsOverM")) {
            const double asec = as * _ec;
            const double aseb = as * _eb;
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
        result.element({0}) = (Hs * ((Mb + Mc) * (Mb + Mc) - Sqq)) / (2. * sqMbMc);
        // Fz
        result.element({1}) = (Hp * ((Mb + Mc) * (Mb + Mc) - Sqq)) / (2. * sqMbMc * (Mb + Mc)) +
                              (Hm * (-(Mb - Mc) * (Mb - Mc) + Sqq)) / (2. * (Mb - Mc) * sqMbMc);
        // Fp
        result.element({2}) = (Hm * (-Mb + Mc)) / (2. * sqMbMc) + (Hp * (Mb + Mc)) / (2. * sqMbMc);
        // Ft
        result.element({3}) = Ht / (2. * sqMbMc);

        result *= h;
    }

    std::unique_ptr<FormFactorBase> FFBtoDBLPRXP::clone(const std::string& label) {
        MAKE_CLONE(FFBtoDBLPRXP, label);
    }

} // namespace Hammer
