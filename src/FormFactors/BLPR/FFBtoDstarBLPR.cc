///
/// @file  FFBtoDstarBLPR.cc
/// @brief \f$ B \rightarrow D^* \f$ BLPR form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/FormFactors/BLPR/FFBtoDstarBLPR.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFBtoDstarBLPR::FFBtoDstarBLPR() {
        // Create tensor rank and dimensions
        IndexList dims = {8};
        string name{"FFBtoDstarBLPR"};

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

    void FFBtoDstarBLPR::defineSettings() {
        //_mFFErrNames = ;
        setPath(getFFErrPrefixGroup().get());
        setUnits("GeV");

        FFBLPRBase::defineSettings(); // common parameters for BLPR are defined in FFBLPRBase

        setInitialized();
    }

    void FFBtoDstarBLPR::evalAtPSPoint(const vector<double>& point, const vector<double>& masses) {
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

        // BLPR expansion parameters
        const double eb = (*getSetting<double>("la")) / (*getSetting<double>("mb") * 2.);
        const double ebReb = (*getSetting<double>("ebR/eb"));
        const double ec = (*getSetting<double>("la")) / (*getSetting<double>("mc") * 2.);
        const double ecRec = (*getSetting<double>("ecR/ec"));
        const double as = (*getSetting<double>("as")) / pi;

        // LO IW function
        const double chi1 = 0;
        auto xi = Xi(w);
        const double h = xi.first / xi.second + 2. * (eb + ec) * chi1;

        // QCD correction functions
        //  const double Cs  = CS(w, zBC);
        const double Cps = _asCorrections.CP(w, _zBC);
        const double Cv1 = _asCorrections.CV1(w, _zBC);
        // const double Cv2 = CV2(w, zBC);
        // const double Cv3 = CV3(w, zBC);
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
        Hps += ec * (_mLi[2](w) + _mLi[3](w) * (w - 1) + _mLi[5](w) - _mLi[6](w) * (w + 1)) +
               eb * (_mLi[1](w) - _mLi[4](w));
        Hv += ec * (_mLi[2](w) - _mLi[5](w)) + eb * (_mLi[1](w) - _mLi[4](w));
        Ha1 += ec * (_mLi[2](w) - _mLi[5](w) * (w - 1) / (w + 1)) + eb * (_mLi[1](w) - _mLi[4](w) * (w - 1) / (w + 1));
        Ha2 += ec * (_mLi[3](w) + _mLi[6](w));
        Ha3 += ec * (_mLi[2](w) - _mLi[3](w) - _mLi[5](w) + _mLi[6](w)) + eb * (_mLi[1](w) - _mLi[4](w));
        Ht1 += ec * (_mLi[2](w)) + eb * (_mLi[1](w));
        Ht2 += ec * (_mLi[5](w)) - eb * (_mLi[4](w));
        Ht3 += ec * (_mLi[6](w) - _mLi[3](w));

        // Upsilon expansion
        const double corrb = eb * (1. - ebReb);
        const double corrc = ec * (1. - ecRec);
        Hps += -(corrc + corrb);
        Hv += -(corrc + corrb);
        Ha1 += -(corrc + corrb) * (w - 1) / (w + 1);
        Ha2 += 2. * corrc / (w + 1.);
        Ha3 += -(corrc * (w - 1.) / (w + 1.) + corrb);
        Ht1 += 0.;
        Ht2 += corrc - corrb;
        Ht3 += 2. * corrc / (w + 1.);

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

    std::unique_ptr<FormFactorBase> FFBtoDstarBLPR::clone(const std::string& label) {
        MAKE_CLONE(FFBtoDstarBLPR, label);
    }

} // namespace Hammer
