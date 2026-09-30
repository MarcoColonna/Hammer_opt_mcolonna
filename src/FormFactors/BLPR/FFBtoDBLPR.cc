///
/// @file  FFBtoDBLPR.cc
/// @brief \f$ B \rightarrow D \f$ BLPR form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/FormFactors/BLPR/FFBtoDBLPR.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFBtoDBLPR::FFBtoDBLPR() {
        // Create tensor rank and dimensions
        IndexList dims = {4};
        string name{"FFBtoDBLPR"};

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

    void FFBtoDBLPR::defineSettings() {
        //_mFFErrNames = ;
        setPath(getFFErrPrefixGroup().get());
        setUnits("GeV");

        FFBLPRBase::defineSettings(); // common parameters for BLPR are defined in FFBLPRBase

        setInitialized();
    }

    void FFBtoDBLPR::evalAtPSPoint(const vector<double>& point, const vector<double>& masses) {
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
        Hs += (ec + eb) * (_mLi[1](w) - _mLi[4](w) * (w - 1) / (w + 1));
        Hp += (ec + eb) * _mLi[1](w);
        Hm += (ec - eb) * _mLi[4](w);
        Ht += (ec + eb) * (_mLi[1](w) - _mLi[4](w));

        // Upsilon expansion
        const double corrb = eb * (1. - ebReb);
        const double corrc = ec * (1. - ecRec);
        Hs += -(corrc + corrb) * (w - 1) / (w + 1);
        Hp += 0.;
        Hm += (corrc - corrb);
        Ht += -(corrc + corrb);

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

    std::unique_ptr<FormFactorBase> FFBtoDBLPR::clone(const std::string& label) {
        MAKE_CLONE(FFBtoDBLPR, label);
    }

} // namespace Hammer
