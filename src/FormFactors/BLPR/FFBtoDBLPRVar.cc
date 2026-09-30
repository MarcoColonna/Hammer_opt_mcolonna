///
/// @file  FFBtoDBLPRVar.cc
/// @brief \f$ B \rightarrow D \f$ BLPRVar form factors with variations
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/FormFactors/BLPR/FFBtoDBLPRVar.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFBtoDBLPRVar::FFBtoDBLPRVar() {
        // Create tensor rank and dimensions
        IndexList dims = {4, 8}; // size is _mFFErrNames + 1 to include central values in zeroth component
        string name{"FFBtoDBLPRVar"};

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

    void FFBtoDBLPRVar::defineSettings() {
        setPath(getFFErrPrefixGroup().get());
        setUnits("GeV");

        FFBLPRVarBase::defineSettings(); // common parameters for BLPRVar are defined in FFBLPRVarBase

        setInitialized();
    }

    void FFBtoDBLPRVar::evalAtPSPoint(const vector<double>& point, const vector<double>& masses) {
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
        const double wSq = w * w;
        const double wm1Sq = (w - 1) * (w - 1);

        const vector<vector<double>>& sliwmat = (*getSetting<vector<vector<double>>>("sliwmatrix"));

        // Linearization of hatted FFs
        vector<vector<double>> FFmat{
            {0, -4 * (eb + ec) * Mb * (-1 + wSq) * sqrC, -4 * (eb + ec) * Mb * wm1Sq * (1 + w) * sqrC,
             12 * (eb + ec) * Mb * (-1 + wSq) * sqrC, -2 * (eb + ec) * Mb * (-1 + w) * sqrC,
             -2 * (eb + ec) * Mb * wm1Sq * sqrC, 0},
            {0, (-4 * (eb + ec) * (-1 + wSq) * sqrC) / (1 + rC), (-4 * (eb + ec) * wm1Sq * (1 + w) * sqrC) / (1 + rC),
             (12 * (eb + ec) * (-1 + wSq) * sqrC) / (1 + rC), (-2 * (eb - ec) * (-1 + w) * sqrC) / (-1 + rC),
             (-2 * (eb - ec) * wm1Sq * sqrC) / (-1 + rC), 0},
            {0, (-2 * (eb + ec) * (-1 + w) * (1 + rC)) / sqrC, (-2 * (eb + ec) * wm1Sq * (1 + rC)) / sqrC,
             (6 * (eb + ec) * (-1 + w) * (1 + rC)) / sqrC, -(((eb - ec) * (-1 + rC)) / sqrC),
             -(((eb - ec) * (-1 + w) * (-1 + rC)) / sqrC), 0},
            {0, (-2 * (eb + ec) * (-1 + w)) / (Mb * sqrC), (-2 * (eb + ec) * wm1Sq) / (Mb * sqrC),
             (6 * (eb + ec) * (-1 + w)) / (Mb * sqrC), -((eb + ec) / (Mb * sqrC)),
             -(((eb + ec) * (-1 + w)) / (Mb * sqrC)), 0}};

        // Linearization of LO IW functions
        const vector<double> xiIWL = XiLinearCoeffs(w);

        for (size_t n = 0; n < 4; ++n) {
            auto ni = static_cast<IndexType>(n);
            for (size_t idx = 0; idx < 7; ++idx) {
                complex<double> entry = 0;
                auto idxi = static_cast<IndexType>(idx + 1u);
                for (size_t idx1 = 0; idx1 < 7; ++idx1) {
                    entry += sliwmat[idx1][idx] * (FFmat[n][idx1] + result.element({ni, 0}) * xiIWL[idx1]);
                }
                result.element({ni, idxi}) = entry;
            }
        }

        result *= h;
        result.toVector();
    }

    std::unique_ptr<FormFactorBase> FFBtoDBLPRVar::clone(const std::string& label) {
        MAKE_CLONE(FFBtoDBLPRVar, label);
    }

} // namespace Hammer
