///
/// @file  FFBtoDstarBLPRVar.cc
/// @brief \f$ B \rightarrow D^* \f$ BLPRVar form factors with variations
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/FormFactors/BLPR/FFBtoDstarBLPRVar.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFBtoDstarBLPRVar::FFBtoDstarBLPRVar() {
        // Create tensor rank and dimensions
        IndexList dims = {8, 8}; // size is _mFFErrNames + 1 to include central values in zeroth component
        string name{"FFBtoDstarBLPRVar"};

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

    void FFBtoDstarBLPRVar::defineSettings() {
        setPath(getFFErrPrefixGroup().get());
        setUnits("GeV");

        FFBLPRVarBase::defineSettings(); // common parameters for BLPRVar are defined in FFBLPRVarBase

        setInitialized();
    }

    void FFBtoDstarBLPRVar::evalAtPSPoint(const vector<double>& point, const vector<double>& masses) {
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

        // L and helper variables
        const double corrb = eb * (1. - ebReb);
        const double corrc = ec * (1. - ecRec);

        // Upsilon expansion
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
        const double sqrC = sqrt(rC);
        const double wSq = w * w;
        const double wm1Sq = (w - 1) * (w - 1);

        const vector<vector<double>>& sliwmat = (*getSetting<vector<vector<double>>>("sliwmatrix"));

        // Linearization of hatted FFs
        vector<vector<double>> FFmat{
            {0, 4 * (eb - ec) * (-1 + w) * sqrC, 4 * (eb - ec) * wm1Sq * sqrC, -4 * (3 * eb - ec) * (-1 + w) * sqrC,
             2 * (eb - ec) * sqrC, 2 * (eb - ec) * (-1 + w) * sqrC, 0},
            {0, -4 * eb * Mb * (-1 + wSq) * sqrC, -4 * eb * Mb * wm1Sq * (1 + w) * sqrC,
             4 * (3 * eb - ec) * Mb * (-1 + wSq) * sqrC, -2 * eb * Mb * (-1 + w) * sqrC, -2 * eb * Mb * wm1Sq * sqrC,
             0},
            {0, (-2 * eb * (-1 + w)) / (Mb * sqrC), (-2 * eb * wm1Sq) / (Mb * sqrC),
             (2 * (3 * eb - ec) * (-1 + w)) / (Mb * sqrC), -(eb / (Mb * sqrC)), -((eb * (-1 + w)) / (Mb * sqrC)), 0},
            {0, (-2 * (eb * (-1 + w) + ec * (1 + rC))) / (Mb * sqrC),
             (-2 * (-1 + w) * (eb * (-1 + w) + ec * (1 + rC))) / (Mb * sqrC),
             (2 * (3 * eb - ec) * (-1 + w)) / (Mb * sqrC), -((eb + ec + eb * w - ec * rC) / (Mb * (1 + w) * sqrC)),
             -(((-1 + w) * (eb + ec + eb * w - ec * rC)) / (Mb * (1 + w) * sqrC)), 0},
            {0, (2 * (ec + eb * (-1 + w) - ec * rC)) / (Mb * sqrC),
             (2 * (-1 + w) * (ec + eb * (-1 + w) - ec * rC)) / (Mb * sqrC),
             (-2 * (3 * eb - ec) * (-1 + w)) / (Mb * sqrC), (eb + ec + eb * w + ec * rC) / (Mb * (1 + w) * sqrC),
             ((-1 + w) * (eb + ec + eb * w + ec * rC)) / (Mb * (1 + w) * sqrC), 0},
            {0, (-2 * ec) / (Mb2 * sqrC), (-2 * ec * (-1 + w)) / (Mb2 * sqrC), 0, -(ec / (Mb2 * (1 + w) * sqrC)),
             -((ec * (-1 + w)) / (Mb2 * (1 + w) * sqrC)), 0},
            {0, (2 * eb * (-1 + w) * (-1 + rC)) / sqrC, (2 * eb * wm1Sq * (-1 + rC)) / sqrC,
             (-2 * (3 * eb - ec) * (-1 + w) * (-1 + rC)) / sqrC, (eb * (1 + rC)) / sqrC,
             (eb * (-1 + w) * (1 + rC)) / sqrC, 0},
            {0, (2 * eb * (-1 + w) * (1 + rC)) / sqrC, (2 * eb * wm1Sq * (1 + rC)) / sqrC,
             (-2 * (3 * eb - ec) * (-1 + w) * (1 + rC)) / sqrC, (eb * (-1 + rC)) / sqrC,
             (eb * (-1 + w) * (-1 + rC)) / sqrC, 0}};

        // Linearization of LO IW functions
        const vector<double> xiIWL = XiLinearCoeffs(w);

        for (size_t n = 0; n < 8; ++n) {
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

    std::unique_ptr<FormFactorBase> FFBtoDstarBLPRVar::clone(const std::string& label) {
        MAKE_CLONE(FFBtoDstarBLPRVar, label);
    }

} // namespace Hammer
