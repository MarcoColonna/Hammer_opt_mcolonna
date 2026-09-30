///
/// @file  FFBtoD1BLRVar.cc
/// @brief \f$ B \rightarrow D_1 \f$ BLR form factors with variations
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/FormFactors/BLR/FFBtoD1BLRVar.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFBtoD1BLRVar::FFBtoD1BLRVar() {
        // Create tensor rank and dimensions
        IndexList dims = {8, 8}; // size is _mFFErrNames + 1 to include central values in zeroth component
        setGroup("BLRVar");              // override BLPR base class FF group setting
        string name{"FFBtoD1BLRVar"};

        setPrefix("BtoD**1");
        _mFFErrLabel = FF_BDSSD1_VAR;
        addProcessSignature(PID::BPLUS, {-PID::DSSD1});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BDSSD1, FF_BDSSD1_VAR})});

        addProcessSignature(PID::BZERO, {PID::DSSD1MINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BDSSD1, FF_BDSSD1_VAR})});

        setPrefix("BstoDs**1");
        _mFFErrLabel = FF_BSDSSDS1_VAR;
        addProcessSignature(PID::BS, {PID::DSSDS1MINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BSDSSDS1, FF_BSDSSDS1_VAR})});

        setSignatureIndex();
    }

    void FFBtoD1BLRVar::defineSettings() {
        setPath(getFFErrPrefixGroup().get());
        setUnits("GeV");

        defineAndAddErrSettings(
            {"delta_t1", "delta_tp", "delta_tau1", "delta_tau2", "delta_eta1", "delta_eta2", "delta_eta3"});

        // 1S scheme: mb1S = 4710, mbBar = 5313, lambda1 = -3 10^5 MeV^2, delta mb-mc = 3400, alpha_s = 26/100
        addSetting<double>("as", 0.26);
        addSetting<double>("mb", 4.710);
        addSetting<double>("mc", 4.710 - 3.400);
        addSetting<double>("t1", 0.7);
        addSetting<double>("tp", -1.6);
        addSetting<double>("tau1", -0.5);
        addSetting<double>("tau2", 2.9);
        addSetting<double>("eta1", 0.);
        addSetting<double>("eta2", 0.);
        addSetting<double>("eta3", 0.);
        addSetting<double>("laB", 0.4);
        addSetting<double>("laP", 0.8);

        // correlation matrix and set zero error eigenvectors
        // Row basis is t1, tp, tau1, tau2, eta1(1), eta2(1), eta3(1)
        vector<vector<double>> sliwmat{{1., 0., 0., 0., 0., 0., 0.}, {0., 1., 0., 0., 0., 0., 0.},
                                       {0., 0., 1., 0., 0., 0., 0.}, {0., 0., 0., 1., 0., 0., 0.},
                                       {0., 0., 0., 0., 1., 0., 0.}, {0., 0., 0., 0., 0., 1., 0.},
                                       {0., 0., 0., 0., 0., 0., 1.}};
        addSetting<vector<vector<double>>>("sliwmatrix", sliwmat);

        setInitialized();
    }

    void FFBtoD1BLRVar::evalAtPSPoint(const vector<double>& point, const vector<double>& masses) {
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

        // BLR parameters
        const double zBC = (*getSetting<double>("mc")) / (*getSetting<double>("mb"));
        const double eB = 1. / (*getSetting<double>("mb") * 2.);
        const double eC = 1. / (*getSetting<double>("mc") * 2.);
        const double as = (*getSetting<double>("as")) / pi;
        const double t1 = (*getSetting<double>("t1"));
        const double tp = (*getSetting<double>("tp"));
        const double tau1 = (*getSetting<double>("tau1"));
        const double tau2 = (*getSetting<double>("tau2"));
        const double eta1 = (*getSetting<double>("eta1"));
        const double eta2 = (*getSetting<double>("eta2"));
        const double eta3 = (*getSetting<double>("eta3"));
        const double laB = (*getSetting<double>("laB"));
        const double laP = (*getSetting<double>("laP"));

        const vector<vector<double>>& sliwmat = (*getSetting<vector<vector<double>>>("sliwmatrix"));

        const double LambdaD32 = -laB + laP * w;
        const double Fb = laB + laP - tau2 - tau1 * (1 + 2 * w);
        const double LOIWtau = t1 + (w - 1) * t1 * tp;

        // QCD correction functions
        const double Cs = _asCorrections.CS(w, zBC);
        // const double Cps = CP(w, zBC);
        const double Cv1 = _asCorrections.CV1(w, zBC);
        const double Cv2 = _asCorrections.CV2(w, zBC);
        const double Cv3 = _asCorrections.CV3(w, zBC);
        const double Ca1 = _asCorrections.CA1(w, zBC);
        // const double Ca2 = CA2(w, zBC);
        // const double Ca3 = CA3(w, zBC);
        const double Ct1 = _asCorrections.CT1(w, zBC);
        const double Ct2 = _asCorrections.CT2(w, zBC);
        const double Ct3 = _asCorrections.CT3(w, zBC);

        // Form factors
        const double Fsc = (-(eC * (4 * LambdaD32 + 2 * (1 + w) * (6 * eta1 + 2 * (w - 1) * eta2 - eta3) -
                                    2 * (w - 1) * ((1 + 2 * w) * tau1 + tau2))) -
                            2 * (1 + w) * (1 + as * Cs) - 2 * (w - 1) * eB * Fb) /
                           sqrt(6.);
        const double Fv1 =
            (-(eC * (4 * (1 + w) * LambdaD32 - (-1 + w * w) * (2 * eta1 + 3 * eta3 + 3 * tau1 - 3 * tau2))) +
             (1 - w * w) * (1 + as * Cv1) - (-1 + w * w) * eB * Fb) /
            sqrt(6.);
        const double Fv2 = (-3 - eC * (10 * eta1 + 4 * (w - 1) * eta2 - 5 * eta3 + (-1 + 4 * w) * tau1 + 5 * tau2) -
                            as * (3 * Cv1 + 2 * (1 + w) * Cv2) - 3 * eB * Fb) /
                           sqrt(6.);
        const double Fv3 = (-2 + w +
                            eC * (4 * LambdaD32 - 2 * (6 + w) * eta1 - 4 * (w - 1) * eta2 - (-2 + 3 * w) * eta3 +
                                  (2 + w) * tau1 + (2 + 3 * w) * tau2) -
                            as * ((2 - w) * Cv1 + 2 * (1 + w) * Cv3) + (2 + w) * eB * Fb) /
                           sqrt(6.);
        const double Fa = (-(eC * (4 * LambdaD32 - (1 + w) * (2 * eta1 + 3 * eta3) - 3 * (w - 1) * (tau1 - tau2))) +
                           (-1 - w) * (1 + as * Ca1) - (w - 1) * eB * Fb) /
                          sqrt(6.);
        const double Ft1 = (-(eC * (4 * LambdaD32 + (1 + w) * (2 * eta1 + 3 * eta3) - 3 * (w - 1) * (tau1 - tau2))) +
                            (1 + w) * (1 + as * (Ct1 + (w - 1) * Ct2)) + (w - 1) * eB * Fb) /
                           sqrt(6.);
        const double Ft2 = (-(eC * (4 * LambdaD32 - (1 + w) * (2 * eta1 + 3 * eta3) - 3 * (w - 1) * (tau1 - tau2))) +
                            (-1 - w) * (1 + as * (Ct1 - (w - 1) * Ct3)) + (w - 1) * eB * Fb) /
                           sqrt(6.);
        const double Ft3 = (3 - eC * (-10 * eta1 - 4 * (w - 1) * eta2 + 5 * eta3 + (-1 + 4 * w) * tau1 + 5 * tau2) +
                            as * (3 * Ct1 - (2 - w) * Ct2 + 3 * Ct3) + 3 * eB * Fb) /
                           sqrt(6.);


        // defined as 1/LOIWtau * D[LOIWtau * f, var]
        const array<double, 7> FscDer = {Fsc / t1,
                                         Fsc / LOIWtau * t1 * (w - 1),
                                         sqrt(2. / 3.) * (eB + eC) * (-1 + w) * (1 + 2 * w),
                                         sqrt(2. / 3.) * (eB + eC) * (-1 + w),
                                         -2 * sqrt(6.) * eC * (1 + w),
                                         -2 * sqrt(2. / 3.) * eC * (-1 + w * w),
                                         sqrt(2. / 3.) * eC * (1 + w)};
        const array<double, 7> Fv1Der = {Fv1 / t1,
                                         Fv1 / LOIWtau * t1 * (w - 1),
                                         ((eB + 3 * eC + 2 * eB * w) * (-1 + w * w)) / sqrt(6.),
                                         ((eB - 3 * eC) * (-1 + w * w)) / sqrt(6.),
                                         sqrt(2. / 3.) * eC * (-1 + w * w),
                                         0.,
                                         sqrt(3. / 2.) * eC * (-1 + w * w)};
        const array<double, 7> Fv2Der = {Fv2 / t1,
                                         Fv2 / LOIWtau * t1 * (w - 1),
                                         (eC - 4 * eC * w + eB * (3 + 6 * w)) / sqrt(6.),
                                         (3 * eB - 5 * eC) / sqrt(6.),
                                         -5 * sqrt(2. / 3.) * eC,
                                         -2 * sqrt(2. / 3.) * eC * (-1 + w),
                                         (5 * eC) / sqrt(6.)};
        const array<double, 7> Fv3Der = {Fv3 / t1,
                                         Fv3 / LOIWtau * t1 * (w - 1),
                                         -(((2 + w) * (eB - eC + 2 * eB * w)) / sqrt(6.)),
                                         (-eB * (2 + w) + eC * (2 + 3 * w)) / sqrt(6.),
                                         -sqrt(2. / 3.) * eC * (6 + w),
                                         -2 * sqrt(2. / 3.) * eC * (-1 + w),
                                         (eC * (2 - 3 * w)) / sqrt(6.)};
        const array<double, 7> FaDer = {Fa / t1,
                                        Fa / LOIWtau * t1 * (w - 1),
                                        ((-1 + w) * (eB + 3 * eC + 2 * eB * w)) / sqrt(6.),
                                        ((eB - 3 * eC) * (-1 + w)) / sqrt(6.),
                                        sqrt(2. / 3.) * eC * (1 + w),
                                        0.,
                                        sqrt(3. / 2.) * eC * (1 + w)};
        const array<double, 7> Ft1Der = {Ft1 / t1,
                                         Ft1 / LOIWtau * t1 * (w - 1),
                                         -(((-1 + w) * (eB - 3 * eC + 2 * eB * w)) / sqrt(6.)),
                                         -(((eB + 3 * eC) * (-1 + w)) / sqrt(6.)),
                                         -sqrt(2. / 3.) * eC * (1 + w),
                                         0.,
                                         -sqrt(3. / 2.) * eC * (1 + w)};
        const array<double, 7> Ft2Der = {Ft2 / t1,
                                         Ft2 / LOIWtau * t1 * (w - 1),
                                         -(((-1 + w) * (eB - 3 * eC + 2 * eB * w)) / sqrt(6.)),
                                         -(((eB + 3 * eC) * (-1 + w)) / sqrt(6.)),
                                         sqrt(2. / 3.) * eC * (1 + w),
                                         0.,
                                         sqrt(3. / 2.) * eC * (1 + w)};
        const array<double, 7> Ft3Der = {Ft3 / t1,
                                         Ft3 / LOIWtau * t1 * (w - 1),
                                         (eC - 4 * eC * w - 3 * eB * (1 + 2 * w)) / sqrt(6.),
                                         -((3 * eB + 5 * eC) / sqrt(6.)),
                                         5 * sqrt(2. / 3.) * eC,
                                         2 * sqrt(2. / 3.) * eC * (-1 + w),
                                         -((5 * eC) / sqrt(6.))};

        // Set elements
        result.element({0, 0}) = Fsc;
        result.element({1, 0}) = Fv1;
        result.element({2, 0}) = Fv2;
        result.element({3, 0}) = Fv3;
        result.element({4, 0}) = Fa;
        result.element({5, 0}) = Ft1;
        result.element({6, 0}) = Ft2;
        result.element({7, 0}) = Ft3;

        for (IndexType i1 = 1; i1 <= 7; ++i1) {
            for (IndexType i2 = 0; i2 < 7; ++i2) {
                result.element({0, i1}) += sliwmat[i2][i1 - 1] * FscDer[i2];
                result.element({1, i1}) += sliwmat[i2][i1 - 1] * Fv1Der[i2];
                result.element({2, i1}) += sliwmat[i2][i1 - 1] * Fv2Der[i2];
                result.element({3, i1}) += sliwmat[i2][i1 - 1] * Fv3Der[i2];
                result.element({4, i1}) += sliwmat[i2][i1 - 1] * FaDer[i2];
                result.element({5, i1}) += sliwmat[i2][i1 - 1] * Ft1Der[i2];
                result.element({6, i1}) += sliwmat[i2][i1 - 1] * Ft2Der[i2];
                result.element({7, i1}) += sliwmat[i2][i1 - 1] * Ft3Der[i2];
            }
        }

        result *= LOIWtau;
        result.toVector();
    }

    std::unique_ptr<FormFactorBase> FFBtoD1BLRVar::clone(const std::string& label) {
        MAKE_CLONE(FFBtoD1BLRVar, label);
    }

} // namespace Hammer
