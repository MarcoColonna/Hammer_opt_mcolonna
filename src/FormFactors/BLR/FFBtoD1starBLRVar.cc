///
/// @file  FFBtoD1starBLRVar.cc
/// @brief \f$ B \rightarrow D_1^* \f$ BLR form factors with variations
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/FormFactors/BLR/FFBtoD1starBLRVar.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFBtoD1starBLRVar::FFBtoD1starBLRVar() {
        // Create tensor rank and dimensions
        IndexList dims = {8, 6}; // size is _mFFErrNames + 1 to include central values in zeroth component
        setGroup("BLRVar");              // override BLPR base class FF group setting
        string name{"FFBtoD1starBLRVar"};

        setPrefix("BtoD**1*");
        _mFFErrLabel = FF_BDSSD1STAR_VAR;
        addProcessSignature(PID::BPLUS, {-PID::DSSD1STAR});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BDSSD1STAR, FF_BDSSD1STAR_VAR})});

        addProcessSignature(PID::BZERO, {PID::DSSD1STARMINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BDSSD1STAR, FF_BDSSD1STAR_VAR})});

        setPrefix("BstoDs**1*");
        _mFFErrLabel = FF_BSDSSDS1STAR_VAR;
        addProcessSignature(PID::BS, {PID::DSSDS1STARMINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BSDSSDS1STAR, FF_BSDSSDS1STAR_VAR})});

        setSignatureIndex();
    }

    void FFBtoD1starBLRVar::defineSettings() {
        setPath(getFFErrPrefixGroup().get());
        setUnits("GeV");

        defineAndAddErrSettings({"delta_zt1", "delta_ztp", "delta_zeta1", "delta_chi1", "delta_chi2"});

        // 1S scheme: mb1S = 4710, mbBar = 5313, lambda1 = -3 10^5 MeV^2, delta mb-mc = 3400, alpha_s = 26/100
        addSetting<double>("as", 0.26);
        addSetting<double>("mb", 4.710);
        addSetting<double>("mc", 4.710 - 3.400);
        addSetting<double>("zt1", 0.7);
        addSetting<double>("ztp", 0.2);
        addSetting<double>("zeta1", 0.6);
        addSetting<double>("chi1", 0.);
        addSetting<double>("chi2", 0.);
        addSetting<double>("laB", 0.4);
        addSetting<double>("laS", 0.76);

        // correlation matrix and set zero error eigenvectors
        // Row basis is zt1, ztp, zeta1, chi1(1), chi2(1)
        vector<vector<double>> sliwmat{{1., 0., 0., 0., 0.},
                                       {0., 1., 0., 0., 0.},
                                       {0., 0., 1., 0., 0.},
                                       {0., 0., 0., 1., 0.},
                                       {0., 0., 0., 0., 1.}};
        addSetting<vector<vector<double>>>("sliwmatrix", sliwmat);

        setInitialized();
    }

    void FFBtoD1starBLRVar::evalAtPSPoint(const vector<double>& point, const vector<double>& masses) {
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
        const double zt1 = (*getSetting<double>("zt1"));
        const double ztp = (*getSetting<double>("ztp"));
        const double zeta1 = (*getSetting<double>("zeta1"));
        const double chi1 = (*getSetting<double>("chi1"));
        const double chi2 = (*getSetting<double>("chi2"));
        const double laB = (*getSetting<double>("laB"));
        const double laS = (*getSetting<double>("laS"));

        const vector<vector<double>>& sliwmat = (*getSetting<vector<vector<double>>>("sliwmatrix"));

        const double LambdaD12 = -laB + laS * w;
        const double Gb = (-(laB * (2 + w)) + laS * (1 + 2 * w)) / (1 + w) - 2 * (w - 1) * zeta1;
        const double LOIWzeta = zt1 + (w - 1) * zt1 * ztp;

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
        const double Gs =
            1 - eC * (LambdaD12 / (1 + w) - 2 * (w - 1) * zeta1 + 2 * chi1 - 2 * (1 + w) * chi2) + as * Cs - eB * Gb;
        const double Gv1 = eC * (LambdaD12 - 2 * (w - 1) * chi1) + (w - 1) * (1 + as * Cv1) - (1 + w) * eB * Gb;
        const double Gv2 = eC * (2 * zeta1 - 2 * chi2) - as * Cv2;
        const double Gv3 =
            -1 - eC * (LambdaD12 / (1 + w) + 2 * zeta1 - 2 * chi1 + 2 * chi2) - as * (Cv1 + Cv3) + eB * Gb;
        const double Ga = 1 + eC * (LambdaD12 / (1 + w) - 2 * chi1) + as * Ca1 - eB * Gb;
        const double Gt1 = -1 + eC * (LambdaD12 / (1 + w) + 2 * chi1) - as * (Ct1 + (w - 1) * Ct2) + eB * Gb;
        const double Gt2 = 1 + eC * (LambdaD12 / (1 + w) - 2 * chi1) + as * (Ct1 - (w - 1) * Ct3) + eB * Gb;
        const double Gt3 = eC * (2 * zeta1 + 2 * chi2) - as * Ct2;


        // defined as 1/LOIWzeta * D[LOIWzeta * f, var]
        const array<double, 5> GsDer = {Gs / zt1, Gs / LOIWzeta * zt1 * (w - 1), 2 * (eB + eC) * (-1 + w), -2 * eC,
                                        2 * eC * (1 + w)};
        const array<double, 5> Gv1Der = {Gv1 / zt1, Gv1 / LOIWzeta * zt1 * (w - 1), 2 * eB * (-1 + w * w),
                                         -2 * eC * (-1 + w), 0};
        const array<double, 5> Gv2Der = {Gv2 / zt1, Gv2 / LOIWzeta * zt1 * (w - 1), 2 * eC, 0, -2 * eC};
        const array<double, 5> Gv3Der = {Gv3 / zt1, Gv3 / LOIWzeta * zt1 * (w - 1), -2 * (eC + eB * (-1 + w)), 2 * eC,
                                         -2 * eC};
        const array<double, 5> GaDer = {Ga / zt1, Ga / LOIWzeta * zt1 * (w - 1), 2 * eB * (-1 + w), -2 * eC, 0};
        const array<double, 5> Gt1Der = {Gt1 / zt1, Gt1 / LOIWzeta * zt1 * (w - 1), -2 * eB * (-1 + w), 2 * eC, 0};
        const array<double, 5> Gt2Der = {Gt2 / zt1, Gt2 / LOIWzeta * zt1 * (w - 1), -2 * eB * (-1 + w), -2 * eC, 0};
        const array<double, 5> Gt3Der = {Gt3 / zt1, Gt3 / LOIWzeta * zt1 * (w - 1), 2 * eC, 0, 2 * eC};

        // Set elements
        result.element({0, 0}) = Gs;
        result.element({1, 0}) = Gv1;
        result.element({2, 0}) = Gv2;
        result.element({3, 0}) = Gv3;
        result.element({4, 0}) = Ga;
        result.element({5, 0}) = Gt1;
        result.element({6, 0}) = Gt2;
        result.element({7, 0}) = Gt3;

        for (IndexType i1 = 1; i1 <= 5; ++i1) {
            for (IndexType i2 = 0; i2 < 5; ++i2) {
                result.element({0, i1}) += sliwmat[i2][i1 - 1] * GsDer[i2];
                result.element({1, i1}) += sliwmat[i2][i1 - 1] * Gv1Der[i2];
                result.element({2, i1}) += sliwmat[i2][i1 - 1] * Gv2Der[i2];
                result.element({3, i1}) += sliwmat[i2][i1 - 1] * Gv3Der[i2];
                result.element({4, i1}) += sliwmat[i2][i1 - 1] * GaDer[i2];
                result.element({5, i1}) += sliwmat[i2][i1 - 1] * Gt1Der[i2];
                result.element({6, i1}) += sliwmat[i2][i1 - 1] * Gt2Der[i2];
                result.element({7, i1}) += sliwmat[i2][i1 - 1] * Gt3Der[i2];
            }
        }

        result *= LOIWzeta;
        result.toVector();
    }

    std::unique_ptr<FormFactorBase> FFBtoD1starBLRVar::clone(const std::string& label) {
        MAKE_CLONE(FFBtoD1starBLRVar, label);
    }

} // namespace Hammer
