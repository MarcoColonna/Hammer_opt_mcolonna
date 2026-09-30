///
/// @file  FFBtoD0starBLRVar.cc
/// @brief \f$ B \rightarrow D_0^* \f$ BLR form factors with variations
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/FormFactors/BLR/FFBtoD0starBLRVar.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFBtoD0starBLRVar::FFBtoD0starBLRVar() {
        // Create tensor rank and dimensions
        IndexList dims = {4, 6}; // size is _mFFErrNames + 1 to include central values in zeroth component
        setGroup("BLRVar");              // override BLPR base class FF group setting
        string name{"FFBtoD0starBLRVar"};

        setPrefix("BtoD**0*");
        _mFFErrLabel = FF_BDSSD0STAR_VAR;
        addProcessSignature(PID::BPLUS, {-PID::DSSD0STAR});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BDSSD0STAR, FF_BDSSD0STAR_VAR})});

        addProcessSignature(PID::BZERO, {PID::DSSD0STARMINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BDSSD0STAR, FF_BDSSD0STAR_VAR})});

        setPrefix("BstoDs**0*");
        _mFFErrLabel = FF_BSDSSDS0STAR_VAR;
        addProcessSignature(PID::BS, {PID::DSSDS0STARMINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BSDSSDS0STAR, FF_BSDSSDS0STAR_VAR})});

        setSignatureIndex();
    }

    void FFBtoD0starBLRVar::defineSettings() {
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

    void FFBtoD0starBLRVar::evalAtPSPoint(const vector<double>& point, const vector<double>& masses) {
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
        //  const double Cs  = CS(w, zBC);
        const double Cps = _asCorrections.CP(w, zBC);
        // const double Cv1 = CV1(w, zBC);
        // const double Cv2 = CV2(w, zBC);
        // const double Cv3 = CV3(w, zBC);
        const double Ca1 = _asCorrections.CA1(w, zBC);
        const double Ca2 = _asCorrections.CA2(w, zBC);
        const double Ca3 = _asCorrections.CA3(w, zBC);
        const double Ct1 = _asCorrections.CT1(w, zBC);
        // const double Ct2 = CT2(w, zBC);
        // const double Ct3 = CT3(w, zBC);

        // Form factors
        const double gps = eC * (3 * LambdaD12 - 2 * (-1 + w * w) * zeta1 + (w - 1) * (6 * chi1 - 2 * (1 + w) * chi2)) +
                           (w - 1) * (1 + as * Cps) - (1 + w) * eB * Gb;
        const double gp =
            -(eC * ((3 * LambdaD12) / (1 + w) - 2 * (w - 1) * zeta1)) + ((w - 1) * as * (Ca2 + Ca3)) / 2. - eB * Gb;
        const double gm = 1 + eC * (6 * chi1 - 2 * (1 + w) * chi2) + as * (Ca1 + ((w - 1) * (Ca2 - Ca3)) / 2.);
        const double gt = 1 + eC * ((3 * LambdaD12) / (1 + w) - 2 * (w - 1) * zeta1 + 6 * chi1 - 2 * (1 + w) * chi2) +
                          as * Ct1 - eB * Gb;

        // defined as 1/LOIWzeta * D[LOIWzeta * f, var]
        const array<double, 5> gpsDer = {gps / zt1, gps / LOIWzeta * zt1 * (w - 1), 2 * (eB - eC) * (-1 + w * w),
                                         6 * eC * (-1 + w), -2 * eC * (-1 + w * w)};
        const array<double, 5> gpDer = {gp / zt1, gp / LOIWzeta * zt1 * (w - 1), 2 * (eB + eC) * (-1 + w), 0., 0.};
        const array<double, 5> gmDer = {gm / zt1, gm / LOIWzeta * zt1 * (w - 1), 0., 6 * eC, -2 * eC * (1 + w)};
        const array<double, 5> gtDer = {gt / zt1, gt / LOIWzeta * zt1 * (w - 1), 2 * (eB - eC) * (-1 + w), 6 * eC,
                                        -2 * eC * (1 + w)};

        // Set elements
        result.element({0, 0}) = gps;
        result.element({1, 0}) = gp;
        result.element({2, 0}) = gm;
        result.element({3, 0}) = gt;

        for (IndexType i1 = 1; i1 <= 5; ++i1) {
            for (IndexType i2 = 0; i2 < 5; ++i2) {
                result.element({0, i1}) += sliwmat[i2][i1 - 1] * gpsDer[i2];
                result.element({1, i1}) += sliwmat[i2][i1 - 1] * gpDer[i2];
                result.element({2, i1}) += sliwmat[i2][i1 - 1] * gmDer[i2];
                result.element({3, i1}) += sliwmat[i2][i1 - 1] * gtDer[i2];
            }
        }

        result *= LOIWzeta;
        result.toVector();
    }

    std::unique_ptr<FormFactorBase> FFBtoD0starBLRVar::clone(const std::string& label) {
        MAKE_CLONE(FFBtoD0starBLRVar, label);
    }

} // namespace Hammer
