///
/// @file  FFBtoDstarBGLVar.cc
/// @brief \f$ B \rightarrow D \f$ BGL form factors with variations
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/FormFactors/BGL/FFBtoDstarBGLVar.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFBtoDstarBGLVar::FFBtoDstarBGLVar() {
        // Create tensor rank and dimensions
        IndexList dims = {8, 11}; // size is _mFFErrNames + 1 to include central values in zeroth component
        setGroup("BGLVar");               // override BGL base class FF group setting
        string name{"FFBtoDstarBGLVar"};

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

    void FFBtoDstarBGLVar::defineSettings() {
        setPath(getFFErrPrefixGroup().get());
        setUnits("GeV");

        defineAndAddErrSettings({"delta_a0", "delta_a1", "delta_a2", "delta_b0", "delta_b1", "delta_b2", "delta_c1",
                                 "delta_c2", "delta_d0", "delta_d1"});

        // values from hep-ph/9705252, u = 0.33
        addSetting<double>("Vcb", 41.5e-3);
        addSetting<double>("Chim", 3.068e-4); // GeV^-2
        addSetting<double>("Chip", 5.280e-4); // GeV^-2
        addSetting<double>("ChimL", 2.466e-3);
        // pole mass from hep-ph/9705252
        addSetting<double>("mb", 4.9);       // GeV
        addSetting<double>("mc", 4.9 - 3.4); // GeV (using mb-mc = 3.4 GeV)


        //        vector<double> avec={0.012,0.7,0.8};
        //        vector<double> bvec={0.01223,-0.054, 0.2};
        //        vector<double> cvec={-0.01,0.12};
        // Using 222 fit from 1902.09553
        vector<double> avec = {0.00038, 0.026905, 0.};
        vector<double> bvec = {0.00055, -0.0020370, 0.};
        vector<double> cvec = {-0.000433, 0.005353};
        // Approximate values from 1707.09509 (adapted by rescaling wrt ChimL value)
        vector<double> dvec = {0.007, -0.036};
        addSetting<vector<double>>("avec", avec);
        addSetting<vector<double>>("bvec", bvec);
        addSetting<vector<double>>("cvec", cvec);
        addSetting<vector<double>>("dvec", dvec);

        vector<double> BcStatesf{6.730, 6.736, 7.135, 7.142}; // GeV
        vector<double> BcStatesg{6.337, 6.899, 7.012, 7.280}; // GeV
        vector<double> BcStatesP1{6.275, 6.842, 7.250};       // GeV
        addSetting<vector<double>>("BcStatesf", BcStatesf);
        addSetting<vector<double>>("BcStatesg", BcStatesg);
        addSetting<vector<double>>("BcStatesP1", BcStatesP1);

        // correlation matrix and set zero error eigenvectors
        // Row basis is a's, b's, c's, then d's!
        vector<vector<double>> abcdmat{
            {1., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 1., 0., 0., 0., 0., 0., 0., 0., 0.},
            {0., 0., 1., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 0., 1., 0., 0., 0., 0., 0., 0.},
            {0., 0., 0., 0., 1., 0., 0., 0., 0., 0.}, {0., 0., 0., 0., 0., 1., 0., 0., 0., 0.},
            {0., 0., 0., 0., 0., 0., 1., 0., 0., 0.}, {0., 0., 0., 0., 0., 0., 0., 1., 0., 0.},
            {0., 0., 0., 0., 0., 0., 0., 0., 1., 0.}, {0., 0., 0., 0., 0., 0., 0., 0., 0., 1.}};
        addSetting<vector<vector<double>>>("abcdmatrix", abcdmat);

        setInitialized();
    }

    void FFBtoDstarBGLVar::evalAtPSPoint(const vector<double>& point, const vector<double>& masses) {
        Tensor& result = getTensor();
        result.clearData();

        if (!isInitialized()) {
            MSG_WARNING("Warning, Settings have not been defined!");
        }

        double Mb = 0.;
        double Mc = 0.;
        double unitres = 1.;
        tie(Mb, Mc, unitres) = getParentDaughterHadMasses(masses);

        const double Sqq = point[0];
        const double Mb2 = Mb * Mb;
        const double Mb3 = Mb2 * Mb;
        const double rC = Mc / Mb;
        const double rC2 = rC * rC;
        const double sqrC = sqrt(rC);

        double w = getW(Sqq, Mb, Mc);
        // safety measure if w==1.0
        if (isZero(w - 1.0)) {
            w += 1e-6;
        }
        const double w2 = w * w;

        const vector<double>& ag = (*getSetting<vector<double>>("avec"));
        const vector<double>& af = (*getSetting<vector<double>>("bvec"));
        const vector<double>& aF1 = (*getSetting<vector<double>>("cvec"));
        const vector<double>& aP1 = (*getSetting<vector<double>>("dvec"));
        const size_t ags = ag.size();
        const size_t afs = af.size();
        const size_t aF1s = aF1.size();
        const size_t aP1s = aP1.size();

        const vector<size_t>& zdeg = {ags, afs, aF1s + 1ul, aP1s};

        const size_t nmax = *max_element(zdeg.begin(), zdeg.end());
        const double z = (sqrt(w + 1) - sqrt2) / (sqrt(w + 1) + sqrt2);
        vector<double> zpow{1.};
        for (size_t n = 1; n < nmax; ++n) {
            zpow.push_back(zpow[n - 1] * z);
        }

        const vector<vector<double>>& abcdmat = (*getSetting<vector<vector<double>>>("abcdmatrix"));

        const double nc = 2.6;
        const double etaEWVcb = 1.0066 * (*getSetting<double>("Vcb"));
        const double chim = (*getSetting<double>("Chim")) / (unitres * unitres);
        const double chip = (*getSetting<double>("Chip")) / (unitres * unitres);
        const double chimL = (*getSetting<double>("ChimL"));

        const double tp = (Mb + Mc) * (Mb + Mc) / Mb2;
        const double tm = (Mb - Mc) * (Mb - Mc) / Mb2;
        const double sqtptm = sqrt(tp - tm);

        vector<double>& BcStatesf = (*getSetting<vector<double>>("BcStatesf"));
        vector<double>& BcStatesg = (*getSetting<vector<double>>("BcStatesg"));
        vector<double>& BcStatesP1 = (*getSetting<vector<double>>("BcStatesP1"));

        const double mb = (*getSetting<double>("mb")) * unitres;
        const double mc = (*getSetting<double>("mc")) * unitres;

        double Pf = 1.;
        for (double state : BcStatesf) {
            double sqtpBc = sqrt(tp - pow(state * unitres / Mb, 2));
            Pf *= ((z - ((sqtpBc - sqtptm) / (sqtpBc + sqtptm))) / (1. - z * ((sqtpBc - sqtptm) / (sqtpBc + sqtptm))));
        }

        double PF1 = Pf;

        double Pg = 1.;
        for (double state : BcStatesg) {
            double sqtpBc = sqrt(tp - pow(state * unitres / Mb, 2));
            Pg *= ((z - ((sqtpBc - sqtptm) / (sqtpBc + sqtptm))) / (1. - z * ((sqtpBc - sqtptm) / (sqtpBc + sqtptm))));
        }

        double PP1 = 1.;
        for (double state : BcStatesP1) {
            double sqtpBc = sqrt(tp - pow(state * unitres / Mb, 2));
            PP1 *= ((z - ((sqtpBc - sqtptm) / (sqtpBc + sqtptm))) / (1. - z * ((sqtpBc - sqtptm) / (sqtpBc + sqtptm))));
        }

        const double phig = sqrt(256. * nc / (3 * pi * chip)) * ((rC2 * pow(1 + z, 2) * pow(1 - z, -0.5)) /
                                                                 pow(((1 + rC) * (1 - z) + 2 * sqrC * (1 + z)), 4));
        const double phif = (1. / Mb2) * sqrt(16. * nc / (3 * pi * chim)) *
                            ((rC * (1 + z) * pow(1 - z, 1.5)) / pow(((1 + rC) * (1 - z) + 2 * sqrC * (1 + z)), 4));
        const double phiF1 = (1. / Mb3) * sqrt(8. * nc / (3 * pi * chim)) *
                             ((rC * (1 + z) * pow(1 - z, 2.5)) / pow(((1 + rC) * (1 - z) + 2 * sqrC * (1 + z)), 5));

        const double phif_0 = 4. * rC * sqrt(nc / chim) / (Mb2 * sqrt(3 * pi) * pow(1 + 2 * sqrC + rC, 4));
        const double phiF1_0 = 2. * sqrt(2 / (3 * pi)) * rC * sqrt(nc / chim) / (Mb3 * pow(1 + 2 * sqrC + rC, 5));

        const double phiP1 = sqrt(nc / (pi * chimL)) * (8. * sqrt2 * rC2 * pow(1 + z, 2) * pow(1 - z, -0.5)) /
                             pow(((1 + rC) * (1 - z) + 2 * sqrC * (1 + z)), 4);

        double g = 0;
        for (size_t n = 0; n < ags; ++n) {
            g += ag[n] * zpow[n];
        }
        g /= (Pg * phig);

        double f = 0;
        for (size_t n = 0; n < afs; ++n) {
            f += af[n] * zpow[n];
        }
        f /= (Pf * phif);

        double F1 = af[0] * (Mb - Mc) * phiF1_0 / phif_0;
        for (size_t n = 0; n < aF1s; ++n) {
            F1 += aF1[n] * zpow[n + 1];
        }
        F1 /= (PF1 * phiF1);

        // From 1707.09509 (sqrC/(1+rC) maps from F2)
        double P1 = 0;
        for (size_t n = 0; n < aP1s; ++n) {
            P1 += aP1[n] * zpow[n];
        }
        P1 *= sqrC / ((1 + rC) * PP1 * phiP1);

        // Mapping to amplitude FF basis
        const double Fpf = (rC - w) / (2. * rC * Mb2 * (w2 - 1));
        const double FpF1 = 1. / (2. * rC * Mb3 * (w2 - 1));
        const double Fmf = (rC + w) / (2. * rC * Mb2 * (w2 - 1));
        const double FmF1 = 1. / (2. * rC * Mb3 * (w2 - 1)) * (rC2 - 1) / (1 + rC2 - 2 * rC * w);
        const double FmP1 = sqrC * (rC + 1) / (Mb * (1 + rC2 - 2 * rC * w));

        const double fm = (Fmf * f + FmF1 * F1 + FmP1 * P1);
        const double fp = (Fpf * f + FpF1 * F1);

        // Pseudoscalar from eqn of motion
        const double fps = -(f + Mb2 * (fp * (1 - rC2) + fm * (1 + rC2 - 2 * rC * w))) / (mb + mc);

        // set elements, mapping to amplitude FF basis
        // Fps (dim 0)
        result.element({0, 0}) = fps;
        // Ff (dim 1)
        result.element({1, 0}) = f;
        // Fg (dim -1) Definition in BGL (hep-ph/9705252) differs from Manohar, Wise or 1610.02045 by factor of 2
        result.element({2, 0}) = g / 2.;
        // Fm (dim -1)
        result.element({3, 0}) = fm;
        // Fp (dim -1)
        result.element({4, 0}) = fp;
        // Fzt (dim -2)
        // result.element({5}) = 0;
        // Fmt (dim 0)
        // result.element({6}) = 0;
        // Fpt (dim 0)
        // result.element({7}) = 0;

        // Now create remaining FF tensor entries
        // Logic: FF^x = a^x_i z^i, a^x_i = M_{ij} delta_j, delta_j = {1, delta_...}
        // n indexes rows of abcdmat; idx2 indexes the columns that contract with delta.
        // g: 1st ag.size() rows; f: 2nd af.size() rows; F1: contribution from 1st f row, plus next aF1.size() rows; P1:
        // last aP1.size() rows Sum over appropriate n and linearly combine to generate idx-th entries for each FF
        for (IndexType idx = 0; idx < static_cast<IndexType>(ags + afs + aF1s + aP1s); ++idx) {
            double fentry = 0.;
            double gentry = 0.;
            double F1entry = abcdmat[ags][idx] * (Mb - Mc) * phiF1_0 / phif_0 / (PF1 * phiF1);
            double P1entry = 0.;
            for (size_t n = 0; n < ags; ++n) {
                gentry += abcdmat[n][idx] * zpow[n] / (Pg * phig);
            }
            for (size_t n = 0; n < afs; ++n) {
                fentry += abcdmat[n + ags][idx] * zpow[n] / (Pf * phif);
            }
            for (size_t n = 0; n < aF1s; ++n) {
                F1entry += abcdmat[n + ags + afs][idx] * zpow[n + 1] / (PF1 * phiF1);
            }
            for (size_t n = 0; n < aP1s; ++n) {
                P1entry += abcdmat[n + ags + afs + aF1s][idx] * zpow[n] * sqrC / (1 + rC) / (PP1 * phiP1);
            }
            double fmentry = Fmf * fentry + FmF1 * F1entry + FmP1 * P1entry;
            double fpentry = Fpf * fentry + FpF1 * F1entry;

            const double fpsentry =
                -(fentry + Mb2 * (fpentry * (1 - rC2) + fmentry * (1 + rC2 - 2 * rC * w))) / (mb + mc);
            const auto idxp1 = static_cast<IndexType>(idx + 1u);
            if (!isZero(fpsentry)) {
                result.element({0, idxp1}) = fpsentry;
            }
            if (!isZero(fentry)) {
                result.element({1, idxp1}) = fentry;
            }
            if (!isZero(gentry)) {
                result.element({2, idxp1}) = gentry / 2.;
            }
            if (!isZero(fmentry)) {
                result.element({3, idxp1}) = fmentry;
            }
            if (!isZero(fpentry)) {
                result.element({4, idxp1}) = fpentry;
            }
        }

        result *= (1. / etaEWVcb);
    }

    unique_ptr<FormFactorBase> FFBtoDstarBGLVar::clone(const string& label) {
        MAKE_CLONE(FFBtoDstarBGLVar, label);
    }

} // namespace Hammer
