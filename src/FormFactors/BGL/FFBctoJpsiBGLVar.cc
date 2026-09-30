///
/// @file  FFBctoJpsiBGLVar.cc
/// @brief \f$ B_c \rightarrow J/\psi \f$ BGL form factors with variations
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/FormFactors/BGL/FFBctoJpsiBGLVar.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFBctoJpsiBGLVar::FFBctoJpsiBGLVar() {
        // Create tensor rank and dimensions
        IndexList dims = {8, 16}; // size is _mFFErrNames + 1 to include central values in zeroth component
        setGroup("BGLVar");               // override BGL base class FF group setting
        string name{"FFBctoJpsiBGLVar"};

        setPrefix("BctoJpsi");
        _mFFErrLabel = FF_BCJPSI_VAR;
        addProcessSignature(PID::BCPLUS, {PID::JPSI});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BCJPSI, FF_BCJPSI_VAR})});

        setSignatureIndex();
    }

    void FFBctoJpsiBGLVar::defineSettings() {
        setPath(getFFErrPrefixGroup().get());
        setUnits("GeV");

        //        defineAndAddErrSettings({"delta_a0","delta_a1","delta_a2","delta_b0","delta_b1","delta_b2","delta_c1","delta_c2",
        //                       "delta_d0","delta_d1","delta_d2"});

        defineAndAddErrSettings({"delta_e0", "delta_e1", "delta_e2", "delta_e3", "delta_e4", "delta_e5", "delta_e6",
                                 "delta_e7", "delta_e8", "delta_e9", "delta_e10", "delta_e11", "delta_e12", "delta_e13",
                                 "delta_e14"});

        // Optimized z-expansion
        addSetting<bool>("OptZ", false);


        //        //Using parameters/fit of 1909.10691. Chi values from private correspondence with 1909.10691 authors
        //        //Uses optimized z expansion
        //        addSetting<double>("Chim", 0.007168); //GeV^-2
        //        addSetting<double>("Chip", 0.012539); //GeV^-2
        //        addSetting<double>("ChimL",0.025042);
        //
        //        vector<double> avec={0.004698, -0.02218, 0.1503};
        //        vector<double> bvec={0.003424, -0.02595, 0.3897};
        //        vector<double> cvec={-0.003164, 0.08731};
        //        vector<double> dvec={0.04011, -0.2132, 0.008157};
        //
        //        vector<double> BcStatesf{6.7305,6.7385,7.1355,7.1435}; //GeV
        //        vector<double> BcStatesg{6.3290,6.8975,7.0065}; //GeV
        //        vector<double> BcStatesP1{6.2749, 6.8710}; //GeV
        //
        //        vector<vector<double>> abcdmat{ {1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
        //                                        {0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
        //                                        {0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0.},
        //                                        {0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0.},
        //                                        {0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0.},
        //                                        {0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0.},
        //                                        {0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0.},
        //                                        {0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0.},
        //                                        {0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0.},
        //                                        {0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0.},
        //                                        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1.}};

        // Using parameters/fit of 2007.06957
        addSetting<double>("Chim", 0.00737058 / (4.78 * 4.78)); // GeV^-2
        addSetting<double>("Chip", 0.0126835 / (4.78 * 4.78));  // GeV^-2
        addSetting<double>("ChimL", 0.0249307);


        vector<double> avec = {0.0257, -0.129, -0.21, -0.09};
        vector<double> bvec = {0.01781, -0.104, 0.30, -0.23};
        vector<double> cvec = {-0.0149, -0.028, 0.21};
        vector<double> dvec = {0.0374, -0.192, -0.34, -0.16};
        addSetting<vector<double>>("avec", avec);
        addSetting<vector<double>>("bvec", bvec);
        addSetting<vector<double>>("cvec", cvec);
        addSetting<vector<double>>("dvec", dvec);

        vector<double> BcStatesf{6.730, 6.736, 7.135, 7.142}; // GeV
        vector<double> BcStatesg{6.337, 6.899, 7.012};        // GeV
        vector<double> BcStatesP1{6.2749, 6.842};             // GeV
        addSetting<vector<double>>("BcStatesf", BcStatesf);
        addSetting<vector<double>>("BcStatesg", BcStatesg);
        addSetting<vector<double>>("BcStatesP1", BcStatesP1);

        // Ratios for t_bd branch point evaluation
        addSetting<double>("Rb*", 5.324 / 6.275);
        addSetting<double>("Rb0", 5.280 / 6.275);
        addSetting<double>("Rd0", 1.865 / 3.096);

        // Covariance eigenvector matrix normalized to sqrt of eigenvalues.
        // Correlations from private correspondence with 2007.06957 authors
        // Row basis is a's, b's, c's, then d's!
        vector<vector<double>> abcdmat{
            {-5.88650e-5, 5.86010e-5, -1.58044e-7, 4.60597e-5, 4.50361e-5, -6.59239e-5, -6.95905e-6, -1.95256e-4,
             6.56077e-5, -2.03392e-4, 3.33453e-6, -2.70740e-4, 1.23183e-3, 6.85636e-6, 0.},
            {-2.13573e-3, 9.40702e-3, -4.84472e-6, -2.28676e-3, -3.59456e-3, 9.87990e-3, 1.03201e-3, 1.17669e-2,
             2.39061e-4, 7.18304e-3, -6.61127e-5, 3.62434e-6, 4.23865e-5, 3.66908e-8, 0.},
            {5.98838e-2, -1.31024e-1, 4.98634e-4, 4.01246e-2, 1.20864e-3, -4.84849e-3, -5.32963e-4, -8.60631e-3,
             1.73284e-3, 5.83595e-3, -3.23448e-5, -2.66721e-6, 9.68678e-6, -2.64376e-8, 0.},
            {-6.96105e-2, 1.47906e-1, -5.88309e-4, -4.53848e-2, 3.09887e-4, -4.96375e-3, -5.47791e-4, -8.26785e-3,
             1.48659e-3, 4.74425e-3, -2.38828e-5, -2.41916e-6, 6.19484e-6, -6.21490e-10, 0.},
            {-4.60543e-6, 1.41196e-5, -8.98518e-6, 3.54436e-6, -4.99904e-6, -8.06371e-6, 8.07436e-5, -4.18273e-6,
             2.43843e-6, -4.59694e-6, -5.25667e-5, 4.56415e-6, 2.44885e-5, -5.91311e-4, 0.},
            {5.24441e-5, -4.13391e-6, -1.57117e-3, 8.80871e-6, 3.78922e-5, 1.17726e-3, -8.17870e-3, -1.01374e-5,
             1.04187e-4, 5.92923e-5, 8.55777e-3, -1.57574e-6, -6.11024e-7, -3.89649e-6, 0.},
            {3.23078e-4, 7.15757e-4, 1.59092e-1, -4.79956e-6, 3.91654e-5, 2.06952e-3, -1.67128e-2, -4.94988e-5,
             -7.86528e-5, -2.52413e-5, -2.31675e-3, 3.06330e-8, 1.45352e-7, -5.58536e-7, 0.},
            {2.60058e-4, 8.23422e-4, 1.79370e-1, -1.56146e-5, -4.37574e-5, -1.83692e-3, 1.47505e-2, 4.35267e-5,
             6.92199e-5, 2.27523e-5, 2.12984e-3, -5.65577e-8, -1.37471e-7, 4.31755e-7, 0.},
            {-1.73212e-4, -5.55322e-4, 7.93676e-6, -2.73502e-4, -1.54895e-3, 1.38471e-3, 1.44363e-4, -4.75921e-4,
             7.02116e-4, -9.48733e-5, -2.13600e-6, 1.68188e-3, -3.19382e-4, -6.00141e-7, 0.},
            {8.33885e-3, 2.09552e-2, -8.94737e-5, 9.56269e-3, 4.12309e-2, -8.68020e-3, -9.50349e-4, 4.73502e-3,
             -1.24498e-3, 8.18800e-4, 1.95465e-5, 1.09959e-4, -8.64251e-6, -1.07985e-7, 0.},
            {6.44976e-2, -1.27610e-1, 4.55625e-4, -9.18163e-2, 4.98466e-3, -5.82995e-4, -6.30203e-5, 4.43582e-4,
             -7.48498e-5, 1.00251e-4, 1.19193e-6, 4.94928e-6, 6.65244e-7, -4.13381e-8, 0.},
            {2.30058e-5, -9.56825e-5, 5.10240e-7, -2.77891e-5, -3.57100e-4, -1.53400e-4, -2.02130e-5, -1.66891e-4,
             -7.01356e-6, -1.50819e-4, 7.33945e-6, 9.14328e-4, 9.54645e-4, 6.11291e-6, 0.},
            {-7.46610e-3, 2.16119e-3, 6.62772e-5, 8.66103e-4, 1.56607e-2, 1.27524e-2, 1.36370e-3, -2.85774e-3,
             9.73306e-3, -1.47416e-3, -9.30442e-5, -1.03248e-4, 1.42499e-5, 1.82884e-7, 0.},
            {2.02434e-1, 5.29407e-2, -3.52061e-4, 1.04593e-3, 9.00007e-3, 1.26944e-2, 1.50970e-3, -6.15327e-3,
             -5.85425e-3, 2.90022e-4, 1.62749e-5, -2.83731e-6, 3.71817e-6, 1.78077e-8, 0.},
            {2.52157e-1, 6.15351e-2, -4.98348e-4, 2.77304e-4, -9.63310e-3, -9.51264e-3, -1.13198e-3, 4.44601e-3,
             5.04974e-3, -3.44721e-4, -1.72538e-5, -4.67612e-6, -2.69712e-6, 2.59282e-9, 0.}};
        addSetting<vector<vector<double>>>("abcdmatrix", abcdmat);

        setInitialized();
    }

    void FFBctoJpsiBGLVar::evalAtPSPoint(const vector<double>& point, const vector<double>& masses) {
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
        const double Mc2 = Mc * Mc;
        const double rC = Mc / Mb;
        const double rC2 = rC * rC;
        const double sqrC = sqrt(rC);

        double w = getW(Sqq, Mb, Mc);
        // safety measure if w==1.0
        if (isZero(w - 1.0)) {
            w += 1e-6;
        }
        const double w2 = w * w;

        const double Rbs = (*getSetting<double>("Rb*"));
        const double Rbz = (*getSetting<double>("Rb0"));
        const double Rdz = (*getSetting<double>("Rd0"));

        const bool OptZ = (*getSetting<bool>("OptZ"));

        // z expansion, using branch points t_bd = (m_B* + m_D)^2 for f,F1,P1 and t_bd = (m_B + m_D)^2 for g.
        // Yields two different z expansions
        // q^2_0 = q_opt^2 for optimized, otherwise q_-^2
        const double sqtbd1 = (Rbs + rC * Rdz);
        const double sqtst1 = sqrt(sqtbd1 * sqtbd1 - rC * rC + 2. * rC * w - 1.);
        const double sqtstm1 = sqrt(sqtbd1 * sqtbd1 - (rC - 1.) * (rC - 1.));
        // const double sqtst01 = sqrt(sqtbd1*sqtstm1);
        const double sqtst01 = OptZ ? sqrt(sqtbd1 * sqtstm1) : sqtstm1;

        const double sqtbd0 = (Rbz + rC * Rdz);
        const double sqtst0 = sqrt(sqtbd0 * sqtbd0 - rC * rC + 2. * rC * w - 1.);
        const double sqtstm0 = sqrt(sqtbd0 * sqtbd0 - (rC - 1.) * (rC - 1.));
        // const double sqtst00 = sqrt(sqtbd0*sqtstm0);
        const double sqtst00 = OptZ ? sqrt(sqtbd0 * sqtstm0) : sqtstm0;

        const double z1 = (sqtst1 - sqtst01) / (sqtst1 + sqtst01);
        const double z1min = (sqtstm1 - sqtst01) / (sqtstm1 + sqtst01);
        const double z0 = (sqtst0 - sqtst00) / (sqtst0 + sqtst00);

        const size_t nmax = 4;
        vector<double> z1pow{1.};
        vector<double> z1minpow{1.};
        vector<double> z0pow{1.};
        for (size_t n = 1; n < nmax; ++n) {
            z1pow.push_back(z1pow[n - 1] * z1);
            z1minpow.push_back(z1minpow[n - 1] * z1min);
            z0pow.push_back(z0pow[n - 1] * z0);
        }

        const vector<double>& ag = (*getSetting<vector<double>>("avec"));
        const vector<double>& af = (*getSetting<vector<double>>("bvec"));
        const vector<double>& aF1 = (*getSetting<vector<double>>("cvec"));
        const vector<double>& aP1 = (*getSetting<vector<double>>("dvec"));

        const vector<vector<double>>& abcdmat = (*getSetting<vector<vector<double>>>("abcdmatrix"));

        const double nc = 1.0;
        const double chim = (*getSetting<double>("Chim")) / (unitres * unitres);
        const double chip = (*getSetting<double>("Chip")) / (unitres * unitres);
        const double chimL = (*getSetting<double>("ChimL"));

        vector<double>& BcStatesf = (*getSetting<vector<double>>("BcStatesf"));
        vector<double>& BcStatesg = (*getSetting<vector<double>>("BcStatesg"));
        vector<double>& BcStatesP1 = (*getSetting<vector<double>>("BcStatesP1"));

        double Pf = 1.;
        for (double state : BcStatesf) {
            double mP = state * unitres;
            double wM = (Mb2 + Mc2 - mP * mP) / (2. * Mb * Mc);
            double sqtst1M = sqrt(sqtbd1 * sqtbd1 - rC * rC + 2. * rC * wM - 1.);
            Pf *= ((z1 - ((sqtst1M - sqtst01) / (sqtst1M + sqtst01))) /
                   (1. - z1 * ((sqtst1M - sqtst01) / (sqtst1M + sqtst01))));
        }
        double PF1 = Pf;

        double Pg = 1.;
        for (double state : BcStatesg) {
            double mP = state * unitres;
            double wM = (Mb2 + Mc2 - mP * mP) / (2. * Mb * Mc);
            double sqtst0M = sqrt(sqtbd0 * sqtbd0 - rC * rC + 2. * rC * wM - 1.);
            Pg *= ((z0 - ((sqtst0M - sqtst00) / (sqtst0M + sqtst00))) /
                   (1. - z0 * ((sqtst0M - sqtst00) / (sqtst0M + sqtst00))));
        }

        double PP1 = 1.;
        for (double state : BcStatesP1) {
            double mP = state * unitres;
            double wM = (Mb2 + Mc2 - mP * mP) / (2. * Mb * Mc);
            double sqtst1M = sqrt(sqtbd1 * sqtbd1 - rC * rC + 2. * rC * wM - 1.);
            PP1 *= ((z1 - ((sqtst1M - sqtst01) / (sqtst1M + sqtst01))) /
                    (1. - z1 * ((sqtst1M - sqtst01) / (sqtst1M + sqtst01))));
        }


        // Outer functions
        const double phig = sqrt(nc / (96 * pi * chip)) * pow(sqtst0 / sqtst00, 0.5) * (sqtst0 + sqtst00) *
                            pow(sqtst0, 1.5) * pow(sqtst0 + sqtstm0, 1.5) / pow(sqtst0 + sqtbd0, 4.);
        const double phif = (1. / Mb2) * sqrt(nc / (24 * pi * chim)) * pow(sqtst1 / sqtst01, 0.5) * (sqtst1 + sqtst01) *
                            pow(sqtst1, 0.5) * pow(sqtst1 + sqtstm1, 0.5) / pow(sqtst1 + sqtbd1, 4.);
        const double phiF1 = phif / (Mb * sqrt2 * (sqtst1 + sqtbd1));
        const double phiP1 = sqrt(nc / (64 * pi * chimL)) * pow(sqtst1 / sqtst01, 0.5) * (sqtst1 + sqtst01) *
                             pow(sqtst1, 1.5) * pow(sqtst1 + sqtstm1, 1.5) / pow(sqtst1 + sqtbd1, 4.);

        // Outer functions at q^2 = q^2_max.
        const double phif_0 = (1. / Mb2) * sqrt(nc / (24 * pi * chim)) * pow(sqtstm1 / sqtst01, 0.5) *
                              (sqtstm1 + sqtst01) * pow(sqtstm1, 0.5) * pow(sqtstm1 + sqtstm1, 0.5) /
                              pow(sqtstm1 + sqtbd1, 4.);
        const double phiF1_0 = phif_0 / (Mb * sqrt2 * (sqtstm1 + sqtbd1));

        double g = 0;
        for (size_t n = 0; n < ag.size(); ++n) {
            g += ag[n] * z0pow[n];
        }
        g /= (Pg * phig);

        double f = 0;
        double f_0 = 0;
        for (size_t n = 0; n < af.size(); ++n) {
            f += af[n] * z1pow[n];
            f_0 += af[n] * z1minpow[n];
        }
        f /= (Pf * phif);

        // Leading F1 coefficient c_0 set by constraint F1(q^2=q^2_max) = (Mb - Mc) f(q^2=q^2_max)
        double F1 = (Mb - Mc) * f_0 * phiF1_0 / phif_0;
        for (size_t n = 0; n < aF1.size(); ++n) {
            F1 -= aF1[n] * z1minpow[n + 1];
        }
        // double F1=af[0]*(Mb-Mc)*phiF1_0/phif_0;
        for (size_t n = 0; n < aF1.size(); ++n) {
            F1 += aF1[n] * z1pow[n + 1];
        }
        F1 /= (PF1 * phiF1);

        // sqrC/(1+rC) rescales outer fn phi_P1 P1 = phi_F2 F2
        double P1 = 0;
        for (size_t n = 0; n < aP1.size(); ++n) {
            P1 += aP1[n] * z1pow[n];
        }
        P1 *= sqrC / ((1 + rC) * PP1 * phiP1);

        // Mapping to amplitude FF basis
        const double Fpf = (rC - w) / (2. * rC * Mb2 * (w2 - 1));
        const double FpF1 = 1. / (2. * rC * Mb3 * (w2 - 1));
        const double Fmf = (rC + w) / (2. * rC * Mb2 * (w2 - 1));
        const double FmF1 = 1. / (2. * rC * Mb3 * (w2 - 1)) * (rC2 - 1) / (1 + rC2 - 2 * rC * w);
        const double FmP1 = sqrC * (rC + 1) / (Mb * (1 + rC2 - 2 * rC * w));

        // set elements, mapping to amplitude FF basis
        // Fs (dim 0)
        // result.element({0}) = 0;
        // Ff (dim 1)
        result.element({1, 0}) = f;
        // Fg (dim -1) Definition in BGL (hep-ph/9705252) differs from Manohar, Wise or 1610.02045 by factor of 2
        result.element({2, 0}) = g / 2.;
        // Fm (dim -1)
        result.element({3, 0}) = (Fmf * f + FmF1 * F1 + FmP1 * P1);
        // Fp (dim -1)
        result.element({4, 0}) = (Fpf * f + FpF1 * F1);
        // Fzt (dim -2)
        // result.element({5}) = 0;
        // Fmt (dim 0)
        // result.element({6}) = 0;
        // Fpt (dim 0)
        // result.element({7}) = 0;


        // Now create remaining FF tensor entries
        const auto Nag = ag.size();
        const auto Naf = af.size();
        const auto NaF1 = aF1.size();
        const auto NaP1 = aP1.size();
        // Logic: FF^x = a^x_i z^i, a^x_i = M_{ij} delta_j, delta_j = {1, delta_...}
        // n indexes rows of abcdmat; idx2 indexes the columns that contract with delta.
        // g: 1st Nag rows; f: 2nd Naf rows; F1: contribution from 1st f row, plus next NaF1 rows; P1: last NaP1 rows
        // Sum over appropriate n and linearly combine to generate idx-th entries for each FF
        for (IndexType idx = 0; idx < static_cast<IndexType>(Nag + Naf + NaF1 + NaP1); ++idx) {
            double fentry = 0.;
            double f_0entry = 0.;
            double gentry = 0.;
            // double F1entry = abcdmat[3][idx]*(Mb-Mc)*phiF1_0/phif_0/(PF1*phiF1);
            double P1entry = 0.;
            for (size_t n = 0; n < Nag; ++n) {
                gentry += abcdmat[n][idx] * z0pow[n] / (Pg * phig);
            }
            for (size_t n = 0; n < Naf; ++n) {
                fentry += abcdmat[n + Nag][idx] * z1pow[n] / (Pf * phif);
                f_0entry += abcdmat[n + Nag][idx] * z1minpow[n];
            }
            double F1entry = (Mb - Mc) * f_0entry * phiF1_0 / phif_0 / (PF1 * phiF1);
            for (size_t n = 0; n < NaF1; ++n) {
                F1entry += abcdmat[n + Nag + Naf][idx] * (z1pow[n + 1] - z1minpow[n + 1]) / (PF1 * phiF1);
            }
            for (size_t n = 0; n < NaP1; ++n) {
                P1entry += abcdmat[n + Nag + Naf + NaF1][idx] * z1pow[n] * sqrC / (1 + rC) / (PP1 * phiP1);
            }
            double Fmentry = Fmf * fentry + FmF1 * F1entry + FmP1 * P1entry;
            double Fpentry = Fpf * fentry + FpF1 * F1entry;

            const auto idxp1 = static_cast<IndexType>(idx + 1u);
            if (!isZero(fentry)) {
                result.element({1, idxp1}) = fentry;
            }
            if (!isZero(gentry)) {
                result.element({2, idxp1}) = gentry / 2.;
            }
            if (!isZero(Fmentry)) {
                result.element({3, idxp1}) = Fmentry;
            }
            if (!isZero(Fpentry)) {
                result.element({4, idxp1}) = Fpentry;
            }
        }
    }

    unique_ptr<FormFactorBase> FFBctoJpsiBGLVar::clone(const string& label) {
        MAKE_CLONE(FFBctoJpsiBGLVar, label);
    }

    void FFBctoJpsiBGLVar::addRefs() const {
        if (!getSettingsHandler()->checkReference("Cohen:2019zev")) {
            string ref =
                "@article{Cohen:2019zev,\n"
                "    author = \"Cohen, Thomas D. and Lamm, Henry and Lebed, Richard F.\",\n"
                "    title = \"{Precision Model-Independent Bounds from Global Analysis of $b \\to c \\ell \\nu$ Form "
                "Factors}\",\n"
                "    eprint = \"1909.10691\",\n"
                "    archivePrefix = \"arXiv\",\n"
                "    primaryClass = \"hep-ph\",\n"
                "    doi = \"10.1103/PhysRevD.100.094503\",\n"
                "    journal = \"Phys. Rev. D\",\n"
                "    volume = \"100\",\n"
                "    number = \"9\",\n"
                "    pages = \"094503\",\n"
                "    year = \"2019\"\n"
                "}\n";
            getSettingsHandler()->addReference("Cohen:2019zev", ref);
        }
    }

} // namespace Hammer
