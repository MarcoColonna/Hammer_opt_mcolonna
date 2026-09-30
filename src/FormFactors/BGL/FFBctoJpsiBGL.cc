///
/// @file  FFBctoJpsiBGL.cc
/// @brief \f$ B_c \rightarrow J/\psi \f$ BGL form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>
#include <tuple>

#include "Hammer/FormFactors/BGL/FFBctoJpsiBGL.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFBctoJpsiBGL::FFBctoJpsiBGL() {
        // Create tensor rank and dimensions
        IndexList dims = {8};
        string name{"FFBctoJpsiBGL"};

        setPrefix("BctoJpsi");
        addProcessSignature(PID::BCPLUS, {PID::JPSI});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BCJPSI})});

        setSignatureIndex();
    }

    void FFBctoJpsiBGL::defineSettings() {
        //_mFFErrNames = ;
        setPath(getFFErrPrefixGroup().get());
        setUnits("GeV");

        // Optimized z-expansion
        addSetting<bool>("OptZ", false);

        //        //Using parameters/fit of 1909.10691. Chi values from private correspondence with 1909.10691 authors.
        //        //Uses optimized z expansion
        //        addSetting<double>("Chim", 0.007168); //GeV^-2
        //        addSetting<double>("Chip", 0.012539); //GeV^-2
        //        addSetting<double>("ChimL",0.025042);

        //        vector<double> avec={0.004698, -0.02218, 0.1503};
        //        vector<double> bvec={0.003424, -0.02595, 0.3897};
        //        vector<double> cvec={-0.003164, 0.08731};
        //        vector<double> dvec={0.04011, -0.2132, 0.008157};

        //        vector<double> BcStatesf{6.7305,6.7385,7.1355,7.1435}; //GeV
        //        vector<double> BcStatesg{6.3290,6.8975,7.0065}; //GeV
        //        vector<double> BcStatesP1{6.2749, 6.8710}; //GeV

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

        setInitialized();
    }

    void FFBctoJpsiBGL::evalAtPSPoint(const vector<double>& point, const vector<double>& masses) {
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
        // for(size_t n = 0; n< aF1.size(); ++n) F1 -= aF1[n] * z1minpow[n+1];
        // double F1=af[0]*(Mb-Mc)*phiF1_0/phif_0;
        for (size_t n = 0; n < aF1.size(); ++n) {
            F1 += aF1[n] * (z1pow[n + 1] - z1minpow[n + 1]);
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
        result.element({1}) = f;
        // Fg (dim -1) Definition in BGL (hep-ph/9705252) differs from Manohar, Wise or 1610.02045 by factor of 2
        result.element({2}) = g / 2.;
        // Fm (dim -1)
        result.element({3}) = (Fmf * f + FmF1 * F1 + FmP1 * P1);
        // Fp (dim -1)
        result.element({4}) = (Fpf * f + FpF1 * F1);
        // Fzt (dim -2)
        // result.element({5}) = 0;
        // Fmt (dim 0)
        // result.element({6}) = 0;
        // Fpt (dim 0)
        // result.element({7}) = 0;
    }

    unique_ptr<FormFactorBase> FFBctoJpsiBGL::clone(const string& label) {
        MAKE_CLONE(FFBctoJpsiBGL, label);
    }

    void FFBctoJpsiBGL::addRefs() const {
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
        if (!getSettingsHandler()->checkReference("Harrison:2020gvo")) {
            string ref =
                "@article{Harrison:2020gvo,\n"
                "    author = \"Harrison, Judd and Davies, Christine T. H. and Lytle, Andrew\",\n"
                "    collaboration = \"HPQCD\",\n"
                "    title = \"{$B_c \\rightarrow J/\\psi$ form factors for the full $q^2$ range from lattice QCD}\",\n"
                "    eprint = \"2007.06957\",\n"
                "    archivePrefix = \"arXiv\",\n"
                "    primaryClass = \"hep-lat\",\n"
                "    doi = \"10.1103/PhysRevD.102.094518\",\n"
                "    journal = \"Phys. Rev. D\",\n"
                "    volume = \"102\",\n"
                "    number = \"9\",\n"
                "    pages = \"094518\",\n"
                "    year = \"2020\"\n"
                "}\n";
            getSettingsHandler()->addReference("Harrison:2020gvo", ref);
        }
    }

} // namespace Hammer
