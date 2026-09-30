///
/// @file  AmplBDstarLepNu.cc
/// @brief \f$ B \rightarrow D^* \tau\nu \f$ amplitude
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/Amplitudes/AmplBDstarLepNu.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;
using namespace complex_literals;

namespace Hammer {

    namespace MD = MultiDimensional;

    AmplBDstarLepNu::AmplBDstarLepNu() {
        // Create tensor rank and dimensions
        IndexList dims{{11, 8, 3, 2, 2, 2}};
        string name{"AmplBDstarLepNu"};
        addProcessSignature(PID::BPLUS, {-PID::DSTAR, PID::NU_TAU, PID::ANTITAU});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCTAUNU, FF_BDSTAR, SPIN_DSTAR, SPIN_NUTAU,
                                                          SPIN_NUTAU_REF, SPIN_TAUP})});

        addProcessSignature(PID::BZERO, {PID::DSTARMINUS, PID::NU_TAU, PID::ANTITAU});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCTAUNU, FF_BDSTAR, SPIN_DSTAR, SPIN_NUTAU,
                                                          SPIN_NUTAU_REF, SPIN_TAUP})});

        addProcessSignature(PID::BPLUS, {-PID::DSTAR, PID::NU_MU, PID::ANTIMUON});
        addTensor(Tensor{name, MD::makeEmptySparse(
                                   dims, {WILSON_BCMUNU, FF_BDSTAR, SPIN_DSTAR, SPIN_NUMU, SPIN_NUMU_REF, SPIN_MUP})});

        addProcessSignature(PID::BZERO, {PID::DSTARMINUS, PID::NU_MU, PID::ANTIMUON});
        addTensor(Tensor{name, MD::makeEmptySparse(
                                   dims, {WILSON_BCMUNU, FF_BDSTAR, SPIN_DSTAR, SPIN_NUMU, SPIN_NUMU_REF, SPIN_MUP})});

        addProcessSignature(PID::BPLUS, {-PID::DSTAR, PID::NU_E, PID::POSITRON});
        addTensor(Tensor{
            name, MD::makeEmptySparse(dims, {WILSON_BCENU, FF_BDSTAR, SPIN_DSTAR, SPIN_NUE, SPIN_NUE_REF, SPIN_EP})});

        addProcessSignature(PID::BZERO, {PID::DSTARMINUS, PID::NU_E, PID::POSITRON});
        addTensor(Tensor{
            name, MD::makeEmptySparse(dims, {WILSON_BCENU, FF_BDSTAR, SPIN_DSTAR, SPIN_NUE, SPIN_NUE_REF, SPIN_EP})});

        // bs -> cs
        addProcessSignature(PID::BS, {PID::DSSTARMINUS, PID::NU_TAU, PID::ANTITAU});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCTAUNU, FF_BSDSSTAR, SPIN_DSSTAR, SPIN_NUTAU,
                                                          SPIN_NUTAU_REF, SPIN_TAUP})});

        addProcessSignature(PID::BS, {PID::DSSTARMINUS, PID::NU_MU, PID::ANTIMUON});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCMUNU, FF_BSDSSTAR, SPIN_DSSTAR, SPIN_NUMU,
                                                          SPIN_NUMU_REF, SPIN_MUP})});

        addProcessSignature(PID::BS, {PID::DSSTARMINUS, PID::NU_E, PID::POSITRON});
        addTensor(Tensor{name, MD::makeEmptySparse(
                                   dims, {WILSON_BCENU, FF_BSDSSTAR, SPIN_DSSTAR, SPIN_NUE, SPIN_NUE_REF, SPIN_EP})});

        // bc -> cc
        addProcessSignature(PID::BCPLUS, {PID::JPSI, PID::NU_TAU, PID::ANTITAU});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCTAUNU, FF_BCJPSI, SPIN_JPSI, SPIN_NUTAU,
                                                          SPIN_NUTAU_REF, SPIN_TAUP})});

        addProcessSignature(PID::BCPLUS, {PID::JPSI, PID::NU_MU, PID::ANTIMUON});
        addTensor(Tensor{name, MD::makeEmptySparse(
                                   dims, {WILSON_BCMUNU, FF_BCJPSI, SPIN_JPSI, SPIN_NUMU, SPIN_NUMU_REF, SPIN_MUP})});

        addProcessSignature(PID::BCPLUS, {PID::JPSI, PID::NU_E, PID::POSITRON});
        addTensor(Tensor{
            name, MD::makeEmptySparse(dims, {WILSON_BCENU, FF_BCJPSI, SPIN_JPSI, SPIN_NUE, SPIN_NUE_REF, SPIN_EP})});
        setSignatureIndex();
    }

    void AmplBDstarLepNu::eval(const Particle& parent, const ParticleList& daughters,
                               const ParticleList& /*references*/) {
        // Momenta
        const FourMomentum& pBmes = parent.momentum();
        const FourMomentum& pDstarmes = daughters[0].momentum();
        const FourMomentum& kNuTau = daughters[1].momentum();
        const FourMomentum& pTau = daughters[2].momentum();

        // kinematic objects
        const double Mb = pBmes.mass();
        const double Mb2 = Mb * Mb;
        const double Mds = pDstarmes.mass();
        const double Mds2 = Mds * Mds;
        const double Mt = pTau.mass();
        const double Sqq = Mb2 + Mds2 - 2. * (pBmes * pDstarmes);
        const double sqSqq = sqrt(Sqq);
        const double Ew = (Mb2 - Mds2 + Sqq) / (2 * Mb);
        // const double Eds = (Mb2 + Mds2 - Sqq) / (2 * Mb);
        const double Pw = sqrt(Ew * Ew - Sqq);
        const double Pw2 = Pw * Pw;
        const double BNuTau = (pBmes * kNuTau);
        const double NuTauQ = (pBmes * kNuTau) - (pDstarmes * kNuTau);
        const double BQ = Mb2 - (pBmes * pDstarmes);

        // Helicity Angles
        const double CosTt = -((Ew * (Sqq * BNuTau - NuTauQ * BQ)) / (sqrt(Ew * Ew - Sqq) * NuTauQ * BQ));
        const double SinTt = sqrt(1. - CosTt * CosTt);
        const double CosTtHalfSq = (1. + CosTt) / 2.;
        const double SinTtHalfSq = (1. - CosTt) / 2.;

        const double TwoGfSq = 2. * GFermi * sqrt(Sqq - Mt * Mt);

        // initialize tensor elements to zero
        Tensor& t = getTensor();
        t.clearData();

        // set non-zero tensor elements
        // NB: Amplitudes include additional Dstar phase, defined wrt Dstar spin, that removes tau phase. Must be
        // include in D* decay amplitudes.
        t.element({0, 1, 0, 0, 0, 0}) = -(Mt * SinTt) / (2. * sqSqq);
        t.element({0, 2, 0, 0, 0, 0}) = -((Mb * Mt * Pw * SinTt) / sqSqq);
        t.element({5, 1, 0, 0, 0, 0}) = (Mt * SinTt) / (2. * sqSqq);
        t.element({5, 2, 0, 0, 0, 0}) = -((Mb * Mt * Pw * SinTt) / sqSqq);
        t.element({7, 1, 0, 0, 0, 0}) = -(Mt * SinTt) / (2. * sqSqq);
        t.element({7, 2, 0, 0, 0, 0}) = -((Mb * Mt * Pw * SinTt) / sqSqq);
        t.element({9, 6, 0, 0, 0, 0}) = 2. * sqSqq * SinTt;
        t.element({9, 7, 0, 0, 0, 0}) = (2. * (Mb2 - Mds2 - 2. * Mb * Pw) * SinTt) / sqSqq;
        t.element({0, 1, 0, 0, 0, 1}) = -CosTtHalfSq;
        t.element({0, 2, 0, 0, 0, 1}) = -2. * Mb * Pw * CosTtHalfSq;
        t.element({5, 1, 0, 0, 0, 1}) = CosTtHalfSq;
        t.element({5, 2, 0, 0, 0, 1}) = -2. * Mb * Pw * CosTtHalfSq;
        t.element({7, 1, 0, 0, 0, 1}) = -CosTtHalfSq;
        t.element({7, 2, 0, 0, 0, 1}) = -2. * Mb * Pw * CosTtHalfSq;
        t.element({9, 6, 0, 0, 0, 1}) = 4. * Mt * CosTtHalfSq;
        t.element({9, 7, 0, 0, 0, 1}) = (4. * Mt * (Mb2 - Mds2 - 2. * Mb * Pw) * CosTtHalfSq) / Sqq;
        t.element({6, 1, 0, 1, 1, 0}) = -SinTtHalfSq;
        t.element({6, 2, 0, 1, 1, 0}) = 2. * Mb * Pw * SinTtHalfSq;
        t.element({8, 1, 0, 1, 1, 0}) = SinTtHalfSq;
        t.element({8, 2, 0, 1, 1, 0}) = 2. * Mb * Pw * SinTtHalfSq;
        t.element({10, 6, 0, 1, 1, 0}) = 4. * Mt * SinTtHalfSq;
        t.element({10, 7, 0, 1, 1, 0}) = (4. * Mt * (Mb2 - Mds2 + 2. * Mb * Pw) * SinTtHalfSq) / Sqq;
        t.element({6, 1, 0, 1, 1, 1}) = -(Mt * SinTt) / (2. * sqSqq);
        t.element({6, 2, 0, 1, 1, 1}) = (Mb * Mt * Pw * SinTt) / sqSqq;
        t.element({8, 1, 0, 1, 1, 1}) = (Mt * SinTt) / (2. * sqSqq);
        t.element({8, 2, 0, 1, 1, 1}) = (Mb * Mt * Pw * SinTt) / sqSqq;
        t.element({10, 6, 0, 1, 1, 1}) = 2. * sqSqq * SinTt;
        t.element({10, 7, 0, 1, 1, 1}) = (2. * (Mb2 - Mds2 + 2. * Mb * Pw) * SinTt) / sqSqq;
        t.element({0, 1, 1, 0, 0, 0}) = -(Mt * (2. * Mb * Pw + (-Mb2 + Mds2 + Sqq) * CosTt)) / (2. * sqrt2 * Mds * Sqq);
        t.element({0, 3, 1, 0, 0, 0}) = -((Mb * Mt * Pw) / (sqrt2 * Mds));
        t.element({0, 4, 1, 0, 0, 0}) = (Mb * Mt * Pw * (-Mb2 + Mds2 + 2. * Mb * Pw * CosTt)) / (sqrt2 * Mds * Sqq);
        t.element({1, 0, 1, 0, 0, 0}) = (Mb * Pw) / (sqrt2 * Mds);
        t.element({3, 0, 1, 0, 0, 0}) = -((Mb * Pw) / (sqrt2 * Mds));
        t.element({5, 1, 1, 0, 0, 0}) = (Mt * (2. * Mb * Pw + (-Mb2 + Mds2 + Sqq) * CosTt)) / (2. * sqrt2 * Mds * Sqq);
        t.element({5, 3, 1, 0, 0, 0}) = (Mb * Mt * Pw) / (sqrt2 * Mds);
        t.element({5, 4, 1, 0, 0, 0}) = -((Mb * Mt * Pw * (-Mb2 + Mds2 + 2. * Mb * Pw * CosTt)) / (sqrt2 * Mds * Sqq));
        t.element({7, 1, 1, 0, 0, 0}) = -(Mt * (2. * Mb * Pw + (-Mb2 + Mds2 + Sqq) * CosTt)) / (2. * sqrt2 * Mds * Sqq);
        t.element({7, 3, 1, 0, 0, 0}) = -((Mb * Mt * Pw) / (sqrt2 * Mds));
        t.element({7, 4, 1, 0, 0, 0}) = (Mb * Mt * Pw * (-Mb2 + Mds2 + 2. * Mb * Pw * CosTt)) / (sqrt2 * Mds * Sqq);
        t.element({9, 5, 1, 0, 0, 0}) = (-4. * sqrt2 * Mb2 * Pw2 * CosTt) / Mds;
        t.element({9, 6, 1, 0, 0, 0}) = (sqrt2 * (-Mb2 + Mds2 + Sqq) * CosTt) / Mds;
        t.element({9, 7, 1, 0, 0, 0}) = (sqrt2 * (-Mb2 - 3. * Mds2 + Sqq) * CosTt) / Mds;
        t.element({0, 1, 1, 0, 0, 1}) = ((-Mb2 + Mds2 + Sqq) * SinTt) / (2. * sqrt2 * Mds * sqSqq);
        t.element({0, 4, 1, 0, 0, 1}) = -((sqrt2 * Mb2 * Pw2 * SinTt) / (Mds * sqSqq));
        t.element({5, 1, 1, 0, 0, 1}) = -((-Mb2 + Mds2 + Sqq) * SinTt) / (2. * sqrt2 * Mds * sqSqq);
        t.element({5, 4, 1, 0, 0, 1}) = (sqrt2 * Mb2 * Pw2 * SinTt) / (Mds * sqSqq);
        t.element({7, 1, 1, 0, 0, 1}) = ((-Mb2 + Mds2 + Sqq) * SinTt) / (2. * sqrt2 * Mds * sqSqq);
        t.element({7, 4, 1, 0, 0, 1}) = -((sqrt2 * Mb2 * Pw2 * SinTt) / (Mds * sqSqq));
        t.element({9, 5, 1, 0, 0, 1}) = (4. * sqrt2 * Mb2 * Mt * Pw2 * SinTt) / (Mds * sqSqq);
        t.element({9, 6, 1, 0, 0, 1}) = -((sqrt2 * Mt * (-Mb2 + Mds2 + Sqq) * SinTt) / (Mds * sqSqq));
        t.element({9, 7, 1, 0, 0, 1}) = (sqrt2 * Mt * (Mb2 + 3. * Mds2 - Sqq) * SinTt) / (Mds * sqSqq);
        t.element({6, 1, 1, 1, 1, 0}) = -((-Mb2 + Mds2 + Sqq) * SinTt) / (2. * sqrt2 * Mds * sqSqq);
        t.element({6, 4, 1, 1, 1, 0}) = (sqrt2 * Mb2 * Pw2 * SinTt) / (Mds * sqSqq);
        t.element({8, 1, 1, 1, 1, 0}) = ((-Mb2 + Mds2 + Sqq) * SinTt) / (2. * sqrt2 * Mds * sqSqq);
        t.element({8, 4, 1, 1, 1, 0}) = -((sqrt2 * Mb2 * Pw2 * SinTt) / (Mds * sqSqq));
        t.element({10, 5, 1, 1, 1, 0}) = (-4. * sqrt2 * Mb2 * Mt * Pw2 * SinTt) / (Mds * sqSqq);
        t.element({10, 6, 1, 1, 1, 0}) = (sqrt2 * Mt * (-Mb2 + Mds2 + Sqq) * SinTt) / (Mds * sqSqq);
        t.element({10, 7, 1, 1, 1, 0}) = (sqrt2 * Mt * (-Mb2 - 3. * Mds2 + Sqq) * SinTt) / (Mds * sqSqq);
        t.element({2, 0, 1, 1, 1, 1}) = -((Mb * Pw) / (sqrt2 * Mds));
        t.element({4, 0, 1, 1, 1, 1}) = (Mb * Pw) / (sqrt2 * Mds);
        t.element({6, 1, 1, 1, 1, 1}) = -(Mt * (2. * Mb * Pw + (-Mb2 + Mds2 + Sqq) * CosTt)) / (2. * sqrt2 * Mds * Sqq);
        t.element({6, 3, 1, 1, 1, 1}) = -((Mb * Mt * Pw) / (sqrt2 * Mds));
        t.element({6, 4, 1, 1, 1, 1}) = (Mb * Mt * Pw * (-Mb2 + Mds2 + 2. * Mb * Pw * CosTt)) / (sqrt2 * Mds * Sqq);
        t.element({8, 1, 1, 1, 1, 1}) = (Mt * (2. * Mb * Pw + (-Mb2 + Mds2 + Sqq) * CosTt)) / (2. * sqrt2 * Mds * Sqq);
        t.element({8, 3, 1, 1, 1, 1}) = (Mb * Mt * Pw) / (sqrt2 * Mds);
        t.element({8, 4, 1, 1, 1, 1}) = -((Mb * Mt * Pw * (-Mb2 + Mds2 + 2. * Mb * Pw * CosTt)) / (sqrt2 * Mds * Sqq));
        t.element({10, 5, 1, 1, 1, 1}) = (-4. * sqrt2 * Mb2 * Pw2 * CosTt) / Mds;
        t.element({10, 6, 1, 1, 1, 1}) = (sqrt2 * (-Mb2 + Mds2 + Sqq) * CosTt) / Mds;
        t.element({10, 7, 1, 1, 1, 1}) = (sqrt2 * (-Mb2 - 3. * Mds2 + Sqq) * CosTt) / Mds;
        t.element({0, 1, 2, 0, 0, 0}) = -(Mt * SinTt) / (2. * sqSqq);
        t.element({0, 2, 2, 0, 0, 0}) = (Mb * Mt * Pw * SinTt) / sqSqq;
        t.element({5, 1, 2, 0, 0, 0}) = (Mt * SinTt) / (2. * sqSqq);
        t.element({5, 2, 2, 0, 0, 0}) = (Mb * Mt * Pw * SinTt) / sqSqq;
        t.element({7, 1, 2, 0, 0, 0}) = -(Mt * SinTt) / (2. * sqSqq);
        t.element({7, 2, 2, 0, 0, 0}) = (Mb * Mt * Pw * SinTt) / sqSqq;
        t.element({9, 6, 2, 0, 0, 0}) = 2. * sqSqq * SinTt;
        t.element({9, 7, 2, 0, 0, 0}) = (2. * (Mb2 - Mds2 + 2. * Mb * Pw) * SinTt) / sqSqq;
        t.element({0, 1, 2, 0, 0, 1}) = SinTtHalfSq;
        t.element({0, 2, 2, 0, 0, 1}) = -2. * Mb * Pw * SinTtHalfSq;
        t.element({5, 1, 2, 0, 0, 1}) = -SinTtHalfSq;
        t.element({5, 2, 2, 0, 0, 1}) = -2. * Mb * Pw * SinTtHalfSq;
        t.element({7, 1, 2, 0, 0, 1}) = SinTtHalfSq;
        t.element({7, 2, 2, 0, 0, 1}) = -2. * Mb * Pw * SinTtHalfSq;
        t.element({9, 6, 2, 0, 0, 1}) = -4. * Mt * SinTtHalfSq;
        t.element({9, 7, 2, 0, 0, 1}) = (4. * Mt * (Mds2 - Mb * (Mb + 2. * Pw)) * SinTtHalfSq) / Sqq;
        t.element({6, 1, 2, 1, 1, 0}) = CosTtHalfSq;
        t.element({6, 2, 2, 1, 1, 0}) = 2. * Mb * Pw * CosTtHalfSq;
        t.element({8, 1, 2, 1, 1, 0}) = -CosTtHalfSq;
        t.element({8, 2, 2, 1, 1, 0}) = 2. * Mb * Pw * CosTtHalfSq;
        t.element({10, 6, 2, 1, 1, 0}) = -4. * Mt * CosTtHalfSq;
        t.element({10, 7, 2, 1, 1, 0}) = (4. * Mt * (-Mb2 + Mds2 + 2. * Mb * Pw) * CosTtHalfSq) / Sqq;
        t.element({6, 1, 2, 1, 1, 1}) = -(Mt * SinTt) / (2. * sqSqq);
        t.element({6, 2, 2, 1, 1, 1}) = -((Mb * Mt * Pw * SinTt) / sqSqq);
        t.element({8, 1, 2, 1, 1, 1}) = (Mt * SinTt) / (2. * sqSqq);
        t.element({8, 2, 2, 1, 1, 1}) = -((Mb * Mt * Pw * SinTt) / sqSqq);
        t.element({10, 6, 2, 1, 1, 1}) = 2. * sqSqq * SinTt;
        t.element({10, 7, 2, 1, 1, 1}) = (2. * (Mb2 - Mds2 - 2. * Mb * Pw) * SinTt) / sqSqq;

        t *= TwoGfSq;
    }

    void AmplBDstarLepNu::addRefs() const {
        if (!getSettingsHandler()->checkReference("Ligeti:2016npd")) {
            string ref =
                "@article{Ligeti:2016npd,\n"
                "      author         = \"Ligeti, Zoltan and Papucci, Michele and Robinson, Dean J.\",\n"
                "      title          = \"{New Physics in the Visible Final States of $B\\to D^{(*)}\\tau\\nu$}\",\n"
                "      journal        = \"JHEP\",\n"
                "      volume         = \"01\",\n"
                "      year           = \"2017\",\n"
                "      pages          = \"083\",\n"
                "      doi            = \"10.1007/JHEP01(2017)083\",\n"
                "      eprint         = \"1610.02045\",\n"
                "      archivePrefix  = \"arXiv\",\n"
                "      primaryClass   = \"hep-ph\",\n"
                "      SLACcitation   = \"%%CITATION = ARXIV:1610.02045;%%\"\n"
                "}\n";
            getSettingsHandler()->addReference("Ligeti:2016npd", ref);
        }
    }


} // namespace Hammer
