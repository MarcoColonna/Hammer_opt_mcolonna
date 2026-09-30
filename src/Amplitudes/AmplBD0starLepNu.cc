///
/// @file  AmplBD0starLepNu.cc
/// @brief \f$ B \rightarrow D_0^* \tau\nu \f$ amplitude
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/Amplitudes/AmplBD0starLepNu.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;
using namespace complex_literals;

namespace Hammer {

    namespace MD = MultiDimensional;

    AmplBD0starLepNu::AmplBD0starLepNu() {
        // Create tensor rank and dimensions
        IndexList dims{{11, 4, 2, 2, 2}};
        string name{"AmplBD0starLepNu"};
        addProcessSignature(PID::BPLUS, {-PID::DSSD0STAR, PID::NU_TAU, PID::ANTITAU});
        addTensor(Tensor{
            name, MD::makeEmptySparse(dims, {WILSON_BCTAUNU, FF_BDSSD0STAR, SPIN_NUTAU, SPIN_NUTAU_REF, SPIN_TAUP})});

        addProcessSignature(PID::BZERO, {PID::DSSD0STARMINUS, PID::NU_TAU, PID::ANTITAU});
        addTensor(Tensor{
            name, MD::makeEmptySparse(dims, {WILSON_BCTAUNU, FF_BDSSD0STAR, SPIN_NUTAU, SPIN_NUTAU_REF, SPIN_TAUP})});

        addProcessSignature(PID::BPLUS, {-PID::DSSD0STAR, PID::NU_MU, PID::ANTIMUON});
        addTensor(Tensor{
            name, MD::makeEmptySparse(dims, {WILSON_BCMUNU, FF_BDSSD0STAR, SPIN_NUMU, SPIN_NUMU_REF, SPIN_MUP})});

        addProcessSignature(PID::BZERO, {PID::DSSD0STARMINUS, PID::NU_MU, PID::ANTIMUON});
        addTensor(Tensor{
            name, MD::makeEmptySparse(dims, {WILSON_BCMUNU, FF_BDSSD0STAR, SPIN_NUMU, SPIN_NUMU_REF, SPIN_MUP})});

        addProcessSignature(PID::BPLUS, {-PID::DSSD0STAR, PID::NU_E, PID::POSITRON});
        addTensor(
            Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCENU, FF_BDSSD0STAR, SPIN_NUE, SPIN_NUE_REF, SPIN_EP})});

        addProcessSignature(PID::BZERO, {PID::DSSD0STARMINUS, PID::NU_E, PID::POSITRON});
        addTensor(
            Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCENU, FF_BDSSD0STAR, SPIN_NUE, SPIN_NUE_REF, SPIN_EP})});

        // bs -> cs
        addProcessSignature(PID::BS, {PID::DSSDS0STARMINUS, PID::NU_TAU, PID::ANTITAU});
        addTensor(Tensor{
            name, MD::makeEmptySparse(dims, {WILSON_BCTAUNU, FF_BSDSSDS0STAR, SPIN_NUTAU, SPIN_NUTAU_REF, SPIN_TAUP})});

        addProcessSignature(PID::BS, {PID::DSSDS0STARMINUS, PID::NU_MU, PID::ANTIMUON});
        addTensor(Tensor{
            name, MD::makeEmptySparse(dims, {WILSON_BCMUNU, FF_BSDSSDS0STAR, SPIN_NUMU, SPIN_NUMU_REF, SPIN_MUP})});

        addProcessSignature(PID::BS, {PID::DSSDS0STARMINUS, PID::NU_E, PID::POSITRON});
        addTensor(
            Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCENU, FF_BSDSSDS0STAR, SPIN_NUE, SPIN_NUE_REF, SPIN_EP})});

        setSignatureIndex();
    }

    void AmplBD0starLepNu::eval(const Particle& parent, const ParticleList& daughters,
                                const ParticleList& /*references*/) {
        // Momenta
        const FourMomentum& pBmes = parent.momentum();
        const FourMomentum& pD0starmes = daughters[0].momentum();
        const FourMomentum& kNuTau = daughters[1].momentum();
        const FourMomentum& pTau = daughters[2].momentum();

        // kinematic objects
        const double Mb = pBmes.mass();
        const double Mb2 = Mb * Mb;
        const double Mc = pD0starmes.mass();
        const double Mc2 = Mc * Mc;
        const double Mt = pTau.mass();
        const double Mt2 = Mt * Mt;
        const double Sqq = Mb2 + Mc2 - 2. * (pBmes * pD0starmes);
        // const double sqSqq = sqrt(Sqq);
        const double Ew = (Mb2 - Mc2 + Sqq) / (2 * Mb);
        const double BNuTau = (pBmes * kNuTau);
        const double NuTauQ = (pBmes * kNuTau) - (pD0starmes * kNuTau);
        const double BQ = Mb2 - (pBmes * pD0starmes);

        const double w = (Mb2 + Mc2 - Sqq) / (2 * Mb * Mc);
        const double rC = Mc / Mb;
        const double rt = Mt / Mb;

        const double mSqq = Sqq / Mb2;
        const double w2m1 = w * w - 1;
        const double Sqw2m1 = sqrt(w2m1);
        // const double SqmSqq = sqrt(mSqq);
        const double Sqw2m1OnmSqq = sqrt(w2m1 / mSqq);
        const double SqmSqqw2m1 = sqrt(mSqq * w2m1);
        // const double w2m132 = pow(w2m1,1.5);

        // Helicity Angles
        const double CosTt = -((Ew * (Sqq * BNuTau - NuTauQ * BQ)) / (sqrt(Ew * Ew - Sqq) * NuTauQ * BQ));
        const double SinTt = sqrt(1. - CosTt * CosTt);
        // const double CosTtHalfSq = (1. + CosTt) / 2.;
        // const double SinTtHalfSq = (1. - CosTt) / 2.;

        const double prefactor = 2 * sqrt2 * GFermi * sqrt(Mb * Mc) * sqrt(Sqq - Mt2);

        // initialize tensor elements to zero
        Tensor& t = getTensor();
        t.clearData();

        // set non-zero tensor elements
        // NB: Amplitudes include additional D0star phase, defined wrt D0star spin!
        t.element({0, 1, 0, 0, 0}) = -(rt * ((rC - 1) * (1 + w) + CosTt * (1 + rC) * Sqw2m1)) / (2. * mSqq);
        t.element({0, 2, 0, 0, 0}) = -(rt * ((1 + rC) * (w - 1) + CosTt * (rC - 1) * Sqw2m1)) / (2. * mSqq);
        t.element({1, 0, 0, 0, 0}) = -0.5;
        t.element({3, 0, 0, 0, 0}) = 0.5;
        t.element({5, 1, 0, 0, 0}) = (rt * ((rC - 1) * (1 + w) + CosTt * (1 + rC) * Sqw2m1)) / (2. * mSqq);
        t.element({5, 2, 0, 0, 0}) = (rt * ((1 + rC) * (w - 1) + CosTt * (rC - 1) * Sqw2m1)) / (2. * mSqq);
        t.element({7, 1, 0, 0, 0}) = -(rt * ((rC - 1) * (1 + w) + CosTt * (1 + rC) * Sqw2m1)) / (2. * mSqq);
        t.element({7, 2, 0, 0, 0}) = -(rt * ((1 + rC) * (w - 1) + CosTt * (rC - 1) * Sqw2m1)) / (2. * mSqq);
        t.element({9, 3, 0, 0, 0}) = 2. * CosTt * Sqw2m1;
        t.element({0, 1, 0, 0, 1}) = ((1 + rC) * SinTt * Sqw2m1OnmSqq) / 2.;
        t.element({0, 2, 0, 0, 1}) = ((rC - 1) * SinTt * Sqw2m1OnmSqq) / 2.;
        t.element({5, 1, 0, 0, 1}) = -((1 + rC) * SinTt * Sqw2m1OnmSqq) / 2.;
        t.element({5, 2, 0, 0, 1}) = -((rC - 1) * SinTt * Sqw2m1OnmSqq) / 2.;
        t.element({7, 1, 0, 0, 1}) = ((1 + rC) * SinTt * Sqw2m1OnmSqq) / 2.;
        t.element({7, 2, 0, 0, 1}) = ((rC - 1) * SinTt * Sqw2m1OnmSqq) / 2.;
        t.element({9, 3, 0, 0, 1}) = -2. * rt * SinTt * Sqw2m1OnmSqq;
        t.element({6, 1, 1, 1, 0}) = -((1 + rC) * SinTt * SqmSqqw2m1) / (2. * mSqq);
        t.element({6, 2, 1, 1, 0}) = -((rC - 1) * SinTt * SqmSqqw2m1) / (2. * mSqq);
        t.element({8, 1, 1, 1, 0}) = ((1 + rC) * SinTt * SqmSqqw2m1) / (2. * mSqq);
        t.element({8, 2, 1, 1, 0}) = ((rC - 1) * SinTt * SqmSqqw2m1) / (2. * mSqq);
        t.element({10, 3, 1, 1, 0}) = 2. * rt * SinTt * Sqw2m1OnmSqq;
        t.element({2, 0, 1, 1, 1}) = 0.5;
        t.element({4, 0, 1, 1, 1}) = -0.5;
        t.element({6, 1, 1, 1, 1}) = -(rt * ((rC - 1) * (1 + w) + CosTt * (1 + rC) * Sqw2m1)) / (2. * mSqq);
        t.element({6, 2, 1, 1, 1}) = -(rt * ((1 + rC) * (w - 1) + CosTt * (rC - 1) * Sqw2m1)) / (2. * mSqq);
        t.element({8, 1, 1, 1, 1}) = (rt * ((rC - 1) * (1 + w) + CosTt * (1 + rC) * Sqw2m1)) / (2. * mSqq);
        t.element({8, 2, 1, 1, 1}) = (rt * ((1 + rC) * (w - 1) + CosTt * (rC - 1) * Sqw2m1)) / (2. * mSqq);
        t.element({10, 3, 1, 1, 1}) = 2. * CosTt * Sqw2m1;

        t *= prefactor;
    }

    void AmplBD0starLepNu::addRefs() const {
        if (!getSettingsHandler()->checkReference("Bernlochner:2017jxt")) {
            string ref =
                "@article{Bernlochner:2017jxt,\n"
                "     author         = \"Bernlochner, Florian U. and Ligeti, Zoltan and Robinson, Dean J.\",\n"
                "     title          = \"{Model independent analysis of semileptonic $B$ decays to $D^{**}$ for "
                "arbitrary new physics}\",\n"
                "     journal        = \"Phys. Rev.\",\n"
                "     volume         = \"D97\",\n"
                "     year           = \"2018\",\n"
                "     number         = \"7\",\n"
                "     pages          = \"075011\",\n"
                "     doi            = \"10.1103/PhysRevD.97.075011\",\n"
                "     eprint         = \"1711.03110\",\n"
                "     archivePrefix  = \"arXiv\",\n"
                "     primaryClass   = \"hep-ph\",\n"
                "     SLACcitation   = \"%%CITATION = ARXIV:1711.03110;%%\"\n"
                "}\n";
            getSettingsHandler()->addReference("Bernlochner:2017jxt", ref);
        }
    }

} // namespace Hammer
