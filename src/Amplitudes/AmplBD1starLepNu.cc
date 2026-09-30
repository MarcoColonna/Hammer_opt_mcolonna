///
/// @file  AmplBD1starLepNu.cc
/// @brief \f$ B \rightarrow D_1^* \tau\nu \f$ amplitude
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/Amplitudes/AmplBD1starLepNu.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;
using namespace complex_literals;

namespace Hammer {

    namespace MD = MultiDimensional;

    AmplBD1starLepNu::AmplBD1starLepNu() {
        // Create tensor rank and dimensions
        IndexList dims{{11, 8, 3, 2, 2, 2}};
        string name{"AmplBD1starLepNu"};
        addProcessSignature(PID::BPLUS, {-PID::DSSD1STAR, PID::NU_TAU, PID::ANTITAU});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCTAUNU, FF_BDSSD1STAR, SPIN_DSSD1STAR, SPIN_NUTAU,
                                                          SPIN_NUTAU_REF, SPIN_TAUP})});

        addProcessSignature(PID::BZERO, {PID::DSSD1STARMINUS, PID::NU_TAU, PID::ANTITAU});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCTAUNU, FF_BDSSD1STAR, SPIN_DSSD1STAR, SPIN_NUTAU,
                                                          SPIN_NUTAU_REF, SPIN_TAUP})});

        addProcessSignature(PID::BPLUS, {-PID::DSSD1STAR, PID::NU_MU, PID::ANTIMUON});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCMUNU, FF_BDSSD1STAR, SPIN_DSSD1STAR, SPIN_NUMU,
                                                          SPIN_NUMU_REF, SPIN_MUP})});

        addProcessSignature(PID::BZERO, {PID::DSSD1STARMINUS, PID::NU_MU, PID::ANTIMUON});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCMUNU, FF_BDSSD1STAR, SPIN_DSSD1STAR, SPIN_NUMU,
                                                          SPIN_NUMU_REF, SPIN_MUP})});

        addProcessSignature(PID::BPLUS, {-PID::DSSD1STAR, PID::NU_E, PID::POSITRON});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCENU, FF_BDSSD1STAR, SPIN_DSSD1STAR, SPIN_NUE,
                                                          SPIN_NUE_REF, SPIN_EP})});

        addProcessSignature(PID::BZERO, {PID::DSSD1STARMINUS, PID::NU_E, PID::POSITRON});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCENU, FF_BDSSD1STAR, SPIN_DSSD1STAR, SPIN_NUE,
                                                          SPIN_NUE_REF, SPIN_EP})});

        // bs -> cs
        addProcessSignature(PID::BS, {PID::DSSDS1STARMINUS, PID::NU_TAU, PID::ANTITAU});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCTAUNU, FF_BSDSSDS1STAR, SPIN_DSSDS1STAR, SPIN_NUTAU,
                                                          SPIN_NUTAU_REF, SPIN_TAUP})});

        addProcessSignature(PID::BS, {PID::DSSDS1STARMINUS, PID::NU_MU, PID::ANTIMUON});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCMUNU, FF_BSDSSDS1STAR, SPIN_DSSDS1STAR, SPIN_NUMU,
                                                          SPIN_NUMU_REF, SPIN_MUP})});

        addProcessSignature(PID::BS, {PID::DSSDS1STARMINUS, PID::NU_E, PID::POSITRON});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCENU, FF_BSDSSDS1STAR, SPIN_DSSDS1STAR, SPIN_NUE,
                                                          SPIN_NUE_REF, SPIN_EP})});

        setSignatureIndex();
    }

    void AmplBD1starLepNu::eval(const Particle& parent, const ParticleList& daughters,
                                const ParticleList& /*references*/) {
        // Momenta
        const FourMomentum& pBmes = parent.momentum();
        const FourMomentum& pD1starmes = daughters[0].momentum();
        const FourMomentum& kNuTau = daughters[1].momentum();
        const FourMomentum& pTau = daughters[2].momentum();

        // kinematic objects
        const double Mb = pBmes.mass();
        const double Mb2 = Mb * Mb;
        const double Mc = pD1starmes.mass();
        const double Mc2 = Mc * Mc;
        const double Mt = pTau.mass();
        const double Mt2 = Mt * Mt;
        const double Sqq = Mb2 + Mc2 - 2. * (pBmes * pD1starmes);
        // const double sqSqq = sqrt(Sqq);
        const double Ew = (Mb2 - Mc2 + Sqq) / (2 * Mb);
        const double BNuTau = (pBmes * kNuTau);
        const double NuTauQ = (pBmes * kNuTau) - (pD1starmes * kNuTau);
        const double BQ = Mb2 - (pBmes * pD1starmes);

        const double w = (Mb2 + Mc2 - Sqq) / (2 * Mb * Mc);
        const double rC = Mc / Mb;
        const double rt = Mt / Mb;

        const double mSqq = Sqq / Mb2;
        const double w2m1 = w * w - 1;
        const double Sqw2m1 = sqrt(w2m1);
        const double SqmSqq = sqrt(mSqq);
        const double Sqw2m1OnmSqq = sqrt(w2m1 / mSqq);
        // const double SqmSqqw2m1 = sqrt(mSqq*w2m1);
        // const double w2m132 = pow(w2m1,1.5);

        // Helicity Angles
        const double CosTt = -((Ew * (Sqq * BNuTau - NuTauQ * BQ)) / (sqrt(Ew * Ew - Sqq) * NuTauQ * BQ));
        const double SinTt = sqrt(1. - CosTt * CosTt);
        const double CosTtHalfSq = (1. + CosTt) / 2.;
        const double SinTtHalfSq = (1. - CosTt) / 2.;

        const double prefactor = 2 * GFermi * sqrt(Mb * Mc) * sqrt(Sqq - Mt2);

        // initialize tensor elements to zero
        Tensor& t = getTensor();
        t.clearData();

        // set non-zero tensor elements
        // NB: Amplitudes include additional D1star phase, defined wrt D1star spin!
        t.element({0, 1, 0, 0, 0, 0}) = (rt * SinTt) / (2. * SqmSqq);
        t.element({0, 4, 0, 0, 0, 0}) = (rt * SinTt * Sqw2m1OnmSqq) / 2.;
        t.element({5, 1, 0, 0, 0, 0}) = (rt * SinTt) / (2. * SqmSqq);
        t.element({5, 4, 0, 0, 0, 0}) = -(rt * SinTt * Sqw2m1OnmSqq) / 2.;
        t.element({7, 1, 0, 0, 0, 0}) = (rt * SinTt) / (2. * SqmSqq);
        t.element({7, 4, 0, 0, 0, 0}) = (rt * SinTt * Sqw2m1OnmSqq) / 2.;
        t.element({9, 5, 0, 0, 0, 0}) = (2. * SinTt * (-1 + rC * (w + Sqw2m1))) / SqmSqq;
        t.element({9, 6, 0, 0, 0, 0}) = (2. * SinTt * (rC - w + Sqw2m1)) / SqmSqq;
        t.element({0, 1, 0, 0, 0, 1}) = CosTtHalfSq;
        t.element({0, 4, 0, 0, 0, 1}) = CosTtHalfSq * Sqw2m1;
        t.element({5, 1, 0, 0, 0, 1}) = CosTtHalfSq;
        t.element({5, 4, 0, 0, 0, 1}) = -(CosTtHalfSq * Sqw2m1);
        t.element({7, 1, 0, 0, 0, 1}) = CosTtHalfSq;
        t.element({7, 4, 0, 0, 0, 1}) = CosTtHalfSq * Sqw2m1;
        t.element({9, 5, 0, 0, 0, 1}) = (4. * CosTtHalfSq * rt * (-1 + rC * (w + Sqw2m1))) / mSqq;
        t.element({9, 6, 0, 0, 0, 1}) = (4. * CosTtHalfSq * rt * (rC - w + Sqw2m1)) / mSqq;
        t.element({6, 1, 0, 1, 1, 0}) = -SinTtHalfSq;
        t.element({6, 4, 0, 1, 1, 0}) = SinTtHalfSq * Sqw2m1;
        t.element({8, 1, 0, 1, 1, 0}) = -SinTtHalfSq;
        t.element({8, 4, 0, 1, 1, 0}) = -(SinTtHalfSq * Sqw2m1);
        t.element({10, 5, 0, 1, 1, 0}) = (4. * rt * SinTtHalfSq * (1 + rC * (-w + Sqw2m1))) / mSqq;
        t.element({10, 6, 0, 1, 1, 0}) = (4. * rt * SinTtHalfSq * (-rC + w + Sqw2m1)) / mSqq;
        t.element({6, 1, 0, 1, 1, 1}) = -(rt * SinTt) / (2. * SqmSqq);
        t.element({6, 4, 0, 1, 1, 1}) = (rt * SinTt * Sqw2m1OnmSqq) / 2.;
        t.element({8, 1, 0, 1, 1, 1}) = -(rt * SinTt) / (2. * SqmSqq);
        t.element({8, 4, 0, 1, 1, 1}) = -(rt * SinTt * Sqw2m1OnmSqq) / 2.;
        t.element({10, 5, 0, 1, 1, 1}) = (2. * SinTt * (1 + rC * (-w + Sqw2m1))) / SqmSqq;
        t.element({10, 6, 0, 1, 1, 1}) = (2. * SinTt * (-rC + w + Sqw2m1)) / SqmSqq;
        t.element({0, 1, 1, 0, 0, 0}) = (rt * (CosTt * (rC - w) + Sqw2m1)) / (sqrt2 * mSqq);
        t.element({0, 2, 1, 0, 0, 0}) = -((rt * ((-1 + rC * w) * Sqw2m1 + CosTt * rC * w2m1)) / (sqrt2 * mSqq));
        t.element({0, 3, 1, 0, 0, 0}) = (rt * ((-rC + w) * Sqw2m1 - CosTt * w2m1)) / (sqrt2 * mSqq);
        t.element({1, 0, 1, 0, 0, 0}) = Sqw2m1 / sqrt2;
        t.element({3, 0, 1, 0, 0, 0}) = Sqw2m1 / sqrt2;
        t.element({5, 1, 1, 0, 0, 0}) = (rt * (CosTt * (rC - w) + Sqw2m1)) / (sqrt2 * mSqq);
        t.element({5, 2, 1, 0, 0, 0}) = -((rt * ((-1 + rC * w) * Sqw2m1 + CosTt * rC * w2m1)) / (sqrt2 * mSqq));
        t.element({5, 3, 1, 0, 0, 0}) = (rt * ((-rC + w) * Sqw2m1 - CosTt * w2m1)) / (sqrt2 * mSqq);
        t.element({7, 1, 1, 0, 0, 0}) = (rt * (CosTt * (rC - w) + Sqw2m1)) / (sqrt2 * mSqq);
        t.element({7, 2, 1, 0, 0, 0}) = -((rt * ((-1 + rC * w) * Sqw2m1 + CosTt * rC * w2m1)) / (sqrt2 * mSqq));
        t.element({7, 3, 1, 0, 0, 0}) = (rt * ((-rC + w) * Sqw2m1 - CosTt * w2m1)) / (sqrt2 * mSqq);
        t.element({9, 5, 1, 0, 0, 0}) = 2. * sqrt2 * CosTt * w;
        t.element({9, 6, 1, 0, 0, 0}) = 2. * sqrt2 * CosTt;
        t.element({9, 7, 1, 0, 0, 0}) = -2. * sqrt2 * CosTt * w2m1;
        t.element({0, 1, 1, 0, 0, 1}) = (SinTt * (-rC + w)) / (sqrt2 * SqmSqq);
        t.element({0, 2, 1, 0, 0, 1}) = (rC * SinTt * w2m1) / (sqrt2 * SqmSqq);
        t.element({0, 3, 1, 0, 0, 1}) = (SinTt * w2m1) / (sqrt2 * SqmSqq);
        t.element({5, 1, 1, 0, 0, 1}) = (SinTt * (-rC + w)) / (sqrt2 * SqmSqq);
        t.element({5, 2, 1, 0, 0, 1}) = (rC * SinTt * w2m1) / (sqrt2 * SqmSqq);
        t.element({5, 3, 1, 0, 0, 1}) = (SinTt * w2m1) / (sqrt2 * SqmSqq);
        t.element({7, 1, 1, 0, 0, 1}) = (SinTt * (-rC + w)) / (sqrt2 * SqmSqq);
        t.element({7, 2, 1, 0, 0, 1}) = (rC * SinTt * w2m1) / (sqrt2 * SqmSqq);
        t.element({7, 3, 1, 0, 0, 1}) = (SinTt * w2m1) / (sqrt2 * SqmSqq);
        t.element({9, 5, 1, 0, 0, 1}) = (-2. * sqrt2 * rt * SinTt * w) / SqmSqq;
        t.element({9, 6, 1, 0, 0, 1}) = (-2. * sqrt2 * rt * SinTt) / SqmSqq;
        t.element({9, 7, 1, 0, 0, 1}) = (2. * sqrt2 * rt * SinTt * w2m1) / SqmSqq;
        t.element({6, 1, 1, 1, 1, 0}) = (SinTt * (-rC + w)) / (sqrt2 * SqmSqq);
        t.element({6, 2, 1, 1, 1, 0}) = (rC * SinTt * w2m1) / (sqrt2 * SqmSqq);
        t.element({6, 3, 1, 1, 1, 0}) = (SinTt * w2m1) / (sqrt2 * SqmSqq);
        t.element({8, 1, 1, 1, 1, 0}) = (SinTt * (-rC + w)) / (sqrt2 * SqmSqq);
        t.element({8, 2, 1, 1, 1, 0}) = (rC * SinTt * w2m1) / (sqrt2 * SqmSqq);
        t.element({8, 3, 1, 1, 1, 0}) = (SinTt * w2m1) / (sqrt2 * SqmSqq);
        t.element({10, 5, 1, 1, 1, 0}) = (-2. * sqrt2 * rt * SinTt * w) / SqmSqq;
        t.element({10, 6, 1, 1, 1, 0}) = (-2. * sqrt2 * rt * SinTt) / SqmSqq;
        t.element({10, 7, 1, 1, 1, 0}) = (2. * sqrt2 * rt * SinTt * w2m1) / SqmSqq;
        t.element({2, 0, 1, 1, 1, 1}) = -(Sqw2m1 / sqrt2);
        t.element({4, 0, 1, 1, 1, 1}) = -(Sqw2m1 / sqrt2);
        t.element({6, 1, 1, 1, 1, 1}) = -((rt * (CosTt * (rC - w) + Sqw2m1)) / (sqrt2 * mSqq));
        t.element({6, 2, 1, 1, 1, 1}) = (rt * ((-1 + rC * w) * Sqw2m1 + CosTt * rC * w2m1)) / (sqrt2 * mSqq);
        t.element({6, 3, 1, 1, 1, 1}) = (rt * ((rC - w) * Sqw2m1 + CosTt * w2m1)) / (sqrt2 * mSqq);
        t.element({8, 1, 1, 1, 1, 1}) = -((rt * (CosTt * (rC - w) + Sqw2m1)) / (sqrt2 * mSqq));
        t.element({8, 2, 1, 1, 1, 1}) = (rt * ((-1 + rC * w) * Sqw2m1 + CosTt * rC * w2m1)) / (sqrt2 * mSqq);
        t.element({8, 3, 1, 1, 1, 1}) = (rt * ((rC - w) * Sqw2m1 + CosTt * w2m1)) / (sqrt2 * mSqq);
        t.element({10, 5, 1, 1, 1, 1}) = -2. * sqrt2 * CosTt * w;
        t.element({10, 6, 1, 1, 1, 1}) = -2. * sqrt2 * CosTt;
        t.element({10, 7, 1, 1, 1, 1}) = 2. * sqrt2 * CosTt * w2m1;
        t.element({0, 1, 2, 0, 0, 0}) = (rt * SinTt) / (2. * SqmSqq);
        t.element({0, 4, 2, 0, 0, 0}) = -(rt * SinTt * Sqw2m1OnmSqq) / 2.;
        t.element({5, 1, 2, 0, 0, 0}) = (rt * SinTt) / (2. * SqmSqq);
        t.element({5, 4, 2, 0, 0, 0}) = (rt * SinTt * Sqw2m1OnmSqq) / 2.;
        t.element({7, 1, 2, 0, 0, 0}) = (rt * SinTt) / (2. * SqmSqq);
        t.element({7, 4, 2, 0, 0, 0}) = -(rt * SinTt * Sqw2m1OnmSqq) / 2.;
        t.element({9, 5, 2, 0, 0, 0}) = (-2. * SinTt * (1 + rC * (-w + Sqw2m1))) / SqmSqq;
        t.element({9, 6, 2, 0, 0, 0}) = (-2. * SinTt * (-rC + w + Sqw2m1)) / SqmSqq;
        t.element({0, 1, 2, 0, 0, 1}) = -SinTtHalfSq;
        t.element({0, 4, 2, 0, 0, 1}) = SinTtHalfSq * Sqw2m1;
        t.element({5, 1, 2, 0, 0, 1}) = -SinTtHalfSq;
        t.element({5, 4, 2, 0, 0, 1}) = -(SinTtHalfSq * Sqw2m1);
        t.element({7, 1, 2, 0, 0, 1}) = -SinTtHalfSq;
        t.element({7, 4, 2, 0, 0, 1}) = SinTtHalfSq * Sqw2m1;
        t.element({9, 5, 2, 0, 0, 1}) = (4. * rt * SinTtHalfSq * (1 + rC * (-w + Sqw2m1))) / mSqq;
        t.element({9, 6, 2, 0, 0, 1}) = (4. * rt * SinTtHalfSq * (-rC + w + Sqw2m1)) / mSqq;
        t.element({6, 1, 2, 1, 1, 0}) = CosTtHalfSq;
        t.element({6, 4, 2, 1, 1, 0}) = CosTtHalfSq * Sqw2m1;
        t.element({8, 1, 2, 1, 1, 0}) = CosTtHalfSq;
        t.element({8, 4, 2, 1, 1, 0}) = -(CosTtHalfSq * Sqw2m1);
        t.element({10, 5, 2, 1, 1, 0}) = (4. * CosTtHalfSq * rt * (-1 + rC * (w + Sqw2m1))) / mSqq;
        t.element({10, 6, 2, 1, 1, 0}) = (4. * CosTtHalfSq * rt * (rC - w + Sqw2m1)) / mSqq;
        t.element({6, 1, 2, 1, 1, 1}) = -(rt * SinTt) / (2. * SqmSqq);
        t.element({6, 4, 2, 1, 1, 1}) = -(rt * SinTt * Sqw2m1OnmSqq) / 2.;
        t.element({8, 1, 2, 1, 1, 1}) = -(rt * SinTt) / (2. * SqmSqq);
        t.element({8, 4, 2, 1, 1, 1}) = (rt * SinTt * Sqw2m1OnmSqq) / 2.;
        t.element({10, 5, 2, 1, 1, 1}) = (-2. * SinTt * (-1 + rC * (w + Sqw2m1))) / SqmSqq;
        t.element({10, 6, 2, 1, 1, 1}) = (-2. * SinTt * (rC - w + Sqw2m1)) / SqmSqq;

        t *= prefactor;
    }

    void AmplBD1starLepNu::addRefs() const {
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
