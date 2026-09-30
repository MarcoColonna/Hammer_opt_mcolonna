///
/// @file  AmplLbLcLepNu.cc
/// @brief \f$ \Lambda_b \rightarrow \Lambda_c \tau\nu \f$ amplitude
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/Amplitudes/AmplLbLcLepNu.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;
using namespace complex_literals;

namespace Hammer {

    namespace MD = MultiDimensional;

    AmplLbLcLepNu::AmplLbLcLepNu() {
        // Create tensor rank and dimensions
        IndexList dims{{11, 12, 2, 2, 2, 2, 2}};
        string name{"AmplLbLcLepNu"};
        addProcessSignature(-PID::LAMBDAB, {PID::LAMBDACMINUS, PID::NU_TAU, PID::ANTITAU});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCTAUNU, FF_LBLC, SPIN_LB, SPIN_LC, SPIN_NUTAU,
                                                          SPIN_NUTAU_REF, SPIN_TAUP})});

        addProcessSignature(-PID::LAMBDAB, {PID::LAMBDACMINUS, PID::NU_MU, PID::ANTIMUON});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCMUNU, FF_LBLC, SPIN_LB, SPIN_LC, SPIN_NUMU,
                                                          SPIN_NUMU_REF, SPIN_MUP})});

        addProcessSignature(-PID::LAMBDAB, {PID::LAMBDACMINUS, PID::NU_E, PID::POSITRON});
        addTensor(Tensor{name, MD::makeEmptySparse(
                                   dims, {WILSON_BCENU, FF_LBLC, SPIN_LB, SPIN_LC, SPIN_NUE, SPIN_NUE_REF, SPIN_EP})});

        setSignatureIndex();
        _multiplicity = 2ul;
    }

    void AmplLbLcLepNu::eval(const Particle& parent, const ParticleList& daughters,
                             const ParticleList& /*references*/) {
        // Momenta
        const FourMomentum& pLbmes = parent.momentum();
        const FourMomentum& pLcmes = daughters[0].momentum();
        const FourMomentum& kNuTau = daughters[1].momentum();
        const FourMomentum& pTau = daughters[2].momentum();

        // kinematic objects
        const double Mb = pLbmes.mass();
        const double Mb2 = Mb * Mb;
        const double Mc = pLcmes.mass();
        const double Mc2 = Mc * Mc;
        const double Mt = pTau.mass();
        const double Sqq = Mb2 + Mc2 - 2. * (pLbmes * pLcmes);
        const double Ew = (Mb2 - Mc2 + Sqq) / (2 * Mb);
        const double BNuTau = (pLbmes * kNuTau);
        const double NuTauQ = (pLbmes * kNuTau) - (pLcmes * kNuTau);
        const double BQ = Mb2 - (pLbmes * pLcmes);

        const double mSqq = Sqq / Mb2;
        const double sqmSqq = sqrt(mSqq);
        const double w = (Mb2 + Mc2 - Sqq) / (2. * Mb * Mc);
        const double rC = Mc / Mb;
        const double rt = Mt / Mb;

        // Helicity Angles
        const double CosTt = -((Ew * (Sqq * BNuTau - NuTauQ * BQ)) / (sqrt(Ew * Ew - Sqq) * NuTauQ * BQ));
        const double SinTt = sqrt(1. - CosTt * CosTt);
        const double CosTtHalfSq = (1. + CosTt) / 2.;
        const double SinTtHalfSq = (1. - CosTt) / 2.;
        // double CscTt = 1. / SinTt;

        // Collective Functions
        const double wp = w + sqrt(w * w - 1);
        const double wm = w - sqrt(w * w - 1);
        const double sqwp = sqrt(wp);
        const double sqwm = sqrt(wm);
        const double WpmP = sqwp + sqwm;
        const double WpmM = sqwp - sqwm;
        const double RCmp = 1 - rC + (1 + rC) * CosTt;
        const double RCmm = 1 - rC - (1 + rC) * CosTt;
        const double RCpm = 1 + rC - (1 - rC) * CosTt;
        const double RCpp = 1 + rC + (1 - rC) * CosTt;
        const double Op = rC - w + sqrt(w * w - 1) * CosTt;
        const double Ot = rC * w - 1 + rC * sqrt(w * w - 1) * CosTt;

        const double TwoSqTwoGfSq = TwoSqTwoGFermi * Mb2 * sqrt(mSqq - pow(rt, 2.)) * sqrt(rC);

        // initialize tensor elements to zero
        Tensor& t = getTensor();
        t.clearData();

        // set non-zero tensor elements
        t.element({0, 2, 0, 0, 0, 0, 0}) = (rt * (RCmp * sqwm + RCmm * sqwp)) / (2. * mSqq);
        t.element({0, 3, 0, 0, 0, 0, 0}) = -(Ot * rt * WpmP) / (2. * mSqq);
        t.element({0, 4, 0, 0, 0, 0, 0}) = -(Op * rt * WpmP) / (2. * mSqq);
        t.element({0, 5, 0, 0, 0, 0, 0}) = (rt * (RCpp * sqwm - RCpm * sqwp)) / (2. * mSqq);
        t.element({0, 6, 0, 0, 0, 0, 0}) = -(Ot * rt * WpmM) / (2. * mSqq);
        t.element({0, 7, 0, 0, 0, 0, 0}) = -(Op * rt * WpmM) / (2. * mSqq);
        t.element({1, 0, 0, 0, 0, 0, 0}) = -WpmP / 2.;
        t.element({1, 1, 0, 0, 0, 0, 0}) = -WpmM / 2.;
        t.element({3, 0, 0, 0, 0, 0, 0}) = -WpmP / 2.;
        t.element({3, 1, 0, 0, 0, 0, 0}) = WpmM / 2.;
        t.element({5, 2, 0, 0, 0, 0, 0}) = (rt * (RCmp * sqwm + RCmm * sqwp)) / (2. * mSqq);
        t.element({5, 3, 0, 0, 0, 0, 0}) = -(Ot * rt * WpmP) / (2. * mSqq);
        t.element({5, 4, 0, 0, 0, 0, 0}) = -(Op * rt * WpmP) / (2. * mSqq);
        t.element({5, 5, 0, 0, 0, 0, 0}) = (-(RCpp * rt * sqwm) + RCpm * rt * sqwp) / (2. * mSqq);
        t.element({5, 6, 0, 0, 0, 0, 0}) = (Ot * rt * WpmM) / (2. * mSqq);
        t.element({5, 7, 0, 0, 0, 0, 0}) = (Op * rt * WpmM) / (2. * mSqq);
        t.element({7, 2, 0, 0, 0, 0, 0}) = (rt * (RCmp * sqwm + RCmm * sqwp)) / (2. * mSqq);
        t.element({7, 3, 0, 0, 0, 0, 0}) = -(Ot * rt * WpmP) / (2. * mSqq);
        t.element({7, 4, 0, 0, 0, 0, 0}) = -(Op * rt * WpmP) / (2. * mSqq);
        t.element({7, 5, 0, 0, 0, 0, 0}) = (rt * (RCpp * sqwm - RCpm * sqwp)) / (2. * mSqq);
        t.element({7, 6, 0, 0, 0, 0, 0}) = -(Ot * rt * WpmM) / (2. * mSqq);
        t.element({7, 7, 0, 0, 0, 0, 0}) = -(Op * rt * WpmM) / (2. * mSqq);
        t.element({9, 8, 0, 0, 0, 0, 0}) = 4. * sqwp * CosTt;
        t.element({9, 9, 0, 0, 0, 0, 0}) = -2. * WpmM * CosTt;
        t.element({9, 10, 0, 0, 0, 0, 0}) = 2. * WpmM * CosTt;
        t.element({9, 11, 0, 0, 0, 0, 0}) = -2. * (1 + w) * WpmM * CosTt;
        t.element({0, 2, 0, 0, 0, 0, 1}) = ((1 + rC) * WpmM * SinTt) / (2. * sqmSqq);
        t.element({0, 3, 0, 0, 0, 0, 1}) = (rC * (1 + w) * WpmM * SinTt) / (2. * sqmSqq);
        t.element({0, 4, 0, 0, 0, 0, 1}) = ((1 + w) * WpmM * SinTt) / (2. * sqmSqq);
        t.element({0, 5, 0, 0, 0, 0, 1}) = ((-1 + rC) * WpmP * SinTt) / (2. * sqmSqq);
        t.element({0, 6, 0, 0, 0, 0, 1}) = (rC * (-1 + w) * WpmP * SinTt) / (2. * sqmSqq);
        t.element({0, 7, 0, 0, 0, 0, 1}) = ((-1 + w) * WpmP * SinTt) / (2. * sqmSqq);
        t.element({5, 2, 0, 0, 0, 0, 1}) = ((1 + rC) * WpmM * SinTt) / (2. * sqmSqq);
        t.element({5, 3, 0, 0, 0, 0, 1}) = (rC * (1 + w) * WpmM * SinTt) / (2. * sqmSqq);
        t.element({5, 4, 0, 0, 0, 0, 1}) = ((1 + w) * WpmM * SinTt) / (2. * sqmSqq);
        t.element({5, 5, 0, 0, 0, 0, 1}) = -((-1 + rC) * WpmP * SinTt) / (2. * sqmSqq);
        t.element({5, 6, 0, 0, 0, 0, 1}) = -(rC * (-1 + w) * WpmP * SinTt) / (2. * sqmSqq);
        t.element({5, 7, 0, 0, 0, 0, 1}) = -((-1 + w) * WpmP * SinTt) / (2. * sqmSqq);
        t.element({7, 2, 0, 0, 0, 0, 1}) = ((1 + rC) * WpmM * SinTt) / (2. * sqmSqq);
        t.element({7, 3, 0, 0, 0, 0, 1}) = (rC * (1 + w) * WpmM * SinTt) / (2. * sqmSqq);
        t.element({7, 4, 0, 0, 0, 0, 1}) = ((1 + w) * WpmM * SinTt) / (2. * sqmSqq);
        t.element({7, 5, 0, 0, 0, 0, 1}) = ((-1 + rC) * WpmP * SinTt) / (2. * sqmSqq);
        t.element({7, 6, 0, 0, 0, 0, 1}) = (rC * (-1 + w) * WpmP * SinTt) / (2. * sqmSqq);
        t.element({7, 7, 0, 0, 0, 0, 1}) = ((-1 + w) * WpmP * SinTt) / (2. * sqmSqq);
        t.element({9, 8, 0, 0, 0, 0, 1}) = -4. * rt * sqrt(wp / mSqq) * SinTt;
        t.element({9, 9, 0, 0, 0, 0, 1}) = (2. * rt * WpmM * SinTt) / sqmSqq;
        t.element({9, 10, 0, 0, 0, 0, 1}) = (-2. * rt * WpmM * SinTt) / sqmSqq;
        t.element({9, 11, 0, 0, 0, 0, 1}) = (2. * rt * (1 + w) * WpmM * SinTt) / sqmSqq;
        t.element({6, 2, 0, 0, 1, 1, 0}) = ((1 + rC) * WpmM * SinTt) / (2. * sqmSqq);
        t.element({6, 3, 0, 0, 1, 1, 0}) = (rC * (1 + w) * WpmM * SinTt) / (2. * sqmSqq);
        t.element({6, 4, 0, 0, 1, 1, 0}) = ((1 + w) * WpmM * SinTt) / (2. * sqmSqq);
        t.element({6, 5, 0, 0, 1, 1, 0}) = -((-1 + rC) * WpmP * SinTt) / (2. * sqmSqq);
        t.element({6, 6, 0, 0, 1, 1, 0}) = -(rC * (-1 + w) * WpmP * SinTt) / (2. * sqmSqq);
        t.element({6, 7, 0, 0, 1, 1, 0}) = -((-1 + w) * WpmP * SinTt) / (2. * sqmSqq);
        t.element({8, 2, 0, 0, 1, 1, 0}) = ((1 + rC) * WpmM * SinTt) / (2. * sqmSqq);
        t.element({8, 3, 0, 0, 1, 1, 0}) = (rC * (1 + w) * WpmM * SinTt) / (2. * sqmSqq);
        t.element({8, 4, 0, 0, 1, 1, 0}) = ((1 + w) * WpmM * SinTt) / (2. * sqmSqq);
        t.element({8, 5, 0, 0, 1, 1, 0}) = ((-1 + rC) * WpmP * SinTt) / (2. * sqmSqq);
        t.element({8, 6, 0, 0, 1, 1, 0}) = (rC * (-1 + w) * WpmP * SinTt) / (2. * sqmSqq);
        t.element({8, 7, 0, 0, 1, 1, 0}) = ((-1 + w) * WpmP * SinTt) / (2. * sqmSqq);
        t.element({10, 8, 0, 0, 1, 1, 0}) = 4. * rt * sqrt(wm / mSqq) * SinTt;
        t.element({10, 9, 0, 0, 1, 1, 0}) = (2. * rt * WpmM * SinTt) / sqmSqq;
        t.element({10, 10, 0, 0, 1, 1, 0}) = (-2. * rt * WpmM * SinTt) / sqmSqq;
        t.element({10, 11, 0, 0, 1, 1, 0}) = (2. * rt * (1 + w) * WpmM * SinTt) / sqmSqq;
        t.element({2, 0, 0, 0, 1, 1, 1}) = WpmP / 2.;
        t.element({2, 1, 0, 0, 1, 1, 1}) = WpmM / 2.;
        t.element({4, 0, 0, 0, 1, 1, 1}) = WpmP / 2.;
        t.element({4, 1, 0, 0, 1, 1, 1}) = -WpmM / 2.;
        t.element({6, 2, 0, 0, 1, 1, 1}) = -(rt * (RCmp * sqwm + RCmm * sqwp)) / (2. * mSqq);
        t.element({6, 3, 0, 0, 1, 1, 1}) = (Ot * rt * WpmP) / (2. * mSqq);
        t.element({6, 4, 0, 0, 1, 1, 1}) = (Op * rt * WpmP) / (2. * mSqq);
        t.element({6, 5, 0, 0, 1, 1, 1}) = (rt * (RCpp * sqwm - RCpm * sqwp)) / (2. * mSqq);
        t.element({6, 6, 0, 0, 1, 1, 1}) = -(Ot * rt * WpmM) / (2. * mSqq);
        t.element({6, 7, 0, 0, 1, 1, 1}) = -(Op * rt * WpmM) / (2. * mSqq);
        t.element({8, 2, 0, 0, 1, 1, 1}) = -(rt * (RCmp * sqwm + RCmm * sqwp)) / (2. * mSqq);
        t.element({8, 3, 0, 0, 1, 1, 1}) = (Ot * rt * WpmP) / (2. * mSqq);
        t.element({8, 4, 0, 0, 1, 1, 1}) = (Op * rt * WpmP) / (2. * mSqq);
        t.element({8, 5, 0, 0, 1, 1, 1}) = (-(RCpp * rt * sqwm) + RCpm * rt * sqwp) / (2. * mSqq);
        t.element({8, 6, 0, 0, 1, 1, 1}) = (Ot * rt * WpmM) / (2. * mSqq);
        t.element({8, 7, 0, 0, 1, 1, 1}) = (Op * rt * WpmM) / (2. * mSqq);
        t.element({10, 8, 0, 0, 1, 1, 1}) = 4. * sqwm * CosTt;
        t.element({10, 9, 0, 0, 1, 1, 1}) = 2. * WpmM * CosTt;
        t.element({10, 10, 0, 0, 1, 1, 1}) = -2. * WpmM * CosTt;
        t.element({10, 11, 0, 0, 1, 1, 1}) = 2. * (1 + w) * WpmM * CosTt;
        t.element({0, 2, 0, 1, 0, 0, 0}) = (rt * WpmM * SinTt) / (2. * sqmSqq);
        t.element({0, 5, 0, 1, 0, 0, 0}) = (rt * WpmP * SinTt) / (2. * sqmSqq);
        t.element({5, 2, 0, 1, 0, 0, 0}) = (rt * WpmM * SinTt) / (2. * sqmSqq);
        t.element({5, 5, 0, 1, 0, 0, 0}) = -(rt * WpmP * SinTt) / (2. * sqmSqq);
        t.element({7, 2, 0, 1, 0, 0, 0}) = (rt * WpmM * SinTt) / (2. * sqmSqq);
        t.element({7, 5, 0, 1, 0, 0, 0}) = (rt * WpmP * SinTt) / (2. * sqmSqq);
        t.element({9, 8, 0, 1, 0, 0, 0}) = (4. * (-rC + wm) * SinTt) / sqrt(mSqq * wm);
        t.element({9, 9, 0, 1, 0, 0, 0}) = (-2. * (-1 + rC * wp) * WpmM * SinTt) / sqmSqq;
        t.element({9, 10, 0, 1, 0, 0, 0}) = (-2. * (rC - wm) * WpmM * SinTt) / sqmSqq;
        t.element({0, 2, 0, 1, 0, 0, 1}) = WpmM * CosTtHalfSq;
        t.element({0, 5, 0, 1, 0, 0, 1}) = WpmP * CosTtHalfSq;
        t.element({5, 2, 0, 1, 0, 0, 1}) = WpmM * CosTtHalfSq;
        t.element({5, 5, 0, 1, 0, 0, 1}) = -(WpmP * CosTtHalfSq);
        t.element({7, 2, 0, 1, 0, 0, 1}) = WpmM * CosTtHalfSq;
        t.element({7, 5, 0, 1, 0, 0, 1}) = WpmP * CosTtHalfSq;
        t.element({9, 8, 0, 1, 0, 0, 1}) = (-8. * rt * sqwm * (-1 + rC * wp) * CosTtHalfSq) / mSqq;
        t.element({9, 9, 0, 1, 0, 0, 1}) = (-4. * rt * (-1 + rC * wp) * WpmM * CosTtHalfSq) / mSqq;
        t.element({9, 10, 0, 1, 0, 0, 1}) = (-4. * rt * (rC - wm) * WpmM * CosTtHalfSq) / mSqq;
        t.element({6, 2, 0, 1, 1, 1, 0}) = -(WpmM * SinTtHalfSq);
        t.element({6, 5, 0, 1, 1, 1, 0}) = WpmP * SinTtHalfSq;
        t.element({8, 2, 0, 1, 1, 1, 0}) = -(WpmM * SinTtHalfSq);
        t.element({8, 5, 0, 1, 1, 1, 0}) = -(WpmP * SinTtHalfSq);
        t.element({10, 8, 0, 1, 1, 1, 0}) = (-8. * rt * (-1 + rC * wm) * sqwp * SinTtHalfSq) / mSqq;
        t.element({10, 9, 0, 1, 1, 1, 0}) = (4. * rt * (-1 + rC * wm) * WpmM * SinTtHalfSq) / mSqq;
        t.element({10, 10, 0, 1, 1, 1, 0}) = (-4. * rt * (-rC + wp) * WpmM * SinTtHalfSq) / mSqq;
        t.element({6, 2, 0, 1, 1, 1, 1}) = -(rt * WpmM * SinTt) / (2. * sqmSqq);
        t.element({6, 5, 0, 1, 1, 1, 1}) = (rt * WpmP * SinTt) / (2. * sqmSqq);
        t.element({8, 2, 0, 1, 1, 1, 1}) = -(rt * WpmM * SinTt) / (2. * sqmSqq);
        t.element({8, 5, 0, 1, 1, 1, 1}) = -(rt * WpmP * SinTt) / (2. * sqmSqq);
        t.element({10, 8, 0, 1, 1, 1, 1}) = (4. * (-rC + wp) * SinTt) / sqrt(mSqq * wp);
        t.element({10, 9, 0, 1, 1, 1, 1}) = (2. * (-1 + rC * wm) * WpmM * SinTt) / sqmSqq;
        t.element({10, 10, 0, 1, 1, 1, 1}) = (-2. * (-rC + wp) * WpmM * SinTt) / sqmSqq;
        t.element({0, 2, 1, 0, 0, 0, 0}) = -(rt * WpmM * SinTt) / (2. * sqmSqq);
        t.element({0, 5, 1, 0, 0, 0, 0}) = (rt * WpmP * SinTt) / (2. * sqmSqq);
        t.element({5, 2, 1, 0, 0, 0, 0}) = -(rt * WpmM * SinTt) / (2. * sqmSqq);
        t.element({5, 5, 1, 0, 0, 0, 0}) = -(rt * WpmP * SinTt) / (2. * sqmSqq);
        t.element({7, 2, 1, 0, 0, 0, 0}) = -(rt * WpmM * SinTt) / (2. * sqmSqq);
        t.element({7, 5, 1, 0, 0, 0, 0}) = (rt * WpmP * SinTt) / (2. * sqmSqq);
        t.element({9, 8, 1, 0, 0, 0, 0}) = (4. * (-rC + wp) * SinTt) / sqrt(mSqq * wp);
        t.element({9, 9, 1, 0, 0, 0, 0}) = (2. * (-1 + rC * wm) * WpmM * SinTt) / sqmSqq;
        t.element({9, 10, 1, 0, 0, 0, 0}) = (-2. * (-rC + wp) * WpmM * SinTt) / sqmSqq;
        t.element({0, 2, 1, 0, 0, 0, 1}) = WpmM * SinTtHalfSq;
        t.element({0, 5, 1, 0, 0, 0, 1}) = -(WpmP * SinTtHalfSq);
        t.element({5, 2, 1, 0, 0, 0, 1}) = WpmM * SinTtHalfSq;
        t.element({5, 5, 1, 0, 0, 0, 1}) = WpmP * SinTtHalfSq;
        t.element({7, 2, 1, 0, 0, 0, 1}) = WpmM * SinTtHalfSq;
        t.element({7, 5, 1, 0, 0, 0, 1}) = -(WpmP * SinTtHalfSq);
        t.element({9, 8, 1, 0, 0, 0, 1}) = (8. * rt * (-1 + rC * wm) * sqwp * SinTtHalfSq) / mSqq;
        t.element({9, 9, 1, 0, 0, 0, 1}) = (-4. * rt * (-1 + rC * wm) * WpmM * SinTtHalfSq) / mSqq;
        t.element({9, 10, 1, 0, 0, 0, 1}) = (-4. * rt * (rC - wp) * WpmM * SinTtHalfSq) / mSqq;
        t.element({6, 2, 1, 0, 1, 1, 0}) = -(WpmM * CosTtHalfSq);
        t.element({6, 5, 1, 0, 1, 1, 0}) = -(WpmP * CosTtHalfSq);
        t.element({8, 2, 1, 0, 1, 1, 0}) = -(WpmM * CosTtHalfSq);
        t.element({8, 5, 1, 0, 1, 1, 0}) = WpmP * CosTtHalfSq;
        t.element({10, 8, 1, 0, 1, 1, 0}) = (8. * rt * sqwm * (-1 + rC * wp) * CosTtHalfSq) / mSqq;
        t.element({10, 9, 1, 0, 1, 1, 0}) = (4. * rt * (-1 + rC * wp) * WpmM * CosTtHalfSq) / mSqq;
        t.element({10, 10, 1, 0, 1, 1, 0}) = (4. * rt * (rC - wm) * WpmM * CosTtHalfSq) / mSqq;
        t.element({6, 2, 1, 0, 1, 1, 1}) = (rt * WpmM * SinTt) / (2. * sqmSqq);
        t.element({6, 5, 1, 0, 1, 1, 1}) = (rt * WpmP * SinTt) / (2. * sqmSqq);
        t.element({8, 2, 1, 0, 1, 1, 1}) = (rt * WpmM * SinTt) / (2. * sqmSqq);
        t.element({8, 5, 1, 0, 1, 1, 1}) = -(rt * WpmP * SinTt) / (2. * sqmSqq);
        t.element({10, 8, 1, 0, 1, 1, 1}) = (4. * (-rC + wm) * SinTt) / sqrt(mSqq * wm);
        t.element({10, 9, 1, 0, 1, 1, 1}) = (-2. * (-1 + rC * wp) * WpmM * SinTt) / sqmSqq;
        t.element({10, 10, 1, 0, 1, 1, 1}) = (-2. * (rC - wm) * WpmM * SinTt) / sqmSqq;
        t.element({0, 2, 1, 1, 0, 0, 0}) = (rt * (RCmp * sqwm + RCmm * sqwp)) / (2. * mSqq);
        t.element({0, 3, 1, 1, 0, 0, 0}) = -(Ot * rt * WpmP) / (2. * mSqq);
        t.element({0, 4, 1, 1, 0, 0, 0}) = -(Op * rt * WpmP) / (2. * mSqq);
        t.element({0, 5, 1, 1, 0, 0, 0}) = (-(RCpp * rt * sqwm) + RCpm * rt * sqwp) / (2. * mSqq);
        t.element({0, 6, 1, 1, 0, 0, 0}) = (Ot * rt * WpmM) / (2. * mSqq);
        t.element({0, 7, 1, 1, 0, 0, 0}) = (Op * rt * WpmM) / (2. * mSqq);
        t.element({1, 0, 1, 1, 0, 0, 0}) = -WpmP / 2.;
        t.element({1, 1, 1, 1, 0, 0, 0}) = WpmM / 2.;
        t.element({3, 0, 1, 1, 0, 0, 0}) = -WpmP / 2.;
        t.element({3, 1, 1, 1, 0, 0, 0}) = -WpmM / 2.;
        t.element({5, 2, 1, 1, 0, 0, 0}) = (rt * (RCmp * sqwm + RCmm * sqwp)) / (2. * mSqq);
        t.element({5, 3, 1, 1, 0, 0, 0}) = -(Ot * rt * WpmP) / (2. * mSqq);
        t.element({5, 4, 1, 1, 0, 0, 0}) = -(Op * rt * WpmP) / (2. * mSqq);
        t.element({5, 5, 1, 1, 0, 0, 0}) = (rt * (RCpp * sqwm - RCpm * sqwp)) / (2. * mSqq);
        t.element({5, 6, 1, 1, 0, 0, 0}) = -(Ot * rt * WpmM) / (2. * mSqq);
        t.element({5, 7, 1, 1, 0, 0, 0}) = -(Op * rt * WpmM) / (2. * mSqq);
        t.element({7, 2, 1, 1, 0, 0, 0}) = (rt * (RCmp * sqwm + RCmm * sqwp)) / (2. * mSqq);
        t.element({7, 3, 1, 1, 0, 0, 0}) = -(Ot * rt * WpmP) / (2. * mSqq);
        t.element({7, 4, 1, 1, 0, 0, 0}) = -(Op * rt * WpmP) / (2. * mSqq);
        t.element({7, 5, 1, 1, 0, 0, 0}) = (-(RCpp * rt * sqwm) + RCpm * rt * sqwp) / (2. * mSqq);
        t.element({7, 6, 1, 1, 0, 0, 0}) = (Ot * rt * WpmM) / (2. * mSqq);
        t.element({7, 7, 1, 1, 0, 0, 0}) = (Op * rt * WpmM) / (2. * mSqq);
        t.element({9, 8, 1, 1, 0, 0, 0}) = -4. * sqwm * CosTt;
        t.element({9, 9, 1, 1, 0, 0, 0}) = -2. * WpmM * CosTt;
        t.element({9, 10, 1, 1, 0, 0, 0}) = 2. * WpmM * CosTt;
        t.element({9, 11, 1, 1, 0, 0, 0}) = -2. * (1 + w) * WpmM * CosTt;
        t.element({0, 2, 1, 1, 0, 0, 1}) = ((1 + rC) * WpmM * SinTt) / (2. * sqmSqq);
        t.element({0, 3, 1, 1, 0, 0, 1}) = (rC * (1 + w) * WpmM * SinTt) / (2. * sqmSqq);
        t.element({0, 4, 1, 1, 0, 0, 1}) = ((1 + w) * WpmM * SinTt) / (2. * sqmSqq);
        t.element({0, 5, 1, 1, 0, 0, 1}) = -((-1 + rC) * WpmP * SinTt) / (2. * sqmSqq);
        t.element({0, 6, 1, 1, 0, 0, 1}) = -(rC * (-1 + w) * WpmP * SinTt) / (2. * sqmSqq);
        t.element({0, 7, 1, 1, 0, 0, 1}) = -((-1 + w) * WpmP * SinTt) / (2. * sqmSqq);
        t.element({5, 2, 1, 1, 0, 0, 1}) = ((1 + rC) * WpmM * SinTt) / (2. * sqmSqq);
        t.element({5, 3, 1, 1, 0, 0, 1}) = (rC * (1 + w) * WpmM * SinTt) / (2. * sqmSqq);
        t.element({5, 4, 1, 1, 0, 0, 1}) = ((1 + w) * WpmM * SinTt) / (2. * sqmSqq);
        t.element({5, 5, 1, 1, 0, 0, 1}) = ((-1 + rC) * WpmP * SinTt) / (2. * sqmSqq);
        t.element({5, 6, 1, 1, 0, 0, 1}) = (rC * (-1 + w) * WpmP * SinTt) / (2. * sqmSqq);
        t.element({5, 7, 1, 1, 0, 0, 1}) = ((-1 + w) * WpmP * SinTt) / (2. * sqmSqq);
        t.element({7, 2, 1, 1, 0, 0, 1}) = ((1 + rC) * WpmM * SinTt) / (2. * sqmSqq);
        t.element({7, 3, 1, 1, 0, 0, 1}) = (rC * (1 + w) * WpmM * SinTt) / (2. * sqmSqq);
        t.element({7, 4, 1, 1, 0, 0, 1}) = ((1 + w) * WpmM * SinTt) / (2. * sqmSqq);
        t.element({7, 5, 1, 1, 0, 0, 1}) = -((-1 + rC) * WpmP * SinTt) / (2. * sqmSqq);
        t.element({7, 6, 1, 1, 0, 0, 1}) = -(rC * (-1 + w) * WpmP * SinTt) / (2. * sqmSqq);
        t.element({7, 7, 1, 1, 0, 0, 1}) = -((-1 + w) * WpmP * SinTt) / (2. * sqmSqq);
        t.element({9, 8, 1, 1, 0, 0, 1}) = 4. * rt * sqrt(wm / mSqq) * SinTt;
        t.element({9, 9, 1, 1, 0, 0, 1}) = (2. * rt * WpmM * SinTt) / sqmSqq;
        t.element({9, 10, 1, 1, 0, 0, 1}) = (-2. * rt * WpmM * SinTt) / sqmSqq;
        t.element({9, 11, 1, 1, 0, 0, 1}) = (2. * rt * (1 + w) * WpmM * SinTt) / sqmSqq;
        t.element({6, 2, 1, 1, 1, 1, 0}) = ((1 + rC) * WpmM * SinTt) / (2. * sqmSqq);
        t.element({6, 3, 1, 1, 1, 1, 0}) = (rC * (1 + w) * WpmM * SinTt) / (2. * sqmSqq);
        t.element({6, 4, 1, 1, 1, 1, 0}) = ((1 + w) * WpmM * SinTt) / (2. * sqmSqq);
        t.element({6, 5, 1, 1, 1, 1, 0}) = ((-1 + rC) * WpmP * SinTt) / (2. * sqmSqq);
        t.element({6, 6, 1, 1, 1, 1, 0}) = (rC * (-1 + w) * WpmP * SinTt) / (2. * sqmSqq);
        t.element({6, 7, 1, 1, 1, 1, 0}) = ((-1 + w) * WpmP * SinTt) / (2. * sqmSqq);
        t.element({8, 2, 1, 1, 1, 1, 0}) = ((1 + rC) * WpmM * SinTt) / (2. * sqmSqq);
        t.element({8, 3, 1, 1, 1, 1, 0}) = (rC * (1 + w) * WpmM * SinTt) / (2. * sqmSqq);
        t.element({8, 4, 1, 1, 1, 1, 0}) = ((1 + w) * WpmM * SinTt) / (2. * sqmSqq);
        t.element({8, 5, 1, 1, 1, 1, 0}) = -((-1 + rC) * WpmP * SinTt) / (2. * sqmSqq);
        t.element({8, 6, 1, 1, 1, 1, 0}) = -(rC * (-1 + w) * WpmP * SinTt) / (2. * sqmSqq);
        t.element({8, 7, 1, 1, 1, 1, 0}) = -((-1 + w) * WpmP * SinTt) / (2. * sqmSqq);
        t.element({10, 8, 1, 1, 1, 1, 0}) = -4. * rt * sqrt(wp / mSqq) * SinTt;
        t.element({10, 9, 1, 1, 1, 1, 0}) = (2. * rt * WpmM * SinTt) / sqmSqq;
        t.element({10, 10, 1, 1, 1, 1, 0}) = (-2. * rt * WpmM * SinTt) / sqmSqq;
        t.element({10, 11, 1, 1, 1, 1, 0}) = (2. * rt * (1 + w) * WpmM * SinTt) / sqmSqq;
        t.element({2, 0, 1, 1, 1, 1, 1}) = WpmP / 2.;
        t.element({2, 1, 1, 1, 1, 1, 1}) = -WpmM / 2.;
        t.element({4, 0, 1, 1, 1, 1, 1}) = WpmP / 2.;
        t.element({4, 1, 1, 1, 1, 1, 1}) = WpmM / 2.;
        t.element({6, 2, 1, 1, 1, 1, 1}) = -(rt * (RCmp * sqwm + RCmm * sqwp)) / (2. * mSqq);
        t.element({6, 3, 1, 1, 1, 1, 1}) = (Ot * rt * WpmP) / (2. * mSqq);
        t.element({6, 4, 1, 1, 1, 1, 1}) = (Op * rt * WpmP) / (2. * mSqq);
        t.element({6, 5, 1, 1, 1, 1, 1}) = (-(RCpp * rt * sqwm) + RCpm * rt * sqwp) / (2. * mSqq);
        t.element({6, 6, 1, 1, 1, 1, 1}) = (Ot * rt * WpmM) / (2. * mSqq);
        t.element({6, 7, 1, 1, 1, 1, 1}) = (Op * rt * WpmM) / (2. * mSqq);
        t.element({8, 2, 1, 1, 1, 1, 1}) = -(rt * (RCmp * sqwm + RCmm * sqwp)) / (2. * mSqq);
        t.element({8, 3, 1, 1, 1, 1, 1}) = (Ot * rt * WpmP) / (2. * mSqq);
        t.element({8, 4, 1, 1, 1, 1, 1}) = (Op * rt * WpmP) / (2. * mSqq);
        t.element({8, 5, 1, 1, 1, 1, 1}) = (rt * (RCpp * sqwm - RCpm * sqwp)) / (2. * mSqq);
        t.element({8, 6, 1, 1, 1, 1, 1}) = -(Ot * rt * WpmM) / (2. * mSqq);
        t.element({8, 7, 1, 1, 1, 1, 1}) = -(Op * rt * WpmM) / (2. * mSqq);
        t.element({10, 8, 1, 1, 1, 1, 1}) = -4. * sqwp * CosTt;
        t.element({10, 9, 1, 1, 1, 1, 1}) = 2. * WpmM * CosTt;
        t.element({10, 10, 1, 1, 1, 1, 1}) = -2. * WpmM * CosTt;
        t.element({10, 11, 1, 1, 1, 1, 1}) = 2. * (1 + w) * WpmM * CosTt;

        t *= TwoSqTwoGfSq;
    }

    void AmplLbLcLepNu::addRefs() const {
        if (!getSettingsHandler()->checkReference("Bernlochner:2018bfn")) {
            string ref =
                "@article{Bernlochner:2018bfn,\n"
                "   author         = \"Bernlochner, Florian U. and Ligeti, Zoltan and Robinson, Dean J. and Sutcliffe, "
                "William L.\",\n"
                "   title          = \"{Precise predictions for $\\Lambda_b \\to \\Lambda_c$ semileptonic decays}\",\n"
                "   journal        = \"Phys. Rev.\",\n"
                "   volume         = \"D99\",\n"
                "   year           = \"2019\",\n"
                "   number         = \"5\",\n"
                "   pages          = \"055008\",\n"
                "   doi            = \"10.1103/PhysRevD.99.055008\",\n"
                "   eprint         = \"1812.07593\",\n"
                "   archivePrefix  = \"arXiv\",\n"
                "   primaryClass   = \"hep-ph\",\n"
                "   SLACcitation   = \"%%CITATION = ARXIV:1812.07593;%%\"\n"
                "}\n";
            getSettingsHandler()->addReference("Bernlochner:2018bfn", ref);
        }
    }

} // namespace Hammer
