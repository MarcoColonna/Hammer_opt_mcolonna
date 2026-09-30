///
/// @file  AmplLbLcstar12LepNu.cc
/// @brief \f$ \Lambda_b \rightarrow L_c^*(2595) \tau\nu \f$ amplitude
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/Amplitudes/AmplLbLcstar12LepNu.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;
using namespace complex_literals;

namespace Hammer {

    namespace MD = MultiDimensional;

    AmplLbLcstar12LepNu::AmplLbLcstar12LepNu() {
        // Create tensor rank and dimensions
        IndexList dims{{11, 12, 2, 2, 2, 2, 2}};
        string name{"AmplLbLcstar12LepNu"};
        addProcessSignature(-PID::LAMBDAB, {PID::LAMBDACSTAR12MINUS, PID::NU_TAU, PID::ANTITAU});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCTAUNU, FF_LBLCSTAR12, SPIN_LB, SPIN_LCSTAR12,
                                                          SPIN_NUTAU, SPIN_NUTAU_REF, SPIN_TAUP})});

        addProcessSignature(-PID::LAMBDAB, {PID::LAMBDACSTAR12MINUS, PID::NU_MU, PID::ANTIMUON});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCMUNU, FF_LBLCSTAR12, SPIN_LB, SPIN_LCSTAR12,
                                                          SPIN_NUMU, SPIN_NUMU_REF, SPIN_MUP})});

        addProcessSignature(-PID::LAMBDAB, {PID::LAMBDACSTAR12MINUS, PID::NU_E, PID::POSITRON});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCENU, FF_LBLCSTAR12, SPIN_LB, SPIN_LCSTAR12, SPIN_NUE,
                                                          SPIN_NUE_REF, SPIN_EP})});

        setSignatureIndex();
        _multiplicity = 2ul;
    }

    void AmplLbLcstar12LepNu::eval(const Particle& parent, const ParticleList& daughters,
                                   const ParticleList& /*references*/) {
        // Momenta
        const FourMomentum& pLbmes = parent.momentum();
        const FourMomentum& pLcstar12mes = daughters[0].momentum();
        const FourMomentum& kNuTau = daughters[1].momentum();
        const FourMomentum& pTau = daughters[2].momentum();

        // kinematic objects
        const double Mb = pLbmes.mass();
        const double Mb2 = Mb * Mb;
        const double Mc = pLcstar12mes.mass();
        const double Mc2 = Mc * Mc;
        const double Mt = pTau.mass();
        const double Sqq = Mb2 + Mc2 - 2. * (pLbmes * pLcstar12mes);
        const double Ew = (Mb2 - Mc2 + Sqq) / (2 * Mb);
        const double BNuTau = (pLbmes * kNuTau);
        const double NuTauQ = (pLbmes * kNuTau) - (pLcstar12mes * kNuTau);
        const double BQ = Mb2 - (pLbmes * pLcstar12mes);

        const double mSqq = Sqq / Mb2;
        const double SqmSqq = sqrt(mSqq);
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
        // const double w2m1 = w*w-1;
        const double Sqw2m1 = sqrt(w * w - 1);
        const double wp = w + Sqw2m1;
        const double wm = w - Sqw2m1;
        const double sqwp = sqrt(wp);
        const double sqwm = sqrt(wm);
        const double WpmP = sqwp + sqwm;
        const double WpmM = sqwp - sqwm;

        const double Rp = (rC + 1) * WpmM + (rC - 1) * WpmP * CosTt;
        const double Rm = (rC - 1) * WpmP + (rC + 1) * WpmM * CosTt;

        const double Op = rC - w + Sqw2m1 * CosTt;
        const double Ot = rC * w - 1 + rC * Sqw2m1 * CosTt;
        // const double Oz = Sqw2m1 + (rC-w)*CosTt;

        const double TwoSqTwoGfSq = TwoSqTwoGFermi * Mb2 * sqrt(mSqq - pow(rt, 2.)) * sqrt(rC);

        // initialize tensor elements to zero
        Tensor& t = getTensor();
        t.clearData();

        // set non-zero tensor elements
        t.element({0, 2, 0, 0, 0, 0, 0}) = (Rp * rt) / (2. * mSqq);
        t.element({0, 3, 0, 0, 0, 0, 0}) = (Ot * rt * WpmM) / (2. * mSqq);
        t.element({0, 4, 0, 0, 0, 0, 0}) = (Op * rt * WpmM) / (2. * mSqq);
        t.element({0, 5, 0, 0, 0, 0, 0}) = (Rm * rt) / (2. * mSqq);
        t.element({0, 6, 0, 0, 0, 0, 0}) = (Ot * rt * WpmP) / (2. * mSqq);
        t.element({0, 7, 0, 0, 0, 0, 0}) = (Op * rt * WpmP) / (2. * mSqq);
        t.element({1, 0, 0, 0, 0, 0, 0}) = -WpmM / 2.;
        t.element({1, 1, 0, 0, 0, 0, 0}) = -WpmP / 2.;
        t.element({3, 0, 0, 0, 0, 0, 0}) = -WpmM / 2.;
        t.element({3, 1, 0, 0, 0, 0, 0}) = WpmP / 2.;
        t.element({5, 2, 0, 0, 0, 0, 0}) = (Rp * rt) / (2. * mSqq);
        t.element({5, 3, 0, 0, 0, 0, 0}) = (Ot * rt * WpmM) / (2. * mSqq);
        t.element({5, 4, 0, 0, 0, 0, 0}) = (Op * rt * WpmM) / (2. * mSqq);
        t.element({5, 5, 0, 0, 0, 0, 0}) = -(Rm * rt) / (2. * mSqq);
        t.element({5, 6, 0, 0, 0, 0, 0}) = -(Ot * rt * WpmP) / (2. * mSqq);
        t.element({5, 7, 0, 0, 0, 0, 0}) = -(Op * rt * WpmP) / (2. * mSqq);
        t.element({7, 2, 0, 0, 0, 0, 0}) = (Rp * rt) / (2. * mSqq);
        t.element({7, 3, 0, 0, 0, 0, 0}) = (Ot * rt * WpmM) / (2. * mSqq);
        t.element({7, 4, 0, 0, 0, 0, 0}) = (Op * rt * WpmM) / (2. * mSqq);
        t.element({7, 5, 0, 0, 0, 0, 0}) = (Rm * rt) / (2. * mSqq);
        t.element({7, 6, 0, 0, 0, 0, 0}) = (Ot * rt * WpmP) / (2. * mSqq);
        t.element({7, 7, 0, 0, 0, 0, 0}) = (Op * rt * WpmP) / (2. * mSqq);
        t.element({9, 8, 0, 0, 0, 0, 0}) = 4. * CosTt * sqwp;
        t.element({9, 9, 0, 0, 0, 0, 0}) = 2. * CosTt * WpmP;
        t.element({9, 10, 0, 0, 0, 0, 0}) = 2. * CosTt * WpmP;
        t.element({9, 11, 0, 0, 0, 0, 0}) = -2. * CosTt * Sqw2m1 * WpmM;
        t.element({0, 2, 0, 0, 0, 0, 1}) = -((-1 + rC) * SinTt * WpmP) / (2. * SqmSqq);
        t.element({0, 3, 0, 0, 0, 0, 1}) = -(rC * SinTt * Sqw2m1 * WpmM) / (2. * SqmSqq);
        t.element({0, 4, 0, 0, 0, 0, 1}) = -(SinTt * Sqw2m1 * WpmM) / (2. * SqmSqq);
        t.element({0, 5, 0, 0, 0, 0, 1}) = -((1 + rC) * SinTt * WpmM) / (2. * SqmSqq);
        t.element({0, 6, 0, 0, 0, 0, 1}) = -(rC * SinTt * Sqw2m1 * WpmP) / (2. * SqmSqq);
        t.element({0, 7, 0, 0, 0, 0, 1}) = -(SinTt * Sqw2m1 * WpmP) / (2. * SqmSqq);
        t.element({5, 2, 0, 0, 0, 0, 1}) = -((-1 + rC) * SinTt * WpmP) / (2. * SqmSqq);
        t.element({5, 3, 0, 0, 0, 0, 1}) = -(rC * SinTt * Sqw2m1 * WpmM) / (2. * SqmSqq);
        t.element({5, 4, 0, 0, 0, 0, 1}) = -(SinTt * Sqw2m1 * WpmM) / (2. * SqmSqq);
        t.element({5, 5, 0, 0, 0, 0, 1}) = ((1 + rC) * SinTt * WpmM) / (2. * SqmSqq);
        t.element({5, 6, 0, 0, 0, 0, 1}) = (rC * SinTt * Sqw2m1 * WpmP) / (2. * SqmSqq);
        t.element({5, 7, 0, 0, 0, 0, 1}) = (SinTt * Sqw2m1 * WpmP) / (2. * SqmSqq);
        t.element({7, 2, 0, 0, 0, 0, 1}) = -((-1 + rC) * SinTt * WpmP) / (2. * SqmSqq);
        t.element({7, 3, 0, 0, 0, 0, 1}) = -(rC * SinTt * Sqw2m1 * WpmM) / (2. * SqmSqq);
        t.element({7, 4, 0, 0, 0, 0, 1}) = -(SinTt * Sqw2m1 * WpmM) / (2. * SqmSqq);
        t.element({7, 5, 0, 0, 0, 0, 1}) = -((1 + rC) * SinTt * WpmM) / (2. * SqmSqq);
        t.element({7, 6, 0, 0, 0, 0, 1}) = -(rC * SinTt * Sqw2m1 * WpmP) / (2. * SqmSqq);
        t.element({7, 7, 0, 0, 0, 0, 1}) = -(SinTt * Sqw2m1 * WpmP) / (2. * SqmSqq);
        t.element({9, 8, 0, 0, 0, 0, 1}) = (-4. * rt * SinTt * sqwp) / SqmSqq;
        t.element({9, 9, 0, 0, 0, 0, 1}) = (-2. * rt * SinTt * WpmP) / SqmSqq;
        t.element({9, 10, 0, 0, 0, 0, 1}) = (-2. * rt * SinTt * WpmP) / SqmSqq;
        t.element({9, 11, 0, 0, 0, 0, 1}) = (2. * rt * SinTt * Sqw2m1 * WpmM) / SqmSqq;
        t.element({6, 2, 0, 0, 1, 1, 0}) = -((-1 + rC) * SinTt * WpmP) / (2. * SqmSqq);
        t.element({6, 3, 0, 0, 1, 1, 0}) = -(rC * SinTt * Sqw2m1 * WpmM) / (2. * SqmSqq);
        t.element({6, 4, 0, 0, 1, 1, 0}) = -(SinTt * Sqw2m1 * WpmM) / (2. * SqmSqq);
        t.element({6, 5, 0, 0, 1, 1, 0}) = ((1 + rC) * SinTt * WpmM) / (2. * SqmSqq);
        t.element({6, 6, 0, 0, 1, 1, 0}) = (rC * SinTt * Sqw2m1 * WpmP) / (2. * SqmSqq);
        t.element({6, 7, 0, 0, 1, 1, 0}) = (SinTt * Sqw2m1 * WpmP) / (2. * SqmSqq);
        t.element({8, 2, 0, 0, 1, 1, 0}) = -((-1 + rC) * SinTt * WpmP) / (2. * SqmSqq);
        t.element({8, 3, 0, 0, 1, 1, 0}) = -(rC * SinTt * Sqw2m1 * WpmM) / (2. * SqmSqq);
        t.element({8, 4, 0, 0, 1, 1, 0}) = -(SinTt * Sqw2m1 * WpmM) / (2. * SqmSqq);
        t.element({8, 5, 0, 0, 1, 1, 0}) = -((1 + rC) * SinTt * WpmM) / (2. * SqmSqq);
        t.element({8, 6, 0, 0, 1, 1, 0}) = -(rC * SinTt * Sqw2m1 * WpmP) / (2. * SqmSqq);
        t.element({8, 7, 0, 0, 1, 1, 0}) = -(SinTt * Sqw2m1 * WpmP) / (2. * SqmSqq);
        t.element({10, 8, 0, 0, 1, 1, 0}) = (-4. * rt * SinTt * sqwm) / SqmSqq;
        t.element({10, 9, 0, 0, 1, 1, 0}) = (-2. * rt * SinTt * WpmP) / SqmSqq;
        t.element({10, 10, 0, 0, 1, 1, 0}) = (-2. * rt * SinTt * WpmP) / SqmSqq;
        t.element({10, 11, 0, 0, 1, 1, 0}) = (2. * rt * SinTt * Sqw2m1 * WpmM) / SqmSqq;
        t.element({2, 0, 0, 0, 1, 1, 1}) = WpmM / 2.;
        t.element({2, 1, 0, 0, 1, 1, 1}) = WpmP / 2.;
        t.element({4, 0, 0, 0, 1, 1, 1}) = WpmM / 2.;
        t.element({4, 1, 0, 0, 1, 1, 1}) = -WpmP / 2.;
        t.element({6, 2, 0, 0, 1, 1, 1}) = -(Rp * rt) / (2. * mSqq);
        t.element({6, 3, 0, 0, 1, 1, 1}) = -(Ot * rt * WpmM) / (2. * mSqq);
        t.element({6, 4, 0, 0, 1, 1, 1}) = -(Op * rt * WpmM) / (2. * mSqq);
        t.element({6, 5, 0, 0, 1, 1, 1}) = (Rm * rt) / (2. * mSqq);
        t.element({6, 6, 0, 0, 1, 1, 1}) = (Ot * rt * WpmP) / (2. * mSqq);
        t.element({6, 7, 0, 0, 1, 1, 1}) = (Op * rt * WpmP) / (2. * mSqq);
        t.element({8, 2, 0, 0, 1, 1, 1}) = -(Rp * rt) / (2. * mSqq);
        t.element({8, 3, 0, 0, 1, 1, 1}) = -(Ot * rt * WpmM) / (2. * mSqq);
        t.element({8, 4, 0, 0, 1, 1, 1}) = -(Op * rt * WpmM) / (2. * mSqq);
        t.element({8, 5, 0, 0, 1, 1, 1}) = -(Rm * rt) / (2. * mSqq);
        t.element({8, 6, 0, 0, 1, 1, 1}) = -(Ot * rt * WpmP) / (2. * mSqq);
        t.element({8, 7, 0, 0, 1, 1, 1}) = -(Op * rt * WpmP) / (2. * mSqq);
        t.element({10, 8, 0, 0, 1, 1, 1}) = -4. * CosTt * sqwm;
        t.element({10, 9, 0, 0, 1, 1, 1}) = -2. * CosTt * WpmP;
        t.element({10, 10, 0, 0, 1, 1, 1}) = -2. * CosTt * WpmP;
        t.element({10, 11, 0, 0, 1, 1, 1}) = 2. * CosTt * Sqw2m1 * WpmM;
        t.element({0, 2, 0, 1, 0, 0, 0}) = -(rt * SinTt * WpmP) / (2. * SqmSqq);
        t.element({0, 5, 0, 1, 0, 0, 0}) = -(rt * SinTt * WpmM) / (2. * SqmSqq);
        t.element({5, 2, 0, 1, 0, 0, 0}) = -(rt * SinTt * WpmP) / (2. * SqmSqq);
        t.element({5, 5, 0, 1, 0, 0, 0}) = (rt * SinTt * WpmM) / (2. * SqmSqq);
        t.element({7, 2, 0, 1, 0, 0, 0}) = -(rt * SinTt * WpmP) / (2. * SqmSqq);
        t.element({7, 5, 0, 1, 0, 0, 0}) = -(rt * SinTt * WpmM) / (2. * SqmSqq);
        t.element({9, 8, 0, 1, 0, 0, 0}) = (-4. * SinTt * sqwm * (rC * wp - 1)) / SqmSqq;
        t.element({9, 9, 0, 1, 0, 0, 0}) = (-2. * SinTt * (rC * wp - 1) * WpmP) / SqmSqq;
        t.element({9, 10, 0, 1, 0, 0, 0}) = (-2. * SinTt * (rC - wm) * WpmP) / SqmSqq;
        t.element({0, 2, 0, 1, 0, 0, 1}) = -(CosTtHalfSq * WpmP);
        t.element({0, 5, 0, 1, 0, 0, 1}) = -(CosTtHalfSq * WpmM);
        t.element({5, 2, 0, 1, 0, 0, 1}) = -(CosTtHalfSq * WpmP);
        t.element({5, 5, 0, 1, 0, 0, 1}) = CosTtHalfSq * WpmM;
        t.element({7, 2, 0, 1, 0, 0, 1}) = -(CosTtHalfSq * WpmP);
        t.element({7, 5, 0, 1, 0, 0, 1}) = -(CosTtHalfSq * WpmM);
        t.element({9, 8, 0, 1, 0, 0, 1}) = (-8 * CosTtHalfSq * rt * sqwm * (rC * wp - 1)) / mSqq;
        t.element({9, 9, 0, 1, 0, 0, 1}) = (-4. * CosTtHalfSq * rt * (rC * wp - 1) * WpmP) / mSqq;
        t.element({9, 10, 0, 1, 0, 0, 1}) = (-4. * CosTtHalfSq * rt * (rC - wm) * WpmP) / mSqq;
        t.element({6, 2, 0, 1, 1, 1, 0}) = SinTtHalfSq * WpmP;
        t.element({6, 5, 0, 1, 1, 1, 0}) = -(SinTtHalfSq * WpmM);
        t.element({8, 2, 0, 1, 1, 1, 0}) = SinTtHalfSq * WpmP;
        t.element({8, 5, 0, 1, 1, 1, 0}) = SinTtHalfSq * WpmM;
        t.element({10, 8, 0, 1, 1, 1, 0}) = (8 * rt * SinTtHalfSq * (rC * wm - 1) * sqwp) / mSqq;
        t.element({10, 9, 0, 1, 1, 1, 0}) = (4. * rt * SinTtHalfSq * (rC * wm - 1) * WpmP) / mSqq;
        t.element({10, 10, 0, 1, 1, 1, 0}) = (4. * rt * SinTtHalfSq * (rC - wp) * WpmP) / mSqq;
        t.element({6, 2, 0, 1, 1, 1, 1}) = (rt * SinTt * WpmP) / (2. * SqmSqq);
        t.element({6, 5, 0, 1, 1, 1, 1}) = -(rt * SinTt * WpmM) / (2. * SqmSqq);
        t.element({8, 2, 0, 1, 1, 1, 1}) = (rt * SinTt * WpmP) / (2. * SqmSqq);
        t.element({8, 5, 0, 1, 1, 1, 1}) = (rt * SinTt * WpmM) / (2. * SqmSqq);
        t.element({10, 8, 0, 1, 1, 1, 1}) = (4. * SinTt * (rC * wm - 1) * sqwp) / SqmSqq;
        t.element({10, 9, 0, 1, 1, 1, 1}) = (2. * SinTt * (rC * wm - 1) * WpmP) / SqmSqq;
        t.element({10, 10, 0, 1, 1, 1, 1}) = (2. * SinTt * (rC - wp) * WpmP) / SqmSqq;
        t.element({0, 2, 1, 0, 0, 0, 0}) = -(rt * SinTt * WpmP) / (2. * SqmSqq);
        t.element({0, 5, 1, 0, 0, 0, 0}) = (rt * SinTt * WpmM) / (2. * SqmSqq);
        t.element({5, 2, 1, 0, 0, 0, 0}) = -(rt * SinTt * WpmP) / (2. * SqmSqq);
        t.element({5, 5, 1, 0, 0, 0, 0}) = -(rt * SinTt * WpmM) / (2. * SqmSqq);
        t.element({7, 2, 1, 0, 0, 0, 0}) = -(rt * SinTt * WpmP) / (2. * SqmSqq);
        t.element({7, 5, 1, 0, 0, 0, 0}) = (rt * SinTt * WpmM) / (2. * SqmSqq);
        t.element({9, 8, 1, 0, 0, 0, 0}) = (-4. * SinTt * (rC * wm - 1) * sqwp) / SqmSqq;
        t.element({9, 9, 1, 0, 0, 0, 0}) = (-2. * SinTt * (rC * wm - 1) * WpmP) / SqmSqq;
        t.element({9, 10, 1, 0, 0, 0, 0}) = (-2. * SinTt * (rC - wp) * WpmP) / SqmSqq;
        t.element({0, 2, 1, 0, 0, 0, 1}) = SinTtHalfSq * WpmP;
        t.element({0, 5, 1, 0, 0, 0, 1}) = -(SinTtHalfSq * WpmM);
        t.element({5, 2, 1, 0, 0, 0, 1}) = SinTtHalfSq * WpmP;
        t.element({5, 5, 1, 0, 0, 0, 1}) = SinTtHalfSq * WpmM;
        t.element({7, 2, 1, 0, 0, 0, 1}) = SinTtHalfSq * WpmP;
        t.element({7, 5, 1, 0, 0, 0, 1}) = -(SinTtHalfSq * WpmM);
        t.element({9, 8, 1, 0, 0, 0, 1}) = (8 * rt * SinTtHalfSq * (rC * wm - 1) * sqwp) / mSqq;
        t.element({9, 9, 1, 0, 0, 0, 1}) = (4. * rt * SinTtHalfSq * (rC * wm - 1) * WpmP) / mSqq;
        t.element({9, 10, 1, 0, 0, 0, 1}) = (4. * rt * SinTtHalfSq * (rC - wp) * WpmP) / mSqq;
        t.element({6, 2, 1, 0, 1, 1, 0}) = -(CosTtHalfSq * WpmP);
        t.element({6, 5, 1, 0, 1, 1, 0}) = -(CosTtHalfSq * WpmM);
        t.element({8, 2, 1, 0, 1, 1, 0}) = -(CosTtHalfSq * WpmP);
        t.element({8, 5, 1, 0, 1, 1, 0}) = CosTtHalfSq * WpmM;
        t.element({10, 8, 1, 0, 1, 1, 0}) = (-8 * CosTtHalfSq * rt * sqwm * (rC * wp - 1)) / mSqq;
        t.element({10, 9, 1, 0, 1, 1, 0}) = (-4. * CosTtHalfSq * rt * (rC * wp - 1) * WpmP) / mSqq;
        t.element({10, 10, 1, 0, 1, 1, 0}) = (-4. * CosTtHalfSq * rt * (rC - wm) * WpmP) / mSqq;
        t.element({6, 2, 1, 0, 1, 1, 1}) = (rt * SinTt * WpmP) / (2. * SqmSqq);
        t.element({6, 5, 1, 0, 1, 1, 1}) = (rt * SinTt * WpmM) / (2. * SqmSqq);
        t.element({8, 2, 1, 0, 1, 1, 1}) = (rt * SinTt * WpmP) / (2. * SqmSqq);
        t.element({8, 5, 1, 0, 1, 1, 1}) = -(rt * SinTt * WpmM) / (2. * SqmSqq);
        t.element({10, 8, 1, 0, 1, 1, 1}) = (4. * SinTt * sqwm * (rC * wp - 1)) / SqmSqq;
        t.element({10, 9, 1, 0, 1, 1, 1}) = (2. * SinTt * (rC * wp - 1) * WpmP) / SqmSqq;
        t.element({10, 10, 1, 0, 1, 1, 1}) = (2. * SinTt * (rC - wm) * WpmP) / SqmSqq;
        t.element({0, 2, 1, 1, 0, 0, 0}) = -(Rp * rt) / (2. * mSqq);
        t.element({0, 3, 1, 1, 0, 0, 0}) = -(Ot * rt * WpmM) / (2. * mSqq);
        t.element({0, 4, 1, 1, 0, 0, 0}) = -(Op * rt * WpmM) / (2. * mSqq);
        t.element({0, 5, 1, 1, 0, 0, 0}) = (Rm * rt) / (2. * mSqq);
        t.element({0, 6, 1, 1, 0, 0, 0}) = (Ot * rt * WpmP) / (2. * mSqq);
        t.element({0, 7, 1, 1, 0, 0, 0}) = (Op * rt * WpmP) / (2. * mSqq);
        t.element({1, 0, 1, 1, 0, 0, 0}) = WpmM / 2.;
        t.element({1, 1, 1, 1, 0, 0, 0}) = -WpmP / 2.;
        t.element({3, 0, 1, 1, 0, 0, 0}) = WpmM / 2.;
        t.element({3, 1, 1, 1, 0, 0, 0}) = WpmP / 2.;
        t.element({5, 2, 1, 1, 0, 0, 0}) = -(Rp * rt) / (2. * mSqq);
        t.element({5, 3, 1, 1, 0, 0, 0}) = -(Ot * rt * WpmM) / (2. * mSqq);
        t.element({5, 4, 1, 1, 0, 0, 0}) = -(Op * rt * WpmM) / (2. * mSqq);
        t.element({5, 5, 1, 1, 0, 0, 0}) = -(Rm * rt) / (2. * mSqq);
        t.element({5, 6, 1, 1, 0, 0, 0}) = -(Ot * rt * WpmP) / (2. * mSqq);
        t.element({5, 7, 1, 1, 0, 0, 0}) = -(Op * rt * WpmP) / (2. * mSqq);
        t.element({7, 2, 1, 1, 0, 0, 0}) = -(Rp * rt) / (2. * mSqq);
        t.element({7, 3, 1, 1, 0, 0, 0}) = -(Ot * rt * WpmM) / (2. * mSqq);
        t.element({7, 4, 1, 1, 0, 0, 0}) = -(Op * rt * WpmM) / (2. * mSqq);
        t.element({7, 5, 1, 1, 0, 0, 0}) = (Rm * rt) / (2. * mSqq);
        t.element({7, 6, 1, 1, 0, 0, 0}) = (Ot * rt * WpmP) / (2. * mSqq);
        t.element({7, 7, 1, 1, 0, 0, 0}) = (Op * rt * WpmP) / (2. * mSqq);
        t.element({9, 8, 1, 1, 0, 0, 0}) = -4. * CosTt * sqwm;
        t.element({9, 9, 1, 1, 0, 0, 0}) = -2. * CosTt * WpmP;
        t.element({9, 10, 1, 1, 0, 0, 0}) = -2. * CosTt * WpmP;
        t.element({9, 11, 1, 1, 0, 0, 0}) = 2. * CosTt * Sqw2m1 * WpmM;
        t.element({0, 2, 1, 1, 0, 0, 1}) = ((-1 + rC) * SinTt * WpmP) / (2. * SqmSqq);
        t.element({0, 3, 1, 1, 0, 0, 1}) = (rC * SinTt * Sqw2m1 * WpmM) / (2. * SqmSqq);
        t.element({0, 4, 1, 1, 0, 0, 1}) = (SinTt * Sqw2m1 * WpmM) / (2. * SqmSqq);
        t.element({0, 5, 1, 1, 0, 0, 1}) = -((1 + rC) * SinTt * WpmM) / (2. * SqmSqq);
        t.element({0, 6, 1, 1, 0, 0, 1}) = -(rC * SinTt * Sqw2m1 * WpmP) / (2. * SqmSqq);
        t.element({0, 7, 1, 1, 0, 0, 1}) = -(SinTt * Sqw2m1 * WpmP) / (2. * SqmSqq);
        t.element({5, 2, 1, 1, 0, 0, 1}) = ((-1 + rC) * SinTt * WpmP) / (2. * SqmSqq);
        t.element({5, 3, 1, 1, 0, 0, 1}) = (rC * SinTt * Sqw2m1 * WpmM) / (2. * SqmSqq);
        t.element({5, 4, 1, 1, 0, 0, 1}) = (SinTt * Sqw2m1 * WpmM) / (2. * SqmSqq);
        t.element({5, 5, 1, 1, 0, 0, 1}) = ((1 + rC) * SinTt * WpmM) / (2. * SqmSqq);
        t.element({5, 6, 1, 1, 0, 0, 1}) = (rC * SinTt * Sqw2m1 * WpmP) / (2. * SqmSqq);
        t.element({5, 7, 1, 1, 0, 0, 1}) = (SinTt * Sqw2m1 * WpmP) / (2. * SqmSqq);
        t.element({7, 2, 1, 1, 0, 0, 1}) = ((-1 + rC) * SinTt * WpmP) / (2. * SqmSqq);
        t.element({7, 3, 1, 1, 0, 0, 1}) = (rC * SinTt * Sqw2m1 * WpmM) / (2. * SqmSqq);
        t.element({7, 4, 1, 1, 0, 0, 1}) = (SinTt * Sqw2m1 * WpmM) / (2. * SqmSqq);
        t.element({7, 5, 1, 1, 0, 0, 1}) = -((1 + rC) * SinTt * WpmM) / (2. * SqmSqq);
        t.element({7, 6, 1, 1, 0, 0, 1}) = -(rC * SinTt * Sqw2m1 * WpmP) / (2. * SqmSqq);
        t.element({7, 7, 1, 1, 0, 0, 1}) = -(SinTt * Sqw2m1 * WpmP) / (2. * SqmSqq);
        t.element({9, 8, 1, 1, 0, 0, 1}) = (4. * rt * SinTt * sqwm) / SqmSqq;
        t.element({9, 9, 1, 1, 0, 0, 1}) = (2. * rt * SinTt * WpmP) / SqmSqq;
        t.element({9, 10, 1, 1, 0, 0, 1}) = (2. * rt * SinTt * WpmP) / SqmSqq;
        t.element({9, 11, 1, 1, 0, 0, 1}) = (-2. * rt * SinTt * Sqw2m1 * WpmM) / SqmSqq;
        t.element({6, 2, 1, 1, 1, 1, 0}) = ((-1 + rC) * SinTt * WpmP) / (2. * SqmSqq);
        t.element({6, 3, 1, 1, 1, 1, 0}) = (rC * SinTt * Sqw2m1 * WpmM) / (2. * SqmSqq);
        t.element({6, 4, 1, 1, 1, 1, 0}) = (SinTt * Sqw2m1 * WpmM) / (2. * SqmSqq);
        t.element({6, 5, 1, 1, 1, 1, 0}) = ((1 + rC) * SinTt * WpmM) / (2. * SqmSqq);
        t.element({6, 6, 1, 1, 1, 1, 0}) = (rC * SinTt * Sqw2m1 * WpmP) / (2. * SqmSqq);
        t.element({6, 7, 1, 1, 1, 1, 0}) = (SinTt * Sqw2m1 * WpmP) / (2. * SqmSqq);
        t.element({8, 2, 1, 1, 1, 1, 0}) = ((-1 + rC) * SinTt * WpmP) / (2. * SqmSqq);
        t.element({8, 3, 1, 1, 1, 1, 0}) = (rC * SinTt * Sqw2m1 * WpmM) / (2. * SqmSqq);
        t.element({8, 4, 1, 1, 1, 1, 0}) = (SinTt * Sqw2m1 * WpmM) / (2. * SqmSqq);
        t.element({8, 5, 1, 1, 1, 1, 0}) = -((1 + rC) * SinTt * WpmM) / (2. * SqmSqq);
        t.element({8, 6, 1, 1, 1, 1, 0}) = -(rC * SinTt * Sqw2m1 * WpmP) / (2. * SqmSqq);
        t.element({8, 7, 1, 1, 1, 1, 0}) = -(SinTt * Sqw2m1 * WpmP) / (2. * SqmSqq);
        t.element({10, 8, 1, 1, 1, 1, 0}) = (4. * rt * SinTt * sqwp) / SqmSqq;
        t.element({10, 9, 1, 1, 1, 1, 0}) = (2. * rt * SinTt * WpmP) / SqmSqq;
        t.element({10, 10, 1, 1, 1, 1, 0}) = (2. * rt * SinTt * WpmP) / SqmSqq;
        t.element({10, 11, 1, 1, 1, 1, 0}) = (-2. * rt * SinTt * Sqw2m1 * WpmM) / SqmSqq;
        t.element({2, 0, 1, 1, 1, 1, 1}) = -WpmM / 2.;
        t.element({2, 1, 1, 1, 1, 1, 1}) = WpmP / 2.;
        t.element({4, 0, 1, 1, 1, 1, 1}) = -WpmM / 2.;
        t.element({4, 1, 1, 1, 1, 1, 1}) = -WpmP / 2.;
        t.element({6, 2, 1, 1, 1, 1, 1}) = (Rp * rt) / (2. * mSqq);
        t.element({6, 3, 1, 1, 1, 1, 1}) = (Ot * rt * WpmM) / (2. * mSqq);
        t.element({6, 4, 1, 1, 1, 1, 1}) = (Op * rt * WpmM) / (2. * mSqq);
        t.element({6, 5, 1, 1, 1, 1, 1}) = (Rm * rt) / (2. * mSqq);
        t.element({6, 6, 1, 1, 1, 1, 1}) = (Ot * rt * WpmP) / (2. * mSqq);
        t.element({6, 7, 1, 1, 1, 1, 1}) = (Op * rt * WpmP) / (2. * mSqq);
        t.element({8, 2, 1, 1, 1, 1, 1}) = (Rp * rt) / (2. * mSqq);
        t.element({8, 3, 1, 1, 1, 1, 1}) = (Ot * rt * WpmM) / (2. * mSqq);
        t.element({8, 4, 1, 1, 1, 1, 1}) = (Op * rt * WpmM) / (2. * mSqq);
        t.element({8, 5, 1, 1, 1, 1, 1}) = -(Rm * rt) / (2. * mSqq);
        t.element({8, 6, 1, 1, 1, 1, 1}) = -(Ot * rt * WpmP) / (2. * mSqq);
        t.element({8, 7, 1, 1, 1, 1, 1}) = -(Op * rt * WpmP) / (2. * mSqq);
        t.element({10, 8, 1, 1, 1, 1, 1}) = 4. * CosTt * sqwp;
        t.element({10, 9, 1, 1, 1, 1, 1}) = 2. * CosTt * WpmP;
        t.element({10, 10, 1, 1, 1, 1, 1}) = 2. * CosTt * WpmP;
        t.element({10, 11, 1, 1, 1, 1, 1}) = -2. * CosTt * Sqw2m1 * WpmM;

        t *= TwoSqTwoGfSq;
    }

    void AmplLbLcstar12LepNu::addRefs() const {
    }

} // namespace Hammer
