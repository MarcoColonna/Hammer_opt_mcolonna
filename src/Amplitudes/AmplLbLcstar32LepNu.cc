///
/// @file  AmplLbLcstar32LepNu.cc
/// @brief \f$ \Lambda_b \rightarrow L_c^*(2625) \tau\nu \f$ amplitude
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/Amplitudes/AmplLbLcstar32LepNu.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;
using namespace complex_literals;

namespace Hammer {

    namespace MD = MultiDimensional;

    AmplLbLcstar32LepNu::AmplLbLcstar32LepNu() {
        // Create tensor rank and dimensions
        IndexList dims{{11, 16, 2, 4, 2, 2, 2}};
        string name{"AmplLbLcstar32LepNu"};
        addProcessSignature(-PID::LAMBDAB, {PID::LAMBDACSTAR32MINUS, PID::NU_TAU, PID::ANTITAU});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCTAUNU, FF_LBLCSTAR32, SPIN_LB, SPIN_LCSTAR32,
                                                          SPIN_NUTAU, SPIN_NUTAU_REF, SPIN_TAUP})});

        addProcessSignature(-PID::LAMBDAB, {PID::LAMBDACSTAR32MINUS, PID::NU_MU, PID::ANTIMUON});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCMUNU, FF_LBLCSTAR32, SPIN_LB, SPIN_LCSTAR32,
                                                          SPIN_NUMU, SPIN_NUMU_REF, SPIN_MUP})});

        addProcessSignature(-PID::LAMBDAB, {PID::LAMBDACSTAR32MINUS, PID::NU_E, PID::POSITRON});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCENU, FF_LBLCSTAR32, SPIN_LB, SPIN_LCSTAR32, SPIN_NUE,
                                                          SPIN_NUE_REF, SPIN_EP})});

        setSignatureIndex();
        _multiplicity = 2ul;
    }

    // NOLINTBEGIN(readability-function-size)
    void AmplLbLcstar32LepNu::eval(const Particle& parent, const ParticleList& daughters,
                                   const ParticleList& /*references*/) {
        // Momenta
        const FourMomentum& pLbmes = parent.momentum();
        const FourMomentum& pLcstar32mes = daughters[0].momentum();
        const FourMomentum& kNuTau = daughters[1].momentum();
        const FourMomentum& pTau = daughters[2].momentum();

        // kinematic objects
        const double Mb = pLbmes.mass();
        const double Mb2 = Mb * Mb;
        const double Mc = pLcstar32mes.mass();
        const double Mc2 = Mc * Mc;
        const double Mt = pTau.mass();
        const double Sqq = Mb2 + Mc2 - 2. * (pLbmes * pLcstar32mes);
        const double Ew = (Mb2 - Mc2 + Sqq) / (2 * Mb);
        const double BNuTau = (pLbmes * kNuTau);
        const double NuTauQ = (pLbmes * kNuTau) - (pLcstar32mes * kNuTau);
        const double BQ = Mb2 - (pLbmes * pLcstar32mes);

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
        const double w2m1 = w * w - 1;
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
        const double Oz = Sqw2m1 + (rC - w) * CosTt;

        const double TwoSqTwoGfSq = TwoSqTwoGFermi * Mb2 * sqrt(mSqq - pow(rt, 2.)) * sqrt(rC);
        const double sqrt6 = sqrt2 * sqrt3;
        const double sqrt23 = sqrt2 / sqrt3;

        // initialize tensor elements to zero
        Tensor& t = getTensor();
        t.clearData();

        // set non-zero tensor elements
        t.element({0, 2, 0, 1, 0, 0, 0}) = -((rt * SinTt * Sqw2m1 * WpmM) / (sqrt6 * SqmSqq));
        t.element({0, 5, 0, 1, 0, 0, 0}) = (rt * SinTt * WpmP) / (2. * sqrt6 * SqmSqq);
        t.element({0, 6, 0, 1, 0, 0, 0}) = -((rt * SinTt * Sqw2m1 * WpmP) / (sqrt6 * SqmSqq));
        t.element({0, 9, 0, 1, 0, 0, 0}) = (rt * SinTt * WpmM) / (2. * sqrt6 * SqmSqq);
        t.element({5, 2, 0, 1, 0, 0, 0}) = -((rt * SinTt * Sqw2m1 * WpmM) / (sqrt6 * SqmSqq));
        t.element({5, 5, 0, 1, 0, 0, 0}) = (rt * SinTt * WpmP) / (2. * sqrt6 * SqmSqq);
        t.element({5, 6, 0, 1, 0, 0, 0}) = (rt * SinTt * Sqw2m1 * WpmP) / (sqrt6 * SqmSqq);
        t.element({5, 9, 0, 1, 0, 0, 0}) = -(rt * SinTt * WpmM) / (2. * sqrt6 * SqmSqq);
        t.element({7, 2, 0, 1, 0, 0, 0}) = -((rt * SinTt * Sqw2m1 * WpmM) / (sqrt6 * SqmSqq));
        t.element({7, 5, 0, 1, 0, 0, 0}) = (rt * SinTt * WpmP) / (2. * sqrt6 * SqmSqq);
        t.element({7, 6, 0, 1, 0, 0, 0}) = -((rt * SinTt * Sqw2m1 * WpmP) / (sqrt6 * SqmSqq));
        t.element({7, 9, 0, 1, 0, 0, 0}) = (rt * SinTt * WpmM) / (2. * sqrt6 * SqmSqq);
        t.element({9, 10, 0, 1, 0, 0, 0}) = (4. * sqrt23 * SinTt * Sqw2m1 * sqwm * (rC * wp - 1)) / SqmSqq;
        t.element({9, 11, 0, 1, 0, 0, 0}) = (2. * sqrt23 * SinTt * Sqw2m1 * (rC * wp - 1) * WpmM) / SqmSqq;
        t.element({9, 12, 0, 1, 0, 0, 0}) = (2. * sqrt23 * SinTt * Sqw2m1 * (rC - wm) * WpmM) / SqmSqq;
        t.element({9, 13, 0, 1, 0, 0, 0}) = (2. * sqrt23 * SinTt * (rC - wm) * sqwm) / SqmSqq;
        t.element({9, 14, 0, 1, 0, 0, 0}) = (sqrt23 * SinTt * (rC - wm) * WpmP) / SqmSqq;
        t.element({0, 2, 0, 1, 0, 0, 1}) = -(sqrt23 * CosTtHalfSq * Sqw2m1 * WpmM);
        t.element({0, 5, 0, 1, 0, 0, 1}) = (CosTtHalfSq * WpmP) / sqrt6;
        t.element({0, 6, 0, 1, 0, 0, 1}) = -(sqrt23 * CosTtHalfSq * Sqw2m1 * WpmP);
        t.element({0, 9, 0, 1, 0, 0, 1}) = (CosTtHalfSq * WpmM) / sqrt6;
        t.element({5, 2, 0, 1, 0, 0, 1}) = -(sqrt23 * CosTtHalfSq * Sqw2m1 * WpmM);
        t.element({5, 5, 0, 1, 0, 0, 1}) = (CosTtHalfSq * WpmP) / sqrt6;
        t.element({5, 6, 0, 1, 0, 0, 1}) = sqrt23 * CosTtHalfSq * Sqw2m1 * WpmP;
        t.element({5, 9, 0, 1, 0, 0, 1}) = -((CosTtHalfSq * WpmM) / sqrt6);
        t.element({7, 2, 0, 1, 0, 0, 1}) = -(sqrt23 * CosTtHalfSq * Sqw2m1 * WpmM);
        t.element({7, 5, 0, 1, 0, 0, 1}) = (CosTtHalfSq * WpmP) / sqrt6;
        t.element({7, 6, 0, 1, 0, 0, 1}) = -(sqrt23 * CosTtHalfSq * Sqw2m1 * WpmP);
        t.element({7, 9, 0, 1, 0, 0, 1}) = (CosTtHalfSq * WpmM) / sqrt6;
        t.element({9, 10, 0, 1, 0, 0, 1}) = (8 * sqrt23 * CosTtHalfSq * rt * Sqw2m1 * sqwm * (rC * wp - 1)) / mSqq;
        t.element({9, 11, 0, 1, 0, 0, 1}) = (4. * sqrt23 * CosTtHalfSq * rt * Sqw2m1 * (rC * wp - 1) * WpmM) / mSqq;
        t.element({9, 12, 0, 1, 0, 0, 1}) = (4. * sqrt23 * CosTtHalfSq * rt * Sqw2m1 * (rC - wm) * WpmM) / mSqq;
        t.element({9, 13, 0, 1, 0, 0, 1}) = (4. * sqrt23 * CosTtHalfSq * rt * (rC - wm) * sqwm) / mSqq;
        t.element({9, 14, 0, 1, 0, 0, 1}) = (2. * sqrt23 * CosTtHalfSq * rt * (rC - wm) * WpmP) / mSqq;
        t.element({6, 2, 0, 1, 1, 1, 0}) = sqrt23 * SinTtHalfSq * Sqw2m1 * WpmM;
        t.element({6, 5, 0, 1, 1, 1, 0}) = -((SinTtHalfSq * WpmP) / sqrt6);
        t.element({6, 6, 0, 1, 1, 1, 0}) = -(sqrt23 * SinTtHalfSq * Sqw2m1 * WpmP);
        t.element({6, 9, 0, 1, 1, 1, 0}) = (SinTtHalfSq * WpmM) / sqrt6;
        t.element({8, 2, 0, 1, 1, 1, 0}) = sqrt23 * SinTtHalfSq * Sqw2m1 * WpmM;
        t.element({8, 5, 0, 1, 1, 1, 0}) = -((SinTtHalfSq * WpmP) / sqrt6);
        t.element({8, 6, 0, 1, 1, 1, 0}) = sqrt23 * SinTtHalfSq * Sqw2m1 * WpmP;
        t.element({8, 9, 0, 1, 1, 1, 0}) = -((SinTtHalfSq * WpmM) / sqrt6);
        t.element({10, 10, 0, 1, 1, 1, 0}) = (8 * sqrt23 * rt * SinTtHalfSq * Sqw2m1 * (rC * wm - 1) * sqwp) / mSqq;
        t.element({10, 11, 0, 1, 1, 1, 0}) = (-4. * sqrt23 * rt * SinTtHalfSq * Sqw2m1 * (rC * wm - 1) * WpmM) / mSqq;
        t.element({10, 12, 0, 1, 1, 1, 0}) = (-4. * sqrt23 * rt * SinTtHalfSq * Sqw2m1 * (rC - wp) * WpmM) / mSqq;
        t.element({10, 13, 0, 1, 1, 1, 0}) = (4. * sqrt23 * rt * SinTtHalfSq * sqwp * (-rC + wp)) / mSqq;
        t.element({10, 14, 0, 1, 1, 1, 0}) = (-2. * sqrt23 * rt * SinTtHalfSq * (rC - wp) * WpmP) / mSqq;
        t.element({6, 2, 0, 1, 1, 1, 1}) = (rt * SinTt * Sqw2m1 * WpmM) / (sqrt6 * SqmSqq);
        t.element({6, 5, 0, 1, 1, 1, 1}) = -(rt * SinTt * WpmP) / (2. * sqrt6 * SqmSqq);
        t.element({6, 6, 0, 1, 1, 1, 1}) = -((rt * SinTt * Sqw2m1 * WpmP) / (sqrt6 * SqmSqq));
        t.element({6, 9, 0, 1, 1, 1, 1}) = (rt * SinTt * WpmM) / (2. * sqrt6 * SqmSqq);
        t.element({8, 2, 0, 1, 1, 1, 1}) = (rt * SinTt * Sqw2m1 * WpmM) / (sqrt6 * SqmSqq);
        t.element({8, 5, 0, 1, 1, 1, 1}) = -(rt * SinTt * WpmP) / (2. * sqrt6 * SqmSqq);
        t.element({8, 6, 0, 1, 1, 1, 1}) = (rt * SinTt * Sqw2m1 * WpmP) / (sqrt6 * SqmSqq);
        t.element({8, 9, 0, 1, 1, 1, 1}) = -(rt * SinTt * WpmM) / (2. * sqrt6 * SqmSqq);
        t.element({10, 10, 0, 1, 1, 1, 1}) = (4. * sqrt23 * SinTt * Sqw2m1 * (rC * wm - 1) * sqwp) / SqmSqq;
        t.element({10, 11, 0, 1, 1, 1, 1}) = (-2. * sqrt23 * SinTt * Sqw2m1 * (rC * wm - 1) * WpmM) / SqmSqq;
        t.element({10, 12, 0, 1, 1, 1, 1}) = (-2. * sqrt23 * SinTt * Sqw2m1 * (rC - wp) * WpmM) / SqmSqq;
        t.element({10, 13, 0, 1, 1, 1, 1}) = (2. * sqrt23 * SinTt * sqwp * (-rC + wp)) / SqmSqq;
        t.element({10, 14, 0, 1, 1, 1, 1}) = -((sqrt23 * SinTt * (rC - wp) * WpmP) / SqmSqq);
        t.element({0, 2, 0, 2, 0, 0, 0}) = (Rm * rt * Sqw2m1) / (sqrt6 * mSqq);
        t.element({0, 3, 0, 2, 0, 0, 0}) = (Ot * rt * Sqw2m1 * WpmP) / (sqrt6 * mSqq);
        t.element({0, 4, 0, 2, 0, 0, 0}) = (Op * rt * Sqw2m1 * WpmP) / (sqrt6 * mSqq);
        t.element({0, 5, 0, 2, 0, 0, 0}) = -((Oz * rt * WpmP) / (sqrt6 * mSqq));
        t.element({0, 6, 0, 2, 0, 0, 0}) = (Rp * rt * Sqw2m1) / (sqrt6 * mSqq);
        t.element({0, 7, 0, 2, 0, 0, 0}) = (Ot * rt * Sqw2m1 * WpmM) / (sqrt6 * mSqq);
        t.element({0, 8, 0, 2, 0, 0, 0}) = (Op * rt * Sqw2m1 * WpmM) / (sqrt6 * mSqq);
        t.element({0, 9, 0, 2, 0, 0, 0}) = -((Oz * rt * WpmM) / (sqrt6 * mSqq));
        t.element({1, 0, 0, 2, 0, 0, 0}) = (Sqw2m1 * WpmP) / sqrt6;
        t.element({1, 1, 0, 2, 0, 0, 0}) = (Sqw2m1 * WpmM) / sqrt6;
        t.element({3, 0, 0, 2, 0, 0, 0}) = (Sqw2m1 * WpmP) / sqrt6;
        t.element({3, 1, 0, 2, 0, 0, 0}) = -((Sqw2m1 * WpmM) / sqrt6);
        t.element({5, 2, 0, 2, 0, 0, 0}) = (Rm * rt * Sqw2m1) / (sqrt6 * mSqq);
        t.element({5, 3, 0, 2, 0, 0, 0}) = (Ot * rt * Sqw2m1 * WpmP) / (sqrt6 * mSqq);
        t.element({5, 4, 0, 2, 0, 0, 0}) = (Op * rt * Sqw2m1 * WpmP) / (sqrt6 * mSqq);
        t.element({5, 5, 0, 2, 0, 0, 0}) = -((Oz * rt * WpmP) / (sqrt6 * mSqq));
        t.element({5, 6, 0, 2, 0, 0, 0}) = -((Rp * rt * Sqw2m1) / (sqrt6 * mSqq));
        t.element({5, 7, 0, 2, 0, 0, 0}) = -((Ot * rt * Sqw2m1 * WpmM) / (sqrt6 * mSqq));
        t.element({5, 8, 0, 2, 0, 0, 0}) = -((Op * rt * Sqw2m1 * WpmM) / (sqrt6 * mSqq));
        t.element({5, 9, 0, 2, 0, 0, 0}) = (Oz * rt * WpmM) / (sqrt6 * mSqq);
        t.element({7, 2, 0, 2, 0, 0, 0}) = (Rm * rt * Sqw2m1) / (sqrt6 * mSqq);
        t.element({7, 3, 0, 2, 0, 0, 0}) = (Ot * rt * Sqw2m1 * WpmP) / (sqrt6 * mSqq);
        t.element({7, 4, 0, 2, 0, 0, 0}) = (Op * rt * Sqw2m1 * WpmP) / (sqrt6 * mSqq);
        t.element({7, 5, 0, 2, 0, 0, 0}) = -((Oz * rt * WpmP) / (sqrt6 * mSqq));
        t.element({7, 6, 0, 2, 0, 0, 0}) = (Rp * rt * Sqw2m1) / (sqrt6 * mSqq);
        t.element({7, 7, 0, 2, 0, 0, 0}) = (Ot * rt * Sqw2m1 * WpmM) / (sqrt6 * mSqq);
        t.element({7, 8, 0, 2, 0, 0, 0}) = (Op * rt * Sqw2m1 * WpmM) / (sqrt6 * mSqq);
        t.element({7, 9, 0, 2, 0, 0, 0}) = -((Oz * rt * WpmM) / (sqrt6 * mSqq));
        t.element({9, 10, 0, 2, 0, 0, 0}) = -4. * sqrt23 * CosTt * Sqw2m1 * sqwp;
        t.element({9, 11, 0, 2, 0, 0, 0}) = 2. * sqrt23 * CosTt * Sqw2m1 * WpmM;
        t.element({9, 12, 0, 2, 0, 0, 0}) = -2. * sqrt23 * CosTt * Sqw2m1 * WpmM;
        t.element({9, 13, 0, 2, 0, 0, 0}) = -4. * sqrt23 * CosTt * sqwm;
        t.element({9, 14, 0, 2, 0, 0, 0}) = -2. * sqrt23 * CosTt * WpmP;
        t.element({9, 15, 0, 2, 0, 0, 0}) = 2. * sqrt23 * CosTt * w2m1 * WpmP;
        t.element({0, 2, 0, 2, 0, 0, 1}) = -(((1 + rC) * SinTt * Sqw2m1 * WpmM) / (sqrt6 * SqmSqq));
        t.element({0, 3, 0, 2, 0, 0, 1}) = -((rC * SinTt * w2m1 * WpmP) / (sqrt6 * SqmSqq));
        t.element({0, 4, 0, 2, 0, 0, 1}) = -((SinTt * w2m1 * WpmP) / (sqrt6 * SqmSqq));
        t.element({0, 5, 0, 2, 0, 0, 1}) = (SinTt * (rC - w) * WpmP) / (sqrt6 * SqmSqq);
        t.element({0, 6, 0, 2, 0, 0, 1}) = -(((-1 + rC) * SinTt * Sqw2m1 * WpmP) / (sqrt6 * SqmSqq));
        t.element({0, 7, 0, 2, 0, 0, 1}) = -((rC * SinTt * w2m1 * WpmM) / (sqrt6 * SqmSqq));
        t.element({0, 8, 0, 2, 0, 0, 1}) = -((SinTt * w2m1 * WpmM) / (sqrt6 * SqmSqq));
        t.element({0, 9, 0, 2, 0, 0, 1}) = (SinTt * (rC - w) * WpmM) / (sqrt6 * SqmSqq);
        t.element({5, 2, 0, 2, 0, 0, 1}) = -(((1 + rC) * SinTt * Sqw2m1 * WpmM) / (sqrt6 * SqmSqq));
        t.element({5, 3, 0, 2, 0, 0, 1}) = -((rC * SinTt * w2m1 * WpmP) / (sqrt6 * SqmSqq));
        t.element({5, 4, 0, 2, 0, 0, 1}) = -((SinTt * w2m1 * WpmP) / (sqrt6 * SqmSqq));
        t.element({5, 5, 0, 2, 0, 0, 1}) = (SinTt * (rC - w) * WpmP) / (sqrt6 * SqmSqq);
        t.element({5, 6, 0, 2, 0, 0, 1}) = ((-1 + rC) * SinTt * Sqw2m1 * WpmP) / (sqrt6 * SqmSqq);
        t.element({5, 7, 0, 2, 0, 0, 1}) = (rC * SinTt * w2m1 * WpmM) / (sqrt6 * SqmSqq);
        t.element({5, 8, 0, 2, 0, 0, 1}) = (SinTt * w2m1 * WpmM) / (sqrt6 * SqmSqq);
        t.element({5, 9, 0, 2, 0, 0, 1}) = -((SinTt * (rC - w) * WpmM) / (sqrt6 * SqmSqq));
        t.element({7, 2, 0, 2, 0, 0, 1}) = -(((1 + rC) * SinTt * Sqw2m1 * WpmM) / (sqrt6 * SqmSqq));
        t.element({7, 3, 0, 2, 0, 0, 1}) = -((rC * SinTt * w2m1 * WpmP) / (sqrt6 * SqmSqq));
        t.element({7, 4, 0, 2, 0, 0, 1}) = -((SinTt * w2m1 * WpmP) / (sqrt6 * SqmSqq));
        t.element({7, 5, 0, 2, 0, 0, 1}) = (SinTt * (rC - w) * WpmP) / (sqrt6 * SqmSqq);
        t.element({7, 6, 0, 2, 0, 0, 1}) = -(((-1 + rC) * SinTt * Sqw2m1 * WpmP) / (sqrt6 * SqmSqq));
        t.element({7, 7, 0, 2, 0, 0, 1}) = -((rC * SinTt * w2m1 * WpmM) / (sqrt6 * SqmSqq));
        t.element({7, 8, 0, 2, 0, 0, 1}) = -((SinTt * w2m1 * WpmM) / (sqrt6 * SqmSqq));
        t.element({7, 9, 0, 2, 0, 0, 1}) = (SinTt * (rC - w) * WpmM) / (sqrt6 * SqmSqq);
        t.element({9, 10, 0, 2, 0, 0, 1}) = (4. * sqrt23 * rt * SinTt * Sqw2m1 * sqwp) / SqmSqq;
        t.element({9, 11, 0, 2, 0, 0, 1}) = (-2. * sqrt23 * rt * SinTt * Sqw2m1 * WpmM) / SqmSqq;
        t.element({9, 12, 0, 2, 0, 0, 1}) = (2. * sqrt23 * rt * SinTt * Sqw2m1 * WpmM) / SqmSqq;
        t.element({9, 13, 0, 2, 0, 0, 1}) = (4. * sqrt23 * rt * SinTt * sqwm) / SqmSqq;
        t.element({9, 14, 0, 2, 0, 0, 1}) = (2. * sqrt23 * rt * SinTt * WpmP) / SqmSqq;
        t.element({9, 15, 0, 2, 0, 0, 1}) = (-2. * sqrt23 * rt * SinTt * w2m1 * WpmP) / SqmSqq;
        t.element({6, 2, 0, 2, 1, 1, 0}) = -(((1 + rC) * SinTt * Sqw2m1 * WpmM) / (sqrt6 * SqmSqq));
        t.element({6, 3, 0, 2, 1, 1, 0}) = -((rC * SinTt * w2m1 * WpmP) / (sqrt6 * SqmSqq));
        t.element({6, 4, 0, 2, 1, 1, 0}) = -((SinTt * w2m1 * WpmP) / (sqrt6 * SqmSqq));
        t.element({6, 5, 0, 2, 1, 1, 0}) = (SinTt * (rC - w) * WpmP) / (sqrt6 * SqmSqq);
        t.element({6, 6, 0, 2, 1, 1, 0}) = ((-1 + rC) * SinTt * Sqw2m1 * WpmP) / (sqrt6 * SqmSqq);
        t.element({6, 7, 0, 2, 1, 1, 0}) = (rC * SinTt * w2m1 * WpmM) / (sqrt6 * SqmSqq);
        t.element({6, 8, 0, 2, 1, 1, 0}) = (SinTt * w2m1 * WpmM) / (sqrt6 * SqmSqq);
        t.element({6, 9, 0, 2, 1, 1, 0}) = -((SinTt * (rC - w) * WpmM) / (sqrt6 * SqmSqq));
        t.element({8, 2, 0, 2, 1, 1, 0}) = -(((1 + rC) * SinTt * Sqw2m1 * WpmM) / (sqrt6 * SqmSqq));
        t.element({8, 3, 0, 2, 1, 1, 0}) = -((rC * SinTt * w2m1 * WpmP) / (sqrt6 * SqmSqq));
        t.element({8, 4, 0, 2, 1, 1, 0}) = -((SinTt * w2m1 * WpmP) / (sqrt6 * SqmSqq));
        t.element({8, 5, 0, 2, 1, 1, 0}) = (SinTt * (rC - w) * WpmP) / (sqrt6 * SqmSqq);
        t.element({8, 6, 0, 2, 1, 1, 0}) = -(((-1 + rC) * SinTt * Sqw2m1 * WpmP) / (sqrt6 * SqmSqq));
        t.element({8, 7, 0, 2, 1, 1, 0}) = -((rC * SinTt * w2m1 * WpmM) / (sqrt6 * SqmSqq));
        t.element({8, 8, 0, 2, 1, 1, 0}) = -((SinTt * w2m1 * WpmM) / (sqrt6 * SqmSqq));
        t.element({8, 9, 0, 2, 1, 1, 0}) = (SinTt * (rC - w) * WpmM) / (sqrt6 * SqmSqq);
        t.element({10, 10, 0, 2, 1, 1, 0}) = (-4. * sqrt23 * rt * SinTt * Sqw2m1 * sqwm) / SqmSqq;
        t.element({10, 11, 0, 2, 1, 1, 0}) = (-2. * sqrt23 * rt * SinTt * Sqw2m1 * WpmM) / SqmSqq;
        t.element({10, 12, 0, 2, 1, 1, 0}) = (2. * sqrt23 * rt * SinTt * Sqw2m1 * WpmM) / SqmSqq;
        t.element({10, 13, 0, 2, 1, 1, 0}) = (4. * sqrt23 * rt * SinTt * sqwp) / SqmSqq;
        t.element({10, 14, 0, 2, 1, 1, 0}) = (2. * sqrt23 * rt * SinTt * WpmP) / SqmSqq;
        t.element({10, 15, 0, 2, 1, 1, 0}) = (-2. * sqrt23 * rt * SinTt * w2m1 * WpmP) / SqmSqq;
        t.element({2, 0, 0, 2, 1, 1, 1}) = -((Sqw2m1 * WpmP) / sqrt6);
        t.element({2, 1, 0, 2, 1, 1, 1}) = -((Sqw2m1 * WpmM) / sqrt6);
        t.element({4, 0, 0, 2, 1, 1, 1}) = -((Sqw2m1 * WpmP) / sqrt6);
        t.element({4, 1, 0, 2, 1, 1, 1}) = (Sqw2m1 * WpmM) / sqrt6;
        t.element({6, 2, 0, 2, 1, 1, 1}) = -((Rm * rt * Sqw2m1) / (sqrt6 * mSqq));
        t.element({6, 3, 0, 2, 1, 1, 1}) = -((Ot * rt * Sqw2m1 * WpmP) / (sqrt6 * mSqq));
        t.element({6, 4, 0, 2, 1, 1, 1}) = -((Op * rt * Sqw2m1 * WpmP) / (sqrt6 * mSqq));
        t.element({6, 5, 0, 2, 1, 1, 1}) = (Oz * rt * WpmP) / (sqrt6 * mSqq);
        t.element({6, 6, 0, 2, 1, 1, 1}) = (Rp * rt * Sqw2m1) / (sqrt6 * mSqq);
        t.element({6, 7, 0, 2, 1, 1, 1}) = (Ot * rt * Sqw2m1 * WpmM) / (sqrt6 * mSqq);
        t.element({6, 8, 0, 2, 1, 1, 1}) = (Op * rt * Sqw2m1 * WpmM) / (sqrt6 * mSqq);
        t.element({6, 9, 0, 2, 1, 1, 1}) = -((Oz * rt * WpmM) / (sqrt6 * mSqq));
        t.element({8, 2, 0, 2, 1, 1, 1}) = -((Rm * rt * Sqw2m1) / (sqrt6 * mSqq));
        t.element({8, 3, 0, 2, 1, 1, 1}) = -((Ot * rt * Sqw2m1 * WpmP) / (sqrt6 * mSqq));
        t.element({8, 4, 0, 2, 1, 1, 1}) = -((Op * rt * Sqw2m1 * WpmP) / (sqrt6 * mSqq));
        t.element({8, 5, 0, 2, 1, 1, 1}) = (Oz * rt * WpmP) / (sqrt6 * mSqq);
        t.element({8, 6, 0, 2, 1, 1, 1}) = -((Rp * rt * Sqw2m1) / (sqrt6 * mSqq));
        t.element({8, 7, 0, 2, 1, 1, 1}) = -((Ot * rt * Sqw2m1 * WpmM) / (sqrt6 * mSqq));
        t.element({8, 8, 0, 2, 1, 1, 1}) = -((Op * rt * Sqw2m1 * WpmM) / (sqrt6 * mSqq));
        t.element({8, 9, 0, 2, 1, 1, 1}) = (Oz * rt * WpmM) / (sqrt6 * mSqq);
        t.element({10, 10, 0, 2, 1, 1, 1}) = -4. * sqrt23 * CosTt * Sqw2m1 * sqwm;
        t.element({10, 11, 0, 2, 1, 1, 1}) = -2. * sqrt23 * CosTt * Sqw2m1 * WpmM;
        t.element({10, 12, 0, 2, 1, 1, 1}) = 2. * sqrt23 * CosTt * Sqw2m1 * WpmM;
        t.element({10, 13, 0, 2, 1, 1, 1}) = 4. * sqrt23 * CosTt * sqwp;
        t.element({10, 14, 0, 2, 1, 1, 1}) = 2. * sqrt23 * CosTt * WpmP;
        t.element({10, 15, 0, 2, 1, 1, 1}) = -2. * sqrt23 * CosTt * w2m1 * WpmP;
        t.element({0, 5, 0, 3, 0, 0, 0}) = -(rt * SinTt * WpmP) / (2. * sqrt2 * SqmSqq);
        t.element({0, 9, 0, 3, 0, 0, 0}) = -(rt * SinTt * WpmM) / (2. * sqrt2 * SqmSqq);
        t.element({5, 5, 0, 3, 0, 0, 0}) = -(rt * SinTt * WpmP) / (2. * sqrt2 * SqmSqq);
        t.element({5, 9, 0, 3, 0, 0, 0}) = (rt * SinTt * WpmM) / (2. * sqrt2 * SqmSqq);
        t.element({7, 5, 0, 3, 0, 0, 0}) = -(rt * SinTt * WpmP) / (2. * sqrt2 * SqmSqq);
        t.element({7, 9, 0, 3, 0, 0, 0}) = -(rt * SinTt * WpmM) / (2. * sqrt2 * SqmSqq);
        t.element({9, 13, 0, 3, 0, 0, 0}) = (-2. * sqrt2 * SinTt * (rC * wm - 1) * sqwp) / SqmSqq;
        t.element({9, 14, 0, 3, 0, 0, 0}) = -((sqrt2 * SinTt * (rC - wp) * WpmP) / SqmSqq);
        t.element({0, 5, 0, 3, 0, 0, 1}) = (SinTtHalfSq * WpmP) / sqrt2;
        t.element({0, 9, 0, 3, 0, 0, 1}) = (SinTtHalfSq * WpmM) / sqrt2;
        t.element({5, 5, 0, 3, 0, 0, 1}) = (SinTtHalfSq * WpmP) / sqrt2;
        t.element({5, 9, 0, 3, 0, 0, 1}) = -((SinTtHalfSq * WpmM) / sqrt2);
        t.element({7, 5, 0, 3, 0, 0, 1}) = (SinTtHalfSq * WpmP) / sqrt2;
        t.element({7, 9, 0, 3, 0, 0, 1}) = (SinTtHalfSq * WpmM) / sqrt2;
        t.element({9, 13, 0, 3, 0, 0, 1}) = (4. * sqrt2 * rt * SinTtHalfSq * (rC * wm - 1) * sqwp) / mSqq;
        t.element({9, 14, 0, 3, 0, 0, 1}) = (2. * sqrt2 * rt * SinTtHalfSq * (rC - wp) * WpmP) / mSqq;
        t.element({6, 5, 0, 3, 1, 1, 0}) = -((CosTtHalfSq * WpmP) / sqrt2);
        t.element({6, 9, 0, 3, 1, 1, 0}) = (CosTtHalfSq * WpmM) / sqrt2;
        t.element({8, 5, 0, 3, 1, 1, 0}) = -((CosTtHalfSq * WpmP) / sqrt2);
        t.element({8, 9, 0, 3, 1, 1, 0}) = -((CosTtHalfSq * WpmM) / sqrt2);
        t.element({10, 13, 0, 3, 1, 1, 0}) = (-4. * sqrt2 * CosTtHalfSq * rt * sqwm * (rC * wp - 1)) / mSqq;
        t.element({10, 14, 0, 3, 1, 1, 0}) = (-2. * sqrt2 * CosTtHalfSq * rt * (rC - wm) * WpmP) / mSqq;
        t.element({6, 5, 0, 3, 1, 1, 1}) = (rt * SinTt * WpmP) / (2. * sqrt2 * SqmSqq);
        t.element({6, 9, 0, 3, 1, 1, 1}) = -(rt * SinTt * WpmM) / (2. * sqrt2 * SqmSqq);
        t.element({8, 5, 0, 3, 1, 1, 1}) = (rt * SinTt * WpmP) / (2. * sqrt2 * SqmSqq);
        t.element({8, 9, 0, 3, 1, 1, 1}) = (rt * SinTt * WpmM) / (2. * sqrt2 * SqmSqq);
        t.element({10, 13, 0, 3, 1, 1, 1}) = (2. * sqrt2 * SinTt * sqwm * (rC * wp - 1)) / SqmSqq;
        t.element({10, 14, 0, 3, 1, 1, 1}) = (sqrt2 * SinTt * (rC - wm) * WpmP) / SqmSqq;
        t.element({0, 5, 1, 0, 0, 0, 0}) = -(rt * SinTt * WpmP) / (2. * sqrt2 * SqmSqq);
        t.element({0, 9, 1, 0, 0, 0, 0}) = (rt * SinTt * WpmM) / (2. * sqrt2 * SqmSqq);
        t.element({5, 5, 1, 0, 0, 0, 0}) = -(rt * SinTt * WpmP) / (2. * sqrt2 * SqmSqq);
        t.element({5, 9, 1, 0, 0, 0, 0}) = -(rt * SinTt * WpmM) / (2. * sqrt2 * SqmSqq);
        t.element({7, 5, 1, 0, 0, 0, 0}) = -(rt * SinTt * WpmP) / (2. * sqrt2 * SqmSqq);
        t.element({7, 9, 1, 0, 0, 0, 0}) = (rt * SinTt * WpmM) / (2. * sqrt2 * SqmSqq);
        t.element({9, 13, 1, 0, 0, 0, 0}) = (-2. * sqrt2 * SinTt * sqwm * (rC * wp - 1)) / SqmSqq;
        t.element({9, 14, 1, 0, 0, 0, 0}) = -((sqrt2 * SinTt * (rC - wm) * WpmP) / SqmSqq);
        t.element({0, 5, 1, 0, 0, 0, 1}) = -((CosTtHalfSq * WpmP) / sqrt2);
        t.element({0, 9, 1, 0, 0, 0, 1}) = (CosTtHalfSq * WpmM) / sqrt2;
        t.element({5, 5, 1, 0, 0, 0, 1}) = -((CosTtHalfSq * WpmP) / sqrt2);
        t.element({5, 9, 1, 0, 0, 0, 1}) = -((CosTtHalfSq * WpmM) / sqrt2);
        t.element({7, 5, 1, 0, 0, 0, 1}) = -((CosTtHalfSq * WpmP) / sqrt2);
        t.element({7, 9, 1, 0, 0, 0, 1}) = (CosTtHalfSq * WpmM) / sqrt2;
        t.element({9, 13, 1, 0, 0, 0, 1}) = (-4. * sqrt2 * CosTtHalfSq * rt * sqwm * (rC * wp - 1)) / mSqq;
        t.element({9, 14, 1, 0, 0, 0, 1}) = (-2. * sqrt2 * CosTtHalfSq * rt * (rC - wm) * WpmP) / mSqq;
        t.element({6, 5, 1, 0, 1, 1, 0}) = (SinTtHalfSq * WpmP) / sqrt2;
        t.element({6, 9, 1, 0, 1, 1, 0}) = (SinTtHalfSq * WpmM) / sqrt2;
        t.element({8, 5, 1, 0, 1, 1, 0}) = (SinTtHalfSq * WpmP) / sqrt2;
        t.element({8, 9, 1, 0, 1, 1, 0}) = -((SinTtHalfSq * WpmM) / sqrt2);
        t.element({10, 13, 1, 0, 1, 1, 0}) = (4. * sqrt2 * rt * SinTtHalfSq * (rC * wm - 1) * sqwp) / mSqq;
        t.element({10, 14, 1, 0, 1, 1, 0}) = (2. * sqrt2 * rt * SinTtHalfSq * (rC - wp) * WpmP) / mSqq;
        t.element({6, 5, 1, 0, 1, 1, 1}) = (rt * SinTt * WpmP) / (2. * sqrt2 * SqmSqq);
        t.element({6, 9, 1, 0, 1, 1, 1}) = (rt * SinTt * WpmM) / (2. * sqrt2 * SqmSqq);
        t.element({8, 5, 1, 0, 1, 1, 1}) = (rt * SinTt * WpmP) / (2. * sqrt2 * SqmSqq);
        t.element({8, 9, 1, 0, 1, 1, 1}) = -(rt * SinTt * WpmM) / (2. * sqrt2 * SqmSqq);
        t.element({10, 13, 1, 0, 1, 1, 1}) = (2. * sqrt2 * SinTt * (rC * wm - 1) * sqwp) / SqmSqq;
        t.element({10, 14, 1, 0, 1, 1, 1}) = (sqrt2 * SinTt * (rC - wp) * WpmP) / SqmSqq;
        t.element({0, 2, 1, 1, 0, 0, 0}) = (Rm * rt * Sqw2m1) / (sqrt6 * mSqq);
        t.element({0, 3, 1, 1, 0, 0, 0}) = (Ot * rt * Sqw2m1 * WpmP) / (sqrt6 * mSqq);
        t.element({0, 4, 1, 1, 0, 0, 0}) = (Op * rt * Sqw2m1 * WpmP) / (sqrt6 * mSqq);
        t.element({0, 5, 1, 1, 0, 0, 0}) = -((Oz * rt * WpmP) / (sqrt6 * mSqq));
        t.element({0, 6, 1, 1, 0, 0, 0}) = -((Rp * rt * Sqw2m1) / (sqrt6 * mSqq));
        t.element({0, 7, 1, 1, 0, 0, 0}) = -((Ot * rt * Sqw2m1 * WpmM) / (sqrt6 * mSqq));
        t.element({0, 8, 1, 1, 0, 0, 0}) = -((Op * rt * Sqw2m1 * WpmM) / (sqrt6 * mSqq));
        t.element({0, 9, 1, 1, 0, 0, 0}) = (Oz * rt * WpmM) / (sqrt6 * mSqq);
        t.element({1, 0, 1, 1, 0, 0, 0}) = (Sqw2m1 * WpmP) / sqrt6;
        t.element({1, 1, 1, 1, 0, 0, 0}) = -((Sqw2m1 * WpmM) / sqrt6);
        t.element({3, 0, 1, 1, 0, 0, 0}) = (Sqw2m1 * WpmP) / sqrt6;
        t.element({3, 1, 1, 1, 0, 0, 0}) = (Sqw2m1 * WpmM) / sqrt6;
        t.element({5, 2, 1, 1, 0, 0, 0}) = (Rm * rt * Sqw2m1) / (sqrt6 * mSqq);
        t.element({5, 3, 1, 1, 0, 0, 0}) = (Ot * rt * Sqw2m1 * WpmP) / (sqrt6 * mSqq);
        t.element({5, 4, 1, 1, 0, 0, 0}) = (Op * rt * Sqw2m1 * WpmP) / (sqrt6 * mSqq);
        t.element({5, 5, 1, 1, 0, 0, 0}) = -((Oz * rt * WpmP) / (sqrt6 * mSqq));
        t.element({5, 6, 1, 1, 0, 0, 0}) = (Rp * rt * Sqw2m1) / (sqrt6 * mSqq);
        t.element({5, 7, 1, 1, 0, 0, 0}) = (Ot * rt * Sqw2m1 * WpmM) / (sqrt6 * mSqq);
        t.element({5, 8, 1, 1, 0, 0, 0}) = (Op * rt * Sqw2m1 * WpmM) / (sqrt6 * mSqq);
        t.element({5, 9, 1, 1, 0, 0, 0}) = -((Oz * rt * WpmM) / (sqrt6 * mSqq));
        t.element({7, 2, 1, 1, 0, 0, 0}) = (Rm * rt * Sqw2m1) / (sqrt6 * mSqq);
        t.element({7, 3, 1, 1, 0, 0, 0}) = (Ot * rt * Sqw2m1 * WpmP) / (sqrt6 * mSqq);
        t.element({7, 4, 1, 1, 0, 0, 0}) = (Op * rt * Sqw2m1 * WpmP) / (sqrt6 * mSqq);
        t.element({7, 5, 1, 1, 0, 0, 0}) = -((Oz * rt * WpmP) / (sqrt6 * mSqq));
        t.element({7, 6, 1, 1, 0, 0, 0}) = -((Rp * rt * Sqw2m1) / (sqrt6 * mSqq));
        t.element({7, 7, 1, 1, 0, 0, 0}) = -((Ot * rt * Sqw2m1 * WpmM) / (sqrt6 * mSqq));
        t.element({7, 8, 1, 1, 0, 0, 0}) = -((Op * rt * Sqw2m1 * WpmM) / (sqrt6 * mSqq));
        t.element({7, 9, 1, 1, 0, 0, 0}) = (Oz * rt * WpmM) / (sqrt6 * mSqq);
        t.element({9, 10, 1, 1, 0, 0, 0}) = 4. * sqrt23 * CosTt * Sqw2m1 * sqwm;
        t.element({9, 11, 1, 1, 0, 0, 0}) = 2. * sqrt23 * CosTt * Sqw2m1 * WpmM;
        t.element({9, 12, 1, 1, 0, 0, 0}) = -2. * sqrt23 * CosTt * Sqw2m1 * WpmM;
        t.element({9, 13, 1, 1, 0, 0, 0}) = -4. * sqrt23 * CosTt * sqwp;
        t.element({9, 14, 1, 1, 0, 0, 0}) = -2. * sqrt23 * CosTt * WpmP;
        t.element({9, 15, 1, 1, 0, 0, 0}) = 2. * sqrt23 * CosTt * w2m1 * WpmP;
        t.element({0, 2, 1, 1, 0, 0, 1}) = -(((1 + rC) * SinTt * Sqw2m1 * WpmM) / (sqrt6 * SqmSqq));
        t.element({0, 3, 1, 1, 0, 0, 1}) = -((rC * SinTt * w2m1 * WpmP) / (sqrt6 * SqmSqq));
        t.element({0, 4, 1, 1, 0, 0, 1}) = -((SinTt * w2m1 * WpmP) / (sqrt6 * SqmSqq));
        t.element({0, 5, 1, 1, 0, 0, 1}) = (SinTt * (rC - w) * WpmP) / (sqrt6 * SqmSqq);
        t.element({0, 6, 1, 1, 0, 0, 1}) = ((-1 + rC) * SinTt * Sqw2m1 * WpmP) / (sqrt6 * SqmSqq);
        t.element({0, 7, 1, 1, 0, 0, 1}) = (rC * SinTt * w2m1 * WpmM) / (sqrt6 * SqmSqq);
        t.element({0, 8, 1, 1, 0, 0, 1}) = (SinTt * w2m1 * WpmM) / (sqrt6 * SqmSqq);
        t.element({0, 9, 1, 1, 0, 0, 1}) = -((SinTt * (rC - w) * WpmM) / (sqrt6 * SqmSqq));
        t.element({5, 2, 1, 1, 0, 0, 1}) = -(((1 + rC) * SinTt * Sqw2m1 * WpmM) / (sqrt6 * SqmSqq));
        t.element({5, 3, 1, 1, 0, 0, 1}) = -((rC * SinTt * w2m1 * WpmP) / (sqrt6 * SqmSqq));
        t.element({5, 4, 1, 1, 0, 0, 1}) = -((SinTt * w2m1 * WpmP) / (sqrt6 * SqmSqq));
        t.element({5, 5, 1, 1, 0, 0, 1}) = (SinTt * (rC - w) * WpmP) / (sqrt6 * SqmSqq);
        t.element({5, 6, 1, 1, 0, 0, 1}) = -(((-1 + rC) * SinTt * Sqw2m1 * WpmP) / (sqrt6 * SqmSqq));
        t.element({5, 7, 1, 1, 0, 0, 1}) = -((rC * SinTt * w2m1 * WpmM) / (sqrt6 * SqmSqq));
        t.element({5, 8, 1, 1, 0, 0, 1}) = -((SinTt * w2m1 * WpmM) / (sqrt6 * SqmSqq));
        t.element({5, 9, 1, 1, 0, 0, 1}) = (SinTt * (rC - w) * WpmM) / (sqrt6 * SqmSqq);
        t.element({7, 2, 1, 1, 0, 0, 1}) = -(((1 + rC) * SinTt * Sqw2m1 * WpmM) / (sqrt6 * SqmSqq));
        t.element({7, 3, 1, 1, 0, 0, 1}) = -((rC * SinTt * w2m1 * WpmP) / (sqrt6 * SqmSqq));
        t.element({7, 4, 1, 1, 0, 0, 1}) = -((SinTt * w2m1 * WpmP) / (sqrt6 * SqmSqq));
        t.element({7, 5, 1, 1, 0, 0, 1}) = (SinTt * (rC - w) * WpmP) / (sqrt6 * SqmSqq);
        t.element({7, 6, 1, 1, 0, 0, 1}) = ((-1 + rC) * SinTt * Sqw2m1 * WpmP) / (sqrt6 * SqmSqq);
        t.element({7, 7, 1, 1, 0, 0, 1}) = (rC * SinTt * w2m1 * WpmM) / (sqrt6 * SqmSqq);
        t.element({7, 8, 1, 1, 0, 0, 1}) = (SinTt * w2m1 * WpmM) / (sqrt6 * SqmSqq);
        t.element({7, 9, 1, 1, 0, 0, 1}) = -((SinTt * (rC - w) * WpmM) / (sqrt6 * SqmSqq));
        t.element({9, 10, 1, 1, 0, 0, 1}) = (-4. * sqrt23 * rt * SinTt * Sqw2m1 * sqwm) / SqmSqq;
        t.element({9, 11, 1, 1, 0, 0, 1}) = (-2. * sqrt23 * rt * SinTt * Sqw2m1 * WpmM) / SqmSqq;
        t.element({9, 12, 1, 1, 0, 0, 1}) = (2. * sqrt23 * rt * SinTt * Sqw2m1 * WpmM) / SqmSqq;
        t.element({9, 13, 1, 1, 0, 0, 1}) = (4. * sqrt23 * rt * SinTt * sqwp) / SqmSqq;
        t.element({9, 14, 1, 1, 0, 0, 1}) = (2. * sqrt23 * rt * SinTt * WpmP) / SqmSqq;
        t.element({9, 15, 1, 1, 0, 0, 1}) = (-2. * sqrt23 * rt * SinTt * w2m1 * WpmP) / SqmSqq;
        t.element({6, 2, 1, 1, 1, 1, 0}) = -(((1 + rC) * SinTt * Sqw2m1 * WpmM) / (sqrt6 * SqmSqq));
        t.element({6, 3, 1, 1, 1, 1, 0}) = -((rC * SinTt * w2m1 * WpmP) / (sqrt6 * SqmSqq));
        t.element({6, 4, 1, 1, 1, 1, 0}) = -((SinTt * w2m1 * WpmP) / (sqrt6 * SqmSqq));
        t.element({6, 5, 1, 1, 1, 1, 0}) = (SinTt * (rC - w) * WpmP) / (sqrt6 * SqmSqq);
        t.element({6, 6, 1, 1, 1, 1, 0}) = -(((-1 + rC) * SinTt * Sqw2m1 * WpmP) / (sqrt6 * SqmSqq));
        t.element({6, 7, 1, 1, 1, 1, 0}) = -((rC * SinTt * w2m1 * WpmM) / (sqrt6 * SqmSqq));
        t.element({6, 8, 1, 1, 1, 1, 0}) = -((SinTt * w2m1 * WpmM) / (sqrt6 * SqmSqq));
        t.element({6, 9, 1, 1, 1, 1, 0}) = (SinTt * (rC - w) * WpmM) / (sqrt6 * SqmSqq);
        t.element({8, 2, 1, 1, 1, 1, 0}) = -(((1 + rC) * SinTt * Sqw2m1 * WpmM) / (sqrt6 * SqmSqq));
        t.element({8, 3, 1, 1, 1, 1, 0}) = -((rC * SinTt * w2m1 * WpmP) / (sqrt6 * SqmSqq));
        t.element({8, 4, 1, 1, 1, 1, 0}) = -((SinTt * w2m1 * WpmP) / (sqrt6 * SqmSqq));
        t.element({8, 5, 1, 1, 1, 1, 0}) = (SinTt * (rC - w) * WpmP) / (sqrt6 * SqmSqq);
        t.element({8, 6, 1, 1, 1, 1, 0}) = ((-1 + rC) * SinTt * Sqw2m1 * WpmP) / (sqrt6 * SqmSqq);
        t.element({8, 7, 1, 1, 1, 1, 0}) = (rC * SinTt * w2m1 * WpmM) / (sqrt6 * SqmSqq);
        t.element({8, 8, 1, 1, 1, 1, 0}) = (SinTt * w2m1 * WpmM) / (sqrt6 * SqmSqq);
        t.element({8, 9, 1, 1, 1, 1, 0}) = -((SinTt * (rC - w) * WpmM) / (sqrt6 * SqmSqq));
        t.element({10, 10, 1, 1, 1, 1, 0}) = (4. * sqrt23 * rt * SinTt * Sqw2m1 * sqwp) / SqmSqq;
        t.element({10, 11, 1, 1, 1, 1, 0}) = (-2. * sqrt23 * rt * SinTt * Sqw2m1 * WpmM) / SqmSqq;
        t.element({10, 12, 1, 1, 1, 1, 0}) = (2. * sqrt23 * rt * SinTt * Sqw2m1 * WpmM) / SqmSqq;
        t.element({10, 13, 1, 1, 1, 1, 0}) = (4. * sqrt23 * rt * SinTt * sqwm) / SqmSqq;
        t.element({10, 14, 1, 1, 1, 1, 0}) = (2. * sqrt23 * rt * SinTt * WpmP) / SqmSqq;
        t.element({10, 15, 1, 1, 1, 1, 0}) = (-2. * sqrt23 * rt * SinTt * w2m1 * WpmP) / SqmSqq;
        t.element({2, 0, 1, 1, 1, 1, 1}) = -((Sqw2m1 * WpmP) / sqrt6);
        t.element({2, 1, 1, 1, 1, 1, 1}) = (Sqw2m1 * WpmM) / sqrt6;
        t.element({4, 0, 1, 1, 1, 1, 1}) = -((Sqw2m1 * WpmP) / sqrt6);
        t.element({4, 1, 1, 1, 1, 1, 1}) = -((Sqw2m1 * WpmM) / sqrt6);
        t.element({6, 2, 1, 1, 1, 1, 1}) = -((Rm * rt * Sqw2m1) / (sqrt6 * mSqq));
        t.element({6, 3, 1, 1, 1, 1, 1}) = -((Ot * rt * Sqw2m1 * WpmP) / (sqrt6 * mSqq));
        t.element({6, 4, 1, 1, 1, 1, 1}) = -((Op * rt * Sqw2m1 * WpmP) / (sqrt6 * mSqq));
        t.element({6, 5, 1, 1, 1, 1, 1}) = (Oz * rt * WpmP) / (sqrt6 * mSqq);
        t.element({6, 6, 1, 1, 1, 1, 1}) = -((Rp * rt * Sqw2m1) / (sqrt6 * mSqq));
        t.element({6, 7, 1, 1, 1, 1, 1}) = -((Ot * rt * Sqw2m1 * WpmM) / (sqrt6 * mSqq));
        t.element({6, 8, 1, 1, 1, 1, 1}) = -((Op * rt * Sqw2m1 * WpmM) / (sqrt6 * mSqq));
        t.element({6, 9, 1, 1, 1, 1, 1}) = (Oz * rt * WpmM) / (sqrt6 * mSqq);
        t.element({8, 2, 1, 1, 1, 1, 1}) = -((Rm * rt * Sqw2m1) / (sqrt6 * mSqq));
        t.element({8, 3, 1, 1, 1, 1, 1}) = -((Ot * rt * Sqw2m1 * WpmP) / (sqrt6 * mSqq));
        t.element({8, 4, 1, 1, 1, 1, 1}) = -((Op * rt * Sqw2m1 * WpmP) / (sqrt6 * mSqq));
        t.element({8, 5, 1, 1, 1, 1, 1}) = (Oz * rt * WpmP) / (sqrt6 * mSqq);
        t.element({8, 6, 1, 1, 1, 1, 1}) = (Rp * rt * Sqw2m1) / (sqrt6 * mSqq);
        t.element({8, 7, 1, 1, 1, 1, 1}) = (Ot * rt * Sqw2m1 * WpmM) / (sqrt6 * mSqq);
        t.element({8, 8, 1, 1, 1, 1, 1}) = (Op * rt * Sqw2m1 * WpmM) / (sqrt6 * mSqq);
        t.element({8, 9, 1, 1, 1, 1, 1}) = -((Oz * rt * WpmM) / (sqrt6 * mSqq));
        t.element({10, 10, 1, 1, 1, 1, 1}) = 4. * sqrt23 * CosTt * Sqw2m1 * sqwp;
        t.element({10, 11, 1, 1, 1, 1, 1}) = -2. * sqrt23 * CosTt * Sqw2m1 * WpmM;
        t.element({10, 12, 1, 1, 1, 1, 1}) = 2. * sqrt23 * CosTt * Sqw2m1 * WpmM;
        t.element({10, 13, 1, 1, 1, 1, 1}) = 4. * sqrt23 * CosTt * sqwm;
        t.element({10, 14, 1, 1, 1, 1, 1}) = 2. * sqrt23 * CosTt * WpmP;
        t.element({10, 15, 1, 1, 1, 1, 1}) = -2. * sqrt23 * CosTt * w2m1 * WpmP;
        t.element({0, 2, 1, 2, 0, 0, 0}) = (rt * SinTt * Sqw2m1 * WpmM) / (sqrt6 * SqmSqq);
        t.element({0, 5, 1, 2, 0, 0, 0}) = -(rt * SinTt * WpmP) / (2. * sqrt6 * SqmSqq);
        t.element({0, 6, 1, 2, 0, 0, 0}) = -((rt * SinTt * Sqw2m1 * WpmP) / (sqrt6 * SqmSqq));
        t.element({0, 9, 1, 2, 0, 0, 0}) = (rt * SinTt * WpmM) / (2. * sqrt6 * SqmSqq);
        t.element({5, 2, 1, 2, 0, 0, 0}) = (rt * SinTt * Sqw2m1 * WpmM) / (sqrt6 * SqmSqq);
        t.element({5, 5, 1, 2, 0, 0, 0}) = -(rt * SinTt * WpmP) / (2. * sqrt6 * SqmSqq);
        t.element({5, 6, 1, 2, 0, 0, 0}) = (rt * SinTt * Sqw2m1 * WpmP) / (sqrt6 * SqmSqq);
        t.element({5, 9, 1, 2, 0, 0, 0}) = -(rt * SinTt * WpmM) / (2. * sqrt6 * SqmSqq);
        t.element({7, 2, 1, 2, 0, 0, 0}) = (rt * SinTt * Sqw2m1 * WpmM) / (sqrt6 * SqmSqq);
        t.element({7, 5, 1, 2, 0, 0, 0}) = -(rt * SinTt * WpmP) / (2. * sqrt6 * SqmSqq);
        t.element({7, 6, 1, 2, 0, 0, 0}) = -((rt * SinTt * Sqw2m1 * WpmP) / (sqrt6 * SqmSqq));
        t.element({7, 9, 1, 2, 0, 0, 0}) = (rt * SinTt * WpmM) / (2. * sqrt6 * SqmSqq);
        t.element({9, 10, 1, 2, 0, 0, 0}) = (4. * sqrt23 * SinTt * Sqw2m1 * (rC * wm - 1) * sqwp) / SqmSqq;
        t.element({9, 11, 1, 2, 0, 0, 0}) = (-2. * sqrt23 * SinTt * Sqw2m1 * (rC * wm - 1) * WpmM) / SqmSqq;
        t.element({9, 12, 1, 2, 0, 0, 0}) = (-2. * sqrt23 * SinTt * Sqw2m1 * (rC - wp) * WpmM) / SqmSqq;
        t.element({9, 13, 1, 2, 0, 0, 0}) = (2. * sqrt23 * SinTt * sqwp * (-rC + wp)) / SqmSqq;
        t.element({9, 14, 1, 2, 0, 0, 0}) = -((sqrt23 * SinTt * (rC - wp) * WpmP) / SqmSqq);
        t.element({0, 2, 1, 2, 0, 0, 1}) = -(sqrt23 * SinTtHalfSq * Sqw2m1 * WpmM);
        t.element({0, 5, 1, 2, 0, 0, 1}) = (SinTtHalfSq * WpmP) / sqrt6;
        t.element({0, 6, 1, 2, 0, 0, 1}) = sqrt23 * SinTtHalfSq * Sqw2m1 * WpmP;
        t.element({0, 9, 1, 2, 0, 0, 1}) = -((SinTtHalfSq * WpmM) / sqrt6);
        t.element({5, 2, 1, 2, 0, 0, 1}) = -(sqrt23 * SinTtHalfSq * Sqw2m1 * WpmM);
        t.element({5, 5, 1, 2, 0, 0, 1}) = (SinTtHalfSq * WpmP) / sqrt6;
        t.element({5, 6, 1, 2, 0, 0, 1}) = -(sqrt23 * SinTtHalfSq * Sqw2m1 * WpmP);
        t.element({5, 9, 1, 2, 0, 0, 1}) = (SinTtHalfSq * WpmM) / sqrt6;
        t.element({7, 2, 1, 2, 0, 0, 1}) = -(sqrt23 * SinTtHalfSq * Sqw2m1 * WpmM);
        t.element({7, 5, 1, 2, 0, 0, 1}) = (SinTtHalfSq * WpmP) / sqrt6;
        t.element({7, 6, 1, 2, 0, 0, 1}) = sqrt23 * SinTtHalfSq * Sqw2m1 * WpmP;
        t.element({7, 9, 1, 2, 0, 0, 1}) = -((SinTtHalfSq * WpmM) / sqrt6);
        t.element({9, 10, 1, 2, 0, 0, 1}) = (-8 * sqrt23 * rt * SinTtHalfSq * Sqw2m1 * (rC * wm - 1) * sqwp) / mSqq;
        t.element({9, 11, 1, 2, 0, 0, 1}) = (4. * sqrt23 * rt * SinTtHalfSq * Sqw2m1 * (rC * wm - 1) * WpmM) / mSqq;
        t.element({9, 12, 1, 2, 0, 0, 1}) = (4. * sqrt23 * rt * SinTtHalfSq * Sqw2m1 * (rC - wp) * WpmM) / mSqq;
        t.element({9, 13, 1, 2, 0, 0, 1}) = (4. * sqrt23 * rt * SinTtHalfSq * (rC - wp) * sqwp) / mSqq;
        t.element({9, 14, 1, 2, 0, 0, 1}) = (2. * sqrt23 * rt * SinTtHalfSq * (rC - wp) * WpmP) / mSqq;
        t.element({6, 2, 1, 2, 1, 1, 0}) = sqrt23 * CosTtHalfSq * Sqw2m1 * WpmM;
        t.element({6, 5, 1, 2, 1, 1, 0}) = -((CosTtHalfSq * WpmP) / sqrt6);
        t.element({6, 6, 1, 2, 1, 1, 0}) = sqrt23 * CosTtHalfSq * Sqw2m1 * WpmP;
        t.element({6, 9, 1, 2, 1, 1, 0}) = -((CosTtHalfSq * WpmM) / sqrt6);
        t.element({8, 2, 1, 2, 1, 1, 0}) = sqrt23 * CosTtHalfSq * Sqw2m1 * WpmM;
        t.element({8, 5, 1, 2, 1, 1, 0}) = -((CosTtHalfSq * WpmP) / sqrt6);
        t.element({8, 6, 1, 2, 1, 1, 0}) = -(sqrt23 * CosTtHalfSq * Sqw2m1 * WpmP);
        t.element({8, 9, 1, 2, 1, 1, 0}) = (CosTtHalfSq * WpmM) / sqrt6;
        t.element({10, 10, 1, 2, 1, 1, 0}) = (-8 * sqrt23 * CosTtHalfSq * rt * Sqw2m1 * sqwm * (rC * wp - 1)) / mSqq;
        t.element({10, 11, 1, 2, 1, 1, 0}) = (-4. * sqrt23 * CosTtHalfSq * rt * Sqw2m1 * (rC * wp - 1) * WpmM) / mSqq;
        t.element({10, 12, 1, 2, 1, 1, 0}) = (-4. * sqrt23 * CosTtHalfSq * rt * Sqw2m1 * (rC - wm) * WpmM) / mSqq;
        t.element({10, 13, 1, 2, 1, 1, 0}) = (-4. * sqrt23 * CosTtHalfSq * rt * (rC - wm) * sqwm) / mSqq;
        t.element({10, 14, 1, 2, 1, 1, 0}) = (-2. * sqrt23 * CosTtHalfSq * rt * (rC - wm) * WpmP) / mSqq;
        t.element({6, 2, 1, 2, 1, 1, 1}) = -((rt * SinTt * Sqw2m1 * WpmM) / (sqrt6 * SqmSqq));
        t.element({6, 5, 1, 2, 1, 1, 1}) = (rt * SinTt * WpmP) / (2. * sqrt6 * SqmSqq);
        t.element({6, 6, 1, 2, 1, 1, 1}) = -((rt * SinTt * Sqw2m1 * WpmP) / (sqrt6 * SqmSqq));
        t.element({6, 9, 1, 2, 1, 1, 1}) = (rt * SinTt * WpmM) / (2. * sqrt6 * SqmSqq);
        t.element({8, 2, 1, 2, 1, 1, 1}) = -((rt * SinTt * Sqw2m1 * WpmM) / (sqrt6 * SqmSqq));
        t.element({8, 5, 1, 2, 1, 1, 1}) = (rt * SinTt * WpmP) / (2. * sqrt6 * SqmSqq);
        t.element({8, 6, 1, 2, 1, 1, 1}) = (rt * SinTt * Sqw2m1 * WpmP) / (sqrt6 * SqmSqq);
        t.element({8, 9, 1, 2, 1, 1, 1}) = -(rt * SinTt * WpmM) / (2. * sqrt6 * SqmSqq);
        t.element({10, 10, 1, 2, 1, 1, 1}) = (4. * sqrt23 * SinTt * Sqw2m1 * sqwm * (rC * wp - 1)) / SqmSqq;
        t.element({10, 11, 1, 2, 1, 1, 1}) = (2. * sqrt23 * SinTt * Sqw2m1 * (rC * wp - 1) * WpmM) / SqmSqq;
        t.element({10, 12, 1, 2, 1, 1, 1}) = (2. * sqrt23 * SinTt * Sqw2m1 * (rC - wm) * WpmM) / SqmSqq;
        t.element({10, 13, 1, 2, 1, 1, 1}) = (2. * sqrt23 * SinTt * (rC - wm) * sqwm) / SqmSqq;
        t.element({10, 14, 1, 2, 1, 1, 1}) = (sqrt23 * SinTt * (rC - wm) * WpmP) / SqmSqq;

        t *= TwoSqTwoGfSq;
    }
    // NOLINTEND(readability-function-size)

    void AmplLbLcstar32LepNu::addRefs() const {
    }

} // namespace Hammer
