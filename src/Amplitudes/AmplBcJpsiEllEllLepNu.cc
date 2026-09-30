///
/// @file  AmplBcJpsiEllEllLepNu.cc
/// @brief \f$ B_c \rightarrow J/\psi \tau\nu, J/\psi \rightarrow \ell \ell \f$ amplitude
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/Amplitudes/AmplBcJpsiEllEllLepNu.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;
using namespace complex_literals;

namespace Hammer {

    namespace MD = MultiDimensional;

    AmplBcJpsiEllEllLepNu::AmplBcJpsiEllEllLepNu() {
        // Create tensor rank and dimensions
        IndexList dims{{11, 8, 2, 2, 2, 2}};
        string name{"AmplBcJpsiEllEllLepNu"};
        // Use ref spin labels for J/psi daughters to avoid index collisions
        addProcessSignature(PID::BCPLUS, {PID::JPSI, PID::NU_TAU, PID::ANTITAU}, {PID::MUON, PID::ANTIMUON});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCTAUNU, FF_BCJPSI, SPIN_MUM_REF, SPIN_NUTAU,
                                                          SPIN_NUTAU_REF, SPIN_TAUP})});

        addProcessSignature(PID::BCPLUS, {PID::JPSI, PID::NU_TAU, PID::ANTITAU}, {PID::ELECTRON, PID::POSITRON});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCTAUNU, FF_BCJPSI, SPIN_EM_REF, SPIN_NUTAU,
                                                          SPIN_NUTAU_REF, SPIN_TAUP})});

        addProcessSignature(PID::BCPLUS, {PID::JPSI, PID::NU_MU, PID::ANTIMUON}, {PID::MUON, PID::ANTIMUON});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCMUNU, FF_BCJPSI, SPIN_MUM_REF, SPIN_NUMU,
                                                          SPIN_NUMU_REF, SPIN_MUP})});

        addProcessSignature(PID::BCPLUS, {PID::JPSI, PID::NU_MU, PID::ANTIMUON}, {PID::ELECTRON, PID::POSITRON});
        addTensor(Tensor{name, MD::makeEmptySparse(
                                   dims, {WILSON_BCMUNU, FF_BCJPSI, SPIN_EM_REF, SPIN_NUMU, SPIN_NUMU_REF, SPIN_MUP})});

        addProcessSignature(PID::BCPLUS, {PID::JPSI, PID::NU_E, PID::POSITRON}, {PID::MUON, PID::ANTIMUON});
        addTensor(Tensor{
            name, MD::makeEmptySparse(dims, {WILSON_BCENU, FF_BCJPSI, SPIN_MUM_REF, SPIN_NUE, SPIN_NUE_REF, SPIN_EP})});

        addProcessSignature(PID::BCPLUS, {PID::JPSI, PID::NU_E, PID::POSITRON}, {PID::ELECTRON, PID::POSITRON});
        addTensor(Tensor{
            name, MD::makeEmptySparse(dims, {WILSON_BCENU, FF_BCJPSI, SPIN_EP_REF, SPIN_NUE, SPIN_NUE_REF, SPIN_EP})});

        setSignatureIndex();
    }

    void AmplBcJpsiEllEllLepNu::eval(const Particle& parent, const ParticleList& daughters,
                                     const ParticleList& /*references*/) {
        // Momenta
        const FourMomentum& pBmes = parent.momentum();
        const FourMomentum& pJpsimes = daughters[0].momentum();
        const FourMomentum& kNuTau = daughters[1].momentum();
        const FourMomentum& pTau = daughters[2].momentum();
        const FourMomentum& pLm = daughters[3].momentum();
        const FourMomentum& pLp = daughters[4].momentum();

        const FourMomentum& Qmes = (pBmes - pJpsimes);
        const FourMomentum& pLmLp = (pLm - pLp);

        // kinematic objects
        const double Mb = pBmes.mass();
        const double Mb2 = Mb * Mb;
        const double Mb3 = Mb2 * Mb;
        const double Mc = pJpsimes.mass();
        const double Mc2 = Mc * Mc;
        const double Mt = pTau.mass();
        const double Sqq = Qmes * Qmes;
        const double sqSqq = sqrt(Sqq);
        const double Ew = (Mb2 - Mc2 + Sqq) / (2 * Mb);
        const double Pw = sqrt(Ew * Ew - Sqq);

        const double rC = Mc / Mb;
        const double rC2 = rC * rC;
        const double w = (Mb2 + Mc2 - Sqq) / (2 * Mb * Mc);
        const double sqw2m1 = sqrt(w * w - 1.);
        const double mSqq = Sqq / Mb2;
        const double sqmSqq = sqSqq / Mb;
        const double rt = Mt / Mb;

        const double BNuTau = pBmes * kNuTau;
        const double NuTauQ = Qmes * kNuTau;
        const double BQ = pBmes * Qmes;
        const double LepLmLp = sqrt(-(pLmLp * pLmLp));
        const double BLep = pBmes * pLmLp;
        const double QLep = Qmes * pLmLp;
        const double NuTauLep = kNuTau * pLmLp;
        const double epsBLmLpNuTau = epsilon(pBmes, pLm, pLp, kNuTau);
        const double TauNuTau = (pTau * kNuTau);

        // Helicity Angles
        const double CosTe = Mc * BLep / (Mb * LepLmLp * Pw);
        const double SinTe = sqrt(1. - CosTe * CosTe);
        const double CosTeHalfSq = (1. + CosTe) / 2.;
        const double SinTeHalfSq = (1. - CosTe) / 2.;

        const double CosTt = -((Ew * (Sqq * BNuTau - NuTauQ * BQ)) / (Pw * NuTauQ * BQ));
        const double SinTt = sqrt(1. - CosTt * CosTt);
        const double CosTtHalfSq = (1. + CosTt) / 2.;
        const double SinTtHalfSq = (1. - CosTt) / 2.;

        const double CosPePt =
            ((NuTauQ * QLep - Sqq * NuTauLep) / (NuTauQ * LepLmLp) - CosTe * CosTt * (Mb2 - Mc2 - BQ) / Mc) /
            (sqSqq * SinTe * SinTt);
        const double SinPePt = 2. * sqSqq * epsBLmLpNuTau / (SinTe * SinTt * Mb * Pw * TauNuTau * LepLmLp);
        const complex<double> ExpIPePt = CosPePt + 1i * SinPePt;
        const complex<double> ExpIPtPe = CosPePt - 1i * SinPePt;

        // Collective Functions
        const double OmZ = SinTe * (CosTt * (rC - w) + sqw2m1);
        const double OmZb = SinTe * (rC2 + 2. * CosTt * rC * sqw2m1 - 1.);
        const complex<double> OmPl = SinTt * (CosTeHalfSq * ExpIPtPe + ExpIPePt * SinTeHalfSq);
        const complex<double> OmMn = SinTt * (CosTeHalfSq * ExpIPtPe - ExpIPePt * SinTeHalfSq);
        const complex<double> OmPlc = SinTt * (CosTeHalfSq * ExpIPePt + ExpIPtPe * SinTeHalfSq);
        const complex<double> OmMnc = SinTt * (CosTeHalfSq * ExpIPePt - ExpIPtPe * SinTeHalfSq);

        const complex<double> SgPl = CosTeHalfSq * CosTtHalfSq * ExpIPtPe + SinTeHalfSq * SinTtHalfSq * ExpIPePt;
        const complex<double> SgMn = CosTeHalfSq * CosTtHalfSq * ExpIPtPe - SinTeHalfSq * SinTtHalfSq * ExpIPePt;
        const complex<double> SgPlc = CosTeHalfSq * CosTtHalfSq * ExpIPePt + SinTeHalfSq * SinTtHalfSq * ExpIPtPe;
        const complex<double> SgMnc = CosTeHalfSq * CosTtHalfSq * ExpIPePt - SinTeHalfSq * SinTtHalfSq * ExpIPtPe;
        const complex<double> SgPlb = CosTtHalfSq * SinTeHalfSq * ExpIPtPe + CosTeHalfSq * SinTtHalfSq * ExpIPePt;
        const complex<double> SgMnb = CosTtHalfSq * SinTeHalfSq * ExpIPtPe - CosTeHalfSq * SinTtHalfSq * ExpIPePt;
        const complex<double> SgPlbc = CosTtHalfSq * SinTeHalfSq * ExpIPePt + CosTeHalfSq * SinTtHalfSq * ExpIPtPe;
        const complex<double> SgMnbc = CosTtHalfSq * SinTeHalfSq * ExpIPePt - CosTeHalfSq * SinTtHalfSq * ExpIPtPe;

        const double prefactor = TwoSqTwoGFermi * Mb * Pw * sqrt(Sqq - Mt * Mt) / sqSqq;
        const double gPsionSqrt2MGam = sqrt(12. * pi / Mc2);

        // initialize tensor elements to zero
        Tensor& t = getTensor();
        t.clearData();

        // set non-zero tensor elements
        t.element({0, 1, 0, 0, 0, 0}) = (rt * (-OmZ - OmMnc * sqmSqq)) / (2. * sqmSqq * sqw2m1);
        t.element({0, 2, 0, 0, 0, 0}) = Mb2 * OmPlc * rC * rt;
        t.element({0, 3, 0, 0, 0, 0}) = -(Mb2 * rt * SinTe * sqmSqq) / 2.;
        t.element({0, 4, 0, 0, 0, 0}) = (Mb2 * OmZb * rt) / (2. * sqmSqq);
        t.element({1, 0, 0, 0, 0, 0}) = (Mb * SinTe * sqmSqq) / 2.;
        t.element({3, 0, 0, 0, 0, 0}) = -(Mb * SinTe * sqmSqq) / 2.;
        t.element({5, 1, 0, 0, 0, 0}) = (rt * (OmZ + OmMnc * sqmSqq)) / (2. * sqmSqq * sqw2m1);
        t.element({5, 2, 0, 0, 0, 0}) = Mb2 * OmPlc * rC * rt;
        t.element({5, 3, 0, 0, 0, 0}) = (Mb2 * rt * SinTe * sqmSqq) / 2.;
        t.element({5, 4, 0, 0, 0, 0}) = -(Mb2 * OmZb * rt) / (2. * sqmSqq);
        t.element({7, 1, 0, 0, 0, 0}) = (rt * (-OmZ - OmMnc * sqmSqq)) / (2. * sqmSqq * sqw2m1);
        t.element({7, 2, 0, 0, 0, 0}) = Mb2 * OmPlc * rC * rt;
        t.element({7, 3, 0, 0, 0, 0}) = -(Mb2 * rt * SinTe * sqmSqq) / 2.;
        t.element({7, 4, 0, 0, 0, 0}) = (Mb2 * OmZb * rt) / (2. * sqmSqq);
        t.element({9, 5, 0, 0, 0, 0}) = -4. * CosTt * Mb3 * rC * SinTe * sqmSqq * sqw2m1;
        t.element({9, 6, 0, 0, 0, 0}) = (Mb * (2. * OmMnc * mSqq + 2. * CosTt * SinTe * sqmSqq * (rC - w))) / sqw2m1;
        t.element({9, 7, 0, 0, 0, 0}) =
            (Mb * (OmMnc * (2 - 2. * rC2) + 4. * OmPlc * rC * sqw2m1 - 2. * CosTt * SinTe * sqmSqq * (rC + w))) /
            sqw2m1;
        t.element({0, 1, 0, 0, 0, 1}) = (2. * SgPlb * sqmSqq + SinTe * SinTt * (rC - w)) / (2. * sqw2m1);
        t.element({0, 2, 0, 0, 0, 1}) = 2. * Mb2 * rC * SgMnb * sqmSqq;
        t.element({0, 4, 0, 0, 0, 1}) = -(Mb2 * rC * SinTe * SinTt * sqw2m1);
        t.element({5, 1, 0, 0, 0, 1}) = (-2. * SgPlb * sqmSqq + SinTe * SinTt * (-rC + w)) / (2. * sqw2m1);
        t.element({5, 2, 0, 0, 0, 1}) = 2. * Mb2 * rC * SgMnb * sqmSqq;
        t.element({5, 4, 0, 0, 0, 1}) = Mb2 * rC * SinTe * SinTt * sqw2m1;
        t.element({7, 1, 0, 0, 0, 1}) = (2. * SgPlb * sqmSqq + SinTe * SinTt * (rC - w)) / (2. * sqw2m1);
        t.element({7, 2, 0, 0, 0, 1}) = 2. * Mb2 * rC * SgMnb * sqmSqq;
        t.element({7, 4, 0, 0, 0, 1}) = -(Mb2 * rC * SinTe * SinTt * sqw2m1);
        t.element({9, 5, 0, 0, 0, 1}) = 4. * Mb3 * rC * rt * SinTe * SinTt * sqw2m1;
        t.element({9, 6, 0, 0, 0, 1}) = (Mb * rt * (-4. * SgPlb * sqmSqq - 2. * SinTe * SinTt * (rC - w))) / sqw2m1;
        t.element({9, 7, 0, 0, 0, 1}) =
            (Mb * rt * (4. * (-1 + rC2) * SgPlb + 8 * rC * SgMnb * sqw2m1 + 2. * SinTe * SinTt * sqmSqq * (rC + w))) /
            (sqmSqq * sqw2m1);
        t.element({6, 1, 0, 1, 1, 0}) = (2. * SgPlc * sqmSqq + SinTe * SinTt * (-rC + w)) / (2. * sqw2m1);
        t.element({6, 2, 0, 1, 1, 0}) = 2. * Mb2 * rC * SgMnc * sqmSqq;
        t.element({6, 4, 0, 1, 1, 0}) = Mb2 * rC * SinTe * SinTt * sqw2m1;
        t.element({8, 1, 0, 1, 1, 0}) = (-2. * SgPlc * sqmSqq + SinTe * SinTt * (rC - w)) / (2. * sqw2m1);
        t.element({8, 2, 0, 1, 1, 0}) = 2. * Mb2 * rC * SgMnc * sqmSqq;
        t.element({8, 4, 0, 1, 1, 0}) = -(Mb2 * rC * SinTe * SinTt * sqw2m1);
        t.element({10, 5, 0, 1, 1, 0}) = -4. * Mb3 * rC * rt * SinTe * SinTt * sqw2m1;
        t.element({10, 6, 0, 1, 1, 0}) = (Mb * rt * (-4. * SgPlc * sqmSqq + 2. * SinTe * SinTt * (rC - w))) / sqw2m1;
        t.element({10, 7, 0, 1, 1, 0}) =
            (Mb * rt * (4. * (-1 + rC2) * SgPlc + 8 * rC * SgMnc * sqw2m1 - 2. * SinTe * SinTt * sqmSqq * (rC + w))) /
            (sqmSqq * sqw2m1);
        t.element({2, 0, 0, 1, 1, 1}) = -(Mb * SinTe * sqmSqq) / 2.;
        t.element({4, 0, 0, 1, 1, 1}) = (Mb * SinTe * sqmSqq) / 2.;
        t.element({6, 1, 0, 1, 1, 1}) = (rt * (-OmZ - OmMnc * sqmSqq)) / (2. * sqmSqq * sqw2m1);
        t.element({6, 2, 0, 1, 1, 1}) = -(Mb2 * OmPlc * rC * rt);
        t.element({6, 3, 0, 1, 1, 1}) = -(Mb2 * rt * SinTe * sqmSqq) / 2.;
        t.element({6, 4, 0, 1, 1, 1}) = (Mb2 * OmZb * rt) / (2. * sqmSqq);
        t.element({8, 1, 0, 1, 1, 1}) = (rt * (OmZ + OmMnc * sqmSqq)) / (2. * sqmSqq * sqw2m1);
        t.element({8, 2, 0, 1, 1, 1}) = -(Mb2 * OmPlc * rC * rt);
        t.element({8, 3, 0, 1, 1, 1}) = (Mb2 * rt * SinTe * sqmSqq) / 2.;
        t.element({8, 4, 0, 1, 1, 1}) = -(Mb2 * OmZb * rt) / (2. * sqmSqq);
        t.element({10, 5, 0, 1, 1, 1}) = -4. * CosTt * Mb3 * rC * SinTe * sqmSqq * sqw2m1;
        t.element({10, 6, 0, 1, 1, 1}) = (Mb * (2. * OmMnc * mSqq + 2. * CosTt * SinTe * sqmSqq * (rC - w))) / sqw2m1;
        t.element({10, 7, 0, 1, 1, 1}) =
            (Mb * (OmMnc * (2 - 2. * rC2) - 4. * OmPlc * rC * sqw2m1 - 2. * CosTt * SinTe * sqmSqq * (rC + w))) /
            sqw2m1;
        t.element({0, 1, 1, 0, 0, 0}) = (rt * (-OmZ - OmMn * sqmSqq)) / (2. * sqmSqq * sqw2m1);
        t.element({0, 2, 1, 0, 0, 0}) = -(Mb2 * OmPl * rC * rt);
        t.element({0, 3, 1, 0, 0, 0}) = -(Mb2 * rt * SinTe * sqmSqq) / 2.;
        t.element({0, 4, 1, 0, 0, 0}) = (Mb2 * OmZb * rt) / (2. * sqmSqq);
        t.element({1, 0, 1, 0, 0, 0}) = (Mb * SinTe * sqmSqq) / 2.;
        t.element({3, 0, 1, 0, 0, 0}) = -(Mb * SinTe * sqmSqq) / 2.;
        t.element({5, 1, 1, 0, 0, 0}) = (rt * (OmZ + OmMn * sqmSqq)) / (2. * sqmSqq * sqw2m1);
        t.element({5, 2, 1, 0, 0, 0}) = -(Mb2 * OmPl * rC * rt);
        t.element({5, 3, 1, 0, 0, 0}) = (Mb2 * rt * SinTe * sqmSqq) / 2.;
        t.element({5, 4, 1, 0, 0, 0}) = -(Mb2 * OmZb * rt) / (2. * sqmSqq);
        t.element({7, 1, 1, 0, 0, 0}) = (rt * (-OmZ - OmMn * sqmSqq)) / (2. * sqmSqq * sqw2m1);
        t.element({7, 2, 1, 0, 0, 0}) = -(Mb2 * OmPl * rC * rt);
        t.element({7, 3, 1, 0, 0, 0}) = -(Mb2 * rt * SinTe * sqmSqq) / 2.;
        t.element({7, 4, 1, 0, 0, 0}) = (Mb2 * OmZb * rt) / (2. * sqmSqq);
        t.element({9, 5, 1, 0, 0, 0}) = -4. * CosTt * Mb3 * rC * SinTe * sqmSqq * sqw2m1;
        t.element({9, 6, 1, 0, 0, 0}) = (Mb * (2. * OmMn * mSqq + 2. * CosTt * SinTe * sqmSqq * (rC - w))) / sqw2m1;
        t.element({9, 7, 1, 0, 0, 0}) =
            (Mb * (OmMn * (2 - 2. * rC2) - 4. * OmPl * rC * sqw2m1 - 2. * CosTt * SinTe * sqmSqq * (rC + w))) / sqw2m1;
        t.element({0, 1, 1, 0, 0, 1}) = (-2. * SgPl * sqmSqq + SinTe * SinTt * (rC - w)) / (2. * sqw2m1);
        t.element({0, 2, 1, 0, 0, 1}) = -2. * Mb2 * rC * SgMn * sqmSqq;
        t.element({0, 4, 1, 0, 0, 1}) = -(Mb2 * rC * SinTe * SinTt * sqw2m1);
        t.element({5, 1, 1, 0, 0, 1}) = (2. * SgPl * sqmSqq + SinTe * SinTt * (-rC + w)) / (2. * sqw2m1);
        t.element({5, 2, 1, 0, 0, 1}) = -2. * Mb2 * rC * SgMn * sqmSqq;
        t.element({5, 4, 1, 0, 0, 1}) = Mb2 * rC * SinTe * SinTt * sqw2m1;
        t.element({7, 1, 1, 0, 0, 1}) = (-2. * SgPl * sqmSqq + SinTe * SinTt * (rC - w)) / (2. * sqw2m1);
        t.element({7, 2, 1, 0, 0, 1}) = -2. * Mb2 * rC * SgMn * sqmSqq;
        t.element({7, 4, 1, 0, 0, 1}) = -(Mb2 * rC * SinTe * SinTt * sqw2m1);
        t.element({9, 5, 1, 0, 0, 1}) = 4. * Mb3 * rC * rt * SinTe * SinTt * sqw2m1;
        t.element({9, 6, 1, 0, 0, 1}) = (Mb * rt * (4. * SgPl * sqmSqq + 2. * SinTe * SinTt * (-rC + w))) / sqw2m1;
        t.element({9, 7, 1, 0, 0, 1}) =
            (Mb * rt * ((4 - 4. * rC2) * SgPl - 8 * rC * SgMn * sqw2m1 + 2. * SinTe * SinTt * sqmSqq * (rC + w))) /
            (sqmSqq * sqw2m1);
        t.element({6, 1, 1, 1, 1, 0}) = (-2. * SgPlbc * sqmSqq + SinTe * SinTt * (-rC + w)) / (2. * sqw2m1);
        t.element({6, 2, 1, 1, 1, 0}) = -2. * Mb2 * rC * SgMnbc * sqmSqq;
        t.element({6, 4, 1, 1, 1, 0}) = Mb2 * rC * SinTe * SinTt * sqw2m1;
        t.element({8, 1, 1, 1, 1, 0}) = (2. * SgPlbc * sqmSqq + SinTe * SinTt * (rC - w)) / (2. * sqw2m1);
        t.element({8, 2, 1, 1, 1, 0}) = -2. * Mb2 * rC * SgMnbc * sqmSqq;
        t.element({8, 4, 1, 1, 1, 0}) = -(Mb2 * rC * SinTe * SinTt * sqw2m1);
        t.element({10, 5, 1, 1, 1, 0}) = -4. * Mb3 * rC * rt * SinTe * SinTt * sqw2m1;
        t.element({10, 6, 1, 1, 1, 0}) = (Mb * rt * (4. * SgPlbc * sqmSqq + 2. * SinTe * SinTt * (rC - w))) / sqw2m1;
        t.element({10, 7, 1, 1, 1, 0}) =
            (Mb * rt * ((4 - 4. * rC2) * SgPlbc - 8 * rC * SgMnbc * sqw2m1 - 2. * SinTe * SinTt * sqmSqq * (rC + w))) /
            (sqmSqq * sqw2m1);
        t.element({2, 0, 1, 1, 1, 1}) = -(Mb * SinTe * sqmSqq) / 2.;
        t.element({4, 0, 1, 1, 1, 1}) = (Mb * SinTe * sqmSqq) / 2.;
        t.element({6, 1, 1, 1, 1, 1}) = (rt * (-OmZ - OmMn * sqmSqq)) / (2. * sqmSqq * sqw2m1);
        t.element({6, 2, 1, 1, 1, 1}) = Mb2 * OmPl * rC * rt;
        t.element({6, 3, 1, 1, 1, 1}) = -(Mb2 * rt * SinTe * sqmSqq) / 2.;
        t.element({6, 4, 1, 1, 1, 1}) = (Mb2 * OmZb * rt) / (2. * sqmSqq);
        t.element({8, 1, 1, 1, 1, 1}) = (rt * (OmZ + OmMn * sqmSqq)) / (2. * sqmSqq * sqw2m1);
        t.element({8, 2, 1, 1, 1, 1}) = Mb2 * OmPl * rC * rt;
        t.element({8, 3, 1, 1, 1, 1}) = (Mb2 * rt * SinTe * sqmSqq) / 2.;
        t.element({8, 4, 1, 1, 1, 1}) = -(Mb2 * OmZb * rt) / (2. * sqmSqq);
        t.element({10, 5, 1, 1, 1, 1}) = -4. * CosTt * Mb3 * rC * SinTe * sqmSqq * sqw2m1;
        t.element({10, 6, 1, 1, 1, 1}) = (Mb * (2. * OmMn * mSqq + 2. * CosTt * SinTe * sqmSqq * (rC - w))) / sqw2m1;
        t.element({10, 7, 1, 1, 1, 1}) =
            (Mb * (OmMn * (2 - 2. * rC2) + 4. * OmPl * rC * sqw2m1 - 2. * CosTt * SinTe * sqmSqq * (rC + w))) / sqw2m1;


        t *= prefactor * gPsionSqrt2MGam;
    }


} // namespace Hammer
