///
/// @file  AmplBRhoPiPiLepNu.cc
/// @brief \f$ B \rightarrow \rho \tau\nu, \rho \rightarrow \pi \pi \f$ amplitude
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/Amplitudes/AmplBRhoPiPiLepNu.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;
using namespace complex_literals;

namespace Hammer {

    namespace MD = MultiDimensional;

    AmplBRhoPiPiLepNu::AmplBRhoPiPiLepNu() {
        // Create tensor rank and dimensions
        IndexList dims{{11, 8, 2, 2, 2}};
        string name{"AmplBRhoPiPiLepNu"};
        addProcessSignature(PID::BPLUS, {PID::RHO0, PID::NU_TAU, PID::ANTITAU}, {PID::PIPLUS, PID::PIMINUS});
        addTensor(
            Tensor{name, MD::makeEmptySparse(dims, {WILSON_BUTAUNU, FF_BRHO, SPIN_NUTAU, SPIN_NUTAU_REF, SPIN_TAUP})});

        addProcessSignature(PID::BZERO, {PID::RHOMINUS, PID::NU_TAU, PID::ANTITAU}, {PID::PIMINUS, PID::PI0});
        addTensor(
            Tensor{name, MD::makeEmptySparse(dims, {WILSON_BUTAUNU, FF_BRHO, SPIN_NUTAU, SPIN_NUTAU_REF, SPIN_TAUP})});

        addProcessSignature(PID::BPLUS, {PID::RHO0, PID::NU_MU, PID::ANTIMUON}, {PID::PIPLUS, PID::PIMINUS});
        addTensor(
            Tensor{name, MD::makeEmptySparse(dims, {WILSON_BUMUNU, FF_BRHO, SPIN_NUMU, SPIN_NUMU_REF, SPIN_MUP})});

        addProcessSignature(PID::BZERO, {PID::RHOMINUS, PID::NU_MU, PID::ANTIMUON}, {PID::PIMINUS, PID::PI0});
        addTensor(
            Tensor{name, MD::makeEmptySparse(dims, {WILSON_BUMUNU, FF_BRHO, SPIN_NUMU, SPIN_NUMU_REF, SPIN_MUP})});

        addProcessSignature(PID::BPLUS, {PID::RHO0, PID::NU_E, PID::POSITRON}, {PID::PIPLUS, PID::PIMINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BUMUNU, FF_BRHO, SPIN_NUE, SPIN_NUE_REF, SPIN_EP})});

        addProcessSignature(PID::BZERO, {PID::RHOMINUS, PID::NU_E, PID::POSITRON}, {PID::PIMINUS, PID::PI0});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BUENU, FF_BRHO, SPIN_NUE, SPIN_NUE_REF, SPIN_EP})});

        setSignatureIndex();
    }

    void AmplBRhoPiPiLepNu::eval(const Particle& parent, const ParticleList& daughters,
                                 const ParticleList& /*references*/) {
        // Momenta
        const FourMomentum& pBmes = parent.momentum();
        const FourMomentum& pRhomes = daughters[0].momentum();
        const FourMomentum& kNuTau = daughters[1].momentum();
        const FourMomentum& pTau = daughters[2].momentum();
        const FourMomentum& pPp = daughters[3].momentum(); // Either pi+ or pi+ or pi-
        const FourMomentum& pPm = daughters[4].momentum(); // Either pi- or pi0 or pi0, respectively

        const FourMomentum& Qmes = (pBmes - pRhomes);

        // combinatoric for neutral unflavored mesons
        const double Cv = (daughters[0].pdgId() == PID::RHO0) ? 1. / sqrt2 : 1.;

        // kinematic objects
        const double Mb = pBmes.mass();
        const double Mb2 = Mb * Mb;
        // const double Mb3 = Mb2*Mb;
        const double Mu = pRhomes.mass();
        const double Mu2 = Mu * Mu;
        const double Mt = pTau.mass();
        const double Sqq = Qmes * Qmes;
        const double sqSqq = sqrt(Sqq);
        const double Ew = (Mb2 - Mu2 + Sqq) / (2 * Mb);
        const double Eu = (Mb2 + Mu2 - Sqq) / (2 * Mb);
        const double Pw = sqrt(Ew * Ew - Sqq);

        const double Mpp = pPp.mass();
        const double Mpm = pPm.mass();
        const double Epp = (Mu2 + Mpp * Mpp - Mpm * Mpm) / (2 * Mu);
        const double Pp = sqrt(Epp * Epp - Mpp * Mpp);

        const double BNuTau = pBmes * kNuTau;
        const double NuTauQ = Qmes * kNuTau;
        const double BQ = pBmes * Qmes;
        // const double epsBPpPmNuTau = epsilon(pBmes, pPm, pPp, kNuTau);
        const double TauNuTau = (pTau * kNuTau);

        // Helicity Angles
        const double CosTt = -((Ew * (Sqq * BNuTau - NuTauQ * BQ)) / (Pw * NuTauQ * BQ));
        const double SinTt = sqrt(1. - CosTt * CosTt);
        const double CosTtHalfSq = (1. + CosTt) / 2.;
        const double SinTtHalfSq = (1. - CosTt) / 2.;

        const double CosTe = (Epp * Eu - Mu * (pBmes * pPp) / Mb) / (Pp * Pw);
        const double SinTe = sqrt(1. - CosTe * CosTe);

        const double CosPePt = (Sqq * (kNuTau * pPp) / NuTauQ - (Qmes * pPp) -
                                Pp * (Mb2 - Mu2 - Sqq) * CosTe * CosTt / (2 * Mu) + Epp * Mb * Pw * CosTt / Mu) /
                               (sqSqq * Pp * SinTe * SinTt);
        const double SinPePt = sqSqq * epsilon(pBmes, pPm, pPp, kNuTau) / (Mb * Pw * Pp * TauNuTau * SinTe * SinTt);
        const complex<double> ExpIPePt = CosPePt + 1i * SinPePt;
        const complex<double> ExpIPtPe = CosPePt - 1i * SinPePt;

        // Collective Functions
        const double OmC = CosPePt * SinTe * SinTt;
        const double OmS = SinPePt * SinTe * SinTt;

        const complex<double> SgPl = SinTe * (CosTtHalfSq * ExpIPtPe + SinTtHalfSq * ExpIPePt);
        const complex<double> SgMn = SinTe * (CosTtHalfSq * ExpIPtPe - SinTtHalfSq * ExpIPePt);
        const complex<double> SgPlc = SinTe * (CosTtHalfSq * ExpIPePt + SinTtHalfSq * ExpIPtPe);
        const complex<double> SgMnc = SinTe * (CosTtHalfSq * ExpIPePt - SinTtHalfSq * ExpIPtPe);

        const double prefactor = 2. * Cv * GFermi * Pp * sqrt(Sqq - Mt * Mt);
        const double gRhoonSqrt2MGam = sqrt(12. * pi / Mu2);

        // initialize tensor elements to zero
        Tensor& t = getTensor();
        t.clearData();

        //Fix sign convention in T2 and T23 definitions, FF indices 6 and 7.
        const double signfix = -1.;
        // set non-zero tensor elements
        t.element({0, 1, 0, 0, 0}) = (-2i * sqrt2 * Mb * Mt * OmS * Pw) / ((Mb + Mu) * sqSqq);
        t.element({0, 2, 0, 0, 0}) = (-2. * sqrt2 * CosTe * Mb * Mt * Pw) / Sqq;
        t.element({0, 3, 0, 0, 0}) = (sqrt2 * Mt * (Mb + Mu) * OmC) / sqSqq;
        t.element({0, 4, 0, 0, 0}) = (8. * sqrt2 * CosTe * CosTt * Mb * Mt * Mu) / Sqq;
        t.element({1, 0, 0, 0, 0}) = (sqrt2 * CosTe * Mb * Pw) / Mu;
        t.element({3, 0, 0, 0, 0}) = -((sqrt2 * CosTe * Mb * Pw) / Mu);
        t.element({5, 1, 0, 0, 0}) = (-2i * sqrt2 * Mb * Mt * OmS * Pw) / ((Mb + Mu) * sqSqq);
        t.element({5, 2, 0, 0, 0}) = (2. * sqrt2 * CosTe * Mb * Mt * Pw) / Sqq;
        t.element({5, 3, 0, 0, 0}) = -((sqrt2 * Mt * (Mb + Mu) * OmC) / sqSqq);
        t.element({5, 4, 0, 0, 0}) = (-8. * sqrt2 * CosTe * CosTt * Mb * Mt * Mu) / Sqq;
        t.element({7, 1, 0, 0, 0}) = (-2i * sqrt2 * Mb * Mt * OmS * Pw) / ((Mb + Mu) * sqSqq);
        t.element({7, 2, 0, 0, 0}) = (-2. * sqrt2 * CosTe * Mb * Mt * Pw) / Sqq;
        t.element({7, 3, 0, 0, 0}) = (sqrt2 * Mt * (Mb + Mu) * OmC) / sqSqq;
        t.element({7, 4, 0, 0, 0}) = (8. * sqrt2 * CosTe * CosTt * Mb * Mt * Mu) / Sqq;
        t.element({9, 5, 0, 0, 0}) = (8i * sqrt2 * Mb * OmS * Pw) / sqSqq;
        t.element({9, 6, 0, 0, 0}) = signfix * (4. * sqrt2 * (-Mb2 + Mu2) * OmC) / sqSqq;
        t.element({9, 7, 0, 0, 0}) = signfix * (-16. * sqrt2 * CosTe * CosTt * Mb * Mu) / (Mb + Mu);
        t.element({0, 1, 0, 0, 1}) = (2. * sqrt2 * Mb * Pw * SgPl) / (Mb + Mu);
        t.element({0, 3, 0, 0, 1}) = sqrt2 * (Mb + Mu) * SgMn;
        t.element({0, 4, 0, 0, 1}) = (-8. * sqrt2 * CosTe * Mb * Mu * SinTt) / sqSqq;
        t.element({5, 1, 0, 0, 1}) = (2. * sqrt2 * Mb * Pw * SgPl) / (Mb + Mu);
        t.element({5, 3, 0, 0, 1}) = -(sqrt2 * (Mb + Mu) * SgMn);
        t.element({5, 4, 0, 0, 1}) = (8. * sqrt2 * CosTe * Mb * Mu * SinTt) / sqSqq;
        t.element({7, 1, 0, 0, 1}) = (2. * sqrt2 * Mb * Pw * SgPl) / (Mb + Mu);
        t.element({7, 3, 0, 0, 1}) = sqrt2 * (Mb + Mu) * SgMn;
        t.element({7, 4, 0, 0, 1}) = (-8. * sqrt2 * CosTe * Mb * Mu * SinTt) / sqSqq;
        t.element({9, 5, 0, 0, 1}) = (-8. * sqrt2 * Mb * Mt * Pw * SgPl) / Sqq;
        t.element({9, 6, 0, 0, 1}) = signfix * (4. * sqrt2 * Mt * (-Mb2 + Mu2) * SgMn) / Sqq;
        t.element({9, 7, 0, 0, 1}) = signfix * (16. * sqrt2 * CosTe * Mb * Mt * Mu * SinTt) / ((Mb + Mu) * sqSqq);
        t.element({6, 1, 1, 1, 0}) = (-2. * sqrt2 * Mb * Pw * SgPlc) / (Mb + Mu);
        t.element({6, 3, 1, 1, 0}) = -(sqrt2 * (Mb + Mu) * SgMnc);
        t.element({6, 4, 1, 1, 0}) = (8. * sqrt2 * CosTe * Mb * Mu * SinTt) / sqSqq;
        t.element({8, 1, 1, 1, 0}) = (-2. * sqrt2 * Mb * Pw * SgPlc) / (Mb + Mu);
        t.element({8, 3, 1, 1, 0}) = sqrt2 * (Mb + Mu) * SgMnc;
        t.element({8, 4, 1, 1, 0}) = (-8. * sqrt2 * CosTe * Mb * Mu * SinTt) / sqSqq;
        t.element({10, 5, 1, 1, 0}) = (8. * sqrt2 * Mb * Mt * Pw * SgPlc) / Sqq;
        t.element({10, 6, 1, 1, 0}) = signfix * (4. * sqrt2 * Mt * (Mb2 - Mu2) * SgMnc) / Sqq;
        t.element({10, 7, 1, 1, 0}) = signfix * (-16. * sqrt2 * CosTe * Mb * Mt * Mu * SinTt) / ((Mb + Mu) * sqSqq);
        t.element({2, 0, 1, 1, 1}) = -((sqrt2 * CosTe * Mb * Pw) / Mu);
        t.element({4, 0, 1, 1, 1}) = (sqrt2 * CosTe * Mb * Pw) / Mu;
        t.element({6, 1, 1, 1, 1}) = (2i * sqrt2 * Mb * Mt * OmS * Pw) / ((Mb + Mu) * sqSqq);
        t.element({6, 2, 1, 1, 1}) = (-2. * sqrt2 * CosTe * Mb * Mt * Pw) / Sqq;
        t.element({6, 3, 1, 1, 1}) = (sqrt2 * Mt * (Mb + Mu) * OmC) / sqSqq;
        t.element({6, 4, 1, 1, 1}) = (8. * sqrt2 * CosTe * CosTt * Mb * Mt * Mu) / Sqq;
        t.element({8, 1, 1, 1, 1}) = (2i * sqrt2 * Mb * Mt * OmS * Pw) / ((Mb + Mu) * sqSqq);
        t.element({8, 2, 1, 1, 1}) = (2. * sqrt2 * CosTe * Mb * Mt * Pw) / Sqq;
        t.element({8, 3, 1, 1, 1}) = -((sqrt2 * Mt * (Mb + Mu) * OmC) / sqSqq);
        t.element({8, 4, 1, 1, 1}) = (-8. * sqrt2 * CosTe * CosTt * Mb * Mt * Mu) / Sqq;
        t.element({10, 5, 1, 1, 1}) = (-8i * sqrt2 * Mb * OmS * Pw) / sqSqq;
        t.element({10, 6, 1, 1, 1}) = signfix * (4. * sqrt2 * (-Mb2 + Mu2) * OmC) / sqSqq;
        t.element({10, 7, 1, 1, 1}) = signfix * (-16. * sqrt2 * CosTe * CosTt * Mb * Mu) / (Mb + Mu);

        t *= prefactor * gRhoonSqrt2MGam;
    }


} // namespace Hammer
