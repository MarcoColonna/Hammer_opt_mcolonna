///
/// @file  AmplBOmega3PiLepNu.cc
/// @brief \f$ B \rightarrow \omega \tau\nu, \omega \rightarrow 3\pi \f$ amplitude
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/Amplitudes/AmplBOmega3PiLepNu.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;
using namespace complex_literals;

namespace Hammer {

    namespace MD = MultiDimensional;

    AmplBOmega3PiLepNu::AmplBOmega3PiLepNu() {
        // Create tensor rank and dimensions
        IndexList dims{{11, 8, 2, 2, 2}};
        string name{"AmplBOmega3PiLepNu"};
        addProcessSignature(PID::BPLUS, {PID::OMEGA, PID::NU_TAU, PID::ANTITAU}, {PID::PIPLUS, PID::PIMINUS, PID::PI0});
        addTensor(Tensor{
            name, MD::makeEmptySparse(dims, {WILSON_BUTAUNU, FF_BOMEGA, SPIN_NUTAU, SPIN_NUTAU_REF, SPIN_TAUP})});

        addProcessSignature(PID::BPLUS, {PID::OMEGA, PID::NU_MU, PID::ANTIMUON}, {PID::PIPLUS, PID::PIMINUS, PID::PI0});
        addTensor(
            Tensor{name, MD::makeEmptySparse(dims, {WILSON_BUMUNU, FF_BOMEGA, SPIN_NUMU, SPIN_NUMU_REF, SPIN_MUP})});

        addProcessSignature(PID::BPLUS, {PID::OMEGA, PID::NU_E, PID::POSITRON}, {PID::PIPLUS, PID::PIMINUS, PID::PI0});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BUMUNU, FF_BOMEGA, SPIN_NUE, SPIN_NUE_REF, SPIN_EP})});

        setSignatureIndex();
    }

    void AmplBOmega3PiLepNu::eval(const Particle& parent, const ParticleList& daughters,
                                  const ParticleList& /*references*/) {
        // Momenta
        const FourMomentum& pBmes = parent.momentum();
        const FourMomentum& pOmegames = daughters[0].momentum();
        const FourMomentum& kNuTau = daughters[1].momentum();
        const FourMomentum& pTau = daughters[2].momentum();
        const FourMomentum& pPp = daughters[3].momentum();
        const FourMomentum& pPm = daughters[4].momentum();
        // const FourMomentum& pPz = daughters[5].momentum();

        const FourMomentum& Qmes = (pBmes - pOmegames);

        // kinematic objects
        const double Mb = pBmes.mass();
        const double Mb2 = Mb * Mb;
        // const double Mb3 = Mb2*Mb;
        const double Mu = pOmegames.mass();
        const double Mu2 = Mu * Mu;
        const double Mt = pTau.mass();
        const double Sqq = Qmes * Qmes;
        const double sqSqq = sqrt(Sqq);
        const double Ew = (Mb2 - Mu2 + Sqq) / (2 * Mb);
        const double Eu = (Mb2 + Mu2 - Sqq) / (2 * Mb);
        const double Pw = sqrt(Ew * Ew - Sqq);

        const double Mpp = pPp.mass();
        const double Mpm = pPm.mass();
        const double Epp = (pPp * pOmegames) / Mu;
        const double Epm = (pPm * pOmegames) / Mu;
        const double Ppp = sqrt(Epp * Epp - Mpp * Mpp);
        const double Ppm = sqrt(Epm * Epm - Mpm * Mpm);

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

        const double CosTep = (Epp * Eu - Mu * (pBmes * pPp) / Mb) / (Ppp * Pw);
        const double SinTep = sqrt(1. - CosTep * CosTep);

        const double CosTem = (Epm * Eu - Mu * (pBmes * pPm) / Mb) / (Ppm * Pw);
        const double SinTem = sqrt(1. - CosTem * CosTem);

        const double CosPepPt = (Sqq * (kNuTau * pPp) / NuTauQ - (Qmes * pPp) -
                                 Ppp * (Mb2 - Mu2 - Sqq) * CosTep * CosTt / (2 * Mu) + Epp * Mb * Pw * CosTt / Mu) /
                                (sqSqq * Ppp * SinTep * SinTt);
        const double SinPepPt =
            sqSqq * epsilon(pBmes, pOmegames, pPp, kNuTau) / (Mb * Pw * Ppp * TauNuTau * SinTep * SinTt);
        const complex<double> ExpIPepPt = CosPepPt + 1i * SinPepPt;
        const complex<double> ExpIPtPep = CosPepPt - 1i * SinPepPt;

        const double CosPemPt = (Sqq * (kNuTau * pPm) / NuTauQ - (Qmes * pPm) -
                                 Ppm * (Mb2 - Mu2 - Sqq) * CosTem * CosTt / (2 * Mu) + Epm * Mb * Pw * CosTt / Mu) /
                                (sqSqq * Ppm * SinTem * SinTt);
        const double SinPemPt =
            sqSqq * epsilon(pBmes, pOmegames, pPm, kNuTau) / (Mb * Pw * Ppm * TauNuTau * SinTem * SinTt);
        const complex<double> ExpIPemPt = CosPemPt + 1i * SinPemPt;
        const complex<double> ExpIPtPem = CosPemPt - 1i * SinPemPt;

        // Collective Functions
        const double Ompm = (SinPemPt * CosPepPt - SinPepPt * CosPemPt) * SinTem * SinTep;
        const double OmCpm = SinTt * (CosPemPt * SinTem * CosTep - CosPepPt * SinTep * CosTem);
        const double OmSpm = SinTt * (SinPemPt * SinTem * CosTep - SinPepPt * SinTep * CosTem);

        const complex<double> SgPlp = SinTep * CosTem * (CosTtHalfSq * ExpIPtPep + SinTtHalfSq * ExpIPepPt);
        const complex<double> SgMnp = SinTep * CosTem * (CosTtHalfSq * ExpIPtPep - SinTtHalfSq * ExpIPepPt);
        const complex<double> SgPlpc = SinTep * CosTem * (CosTtHalfSq * ExpIPepPt + SinTtHalfSq * ExpIPtPep);
        const complex<double> SgMnpc = SinTep * CosTem * (CosTtHalfSq * ExpIPepPt - SinTtHalfSq * ExpIPtPep);
        const complex<double> SgPlm = SinTem * CosTep * (CosTtHalfSq * ExpIPtPem + SinTtHalfSq * ExpIPemPt);
        const complex<double> SgMnm = SinTem * CosTep * (CosTtHalfSq * ExpIPtPem - SinTtHalfSq * ExpIPemPt);
        const complex<double> SgPlmc = SinTem * CosTep * (CosTtHalfSq * ExpIPemPt + SinTtHalfSq * ExpIPtPem);
        const complex<double> SgMnmc = SinTem * CosTep * (CosTtHalfSq * ExpIPemPt - SinTtHalfSq * ExpIPtPem);

        const double prefactor = sqrt2 * GFermi * Ppp * Ppm * sqrt(Sqq - Mt * Mt);
        const double gOmegaonSqrt2MGam = sqrt(12. * pi / Mu2);

        // initialize tensor elements to zero
        Tensor& t = getTensor();
        t.clearData();
        
        //Fix sign convention in T2 and T23 definitions, FF indices 6 and 7.
        const double signfix = -1.;
        // set non-zero tensor elements
        t.element({0, 1, 0, 0, 0}) = -((sqrt2 * Mb * Mt * Mu * OmCpm * Pw) / ((Mb + Mu) * sqSqq));
        t.element({0, 2, 0, 0, 0}) = (1i * sqrt2 * Mb * Mt * Mu * Ompm * Pw) / Sqq;
        t.element({0, 3, 0, 0, 0}) = (1i * Mt * Mu * (Mb + Mu) * OmSpm) / (sqrt2 * sqSqq);
        t.element({0, 4, 0, 0, 0}) = (-4i * sqrt2 * CosTt * Mb * Mt * Mu2 * Ompm) / Sqq;
        t.element({1, 0, 0, 0, 0}) = (-1i * Mb * Ompm * Pw) / sqrt2;
        t.element({3, 0, 0, 0, 0}) = (1i * Mb * Ompm * Pw) / sqrt2;
        t.element({5, 1, 0, 0, 0}) = -((sqrt2 * Mb * Mt * Mu * OmCpm * Pw) / ((Mb + Mu) * sqSqq));
        t.element({5, 2, 0, 0, 0}) = (-1i * sqrt2 * Mb * Mt * Mu * Ompm * Pw) / Sqq;
        t.element({5, 3, 0, 0, 0}) = (-1i * Mt * Mu * (Mb + Mu) * OmSpm) / (sqrt2 * sqSqq);
        t.element({5, 4, 0, 0, 0}) = (4i * sqrt2 * CosTt * Mb * Mt * Mu2 * Ompm) / Sqq;
        t.element({7, 1, 0, 0, 0}) = -((sqrt2 * Mb * Mt * Mu * OmCpm * Pw) / ((Mb + Mu) * sqSqq));
        t.element({7, 2, 0, 0, 0}) = (1i * sqrt2 * Mb * Mt * Mu * Ompm * Pw) / Sqq;
        t.element({7, 3, 0, 0, 0}) = (1i * Mt * Mu * (Mb + Mu) * OmSpm) / (sqrt2 * sqSqq);
        t.element({7, 4, 0, 0, 0}) = (-4i * sqrt2 * CosTt * Mb * Mt * Mu2 * Ompm) / Sqq;
        t.element({9, 5, 0, 0, 0}) = (4. * sqrt2 * Mb * Mu * OmCpm * Pw) / sqSqq;
        t.element({9, 6, 0, 0, 0}) = signfix * (2i * sqrt2 * Mu * (-Mb2 + Mu2) * OmSpm) / sqSqq;
        t.element({9, 7, 0, 0, 0}) = signfix * (8i * sqrt2 * CosTt * Mb * Mu2 * Ompm) / (Mb + Mu);
        t.element({0, 1, 0, 0, 1}) = -((sqrt2 * Mb * Mu * Pw * (SgMnm - SgMnp)) / (Mb + Mu));
        t.element({0, 3, 0, 0, 1}) = -((Mu * (Mb + Mu) * (SgPlm - SgPlp)) / sqrt2);
        t.element({0, 4, 0, 0, 1}) = (4i * sqrt2 * Mb * Mu2 * Ompm * SinTt) / sqSqq;
        t.element({5, 1, 0, 0, 1}) = -((sqrt2 * Mb * Mu * Pw * (SgMnm - SgMnp)) / (Mb + Mu));
        t.element({5, 3, 0, 0, 1}) = (Mu * (Mb + Mu) * (SgPlm - SgPlp)) / sqrt2;
        t.element({5, 4, 0, 0, 1}) = (-4i * sqrt2 * Mb * Mu2 * Ompm * SinTt) / sqSqq;
        t.element({7, 1, 0, 0, 1}) = -((sqrt2 * Mb * Mu * Pw * (SgMnm - SgMnp)) / (Mb + Mu));
        t.element({7, 3, 0, 0, 1}) = -((Mu * (Mb + Mu) * (SgPlm - SgPlp)) / sqrt2);
        t.element({7, 4, 0, 0, 1}) = (4i * sqrt2 * Mb * Mu2 * Ompm * SinTt) / sqSqq;
        t.element({9, 5, 0, 0, 1}) = (4. * sqrt2 * Mb * Mt * Mu * Pw * (SgMnm - SgMnp)) / Sqq;
        t.element({9, 6, 0, 0, 1}) = signfix * (2. * sqrt2 * Mt * Mu * (Mb2 - Mu2) * (SgPlm - SgPlp)) / Sqq;
        t.element({9, 7, 0, 0, 1}) = signfix * (-8i * sqrt2 * Mb * Mt * Mu2 * Ompm * SinTt) / ((Mb + Mu) * sqSqq);
        t.element({6, 1, 1, 1, 0}) = -((sqrt2 * Mb * Mu * Pw * (SgMnmc - SgMnpc)) / (Mb + Mu));
        t.element({6, 3, 1, 1, 0}) = -((Mu * (Mb + Mu) * (SgPlmc - SgPlpc)) / sqrt2);
        t.element({6, 4, 1, 1, 0}) = (-4i * sqrt2 * Mb * Mu2 * Ompm * SinTt) / sqSqq;
        t.element({8, 1, 1, 1, 0}) = -((sqrt2 * Mb * Mu * Pw * (SgMnmc - SgMnpc)) / (Mb + Mu));
        t.element({8, 3, 1, 1, 0}) = (Mu * (Mb + Mu) * (SgPlmc - SgPlpc)) / sqrt2;
        t.element({8, 4, 1, 1, 0}) = (4i * sqrt2 * Mb * Mu2 * Ompm * SinTt) / sqSqq;
        t.element({10, 5, 1, 1, 0}) = (4. * sqrt2 * Mb * Mt * Mu * Pw * (SgMnmc - SgMnpc)) / Sqq;
        t.element({10, 6, 1, 1, 0}) = signfix * (2. * sqrt2 * Mt * Mu * (Mb2 - Mu2) * (SgPlmc - SgPlpc)) / Sqq;
        t.element({10, 7, 1, 1, 0}) = signfix * (8i * sqrt2 * Mb * Mt * Mu2 * Ompm * SinTt) / ((Mb + Mu) * sqSqq);
        t.element({2, 0, 1, 1, 1}) = (1i * Mb * Ompm * Pw) / sqrt2;
        t.element({4, 0, 1, 1, 1}) = (-1i * Mb * Ompm * Pw) / sqrt2;
        t.element({6, 1, 1, 1, 1}) = (sqrt2 * Mb * Mt * Mu * OmCpm * Pw) / ((Mb + Mu) * sqSqq);
        t.element({6, 2, 1, 1, 1}) = (1i * sqrt2 * Mb * Mt * Mu * Ompm * Pw) / Sqq;
        t.element({6, 3, 1, 1, 1}) = (1i * Mt * Mu * (Mb + Mu) * OmSpm) / (sqrt2 * sqSqq);
        t.element({6, 4, 1, 1, 1}) = (-4i * sqrt2 * CosTt * Mb * Mt * Mu2 * Ompm) / Sqq;
        t.element({8, 1, 1, 1, 1}) = (sqrt2 * Mb * Mt * Mu * OmCpm * Pw) / ((Mb + Mu) * sqSqq);
        t.element({8, 2, 1, 1, 1}) = (-1i * sqrt2 * Mb * Mt * Mu * Ompm * Pw) / Sqq;
        t.element({8, 3, 1, 1, 1}) = (-1i * Mt * Mu * (Mb + Mu) * OmSpm) / (sqrt2 * sqSqq);
        t.element({8, 4, 1, 1, 1}) = (4i * sqrt2 * CosTt * Mb * Mt * Mu2 * Ompm) / Sqq;
        t.element({10, 5, 1, 1, 1}) = (-4. * sqrt2 * Mb * Mu * OmCpm * Pw) / sqSqq;
        t.element({10, 6, 1, 1, 1}) = signfix * (2i * sqrt2 * Mu * (-Mb2 + Mu2) * OmSpm) / sqSqq;
        t.element({10, 7, 1, 1, 1}) = signfix * (8i * sqrt2 * CosTt * Mb * Mu2 * Ompm) / (Mb + Mu);
        
        t *= prefactor * gOmegaonSqrt2MGam;
    }


} // namespace Hammer
