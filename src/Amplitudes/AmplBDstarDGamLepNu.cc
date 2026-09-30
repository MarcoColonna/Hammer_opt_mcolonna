///
/// @file  AmplBDstarDGamLepNu.cc
/// @brief \f$ B \rightarrow D^* \tau\nu, D^* \rightarrow D \gamma \f$ amplitude
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/Amplitudes/AmplBDstarDGamLepNu.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;
using namespace complex_literals;

namespace Hammer {

    namespace MD = MultiDimensional;

    AmplBDstarDGamLepNu::AmplBDstarDGamLepNu() {
        // Create tensor rank and dimensions
        IndexList dims{{11, 8, 2, 2, 2, 2}};
        string name{"AmplBDstarDGamLepNu"};
        addProcessSignature(PID::BPLUS, {-PID::DSTAR, PID::NU_TAU, PID::ANTITAU}, {-PID::D0, PID::GAMMA});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCTAUNU, FF_BDSTAR, SPIN_GAMMA, SPIN_NUTAU,
                                                          SPIN_NUTAU_REF, SPIN_TAUP})});

        addProcessSignature(PID::BZERO, {PID::DSTARMINUS, PID::NU_TAU, PID::ANTITAU}, {PID::DMINUS, PID::GAMMA});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCTAUNU, FF_BDSTAR, SPIN_GAMMA, SPIN_NUTAU,
                                                          SPIN_NUTAU_REF, SPIN_TAUP})});

        addProcessSignature(PID::BPLUS, {-PID::DSTAR, PID::NU_MU, PID::ANTIMUON}, {-PID::D0, PID::GAMMA});
        addTensor(Tensor{name, MD::makeEmptySparse(
                                   dims, {WILSON_BCMUNU, FF_BDSTAR, SPIN_GAMMA, SPIN_NUMU, SPIN_NUMU_REF, SPIN_MUP})});

        addProcessSignature(PID::BZERO, {PID::DSTARMINUS, PID::NU_MU, PID::ANTIMUON}, {PID::DMINUS, PID::GAMMA});
        addTensor(Tensor{name, MD::makeEmptySparse(
                                   dims, {WILSON_BCMUNU, FF_BDSTAR, SPIN_GAMMA, SPIN_NUMU, SPIN_NUMU_REF, SPIN_MUP})});

        addProcessSignature(PID::BPLUS, {-PID::DSTAR, PID::NU_E, PID::POSITRON}, {-PID::D0, PID::GAMMA});
        addTensor(Tensor{
            name, MD::makeEmptySparse(dims, {WILSON_BCENU, FF_BDSTAR, SPIN_GAMMA, SPIN_NUE, SPIN_NUE_REF, SPIN_EP})});

        addProcessSignature(PID::BZERO, {PID::DSTARMINUS, PID::NU_E, PID::POSITRON}, {PID::DMINUS, PID::GAMMA});
        addTensor(Tensor{
            name, MD::makeEmptySparse(dims, {WILSON_BCENU, FF_BDSTAR, SPIN_GAMMA, SPIN_NUE, SPIN_NUE_REF, SPIN_EP})});

        // bs -> cs
        addProcessSignature(PID::BS, {PID::DSSTARMINUS, PID::NU_TAU, PID::ANTITAU}, {PID::DSMINUS, PID::GAMMA});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCTAUNU, FF_BSDSSTAR, SPIN_GAMMA, SPIN_NUTAU,
                                                          SPIN_NUTAU_REF, SPIN_TAUP})});

        addProcessSignature(PID::BS, {PID::DSSTARMINUS, PID::NU_MU, PID::ANTIMUON}, {PID::DSMINUS, PID::GAMMA});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCMUNU, FF_BSDSSTAR, SPIN_GAMMA, SPIN_NUMU,
                                                          SPIN_NUMU_REF, SPIN_MUP})});

        addProcessSignature(PID::BS, {PID::DSSTARMINUS, PID::NU_E, PID::POSITRON}, {PID::DSMINUS, PID::GAMMA});
        addTensor(Tensor{
            name, MD::makeEmptySparse(dims, {WILSON_BCENU, FF_BSDSSTAR, SPIN_GAMMA, SPIN_NUE, SPIN_NUE_REF, SPIN_EP})});

        setSignatureIndex();
    }

    void AmplBDstarDGamLepNu::eval(const Particle& parent, const ParticleList& daughters,
                                   const ParticleList& /*references*/) {
        // Momenta
        const FourMomentum& pBmes = parent.momentum();
        const FourMomentum& pDstarmes = daughters[0].momentum();
        const FourMomentum& kNuTau = daughters[1].momentum();
        const FourMomentum& pTau = daughters[2].momentum();
        const FourMomentum& pDmes = daughters[3].momentum();
        const FourMomentum& kPhoton = daughters[4].momentum();

        // kinematic objects
        const double Mb = pBmes.mass();
        const double Mb2 = Mb * Mb;
        const double Mds = pDstarmes.mass();
        const double Mds2 = Mds * Mds;
        const double Mt = pTau.mass();
        const double Md = pDmes.mass();
        const double Md2 = Md * Md;
        const double Sqq = Mb2 + Mds2 - 2. * (pBmes * pDstarmes);
        const double sqSqq = sqrt(Sqq);
        const double Ew = (Mb2 - Mds2 + Sqq) / (2 * Mb);
        // const double Eds = (Mb2 + Mds2 - Sqq) / (2 * Mb);
        const double Pw = sqrt(Ew * Ew - Sqq);
        const double Pw2 = Pw * Pw;
        const double BNuTau = (pBmes * kNuTau);
        const double NuTauQ = (pBmes * kNuTau) - (pDstarmes * kNuTau);
        const double BQ = Mb2 - (pBmes * pDstarmes);
        // const double DstarQ = pDstarmes * pBmes - Mds2;
        const double PhotonD = kPhoton * pDmes;
        const double PhotonQ = (kPhoton * pBmes) - (kPhoton * pDstarmes);
        const double PhotonNuTau = kPhoton * kNuTau;
        const double TauNuTau = (pTau * kNuTau);
        const double epsPhotonBDNuTau = epsilon(kPhoton, pBmes, pDmes, kNuTau);

        // Helicity Angles
        const double CosTe = (Ew * Mb * PhotonD - Sqq * PhotonD - Mds2 * PhotonQ) / (Mb * Pw * PhotonD);
        const double SinTe = sqrt(1. - CosTe * CosTe);
        // const double CosTeHalf = pow((1. + CosTe) / 2., 0.5);
        // const double SinTeHalf = pow((1. - CosTe) / 2., 0.5);
        const double CosTeHalfSq = (1. + CosTe) / 2.;
        const double SinTeHalfSq = (1. - CosTe) / 2.;
        const double CscTe = 1. / SinTe;

        const double CosTt = -((Ew * (Sqq * BNuTau - NuTauQ * BQ)) / (sqrt(Ew * Ew - Sqq) * NuTauQ * BQ));
        const double SinTt = sqrt(1. - CosTt * CosTt);
        const double CscTt = 1. / SinTt;

        // const double CosTtHalf = pow((1. + CosTt) / 2., 0.5);
        // const double SinTtHalf = pow((1. - CosTt) / 2., 0.5);
        const double CosTtHalfSq = (1. + CosTt) / 2.;
        const double SinTtHalfSq = (1. - CosTt) / 2.;

        const double CosPePt =
            (CscTe * CscTt) / (Mds * sqSqq * NuTauQ * PhotonD) *
            (Mds2 * (-PhotonQ * NuTauQ + Sqq * PhotonNuTau) +
             PhotonD * NuTauQ * (-(Mb2 - Mds2) * CosTe * CosTt + BQ * (1. + CosTe * CosTt)) - Sqq * BNuTau * PhotonD);
        const double SinPePt = (Mds * sqSqq * CscTe * CscTt * epsPhotonBDNuTau) / (Mb * Pw * PhotonD * TauNuTau);
        const complex<double> ExpIPePt = CosPePt + 1i * SinPePt;

        // Collective Functions
        const complex<double> Apl = SinTeHalfSq * CosTtHalfSq / ExpIPePt + ExpIPePt * CosTeHalfSq * SinTtHalfSq;
        const complex<double> Aplc = SinTeHalfSq * CosTtHalfSq * ExpIPePt + CosTeHalfSq * SinTtHalfSq / ExpIPePt;
        const complex<double> Ami = SinTeHalfSq * CosTtHalfSq / ExpIPePt - ExpIPePt * CosTeHalfSq * SinTtHalfSq;
        const complex<double> Amic = SinTeHalfSq * CosTtHalfSq * ExpIPePt - CosTeHalfSq * SinTtHalfSq / ExpIPePt;
        const complex<double> Aze = SinTe * SinTt;

        const complex<double> Bpl = SinTt * (CosTeHalfSq * ExpIPePt + SinTeHalfSq / ExpIPePt);
        const complex<double> Bplc = SinTt * (CosTeHalfSq / ExpIPePt + SinTeHalfSq * ExpIPePt);
        const complex<double> Bmi = SinTt * (CosTeHalfSq * ExpIPePt - SinTeHalfSq / ExpIPePt);
        const complex<double> Bmic = SinTt * (CosTeHalfSq / ExpIPePt - SinTeHalfSq * ExpIPePt);
        const complex<double> Bze = SinTe * CosTt;
        const complex<double> Bzeb = SinTe;

        // Parity Conjugate Collective Functions
        const complex<double> APpl = CosTeHalfSq * CosTtHalfSq / ExpIPePt + ExpIPePt * SinTeHalfSq * SinTtHalfSq;
        const complex<double> APplc = CosTeHalfSq * CosTtHalfSq * ExpIPePt + SinTeHalfSq * SinTtHalfSq / ExpIPePt;
        const complex<double> APmi = CosTeHalfSq * CosTtHalfSq / ExpIPePt - ExpIPePt * SinTeHalfSq * SinTtHalfSq;
        const complex<double> APmic = CosTeHalfSq * CosTtHalfSq * ExpIPePt - SinTeHalfSq * SinTtHalfSq / ExpIPePt;
        const complex<double> APze = -Aze;

        const complex<double> BPpl = SinTt * (SinTeHalfSq * ExpIPePt + CosTeHalfSq / ExpIPePt);
        const complex<double> BPplc = SinTt * (SinTeHalfSq / ExpIPePt + CosTeHalfSq * ExpIPePt);
        const complex<double> BPmi = SinTt * (SinTeHalfSq * ExpIPePt - CosTeHalfSq / ExpIPePt);
        const complex<double> BPmic = SinTt * (SinTeHalfSq / ExpIPePt - CosTeHalfSq * ExpIPePt);
        const complex<double> BPze = -Bze;
        const complex<double> BPzeb = -Bzeb;

        const double TwoGfSqMSq = 2 * GFermi * (Mds2 - Md2) * sqrt(Sqq - Mt * Mt);
        const double eMuOnSqrt2MGam = sqrt(48. * pi * Mds2 / pow(Mds2 - Md2, 3.));

        // initialize tensor elements to zero
        Tensor& t = getTensor();

        t.clearData();

        // set non-zero tensor elements
        t.element({6, 1, 0, 1, 1, 0}) = 1i * (Aplc / 2. + (Aze * (-Mb2 + Mds2 + Sqq)) / (8. * Mds * sqSqq));
        t.element({6, 2, 0, 1, 1, 0}) = 1i * Amic * Mb * Pw;
        t.element({6, 4, 0, 1, 1, 0}) = (-0.5i * Aze * Mb2 * Pw2) / (Mds * sqSqq);
        t.element({8, 1, 0, 1, 1, 0}) = -1i * (Aplc / 2. + (Aze * (-Mb2 + Mds2 + Sqq)) / (8. * Mds * sqSqq));
        t.element({8, 2, 0, 1, 1, 0}) = 1i * Amic * Mb * Pw;
        t.element({8, 4, 0, 1, 1, 0}) = (0.5i * Aze * Mb2 * Pw2) / (Mds * sqSqq);
        t.element({10, 5, 0, 1, 1, 0}) = (2i * Aze * Mb2 * Mt * Pw2) / (Mds * sqSqq);
        t.element({10, 6, 0, 1, 1, 0}) = -1i * Mt * (2. * Aplc + (Aze * (-Mb2 + Mds2 + Sqq)) / (2. * Mds * sqSqq));
        t.element({10, 7, 0, 1, 1, 0}) = 1i * Mt *
                                         ((-2. * Aplc * (Mb2 - Mds2)) / Sqq + (4. * Amic * Mb * Pw) / Sqq +
                                          (Aze * (Mb2 + 3. * Mds2 - Sqq)) / (2. * Mds * sqSqq));
        t.element({2, 0, 0, 1, 1, 1}) = (0.25i * Bzeb * Mb * Pw) / Mds;
        t.element({4, 0, 0, 1, 1, 1}) = (-0.25i * Bzeb * Mb * Pw) / Mds;
        t.element({6, 1, 0, 1, 1, 1}) = 1i * Mt *
                                        ((Bzeb * Mb * Pw) / (4. * Mds * Sqq) + Bmic / (4. * sqSqq) +
                                         (Bze * (-Mb2 + Mds2 + Sqq)) / (8. * Mds * Sqq));
        t.element({6, 2, 0, 1, 1, 1}) = (-0.5i * Bplc * Mb * Mt * Pw) / sqSqq;
        t.element({6, 3, 0, 1, 1, 1}) = (0.25i * Bzeb * Mb * Mt * Pw) / Mds;
        t.element({6, 4, 0, 1, 1, 1}) =
            -1i * Mb * Mt * Pw * ((Bzeb * (-Mb2 + Mds2)) / (4. * Mds * Sqq) + (Bze * Mb * Pw) / (2. * Mds * Sqq));
        t.element({8, 1, 0, 1, 1, 1}) = -1i * Mt *
                                        ((Bzeb * Mb * Pw) / (4. * Mds * Sqq) + Bmic / (4. * sqSqq) +
                                         (Bze * (-Mb2 + Mds2 + Sqq)) / (8. * Mds * Sqq));
        t.element({8, 2, 0, 1, 1, 1}) = (-0.5i * Bplc * Mb * Mt * Pw) / sqSqq;
        t.element({8, 3, 0, 1, 1, 1}) = (-0.25i * Bzeb * Mb * Mt * Pw) / Mds;
        t.element({8, 4, 0, 1, 1, 1}) =
            1i * Mb * Mt * Pw * ((Bzeb * (-Mb2 + Mds2)) / (4. * Mds * Sqq) + (Bze * Mb * Pw) / (2. * Mds * Sqq));
        t.element({10, 5, 0, 1, 1, 1}) = (2i * Bze * Mb2 * Pw2) / Mds;
        t.element({10, 6, 0, 1, 1, 1}) = -1i * (Bmic * sqSqq + (Bze * (-Mb2 + Mds2 + Sqq)) / (2. * Mds));
        t.element({10, 7, 0, 1, 1, 1}) = 1i * ((Bze * (Mb2 + 3. * Mds2 - Sqq)) / (2. * Mds) +
                                               (Bmic * (-Mb2 + Mds2)) / sqSqq - (2. * Bplc * Mb * Pw) / sqSqq);
        t.element({0, 1, 1, 0, 0, 0}) =
            -1i * Mt *
            ((Bzeb * Mb * Pw) / (4. * Mds * Sqq) + Bmi / (4. * sqSqq) + (Bze * (-Mb2 + Mds2 + Sqq)) / (8. * Mds * Sqq));
        t.element({0, 2, 1, 0, 0, 0}) = (0.5i * Bpl * Mb * Mt * Pw) / sqSqq;
        t.element({0, 3, 1, 0, 0, 0}) = (-0.25i * Bzeb * Mb * Mt * Pw) / Mds;
        t.element({0, 4, 1, 0, 0, 0}) =
            1i * Mb * Mt * Pw * ((Bzeb * (-Mb2 + Mds2)) / (4. * Mds * Sqq) + (Bze * Mb * Pw) / (2. * Mds * Sqq));
        t.element({1, 0, 1, 0, 0, 0}) = (0.25i * Bzeb * Mb * Pw) / Mds;
        t.element({3, 0, 1, 0, 0, 0}) = (-0.25i * Bzeb * Mb * Pw) / Mds;
        t.element({5, 1, 1, 0, 0, 0}) =
            1i * Mt *
            ((Bzeb * Mb * Pw) / (4. * Mds * Sqq) + Bmi / (4. * sqSqq) + (Bze * (-Mb2 + Mds2 + Sqq)) / (8. * Mds * Sqq));
        t.element({5, 2, 1, 0, 0, 0}) = (0.5i * Bpl * Mb * Mt * Pw) / sqSqq;
        t.element({5, 3, 1, 0, 0, 0}) = (0.25i * Bzeb * Mb * Mt * Pw) / Mds;
        t.element({5, 4, 1, 0, 0, 0}) =
            -1i * Mb * Mt * Pw * ((Bzeb * (-Mb2 + Mds2)) / (4. * Mds * Sqq) + (Bze * Mb * Pw) / (2. * Mds * Sqq));
        t.element({7, 1, 1, 0, 0, 0}) =
            -1i * Mt *
            ((Bzeb * Mb * Pw) / (4. * Mds * Sqq) + Bmi / (4. * sqSqq) + (Bze * (-Mb2 + Mds2 + Sqq)) / (8. * Mds * Sqq));
        t.element({7, 2, 1, 0, 0, 0}) = (0.5i * Bpl * Mb * Mt * Pw) / sqSqq;
        t.element({7, 3, 1, 0, 0, 0}) = (-0.25i * Bzeb * Mb * Mt * Pw) / Mds;
        t.element({7, 4, 1, 0, 0, 0}) =
            1i * Mb * Mt * Pw * ((Bzeb * (-Mb2 + Mds2)) / (4. * Mds * Sqq) + (Bze * Mb * Pw) / (2. * Mds * Sqq));
        t.element({9, 5, 1, 0, 0, 0}) = (-2i * Bze * Mb2 * Pw2) / Mds;
        t.element({9, 6, 1, 0, 0, 0}) = 1i * (Bmi * sqSqq + (Bze * (-Mb2 + Mds2 + Sqq)) / (2. * Mds));
        t.element({9, 7, 1, 0, 0, 0}) = 1i * (-(Bze * (Mb2 + 3. * Mds2 - Sqq)) / (2. * Mds) +
                                              (Bmi * (Mb2 - Mds2)) / sqSqq + (2. * Bpl * Mb * Pw) / sqSqq);
        t.element({0, 1, 1, 0, 0, 1}) = 1i * (Apl / 2. + (Aze * (-Mb2 + Mds2 + Sqq)) / (8. * Mds * sqSqq));
        t.element({0, 2, 1, 0, 0, 1}) = 1i * Ami * Mb * Pw;
        t.element({0, 4, 1, 0, 0, 1}) = (-0.5i * Aze * Mb2 * Pw2) / (Mds * sqSqq);
        t.element({5, 1, 1, 0, 0, 1}) = -1i * (Apl / 2. + (Aze * (-Mb2 + Mds2 + Sqq)) / (8. * Mds * sqSqq));
        t.element({5, 2, 1, 0, 0, 1}) = 1i * Ami * Mb * Pw;
        t.element({5, 4, 1, 0, 0, 1}) = (0.5i * Aze * Mb2 * Pw2) / (Mds * sqSqq);
        t.element({7, 1, 1, 0, 0, 1}) = 1i * (Apl / 2. + (Aze * (-Mb2 + Mds2 + Sqq)) / (8. * Mds * sqSqq));
        t.element({7, 2, 1, 0, 0, 1}) = 1i * Ami * Mb * Pw;
        t.element({7, 4, 1, 0, 0, 1}) = (-0.5i * Aze * Mb2 * Pw2) / (Mds * sqSqq);
        t.element({9, 5, 1, 0, 0, 1}) = (2i * Aze * Mb2 * Mt * Pw2) / (Mds * sqSqq);
        t.element({9, 6, 1, 0, 0, 1}) = -1i * Mt * (2. * Apl + (Aze * (-Mb2 + Mds2 + Sqq)) / (2. * Mds * sqSqq));
        t.element({9, 7, 1, 0, 0, 1}) = 1i * Mt *
                                        ((-2. * Apl * (Mb2 - Mds2)) / Sqq + (4. * Ami * Mb * Pw) / Sqq +
                                         (Aze * (Mb2 + 3. * Mds2 - Sqq)) / (2. * Mds * sqSqq));
        t.element({6, 1, 1, 1, 1, 0}) = 1i * (APplc / 2. + (APze * (-Mb2 + Mds2 + Sqq)) / (8. * Mds * sqSqq));
        t.element({6, 2, 1, 1, 1, 0}) = 1i * APmic * Mb * Pw;
        t.element({6, 4, 1, 1, 1, 0}) = (-0.5i * APze * Mb2 * Pw2) / (Mds * sqSqq);
        t.element({8, 1, 1, 1, 1, 0}) = -1i * (APplc / 2. + (APze * (-Mb2 + Mds2 + Sqq)) / (8. * Mds * sqSqq));
        t.element({8, 2, 1, 1, 1, 0}) = 1i * APmic * Mb * Pw;
        t.element({8, 4, 1, 1, 1, 0}) = (0.5i * APze * Mb2 * Pw2) / (Mds * sqSqq);
        t.element({10, 5, 1, 1, 1, 0}) = (2i * APze * Mb2 * Mt * Pw2) / (Mds * sqSqq);
        t.element({10, 6, 1, 1, 1, 0}) = -1i * Mt * (2. * APplc + (APze * (-Mb2 + Mds2 + Sqq)) / (2. * Mds * sqSqq));
        t.element({10, 7, 1, 1, 1, 0}) = 1i * Mt *
                                         ((-2. * APplc * (Mb2 - Mds2)) / Sqq + (4. * APmic * Mb * Pw) / Sqq +
                                          (APze * (Mb2 + 3. * Mds2 - Sqq)) / (2. * Mds * sqSqq));
        t.element({2, 0, 1, 1, 1, 1}) = (0.25i * BPzeb * Mb * Pw) / Mds;
        t.element({4, 0, 1, 1, 1, 1}) = (-0.25i * BPzeb * Mb * Pw) / Mds;
        t.element({6, 1, 1, 1, 1, 1}) = 1i * Mt *
                                        ((BPzeb * Mb * Pw) / (4. * Mds * Sqq) + BPmic / (4. * sqSqq) +
                                         (BPze * (-Mb2 + Mds2 + Sqq)) / (8. * Mds * Sqq));
        t.element({6, 2, 1, 1, 1, 1}) = (-0.5i * BPplc * Mb * Mt * Pw) / sqSqq;
        t.element({6, 3, 1, 1, 1, 1}) = (0.25i * BPzeb * Mb * Mt * Pw) / Mds;
        t.element({6, 4, 1, 1, 1, 1}) =
            -1i * Mb * Mt * Pw * ((BPzeb * (-Mb2 + Mds2)) / (4. * Mds * Sqq) + (BPze * Mb * Pw) / (2. * Mds * Sqq));
        t.element({8, 1, 1, 1, 1, 1}) = -1i * Mt *
                                        ((BPzeb * Mb * Pw) / (4. * Mds * Sqq) + BPmic / (4. * sqSqq) +
                                         (BPze * (-Mb2 + Mds2 + Sqq)) / (8. * Mds * Sqq));
        t.element({8, 2, 1, 1, 1, 1}) = (-0.5i * BPplc * Mb * Mt * Pw) / sqSqq;
        t.element({8, 3, 1, 1, 1, 1}) = (-0.25i * BPzeb * Mb * Mt * Pw) / Mds;
        t.element({8, 4, 1, 1, 1, 1}) =
            1i * Mb * Mt * Pw * ((BPzeb * (-Mb2 + Mds2)) / (4. * Mds * Sqq) + (BPze * Mb * Pw) / (2. * Mds * Sqq));
        t.element({10, 5, 1, 1, 1, 1}) = (2i * BPze * Mb2 * Pw2) / Mds;
        t.element({10, 6, 1, 1, 1, 1}) = -1i * (BPmic * sqSqq + (BPze * (-Mb2 + Mds2 + Sqq)) / (2. * Mds));
        t.element({10, 7, 1, 1, 1, 1}) = 1i * ((BPze * (Mb2 + 3. * Mds2 - Sqq)) / (2. * Mds) +
                                               (BPmic * (-Mb2 + Mds2)) / sqSqq - (2. * BPplc * Mb * Pw) / sqSqq);
        t.element({0, 1, 0, 0, 0, 0}) = -1i * Mt *
                                        ((BPzeb * Mb * Pw) / (4. * Mds * Sqq) + BPmi / (4. * sqSqq) +
                                         (BPze * (-Mb2 + Mds2 + Sqq)) / (8. * Mds * Sqq));
        t.element({0, 2, 0, 0, 0, 0}) = (0.5i * BPpl * Mb * Mt * Pw) / sqSqq;
        t.element({0, 3, 0, 0, 0, 0}) = (-0.25i * BPzeb * Mb * Mt * Pw) / Mds;
        t.element({0, 4, 0, 0, 0, 0}) =
            1i * Mb * Mt * Pw * ((BPzeb * (-Mb2 + Mds2)) / (4. * Mds * Sqq) + (BPze * Mb * Pw) / (2. * Mds * Sqq));
        t.element({1, 0, 0, 0, 0, 0}) = (0.25i * BPzeb * Mb * Pw) / Mds;
        t.element({3, 0, 0, 0, 0, 0}) = (-0.25i * BPzeb * Mb * Pw) / Mds;
        t.element({5, 1, 0, 0, 0, 0}) = 1i * Mt *
                                        ((BPzeb * Mb * Pw) / (4. * Mds * Sqq) + BPmi / (4. * sqSqq) +
                                         (BPze * (-Mb2 + Mds2 + Sqq)) / (8. * Mds * Sqq));
        t.element({5, 2, 0, 0, 0, 0}) = (0.5i * BPpl * Mb * Mt * Pw) / sqSqq;
        t.element({5, 3, 0, 0, 0, 0}) = (0.25i * BPzeb * Mb * Mt * Pw) / Mds;
        t.element({5, 4, 0, 0, 0, 0}) =
            -1i * Mb * Mt * Pw * ((BPzeb * (-Mb2 + Mds2)) / (4. * Mds * Sqq) + (BPze * Mb * Pw) / (2. * Mds * Sqq));
        t.element({7, 1, 0, 0, 0, 0}) = -1i * Mt *
                                        ((BPzeb * Mb * Pw) / (4. * Mds * Sqq) + BPmi / (4. * sqSqq) +
                                         (BPze * (-Mb2 + Mds2 + Sqq)) / (8. * Mds * Sqq));
        t.element({7, 2, 0, 0, 0, 0}) = (0.5i * BPpl * Mb * Mt * Pw) / sqSqq;
        t.element({7, 3, 0, 0, 0, 0}) = (-0.25i * BPzeb * Mb * Mt * Pw) / Mds;
        t.element({7, 4, 0, 0, 0, 0}) =
            1i * Mb * Mt * Pw * ((BPzeb * (-Mb2 + Mds2)) / (4. * Mds * Sqq) + (BPze * Mb * Pw) / (2. * Mds * Sqq));
        t.element({9, 5, 0, 0, 0, 0}) = (-2i * BPze * Mb2 * Pw2) / Mds;
        t.element({9, 6, 0, 0, 0, 0}) = 1i * (BPmi * sqSqq + (BPze * (-Mb2 + Mds2 + Sqq)) / (2. * Mds));
        t.element({9, 7, 0, 0, 0, 0}) = 1i * (-(BPze * (Mb2 + 3. * Mds2 - Sqq)) / (2. * Mds) +
                                              (BPmi * (Mb2 - Mds2)) / sqSqq + (2. * BPpl * Mb * Pw) / sqSqq);
        t.element({0, 1, 0, 0, 0, 1}) = 1i * (APpl / 2. + (APze * (-Mb2 + Mds2 + Sqq)) / (8. * Mds * sqSqq));
        t.element({0, 2, 0, 0, 0, 1}) = 1i * APmi * Mb * Pw;
        t.element({0, 4, 0, 0, 0, 1}) = (-0.5i * APze * Mb2 * Pw2) / (Mds * sqSqq);
        t.element({5, 1, 0, 0, 0, 1}) = -1i * (APpl / 2. + (APze * (-Mb2 + Mds2 + Sqq)) / (8. * Mds * sqSqq));
        t.element({5, 2, 0, 0, 0, 1}) = 1i * APmi * Mb * Pw;
        t.element({5, 4, 0, 0, 0, 1}) = (0.5i * APze * Mb2 * Pw2) / (Mds * sqSqq);
        t.element({7, 1, 0, 0, 0, 1}) = 1i * (APpl / 2. + (APze * (-Mb2 + Mds2 + Sqq)) / (8. * Mds * sqSqq));
        t.element({7, 2, 0, 0, 0, 1}) = 1i * APmi * Mb * Pw;
        t.element({7, 4, 0, 0, 0, 1}) = (-0.5i * APze * Mb2 * Pw2) / (Mds * sqSqq);
        t.element({9, 5, 0, 0, 0, 1}) = (2i * APze * Mb2 * Mt * Pw2) / (Mds * sqSqq);
        t.element({9, 6, 0, 0, 0, 1}) = -1i * Mt * (2. * APpl + (APze * (-Mb2 + Mds2 + Sqq)) / (2. * Mds * sqSqq));
        t.element({9, 7, 0, 0, 0, 1}) = 1i * Mt *
                                        ((-2. * APpl * (Mb2 - Mds2)) / Sqq + (4. * APmi * Mb * Pw) / Sqq +
                                         (APze * (Mb2 + 3. * Mds2 - Sqq)) / (2. * Mds * sqSqq));

        t *= TwoGfSqMSq * eMuOnSqrt2MGam;
    }

    void AmplBDstarDGamLepNu::addRefs() const {
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
