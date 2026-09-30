///
/// @file  AmplBDstarDPiLepNu.cc
/// @brief \f$ B \rightarrow D^* \tau\nu, D^* \rightarrow D \pi \f$ amplitude
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/Amplitudes/AmplBDstarDPiLepNu.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;
using namespace complex_literals;

namespace Hammer {

    namespace MD = MultiDimensional;

    AmplBDstarDPiLepNu::AmplBDstarDPiLepNu() {
        // Create tensor rank and dimensions
        IndexList dims{{11, 8, 2, 2, 2}};
        string name{"AmplBDstarDPiLepNu"};
        addProcessSignature(PID::BPLUS, {-PID::DSTAR, PID::NU_TAU, PID::ANTITAU}, {-PID::D0, PID::PI0});
        addTensor(Tensor{
            name, MD::makeEmptySparse(dims, {WILSON_BCTAUNU, FF_BDSTAR, SPIN_NUTAU, SPIN_NUTAU_REF, SPIN_TAUP})});

        addProcessSignature(PID::BZERO, {PID::DSTARMINUS, PID::NU_TAU, PID::ANTITAU}, {PID::DMINUS, PID::PI0});
        addTensor(Tensor{
            name, MD::makeEmptySparse(dims, {WILSON_BCTAUNU, FF_BDSTAR, SPIN_NUTAU, SPIN_NUTAU_REF, SPIN_TAUP})});

        addProcessSignature(PID::BZERO, {PID::DSTARMINUS, PID::NU_TAU, PID::ANTITAU}, {-PID::D0, PID::PIMINUS});
        addTensor(Tensor{
            name, MD::makeEmptySparse(dims, {WILSON_BCTAUNU, FF_BDSTAR, SPIN_NUTAU, SPIN_NUTAU_REF, SPIN_TAUP})});

        addProcessSignature(PID::BPLUS, {-PID::DSTAR, PID::NU_MU, PID::ANTIMUON}, {-PID::D0, PID::PI0});
        addTensor(
            Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCMUNU, FF_BDSTAR, SPIN_NUMU, SPIN_NUMU_REF, SPIN_MUP})});

        addProcessSignature(PID::BZERO, {PID::DSTARMINUS, PID::NU_MU, PID::ANTIMUON}, {PID::DMINUS, PID::PI0});
        addTensor(
            Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCMUNU, FF_BDSTAR, SPIN_NUMU, SPIN_NUMU_REF, SPIN_MUP})});

        addProcessSignature(PID::BZERO, {PID::DSTARMINUS, PID::NU_MU, PID::ANTIMUON}, {-PID::D0, PID::PIMINUS});
        addTensor(
            Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCMUNU, FF_BDSTAR, SPIN_NUMU, SPIN_NUMU_REF, SPIN_MUP})});

        addProcessSignature(PID::BPLUS, {-PID::DSTAR, PID::NU_E, PID::POSITRON}, {-PID::D0, PID::PI0});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCENU, FF_BDSTAR, SPIN_NUE, SPIN_NUE_REF, SPIN_EP})});

        addProcessSignature(PID::BZERO, {PID::DSTARMINUS, PID::NU_E, PID::POSITRON}, {PID::DMINUS, PID::PI0});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCENU, FF_BDSTAR, SPIN_NUE, SPIN_NUE_REF, SPIN_EP})});

        addProcessSignature(PID::BZERO, {PID::DSTARMINUS, PID::NU_E, PID::POSITRON}, {-PID::D0, PID::PIMINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCENU, FF_BDSTAR, SPIN_NUE, SPIN_NUE_REF, SPIN_EP})});

        // bs -> cs
        addProcessSignature(PID::BS, {PID::DSSTARMINUS, PID::NU_TAU, PID::ANTITAU}, {PID::DSMINUS, PID::PI0});
        addTensor(Tensor{
            name, MD::makeEmptySparse(dims, {WILSON_BCTAUNU, FF_BSDSSTAR, SPIN_NUTAU, SPIN_NUTAU_REF, SPIN_TAUP})});

        addProcessSignature(PID::BS, {PID::DSSTARMINUS, PID::NU_MU, PID::ANTIMUON}, {PID::DSMINUS, PID::PI0});
        addTensor(
            Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCMUNU, FF_BSDSSTAR, SPIN_NUMU, SPIN_NUMU_REF, SPIN_MUP})});

        addProcessSignature(PID::BS, {PID::DSSTARMINUS, PID::NU_E, PID::POSITRON}, {PID::DSMINUS, PID::PI0});
        addTensor(
            Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCENU, FF_BSDSSTAR, SPIN_NUE, SPIN_NUE_REF, SPIN_EP})});

        setSignatureIndex();
    }

    void AmplBDstarDPiLepNu::eval(const Particle& parent, const ParticleList& daughters,
                                  const ParticleList& /*references*/) {
        // Momenta
        const FourMomentum& pBmes = parent.momentum();
        const FourMomentum& pDstarmes = daughters[0].momentum();
        const FourMomentum& kNuTau = daughters[1].momentum();
        const FourMomentum& pTau = daughters[2].momentum();
        const FourMomentum& pDmes = daughters[3].momentum();
        const FourMomentum& pPion = daughters[4].momentum();

        // kinematic objects
        const double Mb = pBmes.mass();
        const double Mb2 = Mb * Mb;
        const double Mds = pDstarmes.mass();
        const double Mds2 = Mds * Mds;
        const double Mt = pTau.mass();
        const double Md = pDmes.mass();
        const double Md2 = Md * Md;
        const double Mp = pPion.mass();
        const double Mp2 = Mp * Mp;
        const double Sqq = Mb2 + Mds2 - 2. * (pBmes * pDstarmes);
        const double sqSqq = sqrt(Sqq);
        const double Ew = (Mb2 - Mds2 + Sqq) / (2 * Mb);
        const double Eds = (Mb2 + Mds2 - Sqq) / (2 * Mb);
        const double Pw = sqrt(Ew * Ew - Sqq);
        const double Pw2 = Pw * Pw;
        const double BNuTau = (pBmes * kNuTau);
        const double NuTauQ = (pBmes * kNuTau) - (pDstarmes * kNuTau);
        const double BQ = Mb2 - (pBmes * pDstarmes);
        const double Ed = (Mds2 + Md2 - Mp2) / (2. * Mds);
        const double Ep = (Mds2 - Md2 + Mp2) / (2. * Mds);
        const double Pp = sqrt(Ed * Ed - Md2);
        const double DstarQ = pDstarmes * pBmes - Mds2;
        const double PionQ = (pPion * pBmes) - (pPion * pDstarmes);
        // const double PionD = (pPion * pDmes);
        const double PionNuTau = (pPion * kNuTau);
        const double TauNuTau = (pTau * kNuTau);
        const double epsBDPiNuTau = epsilon(pBmes, pDmes, pPion, kNuTau);

        // Helicity Angles
        const double CosTe = (Eds * Ep * Ew + Ep * Pw * Pw - Mds * PionQ) / (Mb * Pp * Pw);
        const double SinTe = sqrt(1. - CosTe * CosTe);
        // double CosTeHalf = pow((1. + CosTe)/2.,0.5);
        // double SinTeHalf = pow((1. - CosTe)/2.,0.5);
        const double CscTe = 1. / SinTe;

        const double CosTt = -((Ew * (Sqq * BNuTau - NuTauQ * BQ)) / (sqrt(Ew * Ew - Sqq) * NuTauQ * BQ));
        const double SinTt = sqrt(1. - CosTt * CosTt);
        const double CscTt = 1. / SinTt;
        // const double CosTtHalf = pow((1. + CosTt) / 2., 0.5);
        // const double SinTtHalf = pow((1. - CosTt) / 2., 0.5);
        const double CosTtHalfSq = (1. + CosTt) / 2.;
        const double SinTtHalfSq = (1. - CosTt) / 2.;

        const double CosPePt = -(CscTe * CscTt) / (Mds * Pp * sqSqq * NuTauQ) *
                               (Ep * (Sqq * BNuTau - NuTauQ * BQ) - Mds * (Sqq * PionNuTau - NuTauQ * PionQ) +
                                Pp * CosTe * CosTt * DstarQ * NuTauQ);
        const double SinPePt = (sqSqq * CscTe * CscTt * epsBDPiNuTau) / (Mb * Pp * Pw * TauNuTau);
        const complex<double> ExpIPePt = CosPePt + 1i * SinPePt;

        // Collective Functions
        const complex<double> Dpl = SinTe * (CosTtHalfSq / ExpIPePt + ExpIPePt * SinTtHalfSq);
        const complex<double> Dplc = SinTe * (CosTtHalfSq * ExpIPePt + SinTtHalfSq / ExpIPePt);
        const complex<double> Dmi = SinTe * (CosTtHalfSq / ExpIPePt - ExpIPePt * SinTtHalfSq);
        const complex<double> Dmic = SinTe * (CosTtHalfSq * ExpIPePt - SinTtHalfSq / ExpIPePt);
        const complex<double> Dze = CosTe * SinTt;

        const double Epl = CosTe * CosTt;
        const double EmiR = SinTe * SinTt * CosPePt;
        const double EmiI = SinTe * SinTt * SinPePt;
        const double Eze = CosTe;

        const double TwoSqTwoGfSqPp = TwoSqTwoGFermi * sqrt(Sqq - Mt * Mt) * Pp;

        const double gPiOnSqrt2MGam = sqrt(3. * pi * Mds / pow(Pp, 3.));

        // initialize tensor elements to zero
        Tensor& t = getTensor();
        t.clearData();

        // set non-zero tensor elements
        t.element({0, 1, 0, 0, 0}) =
            Mt * ((Eze * Mb * Pw) / (Mds * Sqq) - EmiR / sqSqq + (Epl * (-Mb2 + Mds2 + Sqq)) / (2. * Mds * Sqq));
        t.element({0, 2, 0, 0, 0}) = (2i * EmiI * Mb * Mt * Pw) / sqSqq;
        t.element({0, 3, 0, 0, 0}) = (Eze * Mb * Mt * Pw) / Mds;
        t.element({0, 4, 0, 0, 0}) =
            -(Mb * Mt * Pw * ((Eze * (-Mb2 + Mds2)) / (Mds * Sqq) + (2. * Epl * Mb * Pw) / (Mds * Sqq)));
        t.element({1, 0, 0, 0, 0}) = -((Eze * Mb * Pw) / Mds);
        t.element({3, 0, 0, 0, 0}) = (Eze * Mb * Pw) / Mds;
        t.element({5, 1, 0, 0, 0}) =
            -(Mt * ((Eze * Mb * Pw) / (Mds * Sqq) - EmiR / sqSqq + (Epl * (-Mb2 + Mds2 + Sqq)) / (2. * Mds * Sqq)));
        t.element({5, 2, 0, 0, 0}) = (2i * EmiI * Mb * Mt * Pw) / sqSqq;
        t.element({5, 3, 0, 0, 0}) = -((Eze * Mb * Mt * Pw) / Mds);
        t.element({5, 4, 0, 0, 0}) =
            Mb * Mt * Pw * ((Eze * (-Mb2 + Mds2)) / (Mds * Sqq) + (2. * Epl * Mb * Pw) / (Mds * Sqq));
        t.element({7, 1, 0, 0, 0}) =
            Mt * ((Eze * Mb * Pw) / (Mds * Sqq) - EmiR / sqSqq + (Epl * (-Mb2 + Mds2 + Sqq)) / (2. * Mds * Sqq));
        t.element({7, 2, 0, 0, 0}) = (2i * EmiI * Mb * Mt * Pw) / sqSqq;
        t.element({7, 3, 0, 0, 0}) = (Eze * Mb * Mt * Pw) / Mds;
        t.element({7, 4, 0, 0, 0}) =
            -(Mb * Mt * Pw * ((Eze * (-Mb2 + Mds2)) / (Mds * Sqq) + (2. * Epl * Mb * Pw) / (Mds * Sqq)));
        t.element({9, 5, 0, 0, 0}) = (8. * Epl * Mb2 * Pw2) / Mds;
        t.element({9, 6, 0, 0, 0}) = -2. * (-2. * EmiR * sqSqq + (Epl * (-Mb2 + Mds2 + Sqq)) / Mds);
        t.element({9, 7, 0, 0, 0}) = 2. * ((Epl * (Mb2 + 3. * Mds2 - Sqq)) / Mds + (2. * EmiR * (Mb2 - Mds2)) / sqSqq +
                                           (4i * EmiI * Mb * Pw) / sqSqq);
        t.element({0, 1, 0, 0, 1}) = -Dmi - (Dze * (-Mb2 + Mds2 + Sqq)) / (2. * Mds * sqSqq);
        t.element({0, 2, 0, 0, 1}) = -2. * Dpl * Mb * Pw;
        t.element({0, 4, 0, 0, 1}) = (2. * Dze * Mb2 * Pw2) / (Mds * sqSqq);
        t.element({5, 1, 0, 0, 1}) = Dmi + (Dze * (-Mb2 + Mds2 + Sqq)) / (2. * Mds * sqSqq);
        t.element({5, 2, 0, 0, 1}) = -2. * Dpl * Mb * Pw;
        t.element({5, 4, 0, 0, 1}) = (-2. * Dze * Mb2 * Pw2) / (Mds * sqSqq);
        t.element({7, 1, 0, 0, 1}) = -Dmi - (Dze * (-Mb2 + Mds2 + Sqq)) / (2. * Mds * sqSqq);
        t.element({7, 2, 0, 0, 1}) = -2. * Dpl * Mb * Pw;
        t.element({7, 4, 0, 0, 1}) = (2. * Dze * Mb2 * Pw2) / (Mds * sqSqq);
        t.element({9, 5, 0, 0, 1}) = (-8. * Dze * Mb2 * Mt * Pw2) / (Mds * sqSqq);
        t.element({9, 6, 0, 0, 1}) = 2. * Mt * (2. * Dmi + (Dze * (-Mb2 + Mds2 + Sqq)) / (Mds * sqSqq));
        t.element({9, 7, 0, 0, 1}) = -2. * Mt *
                                     ((-2. * Dmi * (Mb2 - Mds2)) / Sqq + (4. * Dpl * Mb * Pw) / Sqq +
                                      (Dze * (Mb2 + 3. * Mds2 - Sqq)) / (Mds * sqSqq));
        t.element({6, 1, 1, 1, 0}) = Dmic + (Dze * (-Mb2 + Mds2 + Sqq)) / (2. * Mds * sqSqq);
        t.element({6, 2, 1, 1, 0}) = 2. * Dplc * Mb * Pw;
        t.element({6, 4, 1, 1, 0}) = (-2. * Dze * Mb2 * Pw2) / (Mds * sqSqq);
        t.element({8, 1, 1, 1, 0}) = -Dmic - (Dze * (-Mb2 + Mds2 + Sqq)) / (2. * Mds * sqSqq);
        t.element({8, 2, 1, 1, 0}) = 2. * Dplc * Mb * Pw;
        t.element({8, 4, 1, 1, 0}) = (2. * Dze * Mb2 * Pw2) / (Mds * sqSqq);
        t.element({10, 5, 1, 1, 0}) = (8. * Dze * Mb2 * Mt * Pw2) / (Mds * sqSqq);
        t.element({10, 6, 1, 1, 0}) = 2. * Mt * (-2. * Dmic - (Dze * (-Mb2 + Mds2 + Sqq)) / (Mds * sqSqq));
        t.element({10, 7, 1, 1, 0}) = 2. * Mt *
                                      ((-2. * Dmic * (Mb2 - Mds2)) / Sqq + (4. * Dplc * Mb * Pw) / Sqq +
                                       (Dze * (Mb2 + 3. * Mds2 - Sqq)) / (Mds * sqSqq));
        t.element({2, 0, 1, 1, 1}) = (Eze * Mb * Pw) / Mds;
        t.element({4, 0, 1, 1, 1}) = -((Eze * Mb * Pw) / Mds);
        t.element({6, 1, 1, 1, 1}) =
            Mt * ((Eze * Mb * Pw) / (Mds * Sqq) - EmiR / sqSqq + (Epl * (-Mb2 + Mds2 + Sqq)) / (2. * Mds * Sqq));
        t.element({6, 2, 1, 1, 1}) = (-2i * EmiI * Mb * Mt * Pw) / sqSqq;
        t.element({6, 3, 1, 1, 1}) = (Eze * Mb * Mt * Pw) / Mds;
        t.element({6, 4, 1, 1, 1}) =
            -(Mb * Mt * Pw * ((Eze * (-Mb2 + Mds2)) / (Mds * Sqq) + (2. * Epl * Mb * Pw) / (Mds * Sqq)));
        t.element({8, 1, 1, 1, 1}) =
            -(Mt * ((Eze * Mb * Pw) / (Mds * Sqq) - EmiR / sqSqq + (Epl * (-Mb2 + Mds2 + Sqq)) / (2. * Mds * Sqq)));
        t.element({8, 2, 1, 1, 1}) = (-2i * EmiI * Mb * Mt * Pw) / sqSqq;
        t.element({8, 3, 1, 1, 1}) = -((Eze * Mb * Mt * Pw) / Mds);
        t.element({8, 4, 1, 1, 1}) =
            Mb * Mt * Pw * ((Eze * (-Mb2 + Mds2)) / (Mds * Sqq) + (2. * Epl * Mb * Pw) / (Mds * Sqq));
        t.element({10, 5, 1, 1, 1}) = (8. * Epl * Mb2 * Pw2) / Mds;
        t.element({10, 6, 1, 1, 1}) = -2. * (-2. * EmiR * sqSqq + (Epl * (-Mb2 + Mds2 + Sqq)) / Mds);
        t.element({10, 7, 1, 1, 1}) = 2. * ((Epl * (Mb2 + 3. * Mds2 - Sqq)) / Mds + (2. * EmiR * (Mb2 - Mds2)) / sqSqq -
                                            (4i * EmiI * Mb * Pw) / sqSqq);

        t *= TwoSqTwoGfSqPp * gPiOnSqrt2MGam;
    }

    void AmplBDstarDPiLepNu::addRefs() const {
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
