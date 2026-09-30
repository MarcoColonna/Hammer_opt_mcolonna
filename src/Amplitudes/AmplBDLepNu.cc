///
/// @file  AmplBDLepNu.cc
/// @brief \f$ B \rightarrow D \tau\nu \f$ amplitude
/// @brief Also: \f$ B \rightarrow \pi \tau\nu \f$ amplitude
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/Amplitudes/AmplBDLepNu.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;
using namespace complex_literals;

namespace Hammer {

    namespace MD = MultiDimensional;

    AmplBDLepNu::AmplBDLepNu() {
        // Create tensor rank and dimensions

        IndexList dims{{11, 4, 2, 2, 2}};
        string name{"AmplBDLepNu"};
        addProcessSignature(PID::BPLUS, {-PID::D0, PID::NU_TAU, PID::ANTITAU});
        addTensor(
            Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCTAUNU, FF_BD, SPIN_NUTAU, SPIN_NUTAU_REF, SPIN_TAUP})});

        addProcessSignature(PID::BZERO, {PID::DMINUS, PID::NU_TAU, PID::ANTITAU});
        addTensor(
            Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCTAUNU, FF_BD, SPIN_NUTAU, SPIN_NUTAU_REF, SPIN_TAUP})});

        addProcessSignature(PID::BPLUS, {-PID::D0, PID::NU_MU, PID::ANTIMUON});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCMUNU, FF_BD, SPIN_NUMU, SPIN_NUMU_REF, SPIN_MUP})});

        addProcessSignature(PID::BZERO, {PID::DMINUS, PID::NU_MU, PID::ANTIMUON});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCMUNU, FF_BD, SPIN_NUMU, SPIN_NUMU_REF, SPIN_MUP})});

        addProcessSignature(PID::BPLUS, {-PID::D0, PID::NU_E, PID::POSITRON});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCENU, FF_BD, SPIN_NUE, SPIN_NUE_REF, SPIN_EP})});

        addProcessSignature(PID::BZERO, {PID::DMINUS, PID::NU_E, PID::POSITRON});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCENU, FF_BD, SPIN_NUE, SPIN_NUE_REF, SPIN_EP})});

        // bs -> cs
        addProcessSignature(PID::BS, {PID::DSMINUS, PID::NU_TAU, PID::ANTITAU});
        addTensor(
            Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCTAUNU, FF_BSDS, SPIN_NUTAU, SPIN_NUTAU_REF, SPIN_TAUP})});

        addProcessSignature(PID::BS, {PID::DSMINUS, PID::NU_MU, PID::ANTIMUON});
        addTensor(
            Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCMUNU, FF_BSDS, SPIN_NUMU, SPIN_NUMU_REF, SPIN_MUP})});

        addProcessSignature(PID::BS, {PID::DSMINUS, PID::NU_E, PID::POSITRON});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCENU, FF_BSDS, SPIN_NUE, SPIN_NUE_REF, SPIN_EP})});

        // b -> u
        addProcessSignature(PID::BPLUS, {PID::PI0, PID::NU_TAU, PID::ANTITAU});
        addTensor(
            Tensor{name, MD::makeEmptySparse(dims, {WILSON_BUTAUNU, FF_BPI, SPIN_NUTAU, SPIN_NUTAU_REF, SPIN_TAUP})});

        addProcessSignature(PID::BZERO, {PID::PIMINUS, PID::NU_TAU, PID::ANTITAU});
        addTensor(
            Tensor{name, MD::makeEmptySparse(dims, {WILSON_BUTAUNU, FF_BPI, SPIN_NUTAU, SPIN_NUTAU_REF, SPIN_TAUP})});

        addProcessSignature(PID::BPLUS, {PID::PI0, PID::NU_MU, PID::ANTIMUON});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BUMUNU, FF_BPI, SPIN_NUMU, SPIN_NUMU_REF, SPIN_MUP})});

        addProcessSignature(PID::BZERO, {PID::PIMINUS, PID::NU_MU, PID::ANTIMUON});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BUMUNU, FF_BPI, SPIN_NUMU, SPIN_NUMU_REF, SPIN_MUP})});

        addProcessSignature(PID::BPLUS, {PID::PI0, PID::NU_E, PID::POSITRON});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BUENU, FF_BPI, SPIN_NUE, SPIN_NUE_REF, SPIN_EP})});

        addProcessSignature(PID::BZERO, {PID::PIMINUS, PID::NU_E, PID::POSITRON});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BUENU, FF_BPI, SPIN_NUE, SPIN_NUE_REF, SPIN_EP})});

        // bs -> us
        addProcessSignature(PID::BS, {PID::KMINUS, PID::NU_TAU, PID::ANTITAU});
        addTensor(
            Tensor{name, MD::makeEmptySparse(dims, {WILSON_BUTAUNU, FF_BSK, SPIN_NUTAU, SPIN_NUTAU_REF, SPIN_TAUP})});

        addProcessSignature(PID::BS, {PID::KMINUS, PID::NU_MU, PID::ANTIMUON});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BUMUNU, FF_BSK, SPIN_NUMU, SPIN_NUMU_REF, SPIN_MUP})});

        addProcessSignature(PID::BS, {PID::KMINUS, PID::NU_E, PID::POSITRON});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BUENU, FF_BSK, SPIN_NUE, SPIN_NUE_REF, SPIN_EP})});


        setSignatureIndex();
    }

    void AmplBDLepNu::eval(const Particle& parent, const ParticleList& daughters, const ParticleList& /*references*/) {
        // Momenta
        const FourMomentum& pBmes = parent.momentum();
        const FourMomentum& pDmes = daughters[0].momentum();
        const FourMomentum& kNuTau = daughters[1].momentum();
        const FourMomentum& pTau = daughters[2].momentum();

        // combinatoric for neutral unflavored mesons
        const double Cv = (daughters[0].pdgId() == PID::PI0) ? 1. / sqrt2 : 1.;

        // kinematic objects
        const double Mb = pBmes.mass();
        const double Mb2 = Mb * Mb;
        const double Md = pDmes.mass();
        const double Md2 = Md * Md;
        const double Mt = pTau.mass();
        const double Sqq = Mb2 + Md2 - 2. * (pBmes * pDmes);
        const double sqSqq = sqrt(Sqq);
        const double Ew = (Mb2 - Md2 + Sqq) / (2 * Mb);
        const double Pw = sqrt(Ew * Ew - Sqq);
        const double BNuTau = (pBmes * kNuTau);
        const double NuTauQ = (pBmes * kNuTau) - (pDmes * kNuTau);
        const double BQ = Mb2 - (pBmes * pDmes);

        // Helicity Angles
        const double CosTt = -((Ew * (Sqq * BNuTau - NuTauQ * BQ)) / (sqrt(Ew * Ew - Sqq) * NuTauQ * BQ));
        const double SinTt = sqrt(1. - CosTt * CosTt);
        // double CscTt = 1. / SinTt;

        const double TwoSqTwoGfSq = TwoSqTwoGFermi * Cv * sqrt(Sqq - Mt * Mt);

        // initialize tensor elements to zero
        Tensor& t = getTensor();
        t.clearData();

        // set non-zero tensor elements
        t.element({0, 1, 0, 0, 0}) = ((-Mb2 + Md2) * Mt) / (2. * Sqq);
        t.element({0, 2, 0, 0, 0}) = (Mb * Mt * Pw * CosTt) / Sqq;
        t.element({1, 0, 0, 0, 0}) = 0.5;
        t.element({3, 0, 0, 0, 0}) = 0.5;
        t.element({5, 1, 0, 0, 0}) = ((-Mb2 + Md2) * Mt) / (2. * Sqq);
        t.element({5, 2, 0, 0, 0}) = (Mb * Mt * Pw * CosTt) / Sqq;
        t.element({7, 1, 0, 0, 0}) = ((-Mb2 + Md2) * Mt) / (2. * Sqq);
        t.element({7, 2, 0, 0, 0}) = (Mb * Mt * Pw * CosTt) / Sqq;
        t.element({9, 3, 0, 0, 0}) = -4. * Mb * Pw * CosTt;
        t.element({0, 2, 0, 0, 1}) = -((Mb * Pw * SinTt) / sqSqq);
        t.element({5, 2, 0, 0, 1}) = -((Mb * Pw * SinTt) / sqSqq);
        t.element({7, 2, 0, 0, 1}) = -((Mb * Pw * SinTt) / sqSqq);
        t.element({9, 3, 0, 0, 1}) = (4. * Mb * Mt * Pw * SinTt) / sqSqq;
        t.element({6, 2, 1, 1, 0}) = -((Mb * Pw * SinTt) / sqSqq);
        t.element({8, 2, 1, 1, 0}) = -((Mb * Pw * SinTt) / sqSqq);
        t.element({10, 3, 1, 1, 0}) = (4. * Mb * Mt * Pw * SinTt) / sqSqq;
        t.element({2, 0, 1, 1, 1}) = -0.5;
        t.element({4, 0, 1, 1, 1}) = -0.5;
        t.element({6, 1, 1, 1, 1}) = ((Mb2 - Md2) * Mt) / (2. * Sqq);
        t.element({6, 2, 1, 1, 1}) = -((Mb * Mt * Pw * CosTt) / Sqq);
        t.element({8, 1, 1, 1, 1}) = ((Mb2 - Md2) * Mt) / (2. * Sqq);
        t.element({8, 2, 1, 1, 1}) = -((Mb * Mt * Pw * CosTt) / Sqq);
        t.element({10, 3, 1, 1, 1}) = 4. * Mb * Pw * CosTt;

        t *= TwoSqTwoGfSq;
    }

    void AmplBDLepNu::addRefs() const {
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
