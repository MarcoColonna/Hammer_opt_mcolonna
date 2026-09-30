///
/// @file  AmplTauEllNuNu.cc
/// @brief \f$ \tau-> \ell\nu\nu \f$ amplitude
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/Amplitudes/AmplTauEllNuNu.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;
using namespace complex_literals;

namespace Hammer {

    namespace MD = MultiDimensional;

    AmplTauEllNuNu::AmplTauEllNuNu() {
        // Create tensor rank and dimensions
        IndexList dims{{2, 2}};
        string name{"AmplTauEllNuNu"};
        addProcessSignature(PID::ANTITAU, {PID::NU_TAUBAR, PID::NU_MU, PID::ANTIMUON});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {SPIN_NUTAU_REF, SPIN_TAUP})});

        addProcessSignature(PID::ANTITAU, {PID::NU_TAUBAR, PID::NU_E, PID::POSITRON});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {SPIN_NUTAU_REF, SPIN_TAUP})});

        setSignatureIndex();
        _multiplicity = 2ul;
    }

    void AmplTauEllNuNu::defineSettings() {
    }

    void AmplTauEllNuNu::eval(const Particle& parent, const ParticleList& daughters, const ParticleList& references) {
        // Momenta
        const FourMomentum& pTau = parent.momentum();
        const FourMomentum& kNuBarTau = daughters[0].momentum();
        const FourMomentum& kNuMuon = daughters[1].momentum();
        const FourMomentum& kMuon = daughters[2].momentum();

        // Siblings for parent tau case
        FourMomentum pDmes{1., 0., 0., 0.};
        FourMomentum kNuTau{1., 0., 0., 1.};
        if (references.size() >= 2) {
            pDmes = references[0].momentum();
            kNuTau = references[1].momentum();
        }
        const FourMomentum& pBmes = pDmes + kNuTau + pTau;

        // kinematic objects
        const double Mb = pBmes.mass();
        const double Mb2 = Mb * Mb;
        const double Md = pDmes.mass();
        const double Md2 = Md * Md;
        const double Mt = pTau.mass();
        const double Mt2 = Mt * Mt;
        const double Sqq = Mb2 + Md2 - 2. * (pBmes * pDmes);
        const double sqSqq = sqrt(Sqq);
        const double Ew = (Mb2 - Md2 + Sqq) / (2 * Mb);
        const double Pw = sqrt(Ew * Ew - Sqq);
        const double BNuTau = (pBmes * kNuTau);
        const double NuTauQ = (pBmes * kNuTau) - (pDmes * kNuTau);
        const double BQ = Mb2 - (pBmes * pDmes);


        const double Sqp = 2. * (kMuon * kNuMuon);
        const double sqSqp = sqrt(Sqp);
        const double BTau = (pBmes * pTau);
        const double MuonNuBarTau = (kMuon * kNuBarTau);
        const double NuMuonNuBarTau = (kNuMuon * kNuBarTau);
        const double NuMuonNuTau = (kNuMuon * kNuTau);
        const double NuTauNuBarTau = (kNuTau * kNuBarTau);
        const double TauNuBarTau = (pTau * kNuBarTau);
        const double TauNuTau = (pTau * kNuTau);
        const double BNuBarTau = (pBmes * kNuBarTau);
        const double epsBDNuTauNuBarTau = epsilon(pBmes, pDmes, kNuTau, kNuBarTau);
        const double epsNuTauNuBarTauEllonNuMuon = epsilon(kNuTau, kNuBarTau, kMuon, kNuMuon);

        const double CosTt = -((Ew * (Sqq * BNuTau - NuTauQ * BQ)) / (sqrt(Ew * Ew - Sqq) * NuTauQ * BQ));
        const double SinTt = sqrt(1. - CosTt * CosTt);
        const double CscTt = 1. / SinTt;

        const double CosTW = ((-(Mt2 * NuTauNuBarTau) + TauNuBarTau * TauNuTau) / (TauNuBarTau * TauNuTau));
        const double SinTW = sqrt(1. - CosTW * CosTW);
        const double CscTW = 1. / SinTW;
        const double CosTWHalf = pow((1. + CosTW) / 2., 0.5);
        const double SinTWHalf = pow((1. - CosTW) / 2., 0.5);
        const double TanTWHalf = SinTW / (CosTW + 1.);

        const double CosTm = (2. * (MuonNuBarTau - NuMuonNuBarTau)) / (Mt2 - Sqp);
        const double SinTm = sqrt(1. - CosTm * CosTm);
        const double CscTm = 1. / SinTm;
        const double CosTmHalf = pow((1. + CosTm) / 2., 0.5);
        // const double SinTmHalf = pow((1. - CosTm) / 2., 0.5);

        const double CosPtPW = (sqSqq * CscTt * CscTW) / (Mb * Mt * Pw * TauNuBarTau * TauNuTau) *
                               (TauNuTau * (Mt2 * BNuBarTau - BTau * TauNuBarTau) -
                                CosTW * TauNuBarTau * (Mt2 * BNuTau - BTau * TauNuTau));
        const double SinPtPW = -((sqSqq * CscTt * epsBDNuTauNuBarTau * TanTWHalf) / (Mb * Mt * Pw * NuTauNuBarTau));
        const complex<double> ExpIPtPW = CosPtPW + 1i * SinPtPW;

        const double CosPmPW = (CscTm * CscTW) / (Mt * sqSqp * TauNuTau) *
                               (Mt2 * (2 * NuMuonNuTau + (CosTm * CosTW - 1) * TauNuTau) +
                                (1 - CosTW) * (1 + CosTm) * TauNuTau * TauNuBarTau);
        const double SinPmPW = ((2. * CscTm * epsNuTauNuBarTauEllonNuMuon * TanTWHalf) / (Mt * sqSqp * NuTauNuBarTau));
        const complex<double> ExpIPmPW = CosPmPW + 1i * SinPmPW;

        const double TwoSqTwoGfSq = TwoSqTwoGFermi * (sqrt(Mt2 - Sqp));

        // initialize tensor elements to zero
        Tensor& t = getTensor();
        t.clearData();

        // set tensor elements
        t.element({0, 0}) = -(Mt * CosTWHalf * SinTm + 2. * ExpIPmPW * sqSqp * pow(CosTmHalf, 2.) * SinTWHalf);
        t.element({0, 1}) =
            (2. * ExpIPmPW * sqSqp * pow(CosTmHalf, 2.) * CosTWHalf - Mt * SinTm * SinTWHalf) / ExpIPtPW;
        t.element({1, 0}) = ExpIPtPW * t.element({0, 0});
        t.element({1, 1}) = ExpIPtPW * t.element({0, 1});

        t *= TwoSqTwoGfSq;
    }

} // namespace Hammer
