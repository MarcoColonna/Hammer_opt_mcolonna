///
/// @file  AmplTau3PiNu.cc
/// @brief \f$ \tau-> \pi\pi\pi\nu \f$ amplitude
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/Amplitudes/AmplTau3PiNu.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;
using namespace complex_literals;

namespace Hammer {

    namespace MD = MultiDimensional;

    AmplTau3PiNu::AmplTau3PiNu() {
        // Create tensor rank and dimensions
        IndexList dims{{2, 2, 3}};
        string name{"AmplTau3PiNu"};
        // Ordering is piminus then piplus as anti-tau parent has pdg < 0
        addProcessSignature(PID::ANTITAU, {PID::PIMINUS, PID::PIPLUS, PID::PIPLUS, PID::NU_TAUBAR});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {SPIN_NUTAU_REF, SPIN_TAUP, FF_TAU3PI})});
        //        ///@todo check the amplitude is correct with this ordering below
        //        The 2 pi^0 amplitude has different symmetry factors
        //        addProcessSignature(PID::ANTITAU, {PID::PIPLUS, PID::PI0, PID::PI0, PID::NU_TAUBAR});
        //        addTensor(Tensor{name, MD::makeEmptySparse(dims, {SPIN_NUTAU_REF, SPIN_TAUP, FF_TAU3PI})});

        setSignatureIndex();
        _multiplicity = 2ul;
    }

    void AmplTau3PiNu::defineSettings() {
    }

    void AmplTau3PiNu::eval(const Particle& parent, const ParticleList& daughters, const ParticleList& references) {
        // Momenta
        const FourMomentum& pTau = parent.momentum();
        const FourMomentum& kNuBarTau = daughters[3].momentum();

        FourMomentum pPiPlus1;
        FourMomentum pPiPlus2;
        FourMomentum pPiMinus;
        if (daughters[0].pdgId() == daughters[1].pdgId()) { // Ordered 211 211 -211
            pPiPlus1 = daughters[0].momentum();
            pPiPlus2 = daughters[1].momentum();
            pPiMinus = daughters[2].momentum();
        } else if (daughters[1].pdgId() == daughters[2].pdgId()) { // 211 -211 -211
            pPiPlus1 = daughters[1].momentum();
            pPiPlus2 = daughters[2].momentum();
            pPiMinus = daughters[0].momentum();
        }


        // Siblings for parent tau case
        FourMomentum pDmes{1., 0., 0., 0.};
        FourMomentum kNuTau{1., 0., 0., 1.};
        if (references.size() >= 2) {
            pDmes = references[0].momentum();
            kNuTau = references[1].momentum();
        }
        const FourMomentum& pBmes = pDmes + kNuTau + pTau;

        const FourMomentum p = pTau - kNuBarTau;
        // Pion subtracted momenta
        const FourMomentum p1 = pPiPlus1 - pPiMinus;
        const FourMomentum p2 = pPiPlus2 - pPiMinus;

        // kinematic objects
        const double Mb = pBmes.mass();
        const double Mb2 = Mb * Mb;
        const double Md = pDmes.mass();
        const double Md2 = Md * Md;
        const double Mt = pTau.mass();
        const double Mt2 = Mt * Mt;
        const double t1 = p1.mass2();
        const double t2 = p2.mass2();
        const double Sqp = p.mass2();
        const double sqSqp = sqrt(Sqp);
        const double E1 = (p * p1) / sqSqp;
        const double E12 = E1 * E1;
        const double sqE12t1 = sqrt(E12 - t1);
        const double E2 = (p * p2) / sqSqp;
        const double E22 = E2 * E2;
        const double sqE22t2 = sqrt(E22 - t2);

        const double Sqq = Mb2 + Md2 - 2. * (pBmes * pDmes);
        const double sqSqq = sqrt(Sqq);
        const double Ew = (Mb2 - Md2 + Sqq) / (2 * Mb);
        const double Pw = sqrt(Ew * Ew - Sqq);
        const double NuTauNuBarTau = (kNuTau * kNuBarTau);
        const double TauNuBarTau = (pTau * kNuBarTau);
        const double TauNuTau = (pTau * kNuTau);
        const double p1Tau = (p1 * pTau);
        const double p2Tau = (p2 * pTau);
        const double p1NuTau = (p1 * kNuTau);
        const double p2NuTau = (p2 * kNuTau);
        const double p1NuBarTau = (p1 * kNuBarTau);
        const double p2NuBarTau = (p2 * kNuBarTau);
        const double p1p = (p1 * p);
        const double p2p = (p2 * p);
        const double BNuTau = (pBmes * kNuTau);
        const double BNuBarTau = (pBmes * kNuBarTau);
        const double NuTauQ = (pBmes * kNuTau) - (pDmes * kNuTau);
        const double BQ = Mb2 - (pBmes * pDmes);
        const double BTau = (pBmes * pTau);
        const double epsBDNuTauNuBarTau = epsilon(pBmes, pDmes, kNuTau, kNuBarTau);

        // Helicity Angles
        const double CosTt = -((Ew * (Sqq * BNuTau - NuTauQ * BQ)) / (sqrt(Ew * Ew - Sqq) * NuTauQ * BQ));
        const double SinTt = sqrt(1. - CosTt * CosTt);
        const double CscTt = 1. / SinTt;

        const double CosTW = ((-(Mt2 * NuTauNuBarTau) + TauNuBarTau * TauNuTau) / (TauNuBarTau * TauNuTau));
        const double SinTW = sqrt(1. - CosTW * CosTW);
        const double CscTW = 1. / SinTW;
        const double CosTWHalf = pow((1. + CosTW) / 2., 0.5);
        const double SinTWHalf = pow((1. - CosTW) / 2., 0.5);
        const double TanTWHalf = SinTW / (CosTW + 1.);

        const double CosT1 = (Sqp * p1NuBarTau - p1p * TauNuBarTau) / (sqSqp * sqE12t1 * TauNuBarTau);
        const double SinT1 = sqrt(1. - CosT1 * CosT1);
        const double CscT1 = 1. / SinT1;

        const double CosT2 = (Sqp * p2NuBarTau - p2p * TauNuBarTau) / (sqSqp * sqE22t2 * TauNuBarTau);
        const double SinT2 = sqrt(1. - CosT2 * CosT2);
        const double CscT2 = 1. / SinT2;

        const double CosP1PW = (CscT1 * CscTW *
                                (-(Mt2 * p1NuTau) + sqSqp * sqE12t1 * CosT1 * CosTW * TauNuTau +
                                 CosTW * p1NuBarTau * TauNuTau + p1Tau * TauNuTau)) /
                               (Mt * sqE12t1 * TauNuTau);
        const double SinP1PW =
            -(Mt * CscT1 * CscTW * epsilon(p1, pTau, kNuTau, p)) / (sqE12t1 * TauNuBarTau * TauNuTau);
        const complex<double> ExpIP1PW = CosP1PW + 1i * SinP1PW;

        const double CosP2PW = (CscT2 * CscTW *
                                (-(Mt2 * p2NuTau) + sqSqp * sqE22t2 * CosT2 * CosTW * TauNuTau +
                                 CosTW * p2NuBarTau * TauNuTau + p2Tau * TauNuTau)) /
                               (Mt * sqE22t2 * TauNuTau);
        const double SinP2PW =
            -(Mt * CscT2 * CscTW * epsilon(p2, pTau, kNuTau, p)) / (sqE22t2 * TauNuBarTau * TauNuTau);
        const complex<double> ExpIP2PW = CosP2PW + 1i * SinP2PW;

        const double CosPtPW = (sqSqq * CscTt * CscTW) / (Mb * Mt * Pw * TauNuBarTau * TauNuTau) *
                               (TauNuTau * (Mt2 * BNuBarTau - BTau * TauNuBarTau) -
                                CosTW * TauNuBarTau * (Mt2 * BNuTau - BTau * TauNuTau));
        const double SinPtPW = -((sqSqq * CscTt * epsBDNuTauNuBarTau * TanTWHalf) / (Mb * Mt * Pw * NuTauNuBarTau));
        const complex<double> ExpIPtPW = CosPtPW + 1i * SinPtPW;

        const double SqTwoGfSq = sqrt2 * GFermi * (sqrt(Mt2 - Sqp));

        //        // Structure Functions
        //        double F1 = 1.;
        //        double F2 = 1.;
        //        double F4 = 1.;

        // initialize tensor elements to zero
        Tensor& t = getTensor();
        t.clearData();

        // set tensor elements
        t.element({0, 0, 0}) = sqE12t1 / sqSqp * (Mt * CosT1 * CosTWHalf - ExpIP1PW * sqSqp * SinT1 * SinTWHalf);
        t.element({0, 0, 1}) = sqE22t2 / sqSqp * (Mt * CosT2 * CosTWHalf - ExpIP2PW * sqSqp * SinT2 * SinTWHalf);
        t.element({0, 0, 2}) = Mt * CosTWHalf;
        t.element({0, 1, 0}) =
            (sqE12t1 / sqSqp * (ExpIP1PW * sqSqp * CosTWHalf * SinT1 + Mt * CosT1 * SinTWHalf)) / ExpIPtPW;
        t.element({0, 1, 1}) =
            (sqE22t2 / sqSqp * (ExpIP2PW * sqSqp * CosTWHalf * SinT2 + Mt * CosT2 * SinTWHalf)) / ExpIPtPW;
        t.element({0, 1, 2}) = (Mt * SinTWHalf) / ExpIPtPW;
        t.element({1, 0, 0}) =
            ExpIPtPW * sqE12t1 / sqSqp * (Mt * CosT1 * CosTWHalf - ExpIP1PW * sqSqp * SinT1 * SinTWHalf);
        t.element({1, 0, 1}) =
            ExpIPtPW * sqE22t2 / sqSqp * (Mt * CosT2 * CosTWHalf - ExpIP2PW * sqSqp * SinT2 * SinTWHalf);
        t.element({1, 0, 2}) = ExpIPtPW * Mt * CosTWHalf;
        t.element({1, 1, 0}) = sqE12t1 / sqSqp * (ExpIP1PW * sqSqp * CosTWHalf * SinT1 + Mt * CosT1 * SinTWHalf);
        t.element({1, 1, 1}) = sqE22t2 / sqSqp * (ExpIP2PW * sqSqp * CosTWHalf * SinT2 + Mt * CosT2 * SinTWHalf);
        t.element({1, 1, 2}) = Mt * SinTWHalf;

        //        t.element({0, 0}) = F4 * Mt * CosTWHalf +
        //                          F1 * sqE12t1 / sqSqp * (Mt * CosT1 * CosTWHalf - ExpIP1PW * sqSqp * SinT1 *
        //                          SinTWHalf) + F2 * sqE22t2 / sqSqp * (Mt * CosT2 * CosTWHalf - ExpIP2PW * sqSqp *
        //                          SinT2 * SinTWHalf);
        //        t.element({0, 1}) =
        //            (F4 * Mt * SinTWHalf) / ExpIPtPW +
        //            (F1 * sqE12t1 / sqSqp * (ExpIP1PW * sqSqp * CosTWHalf * SinT1 + Mt * CosT1 * SinTWHalf)) /
        //            ExpIPtPW + (F2 * sqE22t2 / sqSqp * (ExpIP2PW * sqSqp * CosTWHalf * SinT2 + Mt * CosT2 *
        //            SinTWHalf)) / ExpIPtPW;
        //        t.element({1, 0}) =
        //            ExpIPtPW * F4 * Mt * CosTWHalf +
        //            ExpIPtPW * F1 * sqE12t1 / sqSqp * (Mt * CosT1 * CosTWHalf - ExpIP1PW * sqSqp * SinT1 * SinTWHalf)
        //            + ExpIPtPW * F2 * sqE22t2 / sqSqp * (Mt * CosT2 * CosTWHalf - ExpIP2PW * sqSqp * SinT2 *
        //            SinTWHalf);
        //        t.element({1, 1}) = F4 * Mt * SinTWHalf +
        //                          F1 * sqE12t1 / sqSqp * (ExpIP1PW * sqSqp * CosTWHalf * SinT1 + Mt * CosT1 *
        //                          SinTWHalf) + F2 * sqE22t2 / sqSqp * (ExpIP2PW * sqSqp * CosTWHalf * SinT2 + Mt *
        //                          CosT2 * SinTWHalf);

        t *= SqTwoGfSq;
    }

} // namespace Hammer
