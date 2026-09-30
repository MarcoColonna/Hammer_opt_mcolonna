///
/// @file  AmplD1starDstarDGamPi.cc
/// @brief \f$ D_1^* \rightarrow D^* \pi, D^* \rightarrow D \gamma \f$ amplitude
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/Amplitudes/AmplD1starDstarDGamPi.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;
using namespace complex_literals;

namespace Hammer {

    namespace MD = MultiDimensional;

    AmplD1starDstarDGamPi::AmplD1starDstarDGamPi() {
        // Create tensor rank and dimensions
        IndexList dims = {1, 3, 2};
        string name{"AmplD1starDstarDGamPi"};

        addProcessSignature(-PID::DSSD1STAR, {PID::DSTARMINUS, PID::PIPLUS}, {PID::DMINUS, PID::GAMMA});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD1STAR, SPIN_DSSD1STAR, SPIN_GAMMA})});

        addProcessSignature(-PID::DSSD1STAR, {-PID::DSTAR, PID::PI0}, {-PID::D0, PID::GAMMA});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD1STAR, SPIN_DSSD1STAR, SPIN_GAMMA})});

        addProcessSignature(PID::DSSD1STARMINUS, {-PID::DSTAR, PID::PIMINUS}, {-PID::D0, PID::GAMMA});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD1STAR, SPIN_DSSD1STAR, SPIN_GAMMA})});

        addProcessSignature(PID::DSSD1STARMINUS, {PID::DSTARMINUS, PID::PI0}, {PID::DMINUS, PID::GAMMA});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD1STAR, SPIN_DSSD1STAR, SPIN_GAMMA})});

        // strange
        addProcessSignature(PID::DSSDS1STARMINUS, {PID::DSSTARMINUS, PID::PI0}, {PID::DSMINUS, PID::GAMMA});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSDS1STAR, SPIN_DSSDS1STAR, SPIN_GAMMA})});

        addProcessSignature(PID::DSSDS1STARMINUS, {PID::DSTARMINUS, -PID::K0}, {PID::DMINUS, PID::GAMMA});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSDS1STAR, SPIN_DSSDS1STAR, SPIN_GAMMA})});

        addProcessSignature(PID::DSSDS1STARMINUS, {-PID::DSTAR, PID::KMINUS}, {-PID::D0, PID::GAMMA});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSDS1STAR, SPIN_DSSDS1STAR, SPIN_GAMMA})});

        setSignatureIndex();
        _multiplicity = 3ul;
    }

    void AmplD1starDstarDGamPi::defineSettings() {
    }

    void AmplD1starDstarDGamPi::eval(const Particle& parent, const ParticleList& daughters,
                                     const ParticleList& references) {
        // Momenta
        const FourMomentum& pDssmes = parent.momentum();
        const FourMomentum& pDsmes = daughters[0].momentum();
        const FourMomentum& pPmes = daughters[1].momentum();
        const FourMomentum& pDmes = daughters[2].momentum();
        const FourMomentum& kGmes = daughters[3].momentum();

        // Siblings for parent D** case
        FourMomentum kNuTau{1., 1. / sqrt2, 0., 1. / sqrt2};
        FourMomentum pTau{1., -1. / sqrt2, 0., 1. / sqrt2};
        if (references.size() >= 2) {
            kNuTau = references[0].momentum();
            pTau = references[1].momentum();
        }
        const FourMomentum& Qmes = kNuTau + pTau;
        const FourMomentum& pBmes = pDssmes + pTau + kNuTau;

        // kinematic objects
        const double Mdss = pDssmes.mass();
        const double Mdss2 = Mdss * Mdss;
        const double Mds = pDsmes.mass();
        const double Mds2 = Mds * Mds;
        const double Md = pDmes.mass();
        const double Md2 = Md * Md;
        const double Mp = pPmes.mass();
        const double Mp2 = Mp * Mp;

        const double Eds = (Mdss2 + Mds2 - Mp2) / (2 * Mdss);
        const double Ep = (Mdss2 - Mds2 + Mp2) / (2 * Mdss);
        const double Pp = sqrt(Ep * Ep - Mp2);
        // const double Pp2 = Pp*Pp;

        const double Ed = (Mds2 + Md2) / (2 * Mds);
        const double Eg = (Mds2 - Md2) / (2 * Mds);

        const double Sqq = Qmes * Qmes;
        const double sqSqq = sqrt(Sqq);
        const double Mb = pBmes.mass();
        const double Mb2 = Mb * Mb;
        const double Ew = (Mb2 - Mdss2 + Sqq) / (2 * Mb);
        const double Edss = (Mb2 + Mdss2 - Sqq) / (2 * Mb);
        const double Pw = sqrt(Ew * Ew - Sqq);

        const double BNuTau = pBmes * kNuTau;
        const double NuTauQ = Qmes * kNuTau;
        const double BQ = pBmes * Qmes;

        // Helicity Angles
        const double CosTt = -((Ew * (Sqq * BNuTau - NuTauQ * BQ)) / (Pw * NuTauQ * BQ));
        const double SinTt = sqrt(1. - CosTt * CosTt);

        const double CosTds = (Eds * Edss * Mb - Mdss * (pBmes * pDsmes)) / (Mb * Pp * Pw);
        const double SinTds = sqrt(1. - CosTds * CosTds);
        const double CosTdsHalfSq = (1. + CosTds) / 2.;
        const double SinTdsHalfSq = (1. - CosTds) / 2.;
        //        const double CosTdsHalf = sqrt(CosTdsHalfSq);
        //        const double SinTdsHalf = sqrt(SinTdsHalfSq);
        //        const double CosTdsHalfCu = pow(CosTdsHalfSq, 1.5);
        //        const double SinTdsHalfCu = pow(SinTdsHalfSq, 1.5);

        const double CosTd = (-Ed * Eds * Mdss + Mds * (pDssmes * pDmes)) / (Mdss * Pp * Eg);
        const double SinTd = sqrt(1. - CosTd * CosTd);
        const double CosTdHalfSq = (1. + CosTd) / 2.;
        const double SinTdHalfSq = (1. - CosTd) / 2.;

        const double CosPdsPt =
            (Ep * Sqq * (pDsmes * kNuTau) - Eds * Sqq * (pPmes * kNuTau) +
             NuTauQ * (Pp * (-Ew * Mb + Sqq) * CosTds * CosTt - Ep * (pDsmes * Qmes) + Eds * (pPmes * Qmes))) /
            (SinTds * SinTt * Mdss * Pp * sqSqq * NuTauQ);
        const double SinPdsPt =
            -sqSqq * epsilon(pBmes, pDsmes, pDssmes, kNuTau) / (SinTds * SinTt * Mb * Pp * Pw * NuTauQ);
        const complex<double> ExpIPdsPt = CosPdsPt + 1i * SinPdsPt;
        const complex<double> ExpIPtPds = CosPdsPt - 1i * SinPdsPt;
        //        const complex<double> Exp2IPdsPt = ExpIPdsPt * ExpIPdsPt;
        //        const complex<double> Exp2IPtPds = ExpIPtPds * ExpIPtPds;

        const double CosPdPds = (Eg * Mdss2 * (pBmes * pDmes) - Ed * Mdss2 * (pBmes * kGmes) +
                                 Mb * (Eds * Mdss * Eg * Pw * CosTd * CosTds - Edss * Eg * (pDssmes * pDmes) +
                                       Ed * Edss * (pDssmes * kGmes))) /
                                (SinTds * SinTd * Mb * Mds * Mdss * Eg * Pw);
        const double SinPdPds = -epsilon(pBmes, pDsmes, pDssmes, pDmes) / (SinTds * SinTd * Mb * Pp * Eg * Pw);
        const complex<double> ExpIPdPds = CosPdPds + 1i * SinPdPds;
        const complex<double> ExpIPdsPd = CosPdPds - 1i * SinPdPds;


        const complex<double> prefactor = 1i * Mds * Eg * Mdss * (Mdss - Mds);       // Include prefactor ~ Mdss Epi
        const double eMuOnSqrt2MGam = sqrt(48. * pi * Mds2 / pow(2 * Mds * Eg, 3.)); // for D* decay

        // initialize tensor elements to zero
        Tensor& t = getTensor();
        t.clearData();

        // set non-zero tensor elements
        t.element({0, 0, 0}) = (ExpIPtPds * (-2. * CosTdHalfSq * CosTdsHalfSq * ExpIPdsPd + SinTd * SinTds +
                                             (-1 + CosTd) * ExpIPdPds * SinTdsHalfSq)) /
                               2.;
        t.element({0, 0, 1}) = (ExpIPtPds * ((-1 + CosTd) * CosTdsHalfSq * ExpIPdsPd - SinTd * SinTds -
                                             2. * CosTdHalfSq * ExpIPdPds * SinTdsHalfSq)) /
                               2.;
        t.element({0, 1, 0}) =
            -(2. * CosTds * SinTd + (ExpIPdsPd + CosTd * ExpIPdsPd - 2. * ExpIPdPds * SinTdHalfSq) * SinTds) /
            (2. * sqrt2);
        t.element({0, 1, 1}) =
            (2. * CosTds * SinTd + (1 + CosTd) * ExpIPdPds * SinTds + (-1 + CosTd) * ExpIPdsPd * SinTds) / (2. * sqrt2);
        t.element({0, 2, 0}) = CosTdsHalfSq * ExpIPdPds * ExpIPdsPt * SinTdHalfSq + (ExpIPdsPt * SinTd * SinTds) / 2. +
                               CosTdHalfSq * ExpIPdsPd * ExpIPdsPt * SinTdsHalfSq;
        t.element({0, 2, 1}) = CosTdHalfSq * CosTdsHalfSq * ExpIPdPds * ExpIPdsPt - (ExpIPdsPt * SinTd * SinTds) / 2. +
                               ExpIPdsPd * ExpIPdsPt * SinTdHalfSq * SinTdsHalfSq;

        t *= prefactor * eMuOnSqrt2MGam;
    }


} // namespace Hammer
