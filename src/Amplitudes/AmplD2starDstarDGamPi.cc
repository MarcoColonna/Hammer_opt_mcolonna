///
/// @file  AmplD2starDstarDGamPi.cc
/// @brief \f$ D_2^* \rightarrow D^* \pi, D^* \rightarrow D \gamma \f$ amplitude
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/Amplitudes/AmplD2starDstarDGamPi.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;
using namespace complex_literals;

namespace Hammer {

    namespace MD = MultiDimensional;

    AmplD2starDstarDGamPi::AmplD2starDstarDGamPi() {
        // Create tensor rank and dimensions
        IndexList dims = {1, 5, 2};
        string name{"AmplD2starDstarDGamPi"};

        addProcessSignature(-PID::DSSD2STAR, {PID::DSTARMINUS, PID::PIPLUS}, {PID::DMINUS, PID::GAMMA});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD2STAR, SPIN_DSSD2STAR, SPIN_GAMMA})});

        addProcessSignature(-PID::DSSD2STAR, {-PID::DSTAR, PID::PI0}, {-PID::D0, PID::GAMMA});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD2STAR, SPIN_DSSD2STAR, SPIN_GAMMA})});

        addProcessSignature(PID::DSSD2STARMINUS, {-PID::DSTAR, PID::PIMINUS}, {-PID::D0, PID::GAMMA});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD2STAR, SPIN_DSSD2STAR, SPIN_GAMMA})});

        addProcessSignature(PID::DSSD2STARMINUS, {PID::DSTARMINUS, PID::PI0}, {PID::DMINUS, PID::GAMMA});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD2STAR, SPIN_DSSD2STAR, SPIN_GAMMA})});

        // strange
        addProcessSignature(PID::DSSDS2STARMINUS, {PID::DSSTARMINUS, PID::PI0}, {PID::DSMINUS, PID::GAMMA});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSDS2STAR, SPIN_DSSDS2STAR, SPIN_GAMMA})});

        addProcessSignature(PID::DSSDS2STARMINUS, {PID::DSTARMINUS, -PID::K0}, {PID::DMINUS, PID::GAMMA});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSDS2STAR, SPIN_DSSDS2STAR, SPIN_GAMMA})});

        addProcessSignature(PID::DSSDS2STARMINUS, {-PID::DSTAR, PID::KMINUS}, {-PID::D0, PID::GAMMA});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSDS2STAR, SPIN_DSSDS2STAR, SPIN_GAMMA})});

        setSignatureIndex();
        _multiplicity = 5ul;
    }

    void AmplD2starDstarDGamPi::defineSettings() {
    }

    void AmplD2starDstarDGamPi::eval(const Particle& parent, const ParticleList& daughters,
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
        const double Pp2 = Pp * Pp;

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
        const double CosTdsHalf = sqrt(CosTdsHalfSq);
        const double SinTdsHalf = sqrt(SinTdsHalfSq);
        const double CosTdsHalfCu = pow(CosTdsHalfSq, 1.5);
        const double SinTdsHalfCu = pow(SinTdsHalfSq, 1.5);

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
        const complex<double> Exp2IPdsPt = ExpIPdsPt * ExpIPdsPt;
        const complex<double> Exp2IPtPds = ExpIPtPds * ExpIPtPds;

        const double CosPdPds = (Eg * Mdss2 * (pBmes * pDmes) - Ed * Mdss2 * (pBmes * kGmes) +
                                 Mb * (Eds * Mdss * Eg * Pw * CosTd * CosTds - Edss * Eg * (pDssmes * pDmes) +
                                       Ed * Edss * (pDssmes * kGmes))) /
                                (SinTds * SinTd * Mb * Mds * Mdss * Eg * Pw);
        const double SinPdPds = -epsilon(pBmes, pDsmes, pDssmes, pDmes) / (SinTds * SinTd * Mb * Pp * Eg * Pw);
        const complex<double> ExpIPdPds = CosPdPds + 1i * SinPdPds;
        const complex<double> ExpIPdsPd = CosPdPds - 1i * SinPdPds;


        const complex<double> prefactor = 2i * Mds * Pp2 * Eg / sqrt2;
        const double eMuOnSqrt2MGam = sqrt(48. * pi * Mds2 / pow(2 * Mds * Eg, 3.)); // for D* decay

        // initialize tensor elements to zero
        Tensor& t = getTensor();
        t.clearData();

        // set non-zero tensor elements
        t.element({0, 0, 0}) = -(CosTdHalfSq * CosTdsHalfCu * Exp2IPtPds * ExpIPdsPd * SinTdsHalf) +
                               CosTdsHalf * Exp2IPtPds * ExpIPdPds * SinTdHalfSq * SinTdsHalfCu;
        t.element({0, 0, 1}) = ((-1 + CosTd) * CosTdsHalfCu * Exp2IPtPds * ExpIPdsPd * SinTdsHalf) / 2. +
                               CosTdHalfSq * CosTdsHalf * Exp2IPtPds * ExpIPdPds * SinTdsHalfCu;
        t.element({0, 1, 0}) = (ExpIPtPds * (CosTdHalfSq * (-1 + 2. * CosTds) * CosTdsHalfSq * ExpIPdsPd -
                                             (1 + 2. * CosTds) * ExpIPdPds * SinTdHalfSq * SinTdsHalfSq)) /
                               2.;
        t.element({0, 1, 1}) = (ExpIPtPds * ((-1 + 2. * CosTds) * CosTdsHalfSq * ExpIPdsPd * SinTdHalfSq -
                                             CosTdHalfSq * (1 + 2. * CosTds) * ExpIPdPds * SinTdsHalfSq)) /
                               2.;
        t.element({0, 2, 0}) =
            -(CosTds * CosTdsHalf * (CosTdHalfSq * ExpIPdsPd + ExpIPdPds * SinTdHalfSq) * SinTdsHalf * sqrt(3. / 2));
        t.element({0, 2, 1}) =
            -(CosTds * CosTdsHalf * (CosTdHalfSq * ExpIPdPds + ExpIPdsPd * SinTdHalfSq) * SinTdsHalf * sqrt(3. / 2));
        t.element({0, 3, 0}) = (ExpIPdsPt * ((-1 + 2. * CosTds) * CosTdsHalfSq * ExpIPdPds * SinTdHalfSq -
                                             CosTdHalfSq * (1 + 2. * CosTds) * ExpIPdsPd * SinTdsHalfSq)) /
                               2.;
        t.element({0, 3, 1}) = (ExpIPdsPt * (CosTdHalfSq * (-1 + 2. * CosTds) * CosTdsHalfSq * ExpIPdPds -
                                             (1 + 2. * CosTds) * ExpIPdsPd * SinTdHalfSq * SinTdsHalfSq)) /
                               2.;
        t.element({0, 4, 0}) = ((-1 + CosTd) * CosTdsHalfCu * Exp2IPdsPt * ExpIPdPds * SinTdsHalf) / 2. +
                               CosTdHalfSq * CosTdsHalf * Exp2IPdsPt * ExpIPdsPd * SinTdsHalfCu;
        t.element({0, 4, 1}) = -(CosTdHalfSq * CosTdsHalfCu * Exp2IPdsPt * ExpIPdPds * SinTdsHalf) +
                               CosTdsHalf * Exp2IPdsPt * ExpIPdsPd * SinTdHalfSq * SinTdsHalfCu;

        t *= prefactor * eMuOnSqrt2MGam;
    }


} // namespace Hammer
