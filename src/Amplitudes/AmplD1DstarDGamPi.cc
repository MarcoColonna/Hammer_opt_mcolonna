///
/// @file  AmplD1DstarDGamPi.cc
/// @brief \f$ D_1 \rightarrow D^* \pi, D^* \rightarrow D \gamma \f$ amplitude
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/Amplitudes/AmplD1DstarDGamPi.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;
using namespace complex_literals;

namespace Hammer {

    namespace MD = MultiDimensional;

    AmplD1DstarDGamPi::AmplD1DstarDGamPi() {
        // Create tensor rank and dimensions
        IndexList dims = {2, 3, 2};
        string name{"AmplD1DstarDGamPi"};

        addProcessSignature(-PID::DSSD1, {PID::DSTARMINUS, PID::PIPLUS}, {PID::DMINUS, PID::GAMMA});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD1, SPIN_DSSD1, SPIN_GAMMA})});

        addProcessSignature(-PID::DSSD1, {-PID::DSTAR, PID::PI0}, {-PID::D0, PID::GAMMA});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD1, SPIN_DSSD1, SPIN_GAMMA})});

        addProcessSignature(PID::DSSD1MINUS, {-PID::DSTAR, PID::PIMINUS}, {-PID::D0, PID::GAMMA});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD1, SPIN_DSSD1, SPIN_GAMMA})});

        addProcessSignature(PID::DSSD1MINUS, {PID::DSTARMINUS, PID::PI0}, {PID::DMINUS, PID::GAMMA});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD1, SPIN_DSSD1, SPIN_GAMMA})});

        // strange
        addProcessSignature(PID::DSSDS1MINUS, {PID::DSSTARMINUS, PID::PI0}, {PID::DSMINUS, PID::GAMMA});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSDS1, SPIN_DSSDS1, SPIN_GAMMA})});

        addProcessSignature(PID::DSSDS1MINUS, {PID::DSTARMINUS, -PID::K0}, {PID::DMINUS, PID::GAMMA});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSDS1, SPIN_DSSDS1, SPIN_GAMMA})});

        addProcessSignature(PID::DSSDS1MINUS, {-PID::DSTAR, PID::KMINUS}, {-PID::D0, PID::GAMMA});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSDS1, SPIN_DSSDS1, SPIN_GAMMA})});

        setSignatureIndex();
        _multiplicity = 3ul;
    }

    void AmplD1DstarDGamPi::defineSettings() {
    }

    void AmplD1DstarDGamPi::eval(const Particle& parent, const ParticleList& daughters,
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

        const double CosTd = (-Ed * Eds * Mdss + Mds * (pDssmes * pDmes)) / (Mdss * Pp * Eg);
        const double SinTd = sqrt(1. - CosTd * CosTd);

        const double CosPdsPt =
            (Ep * Sqq * (pDsmes * kNuTau) - Eds * Sqq * (pPmes * kNuTau) +
             NuTauQ * (Pp * (-Ew * Mb + Sqq) * CosTds * CosTt - Ep * (pDsmes * Qmes) + Eds * (pPmes * Qmes))) /
            (SinTds * SinTt * Mdss * Pp * sqSqq * NuTauQ);
        const double SinPdsPt =
            -sqSqq * epsilon(pBmes, pDsmes, pDssmes, kNuTau) / (SinTds * SinTt * Mb * Pp * Pw * NuTauQ);
        const complex<double> ExpIPdsPt = CosPdsPt + 1i * SinPdsPt;
        const complex<double> ExpIPtPds = CosPdsPt - 1i * SinPdsPt;


        const double CosPdPds = (Eg * Mdss2 * (pBmes * pDmes) - Ed * Mdss2 * (pBmes * kGmes) +
                                 Mb * (Eds * Mdss * Eg * Pw * CosTd * CosTds - Edss * Eg * (pDssmes * pDmes) +
                                       Ed * Edss * (pDssmes * kGmes))) /
                                (SinTds * SinTd * Mb * Mds * Mdss * Eg * Pw);
        const double SinPdPds = -epsilon(pBmes, pDsmes, pDssmes, pDmes) / (SinTds * SinTd * Mb * Pp * Eg * Pw);
        const complex<double> ExpIPdPds = CosPdPds + 1i * SinPdPds;
        const complex<double> ExpIPdsPd = CosPdPds - 1i * SinPdPds;


        const complex<double> prefactor = 2i * Mds * Pp2 * Eg / sqrt3;
        const double eMuOnSqrt2MGam = sqrt(48. * pi * Mds2 / pow(2 * Mds * Eg, 3.)); // for D* decay

        // initialize tensor elements to zero
        Tensor& t = getTensor();
        t.clearData();

        // set non-zero tensor elements
        t.element({0, 0, 0}) =
            (3. * ExpIPtPds *
             ((1 + CosTd) * CosTdsHalfSq * ExpIPdsPd - SinTd * SinTds - (-1 + CosTd) * ExpIPdPds * SinTdsHalfSq)) /
            (4. * sqrt2);
        t.element({0, 0, 1}) = (3. * ExpIPtPds *
                                (CosTdsHalfSq * (ExpIPdsPd - CosTd * ExpIPdsPd) + SinTd * SinTds +
                                 (1 + CosTd) * ExpIPdPds * SinTdsHalfSq)) /
                               (4. * sqrt2);
        t.element({0, 1, 0}) =
            (3. * (2. * CosTds * SinTd + ((-1 + CosTd) * ExpIPdPds + (1 + CosTd) * ExpIPdsPd) * SinTds)) / 8.;
        t.element({0, 1, 1}) =
            (-3. * (2. * CosTds * SinTd + ((1 + CosTd) * ExpIPdPds + (-1 + CosTd) * ExpIPdsPd) * SinTds)) / 8.;
        t.element({0, 2, 0}) =
            (3. * ExpIPdsPt *
             ((-1 + CosTd) * CosTdsHalfSq * ExpIPdPds - SinTd * SinTds - (1 + CosTd) * ExpIPdsPd * SinTdsHalfSq)) /
            (4. * sqrt2);
        t.element({0, 2, 1}) =
            (-3. * ExpIPdsPt *
             ((1 + CosTd) * CosTdsHalfSq * ExpIPdPds - SinTd * SinTds - (-1 + CosTd) * ExpIPdsPd * SinTdsHalfSq)) /
            (4. * sqrt2);
        t.element({1, 0, 0}) = (ExpIPtPds * ((1 + CosTd) * CosTdsHalfSq * ExpIPdsPd + 2. * SinTd * SinTds -
                                             (-1 + CosTd) * ExpIPdPds * SinTdsHalfSq)) /
                               (4. * sqrt2);
        t.element({1, 0, 1}) = (ExpIPtPds * (CosTdsHalfSq * (ExpIPdsPd - CosTd * ExpIPdsPd) - 2. * SinTd * SinTds +
                                             (1 + CosTd) * ExpIPdPds * SinTdsHalfSq)) /
                               (4. * sqrt2);
        t.element({1, 1, 0}) =
            (-4. * CosTds * SinTd + ((-1 + CosTd) * ExpIPdPds + (1 + CosTd) * ExpIPdsPd) * SinTds) / 8.;
        t.element({1, 1, 1}) =
            (4. * CosTds * SinTd - ((1 + CosTd) * ExpIPdPds + (-1 + CosTd) * ExpIPdsPd) * SinTds) / 8.;
        t.element({1, 2, 0}) = (ExpIPdsPt * ((-1 + CosTd) * CosTdsHalfSq * ExpIPdPds + 2. * SinTd * SinTds -
                                             (1 + CosTd) * ExpIPdsPd * SinTdsHalfSq)) /
                               (4. * sqrt2);
        t.element({1, 2, 1}) = -(ExpIPdsPt * ((1 + CosTd) * CosTdsHalfSq * ExpIPdPds + 2. * SinTd * SinTds -
                                              (-1 + CosTd) * ExpIPdsPd * SinTdsHalfSq)) /
                               (4. * sqrt2);


        t *= prefactor * eMuOnSqrt2MGam;
    }


} // namespace Hammer
