///
/// @file  AmplD1DstarDPiPi.cc
/// @brief \f$ D_1 \rightarrow D^* \pi, D^* \rightarrow D \pi \f$ amplitude
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/Amplitudes/AmplD1DstarDPiPi.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;
using namespace complex_literals;

namespace Hammer {

    namespace MD = MultiDimensional;

    AmplD1DstarDPiPi::AmplD1DstarDPiPi() {
        // Create tensor rank and dimensions
        IndexList dims = {2, 3};
        string name{"AmplD1DstarDPiPi"};

        addProcessSignature(-PID::DSSD1, {PID::DSTARMINUS, PID::PIPLUS}, {PID::DMINUS, PID::PI0});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD1, SPIN_DSSD1})});

        addProcessSignature(-PID::DSSD1, {PID::DSTARMINUS, PID::PIPLUS}, {-PID::D0, PID::PIMINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD1, SPIN_DSSD1})});

        addProcessSignature(-PID::DSSD1, {-PID::DSTAR, PID::PI0}, {-PID::D0, PID::PI0});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD1, SPIN_DSSD1})});

        addProcessSignature(PID::DSSD1MINUS, {-PID::DSTAR, PID::PIMINUS}, {-PID::D0, PID::PI0});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD1, SPIN_DSSD1})});

        addProcessSignature(PID::DSSD1MINUS, {PID::DSTARMINUS, PID::PI0}, {PID::DMINUS, PID::PI0});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD1, SPIN_DSSD1})});

        addProcessSignature(PID::DSSD1MINUS, {PID::DSTARMINUS, PID::PI0}, {-PID::D0, PID::PIMINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD1, SPIN_DSSD1})});

        // strange
        addProcessSignature(PID::DSSDS1MINUS, {PID::DSSTARMINUS, PID::PI0}, {PID::DSMINUS, PID::PI0});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSDS1, SPIN_DSSDS1})});

        addProcessSignature(PID::DSSDS1MINUS, {PID::DSTARMINUS, -PID::K0}, {PID::DMINUS, PID::PI0});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSDS1, SPIN_DSSDS1})});

        addProcessSignature(PID::DSSDS1MINUS, {PID::DSTARMINUS, -PID::K0}, {-PID::D0, PID::PIMINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSDS1, SPIN_DSSDS1})});

        addProcessSignature(PID::DSSDS1MINUS, {-PID::DSTAR, PID::KMINUS}, {-PID::D0, PID::PI0});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSDS1, SPIN_DSSDS1})});

        setSignatureIndex();
        _multiplicity = 3ul;
    }

    void AmplD1DstarDPiPi::defineSettings() {
    }

    void AmplD1DstarDPiPi::eval(const Particle& parent, const ParticleList& daughters, const ParticleList& references) {
        // Momenta
        const FourMomentum& pDssmes = parent.momentum();
        const FourMomentum& pDsmes = daughters[0].momentum();
        const FourMomentum& pPmes = daughters[1].momentum();
        const FourMomentum& pDmes = daughters[2].momentum();
        const FourMomentum& pPDmes = daughters[3].momentum();

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
        const double MpD2 = pPDmes.mass2();

        const double Eds = (Mdss2 + Mds2 - Mp2) / (2 * Mdss);
        const double Ep = (Mdss2 - Mds2 + Mp2) / (2 * Mdss);
        const double Pp = sqrt(Ep * Ep - Mp2);
        const double Pp2 = Pp * Pp;

        const double Ed = (Mds2 + Md2 - MpD2) / (2 * Mds);
        const double EpD = (Mds2 - Md2 + MpD2) / (2 * Mds);
        const double PpD = sqrt(Ed * Ed - Md2);

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

        const double CosTd = (-Ed * Eds * Mdss + Mds * (pDssmes * pDmes)) / (Mdss * Pp * PpD);
        const double SinTd = sqrt(1. - CosTd * CosTd);

        const double CosPdsPt =
            (Ep * Sqq * (pDsmes * kNuTau) - Eds * Sqq * (pPmes * kNuTau) +
             NuTauQ * (Pp * (-Ew * Mb + Sqq) * CosTds * CosTt - Ep * (pDsmes * Qmes) + Eds * (pPmes * Qmes))) /
            (SinTds * SinTt * Mdss * Pp * sqSqq * NuTauQ);
        const double SinPdsPt =
            -sqSqq * epsilon(pBmes, pDsmes, pDssmes, kNuTau) / (SinTds * SinTt * Mb * Pp * Pw * NuTauQ);
        const complex<double> ExpIPdsPt = CosPdsPt + 1i * SinPdsPt;
        const complex<double> ExpIPtPds = CosPdsPt - 1i * SinPdsPt;


        const double CosPdPds = (EpD * Mdss2 * (pBmes * pDmes) - Ed * Mdss2 * (pBmes * pPDmes) +
                                 Mb * (Eds * Mdss * PpD * Pw * CosTd * CosTds - Edss * EpD * (pDssmes * pDmes) +
                                       Ed * Edss * (pDssmes * pPDmes))) /
                                (SinTds * SinTd * Mb * Mds * Mdss * PpD * Pw);
        const double SinPdPds = -epsilon(pBmes, pDsmes, pDssmes, pDmes) / (SinTds * SinTd * Mb * Pp * PpD * Pw);
        const complex<double> ExpIPdPds = CosPdPds + 1i * SinPdPds;
        const complex<double> ExpIPdsPd = CosPdPds - 1i * SinPdPds;


        const double prefactor = Pp2 * PpD / sqrt3;
        const double gPiOnSqrt2MGam = sqrt(3. * pi * Mds / pow(PpD, 3.)); // for D* decay

        // initialize tensor elements to zero
        Tensor& t = getTensor();
        t.clearData();

        // set non-zero tensor elements
        t.element({0, 0}) =
            -3. * ExpIPtPds * (CosTdsHalfSq * ExpIPdsPd * SinTd + CosTd * SinTds - ExpIPdPds * SinTd * SinTdsHalfSq);
        t.element({0, 1}) = (6 * CosTd * CosTds - 3. * (ExpIPdPds + ExpIPdsPd) * SinTd * SinTds) / sqrt2;
        t.element({0, 2}) =
            -3. * ExpIPdsPt * (CosTdsHalfSq * ExpIPdPds * SinTd + CosTd * SinTds - ExpIPdsPd * SinTd * SinTdsHalfSq);
        t.element({1, 0}) =
            ExpIPtPds * (-(CosTdsHalfSq * ExpIPdsPd * SinTd) + 2. * CosTd * SinTds + ExpIPdPds * SinTd * SinTdsHalfSq);
        t.element({1, 1}) = (-4. * CosTd * CosTds - (ExpIPdPds + ExpIPdsPd) * SinTd * SinTds) / sqrt2;
        t.element({1, 2}) =
            ExpIPdsPt * (-(CosTdsHalfSq * ExpIPdPds * SinTd) + 2. * CosTd * SinTds + ExpIPdsPd * SinTd * SinTdsHalfSq);

        t *= prefactor * gPiOnSqrt2MGam;
    }


} // namespace Hammer
