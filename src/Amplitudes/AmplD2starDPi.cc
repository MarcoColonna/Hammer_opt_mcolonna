///
/// @file  AmplD2starDPi.cc
/// @brief \f$ D_2^* \rightarrow D \pi \f$ amplitude
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/Amplitudes/AmplD2starDPi.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;
using namespace complex_literals;

namespace Hammer {

    namespace MD = MultiDimensional;

    AmplD2starDPi::AmplD2starDPi() {
        // Create tensor rank and dimensions
        IndexList dims = {1, 5};
        string name{"AmplD2starDPi"};

        addProcessSignature(-PID::DSSD2STAR, {PID::DMINUS, PID::PIPLUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD2STAR, SPIN_DSSD2STAR})});

        addProcessSignature(-PID::DSSD2STAR, {-PID::D0, PID::PI0});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD2STAR, SPIN_DSSD2STAR})});

        addProcessSignature(PID::DSSD2STARMINUS, {-PID::D0, PID::PIMINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD2STAR, SPIN_DSSD2STAR})});

        addProcessSignature(PID::DSSD2STARMINUS, {PID::DMINUS, PID::PI0});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD2STAR, SPIN_DSSD2STAR})});

        // strange
        addProcessSignature(PID::DSSDS2STARMINUS, {PID::DSMINUS, PID::PI0});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSDS2STAR, SPIN_DSSDS2STAR})});

        addProcessSignature(PID::DSSDS2STARMINUS, {PID::DMINUS, -PID::K0});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSDS2STAR, SPIN_DSSDS2STAR})});

        addProcessSignature(PID::DSSDS2STARMINUS, {-PID::D0, PID::KMINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSDS2STAR, SPIN_DSSDS2STAR})});

        setSignatureIndex();
        _multiplicity = 5ul;
    }

    void AmplD2starDPi::defineSettings() {
    }

    void AmplD2starDPi::eval(const Particle& parent, const ParticleList& daughters, const ParticleList& references) {
        // Momenta
        const FourMomentum& pDssmes = parent.momentum();
        const FourMomentum& pDmes = daughters[0].momentum();
        const FourMomentum& pPmes = daughters[1].momentum();

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
        const double Md = pDmes.mass();
        const double Md2 = Md * Md;
        const double Mp = pPmes.mass();
        const double Mp2 = Mp * Mp;

        const double Eds = (Mdss2 + Md2 - Mp2) / (2 * Mdss);
        const double Ep = (Mdss2 - Md2 + Mp2) / (2 * Mdss);
        const double Pp = sqrt(Ep * Ep - Mp2);
        const double Pp2 = Pp * Pp;

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

        const double CosTds = (Eds * Edss * Mb - Mdss * (pBmes * pDmes)) / (Mb * Pp * Pw);
        const double SinTds = sqrt(1. - CosTds * CosTds);
        const double CosTds2 = CosTds * CosTds;
        const double SinTds2 = SinTds * SinTds;

        const double CosPdsPt =
            (Ep * Sqq * (pDmes * kNuTau) - Eds * Sqq * (pPmes * kNuTau) +
             NuTauQ * (Pp * (-Ew * Mb + Sqq) * CosTds * CosTt - Ep * (pDmes * Qmes) + Eds * (pPmes * Qmes))) /
            (SinTds * SinTt * Mdss * Pp * sqSqq * NuTauQ);
        const double SinPdsPt =
            -sqSqq * epsilon(pBmes, pDmes, pDssmes, kNuTau) / (SinTds * SinTt * Mb * Pp * Pw * NuTauQ);
        const complex<double> ExpIPdsPt = CosPdsPt + 1i * SinPdsPt;
        const complex<double> ExpIPtPds = CosPdsPt - 1i * SinPdsPt;
        const complex<double> Exp2IPdsPt = ExpIPdsPt * ExpIPdsPt;
        const complex<double> Exp2IPtPds = ExpIPtPds * ExpIPtPds;

        const double prefactor = Pp2;

        // initialize tensor elements to zero
        Tensor& t = getTensor();
        t.clearData();

        // set non-zero tensor elements
        t.element({0, 0}) = (Exp2IPtPds * SinTds2) / 2.;
        t.element({0, 1}) = -(CosTds * ExpIPtPds * SinTds);
        t.element({0, 2}) = (1 - 3 * CosTds2) / sqrt(6);
        t.element({0, 3}) = -(CosTds * ExpIPdsPt * SinTds);
        t.element({0, 4}) = (Exp2IPdsPt * SinTds2) / 2.;

        t *= prefactor;
    }


} // namespace Hammer
