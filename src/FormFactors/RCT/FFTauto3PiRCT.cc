///
/// @file  FFTauto3PiRCT.cc
/// @brief \f$ \tau^+ \rightarrow \bar\nu\pi^+\pi^+\pi^- \f$ form factors
/// see 1203.3955 and 1310.1053
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/FormFactors/RCT/FFTauto3PiRCT.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

// From 1310.1053 and 1203.3955

using namespace std;
using namespace complex_literals;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFTauto3PiRCT::FFTauto3PiRCT() {
        // Create tensor rank and dimensions
        IndexList dims = {3};
        string name{"FFTauto3PiRCT"};
        // Ordering is piminus then piplus as anti-tau parent has pdg < 0

        setPrefix("Tauto3Pi");
        addProcessSignature(PID::ANTITAU, {PID::PIMINUS, PID::PIPLUS, PID::PIPLUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_TAU3PI})});

        setSignatureIndex();
    }

    void FFTauto3PiRCT::defineSettings() {
        //_mFFErrNames = ;
        setPath(getFFErrPrefixGroup().get());
        setUnits("GeV");

        // Fit parameters from 1310.1053 in GeV
        addSetting<double>("alphaS", -8.);
        addSetting<double>("betaS", 9.);
        addSetting<double>("gammaS", 1.2);
        addSetting<double>("deltaS", 0.6);
        addSetting<double>("RS", 1.8); // GeV^-1
        addSetting<double>("MRho", 0.7718);
        addSetting<double>("MRhoP", 1.35);
        //        addSetting<double>("GRho",0.149);
        addSetting<double>("GRhoP", 0.44);
        addSetting<double>("MA", 1.091);
        addSetting<double>("GA", 0.5); // PDG
        addSetting<double>("MS", 0.48);
        addSetting<double>("GS", 0.7);
        addSetting<double>("fPi", 0.0913);
        addSetting<double>("FV", 0.1686);
        addSetting<double>("FA", 0.131);
        addSetting<double>("betaRhoP", -0.318);
        addSetting<double>("Vud", sqrt(1 - 0.2257 * 0.2257));

        setInitialized();
    }

    void FFTauto3PiRCT::eval(const Particle& /*parent*/, const ParticleList& daughters,
                             const ParticleList& /*references*/) {
        // Momenta
        //        const FourMomentum& pTau = parent.momentum();
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

        //        const double Mt = pTau.mass();
        //        const double Mp = pPiMinus.mass();
        // Follows 1203.3955 s1 and s2 convention
        const double s1 = (pPiPlus2 + pPiMinus).mass2();
        const double s2 = (pPiPlus1 + pPiMinus).mass2();
        const double s3 = (pPiPlus1 + pPiPlus2).mass2();
        const double Sqp = (pPiPlus1 + pPiPlus2 + pPiMinus).mass2();

        const double Mp = sqrt((s1 + s2 + s3 - Sqp) / 3.);

        evalAtPSPoint({Sqp, s1, s2}, {Mp});
    }

    void FFTauto3PiRCT::evalAtPSPoint(const vector<double>& point, const vector<double>& masses) {
        Tensor& result = getTensor();
        result.clearData();

        if (!isInitialized()) {
            MSG_WARNING("Warning, Settings have not been defined!");
        }

        double Mp = 0.;
        double unitres = 1.;
        if (!masses.empty()) {
            Mp = masses[0];
            unitres = _units;
        } else {
            Mp = this->masses()[1];
        }
        double Sqp = point[0];
        double s1 = point[1];
        double s2 = point[2];
        double s3 = Sqp + 3 * Mp * Mp - s1 - s2;

        // Fit parameters in GeV
        const double alphaS = (*getSetting<double>("alphaS"));
        const double betaS = (*getSetting<double>("betaS"));
        const double gammaS = (*getSetting<double>("gammaS"));
        const double deltaS = (*getSetting<double>("deltaS"));
        const double RS = (*getSetting<double>("RS")) / unitres;
        const double MRho = unitres * (*getSetting<double>("MRho"));
        const double MRhoP = unitres * (*getSetting<double>("MRhoP"));
        const double GRhoP = unitres * (*getSetting<double>("GRhoP"));
        const double MA = unitres * (*getSetting<double>("MA"));
        const double GA = unitres * (*getSetting<double>("GA"));
        const double MS = unitres * (*getSetting<double>("MS"));
        const double GS = unitres * (*getSetting<double>("GS"));
        const double fPi = unitres * (*getSetting<double>("fPi"));
        const double FV = unitres * (*getSetting<double>("FV"));
        const double FA = unitres * (*getSetting<double>("FA"));
        const double betaRhoP = (*getSetting<double>("betaRhoP"));
        const double Vud = (*getSetting<double>("Vud"));

        const double Mp2 = Mp * Mp;
        const double fPi2 = fPi * fPi;
        const double MRho2 = MRho * MRho;
        const double MRhoP2 = MRhoP * MRhoP;
        const double MA2 = MA * MA;
        const double MS2 = MS * MS;
        const double GV = fPi2 / FV;

        // Widths
        const double sigmapi1 = sqrt(1 - 4 * Mp2 / s1);
        const double sigmapi2 = sqrt(1 - 4 * Mp2 / s2);
        const double sigmapi3 = sqrt(1 - 4 * Mp2 / s3);
        const double MKaon = unitres * 0.494;
        const double MKaon2 = MKaon * MKaon;
        // Following 0911.4436 and tauola code in f3pi rcht.f
        const double GRho1 =
            MRho * s1 / (96 * pi * fPi * fPi) *
            (pow(sigmapi1, 3.) + 0.5 * (1 - 4 * MKaon2 / s1 > 0. ? pow(1 - 4 * MKaon2 / s1, 1.5) : 0.));
        const double GRho2 =
            MRho * s2 / (96 * pi * fPi * fPi) *
            (pow(sigmapi2, 3.) + 0.5 * (1 - 4 * MKaon2 / s2 > 0. ? pow(1 - 4 * MKaon2 / s2, 1.5) : 0.));
        const double GRhoP1 = GRhoP * MRhoP / sqrt(s1) * pow((s1 - 4 * Mp2) / (MRhoP2 - 4 * Mp2), 1.5);
        const double GRhoP2 = GRhoP * MRhoP / sqrt(s2) * pow((s2 - 4 * Mp2) / (MRhoP2 - 4 * Mp2), 1.5);

        const double GS1 = GS * sigmapi1 / sqrt(1 - 4 * Mp2 / MS2);
        const double GS2 = GS * sigmapi2 / sqrt(1 - 4 * Mp2 / MS2);
        const double GAq = GA; // Placeholder. 1310.1053 uses a linear interpolation in s

        // Constituent BWs. Signs of imGamma terms follow tauola code funct_3pi.f and funct_rpt.f, different to
        // 1203.3955, 1310.1053
        const complex<double> BWRho1 = 1. / (s1 - MRho2 + 1i * MRho * GRho1);
        const complex<double> BWRho2 = 1. / (s2 - MRho2 + 1i * MRho * GRho2);
        const complex<double> BWRhoP1 = 1. / (s1 - MRhoP2 + 1i * MRhoP * GRhoP1);
        const complex<double> BWRhoP2 = 1. / (s2 - MRhoP2 + 1i * MRhoP * GRhoP2);
        const complex<double> BWRhoC1 = (BWRho1 + betaRhoP * BWRhoP1) / (1 + betaRhoP);
        const complex<double> BWRhoC2 = (BWRho2 + betaRhoP * BWRhoP2) / (1 + betaRhoP);
        const complex<double> BWS1 = 1. / (MS2 - s1 - 1i * MS * GS1);
        const complex<double> BWS2 = 1. / (MS2 - s2 - 1i * MS * GS2);
        const complex<double> BWA = 1. / (Sqp - MA2 + 1i * MA * GAq);

        // Constituent sigma exponentials
        const double Lambda1 = pow(Sqp - s1 - Mp2, 2.) - 4 * s1 * Mp2;
        const double Lambda2 = pow(Sqp - s2 - Mp2, 2.) - 4 * s2 * Mp2;
        const double FS1 = exp(-Lambda1 * pow(RS, 2.) / (8 * Sqp));
        const double FS2 = exp(-Lambda2 * pow(RS, 2.) / (8 * Sqp));

        // Lambdas & H fn
        const double Lp = fPi2 / (2. * sqrt2 * FA * GV);
        const double Lpp = -(1. - 2. * pow(GV, 2.) / fPi2) * Lp;
        const double Lz = (Lp + Lpp) / 4.;
        const double H1 = -Lz * Mp2 / Sqp + Lp * s1 / Sqp + Lpp;
        const double H2 = -Lz * Mp2 / Sqp + Lp * s2 / Sqp + Lpp;

        // Pieces of structure functions. F1(q^2, s1, s2) <-> F2(q^2, s1, s2) wrt 1203.3955 (Using Z Phys C 56 661
        // (1992) convention for F1 and F2)
        const double F1X = -2. * sqrt2 / 3.;
        const double F2X = -2. * sqrt2 / 3.;

        const complex<double> F1R =
            (3 * s2 * BWRhoC2 - (2 * GV / FV - 1.) * ((2 * Sqp - 2 * s2 - s3) * BWRhoC2 + (s3 - s2) * BWRhoC1) +
             (alphaS * FS2 * BWS2 + betaS * FS1 * BWS1) * MS2) *
            (sqrt2 * FV * GV / (3 * fPi2));
        const complex<double> F2R =
            (3 * s1 * BWRhoC1 - (2 * GV / FV - 1.) * ((2 * Sqp - 2 * s1 - s3) * BWRhoC1 + (s3 - s1) * BWRhoC2) +
             (alphaS * FS1 * BWS1 + betaS * FS2 * BWS2) * MS2) *
            (sqrt2 * FV * GV / (3 * fPi2));

        const complex<double> F1RR = (-(Lp + Lpp) * 3 * s2 * BWRhoC2 + H2 * (2 * Sqp + s2 - s3) * BWRhoC2 +
                                      H1 * (s3 - s2) * BWRhoC1 + (gammaS * BWS2 * FS2 + deltaS * BWS1 * FS1) * MS2) *
                                     (4 * FA * GV / (3 * fPi2) * Sqp * BWA);
        const complex<double> F2RR = (-(Lp + Lpp) * 3 * s1 * BWRhoC1 + H1 * (2 * Sqp + s1 - s3) * BWRhoC1 +
                                      H2 * (s3 - s1) * BWRhoC2 + (gammaS * BWS1 * FS1 + deltaS * BWS2 * FS2) * MS2) *
                                     (4 * FA * GV / (3 * fPi2) * Sqp * BWA);

        const double F4X = 2. * sqrt2 / 3 * Mp2 / (2 * Sqp * (Sqp - Mp2)) * (3 * (s3 - Mp2) - 3 * Sqp);
        const complex<double> F4R = -sqrt2 * pow(GV, 2.) / fPi2 * Mp2 / (Sqp - Mp2) *
                                    (s1 / Sqp * (s3 - s2) * BWRhoC1 + s2 / Sqp * (s3 - s1) * BWRhoC2);

        // Coulomb correction
        const double alphaQED = 1. / 137.036;
        const double v01 = 2 * sigmapi1 / (1 + sigmapi1 * sigmapi1);
        const double v02 = 2 * sigmapi2 / (1 + sigmapi2 * sigmapi2);
        const double v03 = 2 * sigmapi3 / (1 + sigmapi3 * sigmapi3);

        const double cf1 = 2 * pi * alphaQED / v01 / (1. - exp(-2 * pi * alphaQED / v01));
        const double cf2 = 2 * pi * alphaQED / v02 / (1. - exp(-2 * pi * alphaQED / v02));
        const double cf3 = 2 * pi * alphaQED / v03 / (exp(2 * pi * alphaQED / v03) - 1.);

        // Set tensor elements
        result.element({0}) = F1X + F1R + F1RR;
        result.element({1}) = F2X + F2R + F2RR;
        result.element({2}) = F4X + F4R;

        result *= (Vud / fPi) * sqrt(cf1 * cf2 * cf3);
    }

    unique_ptr<FormFactorBase> FFTauto3PiRCT::clone(const string& label) {
        MAKE_CLONE(FFTauto3PiRCT, label);
    }

} // namespace Hammer
