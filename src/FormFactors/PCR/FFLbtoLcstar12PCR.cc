///
/// @file  FFLbtoLcstar12PCR.cc
/// @brief \f$ \Lambda_b \rightarrow Lambda_c^*(2595) \f$ PCR form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/FormFactors/PCR/FFLbtoLcstar12PCR.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFLbtoLcstar12PCR::FFLbtoLcstar12PCR() {
        // Create tensor rank and dimensions
        IndexList dims = {12};
        string name{"FFLbtoLcstar12PCR"};

        setPrefix("LbtoLc*2595");
        addProcessSignature(-PID::LAMBDAB, {PID::LAMBDACSTAR12MINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_LBLCSTAR12})});

        setSignatureIndex();
    }

    void FFLbtoLcstar12PCR::defineSettings() {
        //_mFFErrNames = ;
        setPath(getFFErrPrefixGroup().get());
        setUnits("GeV");

        addSetting<double>("mQ", 5.28);
        addSetting<double>("mq", 1.89);
        addSetting<double>("md", 0.40);
        addSetting<double>("aL", 0.59);
        addSetting<double>("aLp", 0.47);
        addSetting<double>("k2", 1.0);

        setInitialized();
    }

    void FFLbtoLcstar12PCR::eval(const Particle& parent, const ParticleList& daughters,
                                 const ParticleList& /*references*/) {
        // Momenta
        const FourMomentum& pLbmes = parent.momentum();
        const FourMomentum& pLcmes = daughters[0].momentum();
        // const FourMomentum& pTau = daughters[2].momentum();

        // kinematic objects
        const double Mb = pLbmes.mass();
        const double Mc = pLcmes.mass();
        const double Sqq = Mb * Mb + Mc * Mc - 2. * (pLbmes * pLcmes);

        evalAtPSPoint({Sqq}, {Mb, Mc});
    }

    void FFLbtoLcstar12PCR::evalAtPSPoint(const vector<double>& point, const vector<double>& masses) {
        Tensor& result = getTensor();
        result.clearData();

        if (!isInitialized()) {
            MSG_WARNING("Warning, Settings have not been defined!");
        }

        double Mb = 0.;
        double Mc = 0.;
        double unitres = 1.;
        tie(Mb, Mc, unitres) = getParentDaughterHadMasses(masses);


        double Sqq = point[0];

        // const double sqSqq = sqrt(Sqq);
        double w = getW(Sqq, Mb, Mc);
        // safety measure if w==1.0
        if (isZero(w - 1.0)) {
            w += 1e-6;
        }

        // Parameters
        const double mQ = (*getSetting<double>("mQ"));
        const double mq = (*getSetting<double>("mq"));
        const double md = (*getSetting<double>("md"));
        const double aL = (*getSetting<double>("aL"));
        const double aLp = (*getSetting<double>("aLp"));
        const double k2 = (*getSetting<double>("k2"));

        const double aL2 = aL * aL;
        const double aLp2 = aLp * aLp;
        const double aLLp2 = 0.5 * (aL2 + aLp2);
        const double rho2 = 3. * md * md / (2. * k2 * aLLp2);

        const double IH = pow(aL * aLp / aLLp2, 2.5) * exp(-rho2 * (w * w - 1.));

        // set elements F_123 and G_123
        result.element({2}) = aL / 6. * (3. / mq - 1. / mQ);
        result.element({3}) = -(2. * md / aL - aL / (2. * mq) + 2. * md * md * aL / (mQ * aLLp2) -
                                (md * aL / (6. * mq * mQ * aLLp2)) * (3. * aL2 - 2. * aLp2));
        result.element({4}) = 2. * md * md * aL / (mQ * aLLp2);
        result.element({5}) =
            2. * md / aL - aL / (6. * mQ) + (md * aL / (6. * mq * mQ * aLLp2)) * (3. * aL2 - 2. * aLp2);
        result.element({6}) = -2. * md / aL + aL / (2. * mq) + aL / (3. * mQ);
        result.element({7}) = aL / (3. * mQ) * (1. - (md / (2. * mq * aLLp2)) * (3. * aL2 - 2. * aLp2));

        result *= IH;
    }

    std::unique_ptr<FormFactorBase> FFLbtoLcstar12PCR::clone(const std::string& label) {
        MAKE_CLONE(FFLbtoLcstar12PCR, label);
    }

} // namespace Hammer
