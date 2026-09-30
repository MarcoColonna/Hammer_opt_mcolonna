///
/// @file  FFLbtoLcPCR.cc
/// @brief \f$ \Lambda_b \rightarrow Lambda_c \f$ PCR form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/FormFactors/PCR/FFLbtoLcPCR.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFLbtoLcPCR::FFLbtoLcPCR() {
        // Create tensor rank and dimensions
        IndexList dims = {12};
        string name{"FFLbtoLcPCR"};

        setPrefix("LbtoLc");
        addProcessSignature(-PID::LAMBDAB, {PID::LAMBDACMINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_LBLC})});

        setSignatureIndex();
    }

    void FFLbtoLcPCR::defineSettings() {
        //_mFFErrNames = ;
        setPath(getFFErrPrefixGroup().get());
        setUnits("GeV");

        addSetting<double>("mQ", 5.28);
        addSetting<double>("mq", 1.89);
        addSetting<double>("md", 0.40);
        addSetting<double>("aL", 0.59);
        addSetting<double>("aLp", 0.55);
        addSetting<double>("k2", 1.0);

        setInitialized();
    }

    void FFLbtoLcPCR::eval(const Particle& parent, const ParticleList& daughters, const ParticleList& /*references*/) {
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

    void FFLbtoLcPCR::evalAtPSPoint(const vector<double>& point, const vector<double>& masses) {
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

        const double IH = pow(aL * aLp / aLLp2, 1.5) * exp(-rho2 * (w * w - 1.));

        // set elements f_123 and g_123
        result.element({2}) = 1. + (md / aLLp2) * ((aLp2 / mq) + (aL2 / mQ));
        result.element({3}) = -((md / mq) * (aLp2 / aLLp2) - aL2 * aLp2 / (4. * aLLp2 * mq * mQ));
        result.element({4}) = -md * aL2 / (mQ * aLLp2);
        result.element({5}) = 1. - (aL2 * aLp2) / (12. * aLLp2 * mq * mQ);
        result.element({6}) =
            -(md * aLp2 / (mq * aLLp2) + (aL2 * aLp2) / (12. * aLLp2 * mq * mQ) * (1. + 12. * md * md / aLLp2));
        result.element({7}) = md * aL2 / (mQ * aLLp2) + md * md * aL2 * aLp2 / (mq * mQ * aLLp2 * aLLp2);

        result *= IH;
    }

    std::unique_ptr<FormFactorBase> FFLbtoLcPCR::clone(const std::string& label) {
        MAKE_CLONE(FFLbtoLcPCR, label);
    }

} // namespace Hammer
