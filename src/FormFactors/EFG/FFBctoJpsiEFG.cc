///
/// @file  FFBctoJpsiEFG.cc
/// @brief \f$ B_c \rightarrow J/\psi \f$ PCR form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/FormFactors/EFG/FFBctoJpsiEFG.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFBctoJpsiEFG::FFBctoJpsiEFG() {
        // Create tensor rank and dimensions
        IndexList dims = {8};
        string name{"FFBctoJpsiEFG"};

        setPrefix("BctoJpsi");
        addProcessSignature(PID::BCPLUS, {PID::JPSI});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BCJPSI})});

        setSignatureIndex();
    }

    void FFBctoJpsiEFG::defineSettings() {
        //_mFFErrNames = ;
        setPath(getFFErrPrefixGroup().get());
        setUnits("GeV");

        setInitialized();
    }

    void FFBctoJpsiEFG::eval(const Particle& parent, const ParticleList& daughters,
                             const ParticleList& /*references*/) {
        // Momenta
        const FourMomentum& pBcmes = parent.momentum();
        const FourMomentum& pJpsimes = daughters[0].momentum();

        // kinematic objects
        const double Mb = pBcmes.mass();
        const double Mc = pJpsimes.mass();
        const double Sqq = Mb * Mb + Mc * Mc - 2. * (pBcmes * pJpsimes);

        evalAtPSPoint({Sqq}, {Mb, Mc});
    }

    void FFBctoJpsiEFG::evalAtPSPoint(const vector<double>& point, const vector<double>& masses) {
        Tensor& result = getTensor();
        result.clearData();

        if (!isInitialized()) {
            MSG_WARNING("Warning, Settings have not been defined!");
        }

        double Mb = 0.;
        double Mc = 0.;
        double unitres = 1.;
        tie(Mb, Mc, unitres) = getParentDaughterHadMasses(masses);

        const double rC = Mc / Mb;

        const double Sqq = point[0];
        const double t = Sqq / (unitres * unitres); // GeV^2 units

        // Ported from EvtGen BCVFF
        const double Vf = (0.49077824756158533 - 0.0012925655191347828 * t) / (1 - 0.06292520325875656 * t);
        const double A0f = (0.4160345034630221 - 0.0024720095310225023 * t) / (1 - 0.061603451915567785 * t);
        const double A1f = (0.4970212860605933 - 0.0067519730024654745 * t) / (1 - 0.050487026667172176 * t);
        const double A2f = (0.7315284919705497 + 0.0014263826220727142 * t - 0.0006946090066269195 * t * t) /
                           (1 - 0.04885587273651653 * t);

        // Ff
        result.element({1}) = A1f * Mb * (1. + rC);
        // Fg
        result.element({2}) = Vf / (Mb * (1. + rC));
        // Fm
        result.element({3}) = (A2f * (1. - rC) + 2. * A0f * rC - A1f * (1. + rC)) * Mb / Sqq;
        // Fp
        result.element({4}) = -A2f / (Mb * (1. + rC));
    }

    std::unique_ptr<FormFactorBase> FFBctoJpsiEFG::clone(const std::string& label) {
        MAKE_CLONE(FFBctoJpsiEFG, label);
    }

} // namespace Hammer
