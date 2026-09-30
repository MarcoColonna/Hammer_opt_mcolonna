///
/// @file  FFBctoJpsiKiselev.cc
/// @brief \f$ B_c \rightarrow J/\psi \f$ PCR form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/FormFactors/Kiselev/FFBctoJpsiKiselev.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFBctoJpsiKiselev::FFBctoJpsiKiselev() {
        // Create tensor rank and dimensions
        IndexList dims = {8};
        string name{"FFBctoJpsiKiselev"};

        setPrefix("BctoJpsi");
        addProcessSignature(PID::BCPLUS, {PID::JPSI});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BCJPSI})});

        setSignatureIndex();
    }

    void FFBctoJpsiKiselev::defineSettings() {
        //_mFFErrNames = ;
        setPath(getFFErrPrefixGroup().get());
        setUnits("GeV");

        setInitialized();
    }

    void FFBctoJpsiKiselev::eval(const Particle& parent, const ParticleList& daughters,
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

    void FFBctoJpsiKiselev::evalAtPSPoint(const vector<double>& point, const vector<double>& masses) {
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
        const double t = Sqq / (unitres * unitres);    // GeV^2 units
        const double sumMasses = (Mb + Mc) / unitres;  // GeV units
        const double diffMasses = (Mb - Mc) / unitres; // GeV units

        // Ported from EvtGen BCVFF
        double Mpole2 = 4.5 * 4.5;
        double Den = 1. / (1. - t / Mpole2);
        double FV = 0.11 * Den;
        double FAp = -0.074 * Den;
        double FA0 = 5.9 * Den;
        double FAm = 0.12 * Den;

        const double Vf = sumMasses * FV;
        const double A2f = -sumMasses * FAp;
        const double A1f = FA0 / sumMasses;
        const double A0f = (t * FAm + sumMasses * A1f - diffMasses * A2f) / (sumMasses - diffMasses);

        // Ff
        result.element({1}) = A1f * Mb * (1. + rC);
        // Fg
        result.element({2}) = Vf / (Mb * (1. + rC));
        // Fm
        result.element({3}) = (A2f * (1. - rC) + 2. * A0f * rC - A1f * (1. + rC)) * Mb / Sqq;
        // Fp
        result.element({4}) = -A2f / (Mb * (1. + rC));
    }

    std::unique_ptr<FormFactorBase> FFBctoJpsiKiselev::clone(const std::string& label) {
        MAKE_CLONE(FFBctoJpsiKiselev, label);
    }

} // namespace Hammer
