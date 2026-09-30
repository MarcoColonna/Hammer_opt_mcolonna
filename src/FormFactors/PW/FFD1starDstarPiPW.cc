///
/// @file  FFD1starDstarPiPW.cc
/// @brief \f$ D_1^* \rightarrow D^* \pi \f$ partial wave coefficients
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/FormFactors/PW/FFD1starDstarPiPW.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFD1starDstarPiPW::FFD1starDstarPiPW() {
        // Create tensor rank and dimensions
        IndexList dims = {1};
        string name{"FFD1startoDstarPiPW"};

        setPrefix("D**1*toD*Pi");
        addProcessSignature(-PID::DSSD1STAR, {PID::DSTARMINUS, PID::PIPLUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD1STAR})});

        addProcessSignature(-PID::DSSD1STAR, {-PID::DSTAR, PID::PI0});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD1STAR})});

        addProcessSignature(PID::DSSD1STARMINUS, {-PID::DSTAR, PID::PIMINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD1STAR})});

        addProcessSignature(PID::DSSD1STARMINUS, {PID::DSTARMINUS, PID::PI0});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD1STAR})});

        // strange
        setPrefix("Ds**1*toDs*Pi");
        addProcessSignature(PID::DSSDS1STARMINUS, {PID::DSSTARMINUS, PID::PI0});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSDS1STAR})});

        setPrefix("Ds**1*toD*K");
        addProcessSignature(PID::DSSDS1STARMINUS, {PID::DSTARMINUS, -PID::K0});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSDS1STAR})});

        addProcessSignature(PID::DSSDS1STARMINUS, {-PID::DSTAR, PID::KMINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSDS1STAR})});

        setSignatureIndex();
    }

    void FFD1starDstarPiPW::defineSettings() {
        //_mFFErrNames = ;
        setPath(getFFErrPrefixGroup().get());
        setUnits("GeV");

        addSetting<complex<double>>("S", 1.0);

        setInitialized();
    }

    void FFD1starDstarPiPW::eval(const Particle& parent, const ParticleList& daughters,
                                 const ParticleList& /*references*/) {
        // Momenta
        const FourMomentum& pD1 = parent.momentum();
        const FourMomentum& pDs = daughters[0].momentum();
        const FourMomentum& pPi = daughters[1].momentum();

        // kinematic objects
        const double Mdss = pD1.mass();
        const double Mds = pDs.mass();
        const double Mp = pPi.mass();
        // const double Sqq = Mdss*Mdss + Mds*Mds - 2. * (pD1 * pDs);

        evalAtPSPoint({}, {Mdss, Mds, Mp});
    }

    void FFD1starDstarPiPW::evalAtPSPoint(const vector<double>& /*point*/, const vector<double>& /*masses*/) {
        Tensor& result = getTensor();
        result.clearData();

        if (!isInitialized()) {
            MSG_WARNING("Warning, Settings have not been defined!");
        }


        // Conjugation to match phase convention in EvtGen
        const complex<double> S = conj(*getSetting<complex<double>>("S"));

        // set elements S
        result.element({0}) = S;
    }

    std::unique_ptr<FormFactorBase> FFD1starDstarPiPW::clone(const std::string& label) {
        MAKE_CLONE(FFD1starDstarPiPW, label);
    }

} // namespace Hammer
