///
/// @file  FFD1DstarPiPW.cc
/// @brief \f$ D_1 \rightarrow D^* \pi \f$ partial wave coefficients
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/FormFactors/PW/FFD1DstarPiPW.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFD1DstarPiPW::FFD1DstarPiPW() {
        // Create tensor rank and dimensions
        IndexList dims = {2};
        string name{"FFD1toDstarPiPW"};

        setPrefix("D**1toD*Pi");
        addProcessSignature(-PID::DSSD1, {PID::DSTARMINUS, PID::PIPLUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD1})});

        addProcessSignature(-PID::DSSD1, {-PID::DSTAR, PID::PI0});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD1})});

        addProcessSignature(PID::DSSD1MINUS, {-PID::DSTAR, PID::PIMINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD1})});

        addProcessSignature(PID::DSSD1MINUS, {PID::DSTARMINUS, PID::PI0});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD1})});

        // strange
        setPrefix("Ds**1toDs*Pi");
        addProcessSignature(PID::DSSDS1MINUS, {PID::DSSTARMINUS, PID::PI0});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSDS1})});

        setPrefix("Ds**1toD*K");
        addProcessSignature(PID::DSSDS1MINUS, {PID::DSTARMINUS, -PID::K0});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSDS1})});

        addProcessSignature(PID::DSSDS1MINUS, {-PID::DSTAR, PID::KMINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSDS1})});

        setSignatureIndex();
    }

    void FFD1DstarPiPW::defineSettings() {
        //_mFFErrNames = ;
        setPath(getFFErrPrefixGroup().get());
        setUnits("GeV");

        addSetting<complex<double>>("S", 0.0); // GeV^2
        addSetting<complex<double>>("D", 1.0);

        setInitialized();
    }

    void FFD1DstarPiPW::eval(const Particle& parent, const ParticleList& daughters,
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

    void FFD1DstarPiPW::evalAtPSPoint(const vector<double>& /*point*/, const vector<double>& /*masses*/) {
        Tensor& result = getTensor();
        result.clearData();

        if (!isInitialized()) {
            MSG_WARNING("Warning, Settings have not been defined!");
        }


        // Conjugation to match phase convention in EvtGen
        const complex<double> S = conj(*getSetting<complex<double>>("S"));
        const complex<double> D = conj(*getSetting<complex<double>>("D"));

        // set elements S, D
        result.element({0}) = S;
        result.element({1}) = D;
    }

    std::unique_ptr<FormFactorBase> FFD1DstarPiPW::clone(const std::string& label) {
        MAKE_CLONE(FFD1DstarPiPW, label);
    }

} // namespace Hammer
