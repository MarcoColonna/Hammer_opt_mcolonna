///
/// @file  FFD2starDstarPiPW.cc
/// @brief \f$ D_2^* \rightarrow D^* \pi \f$ partial wave coefficients
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/FormFactors/PW/FFD2starDstarPiPW.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFD2starDstarPiPW::FFD2starDstarPiPW() {
        // Create tensor rank and dimensions
        IndexList dims = {1};
        string name{"FFD2startoDstarPiPW"};

        setPrefix("D**2*toD*Pi");
        addProcessSignature(-PID::DSSD2STAR, {PID::DSTARMINUS, PID::PIPLUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD2STAR})});

        addProcessSignature(-PID::DSSD2STAR, {-PID::DSTAR, PID::PI0});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD2STAR})});

        addProcessSignature(PID::DSSD2STARMINUS, {-PID::DSTAR, PID::PIMINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD2STAR})});

        addProcessSignature(PID::DSSD2STARMINUS, {PID::DSTARMINUS, PID::PI0});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD2STAR})});

        // strange
        setPrefix("Ds**2*toDs*Pi");
        addProcessSignature(PID::DSSDS2STARMINUS, {PID::DSSTARMINUS, PID::PI0});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSDS2STAR})});

        setPrefix("Ds**2*toD*K");
        addProcessSignature(PID::DSSDS2STARMINUS, {PID::DSTARMINUS, -PID::K0});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSDS2STAR})});

        addProcessSignature(PID::DSSDS2STARMINUS, {-PID::DSTAR, PID::KMINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSDS2STAR})});

        setSignatureIndex();
    }

    void FFD2starDstarPiPW::defineSettings() {
        //_mFFErrNames = ;
        setPath(getFFErrPrefixGroup().get());
        setUnits("GeV");

        addSetting<complex<double>>("D", 1.0);

        setInitialized();
    }

    void FFD2starDstarPiPW::eval(const Particle& parent, const ParticleList& daughters,
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

    void FFD2starDstarPiPW::evalAtPSPoint(const vector<double>& /*point*/, const vector<double>& /*masses*/) {
        Tensor& result = getTensor();
        result.clearData();

        if (!isInitialized()) {
            MSG_WARNING("Warning, Settings have not been defined!");
        }


        // Conjugation to match phase convention in EvtGen
        const complex<double> D = conj(*getSetting<complex<double>>("D"));

        // set elements D
        result.element({0}) = D;
    }

    std::unique_ptr<FormFactorBase> FFD2starDstarPiPW::clone(const std::string& label) {
        MAKE_CLONE(FFD2starDstarPiPW, label);
    }

} // namespace Hammer
