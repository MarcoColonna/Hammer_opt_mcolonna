///
/// @file  FFD2starDPiPW.cc
/// @brief \f$ D_2^* \rightarrow D \pi \f$ partial wave coefficients
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/FormFactors/PW/FFD2starDPiPW.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFD2starDPiPW::FFD2starDPiPW() {
        // Create tensor rank and dimensions
        IndexList dims = {1};
        string name{"FFD2startoDPiPW"};

        setPrefix("D**2*toDPi");
        addProcessSignature(-PID::DSSD2STAR, {PID::DMINUS, PID::PIPLUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD2STAR})});

        addProcessSignature(-PID::DSSD2STAR, {-PID::D0, PID::PI0});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD2STAR})});

        addProcessSignature(PID::DSSD2STARMINUS, {-PID::D0, PID::PIMINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD2STAR})});

        addProcessSignature(PID::DSSD2STARMINUS, {PID::DMINUS, PID::PI0});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD2STAR})});

        // strange
        setPrefix("Ds**2*toDsPi");
        addProcessSignature(PID::DSSDS2STARMINUS, {PID::DSMINUS, PID::PI0});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSDS2STAR})});

        setPrefix("Ds**2*toDK");
        addProcessSignature(PID::DSSDS2STARMINUS, {PID::DMINUS, -PID::K0});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSDS2STAR})});

        addProcessSignature(PID::DSSDS2STARMINUS, {-PID::D0, PID::KMINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSDS2STAR})});

        setSignatureIndex();
    }

    void FFD2starDPiPW::defineSettings() {
        //_mFFErrNames = ;
        setPath(getFFErrPrefixGroup().get());
        setUnits("GeV");

        addSetting<complex<double>>("D", 1.0);

        setInitialized();
    }

    void FFD2starDPiPW::eval(const Particle& parent, const ParticleList& daughters,
                             const ParticleList& /*references*/) {
        // Momenta
        const FourMomentum& pD1 = parent.momentum();
        const FourMomentum& pD = daughters[0].momentum();
        const FourMomentum& pPi = daughters[1].momentum();

        // kinematic objects
        const double Mdss = pD1.mass();
        const double Md = pD.mass();
        const double Mp = pPi.mass();
        // const double Sqq = Mdss*Mdss + Mds*Mds - 2. * (pD1 * pDs);

        evalAtPSPoint({}, {Mdss, Md, Mp});
    }

    void FFD2starDPiPW::evalAtPSPoint(const vector<double>& /*point*/, const vector<double>& /*masses*/) {
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

    std::unique_ptr<FormFactorBase> FFD2starDPiPW::clone(const std::string& label) {
        MAKE_CLONE(FFD2starDPiPW, label);
    }

} // namespace Hammer
