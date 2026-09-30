///
/// @file  Particle.cc
/// @brief Hammer particle class
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/Particle.hh"
#include "Hammer/Tools/HammerRoot.hh"


using namespace std;

namespace Hammer {

    Particle::Particle() : _code{0}, _momentum(0., 0., 0., 0.) {
    }

    Particle::Particle(const FourMomentum& p, PdgId code) : _code{code}, _momentum{p} {
    }

#ifdef HAVE_ROOT
    Particle Particle::fromRoot(const TLorentzVector& p, PdgId code) {
        return Particle{{p.E(), p.Px(), p.Py(), p.Pz()}, code};
    }

#endif

    Particle& Particle::setMomentum(const FourMomentum& p) {
        _momentum = p;
        return *this;
    }

    Particle& Particle::setPdgId(PdgId code) {
        _code = code;
        return *this;
    }

    PdgId Particle::pdgId() const {
        return _code;
    }

    const FourMomentum& Particle::momentum() const {
        return _momentum;
    }

    FourMomentum& Particle::momentum() {
        return _momentum;
    }

    const FourMomentum& Particle::p() const {
        return momentum();
    }

    FourMomentum& Particle::p() {
        return _momentum;
    }

} // namespace Hammer
