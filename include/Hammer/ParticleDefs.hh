#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
///
/// @file  ParticleDefs.hh
/// @brief Hammer particle class template methods definitions
///
#pragma clang diagnostic pop

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-

namespace Hammer {

#ifdef HAVE_ROOT

    template <class T>
    Particle Particle::fromRoot(const ROOT::Math::LorentzVector<T>& p, PdgId code) {
        return Particle{{p.E(), p.Px(), p.Py(), p.Pz()}, code};
    }

#endif


} // namespace Hammer
