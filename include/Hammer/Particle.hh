///
/// @file  Particle.hh
/// @brief Hammer particle class
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_TOOLS_PARTICLE
#define HAMMER_TOOLS_PARTICLE

#include "Hammer/Particle.fhh"
#include "Hammer/Math/FourMomentum.hh"
#include "Hammer/Tools/Pdg.fhh"

#include "Hammer/Tools/HammerRoot.fhh"

namespace Hammer {

    /// @brief Particle class
    ///
    /// Defines a particle, ...
    ///
    /// @ingroup Core
    class Particle {

    public:

        Particle();

        /// @brief
        /// @param[in] p
        /// @param[in] code
        Particle(const FourMomentum& p, PdgId code);

        Particle(const Particle& other) = default;
        Particle& operator=(const Particle& other) = default;

#ifdef HAVE_ROOT
        static Particle fromRoot(const TLorentzVector& p, PdgId code);

        template <class T>
        static Particle fromRoot(const ROOT::Math::LorentzVector<T>& p, PdgId code);

#endif
        ~Particle() = default;

        /// @brief
        /// @param[in] p
        /// @return
        Particle& setMomentum(const FourMomentum& p);

        /// @brief
        /// @param[in] code
        /// @return
        Particle& setPdgId(PdgId code);

        /// @brief
        /// @return
        [[nodiscard]] PdgId pdgId() const;

        /// @brief
        /// @return
        [[nodiscard]] const FourMomentum& momentum() const;

        /// @brief
        /// @return
        FourMomentum& momentum();

        /// @brief
        /// @return
        [[nodiscard]] const FourMomentum& p() const;

        /// @brief
        /// @return
        FourMomentum& p();

    private:

        PdgId _code;
        FourMomentum _momentum;
    };

} // namespace Hammer

#include "Hammer/ParticleDefs.hh"

#endif
