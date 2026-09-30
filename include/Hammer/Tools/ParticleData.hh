///
/// @file  ParticleData.hh
/// @brief Hammer decay PDG code management
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_TOOLS_ParticleData
#define HAMMER_TOOLS_ParticleData

#ifndef PARTICLEDATA_FRIENDS
#define PARTICLEDATA_FRIENDS typedef bool DummyParticleDataFriend
#endif

#include <vector>

#include "Hammer/Particle.fhh"
#include "Hammer/Tools/Pdg.fhh"
#include "Hammer/IndexTypes.hh"

namespace Hammer {


    /// @brief PDG code process signature class
    ///
    /// Manages the particle signatures (hashes, pdg lists, particle masses) of an object that needs to be associated
    /// to a decay vertex, in order to associate it in a Process at runtime.
    /// Supports multiple signatures to lump together similar decays
    /// (e.g. \f$ B^0 \rightarrow D^- \bar \tau \nu_\tau \f$ with \f$ B^+ \rightarrow D^0 \bar \tau \nu_\tau \f$,
    /// \f$ B^+ \rightarrow D^0 \bar \mu \nu_\mu \f$, etc.). Each signature can be selected by a signature index.
    ///
    /// @ingroup Core
    class ParticleData {

        PARTICLEDATA_FRIENDS;

    public:

        /// default constructor
        ParticleData() = default;

        ParticleData(const ParticleData& other) = default;
        ParticleData& operator=(const ParticleData& other) = default;
        ParticleData(ParticleData&& other) = default;
        ParticleData& operator=(ParticleData&& other) = default;

        virtual ~ParticleData();

    protected:

        /// @brief Decay signature information
        ///
        /// Contains all the necessary information for a decay vertex
        ///
        /// @ingroup Core
        struct Signature {

            Signature() = default;

            PdgId parent; ///< particle ID of the parent particle
            std::vector<PdgId>
                daughters; ///< particle IDs of the (ordered) daughter (and granddaughter) particles.
                           ///< the ordering is: daughters first, then grand-daughters (if available); within each group
                           ///< different type of particles are ordered according to abs(PDG code) and according to
                           ///< ascending PDG code for particles/anti-particles of the same type
            HashId id;     ///< unique ID of the decay process
            HashId hadronicId;          ///< unique ID of the hadronic part of the decay
                                        ///< (parent + hadronic daughters) for form factor associations
            std::vector<double> masses; ///< list of particle mass values parent + daughters (in the same order as the
                                        ///< particle ID vector)

            /// @brief  computes the unique IDs from the PDG codes
            /// @param[in] numDaughters  number of daughters to consider for the hadronic signature calculation
            ///                          (necessary for processes which combines two decays, e.g. \f$ B\rightarrow
            ///                          D^*\ell\nu, D^*->D\pi \f$)
            void calcIds(size_t numDaughters);

            /// @brief  appends the mass of a particle to the mass list
            /// @param[in] id the particle PDG code
            void addMass(PdgId id);

            /// @brief appends the masses of a list of particles to the mass list
            /// @param[in] ids the list of PDG codes of the particles
            void addMasses(const std::vector<PdgId>& ids);
        };

    public:

        /// @brief method to evaluate the object on a specific particle set
        /// @param[in] parent  the parent Particle
        /// @param[in] daughters the daughters (and grand-daughters, if necessary) Particle list
        /// @param[in] references the parent Particle siblings (necessary e.g. for helicity amplitude phase conventions)
        virtual void eval(const Particle& parent, const ParticleList& daughters, const ParticleList& references) = 0;


        /// @brief select a specific signature to be the current signature
        /// @param[in] idx  the signature index
        virtual bool setSignatureIndex(size_t idx = 0);

        /// @brief  returns the unique ID of the current decay signature
        /// @return the UID
        [[nodiscard]] HashId id() const;

        /// @brief  returns the hadronic unique ID (parent + hadronic daughters) of the current decay signature
        /// @return the hadronic UID
        [[nodiscard]] HashId hadronicId() const;

        /// @brief returns the particle masses of the current decay signature
        /// @return  the list of masses
        [[nodiscard]] const std::vector<double>& masses() const;

        /// @brief  returns the number of available signatures
        /// @return  the number of signatures
        [[nodiscard]] size_t numSignatures() const;

    protected:

        /// @brief adds a signature to the list by specifying the particles PDG codes. The internal particle
        ///        ordering is automatically taken into account
        /// @param[in] parent  the parent PDG code
        /// @param[in] daughters the list of daughter PDG codes
        /// @param[in] grandDaughters  the list of subdaughters PDG codes for objects associated to 2-step decays
        void addProcessSignature(PdgId parent, const std::vector<PdgId>& daughters,
                                 const std::vector<PdgId>& grandDaughters = {});


    private:

        static bool checkValidSignature(PdgId parent, const std::vector<PdgId>& daughters,
                                        const std::vector<PdgId>& grandDaughters);

    protected:

        std::vector<Signature> _signatures; ///< the list of signatures
        size_t _signatureIndex{0ul};        ///< the index of the current signature
    };

} // namespace Hammer

#endif
