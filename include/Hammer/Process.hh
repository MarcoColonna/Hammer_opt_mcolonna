///
/// @file  Process.hh
/// @brief Hammer process class
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_PROCESS_HH
#define HAMMER_PROCESS_HH

#include <limits>
#include <map>
#include <string>
#include <vector>
#include <unordered_map>

#include "Hammer/IndexTypes.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/SettingsConsumer.hh"
#include "Hammer/Tools/HammerSerial.fhh"

namespace Hammer {

    class Log;

    /// @brief Decay process class
    ///
    /// Contains the amplitudes, weights and info associated to a decay, ...
    ///
    /// @ingroup Core
    class Process : public SettingsConsumer {

    public:

        Process();

        Process(const Process& other) = default;
        Process& operator=(const Process& other) = default;
        Process(Process&& other) noexcept = default; // NOLINT(bugprone-exception-escape)
        Process& operator=(Process&& other) noexcept = default;

        ~Process() noexcept override = default;

        using TreeMap = std::map<ParticleIndex, ParticleIndices>;
        using const_iterator = TreeMap::const_iterator;

        [[nodiscard]] const_iterator begin() const;

        [[nodiscard]] const_iterator end() const;

        [[nodiscard]] const_iterator find(ParticleIndex particle) const;


        /// @brief
        /// @param[in] p
        /// @return
        ParticleIndex addParticle(const Particle& p);

        /// @brief
        /// @param[in] parent
        /// @param[in] daughters
        void addVertex(ParticleIndex parent, const ParticleIndices& daughters);

        /// @brief
        /// @param[in] parent
        /// @param[in] prune
        void removeVertex(ParticleIndex parent, bool prune = false);

        /// @brief
        /// @param[in] parent
        /// @return
        [[nodiscard]] const ParticleIndices& getDaughtersIds(ParticleIndex parent = 0) const;

        /// @brief
        /// @param[in] daughter
        /// @return
        [[nodiscard]] ParticleIndex getParentId(ParticleIndex daughter) const;

        /// @brief
        /// @param[in] particle
        /// @return
        [[nodiscard]] ParticleIndices getSiblingsIds(ParticleIndex particle = 0) const;

        /// @brief
        /// @param[in] id
        /// @return
        [[nodiscard]] const Particle& getParticle(ParticleIndex id) const;

        /// @brief
        /// @param[in] parent
        /// @param[in] sorted
        /// @return
        [[nodiscard]] ParticleList getDaughters(ParticleIndex parent = 0, bool sorted = false) const;

        /// @brief
        /// @param[in] particle
        /// @param[in] sorted
        /// @return
        [[nodiscard]] ParticleList getSiblings(ParticleIndex particle = 0, bool sorted = false) const;

        /// @brief
        /// @param[in] daughter
        /// @return
        [[nodiscard]] const Particle& getParent(ParticleIndex daughter) const;

        [[nodiscard]] bool isParent(ParticleIndex particle) const;

        [[nodiscard]] ParticleIndex getFirstVertex() const;

        /// @brief
        /// @return
        [[nodiscard]] HashId getId() const;

        [[nodiscard]] const std::set<HashId>& fullId() const;

        /// @brief
        /// @param[in] parent
        /// @param[in] daughters
        /// @return
        [[nodiscard]] std::pair<Particle, ParticleList> getParticlesByVertex(PdgId parent,
                                                                             const std::vector<PdgId>& daughters) const;

        /// @brief
        /// @param[in] vertex
        /// @return
        [[nodiscard]] std::pair<Particle, ParticleList> getParticlesByVertex(const std::string& vertex) const;

        /// @brief
        /// @param[in] vertex
        /// @return
        [[nodiscard]] HashId getVertexId(const std::string& vertex) const;

        bool initialize();

        void disable();

        /// @brief
        /// @param[in] msgwriter
        /// @param[in] msg
        void write(flatbuffers::FlatBufferBuilder* msgwriter, flatbuffers::Offset<Serial::FBProcIDs>* msg) const;

        /// @brief
        /// @param[in] msgreader
        void read(const Serial::FBProcIDs* msgreader);

        [[nodiscard]] size_t numParticles(bool withoutPhotons = false) const;

    private:

        /// @brief purely virtual function for a class to define new settings
        void defineSettings() override;

        /// @brief logging facility
        /// @return   stream to be used for logging
        static Log& getLog();

        // ///@brief
        // ///@TODO: implement or remove
        // void expandSoftPhotonVertices();

        ///@brief
        void pruneSoftPhotons();

        /// @brief
        void calcSignatures();

        [[nodiscard]] ParticleIndex findFirstVertex() const;
        void cacheParticleDependencies();

        [[nodiscard]] ParticleList sortByParent(ParticleIndex parent, ParticleList parts) const;


        ParticleList _particles;
        TreeMap _processTree;

        ProcessUID _hashId{0};
        VertexUIDSet _fullId;
        ProcIdDict<ParticleIndex> _idTree;

        ParticleIndex _firstVertexIdx = std::numeric_limits<size_t>::max();
        std::unordered_map<ParticleIndex, ParticleIndex> _parentDict;
        UniqueParticleIndices _erasedGammas;
        bool _initialized = false;
    };

    inline ParticleIndex parentId(Process::const_iterator it) {
        return it->first;
    }

    inline const ParticleIndices& daughtersId(Process::const_iterator it) {
        return it->second;
    }

} // namespace Hammer

#endif
