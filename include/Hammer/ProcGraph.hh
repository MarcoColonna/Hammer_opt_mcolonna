///
/// @file  ProcGraph.hh
/// @brief Container class for process tree structure and its amplitudes associations
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_PROCESSGRAPH_HH
#define HAMMER_PROCESSGRAPH_HH

#include <map>
#include <vector>
#include <tuple>

#include "Hammer/ProcGraph.fhh"
#include "Hammer/Particle.hh"

namespace Hammer {

    class Log;
    class Process;
    class DictionaryManager;

    template <typename T>
    using EdgeWeightDict = std::map<size_t, T, std::less<>>;

    template <typename T>
    using DepthLevelDict = std::map<size_t, T, std::greater<size_t>>;


    class AmplTriplet {
    public:

        [[nodiscard]] AmplEntry resolvePurePS(bool parentPS, bool daughterPS) const;

        AmplitudeBase* parent;
        AmplitudeBase* daughter;
        AmplitudeBase* edge;

    private:

        /// @brief logging facility
        /// @return   stream to be used for logging
        static Log& getLog();
    };

    struct EdgeEntry {
        ParticleIndex parentIdx;
        ParticleIndex daughterIdx;
        AmplTriplet amplitudes;
        NumDenPair<bool> parentPSFlags;
        NumDenPair<bool> daughterPSFlags;
    };

    struct VertexEntry {
        ParticleIndex parentIdx;
        AmplitudeBase* amplitude;
        NumDenPair<bool> parentPSFlags;
    };

    /// @brief Decay process class
    ///
    /// Contains the amplitudes, weights and info associated to a decay, ...
    ///
    /// @ingroup Core
    class ProcGraph {

    public:

        ProcGraph() = default;

        ProcGraph(const ProcGraph& other) = delete;
        ProcGraph& operator=(const ProcGraph& other) = delete;
        ProcGraph(ProcGraph&& other) = delete;
        ProcGraph& operator=(ProcGraph&& other) = delete;

        ~ProcGraph() noexcept = default;

        using EdgeDictType = EdgeWeightDict<DepthLevelDict<std::vector<EdgeEntry>>>;
        using VertexListType = std::vector<VertexEntry>;
        using SelAmplitudeDict = VertexDict<SelectedAmplEntry>;

        /// @brief
        /// @param[in] dictionaries
        /// @param[in] inputs
        void initialize(const DictionaryManager* dictionaries, const Process* inputs);

        /// @brief
        /// @param[in] ampltype
        /// @param[in] parent
        /// @param[in] daughter
        [[nodiscard]] std::tuple<Particle, ParticleList, ParticleList>
        getParticleVectors(AmplType ampltype, ParticleIndex parent, ParticleIndex daughter) const;

        ///@brief
        /// @param[in] ampl
        /// @param[in] ancestorAmpl
        /// @param[in] isAncestorGrandparent
        [[nodiscard]] static double getMultFactor(AmplEntry ampl, AmplEntry ancestorAmpl, bool isAncestorGrandparent);

        ///@brief
        /// @param[in] ampl
        /// @param[in] parent
        /// @param[in] daughter
        [[nodiscard]] double getMassFactor(AmplEntry ampl, ParticleIndex parent, ParticleIndex daughter) const;

        [[nodiscard]] const SelAmplitudeDict& selectedAmplitudes() const;
        [[nodiscard]] const UniqueParticleIndices& assignedVertices() const;

    protected:

        /// @brief logging facility
        /// @return   stream to be used for logging
        static Log& getLog();

        /// @brief
        void makeDepthMap(ParticleIndex parent, size_t seed);

        void makeVertexEdgeTrees();

        void selectEdgeVertices();

        [[nodiscard]] bool isAssignedVertex(ParticleIndex index) const;

    private:

        VertexDict<size_t> _depthMap;
        EdgeDictType _implementedEdges;
        VertexListType _implementedVertices;

        SelAmplitudeDict _selectedAmplitudes;
        UniqueParticleIndices _assignedVertices;

        const Process* _inputs;
        const DictionaryManager* _dictionaries;
    };

} // namespace Hammer

#endif
