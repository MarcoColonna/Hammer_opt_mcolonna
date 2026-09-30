///
/// @file  ParticleUtils.hh
/// @brief PDG codes to UID functions
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_TOOLS_ParticleUtils
#define HAMMER_TOOLS_ParticleUtils

#include <vector>
#include <functional>

#include "Hammer/IndexTypes.hh"
#include "Hammer/Tools/Pdg.fhh"

namespace Hammer {

    class Particle;

    /// @brief sorting function for PDG ids for computing hashes.
    ///        the codes are ordered if \f$ |a|<|b| \f$ or if \f$ |a|=|b| \f$ and if \f$ a<b \f$
    /// @param[in] a the first PDG code
    /// @param[in] b the second PDG code
    bool pdgSorter(PdgId a, PdgId b);

    /// @brief  compute a unique ID for a given process based on the PDG codes of the parent particle and
    ///         the ordered list of the daughters (and grandaughters if present)
    /// @param[in] parent       the PDG code of the parent particle
    /// @param[in] allDaughters the ordered list of PDG codes of the daughters particles
    /// @return the unique id
    HashId processID(PdgId parent, const std::vector<PdgId>& allDaughters);

    HashId combineProcessIDs(const std::set<HashId>& allIds);

    /// @brief combine list of codes of daughters and grandaughters (for processes which parameterise
    ///        two subsequent decays in a single amplitude) in a single list for computing hashes
    ///        the convention is that daughters are ordered within themselves, grandDaughters are ordered
    ///        within themselves, and then the two groups are concatenated with daughters first
    /// @param[in] daughters
    /// @param[in] subDaughters
    /// @return  the combined list
    std::vector<PdgId> combineDaughters(const std::vector<PdgId>& daughters,
                                        const std::vector<PdgId>& subDaughters = {});

    /// @brief return the PDG codes of the conjugate particles (itself if self-conjugate)
    ///        for all the PDG codes in alist
    /// @param[in] list  the PDG codes
    std::vector<PdgId> flipSigns(const std::vector<PdgId>& list);

    /// @brief return the PDG code of the conjugate particle (itself if self-conjugate)
    /// @param[in] id  the PDG code
    PdgId flipSign(const PdgId& id);

    /// @brief checks whether two particles are ordered according to the PDG code ordering used in computing hashes
    ///        (see `pdgSorter` for more info on the ordering)
    /// @param[in] pdgGetter  the PDG code extractor function
    /// @param[in] a the first particle
    /// @param[in] b the second particle
    bool particlesByPdg(const std::function<PdgId(const Particle&)>& pdgGetter, const Particle& a, const Particle& b);


} // namespace Hammer

#endif
