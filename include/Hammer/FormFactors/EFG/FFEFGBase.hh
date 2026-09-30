///
/// @file  FFEFGBase.hh
/// @brief Hammer base class for EFG form factors.
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_EFG_BASE
#define HAMMER_FF_EFG_BASE

#include "Hammer/FormFactorBase.hh"

namespace Hammer {

    /// @brief Base class for EFG form factors
    /// See Ebert, Faustov, and Galkin, Phys.Rev.D 68 (2003) 094020
    /// Implementation matched to EvtGen
    ///
    /// @ingroup FormFactors
    class FFEFGBase : public FF1to1Base {

    public:

        FFEFGBase();

        FFEFGBase(const FFEFGBase& other) = default;
        FFEFGBase& operator=(const FFEFGBase& other) = delete;
        FFEFGBase(FFEFGBase&& other) = delete;
        FFEFGBase& operator=(FFEFGBase&& other) = delete;
        ~FFEFGBase() override = default;

    protected:

        /// @brief
        /// @param[in] point
        /// @param[in] masses
        /// @return
        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) override = 0;

        /// @brief
        void defineSettings() override = 0;

        void addRefs() const override;
    };

} // namespace Hammer

#endif
