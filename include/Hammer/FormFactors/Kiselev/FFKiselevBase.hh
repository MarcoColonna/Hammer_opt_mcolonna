///
/// @file  FFKiselevBase.hh
/// @brief Hammer base class for Kiselev form factors.
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_Kiselev_BASE
#define HAMMER_FF_Kiselev_BASE

#include "Hammer/FormFactorBase.hh"

namespace Hammer {

    /// @brief Base class for Kiselev form factors
    /// See Kiselev, hep-ph/0211021
    /// Implementation matched to EvtGen
    ///
    /// @ingroup FormFactors
    class FFKiselevBase : public FF1to1Base {

    public:

        FFKiselevBase();

        FFKiselevBase(const FFKiselevBase& other) = default;
        FFKiselevBase& operator=(const FFKiselevBase& other) = delete;
        FFKiselevBase(FFKiselevBase&& other) = delete;
        FFKiselevBase& operator=(FFKiselevBase&& other) = delete;
        ~FFKiselevBase() override = default;

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
