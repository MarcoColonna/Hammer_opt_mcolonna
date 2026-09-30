///
/// @file  FFBCLBase.hh
/// @brief Hammer base class for BCL form factors.
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BCL_BASE
#define HAMMER_FF_BCL_BASE

#include "Hammer/FormFactorBase.hh"

namespace Hammer {

    /// @brief Base class for BCL form factors
    /// Implementation matched to EvtGen
    ///
    /// @ingroup FormFactors
    class FFBCLBase : public FF1to1Base {

    public:

        FFBCLBase();

        FFBCLBase(const FFBCLBase& other) = default;
        FFBCLBase& operator=(const FFBCLBase& other) = delete;
        FFBCLBase(FFBCLBase&& other) = delete;
        FFBCLBase& operator=(FFBCLBase&& other) = delete;
        ~FFBCLBase() override = default;

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
