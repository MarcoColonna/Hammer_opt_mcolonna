///
/// @file  FFPCRBase.hh
/// @brief Hammer base class for PCR form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_PCR_BASE
#define HAMMER_FF_PCR_BASE

#include "Hammer/FormFactorBase.hh"

namespace Hammer {

    /// @brief Base class for PCR form factors
    /// See Pervin, Roberst, and Capstick, Phys. Rev. C72, 035201 (2005)
    /// Implementation matched to EvtGen
    ///
    /// @ingroup FormFactors
    class FFPCRBase : public FF1to1Base {

    public:

        FFPCRBase();

        FFPCRBase(const FFPCRBase& other) = default;
        FFPCRBase& operator=(const FFPCRBase& other) = delete;
        FFPCRBase(FFPCRBase&& other) = delete;
        FFPCRBase& operator=(FFPCRBase&& other) = delete;
        ~FFPCRBase() override = default;

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
