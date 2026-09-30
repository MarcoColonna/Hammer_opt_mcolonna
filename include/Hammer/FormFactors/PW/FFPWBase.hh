///
/// @file  FFPWBase.hh
/// @brief Hammer base class for partial wave coefficients
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_PW_BASE
#define HAMMER_FF_PW_BASE

#include "Hammer/FormFactorBase.hh"

namespace Hammer {

    /// @brief Base class for PW form factors
    /// See Pervin, Roberst, and Capstick, Phys. Rev. C72, 035201 (2005)
    /// Implementation matched to EvtGen
    ///
    /// @ingroup FormFactors
    class FFPWBase : public FormFactorBase {

    public:

        FFPWBase();

        FFPWBase(const FFPWBase& other) = default;
        FFPWBase& operator=(const FFPWBase& other) = delete;
        FFPWBase(FFPWBase&& other) = delete;
        FFPWBase& operator=(FFPWBase&& other) = delete;
        ~FFPWBase() override = default;

    protected:

        /// @brief
        /// @param[in] point
        /// @param[in] masses
        /// @return
        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) override = 0;

        /// @brief
        void defineSettings() override = 0;

        //    protected:
        //        virtual void addRefs() const override;
    };

} // namespace Hammer

#endif
