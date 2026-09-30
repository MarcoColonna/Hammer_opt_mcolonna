///
/// @file  FFCLNBase.hh
/// @brief Hammer base class for CLN form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_CLN_BASE
#define HAMMER_FF_CLN_BASE

#include "Hammer/FormFactorBase.hh"

namespace Hammer {

    /// @brief Base class for CLN form factors
    ///
    /// @ingroup FormFactors
    class FFCLNBase : public FF1to1Base {

    public:

        FFCLNBase();

        FFCLNBase(const FFCLNBase& other) = default;
        FFCLNBase& operator=(const FFCLNBase& other) = delete;
        FFCLNBase(FFCLNBase&& other) = delete;
        FFCLNBase& operator=(FFCLNBase&& other) = delete;
        ~FFCLNBase() override = default;

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
