///
/// @file  FFBLRBase.hh
/// @brief Hammer base class for BLR form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BLR_BASE
#define HAMMER_FF_BLR_BASE

#include "Hammer/FormFactorBase.hh"
#include "Hammer/FormFactors/AsCoeffs.hh"

namespace Hammer {

    /// @brief Base class for BLR form factors
    ///
    /// @ingroup FormFactors
    class FFBLRBase : public FF1to1Base {

    public:

        FFBLRBase();

        FFBLRBase(const FFBLRBase& other) = default;
        FFBLRBase& operator=(const FFBLRBase& other) = delete;
        FFBLRBase(FFBLRBase&& other) = delete;
        FFBLRBase& operator=(FFBLRBase&& other) = delete;
        ~FFBLRBase() override = default;

    protected:

        /// @brief
        /// @param[in] point
        /// @param[in] masses
        /// @return
        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) override = 0;

        /// @brief
        void defineSettings() override = 0;

        void addRefs() const override;

        AsCoeffs _asCorrections;
    };

} // namespace Hammer

#endif
