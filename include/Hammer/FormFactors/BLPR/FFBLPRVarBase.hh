///
/// @file  FFBLPRVarBase.hh
/// @brief Hammer base class for BLPRVar form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BLPRVAR_BASE
#define HAMMER_FF_BLPRVAR_BASE

#include <vector>

#include "Hammer/FormFactors/BLPR/FFBLPRBase.hh"

namespace Hammer {

    /// @brief Base class for BLPRVar form factors
    ///
    /// @ingroup FormFactors
    class FFBLPRVarBase : public FFBLPRBase {

    public:

        FFBLPRVarBase();

        FFBLPRVarBase(const FFBLPRVarBase& other) = default;
        FFBLPRVarBase& operator=(const FFBLPRVarBase& other) = delete;
        FFBLPRVarBase(FFBLPRVarBase&& other) = delete;
        FFBLPRVarBase& operator=(FFBLPRVarBase&& other) = delete;
        ~FFBLPRVarBase() override = default;

    protected:

        /// @brief
        void defineSettings() override;

        std::vector<double> XiLinearCoeffs(double w) const;
    };

} // namespace Hammer

#endif
