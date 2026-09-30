///
/// @file  FFLSPRVarBase.hh
/// @brief Hammer base class for LSPRVar form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_LSPRVAR_BASE
#define HAMMER_FF_LSPRVAR_BASE

#include <vector>

#include "Hammer/FormFactors/LSPR/FFLSPRBase.hh"

namespace Hammer {

    /// @brief Base class for LSPRVar form factors
    ///
    /// @ingroup FormFactors
    class FFLSPRVarBase : public FFLSPRBase {

    public:

        FFLSPRVarBase();

        FFLSPRVarBase(const FFLSPRVarBase& other) = default;
        FFLSPRVarBase& operator=(const FFLSPRVarBase& other) = delete;
        FFLSPRVarBase(FFLSPRVarBase&& other) = delete;
        FFLSPRVarBase& operator=(FFLSPRVarBase&& other) = delete;
        ~FFLSPRVarBase() override = default;

    protected:

        /// @brief
        void defineSettings() override;

        std::vector<double> XiLinearCoeffs(double w) const;
    };

} // namespace Hammer

#endif
