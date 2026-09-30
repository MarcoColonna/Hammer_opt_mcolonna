///
/// @file  FFBLPRXPVarBase.hh
/// @brief Hammer base class for BLPRXPVar form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BLPRXPVAR_BASE
#define HAMMER_FF_BLPRXPVAR_BASE

#ifndef FFBLPRXPVARBASE_FRIENDS
#define FFBLPRXPVARBASE_FRIENDS typedef bool DummyFriend
#endif

#include <array>

#include "Hammer/FormFactors/BLPRXP/FFBLPRXPBase.hh"

namespace Hammer {

    /// @brief Base class for BLPRXPVar form factors
    ///
    /// @ingroup FormFactors
    class FFBLPRXPVarBase : public FFBLPRXPBase {

        FFBLPRXPVARBASE_FRIENDS;

    public:

        FFBLPRXPVarBase();

        FFBLPRXPVarBase(const FFBLPRXPVarBase& other) = default;
        FFBLPRXPVarBase& operator=(const FFBLPRXPVarBase& other) = delete;
        FFBLPRXPVarBase(FFBLPRXPVarBase&& other) = delete;
        FFBLPRXPVarBase& operator=(FFBLPRXPVarBase&& other) = delete;
        ~FFBLPRXPVarBase() override = default;

    protected:

        /// @brief
        void defineSettings() override;

        std::array<double, 10> XiLinearCoeffs(double w) const;
    };

} // namespace Hammer

#endif
