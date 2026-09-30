///
/// @file  FFBtoDBLPRXP.hh
/// @brief \f$ B \rightarrow D \f$ BLPRXP form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BDBLPRXP
#define HAMMER_FF_BDBLPRXP

#include "Hammer/FormFactors/BLPRXP/FFBLPRXPBase.hh"

namespace Hammer {

    class FFBtoDBLPRXP final : public FFBLPRXPBase {

    public:

        FFBtoDBLPRXP();

        FFBtoDBLPRXP(const FFBtoDBLPRXP& other) = default;
        FFBtoDBLPRXP& operator=(const FFBtoDBLPRXP& other) = delete;
        FFBtoDBLPRXP(FFBtoDBLPRXP&& other) = delete;
        FFBtoDBLPRXP& operator=(FFBtoDBLPRXP&& other) = delete;
        ~FFBtoDBLPRXP() final = default;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;
    };

} // namespace Hammer

#endif
