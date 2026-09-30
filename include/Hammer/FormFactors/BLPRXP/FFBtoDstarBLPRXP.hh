///
/// @file  FFBtoDstarBLPRXP.hh
/// @brief \f$ B \rightarrow D^* \f$ BLPRXP form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BDSTARBLPRXP
#define HAMMER_FF_BDSTARBLPRXP

#include "Hammer/FormFactors/BLPRXP/FFBLPRXPBase.hh"

namespace Hammer {

    class FFBtoDstarBLPRXP final : public FFBLPRXPBase {

    public:

        FFBtoDstarBLPRXP();

        FFBtoDstarBLPRXP(const FFBtoDstarBLPRXP& other) = default;
        FFBtoDstarBLPRXP& operator=(const FFBtoDstarBLPRXP& other) = delete;
        FFBtoDstarBLPRXP(FFBtoDstarBLPRXP&& other) = delete;
        FFBtoDstarBLPRXP& operator=(FFBtoDstarBLPRXP&& other) = delete;
        ~FFBtoDstarBLPRXP() final = default;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;
    };

} // namespace Hammer

#endif
