///
/// @file  FFBtoD1BLR.hh
/// @brief \f$ B \rightarrow D_1 \f$ BLR form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BDSSD1BLR
#define HAMMER_FF_BDSSD1BLR

#include "Hammer/FormFactors/BLR/FFBLRBase.hh"

namespace Hammer {

    class FFBtoD1BLR final : public FFBLRBase {

    public:

        FFBtoD1BLR();

        FFBtoD1BLR(const FFBtoD1BLR& other) = default;
        FFBtoD1BLR& operator=(const FFBtoD1BLR& other) = delete;
        FFBtoD1BLR(FFBtoD1BLR&& other) = delete;
        FFBtoD1BLR& operator=(FFBtoD1BLR&& other) = delete;
        ~FFBtoD1BLR() final = default;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;
    };

} // namespace Hammer

#endif
