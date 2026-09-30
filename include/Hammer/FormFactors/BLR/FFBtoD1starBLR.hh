///
/// @file  FFBtoD1starBLR.hh
/// @brief \f$ B \rightarrow D_1^* \f$ BLR form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BDSSD1STARBLR
#define HAMMER_FF_BDSSD1STARBLR

#include "Hammer/FormFactors/BLR/FFBLRBase.hh"

namespace Hammer {

    class FFBtoD1starBLR final : public FFBLRBase {

    public:

        FFBtoD1starBLR();

        FFBtoD1starBLR(const FFBtoD1starBLR& other) = default;
        FFBtoD1starBLR& operator=(const FFBtoD1starBLR& other) = delete;
        FFBtoD1starBLR(FFBtoD1starBLR&& other) = delete;
        FFBtoD1starBLR& operator=(FFBtoD1starBLR&& other) = delete;
        ~FFBtoD1starBLR() final = default;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;
    };

} // namespace Hammer

#endif
