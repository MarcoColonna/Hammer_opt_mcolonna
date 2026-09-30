///
/// @file  FFBtoD1starLLSW.hh
/// @brief \f$ B \rightarrow D_1^* \f$ LLSW form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BDSSD1STARLLSW
#define HAMMER_FF_BDSSD1STARLLSW

#include "Hammer/FormFactors/LLSW/FFLLSWBase.hh"

namespace Hammer {

    class FFBtoD1starLLSW final : public FFLLSWBase {

    public:

        FFBtoD1starLLSW();

        FFBtoD1starLLSW(const FFBtoD1starLLSW& other) = default;
        FFBtoD1starLLSW& operator=(const FFBtoD1starLLSW& other) = delete;
        FFBtoD1starLLSW(FFBtoD1starLLSW&& other) = delete;
        FFBtoD1starLLSW& operator=(FFBtoD1starLLSW&& other) = delete;
        ~FFBtoD1starLLSW() final = default;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;
    };

} // namespace Hammer

#endif
