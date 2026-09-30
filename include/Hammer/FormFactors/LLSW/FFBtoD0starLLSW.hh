///
/// @file  FFBtoD0starLLSW.hh
/// @brief \f$ B \rightarrow D_0^* \f$ LLSW form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BDSSD0STARLLSW
#define HAMMER_FF_BDSSD0STARLLSW

#include "Hammer/FormFactors/LLSW/FFLLSWBase.hh"

namespace Hammer {

    class FFBtoD0starLLSW final : public FFLLSWBase {

    public:

        FFBtoD0starLLSW();

        FFBtoD0starLLSW(const FFBtoD0starLLSW& other) = default;
        FFBtoD0starLLSW& operator=(const FFBtoD0starLLSW& other) = delete;
        FFBtoD0starLLSW(FFBtoD0starLLSW&& other) = delete;
        FFBtoD0starLLSW& operator=(FFBtoD0starLLSW&& other) = delete;
        ~FFBtoD0starLLSW() final = default;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;
    };

} // namespace Hammer

#endif
