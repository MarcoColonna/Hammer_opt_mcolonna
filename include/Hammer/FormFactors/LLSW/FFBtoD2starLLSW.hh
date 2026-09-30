///
/// @file  FFBtoD2starLLSW.hh
/// @brief \f$ B \rightarrow D_2^* \f$ LLSW form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BDSSD2STARLLSW
#define HAMMER_FF_BDSSD2STARLLSW

#include "Hammer/FormFactors/LLSW/FFLLSWBase.hh"

namespace Hammer {

    class FFBtoD2starLLSW final : public FFLLSWBase {

    public:

        FFBtoD2starLLSW();

        FFBtoD2starLLSW(const FFBtoD2starLLSW& other) = default;
        FFBtoD2starLLSW& operator=(const FFBtoD2starLLSW& other) = delete;
        FFBtoD2starLLSW(FFBtoD2starLLSW&& other) = delete;
        FFBtoD2starLLSW& operator=(FFBtoD2starLLSW&& other) = delete;
        ~FFBtoD2starLLSW() final = default;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;
    };

} // namespace Hammer

#endif
