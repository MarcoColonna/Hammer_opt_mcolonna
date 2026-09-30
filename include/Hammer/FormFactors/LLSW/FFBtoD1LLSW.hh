///
/// @file  FFBtoD1LLSW.hh
/// @brief \f$ B \rightarrow D_1 \f$ LLSW form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BDSSD1LLSW
#define HAMMER_FF_BDSSD1LLSW

#include "Hammer/FormFactors/LLSW/FFLLSWBase.hh"

namespace Hammer {

    class FFBtoD1LLSW final : public FFLLSWBase {

    public:

        FFBtoD1LLSW();

        FFBtoD1LLSW(const FFBtoD1LLSW& other) = default;
        FFBtoD1LLSW& operator=(const FFBtoD1LLSW& other) = delete;
        FFBtoD1LLSW(FFBtoD1LLSW&& other) = delete;
        FFBtoD1LLSW& operator=(FFBtoD1LLSW&& other) = delete;
        ~FFBtoD1LLSW() final = default;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;
    };

} // namespace Hammer

#endif
