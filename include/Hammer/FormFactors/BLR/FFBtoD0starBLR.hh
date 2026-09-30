///
/// @file  FFBtoD0starBLR.hh
/// @brief \f$ B \rightarrow D_0^* \f$ BLR form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BDSSD0STARBLR
#define HAMMER_FF_BDSSD0STARBLR

#include "Hammer/FormFactors/BLR/FFBLRBase.hh"

namespace Hammer {

    class FFBtoD0starBLR final : public FFBLRBase {

    public:

        FFBtoD0starBLR();

        FFBtoD0starBLR(const FFBtoD0starBLR& other) = default;
        FFBtoD0starBLR& operator=(const FFBtoD0starBLR& other) = delete;
        FFBtoD0starBLR(FFBtoD0starBLR&& other) = delete;
        FFBtoD0starBLR& operator=(FFBtoD0starBLR&& other) = delete;
        ~FFBtoD0starBLR() final = default;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;
    };

} // namespace Hammer

#endif
