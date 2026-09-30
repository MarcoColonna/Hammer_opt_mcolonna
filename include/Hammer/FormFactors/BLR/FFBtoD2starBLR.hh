///
/// @file  FFBtoD2starBLR.hh
/// @brief \f$ B \rightarrow D_2^* \f$ BLR form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BDSSD2STARBLR
#define HAMMER_FF_BDSSD2STARBLR

#include "Hammer/FormFactors/BLR/FFBLRBase.hh"

namespace Hammer {

    class FFBtoD2starBLR final : public FFBLRBase {

    public:

        FFBtoD2starBLR();

        FFBtoD2starBLR(const FFBtoD2starBLR& other) = default;
        FFBtoD2starBLR& operator=(const FFBtoD2starBLR& other) = delete;
        FFBtoD2starBLR(FFBtoD2starBLR&& other) = delete;
        FFBtoD2starBLR& operator=(FFBtoD2starBLR&& other) = delete;
        ~FFBtoD2starBLR() final = default;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;
    };

} // namespace Hammer

#endif
