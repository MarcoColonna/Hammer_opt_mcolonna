///
/// @file  FFBtoDBLPR.hh
/// @brief \f$ B \rightarrow D \f$ BLPR form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BDBLPR
#define HAMMER_FF_BDBLPR

#include "Hammer/FormFactors/BLPR/FFBLPRBase.hh"

namespace Hammer {

    class FFBtoDBLPR final : public FFBLPRBase {

    public:

        FFBtoDBLPR();

        FFBtoDBLPR(const FFBtoDBLPR& other) = default;
        FFBtoDBLPR& operator=(const FFBtoDBLPR& other) = delete;
        FFBtoDBLPR(FFBtoDBLPR&& other) = delete;
        FFBtoDBLPR& operator=(FFBtoDBLPR&& other) = delete;
        ~FFBtoDBLPR() final = default;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;
    };

} // namespace Hammer

#endif
