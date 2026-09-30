///
/// @file  FFBtoDstarBLPR.hh
/// @brief \f$ B \rightarrow D^* \f$ BLPR form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BDSTARBLPR
#define HAMMER_FF_BDSTARBLPR

#include "Hammer/FormFactors/BLPR/FFBLPRBase.hh"

namespace Hammer {

    class FFBtoDstarBLPR final : public FFBLPRBase {

    public:

        FFBtoDstarBLPR();

        FFBtoDstarBLPR(const FFBtoDstarBLPR& other) = default;
        FFBtoDstarBLPR& operator=(const FFBtoDstarBLPR& other) = delete;
        FFBtoDstarBLPR(FFBtoDstarBLPR&& other) = delete;
        FFBtoDstarBLPR& operator=(FFBtoDstarBLPR&& other) = delete;
        ~FFBtoDstarBLPR() final = default;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;
    };

} // namespace Hammer

#endif
