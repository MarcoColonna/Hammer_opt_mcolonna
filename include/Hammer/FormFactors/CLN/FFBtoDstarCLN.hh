///
/// @file  FFBtoDstarCLN.hh
/// @brief \f$ B \rightarrow D^* \f$ CLN form factors
/// @brief Ported from EvtGen HQET2
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BDSTARCLN
#define HAMMER_FF_BDSTARCLN

#include "Hammer/FormFactors/CLN/FFCLNBase.hh"

namespace Hammer {

    class FFBtoDstarCLN final : public FFCLNBase {

    public:

        FFBtoDstarCLN();

        FFBtoDstarCLN(const FFBtoDstarCLN& other) = default;
        FFBtoDstarCLN& operator=(const FFBtoDstarCLN& other) = delete;
        FFBtoDstarCLN(FFBtoDstarCLN&& other) = delete;
        FFBtoDstarCLN& operator=(FFBtoDstarCLN&& other) = delete;
        ~FFBtoDstarCLN() final = default;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;
    };

} // namespace Hammer

#endif
