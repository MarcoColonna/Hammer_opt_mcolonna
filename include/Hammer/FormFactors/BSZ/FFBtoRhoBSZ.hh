///
/// @file  FFBtoRhoBSZ.hh
/// @brief \f$ B \rightarrow \rho \f$ BSZ form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BRHOBSZ
#define HAMMER_FF_BRHOBSZ

#include "Hammer/FormFactors/BSZ/FFBSZBase.hh"

namespace Hammer {

    class FFBtoRhoBSZ final : public FFBSZBase {

    public:

        FFBtoRhoBSZ();

        FFBtoRhoBSZ(const FFBtoRhoBSZ& other) = default;
        FFBtoRhoBSZ& operator=(const FFBtoRhoBSZ& other) = delete;
        FFBtoRhoBSZ(FFBtoRhoBSZ&& other) = delete;
        FFBtoRhoBSZ& operator=(FFBtoRhoBSZ&& other) = delete;
        ~FFBtoRhoBSZ() final = default;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;
    };

} // namespace Hammer

#endif
