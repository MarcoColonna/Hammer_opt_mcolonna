///
/// @file  FFBtoOmegaBSZ.hh
/// @brief \f$ B \rightarrow \omega \f$ BSZ form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BOMEGABSZ
#define HAMMER_FF_BOMEGABSZ

#include "Hammer/FormFactors/BSZ/FFBSZBase.hh"

namespace Hammer {

    class FFBtoOmegaBSZ final : public FFBSZBase {

    public:

        FFBtoOmegaBSZ();

        FFBtoOmegaBSZ(const FFBtoOmegaBSZ& other) = default;
        FFBtoOmegaBSZ& operator=(const FFBtoOmegaBSZ& other) = delete;
        FFBtoOmegaBSZ(FFBtoOmegaBSZ&& other) = delete;
        FFBtoOmegaBSZ& operator=(FFBtoOmegaBSZ&& other) = delete;
        ~FFBtoOmegaBSZ() final = default;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;
    };

} // namespace Hammer

#endif
