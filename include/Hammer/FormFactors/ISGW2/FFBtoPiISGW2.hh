///
/// @file  FFBtoPiISGW2.hh
/// @brief \f$ B \rightarrow \pi \f$ ISGW2 form factors
/// @brief Ported directly from EvtGen
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BPIISGW2
#define HAMMER_FF_BPIISGW2

#include "Hammer/FormFactors/ISGW2/FFISGW2Base.hh"

namespace Hammer {

    class FFBtoPiISGW2 final : public FFISGW2Base {

    public:

        FFBtoPiISGW2();

        FFBtoPiISGW2(const FFBtoPiISGW2& other) = default;
        FFBtoPiISGW2& operator=(const FFBtoPiISGW2& other) = delete;
        FFBtoPiISGW2(FFBtoPiISGW2&& other) = delete;
        FFBtoPiISGW2& operator=(FFBtoPiISGW2&& other) = delete;
        ~FFBtoPiISGW2() final = default;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;
    };

} // namespace Hammer

#endif
