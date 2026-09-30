///
/// @file  FFBtoD0starISGW2.hh
/// @brief \f$ B \rightarrow D_0^* \f$ ISGW2 form factors
/// @brief Ported directly from EvtGen
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BDSSD0STARISGW2
#define HAMMER_FF_BDSSD0STARISGW2

#include "Hammer/FormFactors/ISGW2/FFISGW2Base.hh"

namespace Hammer {

    class FFBtoD0starISGW2 final : public FFISGW2Base {

    public:

        FFBtoD0starISGW2();

        FFBtoD0starISGW2(const FFBtoD0starISGW2& other) = default;
        FFBtoD0starISGW2& operator=(const FFBtoD0starISGW2& other) = delete;
        FFBtoD0starISGW2(FFBtoD0starISGW2&& other) = delete;
        FFBtoD0starISGW2& operator=(FFBtoD0starISGW2&& other) = delete;
        ~FFBtoD0starISGW2() final = default;


        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;
    };

} // namespace Hammer

#endif
