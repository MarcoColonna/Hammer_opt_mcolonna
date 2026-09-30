///
/// @file  FFBtoDstarBGLVar.hh
/// @brief \f$ B \rightarrow D \f$ BGL form factors with variations
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BDSTARBGLVAR
#define HAMMER_FF_BDSTARBGLVAR

#include "Hammer/FormFactors/BGL/FFBGLBase.hh"

namespace Hammer {

    class FFBtoDstarBGLVar final : public FFBGLBase {

    public:

        FFBtoDstarBGLVar();

        FFBtoDstarBGLVar(const FFBtoDstarBGLVar& other) = default;
        FFBtoDstarBGLVar& operator=(const FFBtoDstarBGLVar& other) = delete;
        FFBtoDstarBGLVar(FFBtoDstarBGLVar&& other) = delete;
        FFBtoDstarBGLVar& operator=(FFBtoDstarBGLVar&& other) = delete;

        ~FFBtoDstarBGLVar() final = default;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;
    };

} // namespace Hammer

#endif
