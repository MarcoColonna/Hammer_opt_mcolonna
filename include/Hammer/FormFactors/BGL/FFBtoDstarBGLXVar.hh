///
/// @file  FFBtoDstarBGLXVar.hh
/// @brief \f$ B \rightarrow D \f$ BGL Extended form factors with variations
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BDSTARBGLXVAR
#define HAMMER_FF_BDSTARBGLXVAR

#include "Hammer/FormFactors/BGL/FFBGLXBase.hh"

namespace Hammer {

    class FFBtoDstarBGLXVar final : public FFBGLXBase {

    public:

        FFBtoDstarBGLXVar();

        FFBtoDstarBGLXVar(const FFBtoDstarBGLXVar& other) = default;
        FFBtoDstarBGLXVar& operator=(const FFBtoDstarBGLXVar& other) = delete;
        FFBtoDstarBGLXVar(FFBtoDstarBGLXVar&& other) = delete;
        FFBtoDstarBGLXVar& operator=(FFBtoDstarBGLXVar&& other) = delete;

        ~FFBtoDstarBGLXVar() final = default;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;
    };

} // namespace Hammer

#endif
