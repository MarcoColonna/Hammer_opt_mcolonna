///
/// @file  FFBtoDBGLVar.hh
/// @brief \f$ B \rightarrow D \f$ BGL form factors with variations
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BDBGLVAR
#define HAMMER_FF_BDBGLVAR

#include "Hammer/FormFactors/BGL/FFBGLBase.hh"

namespace Hammer {

    class FFBtoDBGLVar final : public FFBGLBase {

    public:

        FFBtoDBGLVar();

        FFBtoDBGLVar(const FFBtoDBGLVar& other) = default;
        FFBtoDBGLVar& operator=(const FFBtoDBGLVar& other) = delete;
        FFBtoDBGLVar(FFBtoDBGLVar&& other) = delete;
        FFBtoDBGLVar& operator=(FFBtoDBGLVar&& other) = delete;
        ~FFBtoDBGLVar() final = default;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;
    };

} // namespace Hammer

#endif
