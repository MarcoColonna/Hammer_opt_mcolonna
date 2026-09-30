///
/// @file  FFBtoDBGL.hh
/// @brief \f$ B \rightarrow D \f$ BGL form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BDBGL
#define HAMMER_FF_BDBGL

#include "Hammer/FormFactors/BGL/FFBGLBase.hh"

namespace Hammer {

    class FFBtoDBGL final : public FFBGLBase {

    public:

        FFBtoDBGL();

        FFBtoDBGL(const FFBtoDBGL& other) = default;
        FFBtoDBGL& operator=(const FFBtoDBGL& other) = delete;
        FFBtoDBGL(FFBtoDBGL&& other) = delete;
        FFBtoDBGL& operator=(FFBtoDBGL&& other) = delete;
        ~FFBtoDBGL() final = default;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;
    };

} // namespace Hammer

#endif
