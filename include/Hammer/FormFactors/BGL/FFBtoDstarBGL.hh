///
/// @file  FFBtoDstarBGL.hh
/// @brief \f$ B \rightarrow D \f$ BGL form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BDSTARBGL
#define HAMMER_FF_BDSTARBGL

#include "Hammer/FormFactors/BGL/FFBGLBase.hh"

namespace Hammer {

    class FFBtoDstarBGL final : public FFBGLBase {

    public:

        FFBtoDstarBGL();

        FFBtoDstarBGL(const FFBtoDstarBGL& other) = default;
        FFBtoDstarBGL& operator=(const FFBtoDstarBGL& other) = delete;
        FFBtoDstarBGL(FFBtoDstarBGL&& other) = delete;
        FFBtoDstarBGL& operator=(FFBtoDstarBGL&& other) = delete;

        ~FFBtoDstarBGL() final = default;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;
    };

} // namespace Hammer

#endif
