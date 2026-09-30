///
/// @file  FFBctoJpsiBGL.hh
/// @brief \f$ B_c \rightarrow J/\psi \f$ BGL form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BCJPSIBGL
#define HAMMER_FF_BCJPSIBGL

#include "Hammer/FormFactors/BGL/FFBGLBase.hh"

namespace Hammer {

    class FFBctoJpsiBGL final : public FFBGLBase {

    public:

        FFBctoJpsiBGL();

        FFBctoJpsiBGL(const FFBctoJpsiBGL& other) = default;
        FFBctoJpsiBGL& operator=(const FFBctoJpsiBGL& other) = delete;
        FFBctoJpsiBGL(FFBctoJpsiBGL&& other) = delete;
        FFBctoJpsiBGL& operator=(FFBctoJpsiBGL&& other) = delete;
        ~FFBctoJpsiBGL() final = default;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;

        void addRefs() const final;
    };

} // namespace Hammer

#endif
