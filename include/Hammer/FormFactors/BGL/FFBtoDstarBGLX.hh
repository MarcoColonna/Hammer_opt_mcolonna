///
/// @file  FFBtoDstarBGLX.hh
/// @brief \f$ B \rightarrow D \f$ BGL Extended form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BDSTARBGLX
#define HAMMER_FF_BDSTARBGLX

#include "Hammer/FormFactors/BGL/FFBGLXBase.hh"

namespace Hammer {

    class FFBtoDstarBGLX final : public FFBGLXBase {

    public:

        FFBtoDstarBGLX();

        FFBtoDstarBGLX(const FFBtoDstarBGLX& other) = default;
        FFBtoDstarBGLX& operator=(const FFBtoDstarBGLX& other) = delete;
        FFBtoDstarBGLX(FFBtoDstarBGLX&& other) = delete;
        FFBtoDstarBGLX& operator=(FFBtoDstarBGLX&& other) = delete;

        ~FFBtoDstarBGLX() final = default;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;
    };

} // namespace Hammer

#endif
