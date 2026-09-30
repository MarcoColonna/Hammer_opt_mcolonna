///
/// @file  FFBGLBase.hh
/// @brief Hammer base class for BGL form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BGL_BASE
#define HAMMER_FF_BGL_BASE

#include "Hammer/FormFactorBase.hh"

namespace Hammer {

    /// @brief Base class for BGL form factors
    ///
    /// @ingroup FormFactors
    class FFBGLBase : public FF1to1Base {

    public:

        FFBGLBase();

        FFBGLBase(const FFBGLBase& other) = default;
        FFBGLBase& operator=(const FFBGLBase& other) = delete;
        FFBGLBase(FFBGLBase&& other) = delete;
        FFBGLBase& operator=(FFBGLBase&& other) = delete;
        ~FFBGLBase() override = default;

    protected:

        /// @brief
        /// @param[in] point
        /// @param[in] masses
        /// @return
        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) override = 0;

        /// @brief
        void defineSettings() override = 0;

        void addRefs() const override;
    };

} // namespace Hammer

#endif
