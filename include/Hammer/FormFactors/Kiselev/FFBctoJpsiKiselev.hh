///
/// @file  FFBctoJpsiKiselev.hh
/// @brief \f$ B_c \rightarrow J/\psi \f$ Kiselev form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BCJPSIKiselev
#define HAMMER_FF_BCJPSIKiselev

#include "Hammer/FormFactors/Kiselev/FFKiselevBase.hh"

namespace Hammer {

    class FFBctoJpsiKiselev final : public FFKiselevBase {

    public:

        FFBctoJpsiKiselev();

        FFBctoJpsiKiselev(const FFBctoJpsiKiselev& other) = default;
        FFBctoJpsiKiselev& operator=(const FFBctoJpsiKiselev& other) = delete;
        FFBctoJpsiKiselev(FFBctoJpsiKiselev&& other) = delete;
        FFBctoJpsiKiselev& operator=(FFBctoJpsiKiselev&& other) = delete;
        ~FFBctoJpsiKiselev() final = default;

        void eval(const Particle& parent, const ParticleList& daughters, const ParticleList& references) final;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;
    };

} // namespace Hammer

#endif
