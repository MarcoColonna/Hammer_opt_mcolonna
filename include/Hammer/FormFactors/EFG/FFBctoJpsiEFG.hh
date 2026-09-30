///
/// @file  FFBctoJpsiEFG.hh
/// @brief \f$ B_c \rightarrow J/\psi \f$ EFG form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BCJPSIEFG
#define HAMMER_FF_BCJPSIEFG

#include "Hammer/FormFactors/EFG/FFEFGBase.hh"

namespace Hammer {

    class FFBctoJpsiEFG final : public FFEFGBase {

    public:

        FFBctoJpsiEFG();

        FFBctoJpsiEFG(const FFBctoJpsiEFG& other) = default;
        FFBctoJpsiEFG& operator=(const FFBctoJpsiEFG& other) = delete;
        FFBctoJpsiEFG(FFBctoJpsiEFG&& other) = delete;
        FFBctoJpsiEFG& operator=(FFBctoJpsiEFG&& other) = delete;
        ~FFBctoJpsiEFG() final = default;

        void eval(const Particle& parent, const ParticleList& daughters, const ParticleList& references) final;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;
    };

} // namespace Hammer

#endif
