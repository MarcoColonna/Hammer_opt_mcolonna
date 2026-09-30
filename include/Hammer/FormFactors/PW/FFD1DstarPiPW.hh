///
/// @file  FFD1DstarPiPW.hh
/// @brief \f$ D_1 \rightarrow D^* \pi \f$ partial wave coefficients
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_D1DSTARPIPW
#define HAMMER_FF_D1DSTARPIPW

#include "Hammer/FormFactors/PW/FFPWBase.hh"

namespace Hammer {

    class FFD1DstarPiPW final : public FFPWBase {

    public:

        FFD1DstarPiPW();

        FFD1DstarPiPW(const FFD1DstarPiPW& other) = default;
        FFD1DstarPiPW& operator=(const FFD1DstarPiPW& other) = delete;
        FFD1DstarPiPW(FFD1DstarPiPW&& other) = delete;
        FFD1DstarPiPW& operator=(FFD1DstarPiPW&& other) = delete;
        ~FFD1DstarPiPW() final = default;

        void eval(const Particle& parent, const ParticleList& daughters, const ParticleList& references) final;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;
    };

} // namespace Hammer

#endif
