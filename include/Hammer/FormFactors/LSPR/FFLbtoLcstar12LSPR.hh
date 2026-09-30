///
/// @file  FFLbtoLcstar12LSPR.hh
/// @brief \f$ \Lambda_b \rightarrow Lambda_c^*(2595) \f$ LSPR form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_LBLCSTAR12LSPR
#define HAMMER_FF_LBLCSTAR12LSPR

#include "Hammer/FormFactors/LSPR/FFLSPRBase.hh"

namespace Hammer {

    class FFLbtoLcstar12LSPR final : public FFLSPRBase {

    public:

        FFLbtoLcstar12LSPR();

        FFLbtoLcstar12LSPR(const FFLbtoLcstar12LSPR& other) = default;
        FFLbtoLcstar12LSPR& operator=(const FFLbtoLcstar12LSPR& other) = delete;
        FFLbtoLcstar12LSPR(FFLbtoLcstar12LSPR&& other) = delete;
        FFLbtoLcstar12LSPR& operator=(FFLbtoLcstar12LSPR&& other) = delete;
        ~FFLbtoLcstar12LSPR() final = default;

        void eval(const Particle& parent, const ParticleList& daughters, const ParticleList& references) final;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;
    };

} // namespace Hammer

#endif
