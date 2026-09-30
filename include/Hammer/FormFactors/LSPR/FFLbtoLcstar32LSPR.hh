///
/// @file  FFLbtoLcstar32LSPR.hh
/// @brief \f$ \Lambda_b \rightarrow Lambda_c^*(2625) \f$ LSPR form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_LBLCSTAR32LSPR
#define HAMMER_FF_LBLCSTAR32LSPR

#include "Hammer/FormFactors/LSPR/FFLSPRBase.hh"

namespace Hammer {

    class FFLbtoLcstar32LSPR final : public FFLSPRBase {

    public:

        FFLbtoLcstar32LSPR();

        FFLbtoLcstar32LSPR(const FFLbtoLcstar32LSPR& other) = default;
        FFLbtoLcstar32LSPR& operator=(const FFLbtoLcstar32LSPR& other) = delete;
        FFLbtoLcstar32LSPR(FFLbtoLcstar32LSPR&& other) = delete;
        FFLbtoLcstar32LSPR& operator=(FFLbtoLcstar32LSPR&& other) = delete;
        ~FFLbtoLcstar32LSPR() final = default;

        void eval(const Particle& parent, const ParticleList& daughters, const ParticleList& references) final;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;
    };

} // namespace Hammer

#endif
