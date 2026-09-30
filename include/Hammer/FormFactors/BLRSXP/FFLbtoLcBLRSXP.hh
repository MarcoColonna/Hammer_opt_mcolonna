///
/// @file  FFLbtoLcBLRSXP.hh
/// @brief \f$ \Lambda_b \rightarrow Lambda_c \f$ BLRSXP form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_LBLCBLRSXP
#define HAMMER_FF_LBLCBLRSXP

#include "Hammer/FormFactors/BLRSXP/FFBLRSXPBase.hh"

namespace Hammer {

    class FFLbtoLcBLRSXP final : public FFBLRSXPBase {

    public:

        FFLbtoLcBLRSXP();

        FFLbtoLcBLRSXP(const FFLbtoLcBLRSXP& other) = default;
        FFLbtoLcBLRSXP& operator=(const FFLbtoLcBLRSXP& other) = delete;
        FFLbtoLcBLRSXP(FFLbtoLcBLRSXP&& other) = delete;
        FFLbtoLcBLRSXP& operator=(FFLbtoLcBLRSXP&& other) = delete;
        ~FFLbtoLcBLRSXP() final = default;

        void eval(const Particle& parent, const ParticleList& daughters, const ParticleList& references) final;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;
    };

} // namespace Hammer

#endif
