///
/// @file  FFLbtoLcBLRS.hh
/// @brief \f$ \Lambda_b \rightarrow Lambda_c \f$ BLRS form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_LBLCBLRS
#define HAMMER_FF_LBLCBLRS

#include "Hammer/FormFactors/BLRS/FFBLRSBase.hh"

namespace Hammer {

    class FFLbtoLcBLRS final : public FFBLRSBase {

    public:

        FFLbtoLcBLRS();

        FFLbtoLcBLRS(const FFLbtoLcBLRS& other) = default;
        FFLbtoLcBLRS& operator=(const FFLbtoLcBLRS& other) = delete;
        FFLbtoLcBLRS(FFLbtoLcBLRS&& other) = delete;
        FFLbtoLcBLRS& operator=(FFLbtoLcBLRS&& other) = delete;
        ~FFLbtoLcBLRS() final = default;

        void eval(const Particle& parent, const ParticleList& daughters, const ParticleList& references) final;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;
    };

} // namespace Hammer

#endif
