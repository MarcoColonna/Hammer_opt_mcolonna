///
/// @file  FFLbtoLcstar12PCR.hh
/// @brief \f$ \Lambda_b \rightarrow Lambda_c^*(2595) \f$ PCR form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_LBLCSTAR12PCR
#define HAMMER_FF_LBLCSTAR12PCR

#include "Hammer/FormFactors/PCR/FFPCRBase.hh"

namespace Hammer {

    class FFLbtoLcstar12PCR final : public FFPCRBase {

    public:

        FFLbtoLcstar12PCR();

        FFLbtoLcstar12PCR(const FFLbtoLcstar12PCR& other) = default;
        FFLbtoLcstar12PCR& operator=(const FFLbtoLcstar12PCR& other) = delete;
        FFLbtoLcstar12PCR(FFLbtoLcstar12PCR&& other) = delete;
        FFLbtoLcstar12PCR& operator=(FFLbtoLcstar12PCR&& other) = delete;
        ~FFLbtoLcstar12PCR() final = default;

        void eval(const Particle& parent, const ParticleList& daughters, const ParticleList& references) final;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;
    };

} // namespace Hammer

#endif
