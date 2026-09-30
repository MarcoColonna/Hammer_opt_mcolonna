///
/// @file  FFLbtoLcstar32PCR.hh
/// @brief \f$ \Lambda_b \rightarrow Lambda_c^*(2625) \f$ PCR form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_LBLCSTAR32PCR
#define HAMMER_FF_LBLCSTAR32PCR

#include "Hammer/FormFactors/PCR/FFPCRBase.hh"

namespace Hammer {

    class FFLbtoLcstar32PCR final : public FFPCRBase {

    public:

        FFLbtoLcstar32PCR();

        FFLbtoLcstar32PCR(const FFLbtoLcstar32PCR& other) = default;
        FFLbtoLcstar32PCR& operator=(const FFLbtoLcstar32PCR& other) = delete;
        FFLbtoLcstar32PCR(FFLbtoLcstar32PCR&& other) = delete;
        FFLbtoLcstar32PCR& operator=(FFLbtoLcstar32PCR&& other) = delete;
        ~FFLbtoLcstar32PCR() final = default;

        void eval(const Particle& parent, const ParticleList& daughters, const ParticleList& references) final;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;
    };

} // namespace Hammer

#endif
