///
/// @file  FFLbtoLcPCR.hh
/// @brief \f$ \Lambda_b \rightarrow Lambda_c \f$ PCR form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_LBLCPCR
#define HAMMER_FF_LBLCPCR

#include "Hammer/FormFactors/PCR/FFPCRBase.hh"

namespace Hammer {

    class FFLbtoLcPCR final : public FFPCRBase {

    public:

        FFLbtoLcPCR();

        FFLbtoLcPCR(const FFLbtoLcPCR& other) = default;
        FFLbtoLcPCR& operator=(const FFLbtoLcPCR& other) = delete;
        FFLbtoLcPCR(FFLbtoLcPCR&& other) = delete;
        FFLbtoLcPCR& operator=(FFLbtoLcPCR&& other) = delete;
        ~FFLbtoLcPCR() final = default;

        void eval(const Particle& parent, const ParticleList& daughters, const ParticleList& references) final;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;
    };

} // namespace Hammer

#endif
