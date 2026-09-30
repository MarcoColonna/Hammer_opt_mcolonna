///
/// @file  FFLbtoLcstar32LSPRVar.hh
/// @brief \f$ \Lambda_b \rightarrow Lambda_c^*(2625) \f$ LSPRVar form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_LBLCSTAR32LSPRVAR
#define HAMMER_FF_LBLCSTAR32LSPRVAR

#include "Hammer/FormFactors/LSPR/FFLSPRVarBase.hh"

namespace Hammer {

    class FFLbtoLcstar32LSPRVar final : public FFLSPRVarBase {

    public:

        FFLbtoLcstar32LSPRVar();

        FFLbtoLcstar32LSPRVar(const FFLbtoLcstar32LSPRVar& other) = default;
        FFLbtoLcstar32LSPRVar& operator=(const FFLbtoLcstar32LSPRVar& other) = delete;
        FFLbtoLcstar32LSPRVar(FFLbtoLcstar32LSPRVar&& other) = delete;
        FFLbtoLcstar32LSPRVar& operator=(FFLbtoLcstar32LSPRVar&& other) = delete;
        ~FFLbtoLcstar32LSPRVar() final = default;

        void eval(const Particle& parent, const ParticleList& daughters, const ParticleList& references) final;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;
    };

} // namespace Hammer

#endif
