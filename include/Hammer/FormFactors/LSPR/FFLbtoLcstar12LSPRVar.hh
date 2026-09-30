///
/// @file  FFLbtoLcstar12LSPRVar.hh
/// @brief \f$ \Lambda_b \rightarrow Lambda_c^*(2595) \f$ LSPRVar form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_LBLCSTAR12LSPRVAR
#define HAMMER_FF_LBLCSTAR12LSPRVAR

#include "Hammer/FormFactors/LSPR/FFLSPRVarBase.hh"

namespace Hammer {

    class FFLbtoLcstar12LSPRVar final : public FFLSPRVarBase {

    public:

        FFLbtoLcstar12LSPRVar();

        FFLbtoLcstar12LSPRVar(const FFLbtoLcstar12LSPRVar& other) = default;
        FFLbtoLcstar12LSPRVar& operator=(const FFLbtoLcstar12LSPRVar& other) = delete;
        FFLbtoLcstar12LSPRVar(FFLbtoLcstar12LSPRVar&& other) = delete;
        FFLbtoLcstar12LSPRVar& operator=(FFLbtoLcstar12LSPRVar&& other) = delete;
        ~FFLbtoLcstar12LSPRVar() final = default;

        void eval(const Particle& parent, const ParticleList& daughters, const ParticleList& references) final;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;
    };

} // namespace Hammer

#endif
