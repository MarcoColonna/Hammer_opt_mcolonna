///
/// @file  AmplBToQLepNuBase.hh
/// @brief \f$ b -> c \tau\nu \f$ base amplitude
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_AMPL_BTOQLEPNU_BASE
#define HAMMER_AMPL_BTOQLEPNU_BASE

#include "Hammer/AmplitudeBase.hh"

namespace Hammer {

    class AmplBToQLepNuBase : public AmplitudeBase {

    public:

        AmplBToQLepNuBase();

        ~AmplBToQLepNuBase() override = default;

        void preProcessWCValues(std::vector<std::complex<double>>& data, bool reverse = false) const override;

    protected:

        void defineSettings() override;

        void updateWilsonCoeffLabelPrefix() override;

    private:

        IndexList _perms;
        std::vector<double> _flips;
    };

} // namespace Hammer

#endif
