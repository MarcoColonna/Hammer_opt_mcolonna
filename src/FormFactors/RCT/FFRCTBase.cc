///
/// @file  FFRCTBase.cc
/// @brief Hammer base class for RCT form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/FormFactors/RCT/FFRCTBase.hh"

using namespace std;

namespace Hammer {

    FFRCTBase::FFRCTBase() {
        setGroup("RCT");
    }

    void FFRCTBase::addRefs() const {
        if (!getSettingsHandler()->checkReference("Shekhovtsova:2012ra")) {
            string ref1 =
                "@article{Shekhovtsova:2012ra,\n"
                "      author         = \"Shekhovtsova, O. and Przedzinski, T. and Roig, P. and Was, Z.\",\n"
                "      title          = \"{Resonance chiral Lagrangian currents and $\\tau$ decay Monte Carlo}\",\n"
                "      journal        = \"Phys. Rev.\",\n"
                "      volume         = \"D86\",\n"
                "      year           = \"2012\",\n"
                "      pages          = \"113008\",\n"
                "      doi            = \"10.1103/PhysRevD.86.113008\",\n"
                "      eprint         = \"1203.3955\",\n"
                "      archivePrefix  = \"arXiv\",\n"
                "      primaryClass   = \"hep-ph\",\n"
                "      reportNumber   = \"IFJPAN-IV-2011-6, UAB-FT-695, FTUV-2011-09-29, IFIC-11-53, "
                "CERN-PH-TH-2012-016\",\n"
                "      SLACcitation   = \"%%CITATION = ARXIV:1203.3955;%%\"\n"
                "}\n";
            getSettingsHandler()->addReference("Shekhovtsova:2012ra", ref1);
        }
        if (!getSettingsHandler()->checkReference("Nugent:2013hxa")) {
            string ref2 =
                "@article{Nugent:2013hxa,\n"
                "      author         = \"Nugent, I. M. and Przedzinski, T. and Roig, P. and Shekhovtsova, O. and Was, "
                "Z.\",\n"
                "      title          = \"{Resonance chiral Lagrangian currents and experimental data for "
                "$\\tau^-\\to\\pi^{-}\\pi^{-}\\pi^{+}\\nu_{\\tau}$}\",\n"
                "      journal        = \"Phys. Rev.\",\n"
                "      volume         = \"D88\",\n"
                "      year           = \"2013\",\n"
                "      pages          = \"093012\",\n"
                "      doi            = \"10.1103/PhysRevD.88.093012\",\n"
                "      eprint         = \"1310.1053\",\n"
                "      archivePrefix  = \"arXiv\",\n"
                "      primaryClass   = \"hep-ph\",\n"
                "      reportNumber   = \"IFJPAN-2013-5, UAB-FT-731\",\n"
                "      SLACcitation   = \"%%CITATION = ARXIV:1310.1053;%%\"\n"
                "}\n";
            getSettingsHandler()->addReference("Nugent:2013hxa", ref2);
        }
    }

} // namespace Hammer
