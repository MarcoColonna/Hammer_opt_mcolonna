///
/// @file  FFCLNBase.cc
/// @brief Hammer base class for CLN form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/FormFactors/CLN/FFCLNBase.hh"

using namespace std;

namespace Hammer {

    FFCLNBase::FFCLNBase() {
        setGroup("CLN");
    }

    void FFCLNBase::addRefs() const {
        if (!getSettingsHandler()->checkReference("Caprini:1997mu")) {
            string ref =
                "@article{Caprini:1997mu,\n"
                "      author         = \"Caprini, Irinel and Lellouch, Laurent and Neubert, Matthias\",\n"
                "      title          = \"{Dispersive bounds on the shape of anti-B --> D(*) lepton anti-neutrino "
                "form-factors}\",\n"
                "      journal        = \"Nucl. Phys.\",\n"
                "      volume         = \"B530\",\n"
                "      year           = \"1998\",\n"
                "      pages          = \"153-181\",\n"
                "      doi            = \"10.1016/S0550-3213(98)00350-2\",\n"
                "      eprint         = \"hep-ph/9712417\",\n"
                "      archivePrefix  = \"arXiv\",\n"
                "      primaryClass   = \"hep-ph\",\n"
                "      reportNumber   = \"CERN-TH-97-091, CPT-97-P3480\",\n"
                "      SLACcitation   = \"%%CITATION = HEP-PH/9712417;%%\"\n"
                "}\n";
            getSettingsHandler()->addReference("Caprini:1997mu", ref);
        }
    }

} // namespace Hammer
