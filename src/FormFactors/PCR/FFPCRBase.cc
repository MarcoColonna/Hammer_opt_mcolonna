///
/// @file  FFPCRBase.cc
/// @brief Hammer base class for PCR form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/FormFactors/PCR/FFPCRBase.hh"

using namespace std;

namespace Hammer {

    FFPCRBase::FFPCRBase() {
        setGroup("PCR");
    }

    void FFPCRBase::addRefs() const {
        if (!getSettingsHandler()->checkReference("Pervin:2005ve")) {
            string ref =
                "@article{Pervin:2005ve,\n"
                "      author         = \"Pervin, Muslema and Roberts, Winston and Capstick, Simon\",\n"
                "      title          = \"{Semileptonic decays of heavy lambda baryons in a quark model}\",\n"
                "      journal        = \"Phys. Rev.\",\n"
                "      volume         = \"C72\",\n"
                "      year           = \"2005\",\n"
                "      pages          = \"035201\",\n"
                "      doi            = \"10.1103/PhysRevC.72.035201\",\n"
                "      eprint         = \"nucl-th/0503030\",\n"
                "      archivePrefix  = \"arXiv\",\n"
                "      primaryClass   = \"nucl-th\",\n"
                "      reportNumber   = \"JLAB-THY-05-304\",\n"
                "      SLACcitation   = \"%%CITATION = NUCL-TH/0503030;%%\"\n"
                "}\n";
            getSettingsHandler()->addReference("Pervin:2005ve", ref);
        }
    }

} // namespace Hammer
