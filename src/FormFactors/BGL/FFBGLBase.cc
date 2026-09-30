///
/// @file  FFBGLBase.hh
/// @brief Hammer base class for BGL form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/FormFactors/BGL/FFBGLBase.hh"

using namespace std;

namespace Hammer {

    FFBGLBase::FFBGLBase() {
        setGroup("BGL");
    }

    void FFBGLBase::addRefs() const {
        if (!getSettingsHandler()->checkReference("Boyd:1995sq")) {
            string ref1 =
                "@article{Boyd:1995sq,\n"
                "      author         = \"Boyd, C. Glenn and Grinstein, Benjamin and Lebed, Richard F.\",\n"
                "      title          = \"{Model independent determinations of anti-B ---> D (lepton), D* (lepton) "
                "anti-neutrino form-factors}\",\n"
                "      journal        = \"Nucl. Phys.\",\n"
                "      volume         = \"B461\",\n"
                "      year           = \"1996\",\n"
                "      pages          = \"493-511\",\n"
                "      doi            = \"10.1016/0550-3213(95)00653-2\",\n"
                "      eprint         = \"hep-ph/9508211\",\n"
                "      archivePrefix  = \"arXiv\",\n"
                "      primaryClass   = \"hep-ph\",\n"
                "      reportNumber   = \"UCSD-PTH-95-11\",\n"
                "      SLACcitation   = \"%%CITATION = HEP-PH/9508211;%%\"\n"
                "}\n";
            getSettingsHandler()->addReference("Boyd:1995sq", ref1);
        }
        if (!getSettingsHandler()->checkReference("Boyd:1997kz")) {
            string ref2 =
                "@article{Boyd:1997kz,\n"
                "      author         = \"Boyd, C. Glenn and Grinstein, Benjamin and Lebed, Richard F.\",\n"
                "      title          = \"{Precision corrections to dispersive bounds on form-factors}\",\n"
                "      journal        = \"Phys. Rev.\",\n"
                "      volume         = \"D56\",\n"
                "      year           = \"1997\",\n"
                "      pages          = \"6895-6911\",\n"
                "      doi            = \"10.1103/PhysRevD.56.6895\",\n"
                "      eprint         = \"hep-ph/9705252\",\n"
                "      archivePrefix  = \"arXiv\",\n"
                "      primaryClass   = \"hep-ph\",\n"
                "      reportNumber   = \"CMU-HEP-97-07A, UCSD-PTH-97-12\",\n"
                "      SLACcitation   = \"%%CITATION = HEP-PH/9705252;%%\"\n"
                "}\n";
            getSettingsHandler()->addReference("Boyd:1997kz", ref2);
        }
    }

} // namespace Hammer
