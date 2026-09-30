///
/// @file  FFBSZBase.cc
/// @brief Hammer base class for EFG form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/FormFactors/BSZ/FFBSZBase.hh"

using namespace std;

namespace Hammer {

    FFBSZBase::FFBSZBase() {
        setGroup("BSZ");
    }

    void FFBSZBase::addRefs() const {
        if (!getSettingsHandler()->checkReference("Straub:2015ica")) {
            string ref =
                "@article{Straub:2015ica,\n"
                "    author = \"Bharucha, Aoife and Straub, David M. and Zwicky, Roman\",\n"
                "    title = \"{$B\to V\\ell^+\\ell^-$ in the Standard Model from light-cone sum rules}\",\n"
                "    eprint = \"1503.05534\",\n"
                "    archivePrefix = \"arXiv\",\n"
                "    primaryClass = \"hep-ph\",\n"
                "    reportNumber = \"TUM-HEP-957-14, CP3-Origins-2015-010, DIAS-2015-10\",\n"
                "    doi = \"10.1007/JHEP08(2016)098\",\n"
                "    journal = \"JHEP\",\n"
                "    volume = \"08\",\n"
                "    pages = \"098\",\n"
                "    year = \"2016\"\n"
                "}\n";
            getSettingsHandler()->addReference("Straub:2015ica", ref);
        }
    }
} // namespace Hammer
