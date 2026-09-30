///
/// @file  FFBCLBase.cc
/// @brief Hammer base class for EFG form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/FormFactors/BCL/FFBCLBase.hh"

using namespace std;

namespace Hammer {

    FFBCLBase::FFBCLBase() {
        setGroup("BCL");
    }

    void FFBCLBase::addRefs() const {
        if (!getSettingsHandler()->checkReference("Bourrely:2008za")) {
            string ref =
                "@article{Bourrely:2008za,\n"
                "author = \"Bourrely, Claude and Caprini, Irinel and Lellouch, Laurent\",\n"
                "title = \"{Model-independent description of B ---> pi l nu decays and a determination of |V(ub)|}\",\n"
                "eprint = \"0807.2722\",\n"
                "archivePrefix = \"arXiv\",\n"
                "primaryClass = \"hep-ph\",\n"
                "reportNumber = \"CPT-P36-2007\",\n"
                "doi = \"10.1103/PhysRevD.82.099902\",\n"
                "journal = \"Phys. Rev. D\",\n"
                "volume = \"79\",\n"
                "pages = \"013008\",\n"
                "year = \"2009\",\n"
                "note = \"[Erratum: Phys.Rev.D 82, 099902 (2010)]\"\n"
                "}\n";
            getSettingsHandler()->addReference("Bourrely:2008za", ref);
        }
    }
} // namespace Hammer
