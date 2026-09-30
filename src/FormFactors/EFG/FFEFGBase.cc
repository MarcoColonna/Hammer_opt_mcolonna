///
/// @file  FFEFGBase.cc
/// @brief Hammer base class for EFG form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/FormFactors/EFG/FFEFGBase.hh"

using namespace std;

namespace Hammer {

    FFEFGBase::FFEFGBase() {
        setGroup("EFG");
    }

    void FFEFGBase::addRefs() const {
        if (!getSettingsHandler()->checkReference("Ebert:2003cn")) {
            string ref =
                "@article{Ebert:2003cn,\n"
                "    author = \"Ebert, D. and Faustov, R.N. and Galkin, V.O.\",\n"
                "    archivePrefix = \"arXiv\",\n"
                "    doi = \"10.1103/PhysRevD.68.094020\",\n"
                "    eprint = \"hep-ph/0306306\",\n"
                "    journal = \"Phys. Rev. D\",\n"
                "    pages = \"094020\",\n"
                "    reportNumber = \"HU-EP-03-035\",\n"
                "    title = \"{Weak decays of the $B_c$ meson to charmonium and $D$ mesons in the relativistic quark "
                "model}\",\n"
                "    volume = \"68\",\n"
                "    year = \"2003\"\n"
                "}\n";
            getSettingsHandler()->addReference("Ebert:2003cn", ref);
        }
    }

} // namespace Hammer
