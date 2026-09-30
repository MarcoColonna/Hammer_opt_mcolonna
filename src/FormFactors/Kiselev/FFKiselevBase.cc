///
/// @file  FFKiselevBase.cc
/// @brief Hammer base class for Kiselev form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/FormFactors/Kiselev/FFKiselevBase.hh"

using namespace std;

namespace Hammer {

    FFKiselevBase::FFKiselevBase() {
        setGroup("Kiselev");
    }

    void FFKiselevBase::addRefs() const {
        if (!getSettingsHandler()->checkReference("Kiselev:2002vz")) {
            string ref =
                "@article{Kiselev:2002vz,\n"
                "    author = \"Kiselev, V.V.\",\n"
                "    title = \"{Exclusive decays and lifetime of $B_c$ meson in QCD sum rules}\",\n"
                "    eprint = \"hep-ph/0211021\",\n"
                "    archivePrefix = \"arXiv\",\n"
                "    month = \"11\",\n"
                "    year = \"2002\"\n"
                "}\n";
            getSettingsHandler()->addReference("Kiselev:2002vz", ref);
        }
    }

} // namespace Hammer
