///
/// @file  FFLLSWBase.cc
/// @brief Hammer base class for LLSW form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/FormFactors/LLSW/FFLLSWBase.hh"

using namespace std;

namespace Hammer {

    FFLLSWBase::FFLLSWBase() {
        setGroup("LLSW");
    }

    void FFLLSWBase::addRefs() const {
        if (!getSettingsHandler()->checkReference("Leibovich:1997tu")) {
            string ref1 =
                "@article{Leibovich:1997tu,\n"
                "      author         = \"Leibovich, Adam K. and Ligeti, Zoltan and Stewart, Iain W. and Wise, Mark "
                "B.\",\n"
                "      title          = \"{Model independent results for B ---> D1(2420) lepton anti-neutrino and B "
                "---> D2* (2460) lepton anti-neutrino at order Lambda(QCD) / m(c,b)}\",\n"
                "      journal        = \"Phys. Rev. Lett.\",\n"
                "      volume         = \"78\",\n"
                "      year           = \"1997\",\n"
                "      pages          = \"3995-3998\",\n"
                "      doi            = \"10.1103/PhysRevLett.78.3995\",\n"
                "      eprint         = \"hep-ph/9703213\",\n"
                "      archivePrefix  = \"arXiv\",\n"
                "      primaryClass   = \"hep-ph\",\n"
                "      reportNumber   = \"CALT-68-2102\",\n"
                "      SLACcitation   = \"%%CITATION = HEP-PH/9703213;%%\"\n"
                "}\n";
            getSettingsHandler()->addReference("Leibovich:1997tu", ref1);
        }
        if (!getSettingsHandler()->checkReference("Leibovich:1997em")) {
            string ref2 =
                "@article{Leibovich:1997em,\n"
                "      author         = \"Leibovich, Adam K. and Ligeti, Zoltan and Stewart, Iain W. and Wise, Mark "
                "B.\",\n"
                "      title          = \"{Semileptonic B decays to excited charmed mesons}\",\n"
                "      journal        = \"Phys. Rev.\",\n"
                "      volume         = \"D57\",\n"
                "      year           = \"1998\",\n"
                "      pages          = \"308-330\",\n"
                "      doi            = \"10.1103/PhysRevD.57.308\",\n"
                "      eprint         = \"hep-ph/9705467\",\n"
                "      archivePrefix  = \"arXiv\",\n"
                "      primaryClass   = \"hep-ph\",\n"
                "      reportNumber   = \"CALT-68-2120\",\n"
                "      SLACcitation   = \"%%CITATION = HEP-PH/9705467;%%\"\n"
                "}\n";
            getSettingsHandler()->addReference("Leibovich:1997em", ref2);
        }
    }
} // namespace Hammer
