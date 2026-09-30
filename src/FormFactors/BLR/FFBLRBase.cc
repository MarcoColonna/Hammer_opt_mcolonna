///
/// @file  FFBLRBase.cc
/// @brief Hammer base class for BLR form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/FormFactors/BLR/FFBLRBase.hh"

using namespace std;

namespace Hammer {

    FFBLRBase::FFBLRBase() {
        setGroup("BLR");
    }

    void FFBLRBase::addRefs() const {
        if (!getSettingsHandler()->checkReference("Bernlochner:2016bci")) {
            string ref1 =
                "@article{Bernlochner:2016bci,\n"
                "      author         = \"Bernlochner, Florian U. and Ligeti, Zoltan\",\n"
                "      title          = \"{Semileptonic $B_{(s)}$ decays to excited charmed mesons with $e,\\mu,\\tau$ "
                "and searching for new physics with $R(D^{**})$}\",\n"
                "      journal        = \"Phys. Rev.\",\n"
                "      volume         = \"D95\",\n"
                "      year           = \"2017\",\n"
                "      number         = \"1\",\n"
                "      pages          = \"014022\",\n"
                "      doi            = \"10.1103/PhysRevD.95.014022\",\n"
                "      eprint         = \"1606.09300\",\n"
                "      archivePrefix  = \"arXiv\",\n"
                "      primaryClass   = \"hep-ph\",\n"
                "      SLACcitation   = \"%%CITATION = ARXIV:1606.09300;%%\"\n"
                "}\n";
            getSettingsHandler()->addReference("Bernlochner:2016bci", ref1);
        }
        if (!getSettingsHandler()->checkReference("Bernlochner:2017jxt")) {
            string ref2 =
                "@article{Bernlochner:2017jxt,\n"
                "      author         = \"Bernlochner, Florian U. and Ligeti, Zoltan and Robinson, Dean J.\",\n"
                "      title          = \"{Model independent analysis of semileptonic $B$ decays to $D^{**}$ for "
                "arbitrary new physics}\",\n"
                "      journal        = \"Phys. Rev.\",\n"
                "      volume         = \"D97\",\n"
                "      year           = \"2018\",\n"
                "      number         = \"7\",\n"
                "      pages          = \"075011\",\n"
                "      doi            = \"10.1103/PhysRevD.97.075011\",\n"
                "      eprint         = \"1711.03110\",\n"
                "      archivePrefix  = \"arXiv\",\n"
                "      primaryClass   = \"hep-ph\",\n"
                "      SLACcitation   = \"%%CITATION = ARXIV:1711.03110;%%\"\n"
                "}\n";
            getSettingsHandler()->addReference("Bernlochner:2017jxt", ref2);
        }
    }
} // namespace Hammer
