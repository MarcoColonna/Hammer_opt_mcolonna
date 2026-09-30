///
/// @file  FFLSPRBase.cc
/// @brief Hammer base class for LSPR form factors
//

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/FormFactors/LSPR/FFLSPRBase.hh"

using namespace std;

namespace Hammer {

    FFLSPRBase::FFLSPRBase() {
        setGroup("LSPR");
    }

    void FFLSPRBase::addRefs() const {
        if (!getSettingsHandler()->checkReference("Papucci:2021pmj")) {
            string ref1 =
                "@article{Papucci:2021pmj,\n"
                "      author         = \"Papucci, Michele and Robinson, Dean J.\",\n"
                "      title          = \"{Form factor counting and HQET matching for new physics in "
                "{\\ensuremath{\\Lambda}}b{\\textrightarrow}{\\ensuremath{\\Lambda}}c*l{\\ensuremath{\\nu}}}\",\n"
                "      journal        = \"Phys. Rev. D\",\n"
                "      volume         = \"105\",\n"
                "      number         = \"1\",\n"
                "      year           = \"2022\",\n"
                "      pages          = \"016027\",\n"
                "      doi            = \"10.1103/PhysRevD.105.016027\",\n"
                "      eprint         = \"2105.09330\",\n"
                "      archivePrefix  = \"arXiv\",\n"
                "      primaryClass   = \"hep-ph\",\n"
                "      reportNumber   = \"CALT-2021-020\",\n"
                "}\n";
            getSettingsHandler()->addReference("Papucci:2021pmj", ref1);
        }
    }

    void FFLSPRBase::defineSettings() {

        // 1S scheme: mb1S = 4710 MeV, delta mb-mc = 3400 MeV, alpha_s = 26/100
        addSetting<double>("as", 0.26);
        addSetting<double>("mb1S", 4.710);  // GeV
        addSetting<double>("dmbmc", 3.400); // GeV

        addSetting<double>("s1", 0.97457);
        addSetting<double>("sp1", -1.7743);
        addSetting<double>("s11", 0.90412); // GeV
        addSetting<double>("s1p1", 0.0);    // GeV
        addSetting<double>("pc1", 0.17871); // GeV
        addSetting<double>("pb1", 0.0);     // GeV

        addSetting<double>("barLambda", 0.81);  // GeV
        addSetting<double>("barLambdap", 1.10); // GeV
    }

} // namespace Hammer
