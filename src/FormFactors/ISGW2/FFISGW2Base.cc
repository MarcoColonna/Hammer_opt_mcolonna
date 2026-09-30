///
/// @file  FFISGW2Base.cc
/// @brief Hammer base class for ISGW2 form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/FormFactors/ISGW2/FFISGW2Base.hh"
#include "Hammer/Math/Constants.hh"

using namespace std;

namespace Hammer {

    FFISGW2Base::FFISGW2Base() {
        setGroup("ISGW2");
    }

    double FFISGW2Base::GetGammaji(double z) {
        double value;
        value = 2 + ((2.0 * z) / (1 - z)) * log(z);
        value *= -1.0;
        return value;
    }

    double FFISGW2Base::Getas(double mq1, double mq2) {
        constexpr double lambdaSq = 0.04;
        double Nf = 4;
        if (mq1 < 1.85) {
            Nf = 3.0;
        }
        double value = 12.0 * pi;
        value /= (33.0 - 2.0 * Nf);
        value /= log(mq2 * mq2 / lambdaSq);
        if (mq2 <= 0.6) {
            value = 0.6;
        }
        return value;
    }

    double FFISGW2Base::Getas(double m) {
        return Getas(m, m);
    }

    void FFISGW2Base::addRefs() const {
        if (!getSettingsHandler()->checkReference("Isgur:1988gb")) {
            string ref1 =
                "@article{Isgur:1988gb,\n"
                "      author         = \"Isgur, Nathan and Scora, Daryl and Grinstein, Benjamin and Wise, Mark B.\",\n"
                "      title          = \"{Semileptonic B and D Decays in the Quark Model}\",\n"
                "      journal        = \"Phys. Rev.\",\n"
                "      volume         = \"D39\",\n"
                "      year           = \"1989\",\n"
                "      pages          = \"799-818\",\n"
                "      doi            = \"10.1103/PhysRevD.39.799\",\n"
                "      reportNumber   = \"UTPT-88-12\",\n"
                "      SLACcitation   = \"%%CITATION = PHRVA,D39,799;%%\"\n"
                "}\n";
            getSettingsHandler()->addReference("Isgur:1988gb", ref1);
        }
        if (!getSettingsHandler()->checkReference("Scora:1995ty")) {
            string ref2 =
                "@article{Scora:1995ty,\n"
                "      author         = \"Scora, Daryl and Isgur, Nathan\",\n"
                "      title          = \"{Semileptonic meson decays in the quark model: An update}\",\n"
                "      journal        = \"Phys. Rev.\",\n"
                "      volume         = \"D52\",\n"
                "      year           = \"1995\",\n"
                "      pages          = \"2783-2812\",\n"
                "      doi            = \"10.1103/PhysRevD.52.2783\",\n"
                "      eprint         = \"hep-ph/9503486\",\n"
                "      archivePrefix  = \"arXiv\",\n"
                "      primaryClass   = \"hep-ph\",\n"
                "      reportNumber   = \"CEBAF-TH-94-14\",\n"
                "      SLACcitation   = \"%%CITATION = HEP-PH/9503486;%%\"\n"
                "}\n";
            getSettingsHandler()->addReference("Scora:1995ty", ref2);
        }
    }
} // namespace Hammer
