///
/// @file  Units.cc
/// @brief Unit conversion factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/Math/Units.hh"
#include "Hammer/Exceptions.hh"

using namespace std;

namespace Hammer {

    Units* Units::_theUnits = nullptr;

    Units* Units::getUnitsInstance() {
        if (_theUnits == nullptr) {
            _theUnits = new Units();
            _theUnits->init();
        }
        return _theUnits;
    }

    Units& Units::instance() {
        return *Units::getUnitsInstance();
    }

    void Units::init() {
        _unitfactors.insert({"eV", eV});
        _unitfactors.insert({"ev", eV});
        _unitfactors.insert({"EV", eV});
        _unitfactors.insert({"keV", keV});
        _unitfactors.insert({"KeV", keV});
        _unitfactors.insert({"KEV", keV});
        _unitfactors.insert({"kev", keV});
        _unitfactors.insert({"MeV", MeV});
        _unitfactors.insert({"mev", MeV});
        _unitfactors.insert({"MEV", MeV});
        _unitfactors.insert({"GeV", GeV});
        _unitfactors.insert({"GEV", GeV});
        _unitfactors.insert({"gev", GeV});
        _unitfactors.insert({"TeV", TeV});
        _unitfactors.insert({"TEV", TeV});
        _unitfactors.insert({"tev", TeV});
    }

    double Units::getUnitsRescalingToMC(const string& mcunits, const string& localunits) const {
        auto it = _unitfactors.find(mcunits);
        auto it2 = _unitfactors.find(localunits);
        if (it != _unitfactors.end() && it2 != _unitfactors.end()) {
            return it2->second / it->second;
        }
        if (it == _unitfactors.end()) {
            throw Error("Specified " + mcunits + " units not found. Mars Climate Orbiter would like a word.");
        }
        if (it2 == _unitfactors.end()) {
            throw Error("Specified " + localunits + " units not found. Sometimes the bull wins.");
        }
        return 1.;
    }

} // namespace Hammer
