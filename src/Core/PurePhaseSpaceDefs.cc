///
/// @file  PurePhaseSpaceDefs.cc
/// @brief Container class for pure phase space vertices definitions
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <boost/algorithm/string.hpp>

#include "yaml-cpp/yaml.h"

#include "Hammer/PurePhaseSpaceDefs.hh"
#include "Hammer/Exceptions.hh"
#include "Hammer/Tools/HammerSerial.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Tools/ParticleUtils.hh"

using namespace std;

namespace Hammer {


    Log& PurePhaseSpaceDefs::getLog() {
        return Log::getLog("Hammer.PurePhaseSpaceDefs");
    }

    void PurePhaseSpaceDefs::write(flatbuffers::FlatBufferBuilder* msgwriter,
                                   vector<flatbuffers::Offset<Serial::FBPurePS>>* purePSs) const {
        purePSs->clear();
        auto numlabels = msgwriter->CreateVectorOfStrings(_purePhaseSpaceVerticesNames.numerator);
        vector<uint64_t> ids;
        ids.reserve(_purePhaseSpaceVertices.numerator.size());
        copy(_purePhaseSpaceVertices.numerator.begin(), _purePhaseSpaceVertices.numerator.end(), back_inserter(ids));
        auto numids = msgwriter->CreateVector(ids);
        Serial::FBPurePSBuilder numBuilder{*msgwriter};
        numBuilder.add_labels(numlabels);
        numBuilder.add_ids(numids);
        purePSs->push_back(numBuilder.Finish());
        auto denlabels = msgwriter->CreateVectorOfStrings(_purePhaseSpaceVerticesNames.denominator);
        ids.clear();
        ids.reserve(_purePhaseSpaceVertices.denominator.size());
        copy(_purePhaseSpaceVertices.denominator.begin(), _purePhaseSpaceVertices.denominator.end(),
             back_inserter(ids));
        auto denids = msgwriter->CreateVector(ids);
        Serial::FBPurePSBuilder denBuilder{*msgwriter};
        denBuilder.add_labels(denlabels);
        denBuilder.add_ids(denids);
        purePSs->push_back(denBuilder.Finish());
    }

    bool PurePhaseSpaceDefs::read(const Serial::FBHeader* msgreader, bool merge) {
        if (msgreader != nullptr) {
            const auto* pureps = msgreader->pureps();
            if (pureps == nullptr || pureps->empty()) {
                return false;
            }
            ASSERT(!pureps->empty());
            const auto* ps0 = pureps->Get(0);
            if (ps0 == nullptr || ps0->labels() == nullptr || ps0->ids() == nullptr) {
                return false;
            }
            const auto* labels = ps0->labels();
            const auto* ids = ps0->ids();
            _purePhaseSpaceVerticesNames.numerator.clear();
            _purePhaseSpaceVertices.numerator.clear();
            if (!merge) {
                _purePhaseSpaceVerticesNames.numerator.reserve(labels->size());
                for (unsigned int i = 0; i < labels->size(); ++i) {
                    const auto* lstr = labels->Get(i);
                    if (lstr == nullptr) {
                        return false;
                    }
                    _purePhaseSpaceVerticesNames.numerator.emplace_back(lstr->c_str());
                }
                for (unsigned int i = 0; i < ids->size(); ++i) {
                    _purePhaseSpaceVertices.numerator.insert(ids->Get(i));
                }
            }
            const auto* ps1 = pureps->Get(1);
            if (ps1 == nullptr || ps1->labels() == nullptr || ps1->ids() == nullptr) {
                return false;
            }
            labels = ps1->labels();
            ids = ps1->ids();
            _purePhaseSpaceVerticesNames.denominator.clear();
            _purePhaseSpaceVertices.denominator.clear();
            if (!merge) {
                _purePhaseSpaceVerticesNames.denominator.reserve(labels->size());
                for (unsigned int i = 0; i < labels->size(); ++i) {
                    const auto* lstr = labels->Get(i);
                    if (lstr == nullptr) {
                        return false;
                    }
                    _purePhaseSpaceVerticesNames.denominator.emplace_back(lstr->c_str());
                }
                for (unsigned int i = 0; i < ids->size(); ++i) {
                    _purePhaseSpaceVertices.denominator.insert(ids->Get(i));
                }
            }
            return true;
        }
        return false;
    }

    void PurePhaseSpaceDefs::clearPurePhaseSpaceVertices(WTerm what) {
        switch (what) {
        case WTerm::NUMERATOR:
            _purePhaseSpaceVerticesNames.numerator.clear();
            break;
        case WTerm::DENOMINATOR:
            _purePhaseSpaceVerticesNames.denominator.clear();
            break;
        case WTerm::COMMON:
            _purePhaseSpaceVerticesNames.numerator.clear();
            _purePhaseSpaceVerticesNames.denominator.clear();
            break;
        }
    }

    void PurePhaseSpaceDefs::addPurePhaseSpaceVertices(const set<string>& decays, WTerm what) {
        switch (what) {
        case WTerm::NUMERATOR:
            _purePhaseSpaceVerticesNames.numerator.insert(_purePhaseSpaceVerticesNames.numerator.end(), decays.begin(),
                                                          decays.end());
            break;
        case WTerm::DENOMINATOR:
            _purePhaseSpaceVerticesNames.denominator.insert(_purePhaseSpaceVerticesNames.denominator.end(),
                                                            decays.begin(), decays.end());
            break;
        case WTerm::COMMON:
            _purePhaseSpaceVerticesNames.numerator.insert(_purePhaseSpaceVerticesNames.numerator.end(), decays.begin(),
                                                          decays.end());
            _purePhaseSpaceVerticesNames.denominator.insert(_purePhaseSpaceVerticesNames.denominator.end(),
                                                            decays.begin(), decays.end());
            break;
        }
    }

    void PurePhaseSpaceDefs::addPurePhaseSpaceVertex(const string& decay, WTerm what) {
        switch (what) {
        case WTerm::NUMERATOR:
            _purePhaseSpaceVerticesNames.numerator.push_back(decay);
            break;
        case WTerm::DENOMINATOR:
            _purePhaseSpaceVerticesNames.denominator.push_back(decay);
            break;
        case WTerm::COMMON:
            _purePhaseSpaceVerticesNames.numerator.push_back(decay);
            _purePhaseSpaceVerticesNames.denominator.push_back(decay);
            break;
        }
    }

    NumDenPair<set<HashId>> PurePhaseSpaceDefs::purePhaseSpaceVertices() const {
        NumDenPair<set<HashId>> result{{}, {}};
        const PID& pdg = PID::instance();
        for (const auto& elem : _purePhaseSpaceVerticesNames.numerator) {
            auto temp = pdg.expandToValidVertexUIDs(elem);
            result.numerator.insert(temp.begin(), temp.end());
        }
        for (const auto& elem : _purePhaseSpaceVerticesNames.denominator) {
            auto temp = pdg.expandToValidVertexUIDs(elem);
            result.denominator.insert(temp.begin(), temp.end());
        }
        return result;
    }

    NumDenPair<bool> PurePhaseSpaceDefs::isPurePhaseSpace(PdgId parent, const std::vector<PdgId>& daughters) const {
        auto tmp = combineDaughters(daughters, {});
        HashId id = processID(parent, tmp);
        auto itNum = _purePhaseSpaceVertices.numerator.find(id);
        auto itDen = _purePhaseSpaceVertices.denominator.find(id);
        bool foundnum = itNum != _purePhaseSpaceVertices.numerator.end();
        bool foundden = itDen != _purePhaseSpaceVertices.denominator.end();
        return NumDenPair<bool>{foundnum, foundden};
    }

    void PurePhaseSpaceDefs::init() {
        _purePhaseSpaceVertices = purePhaseSpaceVertices();
    }

    YAML::Emitter& operator<<(YAML::Emitter& out, const PurePhaseSpaceDefs& s) {
        out << YAML::convert<PurePhaseSpaceDefs>::encode(s);
        return out;
    }

} // namespace Hammer

namespace YAML {

    Node convert<::Hammer::PurePhaseSpaceDefs>::encode(const ::Hammer::PurePhaseSpaceDefs& value) {
        Node node;
        if (!value._purePhaseSpaceVerticesNames.numerator.empty()) {
            Node tmpNode;
            for (const auto& elem : value._purePhaseSpaceVerticesNames.numerator) {
                tmpNode.push_back(elem);
            }
            node["Numerator"] = tmpNode;
        }
        if (!value._purePhaseSpaceVerticesNames.denominator.empty()) {
            Node tmpNode;
            for (const auto& elem : value._purePhaseSpaceVerticesNames.denominator) {
                tmpNode.push_back(elem);
            }
            node["Denominator"] = tmpNode;
        }
        return node;
    }

    bool convert<::Hammer::PurePhaseSpaceDefs>::decode(const Node& node, ::Hammer::PurePhaseSpaceDefs& value) {
        if (node.IsMap()) {
            for (auto entry2 : node) {
                auto what = entry2.first.as<string>();
                if (what == "Numerator") {
                    value._purePhaseSpaceVerticesNames.numerator.clear();
                    YAML::Node processes = entry2.second;
                    if (processes.IsSequence()) {
                        for (auto process : processes) {
                            value._purePhaseSpaceVerticesNames.numerator.push_back(process.as<string>());
                        }
                    } else {
                        return false;
                    }
                } else if (what == "Denominator") {
                    value._purePhaseSpaceVerticesNames.denominator.clear();
                    YAML::Node processes = entry2.second;
                    if (processes.IsSequence()) {
                        for (auto process : processes) {
                            value._purePhaseSpaceVerticesNames.denominator.push_back(process.as<string>());
                        }
                    } else {
                        return false;
                    }
                } else {
                    return false;
                }
            }
        } else {
            return false;
        }
        return true;
    }
} // namespace YAML
