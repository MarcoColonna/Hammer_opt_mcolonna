///
/// @file  DictionaryManager.cc
/// @brief Global container class for amplitudes, rates, FFs, data
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <fstream>

#include <boost/algorithm/string.hpp>

#include "yaml-cpp/yaml.h"

#include "Hammer/Config/HammerConfig.hh"
#include "Hammer/DictionaryManager.hh"
#include "Hammer/Math/Tensor.hh"
#include "Hammer/ProcessDefinitions.hh"
#include "Hammer/SchemeDefinitions.hh"
#include "Hammer/SpecializationDefinitions.hh"
#include "Hammer/PurePhaseSpaceDefs.hh"
#include "Hammer/ProvidersRepo.hh"
#include "Hammer/ExternalData.hh"
#include "Hammer/ProcRates.hh"
#include "Hammer/Tools/HammerSerial.hh"
#include "Hammer/AmplitudeBase.hh"
#include "Hammer/FormFactorBase.hh"
#include "Hammer/Tools/Logging.hh"
#include "Hammer/Tools/Utils.hh"
#include "Hammer/Math/Histogram.hh"
#include "Hammer/Math/MultiDim/ScalarContainer.hh"

using namespace std;

namespace Hammer {

    DictionaryManager::DictionaryManager()
        : _schemeDefs{new SchemeDefinitions{}}, _providers{new ProvidersRepo{_schemeDefs.get()}},
          _specDefs{new SpecializationDefinitions{_providers.get()}},
          _external{new ExternalData{_providers.get(), _specDefs.get()}}, _rates{new ProcRates{_external.get()}},
          _purePSDefs{new PurePhaseSpaceDefs{}}, _procDefs{new ProcessDefinitions{}} {
        _providers->initialize();
    }

    void DictionaryManager::setSettingsHandler(SettingsHandler& sh) {
        SettingsConsumer::setSettingsHandler(sh);
        if (_rates) {
            _rates->setSettingsHandler(sh);
        }
        if (_providers) {
            _providers->setSettingsHandler(sh);
        }
        if (_external) {
            _external->setSettingsHandler(sh);
        }
    }

    DictionaryManager::~DictionaryManager() noexcept = default;

    const ProvidersRepo& DictionaryManager::providers() const {
        return *_providers;
    }

    const ExternalData& DictionaryManager::externalData() const {
        return *_external;
    }

    ExternalData& DictionaryManager::externalData() {
        return *_external;
    }

    const ProcRates& DictionaryManager::rates() const {
        return *_rates;
    }

    ProcRates& DictionaryManager::rates() {
        return *_rates;
    }

    const ProcessDefinitions& DictionaryManager::processDefs() const {
        return *_procDefs;
    }

    ProcessDefinitions& DictionaryManager::processDefs() {
        return *_procDefs;
    }

    PurePhaseSpaceDefs& DictionaryManager::purePSDefs() {
        return *_purePSDefs;
    }

    const PurePhaseSpaceDefs& DictionaryManager::purePSDefs() const {
        return *_purePSDefs;
    }

    SchemeDefinitions& DictionaryManager::schemeDefs() {
        return *_schemeDefs;
    }

    const SchemeDefinitions& DictionaryManager::schemeDefs() const {
        return *_schemeDefs;
    }

    SpecializationDefinitions& DictionaryManager::specializationDefs() {
        return *_specDefs;
    }

    const SpecializationDefinitions& DictionaryManager::specializationDefs() const {
        return *_specDefs;
    }


    void DictionaryManager::defineSettings() {
    }

    void DictionaryManager::init() {
        // initSettings();
        //_providers->initialize();
        _rates->init();
        _purePSDefs->init();
        _procDefs->init();
        // do not change the order
        _providers->initPreSchemeDefs();
        _schemeDefs->init(_providers->getFFGroups());
        _providers->initPostSchemeDefs();
        // auto tmp = _schemeDefs->getFFSchemeNames();
        _external->init(_schemeDefs->getFFSchemeNames());
    }

    Log& DictionaryManager::getLog() {
        return Log::getLog("Hammer.DictionaryManager");
    }

    void DictionaryManager::write(flatbuffers::FlatBufferBuilder* msgwriter) const {
        vector<flatbuffers::Offset<Serial::FBFFScheme>> schemes;
        vector<flatbuffers::Offset<Serial::FBPurePS>> pureps;
        vector<flatbuffers::Offset<Serial::FBSetting>> settings;
        vector<flatbuffers::Offset<Serial::FBSpecialization>> specs;
        flatbuffers::Offset<Serial::FBProcDefs> procs;
        getSettingsHandler()->write(msgwriter, &settings);
        _schemeDefs->write(msgwriter, &schemes);
        _purePSDefs->write(msgwriter, &pureps);
        _procDefs->write(msgwriter, &procs);
        _specDefs->write(msgwriter, &specs);
        auto serialschemes = msgwriter->CreateVector(schemes);
        auto serialps = msgwriter->CreateVector(pureps);
        auto serialsettings = msgwriter->CreateVector(settings);
        auto serialspecs = msgwriter->CreateVector(specs);
        Serial::FBHeaderBuilder serialheader{*msgwriter};
        serialheader.add_ffschemes(serialschemes);
        serialheader.add_settings(serialsettings);
        serialheader.add_pureps(serialps);
        serialheader.add_processdefs(procs);
        serialheader.add_specs(serialspecs);
        auto ver = Serial::FBVersion{HAMMER_VERSION_MAJOR, HAMMER_VERSION_MINOR, HAMMER_VERSION_PATCH};
        serialheader.add_version(&ver);
        auto headoffset = serialheader.Finish();
        msgwriter->Finish(headoffset);
    }

    bool DictionaryManager::read(const Serial::FBHeader* msgreader, bool merge) {
        bool result = true;
        result &= _procDefs->read(msgreader->processdefs(), merge);
        result &= _purePSDefs->read(msgreader, merge);
        result &= _schemeDefs->read(msgreader, merge);
        if (flatbuffers::IsFieldPresent(msgreader, Serial::FBHeader::VT_SPECS)) {
            result &= _specDefs->read(msgreader, merge);
        }
        return result;
    }

    void DictionaryManager::parseDecays(const string& data) {
        YAML::Node config = YAML::Load(data);
        processDecays(config);
    }

    void DictionaryManager::readDecays(const string& name) {
        YAML::Node config = YAML::LoadFile(name);
        processDecays(config);
    }

    void DictionaryManager::processDecays(const YAML::Node& config) {
        if (config.IsMap()) {
            if (config["Processes"]) {
                try {
                    *_procDefs = config["Processes"].as<ProcessDefinitions>();
                } catch (YAML::Exception&) {
                    MSG_ERROR("Problem parsing process input. Beware missing spaces after inline ':'.");
                }
            }
            if (config["FormFactors"]) {
                try {
                    *_schemeDefs = config["FormFactors"].as<SchemeDefinitions>();
                } catch (YAML::Exception&) {
                    MSG_ERROR("Problem parsing form factor input. Beware missing spaces after inline ':'.");
                }
            }
            if (config["PurePSVertices"]) {
                try {
                    *_purePSDefs = config["PurePSVertices"].as<PurePhaseSpaceDefs>();
                } catch (YAML::Exception&) {
                    MSG_ERROR("Problem parsing pure PS input. Beware missing spaces after inline ':'.");
                }
            }
            if (config["PurePSVertices"]) {
                try {
                    *_purePSDefs = config["PurePSVertices"].as<PurePhaseSpaceDefs>();
                } catch (YAML::Exception&) {
                    MSG_ERROR("Problem parsing pure PS input. Beware missing spaces after inline ':'.");
                }
            }
            if (config["WCSpecializations"]) {
                try {
                    _specDefs->clear();
                    config["WCSpecializations"] >> *_specDefs;
                } catch (YAML::Exception&) {
                    MSG_ERROR("Problem parsing WC specializations input. Beware missing spaces after inline ':'.");
                }
            }
        } else {
            MSG_ERROR("Problem parsing input. Beware missing spaces after inline ':'.");
        }
    }

    void DictionaryManager::saveDecays(const string& name) {
        ofstream file;
        YAML::Emitter emitter(file);
        file.open(name.c_str());
        if (file.is_open()) {
            emitter << YAML::Comment(version() + string(" header card")) << YAML::Newline;
            emitter << YAML::BeginMap;
            emitter << YAML::Key << "Processes";
            emitter << YAML::Value << *_procDefs;
            emitter << YAML::Key << "FormFactors";
            emitter << YAML::Value << *_schemeDefs;
            auto purePSNode = YAML::convert<PurePhaseSpaceDefs>::encode(*_purePSDefs);
            if (!purePSNode.IsNull()) {
                emitter << YAML::Key << "PurePSVertices";
                emitter << YAML::Value << purePSNode;
            }
            emitter << YAML::Key << "WCSpecializations";
            emitter << YAML::Value << *_specDefs;
            emitter << YAML::EndMap;
            emitter << YAML::Newline;
            file.close();
        }
    }

} // namespace Hammer
