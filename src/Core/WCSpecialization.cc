///
/// @file  SpecializationDef.cc
/// @brief Hammer base amplitude class
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <iterator>

#include "Hammer/WCSpecialization.hh"
#include "Hammer/AmplitudeBase.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/MultiDim/VectorContainer.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"
#include "Hammer/Math/MultiDim/Operations.hh"
#include "Hammer/Tools/HammerSerial.hh"
#include "Hammer/Tools/HammerYaml.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    WCSpecialization::WCSpecialization() : _pad{0ul}, _baseLabel{NONE}, _fullLabel{NONE}, _baseAmpl{nullptr} {
    }

    WCSpecialization::WCSpecialization(const string& prefix, IndexLabel baseLabel, const string& id,
                                       const vector<string>& names)
        : _prefixId{prefix, id}, _pad{0ul}, _coordNames{names}, _baseLabel{baseLabel}, _fullLabel{NONE},
          _baseAmpl{nullptr} {
        if (!_coordNames.empty()) {
            _pad = generateSpecializedPad(id);
            _fullLabel = specializeLabel(baseLabel, _pad);
        }
    }

    WCSpecialization::WCSpecialization(const Serial::FBSpecialization* msgreader)
        : _pad{0ul}, _baseLabel{NONE}, _fullLabel{NONE}, _baseAmpl{nullptr} {
        read(msgreader);
    }

    void WCSpecialization::defineSettings() {
        if (!_coordNames.empty()) {
            string optionPath = _prefixId.get();
            SettingsHandler* settings = getSettingsHandler();
            for (const auto& elem : _coordNames) {
                settings->addSetting<complex<double>>(optionPath, elem, 0.);
            }
        }
    }

    void WCSpecialization::initialize() {
        ASSERT(_baseAmpl != nullptr);
        // set initial origin to current WC settings
        MD::SharedTensorData defaultorigin{};
        _baseAmpl->updateWCTensor(_baseAmpl->getWCVectorFromSettings(WTerm::NUMERATOR), defaultorigin);
        if (isPartialSpecialization()) {
            // initialize rank-2 projector of appropriate size
            vector<map<string, complex<double>>> dummysubspace(_coordNames.size());
            _baseAmpl->createWCProjectionTensor(dummysubspace, getFullLabel(), defaultorigin, _projector);
            ASSERT(_projector->rank() == 2u);
        } else {
            _projector = defaultorigin;
            ASSERT(_projector->rank() == 1u);
        }
    }

    vector<complex<double>>
    WCSpecialization::getWCSpecVectorFromDict(const map<string, complex<double>>& subDict) const {
        if (_coordNames.empty()) {
            MSG_ERROR("Specialization " + _prefixId.get() + " is not a partial specialization!");
            return {};
        }
        vector<complex<double>> result(_coordNames.size() + 1);
        result[0] = 1.;
        //        for (size_t pos = 0ul; pos < _coordNames.size(); ++pos) {
        //            auto it = subDict.find(_coordNames[pos]);
        //            if (it != subDict.end()) {
        //                result[pos + 1ul] = it->second;
        //            }
        //        }
        for (auto it = subDict.begin(); it != subDict.end(); ++it) {
            auto itc = find(_coordNames.begin(), _coordNames.end(), it->first);
            if (itc != _coordNames.end()) {
                auto pos = distance(_coordNames.begin(), itc);
                result[static_cast<size_t>(pos + 1)] = it->second;
            } else {
                MSG_WARNING("'" + it->first + "' does not belong to the specialization " + _prefixId.get() +
                            " set of coordinates. Ignoring.");
            }
        }
        return result;
    }

    vector<complex<double>> WCSpecialization::getWCSpecVectorFromSettings() const {
        if (_coordNames.empty()) {
            MSG_ERROR("Specialization " + _prefixId.get() + " is not a partial specialization!");
            return {};
        }
        WTerm oldwhat = const_cast<WCSpecialization*>(this)->setWeightTerm(WTerm::NUMERATOR);
        vector<complex<double>> result(_coordNames.size() + 1);
        result[0] = 1.;
        for (size_t pos = 0ul; pos < _coordNames.size(); ++pos) {
            result[pos + 1] = *getSetting<complex<double>>(_prefixId.get(), _coordNames[pos]);
        }
        const_cast<WCSpecialization*>(this)->setWeightTerm(oldwhat);
        return result;
    }

    void WCSpecialization::updateWCSpecSettings(const vector<complex<double>>& values) {
        if (_coordNames.empty()) {
            MSG_ERROR("Specialization " + _prefixId.get() + " is not a partial specialization!");
            return;
        }
        this->updateVectorOfSettings(values, _coordNames, _prefixId.get(), WTerm::NUMERATOR);
    }

    void WCSpecialization::updateWCSpecSettings(const map<string, complex<double>>& values) {
        if (_coordNames.empty()) {
            MSG_ERROR("Specialization " + _prefixId.get() + " is not a partial specialization!");
            return;
        }
        this->updateVectorOfSettings(values, _prefixId.get(), WTerm::NUMERATOR);
    }

    vector<complex<double>> WCSpecialization::resetWCSpecialization() {
        vector<complex<double>> result(_coordNames.size() + 1);
        result[0] = 1.;
        if (_coordNames.empty()) {
            MSG_ERROR("Specialization " + _prefixId.get() + " is not a partial specialization!");
        } else {
            for (const auto& _coordName : _coordNames) {
                getSettingsHandler()->resetSetting(_prefixId.get(), _coordName, WTerm::NUMERATOR);
            }
        }
        return result;
    }

    map<string, complex<double>> WCSpecialization::retrieveWCSpecialization() const {
        if (_coordNames.empty()) {
            MSG_ERROR("Specialization " + _prefixId.get() + " is not a partial specialization!");
            return {};
        }
        return this->retrieveVectorOfSettings<complex<double>>(_coordNames, _prefixId.get(), WTerm::NUMERATOR);
    }

    void WCSpecialization::updateWCSpecTensor(vector<complex<double>> values, MD::SharedTensorData& data) const {
        ASSERT(values.size() == _coordNames.size() + 1);
        if (!data || data->rank() == 0) {
            data = MD::SharedTensorData{
                MD::makeVector({static_cast<uint16_t>(values.size())}, {_fullLabel}, values).release()};
        } else {
            for (IndexType i = 0; i < static_cast<IndexType>(values.size()); ++i) {
                data->element({i}) = values[i];
            }
        }
    }

    bool WCSpecialization::isPartialSpecialization() const {
        return !_coordNames.empty();
    }

    const SpecPrefixId& WCSpecialization::getPrefixId() const {
        return _prefixId;
    }

    IndexLabel WCSpecialization::getFullLabel() const {
        return _fullLabel;
    }

    IndexLabel WCSpecialization::getBaseLabel() const {
        return _baseLabel;
    }

    const vector<string>& WCSpecialization::getCoordinates() const {
        return _coordNames;
    }

    uint32_t WCSpecialization::getPad() const {
        return _pad;
    }

    MD::SharedTensorData WCSpecialization::getProjectionTensor() const {
        return _projector;
    }

    void WCSpecialization::setSubspaceOrigin(const vector<complex<double>>& values) {
        if (_baseAmpl != nullptr) {
            //_baseAmpl->updateWCSettings(values, WTerm::NUMERATOR);
            _baseAmpl->updateWCProjectionVector(-1, values, _projector);
        } else {
            MSG_ERROR("Wilson coefficients " + _prefixId.prefix + " or specialization subspace with id " +
                      _prefixId.id + " not found.");
        }
    }

    void WCSpecialization::setSubspaceOrigin(const map<string, complex<double>>& settings) {
        if (_baseAmpl != nullptr) {
            auto values = _baseAmpl->getWCVectorFromDict(settings);
            setSubspaceOrigin(values);
        } else {
            MSG_ERROR("Wilson coefficients " + _prefixId.prefix + " or specialization subspace with id " +
                      _prefixId.id + " not found.");
        }
    }

    void WCSpecialization::setSubspaceBasis(const vector<map<string, complex<double>>>& subspace) {
        if (_baseAmpl != nullptr) {
            if (_coordNames.empty()) {
                MSG_ERROR(
                    "Wilson coefficients " + _prefixId.prefix +
                    " cannot be partially specialized. Subspace coordinates not specified! Alexandretta! Of course!");
                return;
            }
            MD::SharedTensorData data{};
            _baseAmpl->createWCProjectionTensor(subspace, getFullLabel(), _projector, data);
            ASSERT(data->rank() == 2u);
            _projector = data;
            // MSG_INFO(_prefixId.prefix + " Wilson coefficients will be restricted to specialized subspace in
            // weights.");
        } else {
            MSG_ERROR("Wilson coefficients " + _prefixId.prefix + " or specialization subspace with id " +
                      _prefixId.id + " not found.");
        }
    }

    size_t WCSpecialization::coordIndex(const string& name) const {
        auto it = find(_coordNames.begin(), _coordNames.end(), name);
        return static_cast<size_t>(distance(_coordNames.begin(), it));
    }

    void WCSpecialization::setSubspaceVector(const string& coord, const vector<complex<double>>& values) {
        if (_baseAmpl != nullptr) {
            if (!_coordNames.empty()) {
                size_t index = coordIndex(coord);
                if (index < _coordNames.size()) {
                    _baseAmpl->updateWCProjectionVector(static_cast<int>(index), values, _projector);
                } else {
                    MSG_ERROR("Coordinate " + coord + "not found for specialization " + _prefixId.get() +
                              "! Alexandretta! Of course!");
                }
            } else {
                MSG_ERROR("Cannot set a subspace direction for " + _prefixId.prefix +
                          " which is a full specialization!");
            }
        } else {
            MSG_ERROR("Wilson coefficients " + _prefixId.prefix + " or specialization subspace with id " +
                      _prefixId.id + " not found.");
        }
    }

    void WCSpecialization::setSubspaceVector(const string& coord, const map<string, complex<double>>& values) {
        if (_baseAmpl != nullptr) {
            if (!_coordNames.empty()) {
                size_t index = coordIndex(coord);
                if (index < _coordNames.size()) {
                    _baseAmpl->updateWCProjectionVector(static_cast<int>(index), values, _projector);
                } else {
                    MSG_ERROR("Coordinate " + coord + "not found for specialization " + _prefixId.get() +
                              "! Alexandretta! Of course!");
                }
            } else {
                MSG_ERROR("Cannot set a subspace direction for " + _prefixId.prefix +
                          " which is a full specialization!");
            }
        } else {
            MSG_ERROR("Wilson coefficients " + _prefixId.prefix + " or specialization subspace with id " +
                      _prefixId.id + " not found.");
        }
    }

    void WCSpecialization::setAmplitude(AmplitudeBase* amp) {
        _baseAmpl = amp;
        if (_baseAmpl != nullptr) {
            setSettingsHandler(*_baseAmpl);
        }
    }


    void WCSpecialization::write(flatbuffers::FlatBufferBuilder* msgwriter,
                                 flatbuffers::Offset<Serial::FBSpecialization>* msg) const {
        auto serialprefix = msgwriter->CreateString(_prefixId.prefix);
        auto serialid = msgwriter->CreateString(_prefixId.id);
        auto serialcoords = msgwriter->CreateVectorOfStrings(_coordNames);
        auto serialproj = flatbuffers::Offset<Serial::FBSingleTensor>{(_projector->write(msgwriter)).first.o};
        Serial::FBSpecializationBuilder serialSD{*msgwriter};
        serialSD.add_prefix(serialprefix);
        serialSD.add_id(serialid);
        serialSD.add_coords(serialcoords);
        serialSD.add_projector(serialproj);
        *msg = serialSD.Finish();
    }

    void WCSpecialization::read(const Serial::FBSpecialization* msgreader) {
        if (msgreader != nullptr) {
            if (msgreader->prefix() == nullptr || msgreader->id() == nullptr) {
                MSG_ERROR("Missing prefix or id on read. You need a bell and ring it in the morning!");
                return;
            }
            _prefixId.prefix = msgreader->prefix()->c_str();
            _prefixId.id = msgreader->id()->c_str();
            if (msgreader->coords() == nullptr) {
                MSG_ERROR("Missing coords on read. You need a song and sing it in the morning!");
                return;
            }
            const auto* coords = msgreader->coords();
            _coordNames.clear();
            _coordNames.reserve(coords->size());
            for (unsigned int j = 0; j < coords->size(); ++j) {
                const auto* coordStr = coords->Get(j);
                if (coordStr == nullptr) {
                    continue;
                }
                _coordNames.emplace_back(coordStr->c_str());
            }
            if (msgreader->projector() == nullptr) {
                MSG_ERROR("Missing projector on read. Have fun storming the castle!");
                return;
            }
            const auto* data = msgreader->projector();
            if (data->sparse()) {
                _projector.reset(static_cast<MD::IContainer*>(new MD::SparseContainer{data}));
            } else {
                _projector.reset(static_cast<MD::IContainer*>(new MD::VectorContainer{data}));
            }
            _baseLabel = _projector->labels()[0];
            if (_projector->labels().size() == 2) {
                _fullLabel = _projector->labels()[1];
                _pad = getSpecializedLabelComponents(_fullLabel).first;
            } else {
                _fullLabel = NONE;
                _pad = 0ul;
            }
        }
    }

    Log& WCSpecialization::getLog() {
        return Log::getLog("Hammer.WCSpecialization");
    }

    bool operator==(const WCSpecialization& a, const WCSpecialization& b) {
        bool result = (a.getPrefixId().id == b.getPrefixId().id);
        if (result) {
            result &= a.getBaseLabel() == b.getBaseLabel();
        }
        if (result) {
            result &= a.getCoordinates() == b.getCoordinates();
        }
        if (result) {
            auto pa = a.getProjectionTensor();
            auto pb = b.getProjectionTensor();
            if (pa) {
                if (pb) {
                    result &= (*pa) == (*pb);
                } else {
                    result = false;
                }
            } else if (pb) {
                result = false;
            }
        }
        return result;
    }

    bool operator!=(const WCSpecialization& a, const WCSpecialization& b) {
        return !(a == b);
    }

} // namespace Hammer

namespace YAML {

    Node convert<::Hammer::WCSpecialization>::encode(const ::Hammer::WCSpecialization& value) {
        Node names;
        for (const auto& elem : value._coordNames) {
            names.push_back(elem);
        }
        names.SetStyle(YAML::EmitterStyle::Flow);
        Node node;
        node["Coords"] = names;
        bool isPartial = !value._coordNames.empty();
        auto numwcs = value._projector->dims()[0];
        for (uint16_t i = 0; i < static_cast<uint16_t>(value._coordNames.size() + 1ul); ++i) {
            Node tmp;
            vector<complex<double>> tmpvec;
            tmpvec.reserve(numwcs);
            for (uint16_t j = 0; j < numwcs; ++j) {
                tmpvec.push_back(isPartial ? value._projector->element({j, i}) : value._projector->element({j}));
            }
            value._baseAmpl->preProcessWCValues(tmpvec, true); // undo reshuffling
            // TODO: vector of dictionary?
            for (uint16_t j = 0; j < numwcs; ++j) {
                Node tmp2;
                tmp2.push_back(tmpvec[j].real());
                tmp2.push_back(tmpvec[j].imag());
                tmp2.SetStyle(YAML::EmitterStyle::Flow);
                tmp.push_back(tmp2);
            }
            tmp.SetStyle(YAML::EmitterStyle::Flow);
            node[(i == 0 ? "Origin" : value._coordNames[i - 1])] = tmp;
        }
        return node;
    }

    bool convert<::Hammer::WCSpecialization>::decode(const Node& node, ::Hammer::WCSpecialization& value) {
        if (!node.IsMap() || value._baseLabel == 0ul || value._prefixId.id.empty() || value._prefixId.prefix.empty()) {
            return false;
        }
        value._coordNames = node["Coords"].as<vector<string>>();
        const auto& curNode = node["Origin"];
        if (curNode.IsSequence()) {
            vector<complex<double>> data{};
            data.reserve(curNode.size());
            auto ito = data.begin();
            for (YAML::const_iterator itvv = curNode.begin(); itvv != curNode.end(); ++itvv, ++ito) {
                if (!itvv->IsSequence() || itvv->size() != 2) {
                    return false;
                }
                *ito = itvv->as<complex<double>>();
            }
            value.setSubspaceOrigin(data);
        } else if (curNode.IsMap()) {
            map<string, complex<double>> data;
            for (YAML::const_iterator itvv = curNode.begin(); itvv != curNode.end(); ++itvv) {
                if (!itvv->second.IsSequence() || itvv->second.size() != 2) {
                    return false;
                }
                data.insert({itvv->first.as<string>(), itvv->second.as<complex<double>>()});
            }
            value.setSubspaceOrigin(data);
        }
        for (const auto& elem : value._coordNames) {
            const auto& curNode2 = node[elem];
            if (curNode2.IsSequence()) {
                vector<complex<double>> data{};
                data.reserve(curNode2.size());
                auto ito = data.begin();
                for (YAML::const_iterator itvv = curNode2.begin(); itvv != curNode2.end(); ++itvv, ++ito) {
                    if (!itvv->IsSequence() || itvv->size() != 2) {
                        return false;
                    }
                    *ito = itvv->as<complex<double>>();
                }
                value.setSubspaceVector(elem, data);
            } else if (curNode2.IsMap()) {
                map<string, complex<double>> data;
                for (YAML::const_iterator itvv = curNode2.begin(); itvv != curNode2.end(); ++itvv) {
                    if (!itvv->second.IsSequence() || itvv->second.size() != 2) {
                        return false;
                    }
                    data.insert({itvv->first.as<string>(), itvv->second.as<complex<double>>()});
                }
                value.setSubspaceVector(elem, data);
            }
        }
        return true;
    }


} // namespace YAML
