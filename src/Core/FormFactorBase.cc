///
/// @file  FormFactorBase.cc
/// @brief Hammer base form factor class
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/FormFactorBase.hh"
#include "Hammer/Tools/Utils.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/MultiDim/VectorContainer.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"
#include "Hammer/Math/MultiDim/ScalarContainer.hh"
#include "Hammer/Math/Units.hh"
#include "Hammer/Particle.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FormFactorBase::FormFactorBase() : _errPrefixGroup{"", ""}, _mFFErrLabel{NONE} {
    }

    FormFactorBase::FormFactorBase(const FormFactorBase& other)
        : ParticleData(other), SettingsConsumer(other), _errPrefixGroup(other._errPrefixGroup),
          _mFFErrLabel(other._mFFErrLabel), _mFFInfoList{other._mFFInfoList}, _mFFErrNames(other._mFFErrNames),
          _units{other._units} {
        for (const auto& elem : other._tensorList) {
            addTensor(Tensor{elem});
        }
    }

    // FormFactorBase& FormFactorBase::operator=(const FormFactorBase& other) {
    //     ParticleData::operator=(other);
    //     SettingsConsumer::operator=(other);
    //     _errPrefixGroup = other._errPrefixGroup;
    //     _mFFErrLabel = other._mFFErrLabel;
    //     _mFFErrNames = other._mFFErrNames;
    //     _tensorList.clear();
    //     for (auto& elem : other._tensorList) {
    //         addTensor(Tensor{elem});
    //     }
    //     return *this;
    // }

    // FormFactorBase::~FormFactorBase() {
    // }

    vector<double> FormFactorBase::getErrVectorFromDict(const map<std::string, double>& errDict) const {
        vector<double> result(_mFFErrNames->size() + 1);
        result[0] = 1.;
        //        for (size_t pos = 0ul; pos < _mFFErrNames->size(); ++pos) {
        //            auto it = errDict.find(_mFFErrNames->at(pos));
        //            if (it != errDict.end()) {
        //                result[pos + 1] = it->second;
        //            }
        //        }
        for (auto it = errDict.begin(); it != errDict.end(); ++it) {
            auto itf = find(_mFFErrNames->begin(), _mFFErrNames->end(), it->first);
            if (itf != _mFFErrNames->end()) {
                auto pos = distance(_mFFErrNames->begin(), itf);
                result[static_cast<size_t>(pos) + 1ul] = it->second;
            } else {
                MSG_WARNING("'" + it->first + "' does not belong to the parametrization " +
                            getFFErrPrefixGroup().get() + " set of parameters. Ignoring.");
            }
        }
        return result;
    }

    vector<double> FormFactorBase::getErrVectorFromSettings(bool useDefault) const {
        vector<double> result(_mFFErrNames->size() + 1);
        result[0] = 1.;
        for (size_t pos = 0ul; pos < _mFFErrNames->size(); ++pos) {
            if (!useDefault) {
                result[pos + 1] = *getSetting<double>(getFFErrPrefixGroup().get(), _mFFErrNames->at(pos));
            }
        }
        return result;
    }

    void FormFactorBase::updateFFErrSettings(const vector<double>& values) {
        this->updateVectorOfSettings(values, *_mFFErrNames, getFFErrPrefixGroup().get(), WTerm::COMMON);
    }

    void FormFactorBase::updateFFErrSettings(const map<std::string, double>& values) {
        this->updateVectorOfSettings(values, getFFErrPrefixGroup().get(), WTerm::COMMON);
    }

    map<string, double> FormFactorBase::retrieveFFErrSettings() const {
        return this->retrieveVectorOfSettings<double>(*_mFFErrNames, getFFErrPrefixGroup().get(), WTerm::COMMON);
    }

    void FormFactorBase::updateFFErrTensor(const vector<double>& values, MD::SharedTensorData& data) const {
        ASSERT(values.size() == _mFFErrNames->size() + 1);
        if (!data || data->rank() == 0) {
            data = MD::SharedTensorData{
                MD::makeEmptyVector({static_cast<uint16_t>(values.size())}, {_mFFInfoList[_signatureIndex].second})
                    .release()};
        }
        for (IndexType i = 0; i < static_cast<IndexType>(values.size()); ++i) {
            data->element({i}) = values[i];
        }
    }


    const FFPrefixGroup& FormFactorBase::getFFErrPrefixGroup() const {
        return _mFFInfoList[_signatureIndex].first;
    }

    std::pair<FFPrefixGroup, IndexLabel> FormFactorBase::getFFErrInfo() const {
        return _mFFInfoList[_signatureIndex];
    }

    void FormFactorBase::init() {
        // initSettings();
        ///@todo anything else?
    }

    Tensor& FormFactorBase::getTensor() {
        return _tensorList[_signatureIndex];
    }

    const Tensor& FormFactorBase::getTensor() const {
        return _tensorList[_signatureIndex];
    }

    void FormFactorBase::setGroup(const string& name) {
        bool reInitSettings = (_errPrefixGroup.group != name && !_errPrefixGroup.group.empty());
        _errPrefixGroup.group = name;
        for (auto& elem : _mFFInfoList) {
            elem.first.group = name;
        }
        if (reInitSettings && getSettingsHandler() != nullptr) {
            initSettings();
        }
    }

    void FormFactorBase::setPrefix(const string& name) {
        _errPrefixGroup.prefix = name;
    }

    const std::string& FormFactorBase::group() const {
        return _mFFInfoList[_signatureIndex].first.group;
    }

    void FormFactorBase::setUnits(const string& name) {
        addSetting<string>("Units", name);
    }

    void FormFactorBase::calcUnits() {
        string mcunits;
        auto* mc = getSetting<string>("Hammer", "Units");
        if (mc != nullptr) {
            mcunits = *mc;
        } else {
            throw Error("Hammer units not found. Mars Climate Orbiter would like a word.");
        }
        Units& units = Units::instance();
        _units = units.getUnitsRescalingToMC(mcunits, *(getSetting<string>("Units")));
    }

    Tensor FormFactorBase::getFFPSIntegrand(const EvaluationGrid& intPoints) {
        if (!intPoints.empty()) {
            auto newdims = getTensor().dims();
            auto newlabs = getTensor().labels();
            auto oldsize = newdims.size();
            newdims.reserve(2 * oldsize);
            newlabs.reserve(2 * oldsize);
            std::copy_n(newdims.begin(), oldsize, std::back_inserter(newdims));
            transform_n(newlabs.begin(), oldsize, std::back_inserter(newlabs),
                        [](IndexLabel l) -> IndexLabel { return static_cast<IndexLabel>(-l); });
            newlabs.push_back(INTEGRATION_INDEX);
            newdims.push_back(static_cast<IndexType>(intPoints.size()));
            Tensor result{"", MD::makeEmptySparse(newdims, newlabs)};
            for (IndexType i = 0; i < static_cast<IndexType>(intPoints.size()); ++i) {
                evalAtPSPoint(intPoints[i]);
                Tensor t = getTensor();
                t.outerSquare();
                result.addAt(t, INTEGRATION_INDEX, i);
            }
            return result;
        }
        evalAtPSPoint({});
        Tensor t = getTensor();
        t.outerSquare();
        return t;
    }

    Log& FormFactorBase::getLog() {
        return Log::getLog("Hammer.FormFactorBase");
    }

    void FormFactorBase::addTensor(Tensor&& tensor) {
        _tensorList.push_back(std::move(tensor));
    }

    void FormFactorBase::addProcessSignature(PdgId parent, const vector<PdgId>& daughters) {
        ParticleData::addProcessSignature(parent, daughters);
        _mFFInfoList.emplace_back(_errPrefixGroup, _mFFErrLabel);
    }

    void FormFactorBase::defineAndAddErrSettings(const vector<string>& names) {
        auto* preExistingNames = getSetting<vector<string>>("ErrNames");
        if (preExistingNames != nullptr) {
            // pre-existing: may happen for duplicates
            if (names.size() < preExistingNames->size()) {
                MSG_WARNING("FF Error Names renaming: too many names provided. Excess will be ignored.");
            }
            vector<string> newNames = *preExistingNames;
            newNames.resize(names.size());
            for (size_t i = 0; i < newNames.size(); ++i) {
                if (newNames[i].empty()) {
                    newNames[i] = names[i];
                }
            }
            addSetting<vector<string>>("ErrNames", names); // updates the defaults
            for (const auto& elem : names) {               // create the settings for the defaults
                if (getSetting<double>(elem) == nullptr) {
                    addSetting<double>(elem, 0.);
                }
            }
            for (const auto& elem : newNames) { // create the settings for the current ones
                if (getSetting<double>(elem) == nullptr) {
                    addSetting<double>(elem, 0.);
                }
            }
            *preExistingNames = newNames; // updates the values
        } else {
            addSetting<vector<string>>("ErrNames", names);
        }
        bindErrNames();
    }

    void FormFactorBase::bindErrNames() {
        _mFFErrNames = getSetting<vector<string>>("ErrNames");
        for (const auto& elem : *_mFFErrNames) {
            addSetting<double>(elem, 0.);
        }
    }

    void FF1to1Base::eval(const Particle& parent, const ParticleList& daughters, const ParticleList& /*references*/) {

        // Momenta
        const FourMomentum& pParent = parent.momentum();
        const FourMomentum& pDaughter = daughters[0].momentum();
        // const FourMomentum& pTau = daughters[2].momentum();


        // kinematic objects
        const double Mparent = pParent.mass();
        const double Mdaughter = pDaughter.mass();
        // const double Mt = pTau.mass();
        const double Sqq = pow(Mparent, 2.) + pow(Mdaughter, 2.) - 2. * (pParent * pDaughter);

        evalAtPSPoint({Sqq}, {Mparent, Mdaughter});
    }

    tuple<double, double, double> FF1to1Base::getParentDaughterHadMasses(const vector<double>& masses) const {
        double Mparent = 0.;
        double Mdaughter = 0.;
        double unitres = 1.;
        if (masses.size() >= 2) {
            Mparent = masses[0];
            Mdaughter = masses[1];
            unitres = _units;
        } else {
            Mparent = this->masses()[0];
            Mdaughter = this->masses()[1];
        }
        return {Mparent, Mdaughter, unitres};
    }

    double FF1to1Base::getW(double Sqq, double Mparent, double Mdaughter) {
        return (Mparent * Mparent + Mdaughter * Mdaughter - Sqq) / (2. * Mparent * Mdaughter);
    }

} // namespace Hammer
