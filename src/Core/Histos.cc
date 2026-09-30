///
/// @file  Histos.cc
/// @brief Hammer histogram manager
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <functional>

#include "Hammer/Histos.hh"
#include "Hammer/Math/Histogram.hh"
#include "Hammer/DictionaryManager.hh"
#include "Hammer/SchemeDefinitions.hh"
#include "Hammer/ExternalData.hh"
#include "Hammer/Math/MultiDim/BinnedIndexing.hh"
#include "Hammer/Tools/HammerRoot.hh"
#include "Hammer/Tools/HammerSerial.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Tools/Logging.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    Histos::Histos(DictionaryManager* dict) : _dictionaries{dict} {
    }

    Histos::~Histos() noexcept {
        clear();
        _dictionaries = nullptr;
    }

    void Histos::addHistogramDefinition(const string& name, const IndexList& binSizes, bool hasUnderOverFlow,
                                        const MD::BinRangeList& ranges, bool shouldCompress, bool withErrors) {
        if (name == "Total Sum of Weights") { // now defunct as histo is no longer added in initRun
            _histogramDefs.insert(
                {name, HistogramDefinition{hasUnderOverFlow, binSizes, ranges, shouldCompress, withErrors}});

        } else {
            auto it = _histogramDefs.find(name);
            if (it != _histogramDefs.end()) {
                MSG_ERROR(string("A histogram with name '") + name + string("' is already present! Not adding"));
                return;
            }
            auto resdef = _histogramDefs.emplace(
                name, HistogramDefinition{hasUnderOverFlow, binSizes, ranges, shouldCompress, withErrors});
            if (!resdef.second) {
                MSG_ERROR(string("Invalid histogram definition for histogram '") + name +
                          string("'! Not adding. The Sesame Street Count is on his way."));
                return;
            }
        }
        addSetting<bool>("CalcGeneralHistogram_" + name, true);
        auto ith = _histograms.emplace(name, SpecializationDict<SchemeDict<HistogramSet>>{}).first;
        if (_initialized) {
            auto specs = getSpecializationsForHisto(name);
            if (isOn("CalcGeneralHistogram_" + name)) {
                specs.insert(Spec::none());
            }
            SpecializationDict<SchemeDict<HistogramSet>>::iterator ith2;
            for (const auto& elem : specs) {
                ith2 = ith->second.emplace(elem, SchemeDict<HistogramSet>{}).first;
                for (const auto& elem2 : _dictionaries->schemeDefs().getFFSchemeNames()) {
                    ith2->second.emplace(elem2, HistogramSet{});
                }
            }
        }
    }

    void Histos::addHistogramDefinition(const string& name, const MD::BinEdgeList& binEdges, bool hasUnderOverFlow,
                                        bool shouldCompress, bool withErrors) {
        auto it = _histogramDefs.find(name);
        if (it != _histogramDefs.end()) {
            MSG_ERROR(string("A histogram with name '") + name + string("' is already present! Not adding"));
            return;
        }
        auto resdef =
            _histogramDefs.emplace(name, HistogramDefinition{hasUnderOverFlow, binEdges, shouldCompress, withErrors});
        if (!resdef.second) {
            MSG_ERROR(string("Invalid histogram definition for histogram '") + name +
                      string("'! Not adding. The Sesame Street Count is on his way."));
            return;
        }
        addSetting<bool>("CalcGeneralHistogram_" + name, true);
        auto ith = _histograms.emplace(name, SpecializationDict<SchemeDict<HistogramSet>>{}).first;
        if (_initialized) {

            auto specs = getSpecializationsForHisto(name);
            if (isOn("CalcGeneralHistogram_" + name)) {
                specs.insert(Spec::none());
            }
            SpecializationDict<SchemeDict<HistogramSet>>::iterator ith2;
            for (const auto& elem : specs) {
                ith2 = ith->second.emplace(elem, SchemeDict<HistogramSet>{}).first;
                for (const auto& elem2 : _dictionaries->schemeDefs().getFFSchemeNames()) {
                    ith2->second.emplace(elem2, HistogramSet{});
                }
            }
        }
    }

    void Histos::createProjectedHistogram(const string& oldName, const string& newName,
                                          const set<uint16_t>& collapsedIndexPositions) {
        auto itNew = _histogramDefs.find(newName);
        auto itOld = _histogramDefs.find(oldName);
        if (itNew != _histogramDefs.end()) {
            MSG_ERROR(string("Histogram '") + newName + string("' already exists. Not projecting"));
            return;
        }
        if (itOld == _histogramDefs.end()) {
            MSG_ERROR(string("Cannot find histogram '") + oldName + string("'. Not projecting"));
            return;
        }
        MD::BinEdgeList newBinEdges = itOld->second.getIndexing().edges();
        for (auto itPos = collapsedIndexPositions.crbegin(); itPos != collapsedIndexPositions.crend(); ++itPos) {
            newBinEdges.erase(newBinEdges.begin() + *itPos);
        }
        if (newBinEdges.empty()) {
            MSG_ERROR(string("You are trying to project all the dimensions of histogram '") + oldName +
                      string("'! Mascalzone!"));
            return;
        }
        auto resdef = _histogramDefs.emplace(
            newName, HistogramDefinition{itOld->second.getIndexing().hasUnderOverFlow(), newBinEdges,
                                         itOld->second.shouldCompress(), itOld->second.shouldKeepErrors()});
        if (!resdef.second) {
            MSG_ERROR(string("Invalid histogram definition for histogram '") + newName +
                      string("'! Not adding. The Sesame Street Count is on his way."));
            return;
        }
        auto resSpec = _histograms.emplace(newName, SpecializationDict<SchemeDict<HistogramSet>>{});
        auto itHistoOld = _histograms.find(oldName);
        for (auto& elemSpec : itHistoOld->second) {
            auto resScheme = resSpec.first->second.emplace(elemSpec.first, SchemeDict<HistogramSet>{});
            for (auto& elemScheme : elemSpec.second) {
                resScheme.first->second.emplace(elemScheme.first,
                                                HistogramSet{elemScheme.second, resdef.first->first,
                                                             resdef.first->second, collapsedIndexPositions});
            }
        }
        for (const auto& elem : getSpecializationsForHisto(oldName)) {
            addSpecializationForHisto(newName, elem);
        }
    }


    void Histos::removeHistogram(const string& name) {
        auto it = _histogramDefs.find(name);
        if (it == _histogramDefs.end()) {
            MSG_ERROR(string("Cannot find histogram '") + name + string("'. Not removing"));
            return;
        }
        _histogramDefs.erase(it);
        _histograms.erase(name);
        _usedSpecializations.erase(name);
    }


    void Histos::setHistogramCompression(const std::string& name, bool value) {
        auto it = _histogramDefs.find(name);
        if (it == _histogramDefs.end()) {
            MSG_ERROR(string("Cannot find histogram '") + name + string("'. Not setting compression"));
            return;
        }
        it->second.setCompression(value);
    }

    void Histos::setHistogramKeepErrors(const std::string& name, bool value) {
        auto it = _histogramDefs.find(name);
        if (it == _histogramDefs.end()) {
            MSG_ERROR(string("Cannot find histogram '") + name + string("'. Not setting compression"));
            return;
        }
        it->second.setKeepErrors(value);
    }

    void Histos::addHistogramFixedData(const std::string& name, const MD::SharedTensorData& data) {
        auto it = _histogramDefs.find(name);
        if (it == _histogramDefs.end()) {
            MSG_ERROR(string("Cannot find histogram '") + name + string("'. Not adding fixed data."));
            return;
        }
        it->second.addFixedData(data);
    }

    void Histos::resetHistogramFixedData(const std::string& name) {
        auto it = _histogramDefs.find(name);
        if (it == _histogramDefs.end()) {
            MSG_ERROR(string("Cannot find histogram '") + name + string("'. Not resetting fixed data."));
            return;
        }
        it->second.resetFixedData();
    }

    void Histos::init() {
        // initSettings();
        for (auto& elem : _histograms) {
            elem.second.clear();
        }
        auto schemes = _dictionaries->schemeDefs().getFFSchemeNames();
        for (auto& elem : _histograms) {
            auto specs = getSpecializationsForHisto(elem.first);
            if (specs.empty()) {
                addSpecializationForHisto(elem.first, Spec::none());
            }
            if (!isOn("CalcGeneralHistogram_" + elem.first)) {
                removeSpecializationForHisto(elem.first, Spec::none(), false);
                if (getSpecializationsForHisto(elem.first).empty()) {
                    MSG_ERROR("CalcGeneralHistogram set to false with no WC specializations declared in histogram '" +
                              elem.first + "'. The hive is out of balance.");
                }
            } else {
                addSpecializationForHisto(elem.first, Spec::none());
            }
            specs = getSpecializationsForHisto(elem.first);
            for (const auto& specName : specs) {
                auto it2 = elem.second.emplace(specName, SchemeDict<HistogramSet>{}).first;
                for (const auto& schemeName : schemes) {
                    bool compress = _histogramDefs.find(elem.first)->second.shouldCompress();
                    it2->second.emplace(schemeName, HistogramSet{compress});
                }
            }
        }
        for (auto& elem : _histogramDefs) {
            elem.second.setKeepErrors(elem.second.shouldKeepErrors() || *getSetting<bool>("KeepErrors"));
        }
        _initialized = true;
    }

    void Histos::setEventId(EventUID eventId) {
        _currentEventId = std::move(eventId);
    }

    void Histos::clear() {
        _histograms.clear();
        _histogramDefs.clear();
        _currentEventId.clear();
        _initialized = false;
    }

    IndexList Histos::getBinIndices(const string& name, const MD::BinValue& value) const {
        auto it = _histogramDefs.find(name);
        if (it != _histogramDefs.end()) {
            return it->second.getIndexing().valueToPos(value);
        }
        throw Error("Histogram name not found. You can't handle the truth!");
    }

    size_t Histos::size() const {
        return _histogramDefs.size();
    }

    bool Histos::canFill(const string& name, const vector<double>& values) const {
        auto it = _histogramDefs.find(name);
        if (it != _histogramDefs.end()) {
            return it->second.getIndexing().isValid(values);
        }
        MSG_ERROR(string("Can't find name ") + string(name) + string(" in histogram request."));
        return false;
    }

    void Histos::checkExists(const string& name) const {
        auto it = _histograms.find(name);
        if (it == _histograms.end()) {
            MSG_ERROR(string("The histogram '") + string(name) + string("' does not exist. The Dude does not abide."));
        }
    }

    const HistogramSet* Histos::getEntry(const string& name, const string& specialization, const string& scheme) const {
        auto it = _histograms.find(name);
        if (it != _histograms.end()) {
            auto it2 = it->second.find(specialization);
            if (it2 != it->second.end()) {
                auto it3 = it2->second.find(scheme);
                if (it3 != it2->second.end()) {
                    return &it3->second;
                }
                MSG_ERROR(
                    string("Can't find scheme ") + string(scheme) +
                    string(" in histogram request. Returning empty histogram list. Hey Stella!!! (Check addHistogram "
                           "is applied before initRun)"));
                return nullptr;
            }
            MSG_ERROR(string("Can't find specialization ") + string(specialization) +
                      string(" in histogram request. Returning empty histogram list. Hey Stella!!! (Check addHistogram "
                             "is applied before initRun)"));
            return nullptr;
        }
        MSG_ERROR(string("Can't find name ") + string(name) +
                  string(" in histogram request. Returning empty histogram list. Hey Stella!!! (Check addHistogram is "
                         "applied before initRun)"));
        return nullptr;
    }

    HistogramSet* Histos::getEntry(const string& name, const string& specialization, const string& scheme) {
        auto it = _histograms.find(name);
        if (it != _histograms.end()) {
            auto it2 = it->second.find(specialization);
            if (it2 != it->second.end()) {
                auto it3 = it2->second.find(scheme);
                if (it3 != it2->second.end()) {
                    return &it3->second;
                }
                MSG_ERROR(
                    string("Can't find scheme ") + string(scheme) +
                    string(" in histogram request. Returning empty histogram list. Hey Stella!!! (Check addHistogram "
                           "is applied before initRun)"));
                return nullptr;
            }
            MSG_ERROR(string("Can't find specialization ") + string(specialization) +
                      string(" in histogram request. Returning empty histogram list. Hey Stella!!! (Check addHistogram "
                             "is applied before initRun)"));
            return nullptr;
        }
        MSG_ERROR(string("Can't find name ") + string(name) +
                  string(" in histogram request. Returning empty histogram list. Hey Stella!!! (Check addHistogram is "
                         "applied before initRun)"));
        return nullptr;
    }

    EventUIDGroup Histos::getHistogramEventIds(const string& name, const string& specialization,
                                               const string& scheme) const {
        const auto* entry = getEntry(name, specialization, scheme);
        if (entry != nullptr) {
            return entry->getEventIdsInHistogram();
        }
        MSG_ERROR("Wrong name/scheme in histogram request. Returning empty eventID list.");
        return EventUIDGroup{};
    }


    MD::BinEdgeList Histos::getHistogramEdges(const string& name) const {
        auto it = _histogramDefs.find(name);
        if (it != _histogramDefs.end()) {
            return it->second.getIndexing().edges();
        }
        MSG_ERROR(string("Can't find name ") + string(name) +
                  string(" in histogram request. Returning empty bin edges list. Hey Stella!!! (Check addHistogram is "
                         "applied before initRun)"));
        return {};
    }

    IndexList Histos::getHistogramShape(const string& name) const {
        auto it = _histogramDefs.find(name);
        if (it != _histogramDefs.end()) {
            return it->second.getIndexing().dims();
        }
        MSG_ERROR(string("Can't find name ") + string(name) +
                  string(" in histogram request. Returning empty bin edges list. Hey Stella!!! (Check addHistogram is "
                         "applied before initRun)"));
        return {};
    }

    bool Histos::getUnderOverFlows(const string& name) const {
        auto it = _histogramDefs.find(name);
        if (it != _histogramDefs.end()) {
            return it->second.getIndexing().hasUnderOverFlow();
        }
        MSG_ERROR(string("Can't find name ") + string(name) + string(" in histogram request."));
        return false;
    }

    vector<Tensor> Histos::getExternalData(const string& scheme, vector<LabelsList> labels,
                                           const string& specialization) const {
        vector<Tensor> result;
        result.reserve(labels.size());
        transform(labels.begin(), labels.end(), back_inserter(result),
                  [this, &scheme, &specialization](const LabelsList& val) -> const Tensor& {
                      return _dictionaries->externalData().getExternalVectors(scheme, val, specialization);
                  });
        return result;
    }

    EventIdGroupDict<IOHistogram> Histos::getHistograms(const string& name, const string& scheme,
                                                        const string& specialization) const {
        const auto* entry = getEntry(name, specialization, scheme);
        if (entry != nullptr) {
            auto externalData = getExternalData(scheme, entry->getHistogramLabels(), specialization);
            return entry->specializeEventHistograms(externalData);
        }
        MSG_ERROR("Wrong name/scheme in histogram request. Returning empty histogram map");
        return EventIdGroupDict<IOHistogram>{};
    }

    IOHistogram Histos::getHistogram(const string& name, const string& scheme, const string& specialization) const {
        const auto* entry = getEntry(name, specialization, scheme);
        if (entry != nullptr) {
            auto externalData = getExternalData(scheme, entry->getHistogramLabels(), specialization);
            return entry->specializeSumHistogram(externalData);
        }
        MSG_ERROR("Wrong name/scheme in histogram request. Returning empty histogram.");
        return IOHistogram{};
    }

    vector<string> Histos::getHistogramNames() const {
        vector<string> result;
        transform(_histogramDefs.begin(), _histogramDefs.end(), back_inserter(result),
                  [](const pair<string, HistogramDefinition>& val) -> string { return val.first; });
        return result;
    }

    vector<EventUID> Histos::getEventIDRepsForHisto(const string& name, const string& specialization,
                                                    const string& scheme) const {
        const auto* entry = getEntry(name, specialization, scheme);
        if (entry != nullptr) {
            return entry->getEventUIDRepresentatives();
        }
        return {};
    }

    bool Histos::haveSpecializations(const string& name) const {
        return !(_usedSpecializations.find(name) == _usedSpecializations.end());
    }

    set<string> Histos::getSpecializationsForHisto(const string& name) const {
        return getOrDefault(_usedSpecializations, name, {});
    }

    void Histos::addSpecializationForHisto(const string& histoName, const string& specName) {
        _usedSpecializations[histoName].insert(specName);
    }

    void Histos::removeSpecializationForHisto(const string& histoName, const string& specName, bool removeData) {
        auto it = _usedSpecializations.find(histoName);
        if (it != _usedSpecializations.end()) {
            it->second.erase(specName);
            if (it->second.empty()) {
                _usedSpecializations.erase(histoName);
            }
        }
        if (removeData) {
            auto it2 = _histograms.find(histoName);
            if (it2 != _histograms.end()) {
                it2->second.erase(specName);
            }
        }
    }

    void Histos::removeAllSpecializationsForHisto(const string& histoName, bool removeData) {
        _usedSpecializations.erase(histoName);
        if (removeData) {
            auto it = _histograms.find(histoName);
            if (it != _histograms.end()) {
                // wish I had C++20 with erase_if...
                for (auto it2 = it->second.begin(); it2 != it->second.end();) {
                    if (it2->first != Spec::none()) {
                        it2 = it->second.erase(it2);
                    } else {
                        it2++;
                    }
                }
            }
        }
    }

    //    IOHistogram Histos::getHistogram(const string& name, const string& scheme) const {
    //        auto entry = getEntry(name, scheme);
    //        if(entry != nullptr && entry->size() > 0) {
    //            IOHistogram result;
    //            for(auto it = begin(*entry); it != end(*entry); ++it) {
    //                if(it->second) {
    //                    if(result.size() == 0) {
    //                        result = it->second->evalToList(externalData);
    //                    }
    //                    else {
    //                        auto result_temp = it->second->evalToList(externalData);
    //                        transform(result.begin(), result.end(), result_temp.begin(),
    //                                    result.begin(), plus<>());
    //                    }
    //                }
    //            }
    //            return result;
    //        }
    //        MSG_ERROR("Wrong name/scheme in histogram request. Returning empty histogram list");
    //        return IOHistogram{};
    //    }

    bool Histos::isValidHistogram(const string& name, size_t dim) const {
        auto it = _histogramDefs.find(name);
        if (it == _histogramDefs.end()) {
            return false;
        }
        bool valid = true;
        const auto& def = it->second;
        valid &= (def.getIndexing().rank() == dim);
        return valid;
    }

#ifdef HAVE_ROOT

    static unique_ptr<TH1D> createTH1D(const string& name, const vector<vector<double>>& edges) {
        Int_t nbinsX = static_cast<Int_t>(edges[0].size()) - 1;
        auto result = std::make_unique<TH1D>(name.c_str(), name.c_str(), nbinsX, edges[0].data());
        return result;
    }

    static unique_ptr<TH2D> createTH2D(const string& name, const vector<vector<double>>& edges) {
        Int_t nbinsX = static_cast<Int_t>(edges[0].size()) - 1;
        Int_t nbinsY = static_cast<Int_t>(edges[1].size()) - 1;
        auto result =
            std::make_unique<TH2D>(name.c_str(), name.c_str(), nbinsX, edges[0].data(), nbinsY, edges[1].data());
        return result;
    }

    static unique_ptr<TH3D> createTH3D(const string& name, const vector<vector<double>>& edges) {
        Int_t nbinsX = static_cast<Int_t>(edges[0].size()) - 1;
        Int_t nbinsY = static_cast<Int_t>(edges[1].size()) - 1;
        Int_t nbinsZ = static_cast<Int_t>(edges[2].size()) - 1;
        auto result = std::make_unique<TH3D>(name.c_str(), name.c_str(), nbinsX, edges[0].data(), nbinsY,
                                             edges[1].data(), nbinsZ, edges[2].data());
        return result;
    }

    static void evalToTH1D(const IOHistogram& input, const IndexList& dimensions, TH1D& rootHistogram,
                           bool hasUnderOverFlow, bool add = false) {
        ASSERT(dimensions.size() == 1);
        size_t base = hasUnderOverFlow ? 0 : 1;
        ASSERT(static_cast<size_t>(rootHistogram.GetNbinsX()) + 2 * (1 - base) == dimensions[0]);
        double entries = 0.;
        array<double, 4> stats;
        if (add) {
            entries = rootHistogram.GetEntries();
            rootHistogram.GetStats(stats.data());
        }
        for (IndexType i1 = 0; i1 < dimensions[0]; ++i1) {
            double value = input[i1].sumWi;
            double err2 = input[i1].sumWi2;
            stats[0] += value;
            stats[1] += err2;
            if (add) {
                value += rootHistogram.GetBinContent(static_cast<Int_t>(i1 + base));
                err2 += pow(rootHistogram.GetBinError(static_cast<Int_t>(i1 + base)), 2);
            }
            rootHistogram.SetBinContent(static_cast<Int_t>(i1 + base), value);
            rootHistogram.SetBinError(static_cast<Int_t>(i1 + base), sqrt(err2));
            entries += static_cast<double>(input[i1].n);
        }
        rootHistogram.SetEntries(entries);
        rootHistogram.PutStats(stats.data());
    }

    static void evalToTH2D(const IOHistogram& input, const IndexList& dimensions, TH2D& rootHistogram,
                           function<PositionType(const IndexList&)>& calcCoords, bool hasUnderOverFlow,
                           bool add = false) {
        ASSERT(dimensions.size() == 2);
        size_t base = hasUnderOverFlow ? 0 : 1;
        ASSERT(static_cast<size_t>(rootHistogram.GetNbinsX()) + 2 * (1 - base) == dimensions[0]);
        ASSERT(static_cast<size_t>(rootHistogram.GetNbinsY()) + 2 * (1 - base) == dimensions[1]);
        double entries = 0;
        array<double, 7> stats;
        if (add) {
            entries = rootHistogram.GetEntries();
            rootHistogram.GetStats(stats.data());
        }
        for (IndexType i1 = 0; i1 < dimensions[0]; ++i1) {
            for (IndexType i2 = 0; i2 < dimensions[1]; ++i2) {
                double value = input[calcCoords({i1, i2})].sumWi;
                double err2 = input[calcCoords({i1, i2})].sumWi2;
                stats[0] += value;
                stats[1] += err2;
                if (add) {
                    value += rootHistogram.GetBinContent(static_cast<Int_t>(i1 + base), static_cast<Int_t>(i2 + base));
                    err2 +=
                        pow(rootHistogram.GetBinError(static_cast<Int_t>(i1 + base), static_cast<Int_t>(i2 + base)), 2);
                }
                rootHistogram.SetBinContent(static_cast<Int_t>(i1 + base), static_cast<Int_t>(i2 + base), value);
                rootHistogram.SetBinError(static_cast<Int_t>(i1 + base), static_cast<Int_t>(i2 + base), sqrt(err2));
                entries += static_cast<double>(input[calcCoords({i1, i2})].n);
            }
        }
        rootHistogram.SetEntries(entries);
        rootHistogram.PutStats(stats.data());
    }

    static void evalToTH3D(const IOHistogram& input, const IndexList& dimensions, TH3D& rootHistogram,
                           function<PositionType(const IndexList&)>& calcCoords, bool hasUnderOverFlow,
                           bool add = false) {
        ASSERT(dimensions.size() == 3);
        size_t base = hasUnderOverFlow ? 0 : 1;
        ASSERT(static_cast<size_t>(rootHistogram.GetNbinsX()) + 2 * (1 - base) == dimensions[0]);
        ASSERT(static_cast<size_t>(rootHistogram.GetNbinsY()) + 2 * (1 - base) == dimensions[1]);
        ASSERT(static_cast<size_t>(rootHistogram.GetNbinsZ()) + 2 * (1 - base) == dimensions[2]);
        double entries = 0;
        array<double, 11> stats;
        if (add) {
            entries = rootHistogram.GetEntries();
            rootHistogram.GetStats(stats.data());
        }
        for (IndexType i1 = 0; i1 < dimensions[0]; ++i1) {
            for (IndexType i2 = 0; i2 < dimensions[1]; ++i2) {
                for (IndexType i3 = 0; i3 < dimensions[2]; ++i3) {
                    double value = input[calcCoords({i1, i2, i3})].sumWi;
                    double err2 = input[calcCoords({i1, i2, i3})].sumWi2;
                    stats[0] += value;
                    stats[1] += err2;
                    if (add) {
                        value +=
                            rootHistogram.GetBinContent(static_cast<Int_t>(i1 + base), static_cast<Int_t>(i2 + base),
                                                        static_cast<Int_t>(i3 + base));
                        err2 +=
                            pow(rootHistogram.GetBinError(static_cast<Int_t>(i1 + base), static_cast<Int_t>(i2 + base),
                                                          static_cast<Int_t>(i3 + base)),
                                2);
                    }
                    rootHistogram.SetBinContent(static_cast<Int_t>(i1 + base), static_cast<Int_t>(i2 + base),
                                                static_cast<Int_t>(i3 + base), value);
                    rootHistogram.SetBinError(static_cast<Int_t>(i1 + base), static_cast<Int_t>(i2 + base),
                                              static_cast<Int_t>(i3 + base), sqrt(err2));
                    entries += static_cast<double>(input[calcCoords({i1, i2, i3})].n);
                }
            }
        }
        rootHistogram.SetEntries(entries);
        rootHistogram.PutStats(stats.data());
    }

    bool Histos::isValidHistogram(const string& name, const TH1D& h) const {
        auto it = _histogramDefs.find(name);
        if (it == _histogramDefs.end()) {
            return false;
        }
        bool valid = true;
        const auto& def = it->second;
        valid &= (def.getIndexing().rank() == 1);
        if (valid) {
            size_t base = def.getIndexing().hasUnderOverFlow() ? 1 : 0;
            valid &= (static_cast<size_t>(h.GetNbinsX()) + 2 * base == def.getIndexing().dims()[0]);
        }
        return valid;
    }

    bool Histos::isValidHistogram(const string& name, const TH2D& h) const {
        auto it = _histogramDefs.find(name);
        if (it == _histogramDefs.end()) {
            return false;
        }
        bool valid = true;
        const auto& def = it->second;
        valid &= (def.getIndexing().rank() == 2);
        if (valid) {
            size_t base = def.getIndexing().hasUnderOverFlow() ? 1 : 0;
            valid &= (static_cast<size_t>(h.GetNbinsX()) + 2 * base == def.getIndexing().dims()[0]);
            valid &= (static_cast<size_t>(h.GetNbinsY()) + 2 * base == def.getIndexing().dims()[1]);
        }
        return valid;
    }

    bool Histos::isValidHistogram(const string& name, const TH3D& h) const {
        auto it = _histogramDefs.find(name);
        if (it == _histogramDefs.end()) {
            return false;
        }
        bool valid = true;
        const auto& def = it->second;
        valid &= (def.getIndexing().rank() == 3);
        if (valid) {
            size_t base = def.getIndexing().hasUnderOverFlow() ? 1 : 0;
            valid &= (static_cast<size_t>(h.GetNbinsX()) + 2 * base == def.getIndexing().dims()[0]);
            valid &= (static_cast<size_t>(h.GetNbinsY()) + 2 * base == def.getIndexing().dims()[1]);
            valid &= (static_cast<size_t>(h.GetNbinsZ()) + 2 * base == def.getIndexing().dims()[2]);
        }
        return valid;
    }


    EventIdGroupDict<unique_ptr<TH1D>> Histos::getHistograms1D(const string& name, const string& scheme,
                                                               const string& specialization) const {
        if (isValidHistogram(name, 1)) {
            auto hammerHistos = getHistograms(name, scheme, specialization);
            if (!hammerHistos.empty()) {
                const auto& def = _histogramDefs.find(name)->second;
                EventIdGroupDict<unique_ptr<TH1D>> result;
                for (auto& hammerHisto : hammerHistos) {
                    auto th1d = createTH1D(name, def.getIndexing().edges());
                    evalToTH1D(hammerHisto.second, def.getIndexing().dims(), *th1d,
                               def.getIndexing().hasUnderOverFlow());
                    result.insert(make_pair(hammerHisto.first, std::move(th1d)));
                }
                return result;
            }
        }
        MSG_ERROR("Wrong name/scheme/dimensionality in histogram request. Returning empty histogram list");
        return EventIdGroupDict<unique_ptr<TH1D>>{};
    }

    EventIdGroupDict<unique_ptr<TH2D>> Histos::getHistograms2D(const string& name, const string& scheme,
                                                               const string& specialization) const {
        if (isValidHistogram(name, 2)) {
            auto hammerHistos = getHistograms(name, scheme, specialization);
            if (!hammerHistos.empty()) {
                EventIdGroupDict<unique_ptr<TH2D>> result;
                const auto& def = _histogramDefs.find(name)->second;
                function<PositionType(const IndexList&)> fun = [&](const IndexList& val) -> PositionType {
                    return def.getIndexing().indicesToPos(val);
                };
                for (auto& hammerHisto : hammerHistos) {
                    auto th2d = createTH2D(name, def.getIndexing().edges());
                    evalToTH2D(hammerHisto.second, def.getIndexing().dims(), *th2d, fun,
                               def.getIndexing().hasUnderOverFlow());
                    result.insert(make_pair(hammerHisto.first, std::move(th2d)));
                }
                return result;
            }
        }
        MSG_ERROR("Wrong name/scheme/dimensionality in histogram request. Returning empty histogram list");
        return EventIdGroupDict<unique_ptr<TH2D>>{};
    }

    EventIdGroupDict<unique_ptr<TH3D>> Histos::getHistograms3D(const string& name, const string& scheme,
                                                               const string& specialization) const {
        if (isValidHistogram(name, 3)) {
            auto hammerHistos = getHistograms(name, scheme, specialization);
            if (!hammerHistos.empty()) {
                EventIdGroupDict<unique_ptr<TH3D>> result;
                const auto& def = _histogramDefs.find(name)->second;
                function<PositionType(const IndexList&)> fun = [&](const IndexList& val) -> PositionType {
                    return def.getIndexing().indicesToPos(val);
                };
                for (auto& hammerHisto : hammerHistos) {
                    auto th3d = createTH3D(name, def.getIndexing().edges());
                    evalToTH3D(hammerHisto.second, def.getIndexing().dims(), *th3d, fun,
                               def.getIndexing().hasUnderOverFlow());
                    result.insert(make_pair(hammerHisto.first, std::move(th3d)));
                }
                return result;
            }
        }
        MSG_ERROR("Wrong name/scheme/dimensionality in histogram request. Returning empty histogram list");
        return EventIdGroupDict<unique_ptr<TH3D>>{};
    }

    unique_ptr<TH1D> Histos::getHistogram1D(const string& name, const string& scheme,
                                            const string& specialization) const {
        if (isValidHistogram(name, 1)) {
            auto hammerHisto = getHistogram(name, scheme, specialization);
            if (!hammerHisto.empty()) {
                const auto& def = _histogramDefs.find(name)->second;
                auto th1d = createTH1D(name, def.getIndexing().edges());
                evalToTH1D(hammerHisto, def.getIndexing().dims(), *th1d, def.getIndexing().hasUnderOverFlow());
                return th1d;
            }
        }
        MSG_ERROR("Wrong name/scheme/dimensionality in histogram request. Returning empty histogram list");
        return unique_ptr<TH1D>{};
    }

    unique_ptr<TH2D> Histos::getHistogram2D(const string& name, const string& scheme,
                                            const string& specialization) const {
        if (isValidHistogram(name, 2)) {
            auto hammerHisto = getHistogram(name, scheme, specialization);
            if (!hammerHisto.empty()) {
                const auto& def = _histogramDefs.find(name)->second;
                function<PositionType(const IndexList&)> fun = [&](const IndexList& val) -> PositionType {
                    return def.getIndexing().indicesToPos(val);
                };
                auto th2d = createTH2D(name, def.getIndexing().edges());
                evalToTH2D(hammerHisto, def.getIndexing().dims(), *th2d, fun, def.getIndexing().hasUnderOverFlow());
                return th2d;
            }
        }
        MSG_ERROR("Wrong name/scheme/dimensionality in histogram request. Returning empty histogram list");
        return unique_ptr<TH2D>{};
    }

    unique_ptr<TH3D> Histos::getHistogram3D(const string& name, const string& scheme,
                                            const string& specialization) const {
        if (isValidHistogram(name, 3)) {
            auto hammerHisto = getHistogram(name, scheme, specialization);
            if (!hammerHisto.empty()) {
                const auto& def = _histogramDefs.find(name)->second;
                function<PositionType(const IndexList&)> fun = [&](const IndexList& val) -> PositionType {
                    return def.getIndexing().indicesToPos(val);
                };
                auto th3d = createTH3D(name, def.getIndexing().edges());
                evalToTH3D(hammerHisto, def.getIndexing().dims(), *th3d, fun, def.getIndexing().hasUnderOverFlow());
                return th3d;
            }
        }
        MSG_ERROR("Wrong name/scheme/dimensionality in histogram request. Returning empty histogram list");
        return unique_ptr<TH3D>{};
    }

    void Histos::setHistogram1D(const string& name, const string& scheme, TH1D& rootHistogram,
                                const string& specialization) const {
        if (isValidHistogram(name, rootHistogram)) {
            auto hammerHisto = getHistogram(name, scheme, specialization);
            if (!hammerHisto.empty()) {
                const auto& def = _histogramDefs.find(name)->second;
                evalToTH1D(hammerHisto, def.getIndexing().dims(), rootHistogram, def.getIndexing().hasUnderOverFlow());
            } else {
                MSG_ERROR(string("Wrong scheme in setting ROOT histogram '") + name + string("'. Not set."));
            }
        } else {
            MSG_ERROR(string("The Hammer and ROOT histograms corresponding to '") + name +
                      string("' have incompatible layouts. Not set."));
        }
    }

    void Histos::setHistogram2D(const string& name, const string& scheme, TH2D& rootHistogram,
                                const string& specialization) const {
        if (isValidHistogram(name, rootHistogram)) {
            auto hammerHisto = getHistogram(name, scheme, specialization);
            if (!hammerHisto.empty()) {
                const auto& def = _histogramDefs.find(name)->second;
                function<PositionType(const IndexList&)> fun = [&](const IndexList& val) -> PositionType {
                    return def.getIndexing().indicesToPos(val);
                };
                evalToTH2D(hammerHisto, def.getIndexing().dims(), rootHistogram, fun,
                           def.getIndexing().hasUnderOverFlow());
            } else {
                MSG_ERROR(string("Wrong scheme in setting ROOT histogram '") + name + string("'. Not set."));
            }
        } else {
            MSG_ERROR(string("The Hammer and ROOT histograms corresponding to '") + name +
                      string("' have incompatible layouts. Not set."));
        }
    }

    void Histos::setHistogram3D(const string& name, const string& scheme, TH3D& rootHistogram,
                                const string& specialization) const {
        if (isValidHistogram(name, rootHistogram)) {
            auto hammerHisto = getHistogram(name, scheme, specialization);
            if (!hammerHisto.empty()) {
                const auto& def = _histogramDefs.find(name)->second;
                function<PositionType(const IndexList&)> fun = [&](const IndexList& val) -> PositionType {
                    return def.getIndexing().indicesToPos(val);
                };
                evalToTH3D(hammerHisto, def.getIndexing().dims(), rootHistogram, fun,
                           def.getIndexing().hasUnderOverFlow());
            } else {
                MSG_ERROR(string("Wrong scheme in setting ROOT histogram '") + name + string("'. Not set."));
            }
        } else {
            MSG_ERROR(string("The Hammer and ROOT histograms corresponding to '") + name +
                      string("' have incompatible layouts. Not set."));
        }
    }

    void Histos::setHistograms1D(const string& name, const string& scheme,
                                 EventIdGroupDict<unique_ptr<TH1D>>& rootHistograms,
                                 const string& specialization) const {
        if (isValidHistogram(name, *(rootHistograms.begin()->second))) {
            auto hammerHistos = getHistograms(name, scheme, specialization);
            if (!hammerHistos.empty()) {
                const auto& def = _histogramDefs.find(name)->second;
                for (auto& hammerHisto : hammerHistos) {
                    auto itRoot = rootHistograms.find(hammerHisto.first);
                    if (itRoot == rootHistograms.end()) {
                        auto th1d = createTH1D(name, def.getIndexing().edges());
                        auto res = rootHistograms.insert(make_pair(hammerHisto.first, std::move(th1d)));
                        itRoot = res.first;
                    }
                    evalToTH1D(hammerHisto.second, def.getIndexing().dims(), *(itRoot->second),
                               def.getIndexing().hasUnderOverFlow());
                }
            } else {
                MSG_ERROR(string("Wrong scheme in setting ROOT histogram '") + name + string("'. Not set."));
            }
        } else {
            MSG_ERROR(string("The Hammer and ROOT histograms corresponding to '") + name +
                      string("' have incompatible layouts. Not set."));
        }
    }

    void Histos::setHistograms2D(const string& name, const string& scheme,
                                 EventIdGroupDict<unique_ptr<TH2D>>& rootHistograms,
                                 const string& specialization) const {
        if (isValidHistogram(name, *(rootHistograms.begin()->second))) {
            auto hammerHistos = getHistograms(name, scheme, specialization);
            if (!hammerHistos.empty()) {
                const auto& def = _histogramDefs.find(name)->second;
                function<PositionType(const IndexList&)> fun = [&](const IndexList& val) -> PositionType {
                    return def.getIndexing().indicesToPos(val);
                };
                for (auto& hammerHisto : hammerHistos) {
                    auto itRoot = rootHistograms.find(hammerHisto.first);
                    if (itRoot == rootHistograms.end()) {
                        auto th2d = createTH2D(name, def.getIndexing().edges());
                        auto res = rootHistograms.insert(make_pair(hammerHisto.first, std::move(th2d)));
                        itRoot = res.first;
                    }
                    evalToTH2D(hammerHisto.second, def.getIndexing().dims(), *(itRoot->second), fun,
                               def.getIndexing().hasUnderOverFlow());
                }
            } else {
                MSG_ERROR(string("Wrong scheme in setting ROOT histogram '") + name + string("'. Not set."));
            }
        } else {
            MSG_ERROR(string("The Hammer and ROOT histograms corresponding to '") + name +
                      string("' have incompatible layouts. Not set."));
        }
    }

    void Histos::setHistograms3D(const string& name, const string& scheme,
                                 EventIdGroupDict<unique_ptr<TH3D>>& rootHistograms,
                                 const string& specialization) const {
        if (isValidHistogram(name, *(rootHistograms.begin()->second))) {
            auto hammerHistos = getHistograms(name, scheme, specialization);
            if (!hammerHistos.empty()) {
                const auto& def = _histogramDefs.find(name)->second;
                function<PositionType(const IndexList&)> fun = [&](const IndexList& val) -> PositionType {
                    return def.getIndexing().indicesToPos(val);
                };
                for (auto& hammerHisto : hammerHistos) {
                    auto itRoot = rootHistograms.find(hammerHisto.first);
                    if (itRoot == rootHistograms.end()) {
                        auto th3d = createTH3D(name, def.getIndexing().edges());
                        auto res = rootHistograms.insert(make_pair(hammerHisto.first, std::move(th3d)));
                        itRoot = res.first;
                    }
                    evalToTH3D(hammerHisto.second, def.getIndexing().dims(), *(itRoot->second), fun,
                               def.getIndexing().hasUnderOverFlow());
                }
            } else {
                MSG_ERROR(string("Wrong scheme in setting ROOT histogram '") + name + string("'. Not set."));
            }
        } else {
            MSG_ERROR(string("The Hammer and ROOT histograms corresponding to '") + name +
                      string("' have incompatible layouts. Not set."));
        }
    }

#endif

    void Histos::fillHisto(const string& name, const string& spec, const string& scheme, const IndexList& binPosition,
                           Tensor& value, double extraWeight) {
        auto itDef = _histogramDefs.find(name);
        if (itDef == _histogramDefs.end()) {
            MSG_ERROR(string("Unable to find histogram ") + name + string(". Histogram not filled."));
            return;
        }
        LabelsList newlabs = value.labels();
        IndexList newdims = value.dims();
        Tensor temp = itDef->second.getFixedData(_currentEventId, spec, scheme, newlabs);
        bool processed = false;
        if (temp.rank() > 0) {
            temp.dot(value);
            newlabs = temp.labels();
            newdims = temp.dims();
            if (!isZero(extraWeight - 1.0)) {
                temp *= extraWeight;
            }
            processed = true;
        } else {
            if (!isZero(extraWeight - 1.0)) {
                temp = value * extraWeight;
                processed = true;
            }
        }
        auto* entry = getEntry(name, spec, scheme);
        if (entry == nullptr) {
            MSG_ERROR(string("Unable to find histogram ") + name + string(". Histogram not filled."));
            return;
        }
        Histogram* histo = entry->getHistogram(_currentEventId);
        if (histo == nullptr) {
            auto id = entry->addEventId(_currentEventId, newlabs);
            histo = entry->getHistogram(id);
            if (histo == nullptr) {
                histo = entry->addHistogram(id, makeHistogram(name, itDef->second, newdims, newlabs));
            }
        }
        if (processed) {
            histo->fillBin(temp, binPosition);
        } else {
            histo->fillBin(value, binPosition);
        }
    }

    Log& Histos::getLog() {
        return Log::getLog("Hammer.Histos");
    }

    void Histos::defineSettings() {
        setPath("Histos");
        addSetting<bool>("KeepErrors", false);
    }

    bool Histos::writeDefinition(flatbuffers::FlatBufferBuilder* msgwriter, const string& name) const {
        auto ithisto = _histogramDefs.find(name);
        if (ithisto != _histogramDefs.end()) {
            ithisto->second.write(msgwriter, name);
            return true;
        }
        return false;
    }


    bool Histos::writeHistogram(flatbuffers::FlatBufferBuilder* msgwriter, const string& histogramName,
                                const string& specializationName, const string& schemeName,
                                const EventUID& eventIDs) const {
        auto ithisto = _histograms.find(histogramName);
        if (ithisto != _histograms.end()) {
            auto itspec = ithisto->second.find(specializationName);
            if (itspec != ithisto->second.end()) {
                auto itscheme = itspec->second.find(schemeName);
                if (itscheme != itspec->second.end()) {
                    auto serialScheme = msgwriter->CreateString(schemeName);
                    auto serialSpec = msgwriter->CreateString(specializationName);
                    auto res = itscheme->second.write(msgwriter, eventIDs);
                    res->add_scheme(serialScheme);
                    res->add_specialization(serialSpec);
                    auto out = res->Finish();
                    msgwriter->Finish(out);
                    return true;
                }
            }
        }
        return false;
    }


    string Histos::readDefinition(const Serial::FBHistoDefinition* msgreader, bool merge) {
        if (msgreader != nullptr) {
            if (msgreader->name() == nullptr) {
                MSG_ERROR("Missing name in reading Histograms. You need a bell and ring it in the morning.");
                return "";
            }
            auto itdef = _histogramDefs.find(msgreader->name()->c_str());
            if (itdef == _histogramDefs.end()) {
                itdef = _histogramDefs.emplace(msgreader->name()->c_str(), HistogramDefinition{msgreader}).first;
            } else {
                if (merge) {
                    if (!itdef->second.checkDefinition(msgreader)) {
                        MSG_ERROR("Invalid merge for histogram '" + itdef->first + "': definitions do not match.");
                        return "";
                    }
                } else {
                    itdef->second = HistogramDefinition{msgreader};
                }
            }
            addSetting<bool>("CalcGeneralHistogram_" + static_cast<string>(msgreader->name()->c_str()), true);
            return itdef->first;
        }
        MSG_ERROR("Invalid record for histogram definition");
        return "";
    }

    HistoInfo Histos::readHistogram(const Serial::FBHistogram* msgreader, bool merge) {
        if (msgreader != nullptr) {
            if (msgreader->name() == nullptr || msgreader->scheme() == nullptr) {
                MSG_ERROR("Missing name or scheme in reading Histogram. You need a song and sing it in the morning.");
                return HistoInfo{"", "", "", EventUIDGroup{}};
            }
            auto ithisto = _histograms.find(msgreader->name()->c_str());
            if (ithisto == _histograms.end()) {
                auto res =
                    _histograms.emplace(msgreader->name()->c_str(), SpecializationDict<SchemeDict<HistogramSet>>{});
                if (!res.second) {
                    MSG_ERROR(string("Unable to add read histogram '") + msgreader->name()->c_str() + string("'."));
                    return HistoInfo{"", "", "", EventUIDGroup{}};
                }
                ithisto = res.first;
            }
            string specName;
            if (flatbuffers::IsFieldPresent(msgreader, Serial::FBHistogram::VT_SPECIALIZATION)) {
                specName = msgreader->specialization()->c_str();
            } else {
                specName = Spec::none();
            }
            auto itspec = ithisto->second.find(specName);
            if (itspec == ithisto->second.end()) {
                auto res2 = ithisto->second.emplace(specName, SchemeDict<HistogramSet>{});
                if (!res2.second) {
                    MSG_ERROR(string("Unable to add read histogram '") + ithisto->first + string("'."));
                    return HistoInfo{"", "", "", EventUIDGroup{}};
                }
                itspec = res2.first;
            }
            auto itscheme = itspec->second.find(msgreader->scheme()->c_str());
            if (itscheme == itspec->second.end()) {
                auto res3 = itspec->second.emplace(msgreader->scheme()->c_str(), HistogramSet{});
                if (!res3.second) {
                    MSG_ERROR(string("Unable to add read histogram '") + ithisto->first + string("'."));
                    return HistoInfo{"", "", "", EventUIDGroup{}};
                }
                itscheme = res3.first;
            }
            auto itdef = _histogramDefs.find(ithisto->first);
            if (itdef == _histogramDefs.end()) {
                MSG_ERROR(string("Unable to add read histogram '") + ithisto->first +
                          string("': definition not found."));
                return HistoInfo{"", "", "", EventUIDGroup{}};
            }
            auto result = itscheme->second.read(msgreader, itdef->second, merge);
            if (!result.empty()) {
                if (specName != Spec::none()) {
                    _usedSpecializations[msgreader->name()->c_str()].insert(specName);
                }
                return HistoInfo{ithisto->first, itscheme->first, itspec->first, result};
            }
            MSG_ERROR("Invalid merge for histogram '" + ithisto->first + "'");
            return HistoInfo{"", "", "", EventUIDGroup{}};
        }
        MSG_ERROR("Invalid record for histogram");
        return HistoInfo{"", "", "", EventUIDGroup{}};
    }

} // namespace Hammer
