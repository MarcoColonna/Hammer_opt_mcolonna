///
/// @file  HistogramSet.cc
/// @brief Container class for histograms belonging to different event types
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <functional>

#include "Hammer/Math/HistogramSet.hh"
#include "Hammer/Math/Histogram.hh"
#include "Hammer/Math/HistogramDefinition.hh"
#include "Hammer/Tools/HammerSerial.hh"

using namespace std;

namespace Hammer {

    HistogramSet::HistogramSet(bool compressed) : _isCompressed{compressed} {
    }

    HistogramSet::HistogramSet(const HistogramSet& other, const string& newName, const HistogramDefinition& newDef,
                               const std::set<uint16_t>& collapsedIndexPositions)
        : _compressionDict{other._compressionDict}, _compressionReverseDict{other._compressionReverseDict},
          _labelsDict{other._labelsDict}, _labelsReverseDict{other._labelsReverseDict},
          _isCompressed{other._isCompressed} {
        transform(other._data.begin(), other._data.end(), back_inserter(_data),
                  [&](const unique_ptr<Histogram>& original) -> unique_ptr<Histogram> {
                      return original->collapseIntoNew(newName, newDef, collapsedIndexPositions);
                  });
    }

    Histogram* HistogramSet::getHistogram(const EventUID& id) {
        auto it = _compressionDict.find(id);
        if (it != _compressionDict.end()) {
            return _data[it->second].get();
        }
        return nullptr;
    }

    Histogram* HistogramSet::getHistogram(size_t id) {
        if (id < _data.size()) {
            return _data[id].get();
        }
        return nullptr;
    }

    size_t HistogramSet::addEventId(const EventUID& id, const LabelsList& labels) {
        LabelsList tmpLabels = labels;
        sort(tmpLabels.begin(), tmpLabels.end());
        auto it = _compressionDict.find(id);
        if (it != _compressionDict.end()) {
            return it->second;
        }
        size_t pos = _data.size();
        if (_isCompressed) {
            auto itl = _labelsReverseDict.find(tmpLabels);
            if (itl != _labelsReverseDict.end()) {
                pos = itl->second;
                _compressionReverseDict.find(pos)->second.insert(id);
                _compressionDict.insert({id, pos});
                return pos;
            }
            _labelsReverseDict.insert({tmpLabels, pos});
        }
        _labelsDict.insert({pos, tmpLabels});
        _compressionDict.insert({id, pos});
        _compressionReverseDict.insert({pos, {id}});
        return pos;
    }

    Histogram* HistogramSet::addHistogram(size_t id, std::unique_ptr<Histogram> hist) {
        ASSERT(id <= _data.size());
        if (id == _data.size()) {
            _data.push_back(std::move(hist));
        }
        return _data[id].get();
    }

    void HistogramSet::clear() {
        _data.clear();
        _compressionDict.clear();
        _compressionReverseDict.clear();
        _labelsDict.clear();
        _labelsReverseDict.clear();
    }

    EventUIDGroup HistogramSet::read(const Serial::FBHistogram* msgreader, const HistogramDefinition& def, bool merge) {
        const auto* idset = msgreader->idset();
        if (idset == nullptr) {
            return {};
        }
        const auto* idsets = idset->evtids();
        if (idsets == nullptr) {
            return {};
        }
        vector<EventUID> keys;
        keys.reserve(idsets->size());
        for (unsigned int i = 0; i < idsets->size(); ++i) {
            const auto* idsetElem = idsets->Get(i);
            if (idsetElem == nullptr || idsetElem->ids() == nullptr) {
                return {};
            }
            const auto* procids = idsetElem->ids();
            EventUID key;
            for (unsigned int j = 0; j < procids->size(); ++j) {
                key.insert(procids->Get(j));
            }
            keys.push_back(key);
        }
        size_t newK = _data.size();
        set<size_t> candKs;
        auto itNewKeys = partition(keys.begin(), keys.end(), [&](const EventUID& elem) -> bool {
            auto it = _compressionDict.find(elem);
            if (it != _compressionDict.end()) {
                candKs.insert(it->second);
                return true;
            }
            return false;
        });
        if (candKs.size() > 1) {
            throw Error(
                "Event ID grouping in Histogram merging is not consistent.\n"
                "Mahna Mahna\nDo doo be-do-do\nMahna Mahna\nDo do-do do\n"
                "Mahna Mahna\nDo doo be-do-do be-do-do be-do-do be-do-do-doodle do do do-doo do!");
        }
        auto itCompRev = _compressionReverseDict.begin();
        if (!candKs.empty()) {
            newK = *candKs.begin();
            itCompRev = _compressionReverseDict.find(newK);
        } else {
            itCompRev = _compressionReverseDict.insert({newK, set<EventUID>{}}).first;
        }
        for (auto it = itNewKeys; it != keys.end(); ++it) {
            _compressionDict.insert({*it, newK});
            itCompRev->second.insert(*it);
        }
        if (newK == _data.size()) {
            const auto* idlabels = idset->labels();
            if (idlabels == nullptr) {
                return {};
            }
            LabelsList labels;
            for (unsigned int l = 0; l < idlabels->size(); ++l) {
                labels.push_back(static_cast<IndexLabel>(idlabels->Get(l)));
            }
            const auto* idpadpos = idset->labelpadsindices();
            const auto* idpadvals = idset->labelpadsvalues();
            if (idpadpos != nullptr && idpadvals != nullptr) {
                for (unsigned int l = 0; l < idpadpos->size(); ++l) {
                    auto idx = idpadpos->Get(l);
                    uint32_t pad = idpadvals->Get(l);
                    labels[idx] = specializeLabel(labels[idx], pad);
                }
            }
            _labelsDict.insert({newK, labels});
            if (_isCompressed) {
                _labelsReverseDict.insert({labels, newK});
            }
            _data.emplace_back(std::make_unique<Histogram>());
            _data.back()->read(msgreader, def);
        } else {
            if (merge) {
                Histogram h{};
                h.read(msgreader, def);
                *_data[newK] += h;
            } else {
                _data[newK] = std::make_unique<Histogram>();
                _data[newK]->read(msgreader, def);
            }
        }
        EventUIDGroup result{keys.begin(), keys.end()};
        return result;
    }

    unique_ptr<Serial::FBHistogramBuilder> HistogramSet::write(flatbuffers::FlatBufferBuilder* msgwriter,
                                                               const EventUID& id) const {
        auto it = _compressionDict.find(id);
        if (it != _compressionDict.end()) {
            const auto& evtids = _compressionReverseDict.find(it->second)->second;
            vector<flatbuffers::Offset<Serial::FBIdSet>> serialEvtIds;
            for (const auto& elem : evtids) {
                vector<uint64_t> procids;
                procids.reserve(elem.size());
                copy(elem.begin(), elem.end(), back_inserter(procids));
                auto serialProcids = msgwriter->CreateVector(procids);
                Serial::FBIdSetBuilder serialEventId{*msgwriter};
                serialEventId.add_ids(serialProcids);
                auto key = serialEventId.Finish();
                serialEvtIds.push_back(key);
            }
            auto serialIds = msgwriter->CreateVector(serialEvtIds);
            const auto& evtLabels = _labelsDict.find(it->second)->second;
            vector<int16_t> labels;
            vector<uint32_t> padvals;
            vector<uint8_t> padpos;
            for (uint8_t i = 0; i < static_cast<uint8_t>(evtLabels.size()); ++i) {
                auto elem = evtLabels[i];
                if (isSpecializedIndex(elem)) {
                    auto tmplabel = getSpecializedLabelComponents(elem);
                    labels.push_back(static_cast<int16_t>(tmplabel.second));
                    padpos.push_back(i);
                    padvals.push_back(tmplabel.first);
                } else {
                    labels.push_back(static_cast<int16_t>(elem));
                }
            }
            auto serialLabels = msgwriter->CreateVector(labels);
            auto serialpadpos = msgwriter->CreateVector(padpos);
            auto serialpadvals = msgwriter->CreateVector(padvals);
            Serial::FBEvtIdSetBuilder serialEventIdSet{*msgwriter};
            serialEventIdSet.add_evtids(serialIds);
            serialEventIdSet.add_labels(serialLabels);
            serialEventIdSet.add_labelpadsindices(serialpadpos);
            serialEventIdSet.add_labelpadsvalues(serialpadvals);
            auto idset = serialEventIdSet.Finish();
            auto hbuilder = _data[it->second]->write(msgwriter);
            hbuilder->add_idset(idset);
            return hbuilder;
        }
        return nullptr;
    }

    static BinContents operator+(const BinContents& a, const BinContents& b) {
        BinContents result = a;
        result.sumWi += b.sumWi;
        result.sumWi2 += b.sumWi2;
        result.n += b.n;
        return result;
    }

    EventIdGroupDict<IOHistogram> HistogramSet::specializeEventHistograms(const vector<Tensor>& externalData) const {
        EventIdGroupDict<IOHistogram> result;
        vector<IOHistogram> iodata;
        for (size_t i = 0; i < _data.size(); ++i) {
            auto it = _compressionReverseDict.find(i);
            result.emplace(it->second, _data[i]->evalToList(externalData[i]));
        }
        return result;
    }

    IOHistogram HistogramSet::specializeSumHistogram(const vector<Tensor>& externalData) const {
        if (_data.empty()) {
            return {};
        }
        IOHistogram result = _data[0]->evalToList(externalData[0]);
        for (size_t i = 1; i < _data.size(); ++i) {
            auto result_temp = _data[i]->evalToList(externalData[i]);
            transform(result.begin(), result.end(), result_temp.begin(), result.begin(), plus<>());
        }
        return result;
    }

    EventUIDGroup HistogramSet::getEventIdsInHistogram() const {
        EventUIDGroup result;
        for (const auto& elem : _compressionDict) {
            result.insert(elem.first);
        }
        return result;
    }

    vector<LabelsList> HistogramSet::getHistogramLabels() const {
        ASSERT(_labelsDict.size() == _data.size());
        vector<LabelsList> result;
        result.reserve(_labelsDict.size());
        transform(_labelsDict.begin(), _labelsDict.end(), back_inserter(result),
                  [](const pair<size_t, LabelsList>& val) -> LabelsList { return val.second; });
        return result;
    }

    vector<EventUID> HistogramSet::getEventUIDRepresentatives() const {
        vector<EventUID> result;
        result.reserve(_compressionReverseDict.size());
        transform(_compressionReverseDict.begin(), _compressionReverseDict.end(), back_inserter(result),
                  [](const pair<size_t, set<EventUID>>& val) -> EventUID { return *val.second.begin(); });
        return result;
    }

} // namespace Hammer
