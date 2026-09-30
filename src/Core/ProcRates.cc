///
/// @file  ProcRates.cc
/// @brief Container class for process rate tensors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <boost/algorithm/string.hpp>

#include "Hammer/ProcRates.hh"
#include "Hammer/Tools/HammerSerial.hh"
#include "Hammer/ExternalData.hh"
#include "Hammer/Tools/Logging.hh"
#include "Hammer/Tools/Utils.hh"
#include "Hammer/Math/MultiDim/ScalarContainer.hh"

using namespace std;

namespace Hammer {

    ProcRates::ProcRates(const ExternalData* ext) : _external{ext} {
    }

    ProcRates::~ProcRates() noexcept {
        _processRates.clear();
        _external = nullptr;
    }

    void ProcRates::defineSettings() {
    }

    void ProcRates::init() {
        _processRates.clear();
    }

    SpecializationDict<SchemeDict<Tensor>>* ProcRates::getProcessRates(HashId id) {
        auto it = _processRates.find(id);
        if (it == _processRates.end()) {
            auto res = _processRates.insert({id, SpecializationDict<SchemeDict<Tensor>>{}});
            if (res.second) {
                return &(res.first->second);
            }
            return nullptr;
        }
        return &(it->second);
    }

    double ProcRates::getVertexRate(const HashId& id, const string& scheme, const string& spec) const {
        auto it = _processRates.find(id);
        if (it != _processRates.end()) {
            auto itsp = it->second.find(spec);
            if (itsp != it->second.end()) {
                auto itsc = itsp->second.find(scheme);
                if (itsc != itsp->second.end()) {
                    if (!itsc->second.labels().empty()) {
                        auto t = _external->getExternalVectors(scheme, itsc->second.labels());
                        t.dot(itsc->second);
                        if (t.rank() != 0) {
                            MSG_ERROR("Invalid rank at end of evaluation");
                        }
                        return t.element().real();
                    }
                    if (itsc->second.rank() != 0) {
                        MSG_ERROR("Invalid rank at end of evaluation");
                    }
                    return itsc->second.element().real();
                }
                MSG_ERROR("Rate not found for specified specialization '" + spec + "' and FF scheme '" + scheme +
                          "'. Oh, the humanity.");
                return 0.;
            }
            MSG_ERROR("Rate not found for specified specialization '" + spec + "'. Oh, the humanity.");
            return 0.;
        }
        MSG_ERROR("Rate not found for specified vertex ID. Oh, the humanity.");
        return 0.;
    }

    Log& ProcRates::getLog() {
        return Log::getLog("Hammer.ProcRates");
    }

    bool ProcRates::read(const Serial::FBRates* msgreader, bool merge) {
        const auto* rates = msgreader->rates();
        if (rates == nullptr) {
            return false;
        }
        if (!merge) {
            _processRates.clear();
        }
        bool result = true;
        for (unsigned int i = 0; i < rates->size() && result; ++i) {
            auto res = _processRates.insert({rates->Get(i)->id(), SpecializationDict<SchemeDict<Tensor>>{}});
            if (merge && !res.second) {
                res.first = _processRates.find(rates->Get(i)->id());
                if (res.first == _processRates.end()) {
                    return false;
                }
            }
            auto& curProcRateDict = res.first->second;
            const auto* vals = rates->Get(i)->ratevals();
            const auto* chs = rates->Get(i)->ratenames();
            const auto* sps = rates->Get(i)->ratespecs();
            if (vals == nullptr || chs == nullptr) {
                return false;
            }
            for (unsigned int j = 0; j < vals->size() && result; ++j) {
                const auto* spsStr = (sps != nullptr) ? sps->Get(j) : nullptr;
                auto specName = (spsStr != nullptr) ? spsStr->c_str() : Spec::none();
                auto res2 = curProcRateDict.insert({specName, SchemeDict<Tensor>{}});
                if (merge && !res2.second) {
                    res2.first = curProcRateDict.find(specName);
                    if (res2.first == curProcRateDict.end()) {
                        return false;
                    }
                }
                auto itspec = res2.first;
                const auto* chStr = chs->Get(j);
                if (chStr == nullptr) {
                    return false;
                }
                auto itscheme = itspec->second.find(chStr->c_str());
                if (merge && itscheme != itspec->second.end()) {
                    Tensor temp{};
                    temp.read(vals->Get(j));
                    if (!(temp.isEqualTo(itscheme->second))) {
                        MSG_ERROR(
                            "Trying to merge two rates for the same vertex within the same form factor scheme name, '" +
                            string(chStr->c_str()) + "' and same specialization name, '" + string(specName) +
                            "', but rates do not match! I would like to have seen Montana.");
                        return false;
                    }
                }
                auto res3 = itspec->second.insert({chStr->c_str(), Tensor{}});
                res3.first->second.read(vals->Get(j));
                //                if(merge && it != res.first->second.end() && strcmp(chs->Get(j)->c_str(),
                //                "Denominator") != 0 && strcmp(chs->Get(j)->c_str(), "NoFormFactor") != 0) {
                //                    MSG_ERROR("Trying to merge two rates with same form factor scheme name: '" +
                //                              string(chs->Get(j)->c_str()) + "'! I would like to have seen Montana.");
                //                    result = false;
                //
                //                }
                //                else {
                //                    auto res2 = res.first->second.insert({chs->Get(j)->c_str(), Tensor{}});
                //                    res2.first->second.read(vals->Get(j));
                //                }
            }
        }
        return result;
    }

    void ProcRates::write(flatbuffers::FlatBufferBuilder* msgwriter) const {
        vector<flatbuffers::Offset<Serial::FBRate>> rates;
        rates.reserve(_processRates.size());
        for (const auto& elem : _processRates) {
            vector<flatbuffers::Offset<flatbuffers::String>> chs;
            vector<flatbuffers::Offset<flatbuffers::String>> sps;
            vector<flatbuffers::Offset<Serial::FBTensor>> rts;
            for (const auto& elem2 : elem.second) {
                for (const auto& elem3 : elem2.second) {
                    auto serSp = msgwriter->CreateString(elem2.first);
                    auto serCh = msgwriter->CreateString(elem3.first);
                    flatbuffers::Offset<Serial::FBTensor> val;
                    elem3.second.write(&(*msgwriter), &val);
                    sps.push_back(serSp);
                    chs.push_back(serCh);
                    rts.push_back(val);
                }
            }
            auto serialChs = msgwriter->CreateVector(chs);
            auto serialrts = msgwriter->CreateVector(rts);
            auto serialSps = msgwriter->CreateVector(sps);
            Serial::FBRateBuilder serialFF{*msgwriter};
            serialFF.add_id(elem.first);
            serialFF.add_ratevals(serialrts);
            serialFF.add_ratenames(serialChs);
            serialFF.add_ratespecs(serialSps);
            auto resFF = serialFF.Finish();
            rates.push_back(resFF);
        }
        auto serialrates = msgwriter->CreateVector(rates);
        Serial::FBRatesBuilder serialout{*msgwriter};
        serialout.add_rates(serialrates);
        auto headoffset = serialout.Finish();
        msgwriter->Finish(headoffset);
    }


} // namespace Hammer
