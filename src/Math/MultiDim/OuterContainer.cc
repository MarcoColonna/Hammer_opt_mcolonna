///
/// @file  OuterContainer.cc
/// @brief (Sum of) Outer product tensor data container
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cstddef>
#include <numeric>
#include <iterator>

#include "Hammer/Math/MultiDim/OuterContainer.hh"
#include "Hammer/Math/MultiDim/ScalarContainer.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"
#include "Hammer/Math/MultiDim/VectorContainer.hh"
#include "Hammer/Math/MultiDim/Ops/Optimize.hh"
#include "Hammer/Exceptions.hh"
#include "Hammer/Math/MultiDimensional.fhh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Tools/HammerSerial.hh"
#include "Hammer/Tools/Logging.hh"

using namespace std;

namespace Hammer::MultiDimensional {

    OuterContainer::OuterContainer(const OuterContainer& other) : IContainer{other} {
        if (!other._sharedData) {
            auto tmp = other.clone();
            auto* newouter = static_cast<OuterContainer*>(tmp.get());
            _data.swap(newouter->_data);
            _indexing = newouter->_indexing;
            _accessors.swap(newouter->_accessors);
            _sharedData = false;
        } else {
            _data = other._data;
            _indexing = other._indexing;
            _accessors = other._accessors;
            _sharedData = true;
        }
    }

    OuterContainer& OuterContainer::operator=(const OuterContainer& other) {
        if (!other._sharedData) {
            auto tmp = other.clone();
            auto* newouter = static_cast<OuterContainer*>(tmp.get());
            _data.swap(newouter->_data);
            _accessors.swap(newouter->_accessors);
            _indexing = newouter->_indexing;
            _sharedData = false;
        } else {
            _data = other._data;
            _indexing = other._indexing;
            _accessors = other._accessors;
            _sharedData = true;
        }
        return *this;
    }

    OuterContainer::OuterContainer(TensorData left, TensorData right)
        : _indexing{{left->dims(), right->dims()}, {left->labels(), right->labels()}}, _sharedData{false} {
        size_t pos = left->rank();
        vector<pair<SharedTensorData, bool>> vec{{SharedTensorData{left.release()}, false},
                                                 {SharedTensorData{right.release()}, false}};
        ASSERT(vec[0].first.get() != nullptr);
        ASSERT(vec[1].first.get() != nullptr);
        _data.push_back(vec);
        _accessors = {[pos](const IndexList& list, IContainer* elem) -> ElementType {
                          return elem->element(list.begin(), list.begin() + static_cast<ptrdiff_t>(pos));
                      },
                      [pos](const IndexList& list, IContainer* elem) -> ElementType {
                          return elem->element(list.begin() + static_cast<ptrdiff_t>(pos), list.end());
                      }};
    }

    OuterContainer::OuterContainer(TensorData toBeSquared, bool conjugate)
        : _indexing{
              {toBeSquared->dims(), toBeSquared->dims()},
              {toBeSquared->labels(), conjugate ? flipListOfLabels(toBeSquared->labels()) : toBeSquared->labels()}},
          _sharedData{false} {
        size_t pos = toBeSquared->rank();
        auto elem = SharedTensorData{toBeSquared.release()};
        ASSERT(elem.get() != nullptr);
        vector<pair<SharedTensorData, bool>> vec{{elem, false}, {elem, conjugate}};
        _data.push_back(vec);
        _accessors = {[pos](const IndexList& list, IContainer* item) -> ElementType {
                          return item->element(list.begin(), list.begin() + static_cast<ptrdiff_t>(pos));
                      },
                      [pos](const IndexList& list, IContainer* item) -> ElementType {
                          return item->element(list.begin() + static_cast<ptrdiff_t>(pos), list.end());
                      }};
    }

    OuterContainer::OuterContainer(vector<TensorData>&& group) : _sharedData{false} {
        vector<pair<SharedTensorData, bool>> vec;
        std::vector<IndexList> dimensions;
        std::vector<LabelsList> tensLabels;
        dimensions.reserve(group.size());
        tensLabels.reserve(group.size());
        vec.reserve(group.size());
        size_t start = 0ul;
        size_t finish = 0ul;
        for (const auto& elem : group) {
            finish += elem->rank();
            _accessors.emplace_back([start, finish](const IndexList& list, const IContainer* item) -> ElementType {
                return item->element(list.begin() + static_cast<ptrdiff_t>(start),
                                     list.begin() + static_cast<ptrdiff_t>(finish));
            });
            start = finish;
        }
        for (auto& elem : group) {
            dimensions.push_back(elem->dims());
            tensLabels.push_back(elem->labels());
            vec.emplace_back(SharedTensorData{elem.release()}, false);
        }
        _data.emplace_back(vec);
        _indexing = BlockIndexing{dimensions, tensLabels};
    }

    OuterContainer::OuterContainer(SharedTensorData toBeSquared, bool conjugate)
        : _indexing{
              {toBeSquared->dims(), toBeSquared->dims()},
              {toBeSquared->labels(), conjugate ? flipListOfLabels(toBeSquared->labels()) : toBeSquared->labels()}},
          _sharedData{true} {
        vector<pair<SharedTensorData, bool>> vec{{toBeSquared, false}, {toBeSquared, conjugate}};
        _data.push_back(vec);
        size_t pos = toBeSquared->rank();
        _accessors = {[pos](const IndexList& list, IContainer* item) -> ElementType {
                          return item->element(list.begin(), list.begin() + static_cast<ptrdiff_t>(pos));
                      },
                      [pos](const IndexList& list, IContainer* item) -> ElementType {
                          return item->element(list.begin() + static_cast<ptrdiff_t>(pos), list.end());
                      }};
    }

    OuterContainer::OuterContainer(OuterContainer::EntryType&& data) : _sharedData{true} {
        std::vector<IndexList> dimensions;
        std::vector<LabelsList> tensLabels;
        dimensions.reserve(data.size());
        tensLabels.reserve(data.size());
        for (auto& i : data) {
            dimensions.push_back(i.first->dims());
            if (i.second) {
                tensLabels.push_back(flipListOfLabels(i.first->labels()));
            } else {
                tensLabels.push_back(i.first->labels());
            }
        }
        size_t start = 0ul;
        size_t finish = 0ul;
        for (auto& elem : data) {
            finish += elem.first->rank();
            _accessors.emplace_back([start, finish](const IndexList& list, const IContainer* item) -> ElementType {
                return item->element(list.begin() + static_cast<ptrdiff_t>(start),
                                     list.begin() + static_cast<ptrdiff_t>(finish));
            });
            start = finish;
        }
        _data.push_back(std::move(data));
        _indexing = BlockIndexing{dimensions, tensLabels};
    }

    OuterContainer::OuterContainer(const Serial::FBTensorList* input) : _sharedData{false} {
        if (input != nullptr) {
            _sharedData = false;
            const auto* serialvals = input->vals();
            const auto* serialentries = input->entries();
            if (serialvals != nullptr && serialentries != nullptr) {
                vector<SharedTensorData> allEntries;
                allEntries.reserve(serialvals->size());
                for (unsigned int i = 0; i < serialvals->size(); ++i) {
                    const auto* sv = serialvals->Get(i);
                    if (sv == nullptr || sv->vals() == nullptr || sv->dims() == nullptr) {
                        continue;
                    }
                    size_t fills = sv->vals()->size();
                    size_t totals = 1ul;
                    for (unsigned int j = 0; j < sv->dims()->size(); ++j) {
                        totals *= sv->dims()->Get(j);
                    }
                    // IMPORTANT: this assumes tensor was optimized before writing!!
                    if (Ops::shouldBeSparse(fills, totals)) {
                        allEntries.emplace_back(SharedTensorData{new SparseContainer{sv}});
                    } else {
                        allEntries.emplace_back(SharedTensorData{new VectorContainer{sv}});
                    }
                }
                _data.reserve(serialentries->size());
                for (unsigned int i = 0; i < serialentries->size(); ++i) {
                    const auto* serialitem = serialentries->Get(i);
                    if (serialitem->idx() == nullptr || serialitem->conj() == nullptr) {
                        continue;
                    }
                    vector<pair<SharedTensorData, bool>> vec;
                    vec.reserve(serialitem->idx()->size());
                    for (unsigned int j = 0; j < serialitem->idx()->size(); ++j) {
                        vec.emplace_back(allEntries[serialitem->idx()->Get(j)], serialitem->conj()->Get(j));
                    }
                    _data.push_back(vec);
                }
            }
            _accessors.clear();
            if (!_data.empty()) {
                auto entry = _data.begin();
                vector<LabelsList> tmplabels;
                vector<IndexList> tmpdims;
                tmplabels.reserve(entry->size());
                tmpdims.reserve(entry->size());
                size_t start = 0ul;
                size_t finish = 0ul;
                for (auto& it : *entry) {
                    tmplabels.push_back(((it.second) ? flipListOfLabels((it.first)->labels()) : (it.first)->labels()));
                    tmpdims.push_back((it.first)->dims());
                    finish += it.first->rank();
                    _accessors.emplace_back(
                        [start, finish](const IndexList& list, const IContainer* item) -> ElementType {
                            return item->element(list.begin() + static_cast<ptrdiff_t>(start),
                                                 list.begin() + static_cast<ptrdiff_t>(finish));
                        });
                    start = finish;
                }
                _indexing = BlockIndexing{tmpdims, tmplabels};
            }
        } else {
            throw InitializationError("Outer Container null buffer on read: have fun storming the castle!");
        }
    }


    OuterContainer::iterator OuterContainer::begin() {
        return _data.begin();
    }

    OuterContainer::const_iterator OuterContainer::begin() const {
        return _data.begin();
    }

    OuterContainer::iterator OuterContainer::end() {
        return _data.end();
    }

    OuterContainer::const_iterator OuterContainer::end() const {
        return _data.end();
    }


    OuterContainer::ElementType OuterContainer::value(const IndexList& indices) const {
        ASSERT(_indexing.checkValidIndices(indices));
        // auto split = _indexing.splitIndices(indices);
        // for(auto& elem: _data) {
        //     ElementType item = 1.;
        //     auto its = split.begin();
        //     auto it = elem.begin();
        //     for (; it != elem.end(); ++it, ++its) {
        //         auto tmp = it->first->element(*its);
        //         item *= it->second ? conj(tmp) : tmp;
        //     }
        //     result += item;
        // }
        ElementType result = 0;
        for (auto& elem : _data) {
            ElementType item = 1.;
            ASSERT(elem.size() == _accessors.size());
            for (size_t i = 0; i < elem.size(); ++i) {
                auto tmp = _accessors[i](indices, elem[i].first.get());
                item *= elem[i].second ? conj(tmp) : tmp;
            }
            result += item;
        }
        return result;
    }

    OuterContainer::ElementType OuterContainer::value(IndexList::const_iterator first,
                                                      IndexList::const_iterator last) const {
        ASSERT(_indexing.checkValidIndices(first, last));
        // ElementType result = 0;
        // auto split = _indexing.splitIndices(first, last);
        // for (auto& elem : _data) {
        //     ElementType item = 1.;
        //     auto its = split.begin();
        //     for (auto it = elem.begin(); it != elem.end(); ++it) {
        //         auto b = *its;
        //         ++its;
        //         auto e = *its;
        //         IndexList idx{b, e};
        //         auto tmp = it->first->element(idx);
        //         item *= it->second ? conj(tmp) : tmp;
        //     }
        //     result += item;
        // }
        // return result;
        IndexList indices{first, last};
        return value(indices);
    }

    OuterContainer::ElementType OuterContainer::value(const vector<IndexList>& indices) const {
        ASSERT(_indexing.checkValidIndices(indices));
        ElementType result = 0;
        for (auto& elem : _data) {
            auto it = elem.begin();
            auto iti = indices.begin();
            ElementType item = 1.;
            for (; it != elem.end(); ++it, ++iti) {
                auto tmp = it->first->element(*iti);
                item *= it->second ? conj(tmp) : tmp;
            }
            result += item;
        }
        return result;
    }

    size_t OuterContainer::numAddends() const {
        return _data.size();
    }

    size_t OuterContainer::rank() const {
        return _indexing.rank();
    }

    IndexList OuterContainer::dims() const {
        return _indexing.dims();
    }

    LabelsList OuterContainer::labels() const {
        return _indexing.labels();
    }

    size_t OuterContainer::numValues() const {
        return _indexing.numValues();
    }

    IndexType OuterContainer::labelToIndex(IndexLabel label) const {
        return _indexing.labelIndex(label);
    }

    IndexPairList OuterContainer::getSameLabelPairs(const IContainer& other, const UniqueLabelsList& indices) const {
        return _indexing.getSameLabelPairs(other.labels(), indices);
    }

    IndexPairList OuterContainer::getSpinLabelPairs() const {
        return _indexing.getOppositeLabelPairs(_indexing.spinIndices());
    }

    bool OuterContainer::isSameShape(const IContainer& other) const {
        return _indexing.isSameLabelShape(other.labels(), other.dims());
    }

    bool OuterContainer::canAddAt(const IContainer& /*subContainer*/, IndexLabel /*coord*/,
                                  IndexType /*position*/) const {
        return false;
    }

    OuterContainer::reference OuterContainer::element(const IndexList& /*coords*/) {
        throw Error("Cannot access element of OuterContainer by reference");
    }

    OuterContainer::ElementType OuterContainer::element(const IndexList& coords) const {
        return value(coords);
    }

    OuterContainer::reference OuterContainer::element(IndexList::const_iterator /*start*/,
                                                      IndexList::const_iterator /*end*/) {
        throw Error("Cannot access element of OuterContainer by reference");
    }

    OuterContainer::ElementType OuterContainer::element(IndexList::const_iterator start,
                                                        IndexList::const_iterator end) const {
        return value(start, end);
    }

    bool OuterContainer::isDataShared() const {
        return _sharedData;
    }

    size_t OuterContainer::dataSize() const {
        map<IContainer*, size_t> occupancies;
        for (auto& elem : _data) {
            for (auto& elem2 : elem) {
                auto it = occupancies.find(elem2.first.get());
                if (it == occupancies.end()) {
                    occupancies.insert({elem2.first.get(), elem2.first->dataSize()});
                }
            }
        }
        return accumulate(
            occupancies.begin(), occupancies.end(), 0ul,
            [](size_t var, const pair<IContainer*, size_t>& elem) -> size_t { return var + elem.second; });
    }

    size_t OuterContainer::entrySize() const {
        return sizeof(complex<double>);
    }

    LabelSignature OuterContainer::getLabelSignature() const {
        return _indexing.getLabelSignature();
    }


    bool OuterContainer::shouldBeEvaluated() const {
        if (_sharedData) {
            return false;
        }
        map<IContainer*, size_t> occupancies;
        double maxFillingFactor = 0.;
        for (auto& elem : _data) {
            double tmpFillingFactor = 1.;
            for (auto& elem2 : elem) {
                auto it = occupancies.find(elem2.first.get());
                if (it == occupancies.end()) {
                    it = occupancies.insert({elem2.first.get(), elem2.first->dataSize()}).first;
                }
                tmpFillingFactor *= (static_cast<double>(it->second)) /
                                    (static_cast<double>(elem2.first->numValues() * elem2.first->entrySize()));
            }
            maxFillingFactor = max(maxFillingFactor, tmpFillingFactor);
        }
        size_t result =
            accumulate(occupancies.begin(), occupancies.end(), 0ul,
                       [](size_t var, const pair<IContainer*, size_t>& elem) -> size_t { return var + elem.second; });
        return (maxFillingFactor * static_cast<double>(numValues() * (sizeof(complex<double>) + sizeof(PositionType))) <
                static_cast<double>(result));
    }

    bool OuterContainer::compare(const IContainer& other) const {
        const auto* tmp = dynamic_cast<const OuterContainer*>(&other);
        if (tmp == nullptr) {
            return false;
        }
        if (_sharedData != tmp->_sharedData) {
            return false;
        }
        bool dimcheck = _indexing.isSameLabelShape(tmp->_indexing);
        if (!dimcheck) {
            return false;
        }
        if (_data.size() != tmp->_data.size()) {
            return false;
        }
        if (_data.empty()) {
            return true;
        }
        if (_data.begin()->size() != tmp->_data.begin()->size()) {
            return false;
        }
        return equal(
            _data.begin(), _data.end(), tmp->_data.begin(),
            [](const vector<pair<SharedTensorData, bool>>& a, const vector<pair<SharedTensorData, bool>>& b) -> bool {
                bool result = true;
                auto ita = a.begin();
                auto itb = b.begin();
                for (; ita != a.end(); ++ita, ++itb) {
                    result &= (ita->second == itb->second) && ita->first->compare(*(itb->first));
                }
                return result;
            });
    }

    TensorData OuterContainer::clone() const {
        // always do a deep copy even if sharedData is true to avoid changes in contents downstream of dot if dot is
        // semi-trivial. clone should not be triggered when sharedData is true for moving/copying around the
        // outercontainer (i.e. make sure to use standard move/copy semantics)
        auto* t = new OuterContainer{};
        t->_sharedData = false;
        t->_indexing = _indexing;
        map<IContainer*, SharedTensorData> others;
        for (auto& elem : _data) {
            for (auto& elem2 : elem) {
                if (others.find(elem2.first.get()) == others.end()) {
                    others.insert({elem2.first.get(), SharedTensorData{elem2.first->clone().release()}});
                }
            }
        }
        for (auto& elem : _data) {
            vector<pair<SharedTensorData, bool>> tmp;
            transform(elem.begin(), elem.end(), back_inserter(tmp),
                      [others](const pair<SharedTensorData, bool>& in) -> pair<SharedTensorData, bool> {
                          return {others.find(in.first.get())->second, in.second};
                      });
            t->_data.push_back(tmp);
        }
        t->_accessors = _accessors;
        return TensorData{static_cast<IContainer*>(t)};
    }


    void OuterContainer::clear() {
        _data.clear();
    }

    void OuterContainer::swapElement(const IContainer* oldContainer, TensorData newContainer) {
        SharedTensorData newshared{std::move(newContainer)};
        for (auto& elem : _data) {
            for (auto& elem2 : elem) {
                if (elem2.first.get() == oldContainer) {
                    elem2.first = newshared;
                }
            }
        }
    }

    set<IContainer*> OuterContainer::getUniquePtrs(size_t position, bool decoupleConjugates) {
        map<IContainer*, bool> status;
        map<IContainer*, SharedTensorData> others;
        set<IContainer*> result{};
        ASSERT(position < _data.size());
        auto& elem = _data[position];
        if (!decoupleConjugates) {
            for (auto& elem2 : elem) {
                result.insert(elem2.first.get());
            }
        } else {
            set<SharedTensorData> toBeConjugate;
            for (auto& it : elem) {
                auto its = status.find(it.first.get());
                if (its != status.end()) {
                    if (its->second != it.second) {
                        auto ito = others.find(it.first.get());
                        if (ito == others.end()) {
                            SharedTensorData newentry{it.first->clone()};
                            if (it.second) {
                                toBeConjugate.insert(newentry);
                            }
                            others.insert({it.first.get(), newentry});
                            it.first = newentry;
                            result.insert(newentry.get());
                        } else {
                            it.first = ito->second;
                        }
                    }
                } else {
                    status.insert({it.first.get(), it.second});
                    result.insert(it.first.get());
                    if (it.second) {
                        toBeConjugate.insert(it.first);
                    }
                }
            }
            for (const auto& conjElem : toBeConjugate) {
                conjElem->conjugate();
            }
            for (auto& it : elem) {
                it.second = false;
            }
        }
        return result;
    }

    set<IContainer*> OuterContainer::getUniquePtrs(bool decoupleConjugates) {
        map<IContainer*, bool> status;
        set<SharedTensorData> toBeConjugate;
        map<IContainer*, SharedTensorData> others;
        set<IContainer*> result{};
        for (auto& elem : _data) {
            if (!decoupleConjugates) {
                for (auto& elem2 : elem) {
                    result.insert(elem2.first.get());
                }
            } else {
                for (auto& it : elem) {
                    auto its = status.find(it.first.get());
                    if (its != status.end()) {
                        if (its->second != it.second) {
                            auto ito = others.find(it.first.get());
                            if (ito == others.end()) {
                                SharedTensorData newentry{it.first->clone()};
                                if (it.second) {
                                    toBeConjugate.insert(newentry);
                                }
                                others.insert({it.first.get(), newentry});
                                it.first = newentry;
                                result.insert(newentry.get());
                            } else {
                                it.first = ito->second;
                            }
                        }
                    } else {
                        status.insert({it.first.get(), it.second});
                        result.insert(it.first.get());
                        if (it.second) {
                            toBeConjugate.insert(it.first);
                        }
                    }
                }
                for (const auto& conjElem : toBeConjugate) {
                    conjElem->conjugate();
                }
                for (auto& it : elem) {
                    it.second = false;
                }
            }
        }
        return result;
    }

    IContainer& OuterContainer::operator*=(double value) {
        if (_sharedData) {
            MSG_WARNING("Scalar multiplication with shared data: forcing cloning");
            auto tmp = this->clone();
            auto* newouter = static_cast<OuterContainer*>(tmp.get());
            _data.swap(newouter->_data);
            _sharedData = false;
        }
        if (!_data.empty()) {
            set<IContainer*> actualPtrs = getUniquePtrs();
            size_t n = _data.begin()->size();
            ElementType tmpvalue = pow(value, 1. / static_cast<double>(n));
            for (auto* elem : actualPtrs) {
                elem->operator*=(tmpvalue);
            }
        }
        return static_cast<IContainer&>(*this);
    }

    IContainer& OuterContainer::operator*=(ElementType value) {
        if (_sharedData) {
            MSG_WARNING("Scalar multiplication with shared data: forcing cloning");
            auto tmp = this->clone();
            auto* newouter = static_cast<OuterContainer*>(tmp.get());
            _data.swap(newouter->_data);
            _sharedData = false;
        }
        if (!_data.empty()) {
            std::set<IContainer*> actualPtrs = getUniquePtrs(!isZero(value.imag()));
            size_t n = _data.begin()->size();
            ElementType tmpvalue = pow(value, 1. / static_cast<double>(n));
            for (auto* elem : actualPtrs) {
                elem->operator*=(tmpvalue);
            }
        }
        return static_cast<IContainer&>(*this);
    }

    IContainer& OuterContainer::conjugate() {
        for (auto& elem : _data) {
            for (auto& elem2 : elem) {
                elem2.second = !elem2.second;
            }
        }
        _indexing.flipLabels();
        return static_cast<IContainer&>(*this);
    }

    OuterContainer::SerialType OuterContainer::write(flatbuffers::FlatBufferBuilder* msgwriter) const {
        vector<flatbuffers::Offset<Serial::FBSingleTensor>> tensors;
        vector<flatbuffers::Offset<Serial::FBTensorListEntry>> entries;
        entries.reserve(_data.size());
        map<IContainer*, uint16_t> stored;
        uint16_t nextId = 0;
        for (auto& elem : _data) {
            vector<uint16_t> values;
            vector<bool> conjStatus;
            values.reserve(elem.size());
            conjStatus.reserve(elem.size());
            for (auto& it : elem) {
                IContainer* ptr = it.first.get();
                auto itstore = stored.find(ptr);
                if (itstore == stored.end()) {
                    itstore = stored.insert(make_pair(ptr, nextId)).first;
                    ++nextId;
                }
                values.push_back(itstore->second);
                conjStatus.push_back(it.second);
            }
            auto serialvals = msgwriter->CreateVector(values);
            auto serialconj = msgwriter->CreateVector(conjStatus);
            Serial::FBTensorListEntryBuilder serialentry{*msgwriter};
            serialentry.add_idx(serialvals);
            serialentry.add_conj(serialconj);
            auto outEntry = serialentry.Finish();
            entries.push_back(outEntry);
        }
        tensors.resize(stored.size());
        for (auto& elem : stored) {
            auto tmp = elem.first->write(msgwriter);
            flatbuffers::Offset<Serial::FBSingleTensor> offt{tmp.first.o};
            tensors[elem.second] = offt;
        }
        auto serialvals = msgwriter->CreateVector(tensors);
        auto serialentries = msgwriter->CreateVector(entries);
        Serial::FBTensorListBuilder serialtensor{*msgwriter};
        serialtensor.add_vals(serialvals);
        serialtensor.add_entries(serialentries);
        auto out = serialtensor.Finish();
        return make_pair(out.Union(), Serial::FBTensorTypes::FBTensorList);
    }

    Log& OuterContainer::getLog() {
        return Log::getLog("Hammer.OuterContainer");
    }

    const BlockIndexing& OuterContainer::getIndexing() const {
        return _indexing;
    }

    void OuterContainer::swap(DataType values) {
        values.swap(_data);
    }

    void OuterContainer::swapIndexing(BlockIndexing values) {
        _indexing = std::move(values);
    }

    void OuterContainer::reserve(size_t numTerms) {
        size_t newCap = numTerms;
        if (numTerms < _data.capacity()) {
            newCap = max(numTerms * 10ul, 100ul);
        }
        _data.reserve(newCap);
    }

    void OuterContainer::addTerm(vector<pair<SharedTensorData, bool>> tensorsAndConjFlags) {
        if (tensorsAndConjFlags.size() == _indexing.numSubIndexing() || _indexing.numSubIndexing() == 0) {
            _data.push_back(tensorsAndConjFlags);
        } else {
            MSG_ERROR("Inconsistent term: " + to_string(tensorsAndConjFlags.size()) + " vs " +
                      to_string(_indexing.numSubIndexing()) + " with " + to_string(_data.size()) +
                      " terms. Not added.");
        }
    }

    bool OuterContainer::isOuterSquare() const {
        return _data.size() == 1 && _data[0].size() == 2 && _data[0][0].first.get() == _data[0][1].first.get() &&
               _data[0][0].second != _data[0][1].second;
    }

    bool OuterContainer::hasNaNs() const {
        for (auto& elem : _data) {
            for (auto& it : elem) {
                if (it.first->hasNaNs()) {
                    return true;
                }
            }
        }
        return false;
    }

    TensorData makeOuterSquare(const TensorData& base) {
        auto tmp = base->clone();
        auto* result = new OuterContainer{std::move(tmp)};
        return TensorData{static_cast<IContainer*>(result)};
    }

    TensorData makeOuterSquare(TensorData&& base) {
        auto* result = new OuterContainer{std::move(base)};
        return TensorData{static_cast<IContainer*>(result)};
    }

    TensorData combineTensors(const std::vector<TensorData>& base) {
        if (base.empty()) {
            return TensorData{new ScalarContainer()};
        }
        std::vector<TensorData> tmp;
        tmp.reserve(base.size());
        for (const auto& elem : base) {
            tmp.emplace_back(TensorData{elem->clone()});
        }
        return combineTensors(std::move(tmp));
    }

    TensorData combineTensors(std::vector<TensorData>&& base) {
        if (base.empty()) {
            return TensorData{new ScalarContainer()};
        }
        auto* result = new OuterContainer{std::move(base)};
        return TensorData{static_cast<IContainer*>(result)};
    }

    TensorData combineSharedTensors(std::vector<std::pair<SharedTensorData, bool>>&& data) {
        if (data.empty()) {
            return TensorData{new ScalarContainer()};
        }
        auto* result = new OuterContainer{std::move(data)};
        return TensorData{static_cast<IContainer*>(result)};
    }


} // namespace Hammer::MultiDimensional
