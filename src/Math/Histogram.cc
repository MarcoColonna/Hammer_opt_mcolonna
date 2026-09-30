///
/// @file  Histogram.cc
/// @brief Hammer histogram class
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/Math/Histogram.hh"
#include "Hammer/Exceptions.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"
#include "Hammer/Math/MultiDim/ScalarContainer.hh"
#include "Hammer/Math/Tensor.hh"
#include "Hammer/Tools/HammerSerial.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    Histogram::Bin::Bin(const Tensor& t, bool withErrors) : _weight{t}, _withErrors{withErrors} {
        _weight.clearData();
        if (withErrors) {
            _weight2 = t;
            _weight2.outerSquare();
            _weight2.clearData();
        }
    }

    Histogram::Bin::Bin(const Serial::FBTensor* w, const Serial::FBTensor* w2, size_t num) {
        _weight.read(w);
        if (w2 != nullptr) {
            _weight2.read(w2);
            _withErrors = true;
        } else {
            _weight2 = Tensor{};
            _withErrors = false;
        }
        _numEvents = num;
    }

    void Histogram::Bin::clear() {
        _weight = Tensor{};
        _weight2 = Tensor{};
        _numEvents = 0;
    }

    void Histogram::Bin::addWeight(const Tensor& t) {
        _weight += t;
        ++_numEvents;
        if (_withErrors) {
            _weight2 += outerSquare(t);
        }
    }

    const Tensor& Histogram::Bin::getWeight() const {
        return _weight;
    }

    const Tensor& Histogram::Bin::getSquaredWeight() const {
        return _weight2;
    }

    size_t Histogram::Bin::getNumEvents() const {
        return _numEvents;
    }

    Histogram::Bin& Histogram::Bin::operator+=(const Histogram::Bin& other) {
        if (other._numEvents != 0) {
            if (_numEvents == 0) {
                _weight = other._weight;
                if (other._weight2.rank() > 0) {
                    _weight2 = other._weight2;
                }
                _numEvents = other._numEvents;
            } else {
                _weight += other._weight;
                if (_weight2.rank() > 0 && other._weight2.rank() > 0) {
                    _weight2 += other._weight2;
                }
                _numEvents += other._numEvents;
            }
        }
        return (*this);
    }

    Histogram::Histogram() : _name{"Histogram"}, _defaultTensorName{"Bin"} {
    }

    Histogram::~Histogram() {
        clear();
    }

    void Histogram::setName(const string& name) {
        _name = name;
    }

    void Histogram::setDefinition(const HistogramDefinition& def) {
        _definition = &def;
    }

    Histogram::iterator Histogram::begin() {
        return _data.begin();
    }

    Histogram::const_iterator Histogram::begin() const {
        return _data.begin();
    }

    Histogram::iterator Histogram::end() {
        return _data.end();
    }

    Histogram::const_iterator Histogram::end() const {
        return _data.end();
    }

    Histogram::reference Histogram::operator[](PositionType pos) {
        return _data[pos];
    }

    Histogram::const_reference Histogram::operator[](PositionType pos) const {
        return _data[pos];
    }

    size_t Histogram::numValues() const {
        return _definition->getIndexing().numValues();
    }

    size_t Histogram::rank() const {
        return _definition->getIndexing().rank();
    }

    IndexType Histogram::dim(IndexType index) const {
        return _definition->getIndexing().dim(index);
    }

    const IndexList& Histogram::dims() const {
        return _definition->getIndexing().dims();
    }

    const string& Histogram::name() const {
        return _name;
    }

    void Histogram::setDefaultValue(const Tensor& value) {
        if (value.rank() > 0) {
            _defaultTensorIndexing = MD::LabeledIndexing<MD::AlignedIndexing>{value.dims(), value.labels()};
        }
        _defaultTensorName = value.name();
        _initialized = false;
    }

    void Histogram::setDefaultValue(const IndexList& defaultDims, const LabelsList& defaultLabels) {
        if (!defaultDims.empty()) {
            _defaultTensorIndexing = MD::LabeledIndexing<MD::AlignedIndexing>{defaultDims, defaultLabels};
            _defaultTensorName = "";
        }
        _initialized = false;
    }

    void Histogram::clear() {
        _data.clear();
        _definition = nullptr;
        _name = "";
    }

    Log& Histogram::getLog() {
        return Log::getLog("Hammer.Histogram");
    }

    void Histogram::fillBin(const Tensor& value, const IndexList& position) {
        try {
            element(position).addWeight(value);
        } catch (RangeError& e) {
            MSG_ERROR(string("Error setting element: ") + e.what() + string(" Element not set."));
        }
    }

    Histogram::Bin& Histogram::element(const IndexList& indices) {
        if (indices.size() < _definition->getIndexing().rank()) {
            throw RangeError("Not enough indices (" + to_string(indices.size()) + ") for tensor '" + _name +
                             "' of rank " + to_string(_definition->getIndexing().rank()) +
                             ": all your base are belong to us.");
        }
        return _data[_definition->getIndexing().indicesToPos(indices)];
    }

    const Histogram::Bin& Histogram::element(const IndexList& indices) const {
        if (indices.size() < _definition->getIndexing().rank()) {
            throw RangeError("Not enough indices (" + to_string(indices.size()) + ") for tensor '" + _name +
                             "' of rank " + to_string(_definition->getIndexing().rank()) +
                             ": all your base are belong to us.");
        }
        return _data[_definition->getIndexing().indicesToPos(indices)];
    }

    void Histogram::initHistogram() {
        if (!_initialized) {
            clearData();
            _initialized = true;
        }
    }

    void Histogram::clearData() {
        Tensor r{_defaultTensorName, MD::makeEmptyScalar()};
        if (_defaultTensorIndexing.rank() > 0) {
            r = Tensor{_defaultTensorName, MD::makeEmptySparse(_defaultTensorIndexing)};
        }
        vector<Bin> newvec(numValues(), Bin(r, _definition->shouldKeepErrors()));
        _data.swap(newvec);
    }

    unique_ptr<Serial::FBHistogramBuilder> Histogram::write(flatbuffers::FlatBufferBuilder* msgwriter) const {
        vector<flatbuffers::Offset<Serial::FBTensor>> weightbins;
        vector<flatbuffers::Offset<Serial::FBTensor>> weightbins2;
        vector<uint16_t> poslist;
        vector<uint64_t> nEvs;
        weightbins.reserve(_data.size());
        nEvs.reserve(_data.size());
        poslist.reserve(_data.size());
        if (_data.size() > 65535ul) {
            MSG_ERROR("Error saving binary histogram: too many entries!");
        }
        if (_definition->shouldKeepErrors()) {
            weightbins2.reserve(_data.size());
        }
        for (size_t i = 0; i < _data.size(); ++i) {
            if (_data[i].getNumEvents() > 0) {
                flatbuffers::Offset<Serial::FBTensor> val;
                _data[i].getWeight().write(msgwriter, &val);
                weightbins.push_back(val);
                nEvs.push_back(_data[i].getNumEvents());
                poslist.push_back(static_cast<uint16_t>(i));
                if (_definition->shouldKeepErrors()) {
                    flatbuffers::Offset<Serial::FBTensor> val2;
                    _data[i].getSquaredWeight().write(msgwriter, &val2);
                    weightbins2.push_back(val2);
                }
            }
        }
        auto serialW = msgwriter->CreateVector(weightbins);
        auto serialW2 = msgwriter->CreateVector(weightbins2);
        auto serialPos = msgwriter->CreateVector(poslist);
        auto serialNEvs = msgwriter->CreateVector(nEvs);
        auto serialname = msgwriter->CreateString(_name);

        auto serialHisto = std::make_unique<Serial::FBHistogramBuilder>(*msgwriter);
        serialHisto->add_w(serialW);
        serialHisto->add_w2(serialW2);
        serialHisto->add_n(serialNEvs);
        serialHisto->add_pos(serialPos);
        serialHisto->add_name(serialname);
        return serialHisto;
    }


    void Histogram::read(const Serial::FBHistogram* msgreader, const HistogramDefinition& def) {
        _definition = &def;
        if (msgreader != nullptr) {
            const auto* w = msgreader->w();
            if (w == nullptr) {
                MSG_ERROR("Missing w on read. You can't touch this.");
                return;
            }
            const auto* num = msgreader->n();
            if (num == nullptr) {
                MSG_ERROR("Missing n on read. You can't touch this.");
                return;
            }
            const auto* pos = msgreader->pos();
            if (pos == nullptr) {
                MSG_ERROR("Missing pos on read. You can't touch this.");
                return;
            }
            if (msgreader->name() == nullptr) {
                MSG_ERROR("Missing name on read. You can't touch this.");
                return;
            }
            // bool keepErrors = (w2->size() > 0);
            bool keepErrors = _definition->shouldKeepErrors();
            if (keepErrors && msgreader->w2() == nullptr) {
                MSG_ERROR("Missing w2 on read. You can't touch this.");
                return;
            }
            const auto* w2 = msgreader->w2();
            _name = msgreader->name()->c_str();
            _data.clear();
            size_t fullSize = _definition->getIndexing().numValues();
            _data.reserve(fullSize);
            size_t currentIndex = 0ul;
            bool firstread = true;
            Tensor empty;
            for (unsigned int i = 0; i < w->size(); ++i) {
                const auto* elemW = w->Get(i);
                auto elemN = num->Get(i);
                auto elemPos = pos->Get(i);
                const Serial::FBTensor* elemW2 = nullptr;
                if (keepErrors) {
                    elemW2 = w2->Get(i);
                }
                Bin bin(elemW, elemW2, elemN);
                if (firstread) {
                    _defaultTensorName = bin.getWeight().name();
                    LabelsList defaultLabels = bin.getWeight().labels();
                    IndexList defaultDimensions = bin.getWeight().dims();
                    if (!defaultDimensions.empty()) {
                        _defaultTensorIndexing =
                            MD::LabeledIndexing<MD::AlignedIndexing>{defaultDimensions, defaultLabels};
                        empty = Tensor{_defaultTensorName, MD::makeEmptySparse(_defaultTensorIndexing)};
                    } else {
                        empty = Tensor{_defaultTensorName, MD::makeEmptyScalar()};
                    }
                    firstread = false;
                }
                if (elemPos > currentIndex) {
                    _data.insert(_data.end(), elemPos - currentIndex, Bin(empty, keepErrors));
                    currentIndex = elemPos;
                }
                _data.push_back(bin);
                ++currentIndex;
            }
            size_t currentSize = _data.size();
            _data.insert(_data.end(), fullSize - currentSize, Bin(empty, keepErrors));
        }
    }

    Histogram& Histogram::operator+=(const Histogram& other) {
        if (_definition->getIndexing().isSameBinShape(other._definition->getIndexing())) {
            for (size_t i = 0; i < _data.size(); ++i) {
                _data[i] += other._data[i];
            }
            return (*this);
        }
        throw Error("Histogram dimension/sizes/binning do not match");
    }

    const LabelsList& Histogram::getLabelVector() const {
        return _defaultTensorIndexing.labels();
    }

    IOHistogram Histogram::evalToList(const Tensor& externalData) const {
        IOHistogram results;
        // map<size_t, size_t> necessary; //Deals with events weights containing multiple tensors with same indices,
        // auto labsd = _defaultTensorIndexing.labels();
        // multiset<IndexLabel> defaultTensorLabels(labsd.begin(), labsd.end());
        // for(size_t i=0; i< externalData.size(); ++i) {
        //     auto labsi = externalData[i].labels();
        //     set<IndexLabel> externalDataLabels(labsi.begin(), labsi.end());
        //     if(includes(defaultTensorLabels.begin(),defaultTensorLabels.end(),
        //                  externalDataLabels.begin(), externalDataLabels.end())) {
        //         size_t n = static_cast<size_t>(count(defaultTensorLabels.begin(), defaultTensorLabels.end(),
        //         *externalDataLabels.begin())); necessary.insert(make_pair(i,n));
        //     }
        // }
        for (const auto& elem : _data) {
            size_t num = elem.getNumEvents();
            double wVal = 0.;
            double w2Val = 0.;
            if (num > 0) {
                Tensor w = elem.getWeight();
                w.dot(externalData);
                if (w.rank() != 0) {
                    MSG_ERROR("Invalid rank, " + to_string(w.rank()) + ", at end of evaluation.");
                }
                wVal = w.element().real();
                if (_definition->shouldKeepErrors()) {
                    Tensor w2 = elem.getSquaredWeight();
                    auto externalData2 = outerSquare(externalData);
                    w2.dot(externalData2);
                    // w2.dot(externalData).dot(externalData);
                    if (w2.rank() != 0) {
                        MSG_ERROR("Invalid rank at end of evaluation");
                    }
                    w2Val = w2.element().real();
                }
            }
            results.push_back(BinContents{wVal, w2Val, num});
        }
        return results;
    }


    unique_ptr<Histogram> Histogram::collapseIntoNew(const string& newName, const HistogramDefinition& newDef,
                                                     const set<uint16_t>& collapsedIndexPositions) const {
        auto strides = _definition->getIndexing().getBinStrides(collapsedIndexPositions);
        auto result = std::make_unique<Histogram>();
        result->_definition = &newDef;
        result->_defaultTensorIndexing = _defaultTensorIndexing;
        result->_defaultTensorName = _defaultTensorName;
        result->_name = newName;
        result->initHistogram();
        for (size_t i = 0; i < _data.size(); ++i) {
            if (_data[i].getNumEvents() == 0) {
                continue;
            }
            auto pospairs = MD::BinnedIndexing<MD::SequentialIndexing>::splitPosition(i, strides);
            (result->_data)[pospairs.first] += _data[i];
        }
        return result;
    }

    unique_ptr<Histogram> makeHistogram(const string& name, const HistogramDefinition& def,
                                        const Tensor& defaultValue) {
        auto res = std::make_unique<Histogram>();
        res->setName(name);
        res->setDefinition(def);
        res->setDefaultValue(defaultValue);
        res->initHistogram();
        return res;
    }

    unique_ptr<Histogram> makeHistogram(const string& name, const HistogramDefinition& def,
                                        const IndexList& defaultTensorDims, const LabelsList& defaultTensorLabels) {
        auto res = std::make_unique<Histogram>();
        res->setName(name);
        res->setDefinition(def);
        res->setDefaultValue(defaultTensorDims, defaultTensorLabels);
        res->initHistogram();
        return res;
    }

    Histogram operator+(const Histogram& first, const Histogram& second) {
        Histogram result(first);
        result += second;
        return result;
    }

} // namespace Hammer
