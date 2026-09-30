///
/// @file  Tensor.cc
/// @brief Hammer tensor class
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/Math/Tensor.hh"
#include "Hammer/Exceptions.hh"
#include "Hammer/Tools/HammerSerial.hh"
#include "Hammer/Math/MultiDim/IContainer.hh"
#include "Hammer/Math/MultiDim/ScalarContainer.hh"
#include "Hammer/Math/MultiDim/OuterContainer.hh"
#include "Hammer/Math/MultiDim/Operations.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    Tensor::Tensor() : _name{"Tensor"}, _data{MD::makeEmptyScalar()} {
    }

    Tensor::Tensor(string name) : _name{std::move(name)}, _data{MD::makeEmptyScalar()} {
    }

    Tensor::Tensor(string name, MultiDimensional::TensorData container)
        : _name{std::move(name)}, _data{std::move(container)} {
    }

    Tensor::Tensor(std::string name, std::vector<std::pair<MultiDimensional::SharedTensorData, bool>>&& data)
        : _name{std::move(name)}, _data{MD::combineSharedTensors(std::move(data))} {
    }

    Tensor::Tensor(const Tensor& other) : _name{other._name}, _data{other._data->clone()} {
    }

    Tensor& Tensor::operator=(const Tensor& other) {
        if (this == &other) {
            return *this;
        }
        _name = other._name;
        _data = other._data->clone();
        return *this;
    }

    Tensor::Tensor(Tensor&& other) noexcept : _name{std::move(other._name)}, _data{std::move(other._data)} {
    }

    Tensor& Tensor::operator=(Tensor&& other) noexcept {
        _name = std::move(other._name);
        _data = std::move(other._data);
        return *this;
    }

    Tensor::~Tensor() {
        _data.reset();
    }

    complex<double>& Tensor::element(const IndexList& indices) {
        return _data->element(indices);
    }

    complex<double> Tensor::element(const IndexList& indices) const {
        return const_cast<const MD::IContainer*>(_data.get())->element(indices);
    }

    size_t Tensor::rank() const {
        return _data->rank();
    }

    IndexList Tensor::dims() const {
        return _data->dims();
    }

    LabelsList Tensor::labels() const {
        return _data->labels();
    }

    bool Tensor::hasWCLabels() const {
        auto labs = labels();
        return find_if(labs.begin(), labs.end(), [&](IndexLabel& lab) -> bool {
                   return (lab > WC_INDEX_START && lab < WC_INDEX_END);
               }) != labs.end();
    }

    bool Tensor::hasFFLabels() const {
        auto labs = labels();
        return find_if(labs.begin(), labs.end(), [&](IndexLabel& lab) -> bool {
                   return (lab > FF_INDEX_START && lab < FF_INDEX_END);
               }) != labs.end();
    }

    bool Tensor::hasFFVarLabels() const {
        auto labs = labels();
        return find_if(labs.begin(), labs.end(), [&](IndexLabel& lab) -> bool {
                   return (lab > FF_VAR_INDEX_START && lab < FF_VAR_INDEX_END);
               }) != labs.end();
    }

    const string& Tensor::name() const {
        return _name;
    }

    bool Tensor::isEqualTo(const Tensor& other) const {
        return _data->compare(*other._data);
    }

    bool Tensor::isSameLabelShape(const Tensor& other) const {
        return _data->isSameShape(*other._data);
    }

    Log& Tensor::getLog() const {
        return Log::getLog("Hammer.Tensor." + _name);
    }

    Tensor& Tensor::dot(const Tensor& other, const UniqueLabelsList& indices) {
        if (other.rank() == 0) {
            operator*=(other.element());
            return *this;
        }
        if (rank() == 0) {
            complex<double> val = element();
            _data = other._data->clone();
            operator*=(val);
            return *this;
        }
        UniqueLabelsList tmpIndices;
        if (indices.empty()) {
            auto labs = labels();
            tmpIndices.insert(labs.begin(), labs.end());
        } else {
            tmpIndices.insert(indices.begin(), indices.end());
        }
        IndexPairList positions = _data->getSameLabelPairs(*other._data, tmpIndices);
        if (!positions.empty()) {
            _data = MD::calcDot(std::move(_data), *other._data, positions);
        } else {
            MD::TensorData tmp{new MD::OuterContainer{std::move(_data), other._data->clone()}};
            _data.swap(tmp);
        }
        return *this;
    }

    Tensor& Tensor::spinSum() {
        if (rank() == 0) {
            return *this;
        }
        IndexPairList positions = _data->getSpinLabelPairs();
        if (!positions.empty()) {
            _data = MD::calcTrace(std::move(_data), positions);
        }
        return *this;
    }

    Tensor& Tensor::spinAverage() {
        size_t oldDim = _data->numValues();
        spinSum();
        size_t newDim = _data->numValues();
        return operator*=((static_cast<double>(newDim)) / (static_cast<double>(oldDim)));
    }

    Tensor& Tensor::outerSquare() {
        if (rank() == 0) {
            element() = element() * conj(element());
            return *this;
        }
        _data = MD::calcSquare(std::move(_data));
        return *this;
    }

    Tensor& Tensor::toVector() {
        _data = MD::toVector(std::move(_data));
        return *this;
    }

    Tensor& Tensor::toSparse() {
        _data = MD::toSparse(std::move(_data));
        return *this;
    }

    Tensor& Tensor::operator*=(double val) {
        _data->operator*=(val);
        return *this;
    }

    Tensor& Tensor::operator*=(const complex<double> val) {
        _data->operator*=(val);
        return *this;
    }


    Tensor& Tensor::operator+=(const Tensor& other) {
        if (rank() == other.rank()) {
            if (rank() == 0) {
                element() += other.element();
                return (*this);
            }
            if (_data->isSameShape(*other._data)) {
                _data = MD::sum(std::move(_data), *other._data);
                return (*this);
            }
        }
        throw Error(
            "Tensor rank/sizes do not match. Whose gonna do it you, you Lieutenant Weinberg? I have a greater "
            "responsibility than you can possibly fathom.");
    }

    Tensor& Tensor::elementMultiplyBy(const Tensor& other) {
        if (rank() == other.rank()) {
            if (rank() == 0) {
                element() *= other.element();
                return (*this);
            }
            if (_data->isSameShape(*other._data)) {
                _data = MD::elementMultiply(std::move(_data), *other._data);
                return (*this);
            }
        }
        throw Error("Tensor rank/sizes do not match");
    }

    Tensor& Tensor::elementDivideBy(const Tensor& other) {
        if (rank() == other.rank()) {
            if (rank() == 0) {
                element() /= other.element();
                return (*this);
            }
            if (_data->isSameShape(*other._data)) {
                _data = MD::elementDivide(std::move(_data), *other._data);
                return (*this);
            }
        }
        throw Error(
            "Tensor rank/sizes do not match. Whose gonna do it you, you Lieutenant Weinberg? I have a greater "
            "responsibility than you can possibly fathom.");
    }

    void Tensor::clearData() {
        _data->clear();
    }

    void Tensor::write(flatbuffers::FlatBufferBuilder* msgwriter, flatbuffers::Offset<Serial::FBTensor>* msg) const {
        auto serialname = msgwriter->CreateString(_name);
        _data = MD::reOptimize(std::move(_data));
        auto serialvalue = _data->write(msgwriter);
        Serial::FBTensorBuilder serialTensor{*msgwriter};
        serialTensor.add_name(serialname);
        serialTensor.add_data_type(serialvalue.second);
        serialTensor.add_data(serialvalue.first);
        *msg = serialTensor.Finish();
    }


    void Tensor::read(const Serial::FBTensor* msgreader) {
        if (msgreader != nullptr) {
            _data = MD::read(msgreader);
        }
    }

    // Tensor&
    // Tensor::fillTensorElements(const vector<pair<vector<size_t>, complex<double>>>& values) {
    //     for (auto& elem : values) {
    //         element(elem.first) = elem.second;
    //     }
    //     return *this;
    // }

    Tensor& Tensor::addAt(const Tensor& t, IndexLabel coord, IndexType position) {
        if (rank() != 0 && _data->canAddAt(*t._data, coord, position)) {
            IndexType index = _data->labelToIndex(coord);
            _data = MD::addAt(std::move(_data), *t._data, index, position);
            return *this;
        }
        throw IndexLabelError("Unable to addAt");
    }

    bool Tensor::hasNaNs() const {
        return _data->hasNaNs();
    }


    Tensor dot(const Tensor& first, const Tensor& second, const set<IndexLabel>& indices) {
        Tensor result(first);
        result.dot(second, indices);
        return result;
    }

    Tensor spinSum(const Tensor& first) {
        Tensor result(first);
        result.spinSum();
        return result;
    }

    Tensor spinAverage(const Tensor& first) {
        Tensor result(first);
        result.spinAverage();
        return result;
    }

    Tensor operator*(const Tensor& first, double val) {
        Tensor result(first);
        result *= val;
        return result;
    }

    Tensor operator*(double val, const Tensor& first) {
        Tensor result(first);
        result *= val;
        return result;
    }

    Tensor operator*(const complex<double> val, const Tensor& first) {
        Tensor result(first);
        result *= val;
        return result;
    }

    Tensor operator*(const Tensor& first, const complex<double> val) {
        Tensor result(first);
        result *= val;
        return result;
    }

    Tensor operator+(const Tensor& first, const Tensor& second) {
        Tensor result(first);
        result += second;
        return result;
    }

    Tensor outerSquare(const Tensor& first) {
        Tensor result(first);
        result.outerSquare();
        return result;
    }

    Tensor elementMultiply(const Tensor& first, const Tensor& second) {
        Tensor result(first);
        result.elementMultiplyBy(second);
        return result;
    }

    Tensor elementDivide(const Tensor& first, const Tensor& second) {
        Tensor result(first);
        result.elementDivideBy(second);
        return result;
    }

} // namespace Hammer
