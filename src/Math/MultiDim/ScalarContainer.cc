///
/// @file  ScalarContainer.cc
/// @brief Order-0 tensor data container
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/Math/MultiDim/ScalarContainer.hh"
#include "Hammer/Exceptions.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Tools/HammerSerial.hh"

using namespace std;

namespace Hammer::MultiDimensional {

    ScalarContainer::ScalarContainer() : _data{0.} {
    }

    ScalarContainer::ScalarContainer(const Serial::FBComplex* input) : _data{0.} {
        if (input != nullptr) {
            _data = complex<double>{input->re(), input->im()};
        } else {
            throw InitializationError("Scalar Container null buffer on read: have fun storming the castle!");
        }
    }

    size_t ScalarContainer::rank() const {
        return 0ul;
    }

    IndexList ScalarContainer::dims() const {
        return {};
    }

    LabelsList ScalarContainer::labels() const {
        return {};
    }

    size_t ScalarContainer::numValues() const {
        return 1ul;
    }

    size_t ScalarContainer::dataSize() const {
        return sizeof(complex<double>);
    }

    size_t ScalarContainer::entrySize() const {
        return sizeof(complex<double>);
    }

    IndexType ScalarContainer::labelToIndex(IndexLabel /*label*/) const {
        return 0ul;
    }


    IndexPairList ScalarContainer::getSameLabelPairs(const IContainer& /*other*/,
                                                     const UniqueLabelsList& /*indices*/) const {
        return {};
    }

    IndexPairList ScalarContainer::getSpinLabelPairs() const {
        return {};
    }


    bool ScalarContainer::isSameShape(const IContainer& other) const {
        return (other.rank() == 0);
    }

    bool ScalarContainer::canAddAt(const IContainer& /*subContainer*/, IndexLabel /*coord*/,
                                   IndexType /*position*/) const {
        return false;
    }


    ScalarContainer::reference ScalarContainer::element(const IndexList& coords) {
        if (coords.empty()) {
            return _data;
        }
        throw RangeError("Index for multidimensional object out of range (" + to_string(coords.size()) +
                         " indices vs " + "rank 0): " + "all your base are belong to us.");
    }

    ScalarContainer::ElementType ScalarContainer::element(const IndexList& coords) const {
        if (coords.empty()) {
            return _data;
        }
        throw RangeError("Index for multidimensional object out of range (" + to_string(coords.size()) +
                         " indices vs " + "rank 0): " + "all your base are belong to us.");
    }


    ScalarContainer::reference ScalarContainer::element(IndexList::const_iterator start,
                                                        IndexList::const_iterator end) {
        if (start == end) {
            return _data;
        }
        throw RangeError("Index for multidimensional object out of range (" + to_string(distance(start, end)) +
                         " indices vs " + "rank 0): " + "all your base are belong to us.");
    }

    ScalarContainer::ElementType ScalarContainer::element(IndexList::const_iterator start,
                                                          IndexList::const_iterator end) const {
        if (start == end) {
            return _data;
        }
        throw RangeError("Index for multidimensional object out of range (" + to_string(distance(start, end)) +
                         " indices vs " + "rank 0): " + "all your base are belong to us.");
    }

    bool ScalarContainer::compare(const IContainer& other) const {
        if (other.rank() != 0) {
            return false;
        }
        return isZero(_data - other.element());
    }

    TensorData ScalarContainer::clone() const {
        return TensorData{static_cast<IContainer*>(new ScalarContainer{*this})};
    }

    void ScalarContainer::clear() {
        _data = 0.;
    }


    IContainer& ScalarContainer::operator*=(double value) {
        _data *= value;
        return *this;
    }

    IContainer& ScalarContainer::operator*=(ElementType value) {
        _data *= value;
        return *this;
    }

    IContainer& ScalarContainer::conjugate() {
        _data = conj(_data);
        return *this;
    }

    LabelSignature ScalarContainer::getLabelSignature() const {
        return {0, false}; // rank 0 tensor has no labels, so we return a dummy signature
    }


    ScalarContainer::SerialType ScalarContainer::write(flatbuffers::FlatBufferBuilder* msgwriter) const {
        auto resc = Serial::FBComplex{_data.real(), _data.imag()};
        auto out = msgwriter->CreateStruct(resc);
        return make_pair(out.Union(), Serial::FBTensorTypes::FBComplex);
    }

    bool ScalarContainer::hasNaNs() const {
        return is_ieee754_nan(_data);
    }

    TensorData makeEmptyScalar() {
        auto* result = new ScalarContainer{};
        return TensorData{static_cast<IContainer*>(result)};
    }

    TensorData makeScalar(complex<double> value) {
        auto* result = new ScalarContainer{};
        result->element() = value;
        return TensorData{static_cast<IContainer*>(result)};
    }


} // namespace Hammer::MultiDimensional
