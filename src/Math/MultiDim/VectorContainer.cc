///
/// @file  VectorContainer.cc
/// @brief Non-sparse tensor data container
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cstddef>

#include <boost/fusion/adapted/std_pair.hpp>

#include "Hammer/Math/MultiDim/VectorContainer.hh"
#include "Hammer/Exceptions.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Tools/HammerSerial.hh"
#include "Hammer/Tools/Logging.hh"

using namespace std;

namespace Hammer::MultiDimensional {

    VectorContainer::VectorContainer(const IndexList& dimensions, const LabelsList& labels)
        : _indexing{dimensions, labels} {
        clear();
    }

    VectorContainer::VectorContainer(const LabeledIndexing<SequentialIndexing>& indexing) : _indexing{indexing} {
        clear();
    }

    VectorContainer::VectorContainer(const Serial::FBSingleTensor* input) {
        if (input != nullptr) {
            const auto* serialvals = input->vals();
            const auto* serialpos = input->pos();
            const auto* serialdims = input->dims();
            const auto* seriallabels = input->labels();
            if (serialdims != nullptr && seriallabels != nullptr) {
                LabelsList tmplabels;
                tmplabels.reserve(seriallabels->size());
                IndexList tmpdims;
                tmpdims.reserve(serialdims->size());
                for (unsigned int i = 0; i < serialdims->size(); ++i) {
                    tmplabels.push_back(static_cast<IndexLabel>(seriallabels->Get(i)));
                    tmpdims.push_back(serialdims->Get(i));
                }
                if (flatbuffers::IsFieldPresent(input, Serial::FBSingleTensor::VT_LABELPADSINDICES) &&
                    flatbuffers::IsFieldPresent(input, Serial::FBSingleTensor::VT_LABELPADSVALUES)) {
                    const auto* serialpadindices = input->labelpadsindices();
                    const auto* serialpadvalues = input->labelpadsvalues();
                    for (unsigned int i = 0; i < serialpadindices->size(); ++i) {
                        uint32_t pad = serialpadvalues->Get(i);
                        auto j = serialpadindices->Get(i);
                        tmplabels[j] = specializeLabel(tmplabels[j], pad);
                    }
                }
                _indexing = LabeledIndexing<SequentialIndexing>{tmpdims, tmplabels};
            }
            clear();
            if (serialvals != nullptr && serialpos != nullptr) {
                for (unsigned int i = 0; i < serialvals->size(); ++i) {
                    _data[serialpos->Get(i)] = complex<double>{serialvals->Get(i)->re(), serialvals->Get(i)->im()};
                }
            }
        } else {
            throw InitializationError("Vector Container null buffer on read: have fun storming the castle!");
        }
    }

    VectorContainer::ElementType VectorContainer::value(const IndexList& indices) const {
        ASSERT(_indexing.checkValidIndices(indices));
        PositionType pos = _indexing.indicesToPos(indices);
        return _data[pos];
    }

    VectorContainer::ElementType VectorContainer::value(IndexList::const_iterator first,
                                                        IndexList::const_iterator last) const {
        PositionType pos = _indexing.indicesToPos(first, last);
        return _data[pos];
    }

    void VectorContainer::setValue(const IndexList& indices, ElementType value) {
        ASSERT(_indexing.checkValidIndices(indices));
        PositionType pos = _indexing.indicesToPos(indices);
        _data[pos] = value;
    }

    void VectorContainer::setValue(IndexList::const_iterator first, IndexList::const_iterator last, ElementType value) {
        ASSERT(_indexing.checkValidIndices(first, last));
        PositionType pos = _indexing.indicesToPos(first, last);
        _data[pos] = value;
    }

    void VectorContainer::clear() {
        _data.clear();
        DataType(_indexing.numValues()).swap(_data);
    }

    VectorContainer::iterator VectorContainer::begin() {
        return _data.begin();
    }

    VectorContainer::const_iterator VectorContainer::begin() const {
        return _data.begin();
    }

    VectorContainer::iterator VectorContainer::end() {
        return _data.end();
    }

    VectorContainer::const_iterator VectorContainer::end() const {
        return _data.end();
    }

    VectorContainer::reference VectorContainer::operator[](PositionType pos) {
        return _data[pos];
    }

    VectorContainer::const_reference VectorContainer::operator[](PositionType pos) const {
        return _data[pos];
    }

    size_t VectorContainer::rank() const {
        return _indexing.rank();
    }

    IndexList VectorContainer::dims() const {
        return _indexing.dims();
    }

    LabelsList VectorContainer::labels() const {
        return _indexing.labels();
    }

    size_t VectorContainer::numValues() const {
        return _indexing.numValues();
    }

    size_t VectorContainer::dataSize() const {
        return (static_cast<size_t>(
                    count_if(_data.begin(), _data.end(), [](complex<double> val) -> bool { return !isZero(val); })) *
                sizeof(complex<double>));
    }

    size_t VectorContainer::entrySize() const {
        return sizeof(complex<double>);
    }

    IndexType VectorContainer::labelToIndex(IndexLabel label) const {
        return _indexing.labelIndex(label);
    }


    IndexPairList VectorContainer::getSameLabelPairs(const IContainer& other, const UniqueLabelsList& indices) const {
        return _indexing.getSameLabelPairs(other.labels(), indices);
    }

    IndexPairList VectorContainer::getSpinLabelPairs() const {
        return _indexing.getOppositeLabelPairs(_indexing.spinIndices());
    }

    bool VectorContainer::isSameShape(const IContainer& other) const {
        return _indexing.isSameLabelShape(other.labels(), other.dims());
    }

    bool VectorContainer::canAddAt(const IContainer& subContainer, IndexLabel coord, IndexType position) const {
        return _indexing.canAddAt(subContainer.labels(), subContainer.dims(), coord, position);
    }

    VectorContainer::reference VectorContainer::element(const IndexList& coords) {
        ASSERT(_indexing.checkValidIndices(coords));
        PositionType pos = _indexing.indicesToPos(coords);
        return _data[pos];
    }

    VectorContainer::ElementType VectorContainer::element(const IndexList& coords) const {
        return value(coords);
    }


    VectorContainer::reference VectorContainer::element(IndexList::const_iterator start,
                                                        IndexList::const_iterator end) {
        ASSERT(_indexing.checkValidIndices(start, end));
        PositionType pos = _indexing.indicesToPos(start, end);
        return _data[pos];
    }

    VectorContainer::ElementType VectorContainer::element(IndexList::const_iterator start,
                                                          IndexList::const_iterator end) const {
        return value(start, end);
    }

    bool VectorContainer::compare(const IContainer& other) const {
        const auto* tmp = dynamic_cast<const VectorContainer*>(&other);
        if (tmp == nullptr) {
            return false;
        }
        bool dimcheck = _indexing.isSameLabelShape(tmp->_indexing);
        if (!dimcheck) {
            return false;
        }
        if (_data.size() != tmp->_data.size()) {
            return false;
        }
        return equal(_data.begin(), _data.end(), tmp->_data.begin(),
                     [](ElementType a, ElementType b) -> bool { return isZero(a - b); });
    }

    TensorData VectorContainer::clone() const {
        auto* t = new VectorContainer{*this};
        return TensorData{static_cast<IContainer*>(t)};
    }

    ISingleContainer::NonZeroIt VectorContainer::firstNonZero() const {
        auto it = NonZeroIt{new ItSequential{_data.begin(), _data.size()}};
        return it;
    }

    ISingleContainer::NonZeroIt VectorContainer::endNonZero() const {
        return NonZeroIt{new ItSequential{_data.end(), _data.size(), _data.size(),
                                          static_cast<PositionType>(dataSize() / sizeof(complex<double>))}};
    }

    IContainer& VectorContainer::operator*=(double value) {
        for (auto& elem : _data) {
            elem *= value;
        }
        return static_cast<IContainer&>(*this);
    }

    IContainer& VectorContainer::operator*=(ElementType value) {
        for (auto& elem : _data) {
            elem *= value;
        }
        return static_cast<IContainer&>(*this);
    }

    IContainer& VectorContainer::conjugate() {
        for (auto& elem : _data) {
            elem = conj(elem);
        }
        _indexing.flipLabels();
        return *this;
    }

    LabelSignature VectorContainer::getLabelSignature() const {
        return _indexing.getLabelSignature();
    }

    VectorContainer::SerialType VectorContainer::write(flatbuffers::FlatBufferBuilder* msgwriter) const {
        vector<Serial::FBComplex> datalist;
        vector<uint32_t> poslist;
        datalist.reserve(_data.size());
        poslist.reserve(_data.size());
        if (_data.size() > 4294967295ul) {
            MSG_ERROR("Error saving binary tensor: too many entries!");
        }
        // size_t count = 0ul;
        // size_t count2 = 0ul;
        for (size_t i = 0; i < _data.size(); ++i) {
            if (!isZero(_data[i])) {
                // if (!isZero(_data[i].real())) count++;
                // if (!isZero(_data[i].imag())) count++;
                // count2++;
                auto resc = Serial::FBComplex{_data[i].real(), _data[i].imag()};
                datalist.push_back(resc);
                poslist.push_back(static_cast<uint32_t>(i));
            }
        }
        // cout << numValues() << ": " << 1.-(count*(8+4) + 4*4)/(16.*numValues()+4)
        //     << " vs " << 1.-(count2*(16+4)+4*2)/(16.*numValues()+4) << '\n';
        auto serialdata = msgwriter->CreateVectorOfStructs(datalist.data(), datalist.size());
        auto serialpos = msgwriter->CreateVector(poslist);
        IndexList outdims;
        for (auto elem : _indexing.dims()) {
            outdims.push_back(static_cast<uint16_t>(elem));
        }
        vector<int16_t> outlabels;
        vector<uint32_t> padvals;
        vector<uint8_t> padpos;
        for (uint8_t i = 0; i < static_cast<uint8_t>(_indexing.labels().size()); ++i) {
            auto elem = _indexing.labels()[i];
            if (isSpecializedIndex(elem)) {
                auto tmplabel = getSpecializedLabelComponents(elem);
                outlabels.push_back(static_cast<int16_t>(tmplabel.second));
                padpos.push_back(i);
                padvals.push_back(tmplabel.first);
            } else {
                outlabels.push_back(static_cast<int16_t>(elem));
            }
        }
        auto serialdims = msgwriter->CreateVector(outdims);
        auto seriallabels = msgwriter->CreateVector(outlabels);
        auto serialpadpos = msgwriter->CreateVector(padpos);
        auto serialpadvals = msgwriter->CreateVector(padvals);
        Serial::FBSingleTensorBuilder serialtensor{*msgwriter};
        serialtensor.add_sparse(false);
        serialtensor.add_vals(serialdata);
        serialtensor.add_pos(serialpos);
        serialtensor.add_labels(seriallabels);
        serialtensor.add_dims(serialdims);
        serialtensor.add_labelpadsindices(serialpadpos);
        serialtensor.add_labelpadsvalues(serialpadvals);
        auto out = serialtensor.Finish();
        return make_pair(out.Union(), Serial::FBTensorTypes::FBSingleTensor);
    }

    Log& VectorContainer::getLog() {
        return Log::getLog("Hammer.VectorContainer");
    }

    const LabeledIndexing<SequentialIndexing>& VectorContainer::getIndexing() const {
        return _indexing;
    }

    void VectorContainer::swap(vector<complex<double>>& values) {
        ASSERT(_data.size() == values.size());
        values.swap(_data);
    }

    TensorData makeEmptyVector(const IndexList& dimensions, const LabelsList& labels) {
        auto* result = new VectorContainer{dimensions, labels};
        return TensorData{static_cast<IContainer*>(result)};
    }

    TensorData makeEmptyVector(const LabeledIndexing<SequentialIndexing>& indexing) {
        auto* result = new VectorContainer{indexing};
        return TensorData{static_cast<IContainer*>(result)};
    }

    TensorData makeVector(const IndexList& dimensions, const LabelsList& labels, vector<complex<double>> values) {
        auto* result = new VectorContainer{dimensions, labels};
        result->swap(values);
        return TensorData{static_cast<IContainer*>(result)};
    }

    VectorContainer::ItSequential::ItSequential(DataType::const_iterator it, PositionType maxPos, PositionType pos,
                                                PositionType nonZeroPos)
        : _it{it}, _maxPosition{maxPos}, _nonZeroPos{nonZeroPos}, _position{pos} {
        if (_position == 0 && _nonZeroPos == 0 && _maxPosition > 0 && isZero(*_it)) {
            while (_position < _maxPosition && isZero(*_it)) {
                ++_it;
                ++_position;
            }
        }
    }

    VectorContainer::ElementType VectorContainer::ItSequential::value() const {
        return *_it;
    }

    bool VectorContainer::hasNaNs() const {
        return any_of(_data.cbegin(), _data.cend(),
                      [](const ElementType& elem) -> bool { return is_ieee754_nan(elem); });
    }

    PositionType VectorContainer::ItSequential::position() const {
        return _position;
    }

    void VectorContainer::ItSequential::next(int n) {
        if (n > 0) {
            for (size_t i = 0; i < static_cast<size_t>(n) && _position < _maxPosition; ++i) {
                ++_it;
                ++_position;
                while (_position < _maxPosition && isZero(*_it)) {
                    ++_it;
                    ++_position;
                }
                ++_nonZeroPos;
            }
        } else if (n < 0) {
            for (size_t i = 0; i < static_cast<size_t>(-n) && _nonZeroPos > 0; ++i) {
                --_it;
                --_position;
                while (isZero(*_it)) {
                    --_it;
                    --_position;
                }
                --_nonZeroPos;
            }
        }
    }

    bool VectorContainer::ItSequential::isSame(const ItBase& other) const {
        return _it == static_cast<const ItSequential&>(other)._it;
    }

    bool VectorContainer::ItSequential::isAligned() const {
        return false;
    }

    ptrdiff_t VectorContainer::ItSequential::distanceFrom(const ItBase& other) const {
        return static_cast<ptrdiff_t>(_nonZeroPos - static_cast<const ItSequential&>(other)._nonZeroPos);
    }


} // namespace Hammer::MultiDimensional
