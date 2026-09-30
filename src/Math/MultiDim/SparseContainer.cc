///
/// @file  SparseContainer.cc
/// @brief Sparse tensor data container
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <boost/fusion/adapted/std_pair.hpp>

#include "Hammer/Math/MultiDim/SparseContainer.hh"
#include "Hammer/Exceptions.hh"
#include "Hammer/Math/MultiDimensional.fhh"
#include "Hammer/Tools/Utils.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Tools/HammerSerial.hh"
#include "Hammer/Tools/Logging.hh"

using namespace std;

namespace Hammer::MultiDimensional {

    SparseContainer::SparseContainer(const IndexList& dimensions, const LabelsList& labels)
        : _indexing{dimensions, labels} {
        clear();
    }

    SparseContainer::SparseContainer(const LabeledIndexing<AlignedIndexing>& indexing) : _indexing{indexing} {
        clear();
    }

    SparseContainer::SparseContainer(const Serial::FBSingleTensor* input) {
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
                _indexing = LabeledIndexing<AlignedIndexing>{tmpdims, tmplabels};
            }
            if (serialvals != nullptr && serialpos != nullptr) {
                for (unsigned int i = 0; i < serialvals->size(); ++i) {
                    _data[_indexing.posToAlignedPos(serialpos->Get(i))] =
                        complex<double>{serialvals->Get(i)->re(), serialvals->Get(i)->im()};
                }
            }
        } else {
            throw InitializationError("Sparse Container null buffer on read: have fun storming the castle!");
        }
    }

    SparseContainer::ElementType SparseContainer::value(const IndexList& indices) const {
        ASSERT(_indexing.checkValidIndices(indices));
        PositionType pos = _indexing.indicesToPos(indices);
        return getOrDefault(_data, pos, ElementType{});
    }

    SparseContainer::ElementType SparseContainer::value(IndexList::const_iterator first,
                                                        IndexList::const_iterator last) const {
        PositionType pos = _indexing.indicesToPos(first, last);
        return getOrDefault(_data, pos, ElementType{});
    }


    void SparseContainer::setValue(const IndexList& indices, ElementType value) {
        ASSERT(_indexing.checkValidIndices(indices));
        PositionType pos = _indexing.indicesToPos(indices);
        _data[pos] = value;
    }

    void SparseContainer::setValue(IndexList::const_iterator first, IndexList::const_iterator last, ElementType value) {
        ASSERT(_indexing.checkValidIndices(first, last));
        PositionType pos = _indexing.indicesToPos(first, last);
        _data[pos] = value;
    }


    SparseContainer::iterator SparseContainer::begin() {
        return _data.begin();
    }

    SparseContainer::const_iterator SparseContainer::begin() const {
        return _data.begin();
    }

    SparseContainer::iterator SparseContainer::end() {
        return _data.end();
    }

    SparseContainer::const_iterator SparseContainer::end() const {
        return _data.end();
    }

    SparseContainer::iterator SparseContainer::erase(SparseContainer::const_iterator first,
                                                     SparseContainer::const_iterator last) {
        return _data.erase(first, last);
    }

    SparseContainer::reference SparseContainer::operator[](PositionType pos) {
        return _data[pos];
    }

    size_t SparseContainer::rank() const {
        return _indexing.rank();
    }

    IndexList SparseContainer::dims() const {
        return _indexing.dims();
    }

    LabelsList SparseContainer::labels() const {
        return _indexing.labels();
    }

    size_t SparseContainer::numValues() const {
        return _indexing.numValues();
    }

    size_t SparseContainer::dataSize() const {
        return _data.size() * (sizeof(complex<double>) + sizeof(PositionType));
    }

    size_t SparseContainer::entrySize() const {
        return (sizeof(complex<double>) + sizeof(PositionType));
    }

    IndexType SparseContainer::labelToIndex(IndexLabel label) const {
        return _indexing.labelIndex(label);
    }


    IndexPairList SparseContainer::getSameLabelPairs(const IContainer& other, const UniqueLabelsList& indices) const {
        return _indexing.getSameLabelPairs(other.labels(), indices);
    }

    IndexPairList SparseContainer::getSpinLabelPairs() const {
        return _indexing.getOppositeLabelPairs(_indexing.spinIndices());
    }

    bool SparseContainer::isSameShape(const IContainer& other) const {
        return _indexing.isSameLabelShape(other.labels(), other.dims());
    }

    bool SparseContainer::canAddAt(const IContainer& subContainer, IndexLabel coord, IndexType position) const {
        return _indexing.canAddAt(subContainer.labels(), subContainer.dims(), coord, position);
    }

    SparseContainer::reference SparseContainer::element(const IndexList& coords) {
        ASSERT(_indexing.checkValidIndices(coords));
        PositionType pos = _indexing.indicesToPos(coords);
        return _data[pos];
    }

    SparseContainer::ElementType SparseContainer::element(const IndexList& coords) const {
        return value(coords);
    }


    SparseContainer::reference SparseContainer::element(IndexList::const_iterator start,
                                                        IndexList::const_iterator end) {
        ASSERT(_indexing.checkValidIndices(start, end));
        PositionType pos = _indexing.indicesToPos(start, end);
        return _data[pos];
    }

    SparseContainer::ElementType SparseContainer::element(IndexList::const_iterator start,
                                                          IndexList::const_iterator end) const {
        return value(start, end);
    }

    bool SparseContainer::compare(const IContainer& other) const {
        const auto* tmp = dynamic_cast<const SparseContainer*>(&other);
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
                     [](const pair<PositionType, ElementType>& a, const pair<PositionType, ElementType>& b) -> bool {
                         return (a.first == b.first) && isZero(a.second - b.second);
                     });
    }

    TensorData SparseContainer::clone() const {
        auto* t = new SparseContainer{*this};
        return TensorData{static_cast<IContainer*>(t)};
    }


    ISingleContainer::NonZeroIt SparseContainer::firstNonZero() const {
        return NonZeroIt{new ItAligned{_data.begin()}};
    }

    ISingleContainer::NonZeroIt SparseContainer::endNonZero() const {
        return NonZeroIt{new ItAligned{_data.end()}};
    }

    void SparseContainer::clear() {
        _data.clear();
    }

    IContainer& SparseContainer::operator*=(double value) {
        for (auto& elem : _data) {
            elem.second *= value;
        }
        return static_cast<IContainer&>(*this);
    }

    IContainer& SparseContainer::operator*=(ElementType value) {
        for (auto& elem : _data) {
            elem.second *= value;
        }
        return static_cast<IContainer&>(*this);
    }


    IContainer& SparseContainer::conjugate() {
        for (auto& elem : _data) {
            elem.second = conj(elem.second);
        }
        _indexing.flipLabels();
        return *this;
    }

    LabelSignature SparseContainer::getLabelSignature() const {
        return _indexing.getLabelSignature();
    }

    SparseContainer::SerialType SparseContainer::write(flatbuffers::FlatBufferBuilder* msgwriter) const {
        vector<Serial::FBComplex> datalist;
        vector<uint32_t> poslist;
        datalist.reserve(_data.size());
        poslist.reserve(_data.size());
        if (_data.size() > 4294967295ul) {
            MSG_ERROR("Error saving binary tensor: too many entries!");
        }
        // size_t count = 0ul;
        // size_t count2 = 0ul;
        for (auto& elem : _data) {
            if (!isZero(elem.second)) {
                // if (!isZero(_data[i].real())) count++;
                // if (!isZero(_data[i].imag())) count++;
                // count2++;
                auto resc = Serial::FBComplex{elem.second.real(), elem.second.imag()};
                datalist.push_back(resc);
                poslist.push_back(static_cast<uint32_t>(_indexing.alignedPosToPos(elem.first)));
            }
        }
        // std::cout << numValues() << ": " << 1.-(count*(8+4) + 4*4)/(16.*numValues()+4)
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
        serialtensor.add_sparse(true);
        serialtensor.add_vals(serialdata);
        serialtensor.add_pos(serialpos);
        serialtensor.add_labels(seriallabels);
        serialtensor.add_dims(serialdims);
        serialtensor.add_labelpadsindices(serialpadpos);
        serialtensor.add_labelpadsvalues(serialpadvals);
        auto out = serialtensor.Finish();
        return make_pair(out.Union(), Serial::FBTensorTypes::FBSingleTensor);
    }

    Log& SparseContainer::getLog() {
        return Log::getLog("Hammer.SparseContainer");
    }

    const LabeledIndexing<AlignedIndexing>& SparseContainer::getIndexing() const {
        return _indexing;
    }

    bool SparseContainer::hasNaNs() const {
        return any_of(_data.cbegin(), _data.cend(),
                      [](const pair<PositionType, ElementType>& elem) -> bool { return is_ieee754_nan(elem.second); });
    }

    TensorData makeEmptySparse(const IndexList& dimensions, const LabelsList& labels) {
        auto* result = new SparseContainer{dimensions, labels};
        return TensorData{static_cast<IContainer*>(result)};
    }

    TensorData makeEmptySparse(const LabeledIndexing<AlignedIndexing>& indexing) {
        auto* result = new SparseContainer{indexing};
        return TensorData{static_cast<IContainer*>(result)};
    }

    SparseContainer::ItAligned::ItAligned(DataType::const_iterator it) : _it{it} {
    }

    SparseContainer::ElementType SparseContainer::ItAligned::value() const {
        return _it->second;
    }

    PositionType SparseContainer::ItAligned::position() const {
        return _it->first;
    }

    void SparseContainer::ItAligned::next(int n) {
        advance(_it, n);
    }

    bool SparseContainer::ItAligned::isSame(const ItBase& other) const {
        return _it == static_cast<const ItAligned&>(other)._it;
    }

    bool SparseContainer::ItAligned::isAligned() const {
        return true;
    }

    ptrdiff_t SparseContainer::ItAligned::distanceFrom(const ItBase& other) const {
        return distance(static_cast<const ItAligned&>(other)._it, _it);
    }

} // namespace Hammer::MultiDimensional
