///
/// @file  IContainer.hh
/// @brief Interface class for tensor container data structure
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_MATH_IMULTICONTAINER
#define HAMMER_MATH_IMULTICONTAINER

#include <complex>

#include "Hammer/Math/MultiDimensional.fhh"
#include "Hammer/Tools/HammerSerial.fhh"

namespace Hammer::MultiDimensional {

    class IContainer {

    public:

        virtual ~IContainer() = default;
        IContainer() = default;
        IContainer(const IContainer&) = default;
        IContainer(IContainer&&) = default;
        IContainer& operator=(const IContainer&) = default;
        IContainer& operator=(IContainer&&) = default;

        using ElementType = std::complex<double>;
        using reference = ElementType&;
        using const_reference = const ElementType&;

        [[nodiscard]] virtual size_t rank() const = 0;
        [[nodiscard]] virtual IndexList dims() const = 0;
        [[nodiscard]] virtual LabelsList labels() const = 0;
        [[nodiscard]] virtual size_t numValues() const = 0;
        [[nodiscard]] virtual size_t dataSize() const = 0;
        [[nodiscard]] virtual size_t entrySize() const = 0;
        [[nodiscard]] virtual IndexType labelToIndex(IndexLabel label) const = 0;

        [[nodiscard]] virtual IndexPairList getSameLabelPairs(const IContainer& other,
                                                              const UniqueLabelsList& indices) const = 0;
        [[nodiscard]] virtual IndexPairList getSpinLabelPairs() const = 0;

        [[nodiscard]] virtual bool isSameShape(const IContainer& other) const = 0;
        [[nodiscard]] virtual bool canAddAt(const IContainer& subContainer, IndexLabel coord,
                                            IndexType position) const = 0;

        virtual reference element(const IndexList& coords = {}) = 0;
        [[nodiscard]] virtual ElementType element(const IndexList& coords = {}) const = 0;

        virtual reference element(IndexList::const_iterator start, IndexList::const_iterator end) = 0;
        [[nodiscard]] virtual ElementType element(IndexList::const_iterator start,
                                                  IndexList::const_iterator end) const = 0;

        [[nodiscard]] virtual bool compare(const IContainer& other) const = 0;
        [[nodiscard]] virtual TensorData clone() const = 0;
        virtual void clear() = 0;

        virtual IContainer& operator*=(double value) = 0;
        virtual IContainer& operator*=(ElementType value) = 0;

        virtual IContainer& conjugate() = 0;

        [[nodiscard]] virtual bool hasNaNs() const = 0;

        [[nodiscard]] virtual LabelSignature getLabelSignature() const = 0;

        using SerialType = std::pair<flatbuffers::Offset<void>, Serial::FBTensorTypes>;

        virtual SerialType write(flatbuffers::FlatBufferBuilder* msgwriter) const = 0;
    };

} // namespace Hammer::MultiDimensional

#endif
