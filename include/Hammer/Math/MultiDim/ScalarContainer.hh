///
/// @file  ScalarContainer.hh
/// @brief Order-0 tensor data container
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_MATH_MULTIDIM_SCALARCONTAINER
#define HAMMER_MATH_MULTIDIM_SCALARCONTAINER

#include "Hammer/Math/MultiDim/IContainer.hh"

namespace Hammer {

    class Log;

    namespace Serial {

        struct FBComplex;

    }

    namespace MultiDimensional {

        class ScalarContainer final : public IContainer {
        public:

            ScalarContainer();
            explicit ScalarContainer(const Serial::FBComplex* input);

            [[nodiscard]] size_t rank() const override;
            [[nodiscard]] IndexList dims() const override;
            [[nodiscard]] LabelsList labels() const override;
            [[nodiscard]] size_t numValues() const override;
            [[nodiscard]] size_t dataSize() const override;
            [[nodiscard]] size_t entrySize() const override;
            [[nodiscard]] IndexType labelToIndex(IndexLabel label) const override;

            [[nodiscard]] IndexPairList getSameLabelPairs(const IContainer& other,
                                                          const UniqueLabelsList& indices) const override;
            [[nodiscard]] IndexPairList getSpinLabelPairs() const override;

            [[nodiscard]] bool isSameShape(const IContainer& other) const override;
            [[nodiscard]] bool canAddAt(const IContainer& subContainer, IndexLabel coord,
                                        IndexType position) const override;

            reference element(const IndexList& coords = {}) override;
            [[nodiscard]] ElementType element(const IndexList& coords = {}) const override;

            reference element(IndexList::const_iterator start, IndexList::const_iterator end) override;
            [[nodiscard]] ElementType element(IndexList::const_iterator start,
                                              IndexList::const_iterator end) const override;

            [[nodiscard]] bool compare(const IContainer& other) const override;
            [[nodiscard]] TensorData clone() const override;
            void clear() override;

            IContainer& operator*=(double value) override;
            IContainer& operator*=(ElementType value) override;

            IContainer& conjugate() override;

            [[nodiscard]] bool hasNaNs() const override;

            SerialType write(flatbuffers::FlatBufferBuilder* msgwriter) const override;

            [[nodiscard]] LabelSignature getLabelSignature() const override;

        private:

            ElementType _data;
        };


        TensorData makeEmptyScalar();
        TensorData makeScalar(std::complex<double> val);

    } // namespace MultiDimensional

} // namespace Hammer


#endif
