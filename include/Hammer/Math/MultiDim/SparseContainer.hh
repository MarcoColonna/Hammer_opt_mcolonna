///
/// @file  SparseContainer.hh
/// @brief Sparse tensor data container
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_MATH_MULTIDIM_SPARSECONTAINER
#define HAMMER_MATH_MULTIDIM_SPARSECONTAINER

#include <map>
#include <vector>

#include "Hammer/Math/MultiDim/ISingleContainer.hh"
#include "Hammer/Math/MultiDim/AlignedIndexing.hh"
#include "Hammer/Math/MultiDim/LabeledIndexing.hh"

namespace Hammer {

    class Log;

    namespace Serial {

        struct FBSingleTensor;

    }

    namespace MultiDimensional {

        namespace Ops {

            class Sum;
            class Multiply;
            class Divide;
            class Trace;
            class Dot;
            class AddAt;
            class Convert;
            class Compare;

        } // namespace Ops

        class SparseContainer final : public ISingleContainer {
        private:

            using DataType = std::map<PositionType, ElementType>;

        public:

            explicit SparseContainer(const IndexList& dimensions, const LabelsList& labels);
            explicit SparseContainer(const LabeledIndexing<AlignedIndexing>& indexing);
            explicit SparseContainer(const Serial::FBSingleTensor* input);

            [[nodiscard]] ElementType value(const IndexList& indices) const;
            [[nodiscard]] ElementType value(IndexList::const_iterator first, IndexList::const_iterator last) const;
            void setValue(const IndexList& indices, ElementType value = 0.);
            void setValue(IndexList::const_iterator first, IndexList::const_iterator last, ElementType value = 0.);

            using iterator = DataType::iterator;
            using const_iterator = DataType::const_iterator;


            iterator begin();
            [[nodiscard]] const_iterator begin() const;

            iterator end();
            [[nodiscard]] const_iterator end() const;

            reference operator[](PositionType pos);

            iterator erase(const_iterator first, const_iterator last);

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

            SerialType write(flatbuffers::FlatBufferBuilder* msgwriter) const override;

            [[nodiscard]] NonZeroIt firstNonZero() const override;
            [[nodiscard]] NonZeroIt endNonZero() const override;

            [[nodiscard]] bool hasNaNs() const override;

            [[nodiscard]] LabelSignature getLabelSignature() const override;

        protected:

            /// @brief logging facility
            /// @return   stream to be used for logging
            static Log& getLog();

        private:

            class ItAligned : public ItBase {
            public:

                explicit ItAligned(DataType::const_iterator it);

                [[nodiscard]] IContainer::ElementType value() const override;
                [[nodiscard]] PositionType position() const override;
                void next(int n = 1) override;
                [[nodiscard]] bool isSame(const ItBase& other) const override;
                [[nodiscard]] bool isAligned() const override;
                [[nodiscard]] ptrdiff_t distanceFrom(const ItBase& other) const override;

            private:

                friend class SparseContainer;
                DataType::const_iterator _it;
            };

            // all these classes are friends to be able to use getIndexing (not necessary, but cleaner)
            friend class Ops::Sum;
            friend class Ops::Multiply;
            friend class Ops::Divide;
            friend class Ops::Trace;
            friend class Ops::Dot;
            friend class Ops::AddAt;
            friend class Ops::Convert;
            friend class Ops::Compare;

            [[nodiscard]] const LabeledIndexing<AlignedIndexing>& getIndexing() const;

            mutable DataType _data;
            LabeledIndexing<AlignedIndexing> _indexing;
        };

        TensorData makeEmptySparse(const IndexList& dimensions, const LabelsList& labels);
        TensorData makeEmptySparse(const LabeledIndexing<AlignedIndexing>& indexing);

    } // namespace MultiDimensional

} // namespace Hammer


#endif
