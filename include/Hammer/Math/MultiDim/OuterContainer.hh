///
/// @file  OuterContainer.hh
/// @brief (Sum of) Outer product tensor data container
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_MATH_MULTIDIMS_OUTERCONTAINER
#define HAMMER_MATH_MULTIDIMS_OUTERCONTAINER

#include <vector>
#include <functional>

#include "Hammer/Math/MultiDim/IContainer.hh"
#include "Hammer/Math/MultiDim/BlockIndexing.hh"

namespace Hammer {

    class Log;

    namespace Serial {

        struct FBTensorList;

    }

    namespace MultiDimensional {

        namespace Ops {

            class Sum;
            class Multiply;
            class Divide;
            class OuterSquare;
            class Dot;
            class DotOuterOptimizer;
            class Optimize;
            class Trace;
            class Convert;

        } // namespace Ops

        class OuterContainer final : public IContainer {
        private:

            using EntryType = OuterElemIterator::EntryType;
            using DataType = std::vector<EntryType>;

            OuterContainer() = default;

        public:

            explicit OuterContainer(TensorData left, TensorData right);
            explicit OuterContainer(TensorData toBeSquared, bool conjugate = true);
            explicit OuterContainer(std::vector<TensorData>&& group);
            // OuterContainer(std::vector<SharedTensorData>&& group, FlagList&& isConjugate = {});
            explicit OuterContainer(EntryType&& data);
            explicit OuterContainer(SharedTensorData toBeSquared, bool conjugate = true);
            explicit OuterContainer(const Serial::FBTensorList* input);

            OuterContainer(const OuterContainer& other);
            OuterContainer& operator=(const OuterContainer& other);

            OuterContainer(OuterContainer&& other) = default;
            OuterContainer& operator=(OuterContainer&& other) = default;

            ~OuterContainer() override = default;

            [[nodiscard]] ElementType value(const IndexList& indices) const;
            [[nodiscard]] ElementType value(IndexList::const_iterator first, IndexList::const_iterator last) const;
            [[nodiscard]] ElementType value(const std::vector<IndexList>& indices) const;

            [[nodiscard]] size_t numAddends() const;

            using iterator = DataType::iterator;
            using const_iterator = DataType::const_iterator;


            iterator begin();
            [[nodiscard]] const_iterator begin() const;

            iterator end();
            [[nodiscard]] const_iterator end() const;

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

        protected:

            /// @brief logging facility
            /// @return   stream to be used for logging
            static Log& getLog();

            std::set<IContainer*> getUniquePtrs(bool decoupleConjugates = false);
            std::set<IContainer*> getUniquePtrs(size_t pos, bool decoupleConjugates = false);
            void swapElement(const IContainer* oldContainer, TensorData newContainer);
            [[nodiscard]] bool shouldBeEvaluated() const;

        private:

            // all these classes are friends to be able to use getIndexing, swap, swapIndexing, addTerm and reserve
            // (not necessary, but cleaner)
            friend class Ops::Sum;
            friend class Ops::Multiply;
            friend class Ops::Divide;
            friend class Ops::OuterSquare;
            friend class Ops::Dot;
            friend class Ops::DotOuterOptimizer;
            friend class Ops::Optimize;
            friend class Ops::Trace;
            friend class Ops::Convert;

            const BlockIndexing& getIndexing() const;
            void swap(DataType values);
            void swapIndexing(BlockIndexing values);
            void reserve(size_t numTerms);
            void addTerm(std::vector<std::pair<SharedTensorData, bool>> tensorsAndConjFlags);
            [[nodiscard]] bool isDataShared() const;
            [[nodiscard]] bool isOuterSquare() const;

            mutable DataType _data;
            BlockIndexing _indexing;
            bool _sharedData;
            std::vector<std::function<ElementType(const IndexList&, IContainer*)>> _accessors;
        };

        TensorData makeOuterSquare(const TensorData& base);
        TensorData makeOuterSquare(TensorData&& base);
        TensorData combineTensors(const std::vector<TensorData>& base);
        TensorData combineTensors(std::vector<TensorData>&& base);

        TensorData combineSharedTensors(std::vector<std::pair<SharedTensorData, bool>>&& data);


    } // namespace MultiDimensional

} // namespace Hammer

#endif
