///
/// @file  Dot.cc
/// @brief Tensor dot product algorithm
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-

#include <iterator>
#include <set>
#include <numeric>
#include <type_traits>
#include <tuple>

#include <boost/functional/hash.hpp>

#include "Hammer/Math/MultiDim/Ops/Dot.hh"
#include "Hammer/Math/MultiDim/IContainer.hh"
#include "Hammer/Math/MultiDim/Ops/Sum.hh"
#include "Hammer/Math/MultiDim/ScalarContainer.hh"
#include "Hammer/Math/MultiDim/VectorContainer.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"
#include "Hammer/Math/MultiDim/OuterContainer.hh"
#include "Hammer/Math/MultiDim/BruteForceIterator.hh"
#include "Hammer/Math/MultiDim/BlockIndexing.hh"
#include "Hammer/Exceptions.hh"
#include "Hammer/Math/MultiDimensional.fhh"
#include "Hammer/Tools/Utils.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Math/MultiDim/OperationDefs.hh"
#include "Hammer/Math/MultiDim/Ops/DotOuterOptimizer.hh"

using namespace std;

namespace Hammer::MultiDimensional {

    using VTensor = VectorContainer;
    using STensor = SparseContainer;
    using OTensor = OuterContainer;
    using Base = IContainer;

    namespace Ops {

        Dot::Dot(IndexPairList indices, pair<bool, bool> shouldHC) noexcept // NOLINT(bugprone-exception-escape)
            : _indices{std::move(indices)}, _hc{std::move(shouldHC)} {
            for (auto& elem : _indices) {
                _idxLeft.insert(elem.first);
                _idxRight.insert(elem.second);
            }
        }

        Base* Dot::operator()(VTensor& a, const VTensor& b) {
            auto newdimlabs = getNewIndexLabels(a, b);
            auto stridesA = a.getIndexing().getInnerOuterStrides(_indices, b.getIndexing().strides());
            if (newdimlabs.first.empty()) {
                auto newscal = makeEmptyScalar();
                for (size_t i = 0; i < a.numValues(); ++i) {
                    if (isZero(a[i])) {
                        continue;
                    }
                    Base::ElementType firstTerm = _hc.first ? conj(a[i]) : a[i];
                    auto pospairs = LabeledIndexing<SequentialIndexing>::splitPosition(i, stridesA);
                    Base::ElementType secondTerm = _hc.second ? conj(b[pospairs.second]) : b[pospairs.second];
                    newscal->element({}) += firstTerm * secondTerm;
                }
                return newscal.release();
            }
            auto newvect = makeEmptyVector(newdimlabs.first, newdimlabs.second);
            auto* result = static_cast<VTensor*>(newvect.release());
            if (b.rank() == 1) {
                for (size_t i = 0; i < a.numValues(); ++i) {
                    if (isZero(a[i])) {
                        continue;
                    }
                    Base::ElementType firstTerm = _hc.first ? conj(a[i]) : a[i];
                    auto pospairs = LabeledIndexing<SequentialIndexing>::splitPosition(i, stridesA);
                    Base::ElementType secondTerm = _hc.second ? conj(b[pospairs.second]) : b[pospairs.second];
                    (*result)[pospairs.first] += firstTerm * secondTerm;
                }
            } else {
                PositionType reduced = b.getIndexing().reducedNumValues(_indices);
                auto stridesB = b.getIndexing().getOuterStrides2nd(_indices);
                for (size_t i = 0; i < a.numValues(); ++i) {
                    if (isZero(a[i])) {
                        continue;
                    }
                    Base::ElementType firstTerm = _hc.first ? conj(a[i]) : a[i];
                    auto pospairs = LabeledIndexing<SequentialIndexing>::splitPosition(i, stridesA);
                    for (size_t j = 0; j < reduced; ++j) {
                        auto posB = LabeledIndexing<SequentialIndexing>::build2ndPosition(j, pospairs.second, stridesB);
                        Base::ElementType secondTerm = _hc.second ? conj(b[posB]) : b[posB];
                        (*result)[pospairs.first * reduced + j] += firstTerm * secondTerm;
                    }
                }
            }
            return static_cast<Base*>(result);
        }

        Base* Dot::operator()(STensor& a, const STensor& b) {
            auto newdimlabs = getNewIndexLabels(a, b);
            auto leftinfo = a.getIndexing().processShifts(_indices, IndexPairMember::Left);
            auto rightinfo = b.getIndexing().processShifts(_indices, IndexPairMember::Right);
            IndexList inners((a.dims().size() + b.dims().size() - newdimlabs.first.size()) / 2);
            FlagList innerAdds(inners.size(), false);
            if (newdimlabs.first.empty()) {
                auto newscal = makeEmptyScalar();
                for (const auto& elemL : a) {
                    (void) a.getIndexing().splitPosition(elemL.first, get<0>(leftinfo), get<1>(leftinfo), inners,
                                                         innerAdds);
                    Base::ElementType firstTerm = _hc.first ? conj(elemL.second) : elemL.second;
                    for (const auto& elemR : b) {
                        auto tmpRight = b.getIndexing().splitPosition(elemR.first, get<0>(rightinfo), get<1>(rightinfo),
                                                                      inners, innerAdds, true);
                        if (tmpRight == numeric_limits<size_t>::max()) {
                            continue;
                        }
                        Base::ElementType secondTerm = _hc.second ? conj(elemR.second) : elemR.second;
                        newscal->element({}) += firstTerm * secondTerm;
                    }
                }
                return newscal.release();
            }
            auto newsparse = makeEmptySparse(newdimlabs.first, newdimlabs.second);
            auto* result = static_cast<STensor*>(newsparse.release());
            for (const auto& elemL : a) {
                auto tmpLeft =
                    a.getIndexing().splitPosition(elemL.first, get<0>(leftinfo), get<1>(leftinfo), inners, innerAdds);
                Base::ElementType firstTerm = _hc.first ? conj(elemL.second) : elemL.second;
                for (const auto& elemR : b) {
                    auto tmpRight = b.getIndexing().splitPosition(elemR.first, get<0>(rightinfo), get<1>(rightinfo),
                                                                  inners, innerAdds, true);
                    if (tmpRight == numeric_limits<size_t>::max()) {
                        continue;
                    }
                    Base::ElementType secondTerm = _hc.second ? conj(elemR.second) : elemR.second;
                    (*result)[tmpLeft * get<2>(rightinfo) + tmpRight] += firstTerm * secondTerm;
                }
            }
            return static_cast<Base*>(result);
        }

        Base* Dot::operator()(STensor& a, const VTensor& b) {
            auto newdimlabs = getNewIndexLabels(a, b);
            IndexList inners((a.dims().size() + b.dims().size() - newdimlabs.first.size()) / 2);
            FlagList innerAdds(inners.size(), false);
            auto leftinfo = a.getIndexing().processShifts(_indices, IndexPairMember::Left);
            if (newdimlabs.first.empty()) {
                auto newscal = makeEmptyScalar();
                if (b.rank() > 1) {
                    for (const auto& elemL : a) {
                        (void) a.getIndexing().splitPosition(elemL.first, get<0>(leftinfo), get<1>(leftinfo), inners,
                                                             innerAdds);
                        Base::ElementType firstTerm = _hc.first ? conj(elemL.second) : elemL.second;
                        Base::ElementType secondTerm = _hc.second ? conj(b.element(inners)) : b.element(inners);
                        newscal->element({}) += firstTerm * secondTerm;
                    }
                } else {
                    for (const auto& elemL : a) {
                        (void) a.getIndexing().splitPosition(elemL.first, get<0>(leftinfo), get<1>(leftinfo), inners,
                                                             innerAdds);
                        Base::ElementType firstTerm = _hc.first ? conj(elemL.second) : elemL.second;
                        Base::ElementType secondTerm = _hc.second ? conj(b[inners[0]]) : b[inners[0]];
                        newscal->element({}) += firstTerm * secondTerm;
                    }
                }
                return newscal.release();
            }
            auto newsparse = makeEmptySparse(newdimlabs.first, newdimlabs.second);
            auto* result = static_cast<STensor*>(newsparse.release());
            if (b.rank() == _indices.size()) {
                if (b.rank() > 1) {
                    for (const auto& elemL : a) {
                        auto tmpLeft = a.getIndexing().splitPosition(elemL.first, get<0>(leftinfo), get<1>(leftinfo),
                                                                     inners, innerAdds);
                        Base::ElementType firstTerm = _hc.first ? conj(elemL.second) : elemL.second;
                        Base::ElementType secondTerm = _hc.second ? conj(b.element(inners)) : b.element(inners);
                        (*result)[tmpLeft] += firstTerm * secondTerm;
                    }
                } else {
                    for (const auto& elemL : a) {
                        auto tmpLeft = a.getIndexing().splitPosition(elemL.first, get<0>(leftinfo), get<1>(leftinfo),
                                                                     inners, innerAdds);
                        Base::ElementType firstTerm = _hc.first ? conj(elemL.second) : elemL.second;
                        Base::ElementType secondTerm = _hc.second ? conj(b[inners[0]]) : b[inners[0]];
                        (*result)[tmpLeft] += firstTerm * secondTerm;
                    }
                }
            } else {
                auto itBe = b.endNonZero();
                LabeledIndexing<AlignedIndexing> fakeB{b.dims(), b.labels()};
                auto rightinfo = fakeB.processShifts(_indices, IndexPairMember::Right);
                for (const auto& elemL : a) {
                    auto tmpLeft = a.getIndexing().splitPosition(elemL.first, get<0>(leftinfo), get<1>(leftinfo),
                                                                 inners, innerAdds);
                    Base::ElementType firstTerm = _hc.first ? conj(elemL.second) : elemL.second;
                    auto itB = b.firstNonZero();
                    for (; *itB != *itBe; itB->next()) {
                        PositionType alPos = fakeB.posToAlignedPos(itB->position());
                        auto tmpRight =
                            fakeB.splitPosition(alPos, get<0>(rightinfo), get<1>(rightinfo), inners, innerAdds, true);
                        if (tmpRight == numeric_limits<size_t>::max()) {
                            continue;
                        }
                        Base::ElementType secondTerm = _hc.second ? conj(itB->value()) : itB->value();
                        (*result)[tmpLeft * get<2>(rightinfo) + tmpRight] += firstTerm * secondTerm;
                    }
                }
            }
            return static_cast<Base*>(result);
        }

        using SameDotEntryType = pair<map<IndexType, IndexType>, map<IndexType, IndexType>>;
        using SameDotType = vector<map<IndexType, SameDotEntryType>>;


        inline pair<size_t, size_t> hashTempTensors(
            const tuple<OuterElemIterator::EntryType, OuterElemIterator::EntryType, IndexPairList>& dotData) {
            size_t hashA = 0;
            size_t hashB = 0;
            for (const auto& elem : get<0>(dotData)) {
                boost::hash_combine(hashA, elem.first.get());
                boost::hash_combine(hashB, elem.first.get());
                boost::hash_combine(hashA, elem.second);
                boost::hash_combine(hashB, !elem.second);
            }
            for (const auto& elem : get<1>(dotData)) {
                boost::hash_combine(hashA, elem.first.get());
                boost::hash_combine(hashB, elem.first.get());
                boost::hash_combine(hashA, elem.second);
                boost::hash_combine(hashB, !elem.second);
            }
            for (const auto& elem : get<2>(dotData)) {
                boost::hash_combine(hashA, elem.first);
                boost::hash_combine(hashA, elem.second);
                boost::hash_combine(hashB, elem.first);
                boost::hash_combine(hashB, elem.second);
            }
            return {hashA, hashB};
        }

        static pair<bool, bool> isSameDot(const OuterElemIterator::EntryType& a, const OuterElemIterator::EntryType& b,
                                  const DotGroupType& info, const DotGroupType& infoOther) {
            UNUSED(a);
            UNUSED(b);
            UNUSED(info);
            UNUSED(infoOther);
            return {false, false};
        }


        Base* Dot::operator()(OTensor& a, const OTensor& b) {
                // get the chunks
                DotGroupList chunks = partitionContractions(a.getIndexing(), b.getIndexing());
                TensorData fullResult;
                OTensor* oResult = nullptr;
                auto leftinfo = a.getIndexing().processShifts(chunks, IndexPairMember::Left);
                auto rightinfo = b.getIndexing().processShifts(chunks, IndexPairMember::Right);
                for(auto& elemA: a) {
                    for(auto& elemB: b) {
                        // check for repetitions
                        PositionPairList multiplicities(leftinfo.size());
                        FlagList used(leftinfo.size(), false);
                        for(size_t i=0; i < leftinfo.size(); ++i) {
                            if(used[i]) continue;
                            used[i] = true;
                            size_t count = 1;
                            size_t countHc = 0;
                            for (size_t j = i + 1; j < leftinfo.size(); ++j) {
                                auto tmp = isSameDot(elemA, elemB, chunks[i+1], chunks[j+1]);
                                if(tmp.first) {
                                    if(tmp.second) {
                                        ++countHc;
                                    }
                                    else {
                                        ++count;
                                    }
                                    used[j]=true;
                                }
                            }
                            multiplicities[i] = {count, countHc};
                        }
                        Base::ElementType currentWeight = 1.;
                        OuterElemIterator::EntryType currentTerm;
                        for (size_t i = 0; i < leftinfo.size(); ++i) {
                            if (multiplicities[i].first + multiplicities[i].second == 0)
                                continue;
                            // do the dot
                            OuterElemIterator::EntryType leftTensors;
                            leftTensors.reserve(get<0>(chunks[i + 1]).size());
                            transform(get<0>(chunks[i + 1]).begin(), get<0>(chunks[i + 1]).end(), back_inserter(leftTensors),
                                      [&](IndexType idx) -> const pair<SharedTensorData, bool>& { return elemA[idx]; });
                            OuterElemIterator::EntryType rightTensors;
                            rightTensors.reserve(get<1>(chunks[i + 1]).size());
                            transform(get<1>(chunks[i + 1]).begin(), get<1>(chunks[i + 1]).end(), back_inserter(rightTensors),
                                      [&](IndexType idx) -> const pair<SharedTensorData, bool>& { return elemB[idx]; });
                            IndexList inners(get<2>(chunks[i + 1]).size());
                            FlagList innerAdds(inners.size(), false);
                            auto newdimlabs = getNewIndexLabels(a.getIndexing(), b.getIndexing(), chunks[i + 1]);
                            OuterElemIterator itA{leftTensors};
                            OuterElemIterator itAEnd = itA.end();
                            size_t totalRankB = accumulate(rightTensors.begin(), rightTensors.end(), 0ul, [](PositionType tot, const pair<SharedTensorData, bool>& elem) -> PositionType { return tot + elem.first->rank(); });
                            if(totalRankB == inners.size()) {
                                IndexList::iterator itP1, itP2;
                                if (newdimlabs.first.size() == 0) {
                                    Base::ElementType newscal;
                                    for (; itA != itAEnd; ++itA) {
                                        itP1 = inners.begin();
                                        a.getIndexing().splitPosition(itA, chunks[i + 1], get<0>(leftinfo[i]),
                                                                      get<1>(leftinfo[i]), inners, innerAdds);
                                        Base::ElementType firstTerm = _hc.first ? conj(*itA) : *itA;
                                        for(auto& entry: rightTensors) {
                                            itP2 = itP1 + static_cast<ptrdiff_t>(entry.first->rank());
                                            Base::ElementType secondTerm = (!_hc.second != !entry.second)
                                                                                ? conj(entry.first->element(itP1, itP2))
                                                                                : entry.first->element(itP1, itP2);
                                            firstTerm *= secondTerm;
                                            if(isZero(secondTerm)) {
                                                break;
                                            }
                                            itP1 = itP2;
                                        }
                                        newscal += firstTerm;
                                    }
                                    currentWeight *=
                                        pow(newscal, multiplicities[i].first) * pow(conj(newscal), multiplicities[i].second);
                                } else {
                                    auto newsparse = makeEmptySparse(newdimlabs.first, newdimlabs.second);
                                    STensor* result = static_cast<STensor*>(newsparse.get());
                                    for (; itA != itAEnd; ++itA) {
                                        itP1 = inners.begin();
                                        PositionType tmpLeft =
                                            a.getIndexing().splitPosition(itA, chunks[i + 1], get<0>(leftinfo[i]),
                                                                          get<1>(leftinfo[i]), inners, innerAdds);
                                        Base::ElementType firstTerm = _hc.first ? conj(*itA) : *itA;
                                        for (auto& entry : rightTensors) {
                                            itP2 = itP1 + static_cast<ptrdiff_t>(entry.first->rank());
                                            Base::ElementType secondTerm = (!_hc.second != !entry.second)
                                                                               ? conj(entry.first->element(itP1, itP2))
                                                                               : entry.first->element(itP1, itP2);
                                            firstTerm *= secondTerm;
                                            if (isZero(secondTerm)) {
                                                break;
                                            }
                                            itP1 = itP2;
                                        }
                                        (*result)[tmpLeft] += firstTerm;
                                    }
                                    SharedTensorData tmpShared{newsparse.release()};
                                    currentTerm.insert(currentTerm.end(), multiplicities[i].first, {tmpShared, false});
                                    currentTerm.insert(currentTerm.end(), multiplicities[i].second, {tmpShared, true});
                                }
                            }
                            else {
                                OuterElemIterator itBEnd = OuterElemIterator{rightTensors}.end();
                                if (newdimlabs.first.size() == 0) {
                                    Base::ElementType newscal;
                                    for (; itA != itAEnd; ++itA) {
                                        a.getIndexing().splitPosition(itA, chunks[i + 1], get<0>(leftinfo[i]),
                                                                      get<1>(leftinfo[i]), inners, innerAdds);
                                        Base::ElementType firstTerm = _hc.first ? conj(*itA) : *itA;
                                        OuterElemIterator itB{rightTensors};
                                        for (; itB != itBEnd; ++itB) {
                                            PositionType tmpRight = b.getIndexing().splitPosition(
                                                itB, chunks[i + 1], get<0>(rightinfo[i]), get<1>(rightinfo[i]), inners,
                                                innerAdds, true);
                                            if (tmpRight == numeric_limits<size_t>::max())
                                                continue;
                                            Base::ElementType secondTerm = _hc.second ? conj(*itB) : *itB;
                                            newscal += firstTerm * secondTerm;
                                        }
                                    }
                                    currentWeight *= pow(newscal, multiplicities[i].first) *
                                                     pow(conj(newscal), multiplicities[i].second);
                                } else {
                                    auto newsparse = makeEmptySparse(newdimlabs.first, newdimlabs.second);
                                    STensor* result = static_cast<STensor*>(newsparse.get());
                                    for (; itA != itAEnd; ++itA) {
                                        PositionType tmpLeft =
                                            a.getIndexing().splitPosition(itA, chunks[i + 1], get<0>(leftinfo[i]),
                                                                          get<1>(leftinfo[i]), inners, innerAdds);
                                        Base::ElementType firstTerm = _hc.first ? conj(*itA) : *itA;
                                        OuterElemIterator itB{rightTensors};
                                        for (; itB != itBEnd; ++itB) {
                                            PositionType tmpRight = b.getIndexing().splitPosition(
                                                itB, chunks[i + 1], get<0>(rightinfo[i]), get<1>(rightinfo[i]), inners,
                                                innerAdds, true);
                                            if (tmpRight == numeric_limits<size_t>::max())
                                                continue;
                                            Base::ElementType secondTerm = _hc.second ? conj(*itB) : *itB;
                                            (*result)[tmpLeft * get<2>(rightinfo[i]) + tmpRight] +=
                                                firstTerm * secondTerm;
                                        }
                                    }
                                    SharedTensorData tmpShared{newsparse.release()};
                                    currentTerm.insert(currentTerm.end(), multiplicities[i].first, {tmpShared, false});
                                    currentTerm.insert(currentTerm.end(), multiplicities[i].second, {tmpShared, true});
                                }
                            }
                        }
                        // now add those untouched
                        for(auto elem: get<0>(chunks[0])) {
                            currentTerm.insert(currentTerm.end(), {elemA[elem].first, !elemA[elem].second != !_hc.first});
                        }
                        for(auto elem: get<1>(chunks[0])) {
                            currentTerm.insert(currentTerm.end(), {elemB[elem].first, !elemB[elem].second != !_hc.second});
                        }
                        if(currentTerm.size() == 0) {
                            if(fullResult.get() == nullptr) {
                                fullResult = makeScalar(currentWeight.real());
                            }
                            else {
                                fullResult->element({}) += (currentWeight.real());
                            }
                        }
                        else if(currentTerm.size() == 1) {
                            TensorData out = currentTerm[0].first->clone();
                            if(currentTerm[0].second) {
                                out->conjugate();
                            }
                            if (!isZero(currentWeight - 1.)) {
                                out->operator*=(currentWeight);
                            }
                            if(fullResult.get() == nullptr) {
                                fullResult = move(out);
                            }
                            else {
                                Ops::Sum summer{};
                                fullResult = calc2(move(fullResult), *out, summer, "sum_outerdot");
                            }
                        }
                        else {
                            if (!isZero(currentWeight - 1.)) {
                                auto temp = currentTerm.back();
                                currentTerm.pop_back();
                                auto candNew = temp.first->clone();
                                candNew->operator*=(temp.second ? conj(currentWeight) : currentWeight);
                                currentTerm.push_back({SharedTensorData{candNew.release()}, temp.second});
                            }
                            if (oResult != nullptr) {
                                oResult->addTerm(currentTerm);
                            } else {
                                fullResult = combineSharedTensors(move(currentTerm));
                                oResult = static_cast<OTensor*>(fullResult.get());
                            }
                        }
                    }
                }
                return static_cast<Base*>(fullResult.release());
            }


        Base* Dot::operator()(OTensor& a, const STensor& b) { // NOLINT(readability-make-member-function-const)
            DotGroupList chunks = partitionContractions(a.getIndexing(), b.getIndexing());
            ASSERT(get<1>(chunks[0]).empty());
            ASSERT(get<2>(chunks[0]).empty());
            ASSERT(get<1>(chunks[1]).empty());
            TensorData fullResult;
            OTensor* oResult = nullptr;
            auto leftinfo = a.getIndexing().processShifts(chunks, IndexPairMember::Left);
            auto rightinfo = b.getIndexing().processShifts(_indices, IndexPairMember::Right);
            ASSERT(leftinfo.size() == 1);
            for (auto& elemA : a) {
                Base::ElementType currentWeight = 1.;
                OuterElemIterator::EntryType currentTerm;
                OuterElemIterator::EntryType leftTensors;
                leftTensors.reserve(get<0>(chunks[1]).size());
                transform(get<0>(chunks[1]).begin(), get<0>(chunks[1]).end(), back_inserter(leftTensors),
                          [&](IndexType idx) -> const pair<SharedTensorData, bool>& { return elemA[idx]; });
                IndexList inners(get<2>(chunks[1]).size());
                FlagList innerAdds(inners.size(), false);
                auto newdimlabs = getNewIndexLabels(a.getIndexing(), b.getIndexing(), chunks[1]);
                OuterElemIterator itA{leftTensors};
                OuterElemIterator itAEnd = itA.end();
                if (newdimlabs.first.empty()) {
                    Base::ElementType newscal;
                    for (; itA != itAEnd; ++itA) {
                        (void) a.getIndexing().splitPosition(itA, chunks[1], get<0>(leftinfo[0]), get<1>(leftinfo[0]),
                                                             inners, innerAdds);
                        Base::ElementType firstTerm = _hc.first ? conj(*itA) : *itA;
                        for (const auto& elemR : b) {
                            auto tmpRight = b.getIndexing().splitPosition(elemR.first, get<0>(rightinfo),
                                                                          get<1>(rightinfo), inners, innerAdds, true);
                            if (tmpRight == numeric_limits<size_t>::max()) {
                                continue;
                            }
                            Base::ElementType secondTerm = _hc.second ? conj(elemR.second) : elemR.second;
                            newscal += firstTerm * secondTerm;
                        }
                    }
                    currentWeight *= newscal;
                } else {
                    auto newsparse = makeEmptySparse(newdimlabs.first, newdimlabs.second);
                    auto* result = static_cast<STensor*>(newsparse.get());
                    for (; itA != itAEnd; ++itA) {
                        PositionType tmpLeft = a.getIndexing().splitPosition(itA, chunks[1], get<0>(leftinfo[0]),
                                                                             get<1>(leftinfo[0]), inners, innerAdds);
                        Base::ElementType firstTerm = _hc.first ? conj(*itA) : *itA;
                        for (const auto& elemR : b) {
                            auto tmpRight = b.getIndexing().splitPosition(elemR.first, get<0>(rightinfo),
                                                                          get<1>(rightinfo), inners, innerAdds, true);
                            if (tmpRight == numeric_limits<size_t>::max()) {
                                continue;
                            }
                            Base::ElementType secondTerm = _hc.second ? conj(elemR.second) : elemR.second;
                            (*result)[tmpLeft * get<2>(rightinfo) + tmpRight] += firstTerm * secondTerm;
                        }
                    }
                    SharedTensorData tmpShared{newsparse.release()};
                    currentTerm.emplace_back(tmpShared, false);
                }
                for (auto elem : get<0>(chunks[0])) {
                    currentTerm.insert(currentTerm.end(), {elemA[elem].first, !elemA[elem].second != !_hc.first});
                }
                if (currentTerm.empty()) {
                    if (fullResult.get() == nullptr) {
                        fullResult = makeScalar(currentWeight); // real?
                    } else {
                        fullResult->element({}) += currentWeight; // real?
                    }
                } else if (currentTerm.size() == 1) {
                    TensorData out = currentTerm[0].first->clone();
                    if (currentTerm[0].second) {
                        out->conjugate();
                    }
                    if (!isZero(currentWeight - 1.)) {
                        out->operator*=(currentWeight);
                    }
                    if (fullResult.get() == nullptr) {
                        fullResult = std::move(out);
                    } else {
                        Ops::Sum summer{};
                        fullResult = calc2(std::move(fullResult), *out, summer, "sum_outerdot");
                    }
                } else {
                    if (!isZero(currentWeight - 1.)) {
                        auto temp = currentTerm.back();
                        currentTerm.pop_back();
                        auto candNew = temp.first->clone();
                        candNew->operator*=(temp.second ? conj(currentWeight) : currentWeight);
                        currentTerm.emplace_back(SharedTensorData{candNew.release()}, temp.second);
                    }
                    if (oResult != nullptr) {
                        oResult->addTerm(currentTerm);
                    } else {
                        fullResult = combineSharedTensors(std::move(currentTerm));
                        oResult = static_cast<OTensor*>(fullResult.get());
                    }
                }
            }
            return static_cast<Base*>(fullResult.release());
        }


        Base* Dot::operator()(OTensor& a, const Base& b) {
            // inverse ordering: start dotting from rightmost ensures compatibility
            //                   with getNewIndexLabels() for fully collapsed outers
            map<IndexType, IndexPairList, greater<>> edges;
            for (auto& elem : _indices) {
                auto keyvalA = a.getIndexing().getElementIndex(elem.first);
                edges[keyvalA.first].emplace_back(keyvalA.second, elem.second);
            }
            IndexType tmpshift = 0ul;
            for (auto& elem : edges) {
                for (auto& elem2 : elem.second) {
                    elem2.second = static_cast<IndexType>(elem2.second + tmpshift);
                }
                tmpshift = static_cast<IndexType>(tmpshift + a.getIndexing().getSubIndexing(elem.first).rank() -
                                                  elem.second.size());
            }
            if (edges.size() == a.getIndexing().numSubIndexing()) {
                auto newdimlabs = getNewIndexLabels(a, b);
                TensorData result;
                if (newdimlabs.first.empty()) {
                    result = makeEmptyScalar();
                } else {
                    result = makeEmptySparse(newdimlabs.first, newdimlabs.second);
                }
                // Outer is going to be collapsed
                SharedTensorData current_entry;
                Ops::Sum summer{};
                for (auto& elem : a) {
                    bool first = true;
                    for (auto& elem2 : edges) {
                        // dot the subtensor
                        Ops::Dot dotter{elem2.second, {elem[elem2.first].second, false}};
                        if (!first) {
                            current_entry = calc2(elem[elem2.first].first, *current_entry, dotter, "dot_outerdot");
                        } else {
                            first = false;
                            current_entry = calc2(elem[elem2.first].first, b, dotter, "dot_outerdot");
                        }
                    }
                    if (result->rank() == 0) {
                        result->element({}) += current_entry->element({});
                    } else {
                        result = calc2(std::move(result), *current_entry, summer, "sum_outerdot");
                    }
                }
                return result.release();
            }
            OuterContainer::DataType newdata;
            newdata.reserve(a.numAddends());
            vector<pair<SharedTensorData, bool>> current_entries(a.getIndexing().numSubIndexing() - edges.size() + 1);
            for (auto& elem : a) {
                size_t dotIdx = current_entries.size();
                size_t curInIdx = elem.size() - 1;
                size_t curOutIdx = current_entries.size() - 1;
                for (const auto& elem2 : edges) {
                    if (elem2.first != curInIdx) {
                        // copy the pass throughs
                        for (; curInIdx != elem2.first; --curInIdx, --curOutIdx) {
                            current_entries[curOutIdx] = elem[curInIdx];
                        }
                    }
                    // dot the subtensor
                    Ops::Dot dotter{elem2.second, {elem[curInIdx].second, false}};
                    if (dotIdx < current_entries.size()) {
                        current_entries[dotIdx].first =
                            calc2(elem[curInIdx].first, *(current_entries[dotIdx].first), dotter, "dot_outerdot");
                    } else {
                        dotIdx = curOutIdx--;
                        current_entries[dotIdx].first = calc2(elem[curInIdx].first, b, dotter, "dot_outerdot");
                        current_entries[dotIdx].second = false;
                    }
                }
                newdata.push_back(current_entries);
            }
            vector<IndexList> dimlist;
            vector<LabelsList> lablist;
            for (const auto& elem : newdata[0]) {
                dimlist.push_back(elem.first->dims());
                lablist.push_back(elem.first->labels());
            }
            a.swap(newdata);
            a.swapIndexing(BlockIndexing{dimlist, lablist});
            // Rebuild _accessors to match the new sub-tensor layout after partial collapse.
            a._accessors.clear();
            {
                size_t start = 0ul;
                size_t finish = 0ul;
                for (const auto& entry : a._data[0]) {
                    finish += entry.first->rank();
                    a._accessors.emplace_back(start, finish);
                    //a._accessors.emplace_back(
                    //    [start, finish](const IndexList& list, const IContainer* item) -> Base::ElementType {
                    //        return item->element(list.begin() + static_cast<ptrdiff_t>(start),
                    //                             list.begin() + static_cast<ptrdiff_t>(finish));
                    //    });
                    start = finish;
                }
            }
            return static_cast<Base*>(&a);
        }

        Base* Dot::operator()(STensor& a, const OTensor& b) {
            DotGroupList chunks = partitionContractions(a.getIndexing(), b.getIndexing());
            TensorData fullResult;
            OTensor* oResult = nullptr;
            auto leftinfo = a.getIndexing().processShifts(_indices, IndexPairMember::Left);
            auto rightinfo = b.getIndexing().processShifts(chunks, IndexPairMember::Right);
            ASSERT(rightinfo.size() == 1);
            // check if it's a star topology or boomerang
            bool should1ns = true;
            bool shouldboomerang = true;
            for (auto& elem : get<1>(chunks[1])) {
                should1ns &= b.begin()->at(elem).first->rank() == 1;
                shouldboomerang &= b.begin()->at(elem).first->rank() == 2;
            }
            shouldboomerang &= get<1>(chunks[1]).size() == 2 && get<2>(chunks[1]).size() == 2;
            PositionList sortedOrder(_indices.size());
            IndexPairList contractions{};
            if (should1ns) {
                iota(sortedOrder.begin(), sortedOrder.end(), 0);
                sort(sortedOrder.begin(), sortedOrder.end(),
                     [&](size_t x, size_t y) { return _indices[x].first < _indices[y].first; });
                contractions.reserve(_indices.size());
                for (size_t k = 0; k < sortedOrder.size(); ++k) {
                    contractions.emplace_back(_indices[sortedOrder[k]].first, static_cast<IndexType>(k));
                }
            }
            for (const auto& elemB : b) {
                Base::ElementType currentWeight = 1.;
                OuterElemIterator::EntryType currentTerm;
                if (should1ns) { // star topology
                    OuterElemIterator::EntryType rightTensors;
                    rightTensors.reserve(sortedOrder.size());
                    for (auto i : sortedOrder) {
                        rightTensors.push_back(elemB[get<1>(chunks[1])[i]]);
                    }
                    if (_hc.second) {
                        for (auto& elem : rightTensors) {
                            elem.second = !elem.second;
                        }
                    }
                    SharedTensorData tA{&a, [](IContainer*) {}};
                    auto res = SharedTensorData{contractStar({tA, _hc.first}, rightTensors, contractions)};
                    if (res->rank() == 0) {
                        currentWeight *= res->element({});
                    } else {
                        currentTerm.emplace_back(res, false);
                    }
                } else if (shouldboomerang) {
                    // boomerang topology: sparse a contracted with two rank-2 tensors (one index each);
                    // result has same rank as a
                    auto b0_elem = elemB[get<1>(chunks[1])[0]];
                    auto b1_elem = elemB[get<1>(chunks[1])[1]];
                    if (_hc.second) {
                        b0_elem.second = !b0_elem.second;
                        b1_elem.second = !b1_elem.second;
                    }
                    pair<IndexPair, IndexPair> boomContractions = {
                        {_indices[0].first,
                         static_cast<IndexType>(b.getIndexing().getElementIndex(_indices[0].second).second)},
                        {_indices[1].first,
                         static_cast<IndexType>(b.getIndexing().getElementIndex(_indices[1].second).second)}};
                    SharedTensorData tA{&a, [](IContainer*) {}};
                    auto res =
                        SharedTensorData{contractR2Boomerang({tA, _hc.first}, {b0_elem, b1_elem}, boomContractions)};
                    currentTerm.emplace_back(res, false);
                } else {
                    OuterElemIterator::EntryType rightTensors;
                    rightTensors.reserve(get<1>(chunks[1]).size());
                    transform(get<1>(chunks[1]).begin(), get<1>(chunks[1]).end(), back_inserter(rightTensors),
                              [&](IndexType idx) -> const pair<SharedTensorData, bool>& { return elemB[idx]; });
                    IndexList inners(get<2>(chunks[1]).size());
                    FlagList innerAdds(inners.size(), false);
                    auto newdimlabs = getNewIndexLabels(a.getIndexing(), b.getIndexing(), chunks[1]);
                    size_t totalRankB =
                        accumulate(rightTensors.begin(), rightTensors.end(), 0ul,
                                   [](PositionType tot, const pair<SharedTensorData, bool>& elem) -> PositionType {
                                       return tot + elem.first->rank();
                                   });
                    if (totalRankB == inners.size()) {
                        IndexList::iterator itP1;
                        IndexList::iterator itP2;
                        if (newdimlabs.first.empty()) {
                            Base::ElementType newscal;
                            for (const auto& elemL : a) {
                                itP1 = inners.begin();
                                (void) a.getIndexing().splitPosition(elemL.first, get<0>(leftinfo), get<1>(leftinfo),
                                                                     inners, innerAdds);
                                Base::ElementType firstTerm = _hc.first ? conj(elemL.second) : elemL.second;
                                for (auto& entry : rightTensors) {
                                    itP2 = itP1 + static_cast<ptrdiff_t>(entry.first->rank());
                                    Base::ElementType secondTerm = (!_hc.second != !entry.second)
                                                                       ? conj(entry.first->element(itP1, itP2))
                                                                       : entry.first->element(itP1, itP2);
                                    firstTerm *= secondTerm;
                                    if (isZero(secondTerm)) {
                                        break;
                                    }
                                    itP1 = itP2;
                                }
                                newscal += firstTerm;
                            }
                            currentWeight *= newscal;
                        } else {
                            auto newsparse = makeEmptySparse(newdimlabs.first, newdimlabs.second);
                            auto* result = static_cast<STensor*>(newsparse.get());
                            for (const auto& elemL : a) {
                                itP1 = inners.begin();
                                auto tmpLeft = a.getIndexing().splitPosition(elemL.first, get<0>(leftinfo),
                                                                             get<1>(leftinfo), inners, innerAdds);
                                Base::ElementType firstTerm = _hc.first ? conj(elemL.second) : elemL.second;
                                for (auto& entry : rightTensors) {
                                    itP2 = itP1 + static_cast<ptrdiff_t>(entry.first->rank());
                                    Base::ElementType secondTerm = (!_hc.second != !entry.second)
                                                                       ? conj(entry.first->element(itP1, itP2))
                                                                       : entry.first->element(itP1, itP2);
                                    firstTerm *= secondTerm;
                                    if (isZero(secondTerm)) {
                                        break;
                                    }
                                    itP1 = itP2;
                                }
                                (*result)[tmpLeft] += firstTerm;
                            }
                            SharedTensorData tmpShared{newsparse.release()};
                            currentTerm.emplace_back(tmpShared, false);
                        }
                    } else {
                        OuterElemIterator itBEnd = OuterElemIterator{rightTensors}.end();
                        if (newdimlabs.first.empty()) {
                            Base::ElementType newscal;
                            for (const auto& elemL : a) {
                                (void) a.getIndexing().splitPosition(elemL.first, get<0>(leftinfo), get<1>(leftinfo),
                                                                     inners, innerAdds);
                                Base::ElementType firstTerm = _hc.first ? conj(elemL.second) : elemL.second;
                                OuterElemIterator itB{rightTensors};
                                for (; itB != itBEnd; ++itB) {
                                    PositionType tmpRight =
                                        b.getIndexing().splitPosition(itB, chunks[1], get<0>(rightinfo[0]),
                                                                      get<1>(rightinfo[0]), inners, innerAdds, true);
                                    if (tmpRight == numeric_limits<size_t>::max()) {
                                        continue;
                                    }
                                    Base::ElementType secondTerm = _hc.second ? conj(*itB) : *itB;
                                    newscal += firstTerm * secondTerm;
                                }
                            }
                            currentWeight *= newscal;
                        } else {
                            auto newsparse = makeEmptySparse(newdimlabs.first, newdimlabs.second);
                            auto* result = static_cast<STensor*>(newsparse.get());
                            for (const auto& elemL : a) {
                                auto tmpLeft = a.getIndexing().splitPosition(elemL.first, get<0>(leftinfo),
                                                                             get<1>(leftinfo), inners, innerAdds);
                                Base::ElementType firstTerm = _hc.first ? conj(elemL.second) : elemL.second;
                                OuterElemIterator itB{rightTensors};
                                for (; itB != itBEnd; ++itB) {
                                    PositionType tmpRight =
                                        b.getIndexing().splitPosition(itB, chunks[1], get<0>(rightinfo[0]),
                                                                      get<1>(rightinfo[0]), inners, innerAdds, true);
                                    if (tmpRight == numeric_limits<size_t>::max()) {
                                        continue;
                                    }
                                    Base::ElementType secondTerm = _hc.second ? conj(*itB) : *itB;
                                    (*result)[tmpLeft * get<2>(rightinfo[0]) + tmpRight] += firstTerm * secondTerm;
                                }
                            }
                            SharedTensorData tmpShared{newsparse.release()};
                            currentTerm.emplace_back(tmpShared, false);
                        }
                    }
                }
                for (auto elem : get<1>(chunks[0])) {
                    currentTerm.insert(currentTerm.end(), {elemB[elem].first, !elemB[elem].second != !_hc.second});
                }
                if (currentTerm.empty()) {
                    if (fullResult.get() == nullptr) {
                        fullResult = makeScalar(currentWeight); // real?
                    } else {
                        fullResult->element({}) += currentWeight; // real?
                    }
                } else if (currentTerm.size() == 1) {
                    TensorData out = currentTerm[0].first->clone();
                    if (currentTerm[0].second) {
                        out->conjugate();
                    }
                    if (!isZero(currentWeight - 1.)) {
                        out->operator*=(currentWeight);
                    }
                    if (fullResult.get() == nullptr) {
                        fullResult = std::move(out);
                    } else {
                        Ops::Sum summer{};
                        fullResult = calc2(std::move(fullResult), *out, summer, "sum_outerdot");
                    }
                } else {
                    if (!isZero(currentWeight - 1.)) {
                        auto temp = currentTerm.back();
                        currentTerm.pop_back();
                        auto candNew = temp.first->clone();
                        candNew->operator*=(temp.second ? conj(currentWeight) : currentWeight);
                        currentTerm.emplace_back(SharedTensorData{candNew.release()}, temp.second);
                    }
                    if (oResult != nullptr) {
                        oResult->addTerm(currentTerm);
                    } else {
                        fullResult = combineSharedTensors(std::move(currentTerm));
                        oResult = static_cast<OTensor*>(fullResult.get());
                    }
                }
            }
            return static_cast<Base*>(fullResult.release());
        }

        Base* Dot::operator()(VTensor& a, const OTensor& b) {
            // Inlined PartitionContractions. Same structure as for Sparse.Outer: contracted in chunks[1] and
            // pass-throughs in chunks[0]
            DotGroupList chunks(2);
            {
                IndexList rfree(b.getIndexing().numSubIndexing());
                iota(rfree.begin(), rfree.end(), 0);
                get<2>(chunks[1]) = _indices;
                for (const auto& elem : _indices) {
                    auto rloc = b.getIndexing().getElementIndex(elem.second).first;
                    rfree.erase(remove(rfree.begin(), rfree.end(), rloc), rfree.end());
                    get<1>(chunks[1]).push_back(rloc);
                }
                sort(get<2>(chunks[1]).begin(), get<2>(chunks[1]).end(),
                     [](const IndexPair& x, const IndexPair& y) { return x.second < y.second; });
                chunks[0] = DotGroupType{IndexList{}, rfree, {}};
            }

            TensorData fullResult;
            OTensor* oResult = nullptr;
            auto rightinfo = b.getIndexing().processShifts(chunks, IndexPairMember::Right);
            ASSERT(rightinfo.size() == 1);

            // Decide the topology
            bool should1ns = true;
            bool shouldboomerang = true;
            for (auto& elem : get<1>(chunks[1])) {
                should1ns &= b.begin()->at(elem).first->rank() == 1;
                shouldboomerang &= b.begin()->at(elem).first->rank() == 2;
            }
            shouldboomerang &= get<1>(chunks[1]).size() == 2 && get<2>(chunks[1]).size() == 2;

            // build contractions for star topology
            PositionList sortedOrder(_indices.size());
            IndexPairList contractions{};
            if (should1ns) {
                iota(sortedOrder.begin(), sortedOrder.end(), 0);
                sort(sortedOrder.begin(), sortedOrder.end(),
                     [&](size_t x, size_t y) { return _indices[x].first < _indices[y].first; });
                contractions.reserve(_indices.size());
                for (size_t k = 0; k < sortedOrder.size(); ++k) {
                    contractions.emplace_back(_indices[sortedOrder[k]].first, static_cast<IndexType>(k));
                }
            }

            // build strides for the general and boomerang cases
            IndexList innerDims(_indices.size());
            for (size_t k = 0; k < _indices.size(); ++k) {
                innerDims[k] = a.dims()[_indices[k].first];
            }
            SequentialIndexing innerIndexing{innerDims};

            IndexPairList contractionsForStrides;
            contractionsForStrides.reserve(_indices.size());
            for (size_t k = 0; k < _indices.size(); ++k) {
                contractionsForStrides.emplace_back(_indices[k].first, static_cast<IndexType>(k));
            }
            auto stridesA = a.getIndexing().getInnerOuterStrides(contractionsForStrides, innerIndexing.strides());

            IndexList inners(_indices.size());
            FlagList innerAdds(_indices.size(), false);

            for (const auto& elemB : b) {
                Base::ElementType currentWeight = 1.;
                OuterElemIterator::EntryType currentTerm;

                if (should1ns) {
                    OuterElemIterator::EntryType rightTensors;
                    rightTensors.reserve(sortedOrder.size());
                    for (auto i : sortedOrder) {
                        rightTensors.push_back(elemB[get<1>(chunks[1])[i]]);
                    }
                    if (_hc.second) {
                        for (auto& e : rightTensors) {
                            e.second = !e.second;
                        }
                    }
                    SharedTensorData tA{&a, [](IContainer*) {}};
                    auto res = SharedTensorData{contractStar({tA, _hc.first}, rightTensors, contractions)};
                    if (res->rank() == 0) {
                        currentWeight *= res->element({});
                    } else {
                        currentTerm.emplace_back(res, false);
                    }
                } else if (shouldboomerang) {
                    auto b0_elem = elemB[get<1>(chunks[1])[0]];
                    auto b1_elem = elemB[get<1>(chunks[1])[1]];
                    if (_hc.second) {
                        b0_elem.second = !b0_elem.second;
                        b1_elem.second = !b1_elem.second;
                    }
                    pair<IndexPair, IndexPair> boomContractions = {
                        {_indices[0].first,
                         static_cast<IndexType>(b.getIndexing().getElementIndex(_indices[0].second).second)},
                        {_indices[1].first,
                         static_cast<IndexType>(b.getIndexing().getElementIndex(_indices[1].second).second)}};
                    SharedTensorData tA{&a, [](IContainer*) {}};
                    auto res =
                        SharedTensorData{contractR2Boomerang({tA, _hc.first}, {b0_elem, b1_elem}, boomContractions)};
                    currentTerm.emplace_back(res, false);
                } else {
                    OuterElemIterator::EntryType rightTensors;
                    rightTensors.reserve(get<1>(chunks[1]).size());
                    transform(get<1>(chunks[1]).begin(), get<1>(chunks[1]).end(), back_inserter(rightTensors),
                              [&](IndexType idx) -> const pair<SharedTensorData, bool>& { return elemB[idx]; });
                    auto newdimlabs = getNewIndexLabels(a.getIndexing(), b.getIndexing(), chunks[1]);
                    size_t totalRankB =
                        accumulate(rightTensors.begin(), rightTensors.end(), 0ul,
                                   [](PositionType tot, const pair<SharedTensorData, bool>& elem) -> PositionType {
                                       return tot + elem.first->rank();
                                   });
                    if (totalRankB == inners.size()) {
                        IndexList::iterator itP1;
                        IndexList::iterator itP2;
                        if (newdimlabs.first.empty()) {
                            Base::ElementType newscal;
                            for (PositionType i = 0; i < a.numValues(); ++i) {
                                if (isZero(a[i])) {
                                    continue;
                                }
                                itP1 = inners.begin();
                                Base::ElementType firstTerm = _hc.first ? conj(a[i]) : a[i];
                                auto pospairs = LabeledIndexing<SequentialIndexing>::splitPosition(i, stridesA);
                                innerIndexing.posToIndices(pospairs.second, inners);
                                fill(innerAdds.begin(), innerAdds.end(), true);
                                for (auto& entry : rightTensors) {
                                    itP2 = itP1 + static_cast<ptrdiff_t>(entry.first->rank());
                                    Base::ElementType secondTerm = (!_hc.second != !entry.second)
                                                                       ? conj(entry.first->element(itP1, itP2))
                                                                       : entry.first->element(itP1, itP2);
                                    firstTerm *= secondTerm;
                                    if (isZero(secondTerm)) {
                                        break;
                                    }
                                    itP1 = itP2;
                                }
                                newscal += firstTerm;
                            }
                            currentWeight *= newscal;
                        } else {
                            auto newvect = makeEmptyVector(newdimlabs.first, newdimlabs.second);
                            auto* result = static_cast<VTensor*>(newvect.release());
                            for (PositionType i = 0; i < a.numValues(); ++i) {
                                if (isZero(a[i])) {
                                    continue;
                                }
                                itP1 = inners.begin();
                                Base::ElementType firstTerm = _hc.first ? conj(a[i]) : a[i];
                                auto pospairs = LabeledIndexing<SequentialIndexing>::splitPosition(i, stridesA);
                                innerIndexing.posToIndices(pospairs.second, inners);
                                fill(innerAdds.begin(), innerAdds.end(), true);
                                for (auto& entry : rightTensors) {
                                    itP2 = itP1 + static_cast<ptrdiff_t>(entry.first->rank());
                                    Base::ElementType secondTerm = (!_hc.second != !entry.second)
                                                                       ? conj(entry.first->element(itP1, itP2))
                                                                       : entry.first->element(itP1, itP2);
                                    firstTerm *= secondTerm;
                                    if (isZero(secondTerm)) {
                                        break;
                                    }
                                    itP1 = itP2;
                                }
                                (*result)[pospairs.first] += firstTerm;
                            }
                            SharedTensorData tmpShared{result};
                            currentTerm.emplace_back(tmpShared, false);
                        }
                    } else {
                        OuterElemIterator itBEnd = OuterElemIterator{rightTensors}.end();
                        if (newdimlabs.first.empty()) {
                            Base::ElementType newscal;
                            for (PositionType i = 0; i < a.numValues(); ++i) {
                                if (isZero(a[i])) {
                                    continue;
                                }
                                Base::ElementType firstTerm = _hc.first ? conj(a[i]) : a[i];
                                auto pospairs = LabeledIndexing<SequentialIndexing>::splitPosition(i, stridesA);
                                innerIndexing.posToIndices(pospairs.second, inners);
                                fill(innerAdds.begin(), innerAdds.end(), true);
                                OuterElemIterator itB{rightTensors};
                                for (; itB != itBEnd; ++itB) {
                                    PositionType tmpRight =
                                        b.getIndexing().splitPosition(itB, chunks[1], get<0>(rightinfo[0]),
                                                                      get<1>(rightinfo[0]), inners, innerAdds, true);
                                    if (tmpRight == numeric_limits<size_t>::max()) {
                                        continue;
                                    }
                                    Base::ElementType secondTerm = _hc.second ? conj(*itB) : *itB;
                                    newscal += firstTerm * secondTerm;
                                }
                            }
                            currentWeight *= newscal;
                        } else {
                            auto newsparse = makeEmptySparse(newdimlabs.first, newdimlabs.second);
                            auto* result = static_cast<STensor*>(newsparse.get());
                            for (PositionType i = 0; i < a.numValues(); ++i) {
                                if (isZero(a[i])) {
                                    continue;
                                }
                                Base::ElementType firstTerm = _hc.first ? conj(a[i]) : a[i];
                                auto pospairs = LabeledIndexing<SequentialIndexing>::splitPosition(i, stridesA);
                                innerIndexing.posToIndices(pospairs.second, inners);
                                fill(innerAdds.begin(), innerAdds.end(), true);
                                OuterElemIterator itB{rightTensors};
                                for (; itB != itBEnd; ++itB) {
                                    PositionType tmpRight =
                                        b.getIndexing().splitPosition(itB, chunks[1], get<0>(rightinfo[0]),
                                                                      get<1>(rightinfo[0]), inners, innerAdds, true);
                                    if (tmpRight == numeric_limits<size_t>::max()) {
                                        continue;
                                    }
                                    Base::ElementType secondTerm = _hc.second ? conj(*itB) : *itB;
                                    (*result)[pospairs.first * get<2>(rightinfo[0]) + tmpRight] +=
                                        firstTerm * secondTerm;
                                }
                            }
                            SharedTensorData tmpShared{newsparse.release()};
                            currentTerm.emplace_back(tmpShared, false);
                        }
                    }
                }

                for (auto elem : get<1>(chunks[0])) {
                    currentTerm.insert(currentTerm.end(), {elemB[elem].first, !elemB[elem].second != !_hc.second});
                }

                if (currentTerm.empty()) {
                    if (fullResult.get() == nullptr) {
                        fullResult = makeScalar(currentWeight); // real?
                    } else {
                        fullResult->element({}) += currentWeight; // real?
                    }
                } else if (currentTerm.size() == 1) {
                    TensorData out = currentTerm[0].first->clone();
                    if (currentTerm[0].second) {
                        out->conjugate();
                    }
                    if (!isZero(currentWeight - 1.)) {
                        out->operator*=(currentWeight);
                    }
                    if (fullResult.get() == nullptr) {
                        fullResult = std::move(out);
                    } else {
                        Ops::Sum summer{};
                        fullResult = calc2(std::move(fullResult), *out, summer, "sum_outerdot");
                    }
                } else {
                    if (!isZero(currentWeight - 1.)) {
                        auto temp = currentTerm.back();
                        currentTerm.pop_back();
                        auto candNew = temp.first->clone();
                        candNew->operator*=(temp.second ? conj(currentWeight) : currentWeight);
                        currentTerm.emplace_back(SharedTensorData{candNew.release()}, temp.second);
                    }
                    if (oResult != nullptr) {
                        oResult->addTerm(currentTerm);
                    } else {
                        fullResult = combineSharedTensors(std::move(currentTerm));
                        oResult = static_cast<OTensor*>(fullResult.get());
                    }
                }
            }
            return static_cast<Base*>(fullResult.release());
        }

        Base* Dot::operator()(Base& a, const Base& b) {
            auto newdimlabs = getNewIndexLabels(a, b);
            if (newdimlabs.first.empty()) {
                auto newscalar = makeEmptyScalar();
                IContainer::ElementType res = 0.;
                const BruteForceIteratorRange bf{a.dims()};
                for (const auto& elem : bf) {
                    auto aVal = a.element(elem);
                    if (isZero(aVal)) {
                        continue;
                    }
                    if (_hc.first) {
                        aVal = conj(aVal);
                    }
                    IndexList fixed = b.dims();
                    for (auto idx : _indices) {
                        fixed[idx.second] = elem[idx.first];
                    }
                    auto bVal = b.element(fixed);
                    if (isZero(bVal)) {
                        continue;
                    }
                    if (_hc.second) {
                        bVal = conj(bVal);
                    }
                    res += aVal * bVal;
                }
                newscalar->element({}) = res;
                return newscalar.release();
            }
            auto tmp = makeEmptySparse(newdimlabs.first, newdimlabs.second);
            Base* results = tmp.release();
            const BruteForceIteratorRange bf{a.dims()};
            for (const auto& elem : bf) {
                auto aVal = a.element(elem);
                if (isZero(aVal)) {
                    continue;
                }
                if (_hc.first) {
                    aVal = conj(aVal);
                }
                IndexList fixed = b.dims();
                for (auto idx : _indices) {
                    fixed[idx.second] = elem[idx.first];
                }
                const BruteForceIteratorRange bf2{b.dims(), fixed};
                for (const auto& elem2 : bf2) {
                    auto bVal = b.element(elem2);
                    if (isZero(bVal)) {
                        continue;
                    }
                    if (_hc.second) {
                        bVal = conj(bVal);
                    }
                    auto idxRes = combineIndex(elem, elem2);
                    results->element(idxRes) += aVal * bVal;
                }
            }
            return results;
        }

        Base* Dot::error(Base& /*unused*/, const Base& /*unused*/) {
            throw Error("Invalid data types for tensor Dot");
        }

        IndexList Dot::combineIndex(const IndexList& a, const IndexList& b) const {
            IndexList result = a;
            for (auto elem : reverse_range(_idxLeft)) {
                result.erase(result.begin() + elem);
            }
            size_t base = result.size();
            result.insert(result.end(), b.begin(), b.end());
            for (auto elem : reverse_range(_idxRight)) {
                result.erase(result.begin() + static_cast<ptrdiff_t>(base + elem));
            }
            return result;
        }

        pair<IndexList, LabelsList> Dot::getNewIndexLabels(const Base& first, const Base& second) const {
            IndexList resultD = first.dims();
            LabelsList resultL = _hc.first ? flipListOfLabels(first.labels()) : first.labels();
            for (auto elem : reverse_range(_idxLeft)) {
                resultD.erase(resultD.begin() + elem);
                resultL.erase(resultL.begin() + elem);
            }
            size_t base = resultD.size();
            const IndexList& dim2 = second.dims();
            const LabelsList& lab2 = _hc.second ? flipListOfLabels(second.labels()) : second.labels();
            resultD.insert(resultD.end(), dim2.begin(), dim2.end());
            resultL.insert(resultL.end(), lab2.begin(), lab2.end());
            for (auto elem : reverse_range(_idxRight)) {
                resultD.erase(resultD.begin() + static_cast<ptrdiff_t>(base + elem));
                resultL.erase(resultL.begin() + static_cast<ptrdiff_t>(base + elem));
            }
            return make_pair(resultD, resultL);
        }

        SharedTensorData Dot::calcSharedDot(SharedTensorData origin, const IContainer& other,
                                            const IndexPairList& indices, pair<bool, bool> shouldHC) {
            Ops::Dot dotter{indices, shouldHC};
            return calc2(std::move(origin), other, dotter, "dot");
        }

        template <size_t N, typename U, typename... Types>
        enable_if_t<is_convertible_v<multidim_vector<U>, tuple_element_t<N, tuple<Types...>>>, bool>
        matchPartitions(const tuple<Types...>& data, U value) {
            return find(get<N>(data).begin(), get<N>(data).end(), value) != get<N>(data).end();
        }

        template <size_t N, typename U, typename... Types>
        enable_if_t<is_convertible_v<multidim_vector<U>, tuple_element_t<N, tuple<Types...>>>, void>
        addPartitionEntry(tuple<Types...>& data, U value) {
            if (find(get<N>(data).begin(), get<N>(data).end(), value) == get<N>(data).end()) {
                get<N>(data).push_back(value);
            }
        }

        template <size_t N, typename... Types>
        enable_if_t<(N < sizeof...(Types)), void> appendPartitionEntries(const tuple<Types...>& from,
                                                                         tuple<Types...>& to) {
            get<N>(to).insert(get<N>(to).end(), get<N>(from).begin(), get<N>(from).end());
        }


        DotGroupList Dot::partitionContractions(const LabeledIndexing<AlignedIndexing>& /*unused*/,
                                                const BlockIndexing& rhs) const {
            DotGroupList partitions(2);
            IndexList rfree(rhs.numSubIndexing());
            iota(rfree.begin(), rfree.end(), 0);
            get<2>(partitions[1]) = _indices;
            for (const auto& elem : _indices) {
                auto rloc = rhs.getElementIndex(elem.second).first;
                rfree.erase(remove(rfree.begin(), rfree.end(), rloc), rfree.end());
                get<1>(partitions[1]).push_back(rloc);
            }
            // sort contractions
            for (auto& elem : partitions) {
                std::sort(get<2>(elem).begin(), get<2>(elem).end(),
                          [](const IndexPair& a, const IndexPair& b) -> bool { return a.second < b.second; });
            }
            // Add in untouched tensors
            partitions[0] = DotGroupType{IndexList{}, rfree, {}};
            return partitions;
        }

        DotGroupList Dot::partitionContractions(const BlockIndexing& lhs,
                                                const LabeledIndexing<AlignedIndexing>& /*unused*/) const {
            DotGroupList partitions(2);
            IndexList lfree(lhs.numSubIndexing());
            iota(lfree.begin(), lfree.end(), 0);
            get<2>(partitions[1]) = _indices;
            for (const auto& elem : _indices) {
                auto lloc = lhs.getElementIndex(elem.first).first;
                lfree.erase(remove(lfree.begin(), lfree.end(), lloc), lfree.end());
                get<0>(partitions[1]).push_back(lloc);
            }
            // sort contractions
            for (auto& elem : partitions) {
                std::sort(get<2>(elem).begin(), get<2>(elem).end(),
                          [](const IndexPair& a, const IndexPair& b) -> bool { return a.second < b.second; });
            }
            // Add in untouched tensors
            partitions[0] = DotGroupType{lfree, IndexList{}, {}};
            return partitions;
        }

        pair<IndexList, LabelsList> Dot::getNewIndexLabels(const LabeledIndexing<AlignedIndexing>& lhs,
                                                           const BlockIndexing& rhs, const DotGroupType& chunk) {
            IndexList dims;
            LabelsList labels;
            map<IndexType, IndexType> rightPosMaps;
            IndexType offset = 0;
            dims.insert(dims.end(), lhs.dims().begin(), lhs.dims().end());
            labels.insert(labels.end(), lhs.labels().begin(), lhs.labels().end());
            offset = static_cast<IndexType>(offset + lhs.rank());
            for (auto elem : get<1>(chunk)) {
                rightPosMaps.insert({elem, offset});
                dims.insert(dims.end(), rhs.getSubIndexing(elem).dims().begin(), rhs.getSubIndexing(elem).dims().end());
                labels.insert(labels.end(), rhs.getSubIndexing(elem).labels().begin(),
                              rhs.getSubIndexing(elem).labels().end());
                offset = static_cast<IndexType>(offset + rhs.getSubIndexing(elem).rank());
            }
            set<size_t> deletes;
            for (const auto& elem : get<2>(chunk)) {
                auto right = rhs.getElementIndex(elem.second);
                deletes.insert(elem.first);
                deletes.insert(rightPosMaps[right.first] + right.second);
            }
            for (auto elem : reverse_range(deletes)) {
                dims.erase(dims.begin() + static_cast<ptrdiff_t>(elem));
                labels.erase(labels.begin() + static_cast<ptrdiff_t>(elem));
            }
            return {dims, labels};
        }

        pair<IndexList, LabelsList> Dot::getNewIndexLabels(const LabeledIndexing<SequentialIndexing>& lhs,
                                                           const BlockIndexing& rhs, const DotGroupType& chunk) {
            IndexList dims;
            LabelsList labels;
            map<IndexType, IndexType> rightPosMaps;
            IndexType offset = 0;
            dims.insert(dims.end(), lhs.dims().begin(), lhs.dims().end());
            labels.insert(labels.end(), lhs.labels().begin(), lhs.labels().end());
            offset = static_cast<IndexType>(offset + lhs.rank());
            for (auto elem : get<1>(chunk)) {
                rightPosMaps.insert({elem, offset});
                dims.insert(dims.end(), rhs.getSubIndexing(elem).dims().begin(), rhs.getSubIndexing(elem).dims().end());
                labels.insert(labels.end(), rhs.getSubIndexing(elem).labels().begin(),
                              rhs.getSubIndexing(elem).labels().end());
                offset = static_cast<IndexType>(offset + rhs.getSubIndexing(elem).rank());
            }
            set<size_t> deletes;
            for (const auto& elem : get<2>(chunk)) {
                auto right = rhs.getElementIndex(elem.second);
                deletes.insert(elem.first);
                deletes.insert(rightPosMaps[right.first] + right.second);
            }
            for (auto elem : reverse_range(deletes)) {
                dims.erase(dims.begin() + static_cast<ptrdiff_t>(elem));
                labels.erase(labels.begin() + static_cast<ptrdiff_t>(elem));
            }
            return {dims, labels};
        }

        DotGroupList Dot::partitionContractions(const BlockIndexing& lhs,
                                                                   const BlockIndexing& rhs) const {
                DotGroupList partitions;
                FlagList validPartitions;
                IndexList lfree(lhs.numSubIndexing());
                IndexList rfree(rhs.numSubIndexing());
                iota(lfree.begin(), lfree.end(), 0);
                iota(rfree.begin(), rfree.end(), 0);
                //Placeholder for the untouched
                auto frontelem = make_tuple<IndexList, IndexList, IndexPairList>({}, {}, {});
                partitions.push_back(frontelem);
                validPartitions.push_back(true);
                for (auto& elem : _indices) {
                    auto lloc = lhs.getElementIndex(elem.first).first;
                    auto rloc = rhs.getElementIndex(elem.second).first;
                    lfree.erase(remove(lfree.begin(), lfree.end(), lloc), lfree.end());
                    rfree.erase(remove(rfree.begin(), rfree.end(), rloc), rfree.end());
                    PositionList finds{};
                    auto match = [&](size_t pos, IndexType valL, IndexType valR) -> bool {
                        const auto& data = partitions[pos];
                        return matchPartitions<0>(data, valL) || matchPartitions<1>(data, valR);
                    };
                    for (size_t i = 0; i< partitions.size(); ++i) { // Find all matches on the left or the right
                        if (validPartitions[i] && match(i, lloc, rloc)) {
                            finds.push_back(i);
                        }
                    }
                    auto merge = [&](IndexType lidx, IndexType ridx, IndexPair contraction) -> void {
                        auto& data = partitions[finds[0]];
                        addPartitionEntry<0>(data, lidx);
                        addPartitionEntry<1>(data, ridx);
                        addPartitionEntry<2>(data, contraction);
                    };
                    auto merge_range = [&](size_t other) -> void {
                        auto& from = partitions[finds[other]];
                        auto& to = partitions[finds[0]];
                        appendPartitionEntries<0>(from, to);
                        appendPartitionEntries<1>(from, to);
                        appendPartitionEntries<2>(from, to);
                    };
                    if (finds.size() > 0) { // Insert elem into first found match
                        merge(lloc, rloc, elem);
                        for (size_t idx = 1; idx < finds.size(); ++idx) { // merge in other finds into first found match
                            if(validPartitions[finds[idx]]) {
                                merge_range(idx);
                                validPartitions[finds[idx]] = false;
                            }
                        }
                    }
                    else {
                        auto newelem = make_tuple<IndexList, IndexList, IndexPairList>({lloc}, {rloc}, {elem});
                        partitions.push_back(newelem);
                        validPartitions.push_back(true);
                    }
                }
                // eliminate already merged partitions
                for(size_t i = partitions.size(); 0 < i--;) {
                    if(!validPartitions[i]) partitions.erase(partitions.begin() + static_cast<ptrdiff_t>(i));
                }
                // sort contractions
                for(auto& elem: partitions) {
                    std::sort(get<2>(elem).begin(), get<2>(elem).end(),
                              [](const IndexPair& a, const IndexPair& b) -> bool { return a.second < b.second; });
                }
                //Add in untouched tensors
                partitions[0] = DotGroupType{lfree, rfree, {}};
                return partitions;
            }


            pair<IndexList, LabelsList> Dot::getNewIndexLabels(const BlockIndexing& lhs, const BlockIndexing& rhs,
                                                               const DotGroupType& chunk) {
                IndexList dims;
                LabelsList labels;
                map<IndexType, IndexType> leftPosMaps;
                map<IndexType, IndexType> rightPosMaps;
                IndexType offset = 0;
                for (auto elem : get<0>(chunk)) {
                    leftPosMaps.insert({elem, offset});
                    dims.insert(dims.end(), lhs.getSubIndexing(elem).dims().begin(), lhs.getSubIndexing(elem).dims().end());
                    labels.insert(labels.end(), lhs.getSubIndexing(elem).labels().begin(), lhs.getSubIndexing(elem).labels().end());
                    offset = static_cast<IndexType>(offset + lhs.getSubIndexing(elem).rank());
                }
                for (auto elem : get<1>(chunk)) {
                    rightPosMaps.insert({elem, offset});
                    dims.insert(dims.end(), rhs.getSubIndexing(elem).dims().begin(), rhs.getSubIndexing(elem).dims().end());
                    labels.insert(labels.end(), rhs.getSubIndexing(elem).labels().begin(),
                                  rhs.getSubIndexing(elem).labels().end());
                    offset = static_cast<IndexType>(offset + rhs.getSubIndexing(elem).rank());
                }
                set<size_t> deletes;
                for(auto& elem: get<2>(chunk)) {
                    auto left = lhs.getElementIndex(elem.first);
                    auto right = rhs.getElementIndex(elem.second);
                    deletes.insert(leftPosMaps[left.first] + left.second);
                    deletes.insert(rightPosMaps[right.first] + right.second);
                }
                for(auto elem: reverse_range(deletes)) {
                    dims.erase(dims.begin() + static_cast<ptrdiff_t>(elem));
                    labels.erase(labels.begin() + static_cast<ptrdiff_t>(elem));
                }
                return {dims, labels};
            }

        pair<IndexList, LabelsList> Dot::getNewIndexLabels(const BlockIndexing& lhs,
                                                           const LabeledIndexing<AlignedIndexing>& rhs,
                                                           const DotGroupType& chunk) {
            IndexList dims;
            LabelsList labels;
            map<IndexType, IndexType> leftPosMaps;
            IndexType offset = 0;
            for (auto elem : get<0>(chunk)) {
                leftPosMaps.insert({elem, offset});
                dims.insert(dims.end(), lhs.getSubIndexing(elem).dims().begin(), lhs.getSubIndexing(elem).dims().end());
                labels.insert(labels.end(), lhs.getSubIndexing(elem).labels().begin(),
                              lhs.getSubIndexing(elem).labels().end());
                offset = static_cast<IndexType>(offset + lhs.getSubIndexing(elem).rank());
            }
            dims.insert(dims.end(), rhs.dims().begin(), rhs.dims().end());
            labels.insert(labels.end(), rhs.labels().begin(), rhs.labels().end());
            set<size_t> deletes;
            for (const auto& elem : get<2>(chunk)) {
                auto left = lhs.getElementIndex(elem.first);
                deletes.insert(leftPosMaps[left.first] + left.second);
                deletes.insert(offset + elem.second);
            }
            for (auto elem : reverse_range(deletes)) {
                dims.erase(dims.begin() + static_cast<ptrdiff_t>(elem));
                labels.erase(labels.begin() + static_cast<ptrdiff_t>(elem));
            }
            return {dims, labels};
        }


        pair<IndexList, LabelsList> Dot::getNewIndexLabels(const pair<SharedTensorData, bool>& a,
                                                           const IndexPairList& contracted) {
            const IndexList& srcD = a.first->dims();
            const LabelsList srcL = a.second ? flipListOfLabels(a.first->labels()) : a.first->labels();
            const size_t n = srcD.size();
            FlagList skip(n, false);
            for (const auto& elem : contracted) {
                skip[elem.first] = true;
            }
            IndexList resultD;
            LabelsList resultL;
            resultD.reserve(n - contracted.size());
            resultL.reserve(n - contracted.size());
            for (size_t i = 0; i < n; ++i) {
                if (!skip[i]) {
                    resultD.push_back(srcD[i]);
                    resultL.push_back(srcL[i]);
                }
            }
            return make_pair(resultD, resultL);
        }


        Base* Dot::contractSingles(const pair<SharedTensorData, bool>& tensorA,
                                   const pair<SharedTensorData, bool>& tensorB) {
            auto newscal = makeEmptyScalar();
            const auto* as = dynamic_cast<const SparseContainer*>(tensorA.first.get());
            if (as != nullptr) {
                const auto* bs = dynamic_cast<const SparseContainer*>(tensorB.first.get());
                if (bs != nullptr) {
                    // sparse sparse
                    auto itA = as->begin();
                    auto itB = bs->begin();
                    Base::ElementType newTerm = 0.;
                    while (itA != as->end() && itB != bs->end()) {
                        if (itA->first > itB->first) {
                            ++itB;
                        } else if (itB->first > itA->first) {
                            ++itA;
                        } else {
                            newTerm += (tensorA.second ? conj(itA->second) : itA->second) *
                                       (tensorB.second ? conj(itB->second) : itB->second);
                            ++itA;
                            ++itB;
                        }
                    }
                    newscal->element({}) = newTerm;
                } else {
                    const auto* bv = dynamic_cast<const VectorContainer*>(tensorB.first.get());
                    if (bv != nullptr) {
                        // sparse vector
                        for (const auto& elemA : *as) {
                            Base::ElementType newTerm = (*bv)[elemA.first];
                            if (!isZero(newTerm)) {
                                newscal->element({}) += (tensorB.second ? conj(newTerm) : newTerm) *
                                                        (tensorA.second ? conj(elemA.second) : elemA.second);
                            }
                        }
                    } else {
                        throw Error("Dot:invalid tensor type");
                    }
                }
            } else {
                const auto* av = dynamic_cast<const VectorContainer*>(tensorA.first.get());
                if (av == nullptr) {
                    throw Error("Dot:invalid tensor type");
                }
                const auto* bs = dynamic_cast<const SparseContainer*>(tensorB.first.get());
                if (bs != nullptr) {
                    // vector sparse
                    for (const auto& elemB : *bs) {
                        Base::ElementType newTerm = (*av)[elemB.first];
                        if (!isZero(newTerm)) {
                            newscal->element({}) += (tensorA.second ? conj(newTerm) : newTerm) *
                                                    (tensorB.second ? conj(elemB.second) : elemB.second);
                        }
                    }
                } else {
                    const auto* bv = dynamic_cast<const VectorContainer*>(tensorB.first.get());
                    if (bv != nullptr) {
                        // vector vector
                        Base::ElementType newTerm = 0.;
                        for (PositionType i = 0; i < av->dims()[0]; ++i) {
                            if (!isZero((*av)[i]) && !isZero((*bv)[i])) {
                                newTerm += (tensorA.second ? conj((*av)[i]) : (*av)[i]) *
                                           (tensorB.second ? conj((*bv)[i]) : (*bv)[i]);
                            }
                        }
                        newscal->element({}) = newTerm;
                    } else {
                        throw Error("Dot:invalid tensor type");
                    }
                }
            }
            return newscal.release();
        }


        Base* Dot::contractStar(const pair<SharedTensorData, bool>& tensorA,
                                const vector<pair<SharedTensorData, bool>>& tensorsB,
                                const IndexPairList& contractions) {
            auto newdimlabs = getNewIndexLabels(tensorA, contractions);
            const auto* as = dynamic_cast<const SparseContainer*>(tensorA.first.get());
            bool singleEdge = contractions.size() == 1;
            if (as != nullptr) {
                IndexList inners(contractions.size());
                FlagList innerAdds(inners.size(), false);
                auto leftinfo = as->getIndexing().processShifts(contractions, IndexPairMember::Left);
                if (newdimlabs.first.empty()) {
                    auto newscal = makeEmptyScalar();
                    for (const auto& elemL : *as) {
                        (void) as->getIndexing().splitPosition(elemL.first, get<0>(leftinfo), get<1>(leftinfo), inners,
                                                               innerAdds);
                        Base::ElementType newTerm = tensorA.second ? conj(elemL.second) : elemL.second;
                        for (auto idxMap : contractions) {
                            size_t pos = idxMap.second;
                            Base::ElementType secondTermPre =
                                tensorsB[pos].first->element(inners.begin() + static_cast<ptrdiff_t>(pos),
                                                             inners.begin() + static_cast<ptrdiff_t>(pos) + 1);
                            newTerm *= tensorsB[pos].second ? conj(secondTermPre) : secondTermPre;
                        }
                        newscal->element({}) += newTerm;
                    }
                    return newscal.release();
                }
                auto newsparse = makeEmptySparse(newdimlabs.first, newdimlabs.second);
                auto* result = static_cast<STensor*>(newsparse.release());
                for (const auto& elemL : *as) {
                    auto tmpLeft = as->getIndexing().splitPosition(elemL.first, get<0>(leftinfo), get<1>(leftinfo),
                                                                   inners, innerAdds);
                    Base::ElementType newTerm = tensorA.second ? conj(elemL.second) : elemL.second;
                    for (auto idxMap : contractions) {
                        size_t pos = idxMap.second;
                        Base::ElementType secondTermPre =
                            tensorsB[pos].first->element(inners.begin() + static_cast<ptrdiff_t>(pos),
                                                         inners.begin() + static_cast<ptrdiff_t>(pos) + 1);
                        newTerm *= tensorsB[pos].second ? conj(secondTermPre) : secondTermPre;
                    }
                    (*result)[tmpLeft] += newTerm;
                }
                return static_cast<Base*>(result);
            }
            const auto* av = dynamic_cast<const VectorContainer*>(tensorA.first.get());
            if (av != nullptr) {
                IndexList bIndices;
                bIndices.reserve(tensorsB.size());
                transform(tensorsB.begin(), tensorsB.end(), back_inserter(bIndices),
                          [](const pair<SharedTensorData, bool>& elem) -> IndexType { return elem.first->dims()[0]; });
                SequentialIndexing bIndexing{bIndices};
                auto stridesA = av->getIndexing().getInnerOuterStrides(contractions, bIndexing.strides());
                if (newdimlabs.first.empty()) {
                    auto newscal = makeEmptyScalar();

                    for (PositionType i = 0; i < av->numValues(); ++i) {
                        if (isZero((*av)[i])) {
                            continue;
                        }
                        Base::ElementType newTerm = tensorA.second ? conj((*av)[i]) : (*av)[i];
                        auto pospairs = LabeledIndexing<SequentialIndexing>::splitPosition(i, stridesA);
                        if (singleEdge) {
                            newTerm *=
                                tensorsB[0].second
                                    ? conj(tensorsB[0].first->element({static_cast<unsigned short>(pospairs.second)}))
                                    : tensorsB[0].first->element({static_cast<unsigned short>(pospairs.second)});
                        } else {
                            IndexList idxB(contractions.size());
                            bIndexing.posToIndices(pospairs.second, idxB);
                            for (IndexType j = 0; j < static_cast<IndexType>(contractions.size()); ++j) {
                                newTerm *= tensorsB[j].second ? conj(tensorsB[j].first->element({idxB[j]}))
                                                              : tensorsB[j].first->element({idxB[j]});
                            }
                        }
                        newscal->element({}) += newTerm;
                    }
                    return newscal.release();
                }
                auto newvect = makeEmptyVector(newdimlabs.first, newdimlabs.second);
                auto* result = static_cast<VTensor*>(newvect.release());
                for (PositionType i = 0; i < av->numValues(); ++i) {
                    if (isZero((*av)[i])) {
                        continue;
                    }
                    Base::ElementType newTerm = tensorA.second ? conj((*av)[i]) : (*av)[i];
                    auto pospairs = LabeledIndexing<SequentialIndexing>::splitPosition(i, stridesA);
                    if (singleEdge) {
                        newTerm *=
                            tensorsB[0].second
                                ? conj(tensorsB[0].first->element({static_cast<unsigned short>(pospairs.second)}))
                                : tensorsB[0].first->element({static_cast<unsigned short>(pospairs.second)});
                    } else {
                        IndexList idxB(contractions.size());
                        bIndexing.posToIndices(pospairs.second, idxB);
                        for (IndexType j = 0; j < static_cast<IndexType>(contractions.size()); ++j) {
                            newTerm *= tensorsB[j].second ? conj(tensorsB[j].first->element({idxB[j]}))
                                                          : tensorsB[j].first->element({idxB[j]});
                        }
                    }
                    (*result)[pospairs.first] += newTerm;
                }
                return static_cast<Base*>(result);
            }
            throw Error("Dot:invalid tensor type");
        }


        Base* Dot::contractR2Boomerang(const pair<SharedTensorData, bool>& tensorA,
                                       const pair<pair<SharedTensorData, bool>, pair<SharedTensorData, bool>>& tensorsB,
                                       const pair<IndexPair, IndexPair>& contractions) {
            const auto* as = dynamic_cast<const SparseContainer*>(tensorA.first.get());
            const auto* av = (as == nullptr) ? dynamic_cast<const VectorContainer*>(tensorA.first.get()) : nullptr;
            if (as == nullptr && av == nullptr) {
                throw Error("Dot::contractR2Boomerang: tensorA must be sparse or vector. Boom shakalaka!");
            }

            // contracted index in a and the contracted/free indices in each rank-2 B tensor
            IndexType c0 = contractions.first.first;
            IndexType b0c = contractions.first.second;
            auto b0f = static_cast<IndexType>(1u - b0c);
            IndexType c1 = contractions.second.first;
            IndexType b1c = contractions.second.second;
            auto b1f = static_cast<IndexType>(1u - b1c);

            const auto* b0s = dynamic_cast<const SparseContainer*>(tensorsB.first.first.get());
            const auto* b0v =
                (b0s == nullptr) ? dynamic_cast<const VectorContainer*>(tensorsB.first.first.get()) : nullptr;
            const auto* b1s = dynamic_cast<const SparseContainer*>(tensorsB.second.first.get());
            const auto* b1v =
                (b1s == nullptr) ? dynamic_cast<const VectorContainer*>(tensorsB.second.first.get()) : nullptr;
            if ((b0s == nullptr && b0v == nullptr) || (b1s == nullptr && b1v == nullptr)) {
                throw Error("Dot::contractR2Boomerang: tensorsB must be sparse or vector. 6-7!!");
            }

            // result has same rank as a; dims/labels at c0 and c1 are replaced by the free
            // index of B0 and B1 respectively
            IndexList resultDims = tensorA.first->dims();
            LabelsList resultLabels =
                tensorA.second ? flipListOfLabels(tensorA.first->labels()) : tensorA.first->labels();
            LabelsList b0labels = tensorsB.first.second ? flipListOfLabels(tensorsB.first.first->labels())
                                                        : tensorsB.first.first->labels();
            LabelsList b1labels = tensorsB.second.second ? flipListOfLabels(tensorsB.second.first->labels())
                                                         : tensorsB.second.first->labels();
            resultDims[c0] = tensorsB.first.first->dims()[b0f];
            resultDims[c1] = tensorsB.second.first->dims()[b1f];
            resultLabels[c0] = b0labels[b0f];
            resultLabels[c1] = b1labels[b1f];

            auto newsparse = makeEmptySparse(resultDims, resultLabels);
            auto* result = static_cast<STensor*>(newsparse.release());

            IndexType freeD0 = resultDims[c0];
            IndexType freeD1 = resultDims[c1];
            IndexList aIdxList;
            aIdxList.resize(tensorA.first->rank());
            IndexList b0idx(2);
            IndexList b1idx(2);

            auto accum = [&](IndexType f0, Base::ElementType aVal, Base::ElementType val0, IndexType f1,
                             Base::ElementType val1) {
                aIdxList[c0] = f0;
                aIdxList[c1] = f1;
                (*result)[result->getIndexing().indicesToPos(aIdxList.begin(), aIdxList.end())] += aVal * val0 * val1;
            };

            auto doLoop = [&](IndexType v0, IndexType v1, Base::ElementType aVal) {
                auto visitB1 = [&](IndexType f0, Base::ElementType val0) {
                    if (b1s != nullptr) {
                        for (const auto& elemB1 : *b1s) {
                            if (b1s->getIndexing().ithIndexInPos(elemB1.first, b1c) != v1) {
                                continue;
                            }
                            IndexType f1 = b1s->getIndexing().ithIndexInPos(elemB1.first, b1f);
                            accum(f0, aVal, val0, f1, tensorsB.second.second ? conj(elemB1.second) : elemB1.second);
                        }
                    } else {
                        b1idx[b1c] = v1;
                        for (IndexType f1 = 0; f1 < freeD1; ++f1) {
                            b1idx[b1f] = f1;
                            auto raw = (*b1v)[b1v->getIndexing().indicesToPos(b1idx.begin(), b1idx.end())];
                            if (isZero(raw)) {
                                continue;
                            }
                            accum(f0, aVal, val0, f1, tensorsB.second.second ? conj(raw) : raw);
                        }
                    }
                };
                if (b0s != nullptr) {
                    for (const auto& elemB0 : *b0s) {
                        if (b0s->getIndexing().ithIndexInPos(elemB0.first, b0c) != v0) {
                            continue;
                        }
                        IndexType f0 = b0s->getIndexing().ithIndexInPos(elemB0.first, b0f);
                        visitB1(f0, tensorsB.first.second ? conj(elemB0.second) : elemB0.second);
                    }
                } else {
                    b0idx[b0c] = v0;
                    for (IndexType f0 = 0; f0 < freeD0; ++f0) {
                        b0idx[b0f] = f0;
                        auto raw = (*b0v)[b0v->getIndexing().indicesToPos(b0idx.begin(), b0idx.end())];
                        if (isZero(raw)) {
                            continue;
                        }
                        visitB1(f0, tensorsB.first.second ? conj(raw) : raw);
                    }
                }
            };

            if (as != nullptr) {
                for (const auto& elemA : *as) {
                    as->getIndexing().posToIndices(elemA.first, aIdxList);
                    doLoop(aIdxList[c0], aIdxList[c1], tensorA.second ? conj(elemA.second) : elemA.second);
                }
            } else {
                for (PositionType pos = 0; pos < av->numValues(); ++pos) {
                    auto raw = (*av)[pos];
                    if (isZero(raw)) {
                        continue;
                    }
                    av->getIndexing().posToIndices(pos, aIdxList);
                    doLoop(aIdxList[c0], aIdxList[c1], tensorA.second ? conj(raw) : raw);
                }
            }
            return static_cast<Base*>(result);
        }

    } // namespace Ops


} // namespace Hammer::MultiDimensional
