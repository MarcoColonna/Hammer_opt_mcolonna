///
/// @file  Trace.cc
/// @brief Tensor trace algorithm
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/Math/MultiDim/Ops/Trace.hh"
#include "Hammer/Math/MultiDim/Ops/Dot.hh"
#include "Hammer/Math/MultiDim/Ops/Sum.hh"
#include "Hammer/Math/MultiDim/ScalarContainer.hh"
#include "Hammer/Math/MultiDim/VectorContainer.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"
#include "Hammer/Math/MultiDim/OuterContainer.hh"
#include "Hammer/Exceptions.hh"
#include "Hammer/Tools/Utils.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Math/MultiDim/OperationDefs.hh"

using namespace std;

namespace Hammer::MultiDimensional {

    using VTensor = VectorContainer;
    using STensor = SparseContainer;
    using OTensor = OuterContainer;
    using Base = IContainer;

    namespace Ops {

        Trace::Trace(IndexPairList indices) noexcept // NOLINT(bugprone-exception-escape)
            : _indices{std::move(indices)} {
            for (auto& elem : _indices) {
                _idxSet.insert(elem.first);
                _idxSet.insert(elem.second);
            }
        }

        Base* Trace::operator()(STensor& a) {
            auto newdimlabs = getNewIndexLabels(a);
            auto strides = a.getIndexing().processShifts(_indices, IndexPairMember::Both);
            IndexList inners(a.dims().size() - newdimlabs.first.size());
            vector<bool> innerAdds(inners.size(), false);
            if (!newdimlabs.first.empty()) {
                auto newsparse = makeEmptySparse(newdimlabs.first, newdimlabs.second);
                auto* result = static_cast<STensor*>(newsparse.release());
                for (const auto& elem : a) {
                    auto pospairs =
                        a.getIndexing().splitPosition(elem.first, get<0>(strides), get<1>(strides), inners, innerAdds);
                    if (pospairs == numeric_limits<size_t>::max()) {
                        continue;
                    }
                    (*result)[pospairs] += elem.second;
                }
                return static_cast<Base*>(result);
            }
            auto newscalar = makeEmptyScalar();
            IContainer::ElementType res = 0.;
            for (const auto& elem : a) {
                auto pospairs =
                    a.getIndexing().splitPosition(elem.first, get<0>(strides), get<1>(strides), inners, innerAdds);
                if (pospairs == numeric_limits<size_t>::max()) {
                    continue;
                }
                res += elem.second;
            }
            newscalar->element({}) = res;
            return newscalar.release();
        }

        Base* Trace::operator()(VTensor& a) {
            auto newdimlabs = getNewIndexLabels(a);
            const auto& aIndexing = a.getIndexing();
            auto strides = aIndexing.getInnerOuterStrides(_indices, aIndexing.strides(), true);
            if (!newdimlabs.first.empty()) {
                auto newvect = makeEmptyVector(newdimlabs.first, newdimlabs.second);
                auto* result = static_cast<VTensor*>(newvect.release());
                for (size_t i = 0; i < a.numValues(); ++i) {
                    if (isZero(a[i])) {
                        continue;
                    }
                    auto pospairs = LabeledIndexing<SequentialIndexing>::splitPosition(i, strides);
                    if (pospairs.second == 0) {
                        (*result)[pospairs.first] += a[i];
                    }
                }
                return static_cast<Base*>(result);
            }
            auto newscalar = makeEmptyScalar();
            IContainer::ElementType res = 0.;
            for (size_t i = 0; i < a.numValues(); ++i) {
                if (isZero(a[i])) {
                    continue;
                }
                auto pospairs = LabeledIndexing<SequentialIndexing>::splitPosition(i, strides);
                if (pospairs.second == 0) {
                    res += a[i];
                }
            }
            newscalar->element({}) = res;
            return newscalar.release();
        }


        Base* Trace::operator()(OTensor& a) {
            // assume no traces of subtensors
            for (auto& elem : _indices) {
                UNUSED(elem.first);
                ASSERT(a.getIndexing().getElementIndex(elem.first).first !=
                       a.getIndexing().getElementIndex(elem.second).first);
            }
            if (a.getIndexing().numSubIndexing() == 2) {
                // specialize for the standard outersquare case:
                // just one dot and not an outer anymore
                IndexPairList contraction;
                contraction.reserve(_indices.size());
                auto divider = static_cast<IndexType>(a.getIndexing().getSubIndexing(0).rank());
                for (auto& elem : _indices) {
                    if (elem.first < elem.second) {
                        contraction.emplace_back(elem.first, elem.second - divider);
                    } else {
                        contraction.emplace_back(elem.second, elem.first - divider);
                    }
                }
                stable_sort(contraction.begin(), contraction.end(),
                            [](const IndexPair& x, const IndexPair& y) -> bool { return x.second < y.second; });
                auto newdimlabs = getNewIndexLabels(a);
                // Outer is going to be collapsed
                TensorData result;
                if (newdimlabs.first.empty()) {
                    result = makeEmptyScalar();
                } else {
                    result = makeEmptySparse(newdimlabs.first, newdimlabs.second);
                }
                SharedTensorData current_entry;
                Ops::Sum summer{};
                for (auto& elem : a) {
                    Ops::Dot dotter{contraction, {elem[0].second, elem[1].second}};
                    current_entry = calc2(elem[0].first, *elem[1].first, dotter, "dot_outertrace");
                    if (result->rank() == 0) {
                        result->element({}) += current_entry->element({});
                    } else {
                        result = calc2(std::move(result), *current_entry, summer, "sum_outertrace");
                    }
                }
                return result.release();
            }
            /// @todo IMPLEMENT
            throw Error("Unimplemented tr O for general case");
            // // everything's a dot find which ones
            // map<IndexPair, IndexPairList, greater<IndexPair>> contractions;
            // for (auto& elem : _indices) {
            //     auto tmpLeft = a.getIndexing().getElementIndex(elem.first);
            //     auto tmpRight = a.getIndexing().getElementIndex(elem.second);
            //     if (tmpLeft.first > tmpRight.first) {
            //         contractions[{tmpLeft.first, tmpRight.first}].push_back(
            //             {tmpLeft.second, tmpRight.second});
            //     } else {
            //         contractions[{tmpRight.first, tmpLeft.first}].push_back(
            //             {tmpRight.second, tmpLeft.second});
            //     }
            // }
            // // be lazy: do one set of contractions and then call trace again
            // auto it = contractions.begin();
            // for (auto& elem : a) {
            //     Ops::Dot dotter{it->second, {elem.second[it->first.second], elem.second[it->first.first]}};
            //     elem.first[it->first.first] = calc2(elem.first[it->first.second],
            //                                         *elem.first[it->first.first], dotter, "dot_outertrace");
            //     elem.second[it->first.first] = false;
            //     elem.first.erase(elem.first.begin() + it->first.second);
            //     elem.second.erase(elem.second.begin() + it->first.second);
            // }
            // // now fix indexing before recursively calling trace
        }

        Base* Trace::error(Base& /*unused*/) {
            throw Error("Invalid data types for tensor Trace");
        }

        std::pair<IndexList, LabelsList> Trace::getNewIndexLabels(const IContainer& original) const {
            IndexList resultD = original.dims();
            LabelsList resultL = original.labels();
            for (auto elem : reverse_range(_idxSet)) {
                resultD.erase(resultD.begin() + elem);
                resultL.erase(resultL.begin() + elem);
            }
            return make_pair(resultD, resultL);
        }

        IndexList Trace::reducedIndex(const IndexList& a) const {
            IndexList result = a;
            for (auto elem : reverse_range(_idxSet)) {
                result.erase(result.begin() + elem);
            }
            return result;
        }

    } // namespace Ops


} // namespace Hammer::MultiDimensional
