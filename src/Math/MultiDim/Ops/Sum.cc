///
/// @file  Sum.cc
/// @brief Tensor sum algorithm
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <numeric>

#include "Hammer/Math/MultiDim/Ops/Sum.hh"
#include "Hammer/Math/MultiDim/Ops/Convert.hh"
#include "Hammer/Math/MultiDim/VectorContainer.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"
#include "Hammer/Math/MultiDim/OuterContainer.hh"
#include "Hammer/Math/MultiDim/BruteForceIterator.hh"
#include "Hammer/Exceptions.hh"

using namespace std;

namespace Hammer::MultiDimensional {

    using VTensor = VectorContainer;
    using STensor = SparseContainer;
    using OTensor = OuterContainer;
    using Base = IContainer;

    namespace Ops {

        IContainer* Sum::operator()(VTensor& a, const VTensor& b) {
            auto ita = a.begin();
            auto itb = b.begin();
            for (; ita != a.end(); ++itb, ++ita) {
                *ita += *itb;
            }
            return static_cast<IContainer*>(&a);
        }

        IContainer* Sum::operator()(STensor& a, const STensor& b) {
            for (const auto& elem : b) {
                a[elem.first] += elem.second;
            }
            return static_cast<IContainer*>(&a);
        }

        IContainer* Sum::operator()(VTensor& a, const STensor& b) {
            for (const auto& elem : b) {
                a[b.getIndexing().alignedPosToPos(elem.first)] += elem.second;
            }
            return static_cast<IContainer*>(&a);
        }

        IContainer* Sum::operator()(STensor& a, const VTensor& b) {
            auto* tmp = new VTensor{b};
            auto* result = this->operator()(*tmp, a);
            return result;
        }

        IContainer* Sum::operator()(OTensor& a, const OTensor& b) {
            if (a.numAddends() == 0) {
                auto res = b.clone();
                return res.release();
            }
            if (b.numAddends() == 0) {
                return static_cast<IContainer*>(&a);
            }
            size_t numTerms = a.numAddends() + b.numAddends();
            if (numTerms % 10 == 0) {
                if (a.shouldBeEvaluated()) {
                    auto* sparse = dynamic_cast<STensor*>((Convert{true})(a));
                    return (*this)(*sparse, b);
                }
            }
            auto tmp = b;
            a.reserve(numTerms);
            for (auto& elem : tmp._data) {
                a.addTerm(elem);
            }
            return static_cast<IContainer*>(&a);
        }

        IContainer* Sum::operator()(OTensor& a, const VTensor& b) {
            auto* tmp = new VTensor{b};
            auto* result = this->operator()(*tmp, a);
            return result;
        }

        IContainer* Sum::operator()(OTensor& a, const STensor& b) {
            auto* tmp = new STensor{b};
            auto* result = this->operator()(*tmp, a);
            return result;
        }

        IContainer* Sum::operator()(STensor& a, const OTensor& b) {
            if (b.numAddends() == 0) {
                return static_cast<IContainer*>(&a);
            }
            if (a.dataSize() == 0) {
                auto res = b.clone();
                return res.release();
            }
            IndexList chunkIndices(b.begin()->size());
            iota(chunkIndices.begin(), chunkIndices.end(), 0);
            for (const auto& elem : b) {
                OuterElemIterator bi{elem};
                OuterElemIterator ei = bi.end();
                while (bi != ei) {
                    a[b.getIndexing().buildFullPosition(bi, chunkIndices)] += *bi;
                    ++bi;
                }
            }
            return static_cast<IContainer*>(&a);
        }

        IContainer* Sum::operator()(Base& a, const Base& b) {
            BruteForceIterator bf{a.dims()};
            for (auto elem : bf) {
                a.element(elem) += b.element(elem);
            }
            return &a;
        }

        IContainer* Sum::error(Base& /*unused*/, const Base& /*unused*/) {
            throw Error("Invalid data types for tensor sum");
        }

    } // namespace Ops


} // namespace Hammer::MultiDimensional
