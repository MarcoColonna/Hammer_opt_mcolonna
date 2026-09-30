///
/// @file  Divide.cc
/// @brief Tensor element-wise division algorithm
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/Math/MultiDim/Ops/Divide.hh"
#include "Hammer/Math/MultiDim/VectorContainer.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"
#include "Hammer/Math/MultiDim/OuterContainer.hh"
#include "Hammer/Math/MultiDim/BruteForceIterator.hh"
#include "Hammer/Exceptions.hh"
#include "Hammer/Math/Utils.hh"

using namespace std;

namespace Hammer::MultiDimensional {

    using VTensor = VectorContainer;
    using STensor = SparseContainer;
    using OTensor = OuterContainer;
    using Base = IContainer;

    namespace Ops {

        IContainer* Divide::operator()(VTensor& a, const VTensor& b) {
            auto ita = a.begin();
            auto itb = b.begin();
            for (; ita != a.end(); ++itb, ++ita) {
                *ita /= *itb;
            }
            return static_cast<IContainer*>(&a);
        }

        IContainer* Divide::operator()(STensor& a, const STensor& b) {
            auto first1 = a.begin();
            auto first2 = b.begin();
            while (first1 != a.end()) {
                if (first2 == b.end()) {
                    a.erase(first1, a.end());
                }
                if (first1->first < first2->first) {
                    first1->second = isZero(first1->second) ? 0. : numeric_limits<double>::infinity();
                    ++first1;
                } else {
                    if (!(first2->first < first1->first)) {
                        first1->second /= first2->second;
                        ++first1;
                    }
                    ++first2;
                }
            }
            return static_cast<IContainer*>(&a);
        }

        IContainer* Divide::operator()(VTensor& a, const STensor& b) {
            size_t n = 0;
            auto oldsize = a.numValues();
            auto first2 = b.begin();
            while (n < oldsize) {
                if (first2 == b.end()) {
                    for (size_t m = n; m < oldsize; ++m) {
                        a[m] = isZero(a[m]) ? 0. : numeric_limits<double>::infinity();
                    }
                    break;
                }
                if (n < b.getIndexing().alignedPosToPos(first2->first)) {
                    a[n] = isZero(a[n]) ? 0. : numeric_limits<double>::infinity();
                    ++n;
                } else {
                    if (!(b.getIndexing().alignedPosToPos(first2->first) < n)) {
                        a[n] /= first2->second;
                        ++n;
                    }
                    ++first2;
                }
            }
            return static_cast<IContainer*>(&a);
        }

        IContainer* Divide::operator()(STensor& a, const VTensor& b) {
            for (auto& elem : a) {
                elem.second /= b[a.getIndexing().alignedPosToPos(elem.first)];
            }
            return static_cast<IContainer*>(&a);
        }

        IContainer* Divide::operator()(OTensor& a, const Base& b) {
            return error(a, b);
        }

        IContainer* Divide::operator()(Base& a, const Base& b) {
            const BruteForceIteratorRange bf{a.dims()};
            for (const auto& elem : bf) {
                a.element(elem) /= b.element(elem);
            }
            return &a;
        }

        IContainer* Divide::error(Base& /*unused*/, const Base& /*unused*/) {
            throw Error("Invalid data types for tensor Divide");
        }

    } // namespace Ops


} // namespace Hammer::MultiDimensional
