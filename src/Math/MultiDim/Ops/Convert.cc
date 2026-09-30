///
/// @file  Convert.cc
/// @brief Tensor storage type conversion algorithm
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <numeric>

#include "Hammer/Math/MultiDim/Ops/Convert.hh"
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

        Convert::Convert(bool destinationSparse) noexcept : _destIsSparse{destinationSparse} {
        }

        Base* Convert::operator()(VTensor& a) const {
            if (_destIsSparse) {
                return toSparse(a);
            }
            return static_cast<Base*>(&a);
        }

        Base* Convert::operator()(STensor& a) const {
            if (_destIsSparse) {
                return static_cast<Base*>(&a);
            }
            return toVector(a);
        }

        Base* Convert::operator()(OTensor& a) const {
            if (_destIsSparse) {
                return toSparse(a);
            }
            return toVector(a);
        }

        Base* Convert::operator()(Base& a) const {
            return error(a);
        }


        Base* Convert::toVector(STensor& a) {
            auto tmp = makeEmptyVector(a.dims(), a.labels());
            auto* res = static_cast<VTensor*>(tmp.get());
            for (auto elem : a) {
                res->_data[a.getIndexing().alignedPosToPos(elem.first)] = elem.second;
            }
            return tmp.release();
        }

        Base* Convert::toVector(OTensor& a) {
            auto tmp = makeEmptyVector(a.dims(), a.labels());
            auto* res = static_cast<VTensor*>(tmp.get());
            BruteForceIterator bf{a.dims()};
            for (const auto& el: bf) {
                auto val = a.value(el);
                if (!isZero(val)) {
                    res->element(el) = val;
                }
            }
            return tmp.release();
        }

        Base* Convert::toSparse(VTensor& a) {
            auto tmp = makeEmptySparse(a.dims(), a.labels());
            auto* res = static_cast<STensor*>(tmp.get());
            for (PositionType i = 0; i < a.numValues(); ++i) {
                if (!isZero(a[i])) {
                    (*res)[res->getIndexing().posToAlignedPos(i)] = a[i];
                }
            }
            return tmp.release();
        }

        Base* Convert::toSparse(OTensor& a) {
            auto tmp = makeEmptySparse(a.dims(), a.labels());
            auto* res = static_cast<STensor*>(tmp.get());
            IndexList chunkIndices(a.begin()->size());
            iota(chunkIndices.begin(), chunkIndices.end(), 0);
            for (auto& elem : a) {
                OuterElemIterator bi{elem};
                OuterElemIterator ei = bi.end();
                while (bi != ei) {
                    (*res)[a.getIndexing().buildFullPosition(bi, chunkIndices)] += *bi;
                    ++bi;
                }
            }
            return tmp.release();
        }

        Base* Convert::error(Base& /*unused*/) {
            throw Error("Invalid data types for tensor Convert");
        }

    } // namespace Ops


} // namespace Hammer::MultiDimensional
