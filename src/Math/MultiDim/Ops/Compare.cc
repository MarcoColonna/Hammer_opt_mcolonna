///
/// @file  Compare.cc
/// @brief Sub-tensor block insertion algorithm
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/Math/MultiDim/Ops/Compare.hh"
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

        bool Compare::operator()(const VTensor& a, const VTensor& b) {
            return a.compare(b);
        }

        bool Compare::operator()(const VTensor& a, const STensor& b) {
            bool result = checkStructure(a, b);
            auto itb = b.begin();
            for (PositionType i = 0; result && i < a.numValues(); ++i) {
                if (itb->first != b.getIndexing().posToAlignedPos(i)) {
                    result &= isZero(a[i]);
                } else {
                    result &= isZero(a[i] - itb->second);
                    ++itb;
                }
            }
            return result;
        }

        bool Compare::operator()(const STensor& a, const STensor& b) {
            return a.compare(b);
        }

        bool Compare::operator()(const STensor& a, const VTensor& b) {
            return operator()(b, a);
        }

        bool Compare::operator()(const Base& a, const Base& b) {
            bool result = checkStructure(a, b);
            if (a.rank() == 0) {
                result &= isZero(a.element({}) - b.element({}));
            } else {
                BruteForceIteratorRange bf{a.dims()};
                for (auto elem : bf) {
                    result &= isZero(a.element(elem) - b.element(elem));
                    if (!result) {
                        break;
                    }
                }
            }
            return result;
        }

        bool Compare::error(const Base& /*unused*/, const Base& /*unused*/) {
            throw Error("Invalid data types for tensor Compare");
        }

        bool Compare::checkStructure(const Base& a, const Base& b) {
            bool result = (a.rank() == b.rank());
            if (result) {
                result &= a.isSameShape(b);
            }
            if (result) {
                result &= (a.labels() == b.labels());
            }
            return result;
        }

    } // namespace Ops


} // namespace Hammer::MultiDimensional
