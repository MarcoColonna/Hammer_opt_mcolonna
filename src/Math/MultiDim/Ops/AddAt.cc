///
/// @file  AddAt.cc
/// @brief Sub-tensor block insertion algorithm
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/Math/MultiDim/Ops/AddAt.hh"
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

        AddAt::AddAt(IndexType position, IndexType coord) noexcept : _position{position}, _coord{coord} {
        }

        Base* AddAt::operator()(VTensor& a, const VTensor& b) const {
            for (size_t index = 0; index < b.numValues(); ++index) {
                const PositionType newpos = a.getIndexing().extendPosition(index, _position, _coord);
                a[newpos] += b[index];
            }
            return static_cast<Base*>(&a);
        }

        Base* AddAt::operator()(VTensor& a, const STensor& b) const {
            for (auto it = b.begin(); it != b.end(); ++it) {
                const PositionType newpos =
                    a.getIndexing().extendPosition(b.getIndexing().alignedPosToPos(it->first), _position, _coord);
                a[newpos] += it->second;
            }
            return static_cast<Base*>(&a);
        }

        Base* AddAt::operator()(STensor& a, const STensor& b) const {
            for (const auto& it : b) {
                const PositionType newpos = a.getIndexing().extendAlignedPosition(it.first, _position, _coord);
                a[newpos] += it.second;
            }
            return static_cast<Base*>(&a);
        }

        Base* AddAt::operator()(STensor& a, const VTensor& b) const {
            auto tmpStride = _position > 0 ? b.getIndexing().stride(static_cast<IndexType>(_position - 1))
                                           : b.getIndexing().numValues();
            for (size_t index = 0; index < b.numValues(); ++index) {
                const PositionType newpos = a.getIndexing().extendPosition(index, tmpStride, _position, _coord);
                a[newpos] += b[index];
            }
            return static_cast<Base*>(&a);
        }

        Base* AddAt::operator()(OTensor& a, const Base& b) const {
            return error(a, b);
        }

        Base* AddAt::operator()(VTensor& a, const OTensor& b) const {
            BruteForceIterator bf{b.dims()};
            SequentialIndexing seq{b.dims()};
            for (const auto& el: bf) {
                const PositionType newpos = a.getIndexing().extendPosition(seq.indicesToPos(el), _position, _coord);
                a[newpos] += b.value(el);
            }
            return static_cast<Base*>(&a);
        }

        Base* AddAt::operator()(STensor& a, const OTensor& b) const {
            BruteForceIterator bf{b.dims()};
            AlignedIndexing ali{b.dims()};
            for (const auto& el: bf) {
                const PositionType newpos =
                    a.getIndexing().extendAlignedPosition(ali.indicesToPos(el), _position, _coord);
                a[newpos] += b.value(el);
            }
            return static_cast<Base*>(&a);
        }

        Base* AddAt::operator()(Base& a, const Base& b) const {
            return error(a, b);
        }

        Base* AddAt::error(Base& /*unused*/, const Base& /*unused*/) {
            throw Error("Invalid data types for tensor AddAt");
        }
    } // namespace Ops


} // namespace Hammer::MultiDimensional
