///
/// @file  OuterSquare.cc
/// @brief Tensor outer square algorithm
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/Math/MultiDim/Ops/OuterSquare.hh"
#include "Hammer/Math/MultiDim/VectorContainer.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"
#include "Hammer/Math/MultiDim/OuterContainer.hh"
#include "Hammer/Exceptions.hh"

using namespace std;

namespace Hammer::MultiDimensional {

    using VTensor = VectorContainer;
    using STensor = SparseContainer;
    using OTensor = OuterContainer;
    using Base = IContainer;

    namespace Ops {

        Base* OuterSquare::operator()(OTensor& a) {
            TensorData res{a.clone()};
            auto* tmp = static_cast<OuterContainer*>(res.get());
            size_t orig_size = tmp->numAddends();
            auto newvec = tmp->_data[0];
            auto orig_factors = static_cast<ptrdiff_t>(newvec.size());
            for (auto& i : newvec) {
                i.second = !i.second;
            }
            for (size_t i = 0; i < orig_size; ++i) {
                tmp->_data[i].insert(tmp->_data[i].end(), newvec.begin(), newvec.end());
            }
            tmp->_accessors.clear();
            ptrdiff_t start = 0ul;
            ptrdiff_t finish = 0ul;
            for (auto& elem : tmp->_data[0]) {
                finish += static_cast<ptrdiff_t>(elem.first->rank());
                tmp->_accessors.emplace_back(start,finish);
                start = finish;
            }
            vector<IndexList> tmpdims;
            vector<LabelsList> tmplabels;
            for (auto& it : tmp->_data[0]) {
                auto entrydims = (it.first)->dims();
                auto entrylabs = (it.first)->labels();
                if (it.second) {
                    entrylabs = flipListOfLabels(entrylabs);
                }
                tmpdims.push_back(entrydims);
                tmplabels.push_back(entrylabs);
            }
            tmp->_indexing = BlockIndexing{tmpdims, tmplabels};
            for (size_t i = 0; i < orig_size; ++i) {
                for (size_t j = 1; j < orig_size; ++j) {
                    OuterContainer::EntryType newentryvec(tmp->_data[i].begin(), tmp->_data[i].begin() + orig_factors);
                    newentryvec.insert(newentryvec.end(), tmp->_data[j].begin(), tmp->_data[j].begin() + orig_factors);
                    for (auto k = static_cast<size_t>(orig_factors); k < newentryvec.size(); ++k) {
                        newentryvec[k].second = !newentryvec[k].second;
                    }
                    tmp->addTerm(newentryvec);
                }
            }
            return static_cast<Base*>(res.release());
        }

        Base* OuterSquare::operator()(Base& a) {
            TensorData tmp = a.clone();
            auto* res = new OuterContainer{std::move(tmp)};
            return static_cast<Base*>(res);
        }

        Base* OuterSquare::error(Base& /*unused*/) {
            throw Error("Invalid data types for tensor OuterSquare");
        }

    } // namespace Ops


} // namespace Hammer::MultiDimensional
