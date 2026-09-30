///
/// @file  AmplitudeBase.cc
/// @brief Hammer base amplitude class
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/AmplitudeBase.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Math/MultiDim/VectorContainer.hh"
#include "Hammer/Math/MultiDim/ScalarContainer.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    vector<complex<double>> AmplitudeBase::getWCVectorFromDict(const map<string, complex<double>>& wcDict) const {
        vector<complex<double>> result(_mWCNames.size());
        //        for (size_t pos = 0ul; pos < _mWCNames.size(); ++pos) {
        //            auto it = wcDict.find(_mWCNames[pos]);
        //            if (it != wcDict.end()) {
        //                result[pos] = it->second;
        //            }
        //        }
        for (auto it = wcDict.begin(); it != wcDict.end(); ++it) {
            auto itwc = find(_mWCNames.begin(), _mWCNames.end(), it->first);
            if (itwc != _mWCNames.end()) {
                auto pos = distance(_mWCNames.begin(), itwc);
                result[static_cast<size_t>(pos)] = it->second;
            } else {
                MSG_WARNING("'" + it->first + "' does not belong to the " + _mWCPrefix +
                            " set of Wilson coefficients. Ignoring.");
            }
        }
        return result;
    }

    map<string, complex<double>> AmplitudeBase::getDictFromWCVector(const vector<complex<double>>& wcVect) const {
        map<string, complex<double>> result;
        ASSERT(_mWCNames.size() == wcVect.size());
        for (size_t pos = 0ul; pos < _mWCNames.size(); ++pos) {
            if (!isZero(wcVect[pos])) {
                result.insert({_mWCNames[pos], wcVect[pos]});
            }
        }
        return result;
    }

    vector<complex<double>> AmplitudeBase::getWCVectorFromSettings(WTerm what) const {
        WTerm oldwhat = const_cast<AmplitudeBase*>(this)->setWeightTerm(what);
        vector<complex<double>> result(_mWCNames.size());
        for (size_t pos = 0ul; pos < _mWCNames.size(); ++pos) {
            result[pos] = *getSetting<complex<double>>(_mWCPrefix, _mWCNames[pos]);
        }
        const_cast<AmplitudeBase*>(this)->setWeightTerm(oldwhat);
        return result;
    }

    vector<complex<double>> AmplitudeBase::resetWCSettings(WTerm what) {
        vector<complex<double>> result(_mWCNames.size());
        for (size_t pos = 0ul; pos < _mWCNames.size(); ++pos) {
            getSettingsHandler()->resetSetting(_mWCPrefix, _mWCNames[pos], what);
            result[pos] =
                *(getSettingsHandler()->getNamedSettingValue<complex<double>>(_mWCPrefix, _mWCNames[pos], what));
        }
        return result;
    }

    void AmplitudeBase::updateWCSettings(const vector<complex<double>>& values, WTerm what) {
        this->updateVectorOfSettings(values, _mWCNames, _mWCPrefix, what);
    }

    void AmplitudeBase::updateWCSettings(const map<string, complex<double>>& values, WTerm what) {
        this->updateVectorOfSettings(values, _mWCPrefix, what);
    }

    map<string, complex<double>> AmplitudeBase::retrieveWCSettings(WTerm what) const {
        return this->retrieveVectorOfSettings<complex<double>>(_mWCNames, _mWCPrefix, what);
    }

    void AmplitudeBase::updateWCTensor(vector<complex<double>> values, MD::SharedTensorData& data) const {
        ASSERT(values.size() == _mWCNames.size());
        preProcessWCValues(values);
        if (!data || data->rank() == 0) {
            data = MD::SharedTensorData{
                MD::makeVector({static_cast<uint16_t>(values.size())}, {_mWCLabel}, values).release()};
        } else {
            for (IndexType i = 0; i < static_cast<IndexType>(values.size()); ++i) {
                data->element({i}) = values[i];
            }
        }
    }

    pair<string, IndexLabel> AmplitudeBase::getWCInfo() const {
        return make_pair(_mWCPrefix, _mWCLabel);
    }

    void AmplitudeBase::createWCProjectionTensor(const vector<map<string, complex<double>>>& subspace,
                                                 const IndexLabel label, const MD::SharedTensorData& origin,
                                                 MD::SharedTensorData& proj) const {
        // o + sum_j x_j (sum_i wc_i^j * v_i^j)
        vector<IndexType> dims{
            {static_cast<IndexType>(_mWCNames.size()),
             static_cast<IndexType>(subspace.size() + 1)}}; // size is + 1 to include central values in zeroth component
        proj = MD::SharedTensorData{MD::makeEmptySparse(dims, {_mWCLabel, label})};
        auto isProjection = (origin->rank() == 2);
        for (IndexType i = 0; i < origin->dims()[0]; ++i) {
            auto elem = isProjection ? origin->element({i, 0}) : origin->element({i});
            if (!isZero(elem)) {
                proj->element({i, 0}) = elem;
            }
        }
        for (IndexType j = 0; j < static_cast<IndexType>(subspace.size()); ++j) {
            if (subspace[j].empty()) {
                continue;
            }
            auto vec = getWCVectorFromDict(subspace[j]);
            preProcessWCValues(vec);
            for (IndexType i = 0; i < static_cast<IndexType>(vec.size()); ++i) {
                if (!isZero(vec[i])) {
                    proj->element({i, static_cast<IndexType>(j + 1u)}) = vec[i];
                }
            }
        }
    }


    void AmplitudeBase::updateWCProjectionVector(int index, const map<string, complex<double>>& direction,
                                                 MD::SharedTensorData& proj) const {
        auto vec = getWCVectorFromDict(direction);
        updateWCProjectionVector(index, vec, proj);
    }

    void AmplitudeBase::updateWCProjectionVector(int index, const vector<complex<double>>& direction,
                                                 MD::SharedTensorData& proj) const {
        // proj = MD::SharedTensorData{original->clone().release()};
        auto vec = direction;
        preProcessWCValues(vec);
        ASSERT(proj->rank() > 0 && proj->rank() <= 2); // NOLINT(readability-simplify-boolean-expr)
        for (IndexType i = 0; i < proj->dims()[0]; ++i) {
            vector<IndexType> indices =
                (proj->rank() == 2) ? vector<IndexType>{i, static_cast<IndexType>(index + 1)} : vector<IndexType>{i};
            if (!isZero(vec[i] - proj->element(indices))) {
                proj->element(indices) = vec[i];
            }
        }
    }

    void AmplitudeBase::preProcessWCValues(vector<complex<double>>& /*unused*/, bool /*unused*/) const {
    }

    void AmplitudeBase::init() {
        // initSettings();
        ///@todo anything else?
    }

    Tensor& AmplitudeBase::getTensor() {
        return _tensorList[_signatureIndex];
    }

    const Tensor& AmplitudeBase::getTensor() const {
        return _tensorList[_signatureIndex];
    }

    bool AmplitudeBase::setSignatureIndex(size_t idx) {
        bool res = ParticleData::setSignatureIndex(idx);
        if (res) {
            this->updateWilsonCoeffLabelPrefix();
            return true;
        }
        MSG_WARNING("Unable to update the signature index");
        return false;
    }

    void AmplitudeBase::updateWilsonCoeffLabelPrefix() {
        _mWCLabel = NONE;
        _mWCPrefix = "None";
    }

    size_t AmplitudeBase::multiplicityFactor() const {
        return _multiplicity;
    }

    Log& AmplitudeBase::getLog() {
        return Log::getLog("Hammer.AmplitudeBase");
    }

    void AmplitudeBase::addTensor(Tensor&& tensor) {
        _tensorList.push_back(std::move(tensor));
    }


} // namespace Hammer
