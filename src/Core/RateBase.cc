///
/// @file  RateBase.cc
/// @brief Hammer base rate class
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <string>

#include "Hammer/RateBase.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Integrator.hh"

using namespace std;

namespace Hammer {

    RateBase::RateBase() : _nPoints{_integ.getPointsNumber()} {
    }

    // RateBase::RateBase(const RateBase& other)
    //     : ParticleData(other),
    //       SettingsConsumer(other),
    //       _integ(other._integ),
    //       _mPSRanges(other._mPSRanges),
    //       _nPoints(other._nPoints) {
    //     for (auto& elem : other._tensorList) {
    //         addTensor(Tensor{elem});
    //     }
    // }

    // RateBase& RateBase::operator=(const RateBase& other) {
    //     ParticleData::operator=(other);
    //     SettingsConsumer::operator=(other);
    //     _integ = other._integ;
    //     _mPSRanges = other._mPSRanges;
    //     _nPoints = other._nPoints;
    //     _tensorList.clear();
    //     for (auto& elem : other._tensorList) {
    //         addTensor(Tensor{elem});
    //     }
    //     return *this;
    // }


    // RateBase::~RateBase() {
    //     _mPSRanges.clear();
    //     _tensorList.clear();
    // }

    void RateBase::init() {
        // initSettings();
        ///@todo anything else?
    }

    Tensor& RateBase::getTensor() {
        return _tensorList[_signatureIndex];
    }

    string RateBase::getPdgIdString() const {
        auto sig = _signatures[_signatureIndex];
        string ds;
        for (auto& daughter : sig.daughters) {
            ds += " ";
            ds += to_string(daughter);
        }
        return to_string(sig.parent) + " ->" + ds;
    }

    pair<PdgId, vector<PdgId>> RateBase::getPdgIds() const {
        auto sig = _signatures[_signatureIndex];
        return make_pair(sig.parent, sig.daughters);
    }

    const Tensor& RateBase::getTensor() const {
        return _tensorList[_signatureIndex];
    }

    const IntegrationBoundaries& RateBase::getIntegrationBoundaries() const {
        return _mPSRanges[_signatureIndex];
    }

    void RateBase::addIntegrationBoundaries(const IntegrationBoundaries& boundaries) {
        // _mPSRanges.push_back({std::pow(qMin, 2.), std::pow(qMax, 2.)});
        _mPSRanges.push_back(boundaries);
    }

    void RateBase::eval(const Particle& /*parent*/, const ParticleList& /*daughters*/,
                        const ParticleList& /*references*/) {
    }

    void RateBase::defineSettings() {
    }

    const EvaluationGrid& RateBase::getEvaluationPoints() const {
        return _integ.getEvaluationPoints();
    }

    void RateBase::calcTensor() {
        const auto& boundaries = getIntegrationBoundaries();
        if (!boundaries.empty()) {
            _integ.applyRange(boundaries);
            const EvaluationGrid& points = _integ.getEvaluationPoints();
            const EvaluationWeights& weights = _integ.getEvaluationWeights();
            Tensor& t = getTensor();
            MSG_INFO("Integrating rate " + getPdgIdString());
            t.clearData();
            for (IndexType i = 0; i < static_cast<IndexType>(points.size()); ++i) {
                t.addAt(evalAtPSPoint(points[i]) * weights[i], INTEGRATION_INDEX, i);
            }
        } else {
            Tensor& t = getTensor();
            t = evalAtPSPoint({});
        }
    }

    void RateBase::addTensor(Tensor&& tensor) {
        _tensorList.push_back(std::move(tensor));
    }

    Log& RateBase::getLog() {
        return Log::getLog("Hammer.RateBase");
    }

} // namespace Hammer
