///
/// @file  Integrator.cc
/// @brief Hammer numerical integration classes
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/Math/Integrator.hh"
#include "Hammer/Exceptions.hh"

using namespace std;

namespace Hammer {

    double integrate(std::function<double(double)>& func, double low, double high) {

        double value = 0.;
        static constexpr uint8_t N = 19;

        array<double, N> xi;
        array<double, N> wi;

        xi[0] = 0.0;
        wi[0] = 0.16105444984878362;
        xi[1] = 0.16035862813102394;
        wi[1] = 0.15896882598519965;
        xi[2] = -xi[1];
        wi[2] = wi[1];
        xi[3] = 0.3165642266117781;
        wi[3] = 0.1527663007284602;
        xi[4] = -xi[3];
        wi[4] = wi[3];
        xi[5] = 0.46457038687717755;
        wi[5] = 0.1426055640976579;
        xi[6] = -xi[5];
        wi[6] = wi[5];
        xi[7] = 0.6005459122435153;
        wi[7] = 0.1287567549009922;
        xi[8] = -xi[7];
        wi[8] = wi[7];
        xi[9] = 0.7209655043062628;
        wi[9] = 0.11156236184219351;
        xi[10] = -xi[9];
        wi[10] = wi[9];
        xi[11] = 0.8227150489105877;
        wi[11] = 0.09149349482553046;
        xi[12] = -xi[11];
        wi[12] = wi[11];
        xi[13] = 0.9031559802279142;
        wi[13] = 0.06904552773838758;
        xi[14] = -xi[13];
        wi[14] = wi[13];
        xi[15] = 0.9602078316107117;
        wi[15] = 0.044807508218039215;
        xi[16] = -xi[15];
        wi[16] = wi[15];
        xi[17] = 0.9924070086164836;
        wi[17] = 0.019469784466159833;
        xi[18] = -xi[17];
        wi[18] = wi[17];

        for (uint8_t i = 0; i < N; ++i) {
            double yi = (high - low) * xi[i] / 2.0 + (high + low) / 2.0;
            double dGi = func(yi);
            value += wi[i] * dGi;
        }

        value *= (high - low) / 2.0;

        return value;
    }

    Integrator::Integrator()
        : _basePoints{{0., 0.16035862813102394, 0.3165642266117781, 0.46457038687717755, 0.6005459122435153,
                       0.7209655043062628, 0.8227150489105877, 0.9031559802279142, 0.9602078316107117,
                       0.9924070086164836}},
          _baseWeights{{0.16105444984878362, 0.15896882598519965, 0.1527663007284602, 0.1426055640976579,
                        0.1287567549009922, 0.11156236184219351, 0.09149349482553046, 0.06904552773838758,
                        0.044807508218039215, 0.019469784466159833}} {
        for (size_t j = 1; j < 10; ++j) {
            _basePoints.push_back(-_basePoints[j]);
            _baseWeights.push_back(_baseWeights[j]);
        }
    }

    Integrator::~Integrator() {
        _evaluationPoints.clear();
        _evaluationWeights.clear();
        _basePoints.clear();
        _baseWeights.clear();
    }

    const EvaluationGrid& Integrator::getEvaluationPoints() const {
        return _evaluationPoints;
    }

    const EvaluationWeights& Integrator::getEvaluationWeights() const {
        return _evaluationWeights;
    }

    void Integrator::applyRange(IntegrationBoundaries&& boundaries) {
        auto tmpboundaries = std::move(boundaries);
        applyRange(tmpboundaries);
    }

    void Integrator::applyRange(const IntegrationBoundaries& boundaries) {
        auto totN = static_cast<size_t>(pow(_basePoints.size(), boundaries.size()));
        ASSERT(totN > 0);
        ASSERT(totN < numeric_limits<uint16_t>::max());
        ASSERT(!boundaries.empty());
        vector<double> tmp(boundaries.size(), 0.);
        _evaluationPoints.clear();
        _evaluationPoints.reserve(totN);
        _evaluationPoints.insert(_evaluationPoints.end(), totN, tmp);
        _evaluationWeights = EvaluationWeights(totN, 1.);
        size_t outerStride = _evaluationPoints.size();
        for (size_t i = 0; i < boundaries.size(); ++i) {
            outerStride /= (i == 0) ? 1 : _basePoints.size();
            size_t currentStride = outerStride / _basePoints.size();
            for (size_t jOut = 0; jOut < _evaluationPoints.size(); jOut += outerStride) {
                auto coords = _evaluationPoints[jOut];
                coords.resize(i);
                auto range = boundaries[i](coords);
                double sum = (range.second + range.first) / 2.;
                double diff = (range.second - range.first) / 2.;
                for (size_t k = 0; k < _basePoints.size(); ++k) {
                    size_t delta = k * currentStride;
                    for (size_t jIn = 0; jIn < currentStride; ++jIn) {
                        size_t pos = jOut + delta + jIn;
                        _evaluationWeights[pos] *= _baseWeights[k] * diff;
                        _evaluationPoints[pos][i] = sum + diff * _basePoints[k];
                    }
                }
            }
        }
    }

    uint16_t Integrator::getPointsNumber() const {
        return static_cast<uint16_t>(_basePoints.size());
    }

} // namespace Hammer
