///
/// @file  Integrator.hh
/// @brief Hammer numerical integration classes
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_MATH_Integrator_HH
#define HAMMER_MATH_Integrator_HH

#include <functional>
#include <vector>
#include <cstdint>

#include "Hammer/Math/Integrator.fhh"

namespace Hammer {

    /// @brief integration function based on gaussian method
    /// @param[in] func  the function to be integrated
    /// @param[in] low   the lower limit of integration
    /// @param[in] high  the higher limit of integration
    /// @return the integral
    double integrate(std::function<double(double)>& func, double low, double high);


    /// @brief Tensor integration class
    ///
    /// ...
    ///
    /// @ingroup Math
    class Integrator {

    public:

        Integrator();

        Integrator(const Integrator& other) = default;
        Integrator& operator=(const Integrator& other) = default;

        virtual ~Integrator();

        /// @brief returns the position of the evaluation points given an integration range
        /// @return the points where the function (tensor) should be evaluated
        [[nodiscard]] const EvaluationGrid& getEvaluationPoints() const;

        [[nodiscard]] const EvaluationWeights& getEvaluationWeights() const;

        /// @brief apply rescaling factor for the size of the interval from the canonical [-1,1] interval
        ///        to the actual one, \f$ (x_{max}-x_{min})/2 \f$
        /// @param[in] boundaries the lower limit of integration
        void applyRange(const IntegrationBoundaries& boundaries);

        void applyRange(IntegrationBoundaries&& boundaries);

        /// @brief the number of evaluation points
        [[nodiscard]] uint16_t getPointsNumber() const;

    private:

        std::vector<double> _basePoints;  ///< the evaluation points in the range [-1,1]
        std::vector<double> _baseWeights; ///< the evaluation weights for the range [-1,1]
        EvaluationGrid _evaluationPoints;
        EvaluationWeights _evaluationWeights;
    };

} // namespace Hammer

#endif
