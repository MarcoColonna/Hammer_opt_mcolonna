///
/// @file  Differentiator.hh
/// @brief Numerical derivatives class
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_MATH_DIFFERENTIATOR
#define HAMMER_MATH_DIFFERENTIATOR

#include <functional>
#include <limits>
#include <utility>

#include "Hammer/Exceptions.hh"

namespace Hammer {

    static double derEps = 1.E-4;

    template <unsigned int Order>
    double derLeft(const std::function<double(double)>& /*unused*/, double /*unused*/) {
        throw Error("Not implemented;");
    }

    template <unsigned int Order>
    double derRight(const std::function<double(double)>& /*unused*/, double /*unused*/) {
        throw Error("Not implemented;");
    }

    template <unsigned int Order>
    double derCenter(const std::function<double(double)>& /*unused*/, double /*unused*/) {
        throw Error("Not implemented;");
    }

    template <>
    inline double derLeft<1>(const std::function<double(double)>& f, double value) {
        return (4. * f(value - derEps) - f(value - 2. * derEps) - 3. * f(value)) / (-2. * derEps);
    }

    template <>
    inline double derRight<1>(const std::function<double(double)>& f, double value) {
        return (4. * f(value + derEps) - f(value + 2. * derEps) - 3. * f(value)) / (2. * derEps);
    }

    template <>
    inline double derCenter<1>(const std::function<double(double)>& f, double value) {
        return (f(value + derEps) - f(value - derEps)) / (2. * derEps);
    }

    template <>
    inline double derLeft<2>(const std::function<double(double)>& f, double value) {
        return (4. * f(value - 2. * derEps) - f(value - 3. * derEps) - 5. * f(value - derEps) + 2. * f(value)) /
               (derEps * derEps);
    }

    template <>
    inline double derRight<2>(const std::function<double(double)>& f, double value) {
        return (4. * f(value + 2. * derEps) - f(value + 3. * derEps) - 5. * f(value + derEps) + 2. * f(value)) /
               (derEps * derEps);
    }

    template <>
    inline double derCenter<2>(const std::function<double(double)>& f, double value) {
        return (f(value + derEps) + f(value - derEps) - 2. * f(value)) / (derEps * derEps);
    }

    template <unsigned int Order>
    inline double der(const std::function<double(double)>& f, double value,
                      const std::pair<double, double>& boundaries = {-std::numeric_limits<double>::max(),
                                                                     std::numeric_limits<double>::max()}) {
        if (value <= boundaries.first) {
            return derRight<Order>(f, boundaries.first);
        }
        if (value >= boundaries.second) {
            return derLeft<Order>(f, boundaries.second);
        }
        return derCenter<Order>(f, value);
    }

} // namespace Hammer

#endif
