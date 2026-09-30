#ifndef HAMMER_MATH_CONSTANTS
#define HAMMER_MATH_CONSTANTS

///
/// @file  Constants.hh
/// @brief Various numerical constants
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <boost/math/constants/constants.hpp>

#include "Hammer/Math/Units.hh"

namespace Hammer {

    static constexpr double pi = boost::math::constants::pi<double>();
    static constexpr double twoPi = boost::math::constants::two_pi<double>();
    static constexpr double halfPi = boost::math::constants::half_pi<double>();
    static constexpr double pi2 = boost::math::constants::pi_sqr<double>();
    static constexpr double pi3 = boost::math::constants::pi_cubed<double>();
    static constexpr double sqrt2 = boost::math::constants::root_two<double>();
    static constexpr double sqrt3 = boost::math::constants::root_three<double>();

    static constexpr double GFermi = 1.16638e-5 / GeV2;
    static constexpr double TwoSqTwoGFermi = 2. * sqrt2 * GFermi;
    static constexpr double FPion = 93 * MeV;

} // namespace Hammer

#endif
