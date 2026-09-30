///
/// @file  PhaseSpace.hh
/// @brief Phase space integrals
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_MATH_PHASESPACE
#define HAMMER_MATH_PHASESPACE

#include <vector>

namespace Hammer {

    /// @brief
    /// @param[in] mass
    /// @param[in] masses
    /// @return
    double phaseSpaceN(double mass, const std::vector<double>& masses);


    /// @brief
    /// @param[in] mass
    /// @param[in] masses
    /// @return
    double phaseSpaceNBody(double mass, const std::vector<double>& masses);

} // namespace Hammer

#endif
