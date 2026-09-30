///
/// @file  Utils.cc
/// @brief Hammer math utilities class
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/Math/Utils.hh"

using namespace std;

namespace Hammer {

    double compareVals(const double val1, const double val2) {
        if (isZero(val2)) {
            return (val1 + 1.) / (val2 + 1.);
        }
        return val1 / val2;
    }

    complex<double> compareVals(const complex<double>& val1, const complex<double>& val2) {
        if (isZero(val2)) {
            return (val1 + 1.) / (val2 + 1.);
        }
        return val1 / val2;
    }

} // namespace Hammer
