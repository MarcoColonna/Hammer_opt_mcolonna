///
/// @file  FFLSPRVarBase.cc
/// @brief Hammer base class for LSPR form factors
//

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/FormFactors/LSPR/FFLSPRVarBase.hh"

using namespace std;

namespace Hammer {

    FFLSPRVarBase::FFLSPRVarBase() {
        setGroup("LSPRVar");
    }

    void FFLSPRVarBase::defineSettings() {

        defineAndAddErrSettings({"delta_s1", "delta_sp1", "delta_s11", "delta_s1p1", "delta_pc1", "delta_pb1"});

        FFLSPRBase::defineSettings(); // common parameters for LSPR are defined in FFLSPRBase

        // covariance eigenvector matrix
        // Row basis is s1,sp1, s11, s1p1, pc1, pb1
        vector<vector<double>> sliwmat{{1., 0., 0., 0., 0., 0.}, {0., 1., 0., 0., 0., 0.}, {0., 0., 1., 0., 0., 0.},
                                       {0., 0., 0., 1., 0., 0.}, {0., 0., 0., 0., 1., 0.}, {0., 0., 0., 0., 0., 1.}};
        addSetting<vector<vector<double>>>("sliwmatrix", sliwmat);
    }

} // namespace Hammer
