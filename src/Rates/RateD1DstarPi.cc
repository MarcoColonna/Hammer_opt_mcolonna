///
/// @file  RateD1DstarPi.cc
/// @brief \f$ D_1 \rightarrow D^* \pi \f$ total rate
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/Rates/RateD1DstarPi.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"
#include "Hammer/Tools/Pdg.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    RateD1DstarPi::RateD1DstarPi() {
        // Create tensor rank and dimensions
        IndexList dims{{2, 2}};
        string name{"RateD1DstarPi"};
        // auto& pdg = PID::instance();

        addProcessSignature(-PID::DSSD1, {PID::DSTARMINUS, PID::PIPLUS});
        addIntegrationBoundaries({});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD1, FF_DSSD1_HC})});

        addProcessSignature(-PID::DSSD1, {-PID::DSTAR, PID::PI0});
        addIntegrationBoundaries({});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD1, FF_DSSD1_HC})});

        addProcessSignature(PID::DSSD1MINUS, {-PID::DSTAR, PID::PIMINUS});
        addIntegrationBoundaries({});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD1, FF_DSSD1_HC})});

        addProcessSignature(PID::DSSD1MINUS, {PID::DSTARMINUS, PID::PI0});
        addIntegrationBoundaries({});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD1, FF_DSSD1_HC})});

        // strange
        addProcessSignature(PID::DSSDS1MINUS, {PID::DSSTARMINUS, PID::PI0});
        addIntegrationBoundaries({});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSDS1, FF_DSSDS1_HC})});

        addProcessSignature(PID::DSSDS1MINUS, {PID::DSTARMINUS, -PID::K0});
        addIntegrationBoundaries({});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSDS1, FF_DSSDS1_HC})});

        addProcessSignature(PID::DSSDS1MINUS, {-PID::DSTAR, PID::KMINUS});
        addIntegrationBoundaries({});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSDS1, FF_DSSDS1_HC})});

        setSignatureIndex();
    }

    Tensor RateD1DstarPi::evalAtPSPoint(const vector<double>& /*point*/) {
        auto labs = getTensor().labels();
        auto dimensions = getTensor().dims();
        Tensor result{"RateD1DstarPi", MD::makeEmptySparse(dimensions, labs)};

        const double Mdss = masses()[0];
        const double Mdss2 = Mdss * Mdss;
        const double Mds = masses()[1];
        const double Mds2 = Mds * Mds;
        const double Mp = masses()[2];
        const double Mp2 = Mp * Mp;

        const double Ep = (Mdss2 - Mds2 + Mp2) / (2 * Mdss);
        const double Pp = sqrt(Ep * Ep - Mp2);

        const double Pp4 = Pp * Pp * Pp * Pp;

        double RateNorm = Pp / (24. * pi * Mdss2 * FPion * FPion);

        // set non-zero tensor elements, basis S,D
        result.element({0, 0}) = (9 * Pp4) / 2.;
        result.element({1, 1}) = Pp4;

        result *= RateNorm;

        return result;
    }

} // namespace Hammer
