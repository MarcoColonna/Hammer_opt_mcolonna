///
/// @file  RateD1starDstarPi.cc
/// @brief \f$ D_1^* \rightarrow D^* \pi \f$ total rate
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/Rates/RateD1starDstarPi.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"
#include "Hammer/Tools/Pdg.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    RateD1starDstarPi::RateD1starDstarPi() {
        // Create tensor rank and dimensions
        IndexList dims{{1, 1}};
        string name{"RateD1starDstarPi"};
        // auto& pdg = PID::instance();

        addProcessSignature(-PID::DSSD1STAR, {PID::DSTARMINUS, PID::PIPLUS});
        addIntegrationBoundaries({});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD1STAR, FF_DSSD1STAR_HC})});

        addProcessSignature(-PID::DSSD1STAR, {-PID::DSTAR, PID::PI0});
        addIntegrationBoundaries({});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD1STAR, FF_DSSD1STAR_HC})});

        addProcessSignature(PID::DSSD1STARMINUS, {-PID::DSTAR, PID::PIMINUS});
        addIntegrationBoundaries({});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD1STAR, FF_DSSD1STAR_HC})});

        addProcessSignature(PID::DSSD1STARMINUS, {PID::DSTARMINUS, PID::PI0});
        addIntegrationBoundaries({});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD1STAR, FF_DSSD1STAR_HC})});

        // strange
        addProcessSignature(PID::DSSDS1STARMINUS, {PID::DSSTARMINUS, PID::PI0});
        addIntegrationBoundaries({});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSDS1STAR, FF_DSSDS1STAR_HC})});

        addProcessSignature(PID::DSSDS1STARMINUS, {PID::DSTARMINUS, -PID::K0});
        addIntegrationBoundaries({});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSDS1STAR, FF_DSSDS1STAR_HC})});

        addProcessSignature(PID::DSSDS1STARMINUS, {-PID::DSTAR, PID::KMINUS});
        addIntegrationBoundaries({});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSDS1STAR, FF_DSSDS1STAR_HC})});

        setSignatureIndex();
    }

    Tensor RateD1starDstarPi::evalAtPSPoint(const vector<double>& /*point*/) {
        auto labs = getTensor().labels();
        auto dimensions = getTensor().dims();
        Tensor result{"RateD1starDstarPi", MD::makeEmptySparse(dimensions, labs)};

        const double Mdss = masses()[0];
        const double Mdss2 = Mdss * Mdss;
        const double Mds = masses()[1];
        const double Mds2 = Mds * Mds;
        const double Mp = masses()[2];
        const double Mp2 = Mp * Mp;

        const double Ep = (Mdss2 - Mds2 + Mp2) / (2 * Mdss);
        const double Pp = sqrt(Ep * Ep - Mp2);

        double RateNorm = Pp / (24. * pi * Mdss2 * FPion * FPion);

        // set non-zero tensor elements, basis S,...
        result.element({0, 0}) = 3. * Mdss2 * (Mdss - Mds) * (Mdss - Mds); // Include prefactor ~ (Mdss Epi)^2

        result *= RateNorm;

        return result;
    }

} // namespace Hammer
