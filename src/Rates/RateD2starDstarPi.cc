///
/// @file  RateD2starDstarPi.cc
/// @brief \f$ D_2^* \rightarrow D^* \pi \f$ total rate
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/Rates/RateD2starDstarPi.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"
#include "Hammer/Tools/Pdg.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    RateD2starDstarPi::RateD2starDstarPi() {
        // Create tensor rank and dimensions
        IndexList dims{{1, 1}};
        string name{"RateD2starDstarPi"};
        // auto& pdg = PID::instance();

        addProcessSignature(-PID::DSSD2STAR, {PID::DSTARMINUS, PID::PIPLUS});
        addIntegrationBoundaries({});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD2STAR, FF_DSSD2STAR_HC})});

        addProcessSignature(-PID::DSSD2STAR, {-PID::DSTAR, PID::PI0});
        addIntegrationBoundaries({});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD2STAR, FF_DSSD2STAR_HC})});

        addProcessSignature(PID::DSSD2STARMINUS, {-PID::DSTAR, PID::PIMINUS});
        addIntegrationBoundaries({});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD2STAR, FF_DSSD2STAR_HC})});

        addProcessSignature(PID::DSSD2STARMINUS, {PID::DSTARMINUS, PID::PI0});
        addIntegrationBoundaries({});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD2STAR, FF_DSSD2STAR_HC})});

        // strange
        addProcessSignature(PID::DSSDS2STARMINUS, {PID::DSSTARMINUS, PID::PI0});
        addIntegrationBoundaries({});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSDS2STAR, FF_DSSDS2STAR_HC})});

        addProcessSignature(PID::DSSDS2STARMINUS, {PID::DSTARMINUS, -PID::K0});
        addIntegrationBoundaries({});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSDS2STAR, FF_DSSDS2STAR_HC})});

        addProcessSignature(PID::DSSDS2STARMINUS, {-PID::DSTAR, PID::KMINUS});
        addIntegrationBoundaries({});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSDS2STAR, FF_DSSDS2STAR_HC})});

        setSignatureIndex();
    }

    Tensor RateD2starDstarPi::evalAtPSPoint(const vector<double>& /*point*/) {
        auto labs = getTensor().labels();
        auto dimensions = getTensor().dims();
        Tensor result{"RateD2starDstarPi", MD::makeEmptySparse(dimensions, labs)};

        const double Mdss = masses()[0];
        const double Mdss2 = Mdss * Mdss;
        const double Mds = masses()[1];
        const double Mds2 = Mds * Mds;
        const double Mp = masses()[2];
        const double Mp2 = Mp * Mp;

        const double Ep = (Mdss2 - Mds2 + Mp2) / (2 * Mdss);
        const double Pp = sqrt(Ep * Ep - Mp2);

        double RateNorm = Pp / (40. * pi * Mdss2 * FPion * FPion);

        // set non-zero tensor elements, basis S,D
        result.element({0, 0}) = pow(Pp, 4.);

        result *= RateNorm;

        return result;
    }

} // namespace Hammer
