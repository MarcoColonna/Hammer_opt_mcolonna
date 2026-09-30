///
/// @file  RateD2starDPi.cc
/// @brief \f$ D_2^* \rightarrow D \pi \f$ total rate
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/Rates/RateD2starDPi.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"
#include "Hammer/Tools/Pdg.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    RateD2starDPi::RateD2starDPi() {
        // Create tensor rank and dimensions
        IndexList dims{{1, 1}};
        string name{"RateD2starDPi"};
        // auto& pdg = PID::instance();

        addProcessSignature(-PID::DSSD2STAR, {PID::DMINUS, PID::PIPLUS});
        addIntegrationBoundaries({});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD2STAR, FF_DSSD2STAR_HC})});

        addProcessSignature(-PID::DSSD2STAR, {-PID::D0, PID::PI0});
        addIntegrationBoundaries({});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD2STAR, FF_DSSD2STAR_HC})});

        addProcessSignature(PID::DSSD2STARMINUS, {-PID::D0, PID::PIMINUS});
        addIntegrationBoundaries({});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD2STAR, FF_DSSD2STAR_HC})});

        addProcessSignature(PID::DSSD2STARMINUS, {PID::DMINUS, PID::PI0});
        addIntegrationBoundaries({});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSD2STAR, FF_DSSD2STAR_HC})});

        // strange
        addProcessSignature(PID::DSSDS2STARMINUS, {PID::DSMINUS, PID::PI0});
        addIntegrationBoundaries({});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSDS2STAR, FF_DSSDS2STAR_HC})});

        addProcessSignature(PID::DSSDS2STARMINUS, {PID::DMINUS, -PID::K0});
        addIntegrationBoundaries({});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSDS2STAR, FF_DSSDS2STAR_HC})});

        addProcessSignature(PID::DSSDS2STARMINUS, {-PID::D0, PID::KMINUS});
        addIntegrationBoundaries({});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_DSSDS2STAR, FF_DSSDS2STAR_HC})});

        setSignatureIndex();
    }

    Tensor RateD2starDPi::evalAtPSPoint(const vector<double>& /*point*/) {
        auto labs = getTensor().labels();
        auto dimensions = getTensor().dims();
        Tensor result{"RateD2starDPi", MD::makeEmptySparse(dimensions, labs)};

        const double Mdss = masses()[0];
        const double Mdss2 = Mdss * Mdss;
        const double Md = masses()[1];
        const double Md2 = Md * Md;
        const double Mp = masses()[2];
        const double Mp2 = Mp * Mp;

        const double Ep = (Mdss2 - Md2 + Mp2) / (2 * Mdss);
        const double Pp = sqrt(Ep * Ep - Mp2);

        double RateNorm = Pp / (40. * pi * Mdss2 * FPion * FPion);

        // set non-zero tensor elements, basis S,D
        result.element({0, 0}) = 2. / 3. * pow(Pp, 4.);

        result *= RateNorm;

        return result;
    }

} // namespace Hammer
