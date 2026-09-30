///
/// @file  RateBD0starLepNu.cc
/// @brief \f$ B \rightarrow D_0^* \tau\nu \f$ total rate
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/Rates/RateBD0starLepNu.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"
#include "Hammer/Tools/Pdg.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    RateBD0starLepNu::RateBD0starLepNu() {
        // Create tensor rank and dimensions
        IndexList dims{{11, 4, 11, 4, _nPoints}};
        string name{"RateBD0starLepNuQ2"};
        auto& pdg = PID::instance();
        addProcessSignature(PID::BPLUS, {-PID::DSSD0STAR, PID::NU_TAU, PID::ANTITAU});
        addIntegrationBoundaries(
            {PS::makeQ2Function(pdg.getMass(PID::ANTITAU), pdg.getMass(PID::BPLUS) - pdg.getMass(PID::DSSD0STAR))});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCTAUNU, FF_BDSSD0STAR, WILSON_BCTAUNU_HC,
                                                          FF_BDSSD0STAR_HC, INTEGRATION_INDEX})});

        addProcessSignature(PID::BZERO, {PID::DSSD0STARMINUS, PID::NU_TAU, PID::ANTITAU});
        addIntegrationBoundaries({PS::makeQ2Function(pdg.getMass(PID::ANTITAU),
                                                     pdg.getMass(PID::BZERO) - pdg.getMass(PID::DSSD0STARMINUS))});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCTAUNU, FF_BDSSD0STAR, WILSON_BCTAUNU_HC,
                                                          FF_BDSSD0STAR_HC, INTEGRATION_INDEX})});

        addProcessSignature(PID::BPLUS, {-PID::DSSD0STAR, PID::NU_MU, PID::ANTIMUON});
        addIntegrationBoundaries(
            {PS::makeQ2Function(pdg.getMass(PID::ANTIMUON), pdg.getMass(PID::BPLUS) - pdg.getMass(PID::DSSD0STAR))});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCMUNU, FF_BDSSD0STAR, WILSON_BCMUNU_HC,
                                                          FF_BDSSD0STAR_HC, INTEGRATION_INDEX})});

        addProcessSignature(PID::BZERO, {PID::DSSD0STARMINUS, PID::NU_MU, PID::ANTIMUON});
        addIntegrationBoundaries({PS::makeQ2Function(pdg.getMass(PID::ANTIMUON),
                                                     pdg.getMass(PID::BZERO) - pdg.getMass(PID::DSSD0STARMINUS))});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCMUNU, FF_BDSSD0STAR, WILSON_BCMUNU_HC,
                                                          FF_BDSSD0STAR_HC, INTEGRATION_INDEX})});

        addProcessSignature(PID::BPLUS, {-PID::DSSD0STAR, PID::NU_E, PID::POSITRON});
        addIntegrationBoundaries(
            {PS::makeQ2Function(pdg.getMass(PID::POSITRON), pdg.getMass(PID::BPLUS) - pdg.getMass(PID::DSSD0STAR))});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCENU, FF_BDSSD0STAR, WILSON_BCENU_HC,
                                                          FF_BDSSD0STAR_HC, INTEGRATION_INDEX})});

        addProcessSignature(PID::BZERO, {PID::DSSD0STARMINUS, PID::NU_E, PID::POSITRON});
        addIntegrationBoundaries({PS::makeQ2Function(pdg.getMass(PID::POSITRON),
                                                     pdg.getMass(PID::BZERO) - pdg.getMass(PID::DSSD0STARMINUS))});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCENU, FF_BDSSD0STAR, WILSON_BCENU_HC,
                                                          FF_BDSSD0STAR_HC, INTEGRATION_INDEX})});

        // bs -> cs
        addProcessSignature(PID::BS, {PID::DSSDS0STARMINUS, PID::NU_TAU, PID::ANTITAU});
        addIntegrationBoundaries(
            {PS::makeQ2Function(pdg.getMass(PID::ANTITAU), pdg.getMass(PID::BS) - pdg.getMass(PID::DSSDS0STARMINUS))});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCTAUNU, FF_BSDSSDS0STAR, WILSON_BCTAUNU_HC,
                                                          FF_BSDSSDS0STAR_HC, INTEGRATION_INDEX})});

        addProcessSignature(PID::BS, {PID::DSSDS0STARMINUS, PID::NU_MU, PID::ANTIMUON});
        addIntegrationBoundaries(
            {PS::makeQ2Function(pdg.getMass(PID::ANTIMUON), pdg.getMass(PID::BS) - pdg.getMass(PID::DSSDS0STARMINUS))});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCMUNU, FF_BSDSSDS0STAR, WILSON_BCMUNU_HC,
                                                          FF_BSDSSDS0STAR_HC, INTEGRATION_INDEX})});

        addProcessSignature(PID::BS, {PID::DSSDS0STARMINUS, PID::NU_E, PID::POSITRON});
        addIntegrationBoundaries(
            {PS::makeQ2Function(pdg.getMass(PID::POSITRON), pdg.getMass(PID::BS) - pdg.getMass(PID::DSSDS0STARMINUS))});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCENU, FF_BSDSSDS0STAR, WILSON_BCENU_HC,
                                                          FF_BSDSSDS0STAR_HC, INTEGRATION_INDEX})});

        setSignatureIndex();
    }

    Tensor RateBD0starLepNu::evalAtPSPoint(const vector<double>& point) {
        auto labs = getTensor().labels();
        labs.pop_back();
        auto dimensions = getTensor().dims();
        dimensions.pop_back();
        Tensor result{"RateBD0starLepNu", MD::makeEmptySparse(dimensions, labs)};

        const double Mb = masses()[0];
        const double Mc = masses()[1];
        const double Mt = masses()[3];
        const double Mb2 = Mb * Mb;
        const double Mc2 = Mc * Mc;
        const double Mt2 = Mt * Mt;

        // kinematic objects
        const double Sqq = point[0];
        // const double Sqq2 = Sqq*Sqq;
        const double mSqq = Sqq / Mb2;
        const double w = (Mb2 + Mc2 - Sqq) / (2 * Mb * Mc);
        const double rC = Mc / Mb;
        const double rt = Mt / Mb;

        const double mSqqSq = mSqq * mSqq;
        const double wSq = w * w;
        const double rtSq = rt * rt;
        const double rCSq = rC * rC;
        const double rCp1Sq = (rC + 1.) * (rC + 1.);
        const double rCm1Sq = (rC - 1.) * (rC - 1.);
        // const double wp1Sq = (w + 1.)*(w + 1.);
        // const double wm1Sq = (w - 1.)*(w - 1.);
        // const double wSqm1Sq = (wSq - 1.)*(wSq - 1.);
        // const double wp1Cu = wp1Sq*(w + 1.);
        // const double wm1Cu = wm1Sq*(w - 1.);
        // const double wSqm1Cu = wSqm1Sq*(wSq - 1.);

        double RateNorm =
            (GFermi * GFermi * Mb * rCSq * (-Mt2 + Sqq) * (1 + rCSq - rtSq - 2 * rC * w) * sqrt(w * w - 1.)) /
            (32. * pow(pi, 3.) * (1 + rCSq - 2 * rC * w));

        // set non-zero tensor elements
        result.element({0, 1, 0, 1}) =
            ((1 + w) * (mSqq * rCp1Sq * (-1 + w) + rtSq * (1 + 2 * w - 2 * rC * (2 + w) + rCSq * (1 + 2 * w)))) /
            (3. * mSqqSq);
        result.element({0, 1, 0, 2}) = ((-1 + rCSq) * (mSqq + 2 * rtSq) * (-1 + wSq)) / (3. * mSqqSq);
        result.element({0, 1, 1, 0}) = ((-1 + rC) * rt * (1 + w)) / (2. * mSqq);
        result.element({0, 1, 3, 0}) = -((-1 + rC) * rt * (1 + w)) / (2. * mSqq);
        result.element({0, 1, 5, 1}) =
            -((1 + w) * (mSqq * rCp1Sq * (-1 + w) + rtSq * (1 + 2 * w - 2 * rC * (2 + w) + rCSq * (1 + 2 * w)))) /
            (3. * mSqqSq);
        result.element({0, 1, 5, 2}) = -((-1 + rCSq) * (mSqq + 2 * rtSq) * (-1 + wSq)) / (3. * mSqqSq);
        result.element({0, 1, 7, 1}) =
            ((1 + w) * (mSqq * rCp1Sq * (-1 + w) + rtSq * (1 + 2 * w - 2 * rC * (2 + w) + rCSq * (1 + 2 * w)))) /
            (3. * mSqqSq);
        result.element({0, 1, 7, 2}) = ((-1 + rCSq) * (mSqq + 2 * rtSq) * (-1 + wSq)) / (3. * mSqqSq);
        result.element({0, 1, 9, 3}) = (-2 * (1 + rC) * rt * (-1 + wSq)) / mSqq;
        result.element({0, 2, 0, 1}) = ((-1 + rCSq) * (mSqq + 2 * rtSq) * (-1 + wSq)) / (3. * mSqqSq);
        result.element({0, 2, 0, 2}) =
            ((-1 + w) * (mSqq * rCm1Sq * (1 + w) + rtSq * (-1 + 2 * rC * (-2 + w) + 2 * w + rCSq * (-1 + 2 * w)))) /
            (3. * mSqqSq);
        result.element({0, 2, 1, 0}) = ((1 + rC) * rt * (-1 + w)) / (2. * mSqq);
        result.element({0, 2, 3, 0}) = -((1 + rC) * rt * (-1 + w)) / (2. * mSqq);
        result.element({0, 2, 5, 1}) = -((-1 + rCSq) * (mSqq + 2 * rtSq) * (-1 + wSq)) / (3. * mSqqSq);
        result.element({0, 2, 5, 2}) =
            -((-1 + w) * (mSqq * rCm1Sq * (1 + w) + rtSq * (-1 + 2 * rC * (-2 + w) + 2 * w + rCSq * (-1 + 2 * w)))) /
            (3. * mSqqSq);
        result.element({0, 2, 7, 1}) = ((-1 + rCSq) * (mSqq + 2 * rtSq) * (-1 + wSq)) / (3. * mSqqSq);
        result.element({0, 2, 7, 2}) =
            ((-1 + w) * (mSqq * rCm1Sq * (1 + w) + rtSq * (-1 + 2 * rC * (-2 + w) + 2 * w + rCSq * (-1 + 2 * w)))) /
            (3. * mSqqSq);
        result.element({0, 2, 9, 3}) = (-2 * (-1 + rC) * rt * (-1 + wSq)) / mSqq;
        result.element({1, 0, 0, 1}) = ((-1 + rC) * rt * (1 + w)) / (2. * mSqq);
        result.element({1, 0, 0, 2}) = ((1 + rC) * rt * (-1 + w)) / (2. * mSqq);
        result.element({1, 0, 1, 0}) = 0.5;
        result.element({1, 0, 3, 0}) = -0.5;
        result.element({1, 0, 5, 1}) = -((-1 + rC) * rt * (1 + w)) / (2. * mSqq);
        result.element({1, 0, 5, 2}) = -((1 + rC) * rt * (-1 + w)) / (2. * mSqq);
        result.element({1, 0, 7, 1}) = ((-1 + rC) * rt * (1 + w)) / (2. * mSqq);
        result.element({1, 0, 7, 2}) = ((1 + rC) * rt * (-1 + w)) / (2. * mSqq);
        result.element({2, 0, 2, 0}) = 0.5;
        result.element({2, 0, 4, 0}) = -0.5;
        result.element({2, 0, 6, 1}) = -((-1 + rC) * rt * (1 + w)) / (2. * mSqq);
        result.element({2, 0, 6, 2}) = -((1 + rC) * rt * (-1 + w)) / (2. * mSqq);
        result.element({2, 0, 8, 1}) = ((-1 + rC) * rt * (1 + w)) / (2. * mSqq);
        result.element({2, 0, 8, 2}) = ((1 + rC) * rt * (-1 + w)) / (2. * mSqq);
        result.element({3, 0, 0, 1}) = -((-1 + rC) * rt * (1 + w)) / (2. * mSqq);
        result.element({3, 0, 0, 2}) = -((1 + rC) * rt * (-1 + w)) / (2. * mSqq);
        result.element({3, 0, 1, 0}) = -0.5;
        result.element({3, 0, 3, 0}) = 0.5;
        result.element({3, 0, 5, 1}) = ((-1 + rC) * rt * (1 + w)) / (2. * mSqq);
        result.element({3, 0, 5, 2}) = ((1 + rC) * rt * (-1 + w)) / (2. * mSqq);
        result.element({3, 0, 7, 1}) = -((-1 + rC) * rt * (1 + w)) / (2. * mSqq);
        result.element({3, 0, 7, 2}) = -((1 + rC) * rt * (-1 + w)) / (2. * mSqq);
        result.element({4, 0, 2, 0}) = -0.5;
        result.element({4, 0, 4, 0}) = 0.5;
        result.element({4, 0, 6, 1}) = ((-1 + rC) * rt * (1 + w)) / (2. * mSqq);
        result.element({4, 0, 6, 2}) = ((1 + rC) * rt * (-1 + w)) / (2. * mSqq);
        result.element({4, 0, 8, 1}) = -((-1 + rC) * rt * (1 + w)) / (2. * mSqq);
        result.element({4, 0, 8, 2}) = -((1 + rC) * rt * (-1 + w)) / (2. * mSqq);
        result.element({5, 1, 0, 1}) =
            -((1 + w) * (mSqq * rCp1Sq * (-1 + w) + rtSq * (1 + 2 * w - 2 * rC * (2 + w) + rCSq * (1 + 2 * w)))) /
            (3. * mSqqSq);
        result.element({5, 1, 0, 2}) = -((-1 + rCSq) * (mSqq + 2 * rtSq) * (-1 + wSq)) / (3. * mSqqSq);
        result.element({5, 1, 1, 0}) = -((-1 + rC) * rt * (1 + w)) / (2. * mSqq);
        result.element({5, 1, 3, 0}) = ((-1 + rC) * rt * (1 + w)) / (2. * mSqq);
        result.element({5, 1, 5, 1}) =
            ((1 + w) * (mSqq * rCp1Sq * (-1 + w) + rtSq * (1 + 2 * w - 2 * rC * (2 + w) + rCSq * (1 + 2 * w)))) /
            (3. * mSqqSq);
        result.element({5, 1, 5, 2}) = ((-1 + rCSq) * (mSqq + 2 * rtSq) * (-1 + wSq)) / (3. * mSqqSq);
        result.element({5, 1, 7, 1}) =
            -((1 + w) * (mSqq * rCp1Sq * (-1 + w) + rtSq * (1 + 2 * w - 2 * rC * (2 + w) + rCSq * (1 + 2 * w)))) /
            (3. * mSqqSq);
        result.element({5, 1, 7, 2}) = -((-1 + rCSq) * (mSqq + 2 * rtSq) * (-1 + wSq)) / (3. * mSqqSq);
        result.element({5, 1, 9, 3}) = (2 * (1 + rC) * rt * (-1 + wSq)) / mSqq;
        result.element({5, 2, 0, 1}) = -((-1 + rCSq) * (mSqq + 2 * rtSq) * (-1 + wSq)) / (3. * mSqqSq);
        result.element({5, 2, 0, 2}) =
            -((-1 + w) * (mSqq * rCm1Sq * (1 + w) + rtSq * (-1 + 2 * rC * (-2 + w) + 2 * w + rCSq * (-1 + 2 * w)))) /
            (3. * mSqqSq);
        result.element({5, 2, 1, 0}) = -((1 + rC) * rt * (-1 + w)) / (2. * mSqq);
        result.element({5, 2, 3, 0}) = ((1 + rC) * rt * (-1 + w)) / (2. * mSqq);
        result.element({5, 2, 5, 1}) = ((-1 + rCSq) * (mSqq + 2 * rtSq) * (-1 + wSq)) / (3. * mSqqSq);
        result.element({5, 2, 5, 2}) =
            ((-1 + w) * (mSqq * rCm1Sq * (1 + w) + rtSq * (-1 + 2 * rC * (-2 + w) + 2 * w + rCSq * (-1 + 2 * w)))) /
            (3. * mSqqSq);
        result.element({5, 2, 7, 1}) = -((-1 + rCSq) * (mSqq + 2 * rtSq) * (-1 + wSq)) / (3. * mSqqSq);
        result.element({5, 2, 7, 2}) =
            -((-1 + w) * (mSqq * rCm1Sq * (1 + w) + rtSq * (-1 + 2 * rC * (-2 + w) + 2 * w + rCSq * (-1 + 2 * w)))) /
            (3. * mSqqSq);
        result.element({5, 2, 9, 3}) = (2 * (-1 + rC) * rt * (-1 + wSq)) / mSqq;
        result.element({6, 1, 2, 0}) = -((-1 + rC) * rt * (1 + w)) / (2. * mSqq);
        result.element({6, 1, 4, 0}) = ((-1 + rC) * rt * (1 + w)) / (2. * mSqq);
        result.element({6, 1, 6, 1}) =
            ((1 + w) * (mSqq * rCp1Sq * (-1 + w) + rtSq * (1 + 2 * w - 2 * rC * (2 + w) + rCSq * (1 + 2 * w)))) /
            (3. * mSqqSq);
        result.element({6, 1, 6, 2}) = ((-1 + rCSq) * (mSqq + 2 * rtSq) * (-1 + wSq)) / (3. * mSqqSq);
        result.element({6, 1, 8, 1}) =
            -((1 + w) * (mSqq * rCp1Sq * (-1 + w) + rtSq * (1 + 2 * w - 2 * rC * (2 + w) + rCSq * (1 + 2 * w)))) /
            (3. * mSqqSq);
        result.element({6, 1, 8, 2}) = -((-1 + rCSq) * (mSqq + 2 * rtSq) * (-1 + wSq)) / (3. * mSqqSq);
        result.element({6, 1, 10, 3}) = (-2 * (1 + rC) * rt * (-1 + wSq)) / mSqq;
        result.element({6, 2, 2, 0}) = -((1 + rC) * rt * (-1 + w)) / (2. * mSqq);
        result.element({6, 2, 4, 0}) = ((1 + rC) * rt * (-1 + w)) / (2. * mSqq);
        result.element({6, 2, 6, 1}) = ((-1 + rCSq) * (mSqq + 2 * rtSq) * (-1 + wSq)) / (3. * mSqqSq);
        result.element({6, 2, 6, 2}) =
            ((-1 + w) * (mSqq * rCm1Sq * (1 + w) + rtSq * (-1 + 2 * rC * (-2 + w) + 2 * w + rCSq * (-1 + 2 * w)))) /
            (3. * mSqqSq);
        result.element({6, 2, 8, 1}) = -((-1 + rCSq) * (mSqq + 2 * rtSq) * (-1 + wSq)) / (3. * mSqqSq);
        result.element({6, 2, 8, 2}) =
            -((-1 + w) * (mSqq * rCm1Sq * (1 + w) + rtSq * (-1 + 2 * rC * (-2 + w) + 2 * w + rCSq * (-1 + 2 * w)))) /
            (3. * mSqqSq);
        result.element({6, 2, 10, 3}) = (-2 * (-1 + rC) * rt * (-1 + wSq)) / mSqq;
        result.element({7, 1, 0, 1}) =
            ((1 + w) * (mSqq * rCp1Sq * (-1 + w) + rtSq * (1 + 2 * w - 2 * rC * (2 + w) + rCSq * (1 + 2 * w)))) /
            (3. * mSqqSq);
        result.element({7, 1, 0, 2}) = ((-1 + rCSq) * (mSqq + 2 * rtSq) * (-1 + wSq)) / (3. * mSqqSq);
        result.element({7, 1, 1, 0}) = ((-1 + rC) * rt * (1 + w)) / (2. * mSqq);
        result.element({7, 1, 3, 0}) = -((-1 + rC) * rt * (1 + w)) / (2. * mSqq);
        result.element({7, 1, 5, 1}) =
            -((1 + w) * (mSqq * rCp1Sq * (-1 + w) + rtSq * (1 + 2 * w - 2 * rC * (2 + w) + rCSq * (1 + 2 * w)))) /
            (3. * mSqqSq);
        result.element({7, 1, 5, 2}) = -((-1 + rCSq) * (mSqq + 2 * rtSq) * (-1 + wSq)) / (3. * mSqqSq);
        result.element({7, 1, 7, 1}) =
            ((1 + w) * (mSqq * rCp1Sq * (-1 + w) + rtSq * (1 + 2 * w - 2 * rC * (2 + w) + rCSq * (1 + 2 * w)))) /
            (3. * mSqqSq);
        result.element({7, 1, 7, 2}) = ((-1 + rCSq) * (mSqq + 2 * rtSq) * (-1 + wSq)) / (3. * mSqqSq);
        result.element({7, 1, 9, 3}) = (-2 * (1 + rC) * rt * (-1 + wSq)) / mSqq;
        result.element({7, 2, 0, 1}) = ((-1 + rCSq) * (mSqq + 2 * rtSq) * (-1 + wSq)) / (3. * mSqqSq);
        result.element({7, 2, 0, 2}) =
            ((-1 + w) * (mSqq * rCm1Sq * (1 + w) + rtSq * (-1 + 2 * rC * (-2 + w) + 2 * w + rCSq * (-1 + 2 * w)))) /
            (3. * mSqqSq);
        result.element({7, 2, 1, 0}) = ((1 + rC) * rt * (-1 + w)) / (2. * mSqq);
        result.element({7, 2, 3, 0}) = -((1 + rC) * rt * (-1 + w)) / (2. * mSqq);
        result.element({7, 2, 5, 1}) = -((-1 + rCSq) * (mSqq + 2 * rtSq) * (-1 + wSq)) / (3. * mSqqSq);
        result.element({7, 2, 5, 2}) =
            -((-1 + w) * (mSqq * rCm1Sq * (1 + w) + rtSq * (-1 + 2 * rC * (-2 + w) + 2 * w + rCSq * (-1 + 2 * w)))) /
            (3. * mSqqSq);
        result.element({7, 2, 7, 1}) = ((-1 + rCSq) * (mSqq + 2 * rtSq) * (-1 + wSq)) / (3. * mSqqSq);
        result.element({7, 2, 7, 2}) =
            ((-1 + w) * (mSqq * rCm1Sq * (1 + w) + rtSq * (-1 + 2 * rC * (-2 + w) + 2 * w + rCSq * (-1 + 2 * w)))) /
            (3. * mSqqSq);
        result.element({7, 2, 9, 3}) = (-2 * (-1 + rC) * rt * (-1 + wSq)) / mSqq;
        result.element({8, 1, 2, 0}) = ((-1 + rC) * rt * (1 + w)) / (2. * mSqq);
        result.element({8, 1, 4, 0}) = -((-1 + rC) * rt * (1 + w)) / (2. * mSqq);
        result.element({8, 1, 6, 1}) =
            -((1 + w) * (mSqq * rCp1Sq * (-1 + w) + rtSq * (1 + 2 * w - 2 * rC * (2 + w) + rCSq * (1 + 2 * w)))) /
            (3. * mSqqSq);
        result.element({8, 1, 6, 2}) = -((-1 + rCSq) * (mSqq + 2 * rtSq) * (-1 + wSq)) / (3. * mSqqSq);
        result.element({8, 1, 8, 1}) =
            ((1 + w) * (mSqq * rCp1Sq * (-1 + w) + rtSq * (1 + 2 * w - 2 * rC * (2 + w) + rCSq * (1 + 2 * w)))) /
            (3. * mSqqSq);
        result.element({8, 1, 8, 2}) = ((-1 + rCSq) * (mSqq + 2 * rtSq) * (-1 + wSq)) / (3. * mSqqSq);
        result.element({8, 1, 10, 3}) = (2 * (1 + rC) * rt * (-1 + wSq)) / mSqq;
        result.element({8, 2, 2, 0}) = ((1 + rC) * rt * (-1 + w)) / (2. * mSqq);
        result.element({8, 2, 4, 0}) = -((1 + rC) * rt * (-1 + w)) / (2. * mSqq);
        result.element({8, 2, 6, 1}) = -((-1 + rCSq) * (mSqq + 2 * rtSq) * (-1 + wSq)) / (3. * mSqqSq);
        result.element({8, 2, 6, 2}) =
            -((-1 + w) * (mSqq * rCm1Sq * (1 + w) + rtSq * (-1 + 2 * rC * (-2 + w) + 2 * w + rCSq * (-1 + 2 * w)))) /
            (3. * mSqqSq);
        result.element({8, 2, 8, 1}) = ((-1 + rCSq) * (mSqq + 2 * rtSq) * (-1 + wSq)) / (3. * mSqqSq);
        result.element({8, 2, 8, 2}) =
            ((-1 + w) * (mSqq * rCm1Sq * (1 + w) + rtSq * (-1 + 2 * rC * (-2 + w) + 2 * w + rCSq * (-1 + 2 * w)))) /
            (3. * mSqqSq);
        result.element({8, 2, 10, 3}) = (2 * (-1 + rC) * rt * (-1 + wSq)) / mSqq;
        result.element({9, 3, 0, 1}) = (-2 * (1 + rC) * rt * (-1 + wSq)) / mSqq;
        result.element({9, 3, 0, 2}) = (-2 * (-1 + rC) * rt * (-1 + wSq)) / mSqq;
        result.element({9, 3, 5, 1}) = (2 * (1 + rC) * rt * (-1 + wSq)) / mSqq;
        result.element({9, 3, 5, 2}) = (2 * (-1 + rC) * rt * (-1 + wSq)) / mSqq;
        result.element({9, 3, 7, 1}) = (-2 * (1 + rC) * rt * (-1 + wSq)) / mSqq;
        result.element({9, 3, 7, 2}) = (-2 * (-1 + rC) * rt * (-1 + wSq)) / mSqq;
        result.element({9, 3, 9, 3}) = (8 * (mSqq + 2 * rtSq) * (-1 + wSq)) / (3. * mSqq);
        result.element({10, 3, 6, 1}) = (-2 * (1 + rC) * rt * (-1 + wSq)) / mSqq;
        result.element({10, 3, 6, 2}) = (-2 * (-1 + rC) * rt * (-1 + wSq)) / mSqq;
        result.element({10, 3, 8, 1}) = (2 * (1 + rC) * rt * (-1 + wSq)) / mSqq;
        result.element({10, 3, 8, 2}) = (2 * (-1 + rC) * rt * (-1 + wSq)) / mSqq;
        result.element({10, 3, 10, 3}) = (8 * (mSqq + 2 * rtSq) * (-1 + wSq)) / (3. * mSqq);

        result *= RateNorm;

        return result;
    }

} // namespace Hammer
