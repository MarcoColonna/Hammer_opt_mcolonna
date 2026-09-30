///
/// @file  RateBD1starLepNu.cc
/// @brief \f$ B \rightarrow D_1^* \tau\nu \f$ total rate
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/Rates/RateBD1starLepNu.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"
#include "Hammer/Tools/Pdg.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    RateBD1starLepNu::RateBD1starLepNu() {
        // Create tensor rank and dimensions
        IndexList dims{{11, 8, 11, 8, _nPoints}};
        string name{"RateBD1starLepNuQ2"};
        auto& pdg = PID::instance();
        addProcessSignature(PID::BPLUS, {-PID::DSSD1STAR, PID::NU_TAU, PID::ANTITAU});
        addIntegrationBoundaries(
            {PS::makeQ2Function(pdg.getMass(PID::ANTITAU), pdg.getMass(PID::BPLUS) - pdg.getMass(PID::DSSD1STAR))});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCTAUNU, FF_BDSSD1STAR, WILSON_BCTAUNU_HC,
                                                          FF_BDSSD1STAR_HC, INTEGRATION_INDEX})});

        addProcessSignature(PID::BZERO, {PID::DSSD1STARMINUS, PID::NU_TAU, PID::ANTITAU});
        addIntegrationBoundaries({PS::makeQ2Function(pdg.getMass(PID::ANTITAU),
                                                     pdg.getMass(PID::BZERO) - pdg.getMass(PID::DSSD1STARMINUS))});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCTAUNU, FF_BDSSD1STAR, WILSON_BCTAUNU_HC,
                                                          FF_BDSSD1STAR_HC, INTEGRATION_INDEX})});

        addProcessSignature(PID::BPLUS, {-PID::DSSD1STAR, PID::NU_MU, PID::ANTIMUON});
        addIntegrationBoundaries(
            {PS::makeQ2Function(pdg.getMass(PID::ANTIMUON), pdg.getMass(PID::BPLUS) - pdg.getMass(PID::DSSD1STAR))});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCMUNU, FF_BDSSD1STAR, WILSON_BCMUNU_HC,
                                                          FF_BDSSD1STAR_HC, INTEGRATION_INDEX})});

        addProcessSignature(PID::BZERO, {PID::DSSD1STARMINUS, PID::NU_MU, PID::ANTIMUON});
        addIntegrationBoundaries({PS::makeQ2Function(pdg.getMass(PID::ANTIMUON),
                                                     pdg.getMass(PID::BZERO) - pdg.getMass(PID::DSSD1STARMINUS))});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCMUNU, FF_BDSSD1STAR, WILSON_BCMUNU_HC,
                                                          FF_BDSSD1STAR_HC, INTEGRATION_INDEX})});

        addProcessSignature(PID::BPLUS, {-PID::DSSD1STAR, PID::NU_E, PID::POSITRON});
        addIntegrationBoundaries(
            {PS::makeQ2Function(pdg.getMass(PID::POSITRON), pdg.getMass(PID::BPLUS) - pdg.getMass(PID::DSSD1STAR))});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCENU, FF_BDSSD1STAR, WILSON_BCENU_HC,
                                                          FF_BDSSD1STAR_HC, INTEGRATION_INDEX})});

        addProcessSignature(PID::BZERO, {PID::DSSD1STARMINUS, PID::NU_E, PID::POSITRON});
        addIntegrationBoundaries({PS::makeQ2Function(pdg.getMass(PID::POSITRON),
                                                     pdg.getMass(PID::BZERO) - pdg.getMass(PID::DSSD1STARMINUS))});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCENU, FF_BDSSD1STAR, WILSON_BCENU_HC,
                                                          FF_BDSSD1STAR_HC, INTEGRATION_INDEX})});

        // bs -> cs
        addProcessSignature(PID::BS, {PID::DSSDS1STARMINUS, PID::NU_TAU, PID::ANTITAU});
        addIntegrationBoundaries(
            {PS::makeQ2Function(pdg.getMass(PID::ANTITAU), pdg.getMass(PID::BS) - pdg.getMass(PID::DSSDS1STARMINUS))});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCTAUNU, FF_BSDSSDS1STAR, WILSON_BCTAUNU_HC,
                                                          FF_BSDSSDS1STAR_HC, INTEGRATION_INDEX})});

        addProcessSignature(PID::BS, {PID::DSSDS1STARMINUS, PID::NU_MU, PID::ANTIMUON});
        addIntegrationBoundaries(
            {PS::makeQ2Function(pdg.getMass(PID::ANTIMUON), pdg.getMass(PID::BS) - pdg.getMass(PID::DSSDS1STARMINUS))});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCMUNU, FF_BSDSSDS1STAR, WILSON_BCMUNU_HC,
                                                          FF_BSDSSDS1STAR_HC, INTEGRATION_INDEX})});

        addProcessSignature(PID::BS, {PID::DSSDS1STARMINUS, PID::NU_E, PID::POSITRON});
        addIntegrationBoundaries(
            {PS::makeQ2Function(pdg.getMass(PID::POSITRON), pdg.getMass(PID::BS) - pdg.getMass(PID::DSSDS1STARMINUS))});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCENU, FF_BSDSSDS1STAR, WILSON_BCENU_HC,
                                                          FF_BSDSSDS1STAR_HC, INTEGRATION_INDEX})});

        setSignatureIndex();
    }

    Tensor RateBD1starLepNu::evalAtPSPoint(const vector<double>& point) {
        auto labs = getTensor().labels();
        labs.pop_back();
        auto dimensions = getTensor().dims();
        dimensions.pop_back();
        Tensor result{"RateBD1starLepNu", MD::makeEmptySparse(dimensions, labs)};

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
        // const double rCp1Sq = (rC + 1.)*(rC + 1.);
        // const double rCm1Sq = (rC - 1.)*(rC - 1.);
        // const double wp1Sq = (w + 1.)*(w + 1.);
        // const double wm1Sq = (w - 1.)*(w - 1.);
        const double wSqm1Sq = (wSq - 1.) * (wSq - 1.);
        // const double wp1Cu = wp1Sq*(w + 1.);
        // const double wm1Cu = wm1Sq*(w - 1.);
        // const double wSqm1Cu = wSqm1Sq*(wSq - 1.);

        double RateNorm =
            (GFermi * GFermi * Mb * rCSq * (-Mt2 + Sqq) * (1 + rCSq - rtSq - 2 * rC * w) * sqrt(w * w - 1.)) /
            (64. * pow(pi, 3.) * (1 + rCSq - 2 * rC * w));

        // set non-zero tensor elements
        result.element({0, 1, 0, 1}) =
            (4 * mSqqSq + 2 * mSqq * (rCSq + rtSq - 2 * rC * w + wSq) + rtSq * (-3 + rCSq - 2 * rC * w + 4 * wSq)) /
            (3. * mSqqSq);
        result.element({0, 1, 0, 2}) =
            -((-1 + wSq) * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (3. * mSqqSq);
        result.element({0, 1, 0, 3}) = (-2 * (mSqq + 2 * rtSq) * (rC - w) * (-1 + wSq)) / (3. * mSqqSq);
        result.element({0, 1, 1, 0}) = (rt * (-1 + wSq)) / mSqq;
        result.element({0, 1, 3, 0}) = (rt * (-1 + wSq)) / mSqq;
        result.element({0, 1, 5, 1}) =
            (4 * mSqqSq + 2 * mSqq * (rCSq + rtSq - 2 * rC * w + wSq) + rtSq * (-3 + rCSq - 2 * rC * w + 4 * wSq)) /
            (3. * mSqqSq);
        result.element({0, 1, 5, 2}) =
            -((-1 + wSq) * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (3. * mSqqSq);
        result.element({0, 1, 5, 3}) = (-2 * (mSqq + 2 * rtSq) * (rC - w) * (-1 + wSq)) / (3. * mSqqSq);
        result.element({0, 1, 7, 1}) =
            (4 * mSqqSq + 2 * mSqq * (rCSq + rtSq - 2 * rC * w + wSq) + rtSq * (-3 + rCSq - 2 * rC * w + 4 * wSq)) /
            (3. * mSqqSq);
        result.element({0, 1, 7, 2}) =
            -((-1 + wSq) * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (3. * mSqqSq);
        result.element({0, 1, 7, 3}) = (-2 * (mSqq + 2 * rtSq) * (rC - w) * (-1 + wSq)) / (3. * mSqqSq);
        result.element({0, 1, 9, 5}) = (4 * rt * (-2 + 3 * rC * w - wSq)) / mSqq;
        result.element({0, 1, 9, 6}) = (12 * rt * (rC - w)) / mSqq;
        result.element({0, 1, 9, 7}) = (-4 * rt * (rC - w) * (-1 + wSq)) / mSqq;
        result.element({0, 2, 0, 1}) =
            -((-1 + wSq) * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (3. * mSqqSq);
        result.element({0, 2, 0, 2}) =
            ((-1 + wSq) * (2 * mSqq * rCSq * (-1 + wSq) + rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)))) /
            (3. * mSqqSq);
        result.element({0, 2, 0, 3}) =
            ((-1 + wSq) * (2 * mSqq * rC * (-1 + wSq) + rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)))) /
            (3. * mSqqSq);
        result.element({0, 2, 1, 0}) = -((rt * (-1 + rC * w) * (-1 + wSq)) / mSqq);
        result.element({0, 2, 3, 0}) = -((rt * (-1 + rC * w) * (-1 + wSq)) / mSqq);
        result.element({0, 2, 5, 1}) =
            -((-1 + wSq) * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (3. * mSqqSq);
        result.element({0, 2, 5, 2}) =
            ((-1 + wSq) * (2 * mSqq * rCSq * (-1 + wSq) + rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)))) /
            (3. * mSqqSq);
        result.element({0, 2, 5, 3}) =
            ((-1 + wSq) * (2 * mSqq * rC * (-1 + wSq) + rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)))) /
            (3. * mSqqSq);
        result.element({0, 2, 7, 1}) =
            -((-1 + wSq) * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (3. * mSqqSq);
        result.element({0, 2, 7, 2}) =
            ((-1 + wSq) * (2 * mSqq * rCSq * (-1 + wSq) + rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)))) /
            (3. * mSqqSq);
        result.element({0, 2, 7, 3}) =
            ((-1 + wSq) * (2 * mSqq * rC * (-1 + wSq) + rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)))) /
            (3. * mSqqSq);
        result.element({0, 2, 9, 5}) = (-4 * rC * rt * w * (-1 + wSq)) / mSqq;
        result.element({0, 2, 9, 6}) = (-4 * rC * rt * (-1 + wSq)) / mSqq;
        result.element({0, 2, 9, 7}) = (4 * rC * rt * wSqm1Sq) / mSqq;
        result.element({0, 3, 0, 1}) = (-2 * (mSqq + 2 * rtSq) * (rC - w) * (-1 + wSq)) / (3. * mSqqSq);
        result.element({0, 3, 0, 2}) =
            ((-1 + wSq) * (2 * mSqq * rC * (-1 + wSq) + rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)))) /
            (3. * mSqqSq);
        result.element({0, 3, 0, 3}) =
            ((-1 + wSq) * (2 * mSqq * (-1 + wSq) + rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq))) / (3. * mSqqSq);
        result.element({0, 3, 1, 0}) = (rt * (-rC + w) * (-1 + wSq)) / mSqq;
        result.element({0, 3, 3, 0}) = (rt * (-rC + w) * (-1 + wSq)) / mSqq;
        result.element({0, 3, 5, 1}) = (-2 * (mSqq + 2 * rtSq) * (rC - w) * (-1 + wSq)) / (3. * mSqqSq);
        result.element({0, 3, 5, 2}) =
            ((-1 + wSq) * (2 * mSqq * rC * (-1 + wSq) + rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)))) /
            (3. * mSqqSq);
        result.element({0, 3, 5, 3}) =
            ((-1 + wSq) * (2 * mSqq * (-1 + wSq) + rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq))) / (3. * mSqqSq);
        result.element({0, 3, 7, 1}) = (-2 * (mSqq + 2 * rtSq) * (rC - w) * (-1 + wSq)) / (3. * mSqqSq);
        result.element({0, 3, 7, 2}) =
            ((-1 + wSq) * (2 * mSqq * rC * (-1 + wSq) + rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)))) /
            (3. * mSqqSq);
        result.element({0, 3, 7, 3}) =
            ((-1 + wSq) * (2 * mSqq * (-1 + wSq) + rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq))) / (3. * mSqqSq);
        result.element({0, 3, 9, 5}) = (-4 * rt * w * (-1 + wSq)) / mSqq;
        result.element({0, 3, 9, 6}) = (-4 * rt * (-1 + wSq)) / mSqq;
        result.element({0, 3, 9, 7}) = (4 * rt * wSqm1Sq) / mSqq;
        result.element({0, 4, 0, 4}) = (2 * (2 * mSqq + rtSq) * (-1 + wSq)) / (3. * mSqq);
        result.element({0, 4, 5, 4}) = (-2 * (2 * mSqq + rtSq) * (-1 + wSq)) / (3. * mSqq);
        result.element({0, 4, 7, 4}) = (2 * (2 * mSqq + rtSq) * (-1 + wSq)) / (3. * mSqq);
        result.element({0, 4, 9, 5}) = (8 * rC * rt * (-1 + wSq)) / mSqq;
        result.element({0, 4, 9, 6}) = (8 * rt * (-1 + wSq)) / mSqq;
        result.element({1, 0, 0, 1}) = (rt * (-1 + wSq)) / mSqq;
        result.element({1, 0, 0, 2}) = -((rt * (-1 + rC * w) * (-1 + wSq)) / mSqq);
        result.element({1, 0, 0, 3}) = (rt * (-rC + w) * (-1 + wSq)) / mSqq;
        result.element({1, 0, 1, 0}) = -1 + wSq;
        result.element({1, 0, 3, 0}) = -1 + wSq;
        result.element({1, 0, 5, 1}) = (rt * (-1 + wSq)) / mSqq;
        result.element({1, 0, 5, 2}) = -((rt * (-1 + rC * w) * (-1 + wSq)) / mSqq);
        result.element({1, 0, 5, 3}) = (rt * (-rC + w) * (-1 + wSq)) / mSqq;
        result.element({1, 0, 7, 1}) = (rt * (-1 + wSq)) / mSqq;
        result.element({1, 0, 7, 2}) = -((rt * (-1 + rC * w) * (-1 + wSq)) / mSqq);
        result.element({1, 0, 7, 3}) = (rt * (-rC + w) * (-1 + wSq)) / mSqq;
        result.element({2, 0, 2, 0}) = -1 + wSq;
        result.element({2, 0, 4, 0}) = -1 + wSq;
        result.element({2, 0, 6, 1}) = (rt * (-1 + wSq)) / mSqq;
        result.element({2, 0, 6, 2}) = -((rt * (-1 + rC * w) * (-1 + wSq)) / mSqq);
        result.element({2, 0, 6, 3}) = -((rt * (rC - w) * (-1 + wSq)) / mSqq);
        result.element({2, 0, 8, 1}) = (rt * (-1 + wSq)) / mSqq;
        result.element({2, 0, 8, 2}) = -((rt * (-1 + rC * w) * (-1 + wSq)) / mSqq);
        result.element({2, 0, 8, 3}) = -((rt * (rC - w) * (-1 + wSq)) / mSqq);
        result.element({3, 0, 0, 1}) = (rt * (-1 + wSq)) / mSqq;
        result.element({3, 0, 0, 2}) = -((rt * (-1 + rC * w) * (-1 + wSq)) / mSqq);
        result.element({3, 0, 0, 3}) = (rt * (-rC + w) * (-1 + wSq)) / mSqq;
        result.element({3, 0, 1, 0}) = -1 + wSq;
        result.element({3, 0, 3, 0}) = -1 + wSq;
        result.element({3, 0, 5, 1}) = (rt * (-1 + wSq)) / mSqq;
        result.element({3, 0, 5, 2}) = -((rt * (-1 + rC * w) * (-1 + wSq)) / mSqq);
        result.element({3, 0, 5, 3}) = (rt * (-rC + w) * (-1 + wSq)) / mSqq;
        result.element({3, 0, 7, 1}) = (rt * (-1 + wSq)) / mSqq;
        result.element({3, 0, 7, 2}) = -((rt * (-1 + rC * w) * (-1 + wSq)) / mSqq);
        result.element({3, 0, 7, 3}) = (rt * (-rC + w) * (-1 + wSq)) / mSqq;
        result.element({4, 0, 2, 0}) = -1 + wSq;
        result.element({4, 0, 4, 0}) = -1 + wSq;
        result.element({4, 0, 6, 1}) = (rt * (-1 + wSq)) / mSqq;
        result.element({4, 0, 6, 2}) = -((rt * (-1 + rC * w) * (-1 + wSq)) / mSqq);
        result.element({4, 0, 6, 3}) = -((rt * (rC - w) * (-1 + wSq)) / mSqq);
        result.element({4, 0, 8, 1}) = (rt * (-1 + wSq)) / mSqq;
        result.element({4, 0, 8, 2}) = -((rt * (-1 + rC * w) * (-1 + wSq)) / mSqq);
        result.element({4, 0, 8, 3}) = -((rt * (rC - w) * (-1 + wSq)) / mSqq);
        result.element({5, 1, 0, 1}) =
            (4 * mSqqSq + 2 * mSqq * (rCSq + rtSq - 2 * rC * w + wSq) + rtSq * (-3 + rCSq - 2 * rC * w + 4 * wSq)) /
            (3. * mSqqSq);
        result.element({5, 1, 0, 2}) =
            -((-1 + wSq) * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (3. * mSqqSq);
        result.element({5, 1, 0, 3}) = (-2 * (mSqq + 2 * rtSq) * (rC - w) * (-1 + wSq)) / (3. * mSqqSq);
        result.element({5, 1, 1, 0}) = (rt * (-1 + wSq)) / mSqq;
        result.element({5, 1, 3, 0}) = (rt * (-1 + wSq)) / mSqq;
        result.element({5, 1, 5, 1}) =
            (4 * mSqqSq + 2 * mSqq * (rCSq + rtSq - 2 * rC * w + wSq) + rtSq * (-3 + rCSq - 2 * rC * w + 4 * wSq)) /
            (3. * mSqqSq);
        result.element({5, 1, 5, 2}) =
            -((-1 + wSq) * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (3. * mSqqSq);
        result.element({5, 1, 5, 3}) = (-2 * (mSqq + 2 * rtSq) * (rC - w) * (-1 + wSq)) / (3. * mSqqSq);
        result.element({5, 1, 7, 1}) =
            (4 * mSqqSq + 2 * mSqq * (rCSq + rtSq - 2 * rC * w + wSq) + rtSq * (-3 + rCSq - 2 * rC * w + 4 * wSq)) /
            (3. * mSqqSq);
        result.element({5, 1, 7, 2}) =
            -((-1 + wSq) * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (3. * mSqqSq);
        result.element({5, 1, 7, 3}) = (-2 * (mSqq + 2 * rtSq) * (rC - w) * (-1 + wSq)) / (3. * mSqqSq);
        result.element({5, 1, 9, 5}) = (4 * rt * (-2 + 3 * rC * w - wSq)) / mSqq;
        result.element({5, 1, 9, 6}) = (12 * rt * (rC - w)) / mSqq;
        result.element({5, 1, 9, 7}) = (-4 * rt * (rC - w) * (-1 + wSq)) / mSqq;
        result.element({5, 2, 0, 1}) =
            -((-1 + wSq) * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (3. * mSqqSq);
        result.element({5, 2, 0, 2}) =
            ((-1 + wSq) * (2 * mSqq * rCSq * (-1 + wSq) + rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)))) /
            (3. * mSqqSq);
        result.element({5, 2, 0, 3}) =
            ((-1 + wSq) * (2 * mSqq * rC * (-1 + wSq) + rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)))) /
            (3. * mSqqSq);
        result.element({5, 2, 1, 0}) = -((rt * (-1 + rC * w) * (-1 + wSq)) / mSqq);
        result.element({5, 2, 3, 0}) = -((rt * (-1 + rC * w) * (-1 + wSq)) / mSqq);
        result.element({5, 2, 5, 1}) =
            -((-1 + wSq) * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (3. * mSqqSq);
        result.element({5, 2, 5, 2}) =
            ((-1 + wSq) * (2 * mSqq * rCSq * (-1 + wSq) + rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)))) /
            (3. * mSqqSq);
        result.element({5, 2, 5, 3}) =
            ((-1 + wSq) * (2 * mSqq * rC * (-1 + wSq) + rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)))) /
            (3. * mSqqSq);
        result.element({5, 2, 7, 1}) =
            -((-1 + wSq) * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (3. * mSqqSq);
        result.element({5, 2, 7, 2}) =
            ((-1 + wSq) * (2 * mSqq * rCSq * (-1 + wSq) + rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)))) /
            (3. * mSqqSq);
        result.element({5, 2, 7, 3}) =
            ((-1 + wSq) * (2 * mSqq * rC * (-1 + wSq) + rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)))) /
            (3. * mSqqSq);
        result.element({5, 2, 9, 5}) = (-4 * rC * rt * w * (-1 + wSq)) / mSqq;
        result.element({5, 2, 9, 6}) = (-4 * rC * rt * (-1 + wSq)) / mSqq;
        result.element({5, 2, 9, 7}) = (4 * rC * rt * wSqm1Sq) / mSqq;
        result.element({5, 3, 0, 1}) = (-2 * (mSqq + 2 * rtSq) * (rC - w) * (-1 + wSq)) / (3. * mSqqSq);
        result.element({5, 3, 0, 2}) =
            ((-1 + wSq) * (2 * mSqq * rC * (-1 + wSq) + rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)))) /
            (3. * mSqqSq);
        result.element({5, 3, 0, 3}) =
            ((-1 + wSq) * (2 * mSqq * (-1 + wSq) + rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq))) / (3. * mSqqSq);
        result.element({5, 3, 1, 0}) = (rt * (-rC + w) * (-1 + wSq)) / mSqq;
        result.element({5, 3, 3, 0}) = (rt * (-rC + w) * (-1 + wSq)) / mSqq;
        result.element({5, 3, 5, 1}) = (-2 * (mSqq + 2 * rtSq) * (rC - w) * (-1 + wSq)) / (3. * mSqqSq);
        result.element({5, 3, 5, 2}) =
            ((-1 + wSq) * (2 * mSqq * rC * (-1 + wSq) + rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)))) /
            (3. * mSqqSq);
        result.element({5, 3, 5, 3}) =
            ((-1 + wSq) * (2 * mSqq * (-1 + wSq) + rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq))) / (3. * mSqqSq);
        result.element({5, 3, 7, 1}) = (-2 * (mSqq + 2 * rtSq) * (rC - w) * (-1 + wSq)) / (3. * mSqqSq);
        result.element({5, 3, 7, 2}) =
            ((-1 + wSq) * (2 * mSqq * rC * (-1 + wSq) + rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)))) /
            (3. * mSqqSq);
        result.element({5, 3, 7, 3}) =
            ((-1 + wSq) * (2 * mSqq * (-1 + wSq) + rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq))) / (3. * mSqqSq);
        result.element({5, 3, 9, 5}) = (-4 * rt * w * (-1 + wSq)) / mSqq;
        result.element({5, 3, 9, 6}) = (-4 * rt * (-1 + wSq)) / mSqq;
        result.element({5, 3, 9, 7}) = (4 * rt * wSqm1Sq) / mSqq;
        result.element({5, 4, 0, 4}) = (-2 * (2 * mSqq + rtSq) * (-1 + wSq)) / (3. * mSqq);
        result.element({5, 4, 5, 4}) = (2 * (2 * mSqq + rtSq) * (-1 + wSq)) / (3. * mSqq);
        result.element({5, 4, 7, 4}) = (-2 * (2 * mSqq + rtSq) * (-1 + wSq)) / (3. * mSqq);
        result.element({5, 4, 9, 5}) = (-8 * rC * rt * (-1 + wSq)) / mSqq;
        result.element({5, 4, 9, 6}) = (-8 * rt * (-1 + wSq)) / mSqq;
        result.element({6, 1, 2, 0}) = (rt * (-1 + wSq)) / mSqq;
        result.element({6, 1, 4, 0}) = (rt * (-1 + wSq)) / mSqq;
        result.element({6, 1, 6, 1}) =
            (4 * mSqqSq + 2 * mSqq * (rCSq + rtSq - 2 * rC * w + wSq) + rtSq * (-3 + rCSq - 2 * rC * w + 4 * wSq)) /
            (3. * mSqqSq);
        result.element({6, 1, 6, 2}) =
            -((-1 + wSq) * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (3. * mSqqSq);
        result.element({6, 1, 6, 3}) = (-2 * (mSqq + 2 * rtSq) * (rC - w) * (-1 + wSq)) / (3. * mSqqSq);
        result.element({6, 1, 8, 1}) =
            (4 * mSqqSq + 2 * mSqq * (rCSq + rtSq - 2 * rC * w + wSq) + rtSq * (-3 + rCSq - 2 * rC * w + 4 * wSq)) /
            (3. * mSqqSq);
        result.element({6, 1, 8, 2}) =
            -((-1 + wSq) * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (3. * mSqqSq);
        result.element({6, 1, 8, 3}) = (-2 * (mSqq + 2 * rtSq) * (rC - w) * (-1 + wSq)) / (3. * mSqqSq);
        result.element({6, 1, 10, 5}) = (4 * rt * (-2 + 3 * rC * w - wSq)) / mSqq;
        result.element({6, 1, 10, 6}) = (12 * rt * (rC - w)) / mSqq;
        result.element({6, 1, 10, 7}) = (-4 * rt * (rC - w) * (-1 + wSq)) / mSqq;
        result.element({6, 2, 2, 0}) = -((rt * (-1 + rC * w) * (-1 + wSq)) / mSqq);
        result.element({6, 2, 4, 0}) = -((rt * (-1 + rC * w) * (-1 + wSq)) / mSqq);
        result.element({6, 2, 6, 1}) =
            -((-1 + wSq) * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (3. * mSqqSq);
        result.element({6, 2, 6, 2}) =
            ((-1 + wSq) * (2 * mSqq * rCSq * (-1 + wSq) + rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)))) /
            (3. * mSqqSq);
        result.element({6, 2, 6, 3}) =
            ((-1 + wSq) * (2 * mSqq * rC * (-1 + wSq) + rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)))) /
            (3. * mSqqSq);
        result.element({6, 2, 8, 1}) =
            -((-1 + wSq) * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (3. * mSqqSq);
        result.element({6, 2, 8, 2}) =
            ((-1 + wSq) * (2 * mSqq * rCSq * (-1 + wSq) + rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)))) /
            (3. * mSqqSq);
        result.element({6, 2, 8, 3}) =
            ((-1 + wSq) * (2 * mSqq * rC * (-1 + wSq) + rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)))) /
            (3. * mSqqSq);
        result.element({6, 2, 10, 5}) = (-4 * rC * rt * w * (-1 + wSq)) / mSqq;
        result.element({6, 2, 10, 6}) = (-4 * rC * rt * (-1 + wSq)) / mSqq;
        result.element({6, 2, 10, 7}) = (4 * rC * rt * wSqm1Sq) / mSqq;
        result.element({6, 3, 2, 0}) = -((rt * (rC - w) * (-1 + wSq)) / mSqq);
        result.element({6, 3, 4, 0}) = -((rt * (rC - w) * (-1 + wSq)) / mSqq);
        result.element({6, 3, 6, 1}) = (-2 * (mSqq + 2 * rtSq) * (rC - w) * (-1 + wSq)) / (3. * mSqqSq);
        result.element({6, 3, 6, 2}) =
            ((-1 + wSq) * (2 * mSqq * rC * (-1 + wSq) + rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)))) /
            (3. * mSqqSq);
        result.element({6, 3, 6, 3}) =
            ((-1 + wSq) * (2 * mSqq * (-1 + wSq) + rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq))) / (3. * mSqqSq);
        result.element({6, 3, 8, 1}) = (-2 * (mSqq + 2 * rtSq) * (rC - w) * (-1 + wSq)) / (3. * mSqqSq);
        result.element({6, 3, 8, 2}) =
            ((-1 + wSq) * (2 * mSqq * rC * (-1 + wSq) + rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)))) /
            (3. * mSqqSq);
        result.element({6, 3, 8, 3}) =
            ((-1 + wSq) * (2 * mSqq * (-1 + wSq) + rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq))) / (3. * mSqqSq);
        result.element({6, 3, 10, 5}) = (-4 * rt * w * (-1 + wSq)) / mSqq;
        result.element({6, 3, 10, 6}) = (-4 * rt * (-1 + wSq)) / mSqq;
        result.element({6, 3, 10, 7}) = (4 * rt * wSqm1Sq) / mSqq;
        result.element({6, 4, 6, 4}) = (2 * (2 * mSqq + rtSq) * (-1 + wSq)) / (3. * mSqq);
        result.element({6, 4, 8, 4}) = (-2 * (2 * mSqq + rtSq) * (-1 + wSq)) / (3. * mSqq);
        result.element({6, 4, 10, 5}) = (8 * rC * rt * (-1 + wSq)) / mSqq;
        result.element({6, 4, 10, 6}) = (8 * rt * (-1 + wSq)) / mSqq;
        result.element({7, 1, 0, 1}) =
            (4 * mSqqSq + 2 * mSqq * (rCSq + rtSq - 2 * rC * w + wSq) + rtSq * (-3 + rCSq - 2 * rC * w + 4 * wSq)) /
            (3. * mSqqSq);
        result.element({7, 1, 0, 2}) =
            -((-1 + wSq) * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (3. * mSqqSq);
        result.element({7, 1, 0, 3}) = (-2 * (mSqq + 2 * rtSq) * (rC - w) * (-1 + wSq)) / (3. * mSqqSq);
        result.element({7, 1, 1, 0}) = (rt * (-1 + wSq)) / mSqq;
        result.element({7, 1, 3, 0}) = (rt * (-1 + wSq)) / mSqq;
        result.element({7, 1, 5, 1}) =
            (4 * mSqqSq + 2 * mSqq * (rCSq + rtSq - 2 * rC * w + wSq) + rtSq * (-3 + rCSq - 2 * rC * w + 4 * wSq)) /
            (3. * mSqqSq);
        result.element({7, 1, 5, 2}) =
            -((-1 + wSq) * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (3. * mSqqSq);
        result.element({7, 1, 5, 3}) = (-2 * (mSqq + 2 * rtSq) * (rC - w) * (-1 + wSq)) / (3. * mSqqSq);
        result.element({7, 1, 7, 1}) =
            (4 * mSqqSq + 2 * mSqq * (rCSq + rtSq - 2 * rC * w + wSq) + rtSq * (-3 + rCSq - 2 * rC * w + 4 * wSq)) /
            (3. * mSqqSq);
        result.element({7, 1, 7, 2}) =
            -((-1 + wSq) * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (3. * mSqqSq);
        result.element({7, 1, 7, 3}) = (-2 * (mSqq + 2 * rtSq) * (rC - w) * (-1 + wSq)) / (3. * mSqqSq);
        result.element({7, 1, 9, 5}) = (4 * rt * (-2 + 3 * rC * w - wSq)) / mSqq;
        result.element({7, 1, 9, 6}) = (12 * rt * (rC - w)) / mSqq;
        result.element({7, 1, 9, 7}) = (-4 * rt * (rC - w) * (-1 + wSq)) / mSqq;
        result.element({7, 2, 0, 1}) =
            -((-1 + wSq) * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (3. * mSqqSq);
        result.element({7, 2, 0, 2}) =
            ((-1 + wSq) * (2 * mSqq * rCSq * (-1 + wSq) + rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)))) /
            (3. * mSqqSq);
        result.element({7, 2, 0, 3}) =
            ((-1 + wSq) * (2 * mSqq * rC * (-1 + wSq) + rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)))) /
            (3. * mSqqSq);
        result.element({7, 2, 1, 0}) = -((rt * (-1 + rC * w) * (-1 + wSq)) / mSqq);
        result.element({7, 2, 3, 0}) = -((rt * (-1 + rC * w) * (-1 + wSq)) / mSqq);
        result.element({7, 2, 5, 1}) =
            -((-1 + wSq) * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (3. * mSqqSq);
        result.element({7, 2, 5, 2}) =
            ((-1 + wSq) * (2 * mSqq * rCSq * (-1 + wSq) + rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)))) /
            (3. * mSqqSq);
        result.element({7, 2, 5, 3}) =
            ((-1 + wSq) * (2 * mSqq * rC * (-1 + wSq) + rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)))) /
            (3. * mSqqSq);
        result.element({7, 2, 7, 1}) =
            -((-1 + wSq) * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (3. * mSqqSq);
        result.element({7, 2, 7, 2}) =
            ((-1 + wSq) * (2 * mSqq * rCSq * (-1 + wSq) + rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)))) /
            (3. * mSqqSq);
        result.element({7, 2, 7, 3}) =
            ((-1 + wSq) * (2 * mSqq * rC * (-1 + wSq) + rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)))) /
            (3. * mSqqSq);
        result.element({7, 2, 9, 5}) = (-4 * rC * rt * w * (-1 + wSq)) / mSqq;
        result.element({7, 2, 9, 6}) = (-4 * rC * rt * (-1 + wSq)) / mSqq;
        result.element({7, 2, 9, 7}) = (4 * rC * rt * wSqm1Sq) / mSqq;
        result.element({7, 3, 0, 1}) = (-2 * (mSqq + 2 * rtSq) * (rC - w) * (-1 + wSq)) / (3. * mSqqSq);
        result.element({7, 3, 0, 2}) =
            ((-1 + wSq) * (2 * mSqq * rC * (-1 + wSq) + rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)))) /
            (3. * mSqqSq);
        result.element({7, 3, 0, 3}) =
            ((-1 + wSq) * (2 * mSqq * (-1 + wSq) + rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq))) / (3. * mSqqSq);
        result.element({7, 3, 1, 0}) = (rt * (-rC + w) * (-1 + wSq)) / mSqq;
        result.element({7, 3, 3, 0}) = (rt * (-rC + w) * (-1 + wSq)) / mSqq;
        result.element({7, 3, 5, 1}) = (-2 * (mSqq + 2 * rtSq) * (rC - w) * (-1 + wSq)) / (3. * mSqqSq);
        result.element({7, 3, 5, 2}) =
            ((-1 + wSq) * (2 * mSqq * rC * (-1 + wSq) + rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)))) /
            (3. * mSqqSq);
        result.element({7, 3, 5, 3}) =
            ((-1 + wSq) * (2 * mSqq * (-1 + wSq) + rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq))) / (3. * mSqqSq);
        result.element({7, 3, 7, 1}) = (-2 * (mSqq + 2 * rtSq) * (rC - w) * (-1 + wSq)) / (3. * mSqqSq);
        result.element({7, 3, 7, 2}) =
            ((-1 + wSq) * (2 * mSqq * rC * (-1 + wSq) + rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)))) /
            (3. * mSqqSq);
        result.element({7, 3, 7, 3}) =
            ((-1 + wSq) * (2 * mSqq * (-1 + wSq) + rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq))) / (3. * mSqqSq);
        result.element({7, 3, 9, 5}) = (-4 * rt * w * (-1 + wSq)) / mSqq;
        result.element({7, 3, 9, 6}) = (-4 * rt * (-1 + wSq)) / mSqq;
        result.element({7, 3, 9, 7}) = (4 * rt * wSqm1Sq) / mSqq;
        result.element({7, 4, 0, 4}) = (2 * (2 * mSqq + rtSq) * (-1 + wSq)) / (3. * mSqq);
        result.element({7, 4, 5, 4}) = (-2 * (2 * mSqq + rtSq) * (-1 + wSq)) / (3. * mSqq);
        result.element({7, 4, 7, 4}) = (2 * (2 * mSqq + rtSq) * (-1 + wSq)) / (3. * mSqq);
        result.element({7, 4, 9, 5}) = (8 * rC * rt * (-1 + wSq)) / mSqq;
        result.element({7, 4, 9, 6}) = (8 * rt * (-1 + wSq)) / mSqq;
        result.element({8, 1, 2, 0}) = (rt * (-1 + wSq)) / mSqq;
        result.element({8, 1, 4, 0}) = (rt * (-1 + wSq)) / mSqq;
        result.element({8, 1, 6, 1}) =
            (4 * mSqqSq + 2 * mSqq * (rCSq + rtSq - 2 * rC * w + wSq) + rtSq * (-3 + rCSq - 2 * rC * w + 4 * wSq)) /
            (3. * mSqqSq);
        result.element({8, 1, 6, 2}) =
            -((-1 + wSq) * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (3. * mSqqSq);
        result.element({8, 1, 6, 3}) = (-2 * (mSqq + 2 * rtSq) * (rC - w) * (-1 + wSq)) / (3. * mSqqSq);
        result.element({8, 1, 8, 1}) =
            (4 * mSqqSq + 2 * mSqq * (rCSq + rtSq - 2 * rC * w + wSq) + rtSq * (-3 + rCSq - 2 * rC * w + 4 * wSq)) /
            (3. * mSqqSq);
        result.element({8, 1, 8, 2}) =
            -((-1 + wSq) * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (3. * mSqqSq);
        result.element({8, 1, 8, 3}) = (-2 * (mSqq + 2 * rtSq) * (rC - w) * (-1 + wSq)) / (3. * mSqqSq);
        result.element({8, 1, 10, 5}) = (4 * rt * (-2 + 3 * rC * w - wSq)) / mSqq;
        result.element({8, 1, 10, 6}) = (12 * rt * (rC - w)) / mSqq;
        result.element({8, 1, 10, 7}) = (-4 * rt * (rC - w) * (-1 + wSq)) / mSqq;
        result.element({8, 2, 2, 0}) = -((rt * (-1 + rC * w) * (-1 + wSq)) / mSqq);
        result.element({8, 2, 4, 0}) = -((rt * (-1 + rC * w) * (-1 + wSq)) / mSqq);
        result.element({8, 2, 6, 1}) =
            -((-1 + wSq) * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (3. * mSqqSq);
        result.element({8, 2, 6, 2}) =
            ((-1 + wSq) * (2 * mSqq * rCSq * (-1 + wSq) + rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)))) /
            (3. * mSqqSq);
        result.element({8, 2, 6, 3}) =
            ((-1 + wSq) * (2 * mSqq * rC * (-1 + wSq) + rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)))) /
            (3. * mSqqSq);
        result.element({8, 2, 8, 1}) =
            -((-1 + wSq) * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (3. * mSqqSq);
        result.element({8, 2, 8, 2}) =
            ((-1 + wSq) * (2 * mSqq * rCSq * (-1 + wSq) + rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)))) /
            (3. * mSqqSq);
        result.element({8, 2, 8, 3}) =
            ((-1 + wSq) * (2 * mSqq * rC * (-1 + wSq) + rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)))) /
            (3. * mSqqSq);
        result.element({8, 2, 10, 5}) = (-4 * rC * rt * w * (-1 + wSq)) / mSqq;
        result.element({8, 2, 10, 6}) = (-4 * rC * rt * (-1 + wSq)) / mSqq;
        result.element({8, 2, 10, 7}) = (4 * rC * rt * wSqm1Sq) / mSqq;
        result.element({8, 3, 2, 0}) = -((rt * (rC - w) * (-1 + wSq)) / mSqq);
        result.element({8, 3, 4, 0}) = -((rt * (rC - w) * (-1 + wSq)) / mSqq);
        result.element({8, 3, 6, 1}) = (-2 * (mSqq + 2 * rtSq) * (rC - w) * (-1 + wSq)) / (3. * mSqqSq);
        result.element({8, 3, 6, 2}) =
            ((-1 + wSq) * (2 * mSqq * rC * (-1 + wSq) + rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)))) /
            (3. * mSqqSq);
        result.element({8, 3, 6, 3}) =
            ((-1 + wSq) * (2 * mSqq * (-1 + wSq) + rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq))) / (3. * mSqqSq);
        result.element({8, 3, 8, 1}) = (-2 * (mSqq + 2 * rtSq) * (rC - w) * (-1 + wSq)) / (3. * mSqqSq);
        result.element({8, 3, 8, 2}) =
            ((-1 + wSq) * (2 * mSqq * rC * (-1 + wSq) + rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)))) /
            (3. * mSqqSq);
        result.element({8, 3, 8, 3}) =
            ((-1 + wSq) * (2 * mSqq * (-1 + wSq) + rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq))) / (3. * mSqqSq);
        result.element({8, 3, 10, 5}) = (-4 * rt * w * (-1 + wSq)) / mSqq;
        result.element({8, 3, 10, 6}) = (-4 * rt * (-1 + wSq)) / mSqq;
        result.element({8, 3, 10, 7}) = (4 * rt * wSqm1Sq) / mSqq;
        result.element({8, 4, 6, 4}) = (-2 * (2 * mSqq + rtSq) * (-1 + wSq)) / (3. * mSqq);
        result.element({8, 4, 8, 4}) = (2 * (2 * mSqq + rtSq) * (-1 + wSq)) / (3. * mSqq);
        result.element({8, 4, 10, 5}) = (-8 * rC * rt * (-1 + wSq)) / mSqq;
        result.element({8, 4, 10, 6}) = (-8 * rt * (-1 + wSq)) / mSqq;
        result.element({9, 5, 0, 1}) = (4 * rt * (-2 + 3 * rC * w - wSq)) / mSqq;
        result.element({9, 5, 0, 2}) = (-4 * rC * rt * w * (-1 + wSq)) / mSqq;
        result.element({9, 5, 0, 3}) = (-4 * rt * w * (-1 + wSq)) / mSqq;
        result.element({9, 5, 0, 4}) = (8 * rC * rt * (-1 + wSq)) / mSqq;
        result.element({9, 5, 5, 1}) = (4 * rt * (-2 + 3 * rC * w - wSq)) / mSqq;
        result.element({9, 5, 5, 2}) = (-4 * rC * rt * w * (-1 + wSq)) / mSqq;
        result.element({9, 5, 5, 3}) = (-4 * rt * w * (-1 + wSq)) / mSqq;
        result.element({9, 5, 5, 4}) = (-8 * rC * rt * (-1 + wSq)) / mSqq;
        result.element({9, 5, 7, 1}) = (4 * rt * (-2 + 3 * rC * w - wSq)) / mSqq;
        result.element({9, 5, 7, 2}) = (-4 * rC * rt * w * (-1 + wSq)) / mSqq;
        result.element({9, 5, 7, 3}) = (-4 * rt * w * (-1 + wSq)) / mSqq;
        result.element({9, 5, 7, 4}) = (8 * rC * rt * (-1 + wSq)) / mSqq;
        result.element({9, 5, 9, 5}) = (16 * (mSqqSq * wSq + 4 * rtSq * (1 - 2 * rC * w + rCSq * (-1 + 2 * wSq)) +
                                              2 * mSqq * (1 - 2 * rC * w + rtSq * wSq + rCSq * (-1 + 2 * wSq)))) /
                                       (3. * mSqqSq);
        result.element({9, 5, 9, 6}) =
            (16 * (mSqq + 2 * rtSq) * (-4 * rC + (2 + mSqq) * w + 2 * rCSq * w)) / (3. * mSqqSq);
        result.element({9, 5, 9, 7}) = (-16 * (mSqq + 2 * rtSq) * w * (-1 + wSq)) / (3. * mSqq);
        result.element({9, 6, 0, 1}) = (12 * rt * (rC - w)) / mSqq;
        result.element({9, 6, 0, 2}) = (-4 * rC * rt * (-1 + wSq)) / mSqq;
        result.element({9, 6, 0, 3}) = (-4 * rt * (-1 + wSq)) / mSqq;
        result.element({9, 6, 0, 4}) = (8 * rt * (-1 + wSq)) / mSqq;
        result.element({9, 6, 5, 1}) = (12 * rt * (rC - w)) / mSqq;
        result.element({9, 6, 5, 2}) = (-4 * rC * rt * (-1 + wSq)) / mSqq;
        result.element({9, 6, 5, 3}) = (-4 * rt * (-1 + wSq)) / mSqq;
        result.element({9, 6, 5, 4}) = (-8 * rt * (-1 + wSq)) / mSqq;
        result.element({9, 6, 7, 1}) = (12 * rt * (rC - w)) / mSqq;
        result.element({9, 6, 7, 2}) = (-4 * rC * rt * (-1 + wSq)) / mSqq;
        result.element({9, 6, 7, 3}) = (-4 * rt * (-1 + wSq)) / mSqq;
        result.element({9, 6, 7, 4}) = (8 * rt * (-1 + wSq)) / mSqq;
        result.element({9, 6, 9, 5}) =
            (16 * (mSqq + 2 * rtSq) * (-4 * rC + (2 + mSqq) * w + 2 * rCSq * w)) / (3. * mSqqSq);
        result.element({9, 6, 9, 6}) = (16 * (mSqqSq + 4 * rtSq * (-1 + rCSq - 2 * rC * w + 2 * wSq) +
                                              2 * mSqq * (-1 + rCSq + rtSq - 2 * rC * w + 2 * wSq))) /
                                       (3. * mSqqSq);
        result.element({9, 6, 9, 7}) = (-16 * (mSqq + 2 * rtSq) * (-1 + wSq)) / (3. * mSqq);
        result.element({9, 7, 0, 1}) = (-4 * rt * (rC - w) * (-1 + wSq)) / mSqq;
        result.element({9, 7, 0, 2}) = (4 * rC * rt * wSqm1Sq) / mSqq;
        result.element({9, 7, 0, 3}) = (4 * rt * wSqm1Sq) / mSqq;
        result.element({9, 7, 5, 1}) = (-4 * rt * (rC - w) * (-1 + wSq)) / mSqq;
        result.element({9, 7, 5, 2}) = (4 * rC * rt * wSqm1Sq) / mSqq;
        result.element({9, 7, 5, 3}) = (4 * rt * wSqm1Sq) / mSqq;
        result.element({9, 7, 7, 1}) = (-4 * rt * (rC - w) * (-1 + wSq)) / mSqq;
        result.element({9, 7, 7, 2}) = (4 * rC * rt * wSqm1Sq) / mSqq;
        result.element({9, 7, 7, 3}) = (4 * rt * wSqm1Sq) / mSqq;
        result.element({9, 7, 9, 5}) = (-16 * (mSqq + 2 * rtSq) * w * (-1 + wSq)) / (3. * mSqq);
        result.element({9, 7, 9, 6}) = (-16 * (mSqq + 2 * rtSq) * (-1 + wSq)) / (3. * mSqq);
        result.element({9, 7, 9, 7}) = (16 * (mSqq + 2 * rtSq) * wSqm1Sq) / (3. * mSqq);
        result.element({10, 5, 6, 1}) = (4 * rt * (-2 + 3 * rC * w - wSq)) / mSqq;
        result.element({10, 5, 6, 2}) = (-4 * rC * rt * w * (-1 + wSq)) / mSqq;
        result.element({10, 5, 6, 3}) = (-4 * rt * w * (-1 + wSq)) / mSqq;
        result.element({10, 5, 6, 4}) = (8 * rC * rt * (-1 + wSq)) / mSqq;
        result.element({10, 5, 8, 1}) = (4 * rt * (-2 + 3 * rC * w - wSq)) / mSqq;
        result.element({10, 5, 8, 2}) = (-4 * rC * rt * w * (-1 + wSq)) / mSqq;
        result.element({10, 5, 8, 3}) = (-4 * rt * w * (-1 + wSq)) / mSqq;
        result.element({10, 5, 8, 4}) = (-8 * rC * rt * (-1 + wSq)) / mSqq;
        result.element({10, 5, 10, 5}) = (16 * (mSqqSq * wSq + 4 * rtSq * (1 - 2 * rC * w + rCSq * (-1 + 2 * wSq)) +
                                                2 * mSqq * (1 - 2 * rC * w + rtSq * wSq + rCSq * (-1 + 2 * wSq)))) /
                                         (3. * mSqqSq);
        result.element({10, 5, 10, 6}) =
            (16 * (mSqq + 2 * rtSq) * (-4 * rC + (2 + mSqq) * w + 2 * rCSq * w)) / (3. * mSqqSq);
        result.element({10, 5, 10, 7}) = (-16 * (mSqq + 2 * rtSq) * w * (-1 + wSq)) / (3. * mSqq);
        result.element({10, 6, 6, 1}) = (12 * rt * (rC - w)) / mSqq;
        result.element({10, 6, 6, 2}) = (-4 * rC * rt * (-1 + wSq)) / mSqq;
        result.element({10, 6, 6, 3}) = (-4 * rt * (-1 + wSq)) / mSqq;
        result.element({10, 6, 6, 4}) = (8 * rt * (-1 + wSq)) / mSqq;
        result.element({10, 6, 8, 1}) = (12 * rt * (rC - w)) / mSqq;
        result.element({10, 6, 8, 2}) = (-4 * rC * rt * (-1 + wSq)) / mSqq;
        result.element({10, 6, 8, 3}) = (-4 * rt * (-1 + wSq)) / mSqq;
        result.element({10, 6, 8, 4}) = (-8 * rt * (-1 + wSq)) / mSqq;
        result.element({10, 6, 10, 5}) =
            (16 * (mSqq + 2 * rtSq) * (-4 * rC + (2 + mSqq) * w + 2 * rCSq * w)) / (3. * mSqqSq);
        result.element({10, 6, 10, 6}) = (16 * (mSqqSq + 4 * rtSq * (-1 + rCSq - 2 * rC * w + 2 * wSq) +
                                                2 * mSqq * (-1 + rCSq + rtSq - 2 * rC * w + 2 * wSq))) /
                                         (3. * mSqqSq);
        result.element({10, 6, 10, 7}) = (-16 * (mSqq + 2 * rtSq) * (-1 + wSq)) / (3. * mSqq);
        result.element({10, 7, 6, 1}) = (-4 * rt * (rC - w) * (-1 + wSq)) / mSqq;
        result.element({10, 7, 6, 2}) = (4 * rC * rt * wSqm1Sq) / mSqq;
        result.element({10, 7, 6, 3}) = (4 * rt * wSqm1Sq) / mSqq;
        result.element({10, 7, 8, 1}) = (-4 * rt * (rC - w) * (-1 + wSq)) / mSqq;
        result.element({10, 7, 8, 2}) = (4 * rC * rt * wSqm1Sq) / mSqq;
        result.element({10, 7, 8, 3}) = (4 * rt * wSqm1Sq) / mSqq;
        result.element({10, 7, 10, 5}) = (-16 * (mSqq + 2 * rtSq) * w * (-1 + wSq)) / (3. * mSqq);
        result.element({10, 7, 10, 6}) = (-16 * (mSqq + 2 * rtSq) * (-1 + wSq)) / (3. * mSqq);
        result.element({10, 7, 10, 7}) = (16 * (mSqq + 2 * rtSq) * wSqm1Sq) / (3. * mSqq);

        result *= RateNorm;

        return result;
    }

} // namespace Hammer
