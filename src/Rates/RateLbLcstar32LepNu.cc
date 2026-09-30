///
/// @file  RateLbLcstar32LepNu.cc
/// @brief \f$ \Lambda_b \rightarrow L_c^*(2625) \tau\nu \f$ total rate
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/Rates/RateLbLcstar32LepNu.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    RateLbLcstar32LepNu::RateLbLcstar32LepNu() {
        // Create tensor rank and dimensions
        IndexList dims{{11, 16, 11, 16, _nPoints}};
        string name{"RateLbLcstar32LepNuQ2"};
        auto& pdg = PID::instance();

        addProcessSignature(-PID::LAMBDAB, {PID::LAMBDACSTAR32MINUS, PID::NU_TAU, PID::ANTITAU});
        addIntegrationBoundaries({PS::makeQ2Function(
            pdg.getMass(PID::ANTITAU), pdg.getMass(PID::LAMBDAB) - pdg.getMass(PID::LAMBDACSTAR32MINUS))});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCTAUNU, FF_LBLCSTAR32, WILSON_BCTAUNU_HC,
                                                          FF_LBLCSTAR32_HC, INTEGRATION_INDEX})});

        addProcessSignature(-PID::LAMBDAB, {PID::LAMBDACSTAR32MINUS, PID::NU_MU, PID::ANTIMUON});
        addIntegrationBoundaries({PS::makeQ2Function(
            pdg.getMass(PID::ANTIMUON), pdg.getMass(PID::LAMBDAB) - pdg.getMass(PID::LAMBDACSTAR32MINUS))});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCMUNU, FF_LBLCSTAR32, WILSON_BCMUNU_HC,
                                                          FF_LBLCSTAR32_HC, INTEGRATION_INDEX})});

        addProcessSignature(-PID::LAMBDAB, {PID::LAMBDACSTAR32MINUS, PID::NU_E, PID::POSITRON});
        addIntegrationBoundaries({PS::makeQ2Function(
            pdg.getMass(PID::POSITRON), pdg.getMass(PID::LAMBDAB) - pdg.getMass(PID::LAMBDACSTAR32MINUS))});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCENU, FF_LBLCSTAR32, WILSON_BCENU_HC,
                                                          FF_LBLCSTAR32_HC, INTEGRATION_INDEX})});

        setSignatureIndex();
    }

    // NOLINTBEGIN(readability-function-size)
    Tensor RateLbLcstar32LepNu::evalAtPSPoint(const vector<double>& point) {
        auto labs = getTensor().labels();
        labs.pop_back();
        auto dimensions = getTensor().dims();
        dimensions.pop_back();
        Tensor result{"RateLbLcstar32LepNu", MD::makeEmptySparse(dimensions, labs)};

        const double Mb = masses()[0];
        const double Mc = masses()[1];
        const double Mt = masses()[3];

        // kinematic objects
        const double Sqq = point[0];
        const double mSqq = Sqq / (Mb * Mb);
        const double w = (Mb * Mb + Mc * Mc - Sqq) / (2. * Mb * Mc);
        const double rC = Mc / Mb;
        const double rt = Mt / Mb;

        const double mSqq2 = mSqq * mSqq;
        const double wSq = w * w;
        const double w2m1 = w * w - 1;
        const double wm1Sq = (w - 1) * (w - 1);
        const double wp1Sq = (w + 1) * (w + 1);
        const double rCSq = rC * rC;
        const double rtSq = rt * rt;
        const double wp1Cu = wp1Sq * (w + 1);
        const double w2m1Sq = w2m1 * w2m1;
        const double wCu = wSq * w;

        const double RateNorm =
            (pow(GFermi, 2.) * pow(Mb, 3.) * pow(rC, 2.) * pow(mSqq - rtSq, 2.) * sqrt(-1 + wSq)) / (32. * mSqq * pi3);

        // set non-zero tensor elements
        result.element({0, 2, 0, 2}) = (4 *
                                        (2 * mSqq2 * (-1 + w) + mSqq * (1 + 2 * rC + rCSq + rtSq) * (-1 + w) +
                                         rtSq * (1 + 2 * w - 2 * rC * (2 + w) + rCSq * (1 + 2 * w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({0, 2, 0, 3}) =
            (2 * (-1 + w) * wp1Sq *
             (2 * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2 * rC * (2 + w) + rCSq * (-1 + 4 * w)))) /
            (9. * mSqq2);
        result.element({0, 2, 0, 4}) =
            (2 * (-1 + w) * wp1Sq *
             (2 * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3 * rCSq + 4 * w - 2 * rC * (2 + w)))) /
            (9. * mSqq2);
        result.element({0, 2, 0, 5}) = (-2 *
                                        (2 * mSqq2 + mSqq * (2 * rCSq + rtSq - 2 * rC * (-1 + w) - 2 * w) +
                                         rtSq * (-3 + rCSq - 4 * w + 2 * rC * (2 + w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({0, 2, 1, 2}) = (4 *
                                        (2 * mSqq2 * (-1 + w) + mSqq * (1 + 2 * rC + rCSq + rtSq) * (-1 + w) +
                                         rtSq * (1 + 2 * w - 2 * rC * (2 + w) + rCSq * (1 + 2 * w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({0, 2, 1, 3}) =
            (2 * (-1 + w) * wp1Sq *
             (2 * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2 * rC * (2 + w) + rCSq * (-1 + 4 * w)))) /
            (9. * mSqq2);
        result.element({0, 2, 1, 4}) =
            (2 * (-1 + w) * wp1Sq *
             (2 * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3 * rCSq + 4 * w - 2 * rC * (2 + w)))) /
            (9. * mSqq2);
        result.element({0, 2, 1, 5}) = (-2 *
                                        (2 * mSqq2 + mSqq * (2 * rCSq + rtSq - 2 * rC * (-1 + w) - 2 * w) +
                                         rtSq * (-3 + rCSq - 4 * w + 2 * rC * (2 + w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({0, 2, 2, 2}) = (4 *
                                        (2 * mSqq2 * (-1 + w) + mSqq * (1 + 2 * rC + rCSq + rtSq) * (-1 + w) +
                                         rtSq * (1 + 2 * w - 2 * rC * (2 + w) + rCSq * (1 + 2 * w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({0, 2, 2, 3}) =
            (2 * (-1 + w) * wp1Sq *
             (2 * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2 * rC * (2 + w) + rCSq * (-1 + 4 * w)))) /
            (9. * mSqq2);
        result.element({0, 2, 2, 4}) =
            (2 * (-1 + w) * wp1Sq *
             (2 * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3 * rCSq + 4 * w - 2 * rC * (2 + w)))) /
            (9. * mSqq2);
        result.element({0, 2, 2, 5}) = (-2 *
                                        (2 * mSqq2 + mSqq * (2 * rCSq + rtSq - 2 * rC * (-1 + w) - 2 * w) +
                                         rtSq * (-3 + rCSq - 4 * w + 2 * rC * (2 + w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({0, 2, 3, 0}) = (2 * (-1 + rC) * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({0, 2, 4, 0}) = (2 * (-1 + rC) * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({0, 2, 5, 10}) = (-8 * (1 + rC) * rt * wm1Sq * (1 + w)) / mSqq;
        result.element({0, 2, 5, 11}) = (-8 * rt * wm1Sq * (1 + w) * (-3 + rC * (-1 + 2 * w))) / (3. * mSqq);
        result.element({0, 2, 5, 12}) = (-8 * rt * (1 + 3 * rC - 2 * w) * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({0, 2, 5, 13}) = (-16 * rt * (1 + rC - w) * w2m1) / (3. * mSqq);
        result.element({0, 2, 5, 14}) = (-8 * rt * (1 + 2 * rC - w) * w2m1) / (3. * mSqq);
        result.element({0, 2, 5, 15}) = (8 * (1 + rC) * rt * w2m1Sq) / (3. * mSqq);
        result.element({0, 3, 0, 2}) =
            (2 * (-1 + w) * wp1Sq *
             (2 * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2 * rC * (2 + w) + rCSq * (-1 + 4 * w)))) /
            (9. * mSqq2);
        result.element({0, 3, 0, 3}) =
            (2 * (-1 + w) * wp1Sq * (rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)) + 2 * mSqq * rCSq * w2m1)) /
            (9. * mSqq2);
        result.element({0, 3, 0, 4}) =
            (2 * (-1 + w) * wp1Sq * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) /
            (9. * mSqq2);
        result.element({0, 3, 0, 5}) =
            (-2 * (-1 + w) * wp1Sq * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (9. * mSqq2);
        result.element({0, 3, 1, 2}) =
            (2 * (-1 + w) * wp1Sq *
             (2 * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2 * rC * (2 + w) + rCSq * (-1 + 4 * w)))) /
            (9. * mSqq2);
        result.element({0, 3, 1, 3}) =
            (2 * (-1 + w) * wp1Sq * (rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)) + 2 * mSqq * rCSq * w2m1)) /
            (9. * mSqq2);
        result.element({0, 3, 1, 4}) =
            (2 * (-1 + w) * wp1Sq * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) /
            (9. * mSqq2);
        result.element({0, 3, 1, 5}) =
            (-2 * (-1 + w) * wp1Sq * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (9. * mSqq2);
        result.element({0, 3, 2, 2}) =
            (2 * (-1 + w) * wp1Sq *
             (2 * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2 * rC * (2 + w) + rCSq * (-1 + 4 * w)))) /
            (9. * mSqq2);
        result.element({0, 3, 2, 3}) =
            (2 * (-1 + w) * wp1Sq * (rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)) + 2 * mSqq * rCSq * w2m1)) /
            (9. * mSqq2);
        result.element({0, 3, 2, 4}) =
            (2 * (-1 + w) * wp1Sq * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) /
            (9. * mSqq2);
        result.element({0, 3, 2, 5}) =
            (-2 * (-1 + w) * wp1Sq * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (9. * mSqq2);
        result.element({0, 3, 3, 0}) = (2 * rt * (-1 + w) * wp1Sq * (-1 + rC * w)) / (3. * mSqq);
        result.element({0, 3, 4, 0}) = (2 * rt * (-1 + w) * wp1Sq * (-1 + rC * w)) / (3. * mSqq);
        result.element({0, 3, 5, 10}) = (-8 * rC * rt * w2m1Sq) / (3. * mSqq);
        result.element({0, 3, 5, 11}) = (8 * rC * rt * w2m1Sq) / (3. * mSqq);
        result.element({0, 3, 5, 12}) = (-8 * rC * rt * w2m1Sq) / (3. * mSqq);
        result.element({0, 3, 5, 13}) = (-8 * rC * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({0, 3, 5, 14}) = (-8 * rC * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({0, 3, 5, 15}) = (8 * rC * rt * wm1Sq * wp1Cu) / (3. * mSqq);
        result.element({0, 4, 0, 2}) =
            (2 * (-1 + w) * wp1Sq *
             (2 * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3 * rCSq + 4 * w - 2 * rC * (2 + w)))) /
            (9. * mSqq2);
        result.element({0, 4, 0, 3}) =
            (2 * (-1 + w) * wp1Sq * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) /
            (9. * mSqq2);
        result.element({0, 4, 0, 4}) =
            (2 * (-1 + w) * wp1Sq * (rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq) + 2 * mSqq * w2m1)) / (9. * mSqq2);
        result.element({0, 4, 0, 5}) = (-4 * (mSqq + 2 * rtSq) * (rC - w) * (-1 + w) * wp1Sq) / (9. * mSqq2);
        result.element({0, 4, 1, 2}) =
            (2 * (-1 + w) * wp1Sq *
             (2 * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3 * rCSq + 4 * w - 2 * rC * (2 + w)))) /
            (9. * mSqq2);
        result.element({0, 4, 1, 3}) =
            (2 * (-1 + w) * wp1Sq * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) /
            (9. * mSqq2);
        result.element({0, 4, 1, 4}) =
            (2 * (-1 + w) * wp1Sq * (rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq) + 2 * mSqq * w2m1)) / (9. * mSqq2);
        result.element({0, 4, 1, 5}) = (-4 * (mSqq + 2 * rtSq) * (rC - w) * (-1 + w) * wp1Sq) / (9. * mSqq2);
        result.element({0, 4, 2, 2}) =
            (2 * (-1 + w) * wp1Sq *
             (2 * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3 * rCSq + 4 * w - 2 * rC * (2 + w)))) /
            (9. * mSqq2);
        result.element({0, 4, 2, 3}) =
            (2 * (-1 + w) * wp1Sq * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) /
            (9. * mSqq2);
        result.element({0, 4, 2, 4}) =
            (2 * (-1 + w) * wp1Sq * (rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq) + 2 * mSqq * w2m1)) / (9. * mSqq2);
        result.element({0, 4, 2, 5}) = (-4 * (mSqq + 2 * rtSq) * (rC - w) * (-1 + w) * wp1Sq) / (9. * mSqq2);
        result.element({0, 4, 3, 0}) = (2 * rt * (rC - w) * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({0, 4, 4, 0}) = (2 * rt * (rC - w) * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({0, 4, 5, 10}) = (-8 * rt * w2m1Sq) / (3. * mSqq);
        result.element({0, 4, 5, 11}) = (8 * rt * w2m1Sq) / (3. * mSqq);
        result.element({0, 4, 5, 12}) = (-8 * rt * w2m1Sq) / (3. * mSqq);
        result.element({0, 4, 5, 13}) = (-8 * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({0, 4, 5, 14}) = (-8 * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({0, 4, 5, 15}) = (8 * rt * wm1Sq * wp1Cu) / (3. * mSqq);
        result.element({0, 5, 0, 2}) = (-2 *
                                        (2 * mSqq2 + mSqq * (2 * rCSq + rtSq - 2 * rC * (-1 + w) - 2 * w) +
                                         rtSq * (-3 + rCSq - 4 * w + 2 * rC * (2 + w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({0, 5, 0, 3}) =
            (-2 * (-1 + w) * wp1Sq * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (9. * mSqq2);
        result.element({0, 5, 0, 4}) = (-4 * (mSqq + 2 * rtSq) * (rC - w) * (-1 + w) * wp1Sq) / (9. * mSqq2);
        result.element({0, 5, 0, 5}) =
            (2 * (1 + w) *
             (4 * mSqq2 + 2 * mSqq * (rCSq + rtSq - 2 * rC * w + wSq) + rtSq * (-3 + rCSq - 2 * rC * w + 4 * wSq))) /
            (9. * mSqq2);
        result.element({0, 5, 1, 2}) = (-2 *
                                        (2 * mSqq2 + mSqq * (2 * rCSq + rtSq - 2 * rC * (-1 + w) - 2 * w) +
                                         rtSq * (-3 + rCSq - 4 * w + 2 * rC * (2 + w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({0, 5, 1, 3}) =
            (-2 * (-1 + w) * wp1Sq * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (9. * mSqq2);
        result.element({0, 5, 1, 4}) = (-4 * (mSqq + 2 * rtSq) * (rC - w) * (-1 + w) * wp1Sq) / (9. * mSqq2);
        result.element({0, 5, 1, 5}) =
            (2 * (1 + w) *
             (4 * mSqq2 + 2 * mSqq * (rCSq + rtSq - 2 * rC * w + wSq) + rtSq * (-3 + rCSq - 2 * rC * w + 4 * wSq))) /
            (9. * mSqq2);
        result.element({0, 5, 2, 2}) = (-2 *
                                        (2 * mSqq2 + mSqq * (2 * rCSq + rtSq - 2 * rC * (-1 + w) - 2 * w) +
                                         rtSq * (-3 + rCSq - 4 * w + 2 * rC * (2 + w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({0, 5, 2, 3}) =
            (-2 * (-1 + w) * wp1Sq * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (9. * mSqq2);
        result.element({0, 5, 2, 4}) = (-4 * (mSqq + 2 * rtSq) * (rC - w) * (-1 + w) * wp1Sq) / (9. * mSqq2);
        result.element({0, 5, 2, 5}) =
            (2 * (1 + w) *
             (4 * mSqq2 + 2 * mSqq * (rCSq + rtSq - 2 * rC * w + wSq) + rtSq * (-3 + rCSq - 2 * rC * w + 4 * wSq))) /
            (9. * mSqq2);
        result.element({0, 5, 3, 0}) = (-2 * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({0, 5, 4, 0}) = (-2 * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({0, 5, 5, 10}) = (8 * rt * (1 + 2 * rC - w) * w2m1) / (3. * mSqq);
        result.element({0, 5, 5, 11}) = (8 * (1 + rC) * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({0, 5, 5, 12}) = (16 * rt * (rC - w) * w2m1) / (3. * mSqq);
        result.element({0, 5, 5, 13}) = (8 * rt * (-1 + 3 * rC - 2 * w) * (1 + w)) / (3. * mSqq);
        result.element({0, 5, 5, 14}) = (8 * rt * (rC - w) * (1 + w)) / mSqq;
        result.element({0, 5, 5, 15}) = (-8 * rt * (rC - w) * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({0, 6, 0, 6}) = (4 *
                                        (2 * mSqq2 * (1 + w) + mSqq * (1 - 2 * rC + rCSq + rtSq) * (1 + w) +
                                         rtSq * (-1 + 2 * rC * (-2 + w) + 2 * w + rCSq * (-1 + 2 * w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({0, 6, 0, 7}) =
            (2 * wm1Sq * (1 + w) *
             (2 * mSqq * (-1 + rC) * rC * (1 + w) + rtSq * (-3 + 2 * rC * (-2 + w) + rCSq * (1 + 4 * w)))) /
            (9. * mSqq2);
        result.element({0, 6, 0, 8}) =
            (2 * wm1Sq * (1 + w) *
             (rtSq * (-1 + 3 * rCSq - 2 * rC * (-2 + w) - 4 * w) + 2 * mSqq * (-1 + rC) * (1 + w))) /
            (9. * mSqq2);
        result.element({0, 6, 0, 9}) = (-2 *
                                        (2 * mSqq2 + rtSq * (-3 + rCSq + 2 * rC * (-2 + w) + 4 * w) +
                                         mSqq * (2 * rCSq + rtSq + 2 * w - 2 * rC * (1 + w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({0, 6, 1, 6}) = (4 *
                                        (2 * mSqq2 * (1 + w) + mSqq * (1 - 2 * rC + rCSq + rtSq) * (1 + w) +
                                         rtSq * (-1 + 2 * rC * (-2 + w) + 2 * w + rCSq * (-1 + 2 * w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({0, 6, 1, 7}) =
            (2 * wm1Sq * (1 + w) *
             (2 * mSqq * (-1 + rC) * rC * (1 + w) + rtSq * (-3 + 2 * rC * (-2 + w) + rCSq * (1 + 4 * w)))) /
            (9. * mSqq2);
        result.element({0, 6, 1, 8}) =
            (2 * wm1Sq * (1 + w) *
             (rtSq * (-1 + 3 * rCSq - 2 * rC * (-2 + w) - 4 * w) + 2 * mSqq * (-1 + rC) * (1 + w))) /
            (9. * mSqq2);
        result.element({0, 6, 1, 9}) = (-2 *
                                        (2 * mSqq2 + rtSq * (-3 + rCSq + 2 * rC * (-2 + w) + 4 * w) +
                                         mSqq * (2 * rCSq + rtSq + 2 * w - 2 * rC * (1 + w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({0, 6, 2, 6}) = (-4 *
                                        (2 * mSqq2 * (1 + w) + mSqq * (1 - 2 * rC + rCSq + rtSq) * (1 + w) +
                                         rtSq * (-1 + 2 * rC * (-2 + w) + 2 * w + rCSq * (-1 + 2 * w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({0, 6, 2, 7}) =
            (-2 * wm1Sq * (1 + w) *
             (2 * mSqq * (-1 + rC) * rC * (1 + w) + rtSq * (-3 + 2 * rC * (-2 + w) + rCSq * (1 + 4 * w)))) /
            (9. * mSqq2);
        result.element({0, 6, 2, 8}) =
            (-2 * wm1Sq * (1 + w) *
             (rtSq * (-1 + 3 * rCSq - 2 * rC * (-2 + w) - 4 * w) + 2 * mSqq * (-1 + rC) * (1 + w))) /
            (9. * mSqq2);
        result.element({0, 6, 2, 9}) = (2 *
                                        (2 * mSqq2 + rtSq * (-3 + rCSq + 2 * rC * (-2 + w) + 4 * w) +
                                         mSqq * (2 * rCSq + rtSq + 2 * w - 2 * rC * (1 + w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({0, 6, 3, 1}) = (-2 * (1 + rC) * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({0, 6, 4, 1}) = (2 * (1 + rC) * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({0, 6, 5, 10}) = (-8 * (-1 + rC) * rt * (-1 + w) * wp1Sq) / mSqq;
        result.element({0, 6, 5, 11}) = (-16 * rC * rt * w2m1Sq) / (3. * mSqq);
        result.element({0, 6, 5, 12}) = (-16 * rt * w2m1Sq) / (3. * mSqq);
        result.element({0, 6, 5, 13}) = (16 * rt * (-1 + rC - w) * w2m1) / (3. * mSqq);
        result.element({0, 6, 5, 14}) = (-8 * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({0, 7, 0, 6}) =
            (2 * wm1Sq * (1 + w) *
             (2 * mSqq * (-1 + rC) * rC * (1 + w) + rtSq * (-3 + 2 * rC * (-2 + w) + rCSq * (1 + 4 * w)))) /
            (9. * mSqq2);
        result.element({0, 7, 0, 7}) =
            (2 * wm1Sq * (1 + w) * (rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)) + 2 * mSqq * rCSq * w2m1)) /
            (9. * mSqq2);
        result.element({0, 7, 0, 8}) =
            (2 * wm1Sq * (1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) /
            (9. * mSqq2);
        result.element({0, 7, 0, 9}) =
            (-2 * wm1Sq * (1 + w) * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (9. * mSqq2);
        result.element({0, 7, 1, 6}) =
            (2 * wm1Sq * (1 + w) *
             (2 * mSqq * (-1 + rC) * rC * (1 + w) + rtSq * (-3 + 2 * rC * (-2 + w) + rCSq * (1 + 4 * w)))) /
            (9. * mSqq2);
        result.element({0, 7, 1, 7}) =
            (2 * wm1Sq * (1 + w) * (rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)) + 2 * mSqq * rCSq * w2m1)) /
            (9. * mSqq2);
        result.element({0, 7, 1, 8}) =
            (2 * wm1Sq * (1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) /
            (9. * mSqq2);
        result.element({0, 7, 1, 9}) =
            (-2 * wm1Sq * (1 + w) * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (9. * mSqq2);
        result.element({0, 7, 2, 6}) =
            (-2 * wm1Sq * (1 + w) *
             (2 * mSqq * (-1 + rC) * rC * (1 + w) + rtSq * (-3 + 2 * rC * (-2 + w) + rCSq * (1 + 4 * w)))) /
            (9. * mSqq2);
        result.element({0, 7, 2, 7}) =
            (-2 * wm1Sq * (1 + w) * (rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)) + 2 * mSqq * rCSq * w2m1)) /
            (9. * mSqq2);
        result.element({0, 7, 2, 8}) =
            (-2 * wm1Sq * (1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) /
            (9. * mSqq2);
        result.element({0, 7, 2, 9}) =
            (2 * wm1Sq * (1 + w) * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (9. * mSqq2);
        result.element({0, 7, 3, 1}) = (-2 * rt * wm1Sq * (1 + w) * (-1 + rC * w)) / (3. * mSqq);
        result.element({0, 7, 4, 1}) = (2 * rt * wm1Sq * (1 + w) * (-1 + rC * w)) / (3. * mSqq);
        result.element({0, 7, 5, 10}) = (-8 * rC * rt * w2m1Sq) / (3. * mSqq);
        result.element({0, 7, 5, 13}) = (8 * rC * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({0, 8, 0, 6}) =
            (2 * wm1Sq * (1 + w) *
             (rtSq * (-1 + 3 * rCSq - 2 * rC * (-2 + w) - 4 * w) + 2 * mSqq * (-1 + rC) * (1 + w))) /
            (9. * mSqq2);
        result.element({0, 8, 0, 7}) =
            (2 * wm1Sq * (1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) /
            (9. * mSqq2);
        result.element({0, 8, 0, 8}) =
            (2 * wm1Sq * (1 + w) * (rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq) + 2 * mSqq * w2m1)) / (9. * mSqq2);
        result.element({0, 8, 0, 9}) = (-4 * (mSqq + 2 * rtSq) * (rC - w) * wm1Sq * (1 + w)) / (9. * mSqq2);
        result.element({0, 8, 1, 6}) =
            (2 * wm1Sq * (1 + w) *
             (rtSq * (-1 + 3 * rCSq - 2 * rC * (-2 + w) - 4 * w) + 2 * mSqq * (-1 + rC) * (1 + w))) /
            (9. * mSqq2);
        result.element({0, 8, 1, 7}) =
            (2 * wm1Sq * (1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) /
            (9. * mSqq2);
        result.element({0, 8, 1, 8}) =
            (2 * wm1Sq * (1 + w) * (rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq) + 2 * mSqq * w2m1)) / (9. * mSqq2);
        result.element({0, 8, 1, 9}) = (-4 * (mSqq + 2 * rtSq) * (rC - w) * wm1Sq * (1 + w)) / (9. * mSqq2);
        result.element({0, 8, 2, 6}) =
            (-2 * wm1Sq * (1 + w) *
             (rtSq * (-1 + 3 * rCSq - 2 * rC * (-2 + w) - 4 * w) + 2 * mSqq * (-1 + rC) * (1 + w))) /
            (9. * mSqq2);
        result.element({0, 8, 2, 7}) =
            (-2 * wm1Sq * (1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) /
            (9. * mSqq2);
        result.element({0, 8, 2, 8}) =
            (-2 * wm1Sq * (1 + w) * (rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq) + 2 * mSqq * w2m1)) / (9. * mSqq2);
        result.element({0, 8, 2, 9}) = (4 * (mSqq + 2 * rtSq) * (rC - w) * wm1Sq * (1 + w)) / (9. * mSqq2);
        result.element({0, 8, 3, 1}) = (-2 * rt * (rC - w) * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({0, 8, 4, 1}) = (2 * rt * (rC - w) * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({0, 8, 5, 10}) = (-8 * rt * w2m1Sq) / (3. * mSqq);
        result.element({0, 8, 5, 13}) = (8 * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({0, 9, 0, 6}) = (-2 *
                                        (2 * mSqq2 + rtSq * (-3 + rCSq + 2 * rC * (-2 + w) + 4 * w) +
                                         mSqq * (2 * rCSq + rtSq + 2 * w - 2 * rC * (1 + w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({0, 9, 0, 7}) =
            (-2 * wm1Sq * (1 + w) * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (9. * mSqq2);
        result.element({0, 9, 0, 8}) = (-4 * (mSqq + 2 * rtSq) * (rC - w) * wm1Sq * (1 + w)) / (9. * mSqq2);
        result.element({0, 9, 0, 9}) =
            (2 * (-1 + w) *
             (4 * mSqq2 + 2 * mSqq * (rCSq + rtSq - 2 * rC * w + wSq) + rtSq * (-3 + rCSq - 2 * rC * w + 4 * wSq))) /
            (9. * mSqq2);
        result.element({0, 9, 1, 6}) = (-2 *
                                        (2 * mSqq2 + rtSq * (-3 + rCSq + 2 * rC * (-2 + w) + 4 * w) +
                                         mSqq * (2 * rCSq + rtSq + 2 * w - 2 * rC * (1 + w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({0, 9, 1, 7}) =
            (-2 * wm1Sq * (1 + w) * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (9. * mSqq2);
        result.element({0, 9, 1, 8}) = (-4 * (mSqq + 2 * rtSq) * (rC - w) * wm1Sq * (1 + w)) / (9. * mSqq2);
        result.element({0, 9, 1, 9}) =
            (2 * (-1 + w) *
             (4 * mSqq2 + 2 * mSqq * (rCSq + rtSq - 2 * rC * w + wSq) + rtSq * (-3 + rCSq - 2 * rC * w + 4 * wSq))) /
            (9. * mSqq2);
        result.element({0, 9, 2, 6}) = (2 *
                                        (2 * mSqq2 + rtSq * (-3 + rCSq + 2 * rC * (-2 + w) + 4 * w) +
                                         mSqq * (2 * rCSq + rtSq + 2 * w - 2 * rC * (1 + w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({0, 9, 2, 7}) =
            (2 * wm1Sq * (1 + w) * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (9. * mSqq2);
        result.element({0, 9, 2, 8}) = (4 * (mSqq + 2 * rtSq) * (rC - w) * wm1Sq * (1 + w)) / (9. * mSqq2);
        result.element({0, 9, 2, 9}) =
            (-2 * (-1 + w) *
             (4 * mSqq2 + 2 * mSqq * (rCSq + rtSq - 2 * rC * w + wSq) + rtSq * (-3 + rCSq - 2 * rC * w + 4 * wSq))) /
            (9. * mSqq2);
        result.element({0, 9, 3, 1}) = (2 * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({0, 9, 4, 1}) = (-2 * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({0, 9, 5, 10}) = (-8 * rt * (1 - 2 * rC + w) * w2m1) / (3. * mSqq);
        result.element({0, 9, 5, 11}) = (8 * rC * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({0, 9, 5, 12}) = (8 * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({0, 9, 5, 13}) = (-8 * rt * (1 + 3 * rC - 2 * w) * (-1 + w)) / (3. * mSqq);
        result.element({0, 9, 5, 14}) = (-8 * rt * w2m1) / (3. * mSqq);
        result.element({1, 2, 0, 2}) = (4 *
                                        (2 * mSqq2 * (-1 + w) + mSqq * (1 + 2 * rC + rCSq + rtSq) * (-1 + w) +
                                         rtSq * (1 + 2 * w - 2 * rC * (2 + w) + rCSq * (1 + 2 * w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({1, 2, 0, 3}) =
            (2 * (-1 + w) * wp1Sq *
             (2 * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2 * rC * (2 + w) + rCSq * (-1 + 4 * w)))) /
            (9. * mSqq2);
        result.element({1, 2, 0, 4}) =
            (2 * (-1 + w) * wp1Sq *
             (2 * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3 * rCSq + 4 * w - 2 * rC * (2 + w)))) /
            (9. * mSqq2);
        result.element({1, 2, 0, 5}) = (-2 *
                                        (2 * mSqq2 + mSqq * (2 * rCSq + rtSq - 2 * rC * (-1 + w) - 2 * w) +
                                         rtSq * (-3 + rCSq - 4 * w + 2 * rC * (2 + w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({1, 2, 1, 2}) = (4 *
                                        (2 * mSqq2 * (-1 + w) + mSqq * (1 + 2 * rC + rCSq + rtSq) * (-1 + w) +
                                         rtSq * (1 + 2 * w - 2 * rC * (2 + w) + rCSq * (1 + 2 * w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({1, 2, 1, 3}) =
            (2 * (-1 + w) * wp1Sq *
             (2 * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2 * rC * (2 + w) + rCSq * (-1 + 4 * w)))) /
            (9. * mSqq2);
        result.element({1, 2, 1, 4}) =
            (2 * (-1 + w) * wp1Sq *
             (2 * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3 * rCSq + 4 * w - 2 * rC * (2 + w)))) /
            (9. * mSqq2);
        result.element({1, 2, 1, 5}) = (-2 *
                                        (2 * mSqq2 + mSqq * (2 * rCSq + rtSq - 2 * rC * (-1 + w) - 2 * w) +
                                         rtSq * (-3 + rCSq - 4 * w + 2 * rC * (2 + w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({1, 2, 2, 2}) = (4 *
                                        (2 * mSqq2 * (-1 + w) + mSqq * (1 + 2 * rC + rCSq + rtSq) * (-1 + w) +
                                         rtSq * (1 + 2 * w - 2 * rC * (2 + w) + rCSq * (1 + 2 * w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({1, 2, 2, 3}) =
            (2 * (-1 + w) * wp1Sq *
             (2 * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2 * rC * (2 + w) + rCSq * (-1 + 4 * w)))) /
            (9. * mSqq2);
        result.element({1, 2, 2, 4}) =
            (2 * (-1 + w) * wp1Sq *
             (2 * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3 * rCSq + 4 * w - 2 * rC * (2 + w)))) /
            (9. * mSqq2);
        result.element({1, 2, 2, 5}) = (-2 *
                                        (2 * mSqq2 + mSqq * (2 * rCSq + rtSq - 2 * rC * (-1 + w) - 2 * w) +
                                         rtSq * (-3 + rCSq - 4 * w + 2 * rC * (2 + w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({1, 2, 3, 0}) = (2 * (-1 + rC) * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({1, 2, 4, 0}) = (2 * (-1 + rC) * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({1, 2, 5, 10}) = (-8 * (1 + rC) * rt * wm1Sq * (1 + w)) / mSqq;
        result.element({1, 2, 5, 11}) = (-8 * rt * wm1Sq * (1 + w) * (-3 + rC * (-1 + 2 * w))) / (3. * mSqq);
        result.element({1, 2, 5, 12}) = (-8 * rt * (1 + 3 * rC - 2 * w) * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({1, 2, 5, 13}) = (-16 * rt * (1 + rC - w) * w2m1) / (3. * mSqq);
        result.element({1, 2, 5, 14}) = (-8 * rt * (1 + 2 * rC - w) * w2m1) / (3. * mSqq);
        result.element({1, 2, 5, 15}) = (8 * (1 + rC) * rt * w2m1Sq) / (3. * mSqq);
        result.element({1, 3, 0, 2}) =
            (2 * (-1 + w) * wp1Sq *
             (2 * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2 * rC * (2 + w) + rCSq * (-1 + 4 * w)))) /
            (9. * mSqq2);
        result.element({1, 3, 0, 3}) =
            (2 * (-1 + w) * wp1Sq * (rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)) + 2 * mSqq * rCSq * w2m1)) /
            (9. * mSqq2);
        result.element({1, 3, 0, 4}) =
            (2 * (-1 + w) * wp1Sq * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) /
            (9. * mSqq2);
        result.element({1, 3, 0, 5}) =
            (-2 * (-1 + w) * wp1Sq * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (9. * mSqq2);
        result.element({1, 3, 1, 2}) =
            (2 * (-1 + w) * wp1Sq *
             (2 * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2 * rC * (2 + w) + rCSq * (-1 + 4 * w)))) /
            (9. * mSqq2);
        result.element({1, 3, 1, 3}) =
            (2 * (-1 + w) * wp1Sq * (rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)) + 2 * mSqq * rCSq * w2m1)) /
            (9. * mSqq2);
        result.element({1, 3, 1, 4}) =
            (2 * (-1 + w) * wp1Sq * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) /
            (9. * mSqq2);
        result.element({1, 3, 1, 5}) =
            (-2 * (-1 + w) * wp1Sq * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (9. * mSqq2);
        result.element({1, 3, 2, 2}) =
            (2 * (-1 + w) * wp1Sq *
             (2 * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2 * rC * (2 + w) + rCSq * (-1 + 4 * w)))) /
            (9. * mSqq2);
        result.element({1, 3, 2, 3}) =
            (2 * (-1 + w) * wp1Sq * (rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)) + 2 * mSqq * rCSq * w2m1)) /
            (9. * mSqq2);
        result.element({1, 3, 2, 4}) =
            (2 * (-1 + w) * wp1Sq * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) /
            (9. * mSqq2);
        result.element({1, 3, 2, 5}) =
            (-2 * (-1 + w) * wp1Sq * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (9. * mSqq2);
        result.element({1, 3, 3, 0}) = (2 * rt * (-1 + w) * wp1Sq * (-1 + rC * w)) / (3. * mSqq);
        result.element({1, 3, 4, 0}) = (2 * rt * (-1 + w) * wp1Sq * (-1 + rC * w)) / (3. * mSqq);
        result.element({1, 3, 5, 10}) = (-8 * rC * rt * w2m1Sq) / (3. * mSqq);
        result.element({1, 3, 5, 11}) = (8 * rC * rt * w2m1Sq) / (3. * mSqq);
        result.element({1, 3, 5, 12}) = (-8 * rC * rt * w2m1Sq) / (3. * mSqq);
        result.element({1, 3, 5, 13}) = (-8 * rC * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({1, 3, 5, 14}) = (-8 * rC * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({1, 3, 5, 15}) = (8 * rC * rt * wm1Sq * wp1Cu) / (3. * mSqq);
        result.element({1, 4, 0, 2}) =
            (2 * (-1 + w) * wp1Sq *
             (2 * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3 * rCSq + 4 * w - 2 * rC * (2 + w)))) /
            (9. * mSqq2);
        result.element({1, 4, 0, 3}) =
            (2 * (-1 + w) * wp1Sq * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) /
            (9. * mSqq2);
        result.element({1, 4, 0, 4}) =
            (2 * (-1 + w) * wp1Sq * (rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq) + 2 * mSqq * w2m1)) / (9. * mSqq2);
        result.element({1, 4, 0, 5}) = (-4 * (mSqq + 2 * rtSq) * (rC - w) * (-1 + w) * wp1Sq) / (9. * mSqq2);
        result.element({1, 4, 1, 2}) =
            (2 * (-1 + w) * wp1Sq *
             (2 * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3 * rCSq + 4 * w - 2 * rC * (2 + w)))) /
            (9. * mSqq2);
        result.element({1, 4, 1, 3}) =
            (2 * (-1 + w) * wp1Sq * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) /
            (9. * mSqq2);
        result.element({1, 4, 1, 4}) =
            (2 * (-1 + w) * wp1Sq * (rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq) + 2 * mSqq * w2m1)) / (9. * mSqq2);
        result.element({1, 4, 1, 5}) = (-4 * (mSqq + 2 * rtSq) * (rC - w) * (-1 + w) * wp1Sq) / (9. * mSqq2);
        result.element({1, 4, 2, 2}) =
            (2 * (-1 + w) * wp1Sq *
             (2 * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3 * rCSq + 4 * w - 2 * rC * (2 + w)))) /
            (9. * mSqq2);
        result.element({1, 4, 2, 3}) =
            (2 * (-1 + w) * wp1Sq * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) /
            (9. * mSqq2);
        result.element({1, 4, 2, 4}) =
            (2 * (-1 + w) * wp1Sq * (rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq) + 2 * mSqq * w2m1)) / (9. * mSqq2);
        result.element({1, 4, 2, 5}) = (-4 * (mSqq + 2 * rtSq) * (rC - w) * (-1 + w) * wp1Sq) / (9. * mSqq2);
        result.element({1, 4, 3, 0}) = (2 * rt * (rC - w) * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({1, 4, 4, 0}) = (2 * rt * (rC - w) * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({1, 4, 5, 10}) = (-8 * rt * w2m1Sq) / (3. * mSqq);
        result.element({1, 4, 5, 11}) = (8 * rt * w2m1Sq) / (3. * mSqq);
        result.element({1, 4, 5, 12}) = (-8 * rt * w2m1Sq) / (3. * mSqq);
        result.element({1, 4, 5, 13}) = (-8 * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({1, 4, 5, 14}) = (-8 * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({1, 4, 5, 15}) = (8 * rt * wm1Sq * wp1Cu) / (3. * mSqq);
        result.element({1, 5, 0, 2}) = (-2 *
                                        (2 * mSqq2 + mSqq * (2 * rCSq + rtSq - 2 * rC * (-1 + w) - 2 * w) +
                                         rtSq * (-3 + rCSq - 4 * w + 2 * rC * (2 + w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({1, 5, 0, 3}) =
            (-2 * (-1 + w) * wp1Sq * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (9. * mSqq2);
        result.element({1, 5, 0, 4}) = (-4 * (mSqq + 2 * rtSq) * (rC - w) * (-1 + w) * wp1Sq) / (9. * mSqq2);
        result.element({1, 5, 0, 5}) =
            (2 * (1 + w) *
             (4 * mSqq2 + 2 * mSqq * (rCSq + rtSq - 2 * rC * w + wSq) + rtSq * (-3 + rCSq - 2 * rC * w + 4 * wSq))) /
            (9. * mSqq2);
        result.element({1, 5, 1, 2}) = (-2 *
                                        (2 * mSqq2 + mSqq * (2 * rCSq + rtSq - 2 * rC * (-1 + w) - 2 * w) +
                                         rtSq * (-3 + rCSq - 4 * w + 2 * rC * (2 + w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({1, 5, 1, 3}) =
            (-2 * (-1 + w) * wp1Sq * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (9. * mSqq2);
        result.element({1, 5, 1, 4}) = (-4 * (mSqq + 2 * rtSq) * (rC - w) * (-1 + w) * wp1Sq) / (9. * mSqq2);
        result.element({1, 5, 1, 5}) =
            (2 * (1 + w) *
             (4 * mSqq2 + 2 * mSqq * (rCSq + rtSq - 2 * rC * w + wSq) + rtSq * (-3 + rCSq - 2 * rC * w + 4 * wSq))) /
            (9. * mSqq2);
        result.element({1, 5, 2, 2}) = (-2 *
                                        (2 * mSqq2 + mSqq * (2 * rCSq + rtSq - 2 * rC * (-1 + w) - 2 * w) +
                                         rtSq * (-3 + rCSq - 4 * w + 2 * rC * (2 + w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({1, 5, 2, 3}) =
            (-2 * (-1 + w) * wp1Sq * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (9. * mSqq2);
        result.element({1, 5, 2, 4}) = (-4 * (mSqq + 2 * rtSq) * (rC - w) * (-1 + w) * wp1Sq) / (9. * mSqq2);
        result.element({1, 5, 2, 5}) =
            (2 * (1 + w) *
             (4 * mSqq2 + 2 * mSqq * (rCSq + rtSq - 2 * rC * w + wSq) + rtSq * (-3 + rCSq - 2 * rC * w + 4 * wSq))) /
            (9. * mSqq2);
        result.element({1, 5, 3, 0}) = (-2 * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({1, 5, 4, 0}) = (-2 * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({1, 5, 5, 10}) = (8 * rt * (1 + 2 * rC - w) * w2m1) / (3. * mSqq);
        result.element({1, 5, 5, 11}) = (8 * (1 + rC) * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({1, 5, 5, 12}) = (16 * rt * (rC - w) * w2m1) / (3. * mSqq);
        result.element({1, 5, 5, 13}) = (8 * rt * (-1 + 3 * rC - 2 * w) * (1 + w)) / (3. * mSqq);
        result.element({1, 5, 5, 14}) = (8 * rt * (rC - w) * (1 + w)) / mSqq;
        result.element({1, 5, 5, 15}) = (-8 * rt * (rC - w) * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({1, 6, 0, 6}) = (4 *
                                        (2 * mSqq2 * (1 + w) + mSqq * (1 - 2 * rC + rCSq + rtSq) * (1 + w) +
                                         rtSq * (-1 + 2 * rC * (-2 + w) + 2 * w + rCSq * (-1 + 2 * w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({1, 6, 0, 7}) =
            (2 * wm1Sq * (1 + w) *
             (2 * mSqq * (-1 + rC) * rC * (1 + w) + rtSq * (-3 + 2 * rC * (-2 + w) + rCSq * (1 + 4 * w)))) /
            (9. * mSqq2);
        result.element({1, 6, 0, 8}) =
            (2 * wm1Sq * (1 + w) *
             (rtSq * (-1 + 3 * rCSq - 2 * rC * (-2 + w) - 4 * w) + 2 * mSqq * (-1 + rC) * (1 + w))) /
            (9. * mSqq2);
        result.element({1, 6, 0, 9}) = (-2 *
                                        (2 * mSqq2 + rtSq * (-3 + rCSq + 2 * rC * (-2 + w) + 4 * w) +
                                         mSqq * (2 * rCSq + rtSq + 2 * w - 2 * rC * (1 + w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({1, 6, 1, 6}) = (4 *
                                        (2 * mSqq2 * (1 + w) + mSqq * (1 - 2 * rC + rCSq + rtSq) * (1 + w) +
                                         rtSq * (-1 + 2 * rC * (-2 + w) + 2 * w + rCSq * (-1 + 2 * w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({1, 6, 1, 7}) =
            (2 * wm1Sq * (1 + w) *
             (2 * mSqq * (-1 + rC) * rC * (1 + w) + rtSq * (-3 + 2 * rC * (-2 + w) + rCSq * (1 + 4 * w)))) /
            (9. * mSqq2);
        result.element({1, 6, 1, 8}) =
            (2 * wm1Sq * (1 + w) *
             (rtSq * (-1 + 3 * rCSq - 2 * rC * (-2 + w) - 4 * w) + 2 * mSqq * (-1 + rC) * (1 + w))) /
            (9. * mSqq2);
        result.element({1, 6, 1, 9}) = (-2 *
                                        (2 * mSqq2 + rtSq * (-3 + rCSq + 2 * rC * (-2 + w) + 4 * w) +
                                         mSqq * (2 * rCSq + rtSq + 2 * w - 2 * rC * (1 + w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({1, 6, 2, 6}) = (-4 *
                                        (2 * mSqq2 * (1 + w) + mSqq * (1 - 2 * rC + rCSq + rtSq) * (1 + w) +
                                         rtSq * (-1 + 2 * rC * (-2 + w) + 2 * w + rCSq * (-1 + 2 * w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({1, 6, 2, 7}) =
            (-2 * wm1Sq * (1 + w) *
             (2 * mSqq * (-1 + rC) * rC * (1 + w) + rtSq * (-3 + 2 * rC * (-2 + w) + rCSq * (1 + 4 * w)))) /
            (9. * mSqq2);
        result.element({1, 6, 2, 8}) =
            (-2 * wm1Sq * (1 + w) *
             (rtSq * (-1 + 3 * rCSq - 2 * rC * (-2 + w) - 4 * w) + 2 * mSqq * (-1 + rC) * (1 + w))) /
            (9. * mSqq2);
        result.element({1, 6, 2, 9}) = (2 *
                                        (2 * mSqq2 + rtSq * (-3 + rCSq + 2 * rC * (-2 + w) + 4 * w) +
                                         mSqq * (2 * rCSq + rtSq + 2 * w - 2 * rC * (1 + w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({1, 6, 3, 1}) = (-2 * (1 + rC) * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({1, 6, 4, 1}) = (2 * (1 + rC) * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({1, 6, 5, 10}) = (-8 * (-1 + rC) * rt * (-1 + w) * wp1Sq) / mSqq;
        result.element({1, 6, 5, 11}) = (-16 * rC * rt * w2m1Sq) / (3. * mSqq);
        result.element({1, 6, 5, 12}) = (-16 * rt * w2m1Sq) / (3. * mSqq);
        result.element({1, 6, 5, 13}) = (16 * rt * (-1 + rC - w) * w2m1) / (3. * mSqq);
        result.element({1, 6, 5, 14}) = (-8 * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({1, 7, 0, 6}) =
            (2 * wm1Sq * (1 + w) *
             (2 * mSqq * (-1 + rC) * rC * (1 + w) + rtSq * (-3 + 2 * rC * (-2 + w) + rCSq * (1 + 4 * w)))) /
            (9. * mSqq2);
        result.element({1, 7, 0, 7}) =
            (2 * wm1Sq * (1 + w) * (rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)) + 2 * mSqq * rCSq * w2m1)) /
            (9. * mSqq2);
        result.element({1, 7, 0, 8}) =
            (2 * wm1Sq * (1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) /
            (9. * mSqq2);
        result.element({1, 7, 0, 9}) =
            (-2 * wm1Sq * (1 + w) * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (9. * mSqq2);
        result.element({1, 7, 1, 6}) =
            (2 * wm1Sq * (1 + w) *
             (2 * mSqq * (-1 + rC) * rC * (1 + w) + rtSq * (-3 + 2 * rC * (-2 + w) + rCSq * (1 + 4 * w)))) /
            (9. * mSqq2);
        result.element({1, 7, 1, 7}) =
            (2 * wm1Sq * (1 + w) * (rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)) + 2 * mSqq * rCSq * w2m1)) /
            (9. * mSqq2);
        result.element({1, 7, 1, 8}) =
            (2 * wm1Sq * (1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) /
            (9. * mSqq2);
        result.element({1, 7, 1, 9}) =
            (-2 * wm1Sq * (1 + w) * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (9. * mSqq2);
        result.element({1, 7, 2, 6}) =
            (-2 * wm1Sq * (1 + w) *
             (2 * mSqq * (-1 + rC) * rC * (1 + w) + rtSq * (-3 + 2 * rC * (-2 + w) + rCSq * (1 + 4 * w)))) /
            (9. * mSqq2);
        result.element({1, 7, 2, 7}) =
            (-2 * wm1Sq * (1 + w) * (rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)) + 2 * mSqq * rCSq * w2m1)) /
            (9. * mSqq2);
        result.element({1, 7, 2, 8}) =
            (-2 * wm1Sq * (1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) /
            (9. * mSqq2);
        result.element({1, 7, 2, 9}) =
            (2 * wm1Sq * (1 + w) * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (9. * mSqq2);
        result.element({1, 7, 3, 1}) = (-2 * rt * wm1Sq * (1 + w) * (-1 + rC * w)) / (3. * mSqq);
        result.element({1, 7, 4, 1}) = (2 * rt * wm1Sq * (1 + w) * (-1 + rC * w)) / (3. * mSqq);
        result.element({1, 7, 5, 10}) = (-8 * rC * rt * w2m1Sq) / (3. * mSqq);
        result.element({1, 7, 5, 13}) = (8 * rC * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({1, 8, 0, 6}) =
            (2 * wm1Sq * (1 + w) *
             (rtSq * (-1 + 3 * rCSq - 2 * rC * (-2 + w) - 4 * w) + 2 * mSqq * (-1 + rC) * (1 + w))) /
            (9. * mSqq2);
        result.element({1, 8, 0, 7}) =
            (2 * wm1Sq * (1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) /
            (9. * mSqq2);
        result.element({1, 8, 0, 8}) =
            (2 * wm1Sq * (1 + w) * (rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq) + 2 * mSqq * w2m1)) / (9. * mSqq2);
        result.element({1, 8, 0, 9}) = (-4 * (mSqq + 2 * rtSq) * (rC - w) * wm1Sq * (1 + w)) / (9. * mSqq2);
        result.element({1, 8, 1, 6}) =
            (2 * wm1Sq * (1 + w) *
             (rtSq * (-1 + 3 * rCSq - 2 * rC * (-2 + w) - 4 * w) + 2 * mSqq * (-1 + rC) * (1 + w))) /
            (9. * mSqq2);
        result.element({1, 8, 1, 7}) =
            (2 * wm1Sq * (1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) /
            (9. * mSqq2);
        result.element({1, 8, 1, 8}) =
            (2 * wm1Sq * (1 + w) * (rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq) + 2 * mSqq * w2m1)) / (9. * mSqq2);
        result.element({1, 8, 1, 9}) = (-4 * (mSqq + 2 * rtSq) * (rC - w) * wm1Sq * (1 + w)) / (9. * mSqq2);
        result.element({1, 8, 2, 6}) =
            (-2 * wm1Sq * (1 + w) *
             (rtSq * (-1 + 3 * rCSq - 2 * rC * (-2 + w) - 4 * w) + 2 * mSqq * (-1 + rC) * (1 + w))) /
            (9. * mSqq2);
        result.element({1, 8, 2, 7}) =
            (-2 * wm1Sq * (1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) /
            (9. * mSqq2);
        result.element({1, 8, 2, 8}) =
            (-2 * wm1Sq * (1 + w) * (rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq) + 2 * mSqq * w2m1)) / (9. * mSqq2);
        result.element({1, 8, 2, 9}) = (4 * (mSqq + 2 * rtSq) * (rC - w) * wm1Sq * (1 + w)) / (9. * mSqq2);
        result.element({1, 8, 3, 1}) = (-2 * rt * (rC - w) * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({1, 8, 4, 1}) = (2 * rt * (rC - w) * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({1, 8, 5, 10}) = (-8 * rt * w2m1Sq) / (3. * mSqq);
        result.element({1, 8, 5, 13}) = (8 * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({1, 9, 0, 6}) = (-2 *
                                        (2 * mSqq2 + rtSq * (-3 + rCSq + 2 * rC * (-2 + w) + 4 * w) +
                                         mSqq * (2 * rCSq + rtSq + 2 * w - 2 * rC * (1 + w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({1, 9, 0, 7}) =
            (-2 * wm1Sq * (1 + w) * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (9. * mSqq2);
        result.element({1, 9, 0, 8}) = (-4 * (mSqq + 2 * rtSq) * (rC - w) * wm1Sq * (1 + w)) / (9. * mSqq2);
        result.element({1, 9, 0, 9}) =
            (2 * (-1 + w) *
             (4 * mSqq2 + 2 * mSqq * (rCSq + rtSq - 2 * rC * w + wSq) + rtSq * (-3 + rCSq - 2 * rC * w + 4 * wSq))) /
            (9. * mSqq2);
        result.element({1, 9, 1, 6}) = (-2 *
                                        (2 * mSqq2 + rtSq * (-3 + rCSq + 2 * rC * (-2 + w) + 4 * w) +
                                         mSqq * (2 * rCSq + rtSq + 2 * w - 2 * rC * (1 + w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({1, 9, 1, 7}) =
            (-2 * wm1Sq * (1 + w) * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (9. * mSqq2);
        result.element({1, 9, 1, 8}) = (-4 * (mSqq + 2 * rtSq) * (rC - w) * wm1Sq * (1 + w)) / (9. * mSqq2);
        result.element({1, 9, 1, 9}) =
            (2 * (-1 + w) *
             (4 * mSqq2 + 2 * mSqq * (rCSq + rtSq - 2 * rC * w + wSq) + rtSq * (-3 + rCSq - 2 * rC * w + 4 * wSq))) /
            (9. * mSqq2);
        result.element({1, 9, 2, 6}) = (2 *
                                        (2 * mSqq2 + rtSq * (-3 + rCSq + 2 * rC * (-2 + w) + 4 * w) +
                                         mSqq * (2 * rCSq + rtSq + 2 * w - 2 * rC * (1 + w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({1, 9, 2, 7}) =
            (2 * wm1Sq * (1 + w) * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (9. * mSqq2);
        result.element({1, 9, 2, 8}) = (4 * (mSqq + 2 * rtSq) * (rC - w) * wm1Sq * (1 + w)) / (9. * mSqq2);
        result.element({1, 9, 2, 9}) =
            (-2 * (-1 + w) *
             (4 * mSqq2 + 2 * mSqq * (rCSq + rtSq - 2 * rC * w + wSq) + rtSq * (-3 + rCSq - 2 * rC * w + 4 * wSq))) /
            (9. * mSqq2);
        result.element({1, 9, 3, 1}) = (2 * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({1, 9, 4, 1}) = (-2 * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({1, 9, 5, 10}) = (-8 * rt * (1 - 2 * rC + w) * w2m1) / (3. * mSqq);
        result.element({1, 9, 5, 11}) = (8 * rC * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({1, 9, 5, 12}) = (8 * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({1, 9, 5, 13}) = (-8 * rt * (1 + 3 * rC - 2 * w) * (-1 + w)) / (3. * mSqq);
        result.element({1, 9, 5, 14}) = (-8 * rt * w2m1) / (3. * mSqq);
        result.element({2, 2, 0, 2}) = (4 *
                                        (2 * mSqq2 * (-1 + w) + mSqq * (1 + 2 * rC + rCSq + rtSq) * (-1 + w) +
                                         rtSq * (1 + 2 * w - 2 * rC * (2 + w) + rCSq * (1 + 2 * w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({2, 2, 0, 3}) =
            (2 * (-1 + w) * wp1Sq *
             (2 * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2 * rC * (2 + w) + rCSq * (-1 + 4 * w)))) /
            (9. * mSqq2);
        result.element({2, 2, 0, 4}) =
            (2 * (-1 + w) * wp1Sq *
             (2 * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3 * rCSq + 4 * w - 2 * rC * (2 + w)))) /
            (9. * mSqq2);
        result.element({2, 2, 0, 5}) = (-2 *
                                        (2 * mSqq2 + mSqq * (2 * rCSq + rtSq - 2 * rC * (-1 + w) - 2 * w) +
                                         rtSq * (-3 + rCSq - 4 * w + 2 * rC * (2 + w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({2, 2, 1, 2}) = (4 *
                                        (2 * mSqq2 * (-1 + w) + mSqq * (1 + 2 * rC + rCSq + rtSq) * (-1 + w) +
                                         rtSq * (1 + 2 * w - 2 * rC * (2 + w) + rCSq * (1 + 2 * w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({2, 2, 1, 3}) =
            (2 * (-1 + w) * wp1Sq *
             (2 * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2 * rC * (2 + w) + rCSq * (-1 + 4 * w)))) /
            (9. * mSqq2);
        result.element({2, 2, 1, 4}) =
            (2 * (-1 + w) * wp1Sq *
             (2 * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3 * rCSq + 4 * w - 2 * rC * (2 + w)))) /
            (9. * mSqq2);
        result.element({2, 2, 1, 5}) = (-2 *
                                        (2 * mSqq2 + mSqq * (2 * rCSq + rtSq - 2 * rC * (-1 + w) - 2 * w) +
                                         rtSq * (-3 + rCSq - 4 * w + 2 * rC * (2 + w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({2, 2, 2, 2}) = (4 *
                                        (2 * mSqq2 * (-1 + w) + mSqq * (1 + 2 * rC + rCSq + rtSq) * (-1 + w) +
                                         rtSq * (1 + 2 * w - 2 * rC * (2 + w) + rCSq * (1 + 2 * w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({2, 2, 2, 3}) =
            (2 * (-1 + w) * wp1Sq *
             (2 * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2 * rC * (2 + w) + rCSq * (-1 + 4 * w)))) /
            (9. * mSqq2);
        result.element({2, 2, 2, 4}) =
            (2 * (-1 + w) * wp1Sq *
             (2 * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3 * rCSq + 4 * w - 2 * rC * (2 + w)))) /
            (9. * mSqq2);
        result.element({2, 2, 2, 5}) = (-2 *
                                        (2 * mSqq2 + mSqq * (2 * rCSq + rtSq - 2 * rC * (-1 + w) - 2 * w) +
                                         rtSq * (-3 + rCSq - 4 * w + 2 * rC * (2 + w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({2, 2, 3, 0}) = (2 * (-1 + rC) * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({2, 2, 4, 0}) = (2 * (-1 + rC) * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({2, 2, 5, 10}) = (-8 * (1 + rC) * rt * wm1Sq * (1 + w)) / mSqq;
        result.element({2, 2, 5, 11}) = (-8 * rt * wm1Sq * (1 + w) * (-3 + rC * (-1 + 2 * w))) / (3. * mSqq);
        result.element({2, 2, 5, 12}) = (-8 * rt * (1 + 3 * rC - 2 * w) * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({2, 2, 5, 13}) = (-16 * rt * (1 + rC - w) * w2m1) / (3. * mSqq);
        result.element({2, 2, 5, 14}) = (-8 * rt * (1 + 2 * rC - w) * w2m1) / (3. * mSqq);
        result.element({2, 2, 5, 15}) = (8 * (1 + rC) * rt * w2m1Sq) / (3. * mSqq);
        result.element({2, 3, 0, 2}) =
            (2 * (-1 + w) * wp1Sq *
             (2 * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2 * rC * (2 + w) + rCSq * (-1 + 4 * w)))) /
            (9. * mSqq2);
        result.element({2, 3, 0, 3}) =
            (2 * (-1 + w) * wp1Sq * (rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)) + 2 * mSqq * rCSq * w2m1)) /
            (9. * mSqq2);
        result.element({2, 3, 0, 4}) =
            (2 * (-1 + w) * wp1Sq * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) /
            (9. * mSqq2);
        result.element({2, 3, 0, 5}) =
            (-2 * (-1 + w) * wp1Sq * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (9. * mSqq2);
        result.element({2, 3, 1, 2}) =
            (2 * (-1 + w) * wp1Sq *
             (2 * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2 * rC * (2 + w) + rCSq * (-1 + 4 * w)))) /
            (9. * mSqq2);
        result.element({2, 3, 1, 3}) =
            (2 * (-1 + w) * wp1Sq * (rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)) + 2 * mSqq * rCSq * w2m1)) /
            (9. * mSqq2);
        result.element({2, 3, 1, 4}) =
            (2 * (-1 + w) * wp1Sq * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) /
            (9. * mSqq2);
        result.element({2, 3, 1, 5}) =
            (-2 * (-1 + w) * wp1Sq * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (9. * mSqq2);
        result.element({2, 3, 2, 2}) =
            (2 * (-1 + w) * wp1Sq *
             (2 * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2 * rC * (2 + w) + rCSq * (-1 + 4 * w)))) /
            (9. * mSqq2);
        result.element({2, 3, 2, 3}) =
            (2 * (-1 + w) * wp1Sq * (rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)) + 2 * mSqq * rCSq * w2m1)) /
            (9. * mSqq2);
        result.element({2, 3, 2, 4}) =
            (2 * (-1 + w) * wp1Sq * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) /
            (9. * mSqq2);
        result.element({2, 3, 2, 5}) =
            (-2 * (-1 + w) * wp1Sq * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (9. * mSqq2);
        result.element({2, 3, 3, 0}) = (2 * rt * (-1 + w) * wp1Sq * (-1 + rC * w)) / (3. * mSqq);
        result.element({2, 3, 4, 0}) = (2 * rt * (-1 + w) * wp1Sq * (-1 + rC * w)) / (3. * mSqq);
        result.element({2, 3, 5, 10}) = (-8 * rC * rt * w2m1Sq) / (3. * mSqq);
        result.element({2, 3, 5, 11}) = (8 * rC * rt * w2m1Sq) / (3. * mSqq);
        result.element({2, 3, 5, 12}) = (-8 * rC * rt * w2m1Sq) / (3. * mSqq);
        result.element({2, 3, 5, 13}) = (-8 * rC * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({2, 3, 5, 14}) = (-8 * rC * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({2, 3, 5, 15}) = (8 * rC * rt * wm1Sq * wp1Cu) / (3. * mSqq);
        result.element({2, 4, 0, 2}) =
            (2 * (-1 + w) * wp1Sq *
             (2 * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3 * rCSq + 4 * w - 2 * rC * (2 + w)))) /
            (9. * mSqq2);
        result.element({2, 4, 0, 3}) =
            (2 * (-1 + w) * wp1Sq * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) /
            (9. * mSqq2);
        result.element({2, 4, 0, 4}) =
            (2 * (-1 + w) * wp1Sq * (rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq) + 2 * mSqq * w2m1)) / (9. * mSqq2);
        result.element({2, 4, 0, 5}) = (-4 * (mSqq + 2 * rtSq) * (rC - w) * (-1 + w) * wp1Sq) / (9. * mSqq2);
        result.element({2, 4, 1, 2}) =
            (2 * (-1 + w) * wp1Sq *
             (2 * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3 * rCSq + 4 * w - 2 * rC * (2 + w)))) /
            (9. * mSqq2);
        result.element({2, 4, 1, 3}) =
            (2 * (-1 + w) * wp1Sq * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) /
            (9. * mSqq2);
        result.element({2, 4, 1, 4}) =
            (2 * (-1 + w) * wp1Sq * (rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq) + 2 * mSqq * w2m1)) / (9. * mSqq2);
        result.element({2, 4, 1, 5}) = (-4 * (mSqq + 2 * rtSq) * (rC - w) * (-1 + w) * wp1Sq) / (9. * mSqq2);
        result.element({2, 4, 2, 2}) =
            (2 * (-1 + w) * wp1Sq *
             (2 * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3 * rCSq + 4 * w - 2 * rC * (2 + w)))) /
            (9. * mSqq2);
        result.element({2, 4, 2, 3}) =
            (2 * (-1 + w) * wp1Sq * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) /
            (9. * mSqq2);
        result.element({2, 4, 2, 4}) =
            (2 * (-1 + w) * wp1Sq * (rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq) + 2 * mSqq * w2m1)) / (9. * mSqq2);
        result.element({2, 4, 2, 5}) = (-4 * (mSqq + 2 * rtSq) * (rC - w) * (-1 + w) * wp1Sq) / (9. * mSqq2);
        result.element({2, 4, 3, 0}) = (2 * rt * (rC - w) * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({2, 4, 4, 0}) = (2 * rt * (rC - w) * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({2, 4, 5, 10}) = (-8 * rt * w2m1Sq) / (3. * mSqq);
        result.element({2, 4, 5, 11}) = (8 * rt * w2m1Sq) / (3. * mSqq);
        result.element({2, 4, 5, 12}) = (-8 * rt * w2m1Sq) / (3. * mSqq);
        result.element({2, 4, 5, 13}) = (-8 * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({2, 4, 5, 14}) = (-8 * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({2, 4, 5, 15}) = (8 * rt * wm1Sq * wp1Cu) / (3. * mSqq);
        result.element({2, 5, 0, 2}) = (-2 *
                                        (2 * mSqq2 + mSqq * (2 * rCSq + rtSq - 2 * rC * (-1 + w) - 2 * w) +
                                         rtSq * (-3 + rCSq - 4 * w + 2 * rC * (2 + w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({2, 5, 0, 3}) =
            (-2 * (-1 + w) * wp1Sq * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (9. * mSqq2);
        result.element({2, 5, 0, 4}) = (-4 * (mSqq + 2 * rtSq) * (rC - w) * (-1 + w) * wp1Sq) / (9. * mSqq2);
        result.element({2, 5, 0, 5}) =
            (2 * (1 + w) *
             (4 * mSqq2 + 2 * mSqq * (rCSq + rtSq - 2 * rC * w + wSq) + rtSq * (-3 + rCSq - 2 * rC * w + 4 * wSq))) /
            (9. * mSqq2);
        result.element({2, 5, 1, 2}) = (-2 *
                                        (2 * mSqq2 + mSqq * (2 * rCSq + rtSq - 2 * rC * (-1 + w) - 2 * w) +
                                         rtSq * (-3 + rCSq - 4 * w + 2 * rC * (2 + w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({2, 5, 1, 3}) =
            (-2 * (-1 + w) * wp1Sq * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (9. * mSqq2);
        result.element({2, 5, 1, 4}) = (-4 * (mSqq + 2 * rtSq) * (rC - w) * (-1 + w) * wp1Sq) / (9. * mSqq2);
        result.element({2, 5, 1, 5}) =
            (2 * (1 + w) *
             (4 * mSqq2 + 2 * mSqq * (rCSq + rtSq - 2 * rC * w + wSq) + rtSq * (-3 + rCSq - 2 * rC * w + 4 * wSq))) /
            (9. * mSqq2);
        result.element({2, 5, 2, 2}) = (-2 *
                                        (2 * mSqq2 + mSqq * (2 * rCSq + rtSq - 2 * rC * (-1 + w) - 2 * w) +
                                         rtSq * (-3 + rCSq - 4 * w + 2 * rC * (2 + w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({2, 5, 2, 3}) =
            (-2 * (-1 + w) * wp1Sq * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (9. * mSqq2);
        result.element({2, 5, 2, 4}) = (-4 * (mSqq + 2 * rtSq) * (rC - w) * (-1 + w) * wp1Sq) / (9. * mSqq2);
        result.element({2, 5, 2, 5}) =
            (2 * (1 + w) *
             (4 * mSqq2 + 2 * mSqq * (rCSq + rtSq - 2 * rC * w + wSq) + rtSq * (-3 + rCSq - 2 * rC * w + 4 * wSq))) /
            (9. * mSqq2);
        result.element({2, 5, 3, 0}) = (-2 * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({2, 5, 4, 0}) = (-2 * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({2, 5, 5, 10}) = (8 * rt * (1 + 2 * rC - w) * w2m1) / (3. * mSqq);
        result.element({2, 5, 5, 11}) = (8 * (1 + rC) * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({2, 5, 5, 12}) = (16 * rt * (rC - w) * w2m1) / (3. * mSqq);
        result.element({2, 5, 5, 13}) = (8 * rt * (-1 + 3 * rC - 2 * w) * (1 + w)) / (3. * mSqq);
        result.element({2, 5, 5, 14}) = (8 * rt * (rC - w) * (1 + w)) / mSqq;
        result.element({2, 5, 5, 15}) = (-8 * rt * (rC - w) * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({2, 6, 0, 6}) = (-4 *
                                        (2 * mSqq2 * (1 + w) + mSqq * (1 - 2 * rC + rCSq + rtSq) * (1 + w) +
                                         rtSq * (-1 + 2 * rC * (-2 + w) + 2 * w + rCSq * (-1 + 2 * w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({2, 6, 0, 7}) =
            (-2 * wm1Sq * (1 + w) *
             (2 * mSqq * (-1 + rC) * rC * (1 + w) + rtSq * (-3 + 2 * rC * (-2 + w) + rCSq * (1 + 4 * w)))) /
            (9. * mSqq2);
        result.element({2, 6, 0, 8}) =
            (-2 * wm1Sq * (1 + w) *
             (rtSq * (-1 + 3 * rCSq - 2 * rC * (-2 + w) - 4 * w) + 2 * mSqq * (-1 + rC) * (1 + w))) /
            (9. * mSqq2);
        result.element({2, 6, 0, 9}) = (2 *
                                        (2 * mSqq2 + rtSq * (-3 + rCSq + 2 * rC * (-2 + w) + 4 * w) +
                                         mSqq * (2 * rCSq + rtSq + 2 * w - 2 * rC * (1 + w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({2, 6, 1, 6}) = (-4 *
                                        (2 * mSqq2 * (1 + w) + mSqq * (1 - 2 * rC + rCSq + rtSq) * (1 + w) +
                                         rtSq * (-1 + 2 * rC * (-2 + w) + 2 * w + rCSq * (-1 + 2 * w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({2, 6, 1, 7}) =
            (-2 * wm1Sq * (1 + w) *
             (2 * mSqq * (-1 + rC) * rC * (1 + w) + rtSq * (-3 + 2 * rC * (-2 + w) + rCSq * (1 + 4 * w)))) /
            (9. * mSqq2);
        result.element({2, 6, 1, 8}) =
            (-2 * wm1Sq * (1 + w) *
             (rtSq * (-1 + 3 * rCSq - 2 * rC * (-2 + w) - 4 * w) + 2 * mSqq * (-1 + rC) * (1 + w))) /
            (9. * mSqq2);
        result.element({2, 6, 1, 9}) = (2 *
                                        (2 * mSqq2 + rtSq * (-3 + rCSq + 2 * rC * (-2 + w) + 4 * w) +
                                         mSqq * (2 * rCSq + rtSq + 2 * w - 2 * rC * (1 + w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({2, 6, 2, 6}) = (4 *
                                        (2 * mSqq2 * (1 + w) + mSqq * (1 - 2 * rC + rCSq + rtSq) * (1 + w) +
                                         rtSq * (-1 + 2 * rC * (-2 + w) + 2 * w + rCSq * (-1 + 2 * w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({2, 6, 2, 7}) =
            (2 * wm1Sq * (1 + w) *
             (2 * mSqq * (-1 + rC) * rC * (1 + w) + rtSq * (-3 + 2 * rC * (-2 + w) + rCSq * (1 + 4 * w)))) /
            (9. * mSqq2);
        result.element({2, 6, 2, 8}) =
            (2 * wm1Sq * (1 + w) *
             (rtSq * (-1 + 3 * rCSq - 2 * rC * (-2 + w) - 4 * w) + 2 * mSqq * (-1 + rC) * (1 + w))) /
            (9. * mSqq2);
        result.element({2, 6, 2, 9}) = (-2 *
                                        (2 * mSqq2 + rtSq * (-3 + rCSq + 2 * rC * (-2 + w) + 4 * w) +
                                         mSqq * (2 * rCSq + rtSq + 2 * w - 2 * rC * (1 + w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({2, 6, 3, 1}) = (2 * (1 + rC) * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({2, 6, 4, 1}) = (-2 * (1 + rC) * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({2, 6, 5, 10}) = (8 * (-1 + rC) * rt * (-1 + w) * wp1Sq) / mSqq;
        result.element({2, 6, 5, 11}) = (16 * rC * rt * w2m1Sq) / (3. * mSqq);
        result.element({2, 6, 5, 12}) = (16 * rt * w2m1Sq) / (3. * mSqq);
        result.element({2, 6, 5, 13}) = (8 * rt * (2 - 2 * rC + 2 * w) * w2m1) / (3. * mSqq);
        result.element({2, 6, 5, 14}) = (8 * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({2, 7, 0, 6}) =
            (-2 * wm1Sq * (1 + w) *
             (2 * mSqq * (-1 + rC) * rC * (1 + w) + rtSq * (-3 + 2 * rC * (-2 + w) + rCSq * (1 + 4 * w)))) /
            (9. * mSqq2);
        result.element({2, 7, 0, 7}) =
            (-2 * wm1Sq * (1 + w) * (rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)) + 2 * mSqq * rCSq * w2m1)) /
            (9. * mSqq2);
        result.element({2, 7, 0, 8}) =
            (-2 * wm1Sq * (1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) /
            (9. * mSqq2);
        result.element({2, 7, 0, 9}) =
            (2 * wm1Sq * (1 + w) * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (9. * mSqq2);
        result.element({2, 7, 1, 6}) =
            (-2 * wm1Sq * (1 + w) *
             (2 * mSqq * (-1 + rC) * rC * (1 + w) + rtSq * (-3 + 2 * rC * (-2 + w) + rCSq * (1 + 4 * w)))) /
            (9. * mSqq2);
        result.element({2, 7, 1, 7}) =
            (-2 * wm1Sq * (1 + w) * (rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)) + 2 * mSqq * rCSq * w2m1)) /
            (9. * mSqq2);
        result.element({2, 7, 1, 8}) =
            (-2 * wm1Sq * (1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) /
            (9. * mSqq2);
        result.element({2, 7, 1, 9}) =
            (2 * wm1Sq * (1 + w) * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (9. * mSqq2);
        result.element({2, 7, 2, 6}) =
            (2 * wm1Sq * (1 + w) *
             (2 * mSqq * (-1 + rC) * rC * (1 + w) + rtSq * (-3 + 2 * rC * (-2 + w) + rCSq * (1 + 4 * w)))) /
            (9. * mSqq2);
        result.element({2, 7, 2, 7}) =
            (2 * wm1Sq * (1 + w) * (rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)) + 2 * mSqq * rCSq * w2m1)) /
            (9. * mSqq2);
        result.element({2, 7, 2, 8}) =
            (2 * wm1Sq * (1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) /
            (9. * mSqq2);
        result.element({2, 7, 2, 9}) =
            (-2 * wm1Sq * (1 + w) * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (9. * mSqq2);
        result.element({2, 7, 3, 1}) = (2 * rt * wm1Sq * (1 + w) * (-1 + rC * w)) / (3. * mSqq);
        result.element({2, 7, 4, 1}) = (-2 * rt * wm1Sq * (1 + w) * (-1 + rC * w)) / (3. * mSqq);
        result.element({2, 7, 5, 10}) = (8 * rC * rt * w2m1Sq) / (3. * mSqq);
        result.element({2, 7, 5, 13}) = (-8 * rC * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({2, 8, 0, 6}) =
            (-2 * wm1Sq * (1 + w) *
             (rtSq * (-1 + 3 * rCSq - 2 * rC * (-2 + w) - 4 * w) + 2 * mSqq * (-1 + rC) * (1 + w))) /
            (9. * mSqq2);
        result.element({2, 8, 0, 7}) =
            (-2 * wm1Sq * (1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) /
            (9. * mSqq2);
        result.element({2, 8, 0, 8}) =
            (-2 * wm1Sq * (1 + w) * (rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq) + 2 * mSqq * w2m1)) / (9. * mSqq2);
        result.element({2, 8, 0, 9}) = (4 * (mSqq + 2 * rtSq) * (rC - w) * wm1Sq * (1 + w)) / (9. * mSqq2);
        result.element({2, 8, 1, 6}) =
            (-2 * wm1Sq * (1 + w) *
             (rtSq * (-1 + 3 * rCSq - 2 * rC * (-2 + w) - 4 * w) + 2 * mSqq * (-1 + rC) * (1 + w))) /
            (9. * mSqq2);
        result.element({2, 8, 1, 7}) =
            (-2 * wm1Sq * (1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) /
            (9. * mSqq2);
        result.element({2, 8, 1, 8}) =
            (-2 * wm1Sq * (1 + w) * (rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq) + 2 * mSqq * w2m1)) / (9. * mSqq2);
        result.element({2, 8, 1, 9}) = (4 * (mSqq + 2 * rtSq) * (rC - w) * wm1Sq * (1 + w)) / (9. * mSqq2);
        result.element({2, 8, 2, 6}) =
            (2 * wm1Sq * (1 + w) *
             (rtSq * (-1 + 3 * rCSq - 2 * rC * (-2 + w) - 4 * w) + 2 * mSqq * (-1 + rC) * (1 + w))) /
            (9. * mSqq2);
        result.element({2, 8, 2, 7}) =
            (2 * wm1Sq * (1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) /
            (9. * mSqq2);
        result.element({2, 8, 2, 8}) =
            (2 * wm1Sq * (1 + w) * (rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq) + 2 * mSqq * w2m1)) / (9. * mSqq2);
        result.element({2, 8, 2, 9}) = (-4 * (mSqq + 2 * rtSq) * (rC - w) * wm1Sq * (1 + w)) / (9. * mSqq2);
        result.element({2, 8, 3, 1}) = (2 * rt * (rC - w) * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({2, 8, 4, 1}) = (-2 * rt * (rC - w) * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({2, 8, 5, 10}) = (8 * rt * w2m1Sq) / (3. * mSqq);
        result.element({2, 8, 5, 13}) = (-8 * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({2, 9, 0, 6}) = (2 *
                                        (2 * mSqq2 + rtSq * (-3 + rCSq + 2 * rC * (-2 + w) + 4 * w) +
                                         mSqq * (2 * rCSq + rtSq + 2 * w - 2 * rC * (1 + w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({2, 9, 0, 7}) =
            (2 * wm1Sq * (1 + w) * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (9. * mSqq2);
        result.element({2, 9, 0, 8}) = (4 * (mSqq + 2 * rtSq) * (rC - w) * wm1Sq * (1 + w)) / (9. * mSqq2);
        result.element({2, 9, 0, 9}) =
            (-2 * (-1 + w) *
             (4 * mSqq2 + 2 * mSqq * (rCSq + rtSq - 2 * rC * w + wSq) + rtSq * (-3 + rCSq - 2 * rC * w + 4 * wSq))) /
            (9. * mSqq2);
        result.element({2, 9, 1, 6}) = (2 *
                                        (2 * mSqq2 + rtSq * (-3 + rCSq + 2 * rC * (-2 + w) + 4 * w) +
                                         mSqq * (2 * rCSq + rtSq + 2 * w - 2 * rC * (1 + w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({2, 9, 1, 7}) =
            (2 * wm1Sq * (1 + w) * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (9. * mSqq2);
        result.element({2, 9, 1, 8}) = (4 * (mSqq + 2 * rtSq) * (rC - w) * wm1Sq * (1 + w)) / (9. * mSqq2);
        result.element({2, 9, 1, 9}) =
            (-2 * (-1 + w) *
             (4 * mSqq2 + 2 * mSqq * (rCSq + rtSq - 2 * rC * w + wSq) + rtSq * (-3 + rCSq - 2 * rC * w + 4 * wSq))) /
            (9. * mSqq2);
        result.element({2, 9, 2, 6}) = (-2 *
                                        (2 * mSqq2 + rtSq * (-3 + rCSq + 2 * rC * (-2 + w) + 4 * w) +
                                         mSqq * (2 * rCSq + rtSq + 2 * w - 2 * rC * (1 + w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({2, 9, 2, 7}) =
            (-2 * wm1Sq * (1 + w) * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (9. * mSqq2);
        result.element({2, 9, 2, 8}) = (-4 * (mSqq + 2 * rtSq) * (rC - w) * wm1Sq * (1 + w)) / (9. * mSqq2);
        result.element({2, 9, 2, 9}) =
            (2 * (-1 + w) *
             (4 * mSqq2 + 2 * mSqq * (rCSq + rtSq - 2 * rC * w + wSq) + rtSq * (-3 + rCSq - 2 * rC * w + 4 * wSq))) /
            (9. * mSqq2);
        result.element({2, 9, 3, 1}) = (-2 * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({2, 9, 4, 1}) = (2 * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({2, 9, 5, 10}) = (8 * rt * (1 - 2 * rC + w) * w2m1) / (3. * mSqq);
        result.element({2, 9, 5, 11}) = (-8 * rC * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({2, 9, 5, 12}) = (-8 * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({2, 9, 5, 13}) = (8 * rt * (1 + 3 * rC - 2 * w) * (-1 + w)) / (3. * mSqq);
        result.element({2, 9, 5, 14}) = (8 * rt * w2m1) / (3. * mSqq);
        result.element({3, 0, 0, 2}) = (2 * (-1 + rC) * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({3, 0, 0, 3}) = (2 * rt * (-1 + w) * wp1Sq * (-1 + rC * w)) / (3. * mSqq);
        result.element({3, 0, 0, 4}) = (2 * rt * (rC - w) * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({3, 0, 0, 5}) = (-2 * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({3, 0, 1, 2}) = (2 * (-1 + rC) * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({3, 0, 1, 3}) = (2 * rt * (-1 + w) * wp1Sq * (-1 + rC * w)) / (3. * mSqq);
        result.element({3, 0, 1, 4}) = (2 * rt * (rC - w) * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({3, 0, 1, 5}) = (-2 * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({3, 0, 2, 2}) = (2 * (-1 + rC) * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({3, 0, 2, 3}) = (2 * rt * (-1 + w) * wp1Sq * (-1 + rC * w)) / (3. * mSqq);
        result.element({3, 0, 2, 4}) = (2 * rt * (rC - w) * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({3, 0, 2, 5}) = (-2 * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({3, 0, 3, 0}) = (2 * (-1 + w) * wp1Sq) / 3.;
        result.element({3, 0, 4, 0}) = (2 * (-1 + w) * wp1Sq) / 3.;
        result.element({3, 1, 0, 6}) = (-2 * (1 + rC) * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({3, 1, 0, 7}) = (-2 * rt * wm1Sq * (1 + w) * (-1 + rC * w)) / (3. * mSqq);
        result.element({3, 1, 0, 8}) = (-2 * rt * (rC - w) * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({3, 1, 0, 9}) = (2 * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({3, 1, 1, 6}) = (-2 * (1 + rC) * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({3, 1, 1, 7}) = (-2 * rt * wm1Sq * (1 + w) * (-1 + rC * w)) / (3. * mSqq);
        result.element({3, 1, 1, 8}) = (-2 * rt * (rC - w) * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({3, 1, 1, 9}) = (2 * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({3, 1, 2, 6}) = (2 * (1 + rC) * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({3, 1, 2, 7}) = (2 * rt * wm1Sq * (1 + w) * (-1 + rC * w)) / (3. * mSqq);
        result.element({3, 1, 2, 8}) = (2 * rt * (rC - w) * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({3, 1, 2, 9}) = (-2 * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({3, 1, 3, 1}) = (2 * wm1Sq * (1 + w)) / 3.;
        result.element({3, 1, 4, 1}) = (-2 * wm1Sq * (1 + w)) / 3.;
        result.element({4, 0, 0, 2}) = (2 * (-1 + rC) * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({4, 0, 0, 3}) = (2 * rt * (-1 + w) * wp1Sq * (-1 + rC * w)) / (3. * mSqq);
        result.element({4, 0, 0, 4}) = (2 * rt * (rC - w) * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({4, 0, 0, 5}) = (-2 * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({4, 0, 1, 2}) = (2 * (-1 + rC) * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({4, 0, 1, 3}) = (2 * rt * (-1 + w) * wp1Sq * (-1 + rC * w)) / (3. * mSqq);
        result.element({4, 0, 1, 4}) = (2 * rt * (rC - w) * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({4, 0, 1, 5}) = (-2 * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({4, 0, 2, 2}) = (2 * (-1 + rC) * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({4, 0, 2, 3}) = (2 * rt * (-1 + w) * wp1Sq * (-1 + rC * w)) / (3. * mSqq);
        result.element({4, 0, 2, 4}) = (2 * rt * (rC - w) * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({4, 0, 2, 5}) = (-2 * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({4, 0, 3, 0}) = (2 * (-1 + w) * wp1Sq) / 3.;
        result.element({4, 0, 4, 0}) = (2 * (-1 + w) * wp1Sq) / 3.;
        result.element({4, 1, 0, 6}) = (2 * (1 + rC) * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({4, 1, 0, 7}) = (2 * rt * wm1Sq * (1 + w) * (-1 + rC * w)) / (3. * mSqq);
        result.element({4, 1, 0, 8}) = (2 * rt * (rC - w) * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({4, 1, 0, 9}) = (-2 * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({4, 1, 1, 6}) = (2 * (1 + rC) * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({4, 1, 1, 7}) = (2 * rt * wm1Sq * (1 + w) * (-1 + rC * w)) / (3. * mSqq);
        result.element({4, 1, 1, 8}) = (2 * rt * (rC - w) * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({4, 1, 1, 9}) = (-2 * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({4, 1, 2, 6}) = (-2 * (1 + rC) * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({4, 1, 2, 7}) = (-2 * rt * wm1Sq * (1 + w) * (-1 + rC * w)) / (3. * mSqq);
        result.element({4, 1, 2, 8}) = (-2 * rt * (rC - w) * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({4, 1, 2, 9}) = (2 * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({4, 1, 3, 1}) = (-2 * wm1Sq * (1 + w)) / 3.;
        result.element({4, 1, 4, 1}) = (2 * wm1Sq * (1 + w)) / 3.;
        result.element({5, 10, 0, 2}) = (-8 * (1 + rC) * rt * wm1Sq * (1 + w)) / mSqq;
        result.element({5, 10, 0, 3}) = (-8 * rC * rt * w2m1Sq) / (3. * mSqq);
        result.element({5, 10, 0, 4}) = (-8 * rt * w2m1Sq) / (3. * mSqq);
        result.element({5, 10, 0, 5}) = (8 * rt * (1 + 2 * rC - w) * w2m1) / (3. * mSqq);
        result.element({5, 10, 0, 6}) = (-8 * (-1 + rC) * rt * (-1 + w) * wp1Sq) / mSqq;
        result.element({5, 10, 0, 7}) = (-8 * rC * rt * w2m1Sq) / (3. * mSqq);
        result.element({5, 10, 0, 8}) = (-8 * rt * w2m1Sq) / (3. * mSqq);
        result.element({5, 10, 0, 9}) = (-8 * rt * (1 - 2 * rC + w) * w2m1) / (3. * mSqq);
        result.element({5, 10, 1, 2}) = (-8 * (1 + rC) * rt * wm1Sq * (1 + w)) / mSqq;
        result.element({5, 10, 1, 3}) = (-8 * rC * rt * w2m1Sq) / (3. * mSqq);
        result.element({5, 10, 1, 4}) = (-8 * rt * w2m1Sq) / (3. * mSqq);
        result.element({5, 10, 1, 5}) = (8 * rt * (1 + 2 * rC - w) * w2m1) / (3. * mSqq);
        result.element({5, 10, 1, 6}) = (-8 * (-1 + rC) * rt * (-1 + w) * wp1Sq) / mSqq;
        result.element({5, 10, 1, 7}) = (-8 * rC * rt * w2m1Sq) / (3. * mSqq);
        result.element({5, 10, 1, 8}) = (-8 * rt * w2m1Sq) / (3. * mSqq);
        result.element({5, 10, 1, 9}) = (-8 * rt * (1 - 2 * rC + w) * w2m1) / (3. * mSqq);
        result.element({5, 10, 2, 2}) = (-8 * (1 + rC) * rt * wm1Sq * (1 + w)) / mSqq;
        result.element({5, 10, 2, 3}) = (-8 * rC * rt * w2m1Sq) / (3. * mSqq);
        result.element({5, 10, 2, 4}) = (-8 * rt * w2m1Sq) / (3. * mSqq);
        result.element({5, 10, 2, 5}) = (8 * rt * (1 + 2 * rC - w) * w2m1) / (3. * mSqq);
        result.element({5, 10, 2, 6}) = (8 * (-1 + rC) * rt * (-1 + w) * wp1Sq) / mSqq;
        result.element({5, 10, 2, 7}) = (8 * rC * rt * w2m1Sq) / (3. * mSqq);
        result.element({5, 10, 2, 8}) = (8 * rt * w2m1Sq) / (3. * mSqq);
        result.element({5, 10, 2, 9}) = (8 * rt * (1 - 2 * rC + w) * w2m1) / (3. * mSqq);
        result.element({5, 10, 5, 10}) =
            (64 * (mSqq + 2 * rtSq) * (-4 * rC + (2 + mSqq) * w + 2 * rCSq * w) * w2m1) / (9. * mSqq2);
        result.element({5, 10, 5, 11}) =
            (-32 * (mSqq + 2 * rtSq) * wm1Sq * (1 + w) * (2 + mSqq + 4 * rC - 2 * rCSq * (1 + 2 * w))) / (9. * mSqq2);
        result.element({5, 10, 5, 12}) =
            (32 * (mSqq + 2 * rtSq) * (-2 + mSqq + 4 * rC + 2 * rCSq - 4 * w) * wm1Sq * (1 + w)) / (9. * mSqq2);
        result.element({5, 10, 5, 13}) = (128 * (mSqq + 2 * rtSq) * (rC - w) * w2m1) / (9. * mSqq2);
        result.element({5, 10, 5, 14}) =
            (32 * (mSqq + 2 * rtSq) * (-1 + mSqq + 2 * rC + rCSq - 2 * w) * w2m1) / (9. * mSqq2);
        result.element({5, 10, 5, 15}) = (-32 * (mSqq + 2 * rtSq) * w2m1Sq) / (9. * mSqq);
        result.element({5, 11, 0, 2}) = (-8 * rt * wm1Sq * (1 + w) * (-3 + rC * (-1 + 2 * w))) / (3. * mSqq);
        result.element({5, 11, 0, 3}) = (8 * rC * rt * w2m1Sq) / (3. * mSqq);
        result.element({5, 11, 0, 4}) = (8 * rt * w2m1Sq) / (3. * mSqq);
        result.element({5, 11, 0, 5}) = (8 * (1 + rC) * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({5, 11, 0, 6}) = (-16 * rC * rt * w2m1Sq) / (3. * mSqq);
        result.element({5, 11, 0, 9}) = (8 * rC * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({5, 11, 1, 2}) = (-8 * rt * wm1Sq * (1 + w) * (-3 + rC * (-1 + 2 * w))) / (3. * mSqq);
        result.element({5, 11, 1, 3}) = (8 * rC * rt * w2m1Sq) / (3. * mSqq);
        result.element({5, 11, 1, 4}) = (8 * rt * w2m1Sq) / (3. * mSqq);
        result.element({5, 11, 1, 5}) = (8 * (1 + rC) * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({5, 11, 1, 6}) = (-16 * rC * rt * w2m1Sq) / (3. * mSqq);
        result.element({5, 11, 1, 9}) = (8 * rC * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({5, 11, 2, 2}) = (-8 * rt * wm1Sq * (1 + w) * (-3 + rC * (-1 + 2 * w))) / (3. * mSqq);
        result.element({5, 11, 2, 3}) = (8 * rC * rt * w2m1Sq) / (3. * mSqq);
        result.element({5, 11, 2, 4}) = (8 * rt * w2m1Sq) / (3. * mSqq);
        result.element({5, 11, 2, 5}) = (8 * (1 + rC) * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({5, 11, 2, 6}) = (16 * rC * rt * w2m1Sq) / (3. * mSqq);
        result.element({5, 11, 2, 9}) = (-8 * rC * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({5, 11, 5, 10}) =
            (-32 * (mSqq + 2 * rtSq) * wm1Sq * (1 + w) * (2 + mSqq + 4 * rC - 2 * rCSq * (1 + 2 * w))) / (9. * mSqq2);
        result.element({5, 11, 5, 11}) = (32 * wm1Sq * (1 + w) *
                                          (mSqq2 + 4 * rtSq * (1 - 2 * rC * w + rCSq * (-1 + 2 * wSq)) +
                                           2 * mSqq * (1 + rtSq - 2 * rC * w + rCSq * (-1 + 2 * wSq)))) /
                                         (9. * mSqq2);
        result.element({5, 11, 5, 12}) =
            (-32 * (mSqq + 2 * rtSq) * wm1Sq * (1 + w) * (mSqq - 2 * (-2 * rC + w + rCSq * w))) / (9. * mSqq2);
        result.element({5, 11, 5, 13}) =
            (-32 * (mSqq + 2 * rtSq) * (1 + mSqq + 2 * rC - rCSq - 2 * w) * w2m1) / (9. * mSqq2);
        result.element({5, 11, 5, 14}) =
            (-32 * (mSqq + 2 * rtSq) * (mSqq + 2 * rC - w - rCSq * w) * w2m1) / (9. * mSqq2);
        result.element({5, 11, 5, 15}) = (32 * (mSqq + 2 * rtSq) * w2m1Sq) / (9. * mSqq);
        result.element({5, 12, 0, 2}) = (-8 * rt * (1 + 3 * rC - 2 * w) * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({5, 12, 0, 3}) = (-8 * rC * rt * w2m1Sq) / (3. * mSqq);
        result.element({5, 12, 0, 4}) = (-8 * rt * w2m1Sq) / (3. * mSqq);
        result.element({5, 12, 0, 5}) = (16 * rt * (rC - w) * w2m1) / (3. * mSqq);
        result.element({5, 12, 0, 6}) = (-16 * rt * w2m1Sq) / (3. * mSqq);
        result.element({5, 12, 0, 9}) = (8 * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({5, 12, 1, 2}) = (-8 * rt * (1 + 3 * rC - 2 * w) * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({5, 12, 1, 3}) = (-8 * rC * rt * w2m1Sq) / (3. * mSqq);
        result.element({5, 12, 1, 4}) = (-8 * rt * w2m1Sq) / (3. * mSqq);
        result.element({5, 12, 1, 5}) = (16 * rt * (rC - w) * w2m1) / (3. * mSqq);
        result.element({5, 12, 1, 6}) = (-16 * rt * w2m1Sq) / (3. * mSqq);
        result.element({5, 12, 1, 9}) = (8 * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({5, 12, 2, 2}) = (-8 * rt * (1 + 3 * rC - 2 * w) * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({5, 12, 2, 3}) = (-8 * rC * rt * w2m1Sq) / (3. * mSqq);
        result.element({5, 12, 2, 4}) = (-8 * rt * w2m1Sq) / (3. * mSqq);
        result.element({5, 12, 2, 5}) = (16 * rt * (rC - w) * w2m1) / (3. * mSqq);
        result.element({5, 12, 2, 6}) = (16 * rt * w2m1Sq) / (3. * mSqq);
        result.element({5, 12, 2, 9}) = (-8 * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({5, 12, 5, 10}) =
            (32 * (mSqq + 2 * rtSq) * (-2 + mSqq + 4 * rC + 2 * rCSq - 4 * w) * wm1Sq * (1 + w)) / (9. * mSqq2);
        result.element({5, 12, 5, 11}) =
            (-32 * (mSqq + 2 * rtSq) * wm1Sq * (1 + w) * (mSqq - 2 * (-2 * rC + w + rCSq * w))) / (9. * mSqq2);
        result.element({5, 12, 5, 12}) = (32 * wm1Sq * (1 + w) *
                                          (mSqq2 + 4 * rtSq * (-1 + rCSq - 2 * rC * w + 2 * wSq) +
                                           2 * mSqq * (-1 + rCSq + rtSq - 2 * rC * w + 2 * wSq))) /
                                         (9. * mSqq2);
        result.element({5, 12, 5, 13}) =
            (32 * (mSqq + 2 * rtSq) * (-1 + mSqq + rCSq + rC * (2 - 4 * w) - 2 * w + 4 * wSq) * w2m1) / (9. * mSqq2);
        result.element({5, 12, 5, 14}) = (32 *
                                          (mSqq2 + 2 * rtSq * (-1 + rCSq - 2 * rC * w + 2 * wSq) +
                                           mSqq * (-1 + rCSq + 2 * rtSq - 2 * rC * w + 2 * wSq)) *
                                          w2m1) /
                                         (9. * mSqq2);
        result.element({5, 12, 5, 15}) = (-32 * (mSqq + 2 * rtSq) * w2m1Sq) / (9. * mSqq);
        result.element({5, 13, 0, 2}) = (-16 * rt * (1 + rC - w) * w2m1) / (3. * mSqq);
        result.element({5, 13, 0, 3}) = (-8 * rC * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({5, 13, 0, 4}) = (-8 * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({5, 13, 0, 5}) = (8 * rt * (-1 + 3 * rC - 2 * w) * (1 + w)) / (3. * mSqq);
        result.element({5, 13, 0, 6}) = (16 * rt * (-1 + rC - w) * w2m1) / (3. * mSqq);
        result.element({5, 13, 0, 7}) = (8 * rC * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({5, 13, 0, 8}) = (8 * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({5, 13, 0, 9}) = (-8 * rt * (1 + 3 * rC - 2 * w) * (-1 + w)) / (3. * mSqq);
        result.element({5, 13, 1, 2}) = (-16 * rt * (1 + rC - w) * w2m1) / (3. * mSqq);
        result.element({5, 13, 1, 3}) = (-8 * rC * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({5, 13, 1, 4}) = (-8 * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({5, 13, 1, 5}) = (8 * rt * (-1 + 3 * rC - 2 * w) * (1 + w)) / (3. * mSqq);
        result.element({5, 13, 1, 6}) = (16 * rt * (-1 + rC - w) * w2m1) / (3. * mSqq);
        result.element({5, 13, 1, 7}) = (8 * rC * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({5, 13, 1, 8}) = (8 * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({5, 13, 1, 9}) = (-8 * rt * (1 + 3 * rC - 2 * w) * (-1 + w)) / (3. * mSqq);
        result.element({5, 13, 2, 2}) = (-16 * rt * (1 + rC - w) * w2m1) / (3. * mSqq);
        result.element({5, 13, 2, 3}) = (-8 * rC * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({5, 13, 2, 4}) = (-8 * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({5, 13, 2, 5}) = (8 * rt * (-1 + 3 * rC - 2 * w) * (1 + w)) / (3. * mSqq);
        result.element({5, 13, 2, 6}) = (8 * rt * (2 - 2 * rC + 2 * w) * w2m1) / (3. * mSqq);
        result.element({5, 13, 2, 7}) = (-8 * rC * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({5, 13, 2, 8}) = (-8 * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({5, 13, 2, 9}) = (8 * rt * (1 + 3 * rC - 2 * w) * (-1 + w)) / (3. * mSqq);
        result.element({5, 13, 5, 10}) = (128 * (mSqq + 2 * rtSq) * (rC - w) * w2m1) / (9. * mSqq2);
        result.element({5, 13, 5, 11}) =
            (-32 * (mSqq + 2 * rtSq) * (1 + mSqq + 2 * rC - rCSq - 2 * w) * w2m1) / (9. * mSqq2);
        result.element({5, 13, 5, 12}) =
            (32 * (mSqq + 2 * rtSq) * (-1 + mSqq + rCSq + rC * (2 - 4 * w) - 2 * w + 4 * wSq) * w2m1) / (9. * mSqq2);
        result.element({5, 13, 5, 13}) = (64 * (mSqq2 * w + 4 * rtSq * (rCSq * w + wCu - rC * (1 + wSq)) +
                                                2 * mSqq * (rCSq * w - rC * (1 + wSq) + w * (rtSq + wSq)))) /
                                         (9. * mSqq2);
        result.element({5, 13, 5, 14}) =
            (32 * (mSqq + 2 * rtSq) * (1 + w) * (mSqq + 2 * (rCSq + w - rC * (1 + w) + w2m1))) / (9. * mSqq2);
        result.element({5, 13, 5, 15}) = (-32 * (mSqq + 2 * rtSq) * (-1 + w) * wp1Sq) / (9. * mSqq);
        result.element({5, 14, 0, 2}) = (-8 * rt * (1 + 2 * rC - w) * w2m1) / (3. * mSqq);
        result.element({5, 14, 0, 3}) = (-8 * rC * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({5, 14, 0, 4}) = (-8 * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({5, 14, 0, 5}) = (8 * rt * (rC - w) * (1 + w)) / mSqq;
        result.element({5, 14, 0, 6}) = (-8 * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({5, 14, 0, 9}) = (-8 * rt * w2m1) / (3. * mSqq);
        result.element({5, 14, 1, 2}) = (-8 * rt * (1 + 2 * rC - w) * w2m1) / (3. * mSqq);
        result.element({5, 14, 1, 3}) = (-8 * rC * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({5, 14, 1, 4}) = (-8 * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({5, 14, 1, 5}) = (8 * rt * (rC - w) * (1 + w)) / mSqq;
        result.element({5, 14, 1, 6}) = (-8 * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({5, 14, 1, 9}) = (-8 * rt * w2m1) / (3. * mSqq);
        result.element({5, 14, 2, 2}) = (-8 * rt * (1 + 2 * rC - w) * w2m1) / (3. * mSqq);
        result.element({5, 14, 2, 3}) = (-8 * rC * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({5, 14, 2, 4}) = (-8 * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({5, 14, 2, 5}) = (8 * rt * (rC - w) * (1 + w)) / mSqq;
        result.element({5, 14, 2, 6}) = (8 * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({5, 14, 2, 9}) = (8 * rt * w2m1) / (3. * mSqq);
        result.element({5, 14, 5, 10}) =
            (32 * (mSqq + 2 * rtSq) * (-1 + mSqq + 2 * rC + rCSq - 2 * w) * w2m1) / (9. * mSqq2);
        result.element({5, 14, 5, 11}) =
            (-32 * (mSqq + 2 * rtSq) * (mSqq + 2 * rC - w - rCSq * w) * w2m1) / (9. * mSqq2);
        result.element({5, 14, 5, 12}) = (32 *
                                          (mSqq2 + 2 * rtSq * (-1 + rCSq - 2 * rC * w + 2 * wSq) +
                                           mSqq * (-1 + rCSq + 2 * rtSq - 2 * rC * w + 2 * wSq)) *
                                          w2m1) /
                                         (9. * mSqq2);
        result.element({5, 14, 5, 13}) =
            (32 * (mSqq + 2 * rtSq) * (1 + w) * (mSqq + 2 * (rCSq + w - rC * (1 + w) + w2m1))) / (9. * mSqq2);
        result.element({5, 14, 5, 14}) = (32 * (1 + w) *
                                          (mSqq2 + 4 * rtSq * (-1 + rCSq - 2 * rC * w + 2 * wSq) +
                                           2 * mSqq * (-1 + rCSq + rtSq - 2 * rC * w + 2 * wSq))) /
                                         (9. * mSqq2);
        result.element({5, 14, 5, 15}) = (-32 * (mSqq + 2 * rtSq) * (-1 + w) * wp1Sq) / (9. * mSqq);
        result.element({5, 15, 0, 2}) = (8 * (1 + rC) * rt * w2m1Sq) / (3. * mSqq);
        result.element({5, 15, 0, 3}) = (8 * rC * rt * wm1Sq * wp1Cu) / (3. * mSqq);
        result.element({5, 15, 0, 4}) = (8 * rt * wm1Sq * wp1Cu) / (3. * mSqq);
        result.element({5, 15, 0, 5}) = (-8 * rt * (rC - w) * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({5, 15, 1, 2}) = (8 * (1 + rC) * rt * w2m1Sq) / (3. * mSqq);
        result.element({5, 15, 1, 3}) = (8 * rC * rt * wm1Sq * wp1Cu) / (3. * mSqq);
        result.element({5, 15, 1, 4}) = (8 * rt * wm1Sq * wp1Cu) / (3. * mSqq);
        result.element({5, 15, 1, 5}) = (-8 * rt * (rC - w) * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({5, 15, 2, 2}) = (8 * (1 + rC) * rt * w2m1Sq) / (3. * mSqq);
        result.element({5, 15, 2, 3}) = (8 * rC * rt * wm1Sq * wp1Cu) / (3. * mSqq);
        result.element({5, 15, 2, 4}) = (8 * rt * wm1Sq * wp1Cu) / (3. * mSqq);
        result.element({5, 15, 2, 5}) = (-8 * rt * (rC - w) * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({5, 15, 5, 10}) = (-32 * (mSqq + 2 * rtSq) * w2m1Sq) / (9. * mSqq);
        result.element({5, 15, 5, 11}) = (32 * (mSqq + 2 * rtSq) * w2m1Sq) / (9. * mSqq);
        result.element({5, 15, 5, 12}) = (-32 * (mSqq + 2 * rtSq) * w2m1Sq) / (9. * mSqq);
        result.element({5, 15, 5, 13}) = (-32 * (mSqq + 2 * rtSq) * (-1 + w) * wp1Sq) / (9. * mSqq);
        result.element({5, 15, 5, 14}) = (-32 * (mSqq + 2 * rtSq) * (-1 + w) * wp1Sq) / (9. * mSqq);
        result.element({5, 15, 5, 15}) = (32 * (mSqq + 2 * rtSq) * wm1Sq * wp1Cu) / (9. * mSqq);
        result.element({6, 2, 6, 2}) = (4 *
                                        (2 * mSqq2 * (-1 + w) + mSqq * (1 + 2 * rC + rCSq + rtSq) * (-1 + w) +
                                         rtSq * (1 + 2 * w - 2 * rC * (2 + w) + rCSq * (1 + 2 * w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({6, 2, 6, 3}) =
            (2 * (-1 + w) * wp1Sq *
             (2 * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2 * rC * (2 + w) + rCSq * (-1 + 4 * w)))) /
            (9. * mSqq2);
        result.element({6, 2, 6, 4}) =
            (2 * (-1 + w) * wp1Sq *
             (2 * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3 * rCSq + 4 * w - 2 * rC * (2 + w)))) /
            (9. * mSqq2);
        result.element({6, 2, 6, 5}) = (-2 *
                                        (2 * mSqq2 + mSqq * (2 * rCSq + rtSq - 2 * rC * (-1 + w) - 2 * w) +
                                         rtSq * (-3 + rCSq - 4 * w + 2 * rC * (2 + w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({6, 2, 7, 2}) = (4 *
                                        (2 * mSqq2 * (-1 + w) + mSqq * (1 + 2 * rC + rCSq + rtSq) * (-1 + w) +
                                         rtSq * (1 + 2 * w - 2 * rC * (2 + w) + rCSq * (1 + 2 * w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({6, 2, 7, 3}) =
            (2 * (-1 + w) * wp1Sq *
             (2 * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2 * rC * (2 + w) + rCSq * (-1 + 4 * w)))) /
            (9. * mSqq2);
        result.element({6, 2, 7, 4}) =
            (2 * (-1 + w) * wp1Sq *
             (2 * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3 * rCSq + 4 * w - 2 * rC * (2 + w)))) /
            (9. * mSqq2);
        result.element({6, 2, 7, 5}) = (-2 *
                                        (2 * mSqq2 + mSqq * (2 * rCSq + rtSq - 2 * rC * (-1 + w) - 2 * w) +
                                         rtSq * (-3 + rCSq - 4 * w + 2 * rC * (2 + w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({6, 2, 8, 0}) = (2 * (-1 + rC) * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({6, 2, 9, 0}) = (2 * (-1 + rC) * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({6, 2, 10, 10}) = (-8 * (1 + rC) * rt * wm1Sq * (1 + w)) / mSqq;
        result.element({6, 2, 10, 11}) = (-8 * rt * wm1Sq * (1 + w) * (-3 + rC * (-1 + 2 * w))) / (3. * mSqq);
        result.element({6, 2, 10, 12}) = (-8 * rt * (1 + 3 * rC - 2 * w) * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({6, 2, 10, 13}) = (-16 * rt * (1 + rC - w) * w2m1) / (3. * mSqq);
        result.element({6, 2, 10, 14}) = (-8 * rt * (1 + 2 * rC - w) * w2m1) / (3. * mSqq);
        result.element({6, 2, 10, 15}) = (8 * (1 + rC) * rt * w2m1Sq) / (3. * mSqq);
        result.element({6, 3, 6, 2}) =
            (2 * (-1 + w) * wp1Sq *
             (2 * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2 * rC * (2 + w) + rCSq * (-1 + 4 * w)))) /
            (9. * mSqq2);
        result.element({6, 3, 6, 3}) =
            (2 * (-1 + w) * wp1Sq * (rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)) + 2 * mSqq * rCSq * w2m1)) /
            (9. * mSqq2);
        result.element({6, 3, 6, 4}) =
            (2 * (-1 + w) * wp1Sq * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) /
            (9. * mSqq2);
        result.element({6, 3, 6, 5}) =
            (-2 * (-1 + w) * wp1Sq * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (9. * mSqq2);
        result.element({6, 3, 7, 2}) =
            (2 * (-1 + w) * wp1Sq *
             (2 * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2 * rC * (2 + w) + rCSq * (-1 + 4 * w)))) /
            (9. * mSqq2);
        result.element({6, 3, 7, 3}) =
            (2 * (-1 + w) * wp1Sq * (rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)) + 2 * mSqq * rCSq * w2m1)) /
            (9. * mSqq2);
        result.element({6, 3, 7, 4}) =
            (2 * (-1 + w) * wp1Sq * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) /
            (9. * mSqq2);
        result.element({6, 3, 7, 5}) =
            (-2 * (-1 + w) * wp1Sq * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (9. * mSqq2);
        result.element({6, 3, 8, 0}) = (2 * rt * (-1 + w) * wp1Sq * (-1 + rC * w)) / (3. * mSqq);
        result.element({6, 3, 9, 0}) = (2 * rt * (-1 + w) * wp1Sq * (-1 + rC * w)) / (3. * mSqq);
        result.element({6, 3, 10, 10}) = (-8 * rC * rt * w2m1Sq) / (3. * mSqq);
        result.element({6, 3, 10, 11}) = (8 * rC * rt * w2m1Sq) / (3. * mSqq);
        result.element({6, 3, 10, 12}) = (-8 * rC * rt * w2m1Sq) / (3. * mSqq);
        result.element({6, 3, 10, 13}) = (-8 * rC * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({6, 3, 10, 14}) = (-8 * rC * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({6, 3, 10, 15}) = (8 * rC * rt * wm1Sq * wp1Cu) / (3. * mSqq);
        result.element({6, 4, 6, 2}) =
            (2 * (-1 + w) * wp1Sq *
             (2 * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3 * rCSq + 4 * w - 2 * rC * (2 + w)))) /
            (9. * mSqq2);
        result.element({6, 4, 6, 3}) =
            (2 * (-1 + w) * wp1Sq * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) /
            (9. * mSqq2);
        result.element({6, 4, 6, 4}) =
            (2 * (-1 + w) * wp1Sq * (rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq) + 2 * mSqq * w2m1)) / (9. * mSqq2);
        result.element({6, 4, 6, 5}) = (-4 * (mSqq + 2 * rtSq) * (rC - w) * (-1 + w) * wp1Sq) / (9. * mSqq2);
        result.element({6, 4, 7, 2}) =
            (2 * (-1 + w) * wp1Sq *
             (2 * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3 * rCSq + 4 * w - 2 * rC * (2 + w)))) /
            (9. * mSqq2);
        result.element({6, 4, 7, 3}) =
            (2 * (-1 + w) * wp1Sq * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) /
            (9. * mSqq2);
        result.element({6, 4, 7, 4}) =
            (2 * (-1 + w) * wp1Sq * (rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq) + 2 * mSqq * w2m1)) / (9. * mSqq2);
        result.element({6, 4, 7, 5}) = (-4 * (mSqq + 2 * rtSq) * (rC - w) * (-1 + w) * wp1Sq) / (9. * mSqq2);
        result.element({6, 4, 8, 0}) = (2 * rt * (rC - w) * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({6, 4, 9, 0}) = (2 * rt * (rC - w) * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({6, 4, 10, 10}) = (-8 * rt * w2m1Sq) / (3. * mSqq);
        result.element({6, 4, 10, 11}) = (8 * rt * w2m1Sq) / (3. * mSqq);
        result.element({6, 4, 10, 12}) = (-8 * rt * w2m1Sq) / (3. * mSqq);
        result.element({6, 4, 10, 13}) = (-8 * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({6, 4, 10, 14}) = (-8 * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({6, 4, 10, 15}) = (8 * rt * wm1Sq * wp1Cu) / (3. * mSqq);
        result.element({6, 5, 6, 2}) = (-2 *
                                        (2 * mSqq2 + mSqq * (2 * rCSq + rtSq - 2 * rC * (-1 + w) - 2 * w) +
                                         rtSq * (-3 + rCSq - 4 * w + 2 * rC * (2 + w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({6, 5, 6, 3}) =
            (-2 * (-1 + w) * wp1Sq * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (9. * mSqq2);
        result.element({6, 5, 6, 4}) = (-4 * (mSqq + 2 * rtSq) * (rC - w) * (-1 + w) * wp1Sq) / (9. * mSqq2);
        result.element({6, 5, 6, 5}) =
            (2 * (1 + w) *
             (4 * mSqq2 + 2 * mSqq * (rCSq + rtSq - 2 * rC * w + wSq) + rtSq * (-3 + rCSq - 2 * rC * w + 4 * wSq))) /
            (9. * mSqq2);
        result.element({6, 5, 7, 2}) = (-2 *
                                        (2 * mSqq2 + mSqq * (2 * rCSq + rtSq - 2 * rC * (-1 + w) - 2 * w) +
                                         rtSq * (-3 + rCSq - 4 * w + 2 * rC * (2 + w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({6, 5, 7, 3}) =
            (-2 * (-1 + w) * wp1Sq * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (9. * mSqq2);
        result.element({6, 5, 7, 4}) = (-4 * (mSqq + 2 * rtSq) * (rC - w) * (-1 + w) * wp1Sq) / (9. * mSqq2);
        result.element({6, 5, 7, 5}) =
            (2 * (1 + w) *
             (4 * mSqq2 + 2 * mSqq * (rCSq + rtSq - 2 * rC * w + wSq) + rtSq * (-3 + rCSq - 2 * rC * w + 4 * wSq))) /
            (9. * mSqq2);
        result.element({6, 5, 8, 0}) = (-2 * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({6, 5, 9, 0}) = (-2 * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({6, 5, 10, 10}) = (8 * rt * (1 + 2 * rC - w) * w2m1) / (3. * mSqq);
        result.element({6, 5, 10, 11}) = (8 * (1 + rC) * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({6, 5, 10, 12}) = (16 * rt * (rC - w) * w2m1) / (3. * mSqq);
        result.element({6, 5, 10, 13}) = (8 * rt * (-1 + 3 * rC - 2 * w) * (1 + w)) / (3. * mSqq);
        result.element({6, 5, 10, 14}) = (8 * rt * (rC - w) * (1 + w)) / mSqq;
        result.element({6, 5, 10, 15}) = (-8 * rt * (rC - w) * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({6, 6, 6, 6}) = (4 *
                                        (2 * mSqq2 * (1 + w) + mSqq * (1 - 2 * rC + rCSq + rtSq) * (1 + w) +
                                         rtSq * (-1 + 2 * rC * (-2 + w) + 2 * w + rCSq * (-1 + 2 * w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({6, 6, 6, 7}) =
            (2 * wm1Sq * (1 + w) *
             (2 * mSqq * (-1 + rC) * rC * (1 + w) + rtSq * (-3 + 2 * rC * (-2 + w) + rCSq * (1 + 4 * w)))) /
            (9. * mSqq2);
        result.element({6, 6, 6, 8}) =
            (2 * wm1Sq * (1 + w) *
             (rtSq * (-1 + 3 * rCSq - 2 * rC * (-2 + w) - 4 * w) + 2 * mSqq * (-1 + rC) * (1 + w))) /
            (9. * mSqq2);
        result.element({6, 6, 6, 9}) = (-2 *
                                        (2 * mSqq2 + rtSq * (-3 + rCSq + 2 * rC * (-2 + w) + 4 * w) +
                                         mSqq * (2 * rCSq + rtSq + 2 * w - 2 * rC * (1 + w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({6, 6, 7, 6}) = (-4 *
                                        (2 * mSqq2 * (1 + w) + mSqq * (1 - 2 * rC + rCSq + rtSq) * (1 + w) +
                                         rtSq * (-1 + 2 * rC * (-2 + w) + 2 * w + rCSq * (-1 + 2 * w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({6, 6, 7, 7}) =
            (-2 * wm1Sq * (1 + w) *
             (2 * mSqq * (-1 + rC) * rC * (1 + w) + rtSq * (-3 + 2 * rC * (-2 + w) + rCSq * (1 + 4 * w)))) /
            (9. * mSqq2);
        result.element({6, 6, 7, 8}) =
            (-2 * wm1Sq * (1 + w) *
             (rtSq * (-1 + 3 * rCSq - 2 * rC * (-2 + w) - 4 * w) + 2 * mSqq * (-1 + rC) * (1 + w))) /
            (9. * mSqq2);
        result.element({6, 6, 7, 9}) = (2 *
                                        (2 * mSqq2 + rtSq * (-3 + rCSq + 2 * rC * (-2 + w) + 4 * w) +
                                         mSqq * (2 * rCSq + rtSq + 2 * w - 2 * rC * (1 + w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({6, 6, 8, 1}) = (-2 * (1 + rC) * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({6, 6, 9, 1}) = (2 * (1 + rC) * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({6, 6, 10, 10}) = (8 * (-1 + rC) * rt * (-1 + w) * wp1Sq) / mSqq;
        result.element({6, 6, 10, 11}) = (16 * rC * rt * w2m1Sq) / (3. * mSqq);
        result.element({6, 6, 10, 12}) = (16 * rt * w2m1Sq) / (3. * mSqq);
        result.element({6, 6, 10, 13}) = (8 * rt * (2 - 2 * rC + 2 * w) * w2m1) / (3. * mSqq);
        result.element({6, 6, 10, 14}) = (8 * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({6, 7, 6, 6}) =
            (2 * wm1Sq * (1 + w) *
             (2 * mSqq * (-1 + rC) * rC * (1 + w) + rtSq * (-3 + 2 * rC * (-2 + w) + rCSq * (1 + 4 * w)))) /
            (9. * mSqq2);
        result.element({6, 7, 6, 7}) =
            (2 * wm1Sq * (1 + w) * (rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)) + 2 * mSqq * rCSq * w2m1)) /
            (9. * mSqq2);
        result.element({6, 7, 6, 8}) =
            (2 * wm1Sq * (1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) /
            (9. * mSqq2);
        result.element({6, 7, 6, 9}) =
            (-2 * wm1Sq * (1 + w) * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (9. * mSqq2);
        result.element({6, 7, 7, 6}) =
            (-2 * wm1Sq * (1 + w) *
             (2 * mSqq * (-1 + rC) * rC * (1 + w) + rtSq * (-3 + 2 * rC * (-2 + w) + rCSq * (1 + 4 * w)))) /
            (9. * mSqq2);
        result.element({6, 7, 7, 7}) =
            (-2 * wm1Sq * (1 + w) * (rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)) + 2 * mSqq * rCSq * w2m1)) /
            (9. * mSqq2);
        result.element({6, 7, 7, 8}) =
            (-2 * wm1Sq * (1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) /
            (9. * mSqq2);
        result.element({6, 7, 7, 9}) =
            (2 * wm1Sq * (1 + w) * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (9. * mSqq2);
        result.element({6, 7, 8, 1}) = (-2 * rt * wm1Sq * (1 + w) * (-1 + rC * w)) / (3. * mSqq);
        result.element({6, 7, 9, 1}) = (2 * rt * wm1Sq * (1 + w) * (-1 + rC * w)) / (3. * mSqq);
        result.element({6, 7, 10, 10}) = (8 * rC * rt * w2m1Sq) / (3. * mSqq);
        result.element({6, 7, 10, 13}) = (-8 * rC * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({6, 8, 6, 6}) =
            (2 * wm1Sq * (1 + w) *
             (rtSq * (-1 + 3 * rCSq - 2 * rC * (-2 + w) - 4 * w) + 2 * mSqq * (-1 + rC) * (1 + w))) /
            (9. * mSqq2);
        result.element({6, 8, 6, 7}) =
            (2 * wm1Sq * (1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) /
            (9. * mSqq2);
        result.element({6, 8, 6, 8}) =
            (2 * wm1Sq * (1 + w) * (rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq) + 2 * mSqq * w2m1)) / (9. * mSqq2);
        result.element({6, 8, 6, 9}) = (-4 * (mSqq + 2 * rtSq) * (rC - w) * wm1Sq * (1 + w)) / (9. * mSqq2);
        result.element({6, 8, 7, 6}) =
            (-2 * wm1Sq * (1 + w) *
             (rtSq * (-1 + 3 * rCSq - 2 * rC * (-2 + w) - 4 * w) + 2 * mSqq * (-1 + rC) * (1 + w))) /
            (9. * mSqq2);
        result.element({6, 8, 7, 7}) =
            (-2 * wm1Sq * (1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) /
            (9. * mSqq2);
        result.element({6, 8, 7, 8}) =
            (-2 * wm1Sq * (1 + w) * (rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq) + 2 * mSqq * w2m1)) / (9. * mSqq2);
        result.element({6, 8, 7, 9}) = (4 * (mSqq + 2 * rtSq) * (rC - w) * wm1Sq * (1 + w)) / (9. * mSqq2);
        result.element({6, 8, 8, 1}) = (-2 * rt * (rC - w) * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({6, 8, 9, 1}) = (2 * rt * (rC - w) * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({6, 8, 10, 10}) = (8 * rt * w2m1Sq) / (3. * mSqq);
        result.element({6, 8, 10, 13}) = (-8 * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({6, 9, 6, 6}) = (-2 *
                                        (2 * mSqq2 + rtSq * (-3 + rCSq + 2 * rC * (-2 + w) + 4 * w) +
                                         mSqq * (2 * rCSq + rtSq + 2 * w - 2 * rC * (1 + w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({6, 9, 6, 7}) =
            (-2 * wm1Sq * (1 + w) * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (9. * mSqq2);
        result.element({6, 9, 6, 8}) = (-4 * (mSqq + 2 * rtSq) * (rC - w) * wm1Sq * (1 + w)) / (9. * mSqq2);
        result.element({6, 9, 6, 9}) =
            (2 * (-1 + w) *
             (4 * mSqq2 + 2 * mSqq * (rCSq + rtSq - 2 * rC * w + wSq) + rtSq * (-3 + rCSq - 2 * rC * w + 4 * wSq))) /
            (9. * mSqq2);
        result.element({6, 9, 7, 6}) = (2 *
                                        (2 * mSqq2 + rtSq * (-3 + rCSq + 2 * rC * (-2 + w) + 4 * w) +
                                         mSqq * (2 * rCSq + rtSq + 2 * w - 2 * rC * (1 + w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({6, 9, 7, 7}) =
            (2 * wm1Sq * (1 + w) * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (9. * mSqq2);
        result.element({6, 9, 7, 8}) = (4 * (mSqq + 2 * rtSq) * (rC - w) * wm1Sq * (1 + w)) / (9. * mSqq2);
        result.element({6, 9, 7, 9}) =
            (-2 * (-1 + w) *
             (4 * mSqq2 + 2 * mSqq * (rCSq + rtSq - 2 * rC * w + wSq) + rtSq * (-3 + rCSq - 2 * rC * w + 4 * wSq))) /
            (9. * mSqq2);
        result.element({6, 9, 8, 1}) = (2 * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({6, 9, 9, 1}) = (-2 * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({6, 9, 10, 10}) = (8 * rt * (1 - 2 * rC + w) * w2m1) / (3. * mSqq);
        result.element({6, 9, 10, 11}) = (-8 * rC * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({6, 9, 10, 12}) = (-8 * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({6, 9, 10, 13}) = (8 * rt * (1 + 3 * rC - 2 * w) * (-1 + w)) / (3. * mSqq);
        result.element({6, 9, 10, 14}) = (8 * rt * w2m1) / (3. * mSqq);
        result.element({7, 2, 6, 2}) = (4 *
                                        (2 * mSqq2 * (-1 + w) + mSqq * (1 + 2 * rC + rCSq + rtSq) * (-1 + w) +
                                         rtSq * (1 + 2 * w - 2 * rC * (2 + w) + rCSq * (1 + 2 * w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({7, 2, 6, 3}) =
            (2 * (-1 + w) * wp1Sq *
             (2 * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2 * rC * (2 + w) + rCSq * (-1 + 4 * w)))) /
            (9. * mSqq2);
        result.element({7, 2, 6, 4}) =
            (2 * (-1 + w) * wp1Sq *
             (2 * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3 * rCSq + 4 * w - 2 * rC * (2 + w)))) /
            (9. * mSqq2);
        result.element({7, 2, 6, 5}) = (-2 *
                                        (2 * mSqq2 + mSqq * (2 * rCSq + rtSq - 2 * rC * (-1 + w) - 2 * w) +
                                         rtSq * (-3 + rCSq - 4 * w + 2 * rC * (2 + w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({7, 2, 7, 2}) = (4 *
                                        (2 * mSqq2 * (-1 + w) + mSqq * (1 + 2 * rC + rCSq + rtSq) * (-1 + w) +
                                         rtSq * (1 + 2 * w - 2 * rC * (2 + w) + rCSq * (1 + 2 * w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({7, 2, 7, 3}) =
            (2 * (-1 + w) * wp1Sq *
             (2 * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2 * rC * (2 + w) + rCSq * (-1 + 4 * w)))) /
            (9. * mSqq2);
        result.element({7, 2, 7, 4}) =
            (2 * (-1 + w) * wp1Sq *
             (2 * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3 * rCSq + 4 * w - 2 * rC * (2 + w)))) /
            (9. * mSqq2);
        result.element({7, 2, 7, 5}) = (-2 *
                                        (2 * mSqq2 + mSqq * (2 * rCSq + rtSq - 2 * rC * (-1 + w) - 2 * w) +
                                         rtSq * (-3 + rCSq - 4 * w + 2 * rC * (2 + w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({7, 2, 8, 0}) = (2 * (-1 + rC) * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({7, 2, 9, 0}) = (2 * (-1 + rC) * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({7, 2, 10, 10}) = (-8 * (1 + rC) * rt * wm1Sq * (1 + w)) / mSqq;
        result.element({7, 2, 10, 11}) = (-8 * rt * wm1Sq * (1 + w) * (-3 + rC * (-1 + 2 * w))) / (3. * mSqq);
        result.element({7, 2, 10, 12}) = (-8 * rt * (1 + 3 * rC - 2 * w) * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({7, 2, 10, 13}) = (-16 * rt * (1 + rC - w) * w2m1) / (3. * mSqq);
        result.element({7, 2, 10, 14}) = (-8 * rt * (1 + 2 * rC - w) * w2m1) / (3. * mSqq);
        result.element({7, 2, 10, 15}) = (8 * (1 + rC) * rt * w2m1Sq) / (3. * mSqq);
        result.element({7, 3, 6, 2}) =
            (2 * (-1 + w) * wp1Sq *
             (2 * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2 * rC * (2 + w) + rCSq * (-1 + 4 * w)))) /
            (9. * mSqq2);
        result.element({7, 3, 6, 3}) =
            (2 * (-1 + w) * wp1Sq * (rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)) + 2 * mSqq * rCSq * w2m1)) /
            (9. * mSqq2);
        result.element({7, 3, 6, 4}) =
            (2 * (-1 + w) * wp1Sq * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) /
            (9. * mSqq2);
        result.element({7, 3, 6, 5}) =
            (-2 * (-1 + w) * wp1Sq * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (9. * mSqq2);
        result.element({7, 3, 7, 2}) =
            (2 * (-1 + w) * wp1Sq *
             (2 * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2 * rC * (2 + w) + rCSq * (-1 + 4 * w)))) /
            (9. * mSqq2);
        result.element({7, 3, 7, 3}) =
            (2 * (-1 + w) * wp1Sq * (rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)) + 2 * mSqq * rCSq * w2m1)) /
            (9. * mSqq2);
        result.element({7, 3, 7, 4}) =
            (2 * (-1 + w) * wp1Sq * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) /
            (9. * mSqq2);
        result.element({7, 3, 7, 5}) =
            (-2 * (-1 + w) * wp1Sq * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (9. * mSqq2);
        result.element({7, 3, 8, 0}) = (2 * rt * (-1 + w) * wp1Sq * (-1 + rC * w)) / (3. * mSqq);
        result.element({7, 3, 9, 0}) = (2 * rt * (-1 + w) * wp1Sq * (-1 + rC * w)) / (3. * mSqq);
        result.element({7, 3, 10, 10}) = (-8 * rC * rt * w2m1Sq) / (3. * mSqq);
        result.element({7, 3, 10, 11}) = (8 * rC * rt * w2m1Sq) / (3. * mSqq);
        result.element({7, 3, 10, 12}) = (-8 * rC * rt * w2m1Sq) / (3. * mSqq);
        result.element({7, 3, 10, 13}) = (-8 * rC * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({7, 3, 10, 14}) = (-8 * rC * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({7, 3, 10, 15}) = (8 * rC * rt * wm1Sq * wp1Cu) / (3. * mSqq);
        result.element({7, 4, 6, 2}) =
            (2 * (-1 + w) * wp1Sq *
             (2 * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3 * rCSq + 4 * w - 2 * rC * (2 + w)))) /
            (9. * mSqq2);
        result.element({7, 4, 6, 3}) =
            (2 * (-1 + w) * wp1Sq * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) /
            (9. * mSqq2);
        result.element({7, 4, 6, 4}) =
            (2 * (-1 + w) * wp1Sq * (rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq) + 2 * mSqq * w2m1)) / (9. * mSqq2);
        result.element({7, 4, 6, 5}) = (-4 * (mSqq + 2 * rtSq) * (rC - w) * (-1 + w) * wp1Sq) / (9. * mSqq2);
        result.element({7, 4, 7, 2}) =
            (2 * (-1 + w) * wp1Sq *
             (2 * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3 * rCSq + 4 * w - 2 * rC * (2 + w)))) /
            (9. * mSqq2);
        result.element({7, 4, 7, 3}) =
            (2 * (-1 + w) * wp1Sq * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) /
            (9. * mSqq2);
        result.element({7, 4, 7, 4}) =
            (2 * (-1 + w) * wp1Sq * (rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq) + 2 * mSqq * w2m1)) / (9. * mSqq2);
        result.element({7, 4, 7, 5}) = (-4 * (mSqq + 2 * rtSq) * (rC - w) * (-1 + w) * wp1Sq) / (9. * mSqq2);
        result.element({7, 4, 8, 0}) = (2 * rt * (rC - w) * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({7, 4, 9, 0}) = (2 * rt * (rC - w) * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({7, 4, 10, 10}) = (-8 * rt * w2m1Sq) / (3. * mSqq);
        result.element({7, 4, 10, 11}) = (8 * rt * w2m1Sq) / (3. * mSqq);
        result.element({7, 4, 10, 12}) = (-8 * rt * w2m1Sq) / (3. * mSqq);
        result.element({7, 4, 10, 13}) = (-8 * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({7, 4, 10, 14}) = (-8 * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({7, 4, 10, 15}) = (8 * rt * wm1Sq * wp1Cu) / (3. * mSqq);
        result.element({7, 5, 6, 2}) = (-2 *
                                        (2 * mSqq2 + mSqq * (2 * rCSq + rtSq - 2 * rC * (-1 + w) - 2 * w) +
                                         rtSq * (-3 + rCSq - 4 * w + 2 * rC * (2 + w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({7, 5, 6, 3}) =
            (-2 * (-1 + w) * wp1Sq * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (9. * mSqq2);
        result.element({7, 5, 6, 4}) = (-4 * (mSqq + 2 * rtSq) * (rC - w) * (-1 + w) * wp1Sq) / (9. * mSqq2);
        result.element({7, 5, 6, 5}) =
            (2 * (1 + w) *
             (4 * mSqq2 + 2 * mSqq * (rCSq + rtSq - 2 * rC * w + wSq) + rtSq * (-3 + rCSq - 2 * rC * w + 4 * wSq))) /
            (9. * mSqq2);
        result.element({7, 5, 7, 2}) = (-2 *
                                        (2 * mSqq2 + mSqq * (2 * rCSq + rtSq - 2 * rC * (-1 + w) - 2 * w) +
                                         rtSq * (-3 + rCSq - 4 * w + 2 * rC * (2 + w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({7, 5, 7, 3}) =
            (-2 * (-1 + w) * wp1Sq * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (9. * mSqq2);
        result.element({7, 5, 7, 4}) = (-4 * (mSqq + 2 * rtSq) * (rC - w) * (-1 + w) * wp1Sq) / (9. * mSqq2);
        result.element({7, 5, 7, 5}) =
            (2 * (1 + w) *
             (4 * mSqq2 + 2 * mSqq * (rCSq + rtSq - 2 * rC * w + wSq) + rtSq * (-3 + rCSq - 2 * rC * w + 4 * wSq))) /
            (9. * mSqq2);
        result.element({7, 5, 8, 0}) = (-2 * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({7, 5, 9, 0}) = (-2 * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({7, 5, 10, 10}) = (8 * rt * (1 + 2 * rC - w) * w2m1) / (3. * mSqq);
        result.element({7, 5, 10, 11}) = (8 * (1 + rC) * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({7, 5, 10, 12}) = (16 * rt * (rC - w) * w2m1) / (3. * mSqq);
        result.element({7, 5, 10, 13}) = (8 * rt * (-1 + 3 * rC - 2 * w) * (1 + w)) / (3. * mSqq);
        result.element({7, 5, 10, 14}) = (8 * rt * (rC - w) * (1 + w)) / mSqq;
        result.element({7, 5, 10, 15}) = (-8 * rt * (rC - w) * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({7, 6, 6, 6}) = (-4 *
                                        (2 * mSqq2 * (1 + w) + mSqq * (1 - 2 * rC + rCSq + rtSq) * (1 + w) +
                                         rtSq * (-1 + 2 * rC * (-2 + w) + 2 * w + rCSq * (-1 + 2 * w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({7, 6, 6, 7}) =
            (-2 * wm1Sq * (1 + w) *
             (2 * mSqq * (-1 + rC) * rC * (1 + w) + rtSq * (-3 + 2 * rC * (-2 + w) + rCSq * (1 + 4 * w)))) /
            (9. * mSqq2);
        result.element({7, 6, 6, 8}) =
            (-2 * wm1Sq * (1 + w) *
             (rtSq * (-1 + 3 * rCSq - 2 * rC * (-2 + w) - 4 * w) + 2 * mSqq * (-1 + rC) * (1 + w))) /
            (9. * mSqq2);
        result.element({7, 6, 6, 9}) = (2 *
                                        (2 * mSqq2 + rtSq * (-3 + rCSq + 2 * rC * (-2 + w) + 4 * w) +
                                         mSqq * (2 * rCSq + rtSq + 2 * w - 2 * rC * (1 + w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({7, 6, 7, 6}) = (4 *
                                        (2 * mSqq2 * (1 + w) + mSqq * (1 - 2 * rC + rCSq + rtSq) * (1 + w) +
                                         rtSq * (-1 + 2 * rC * (-2 + w) + 2 * w + rCSq * (-1 + 2 * w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({7, 6, 7, 7}) =
            (2 * wm1Sq * (1 + w) *
             (2 * mSqq * (-1 + rC) * rC * (1 + w) + rtSq * (-3 + 2 * rC * (-2 + w) + rCSq * (1 + 4 * w)))) /
            (9. * mSqq2);
        result.element({7, 6, 7, 8}) =
            (2 * wm1Sq * (1 + w) *
             (rtSq * (-1 + 3 * rCSq - 2 * rC * (-2 + w) - 4 * w) + 2 * mSqq * (-1 + rC) * (1 + w))) /
            (9. * mSqq2);
        result.element({7, 6, 7, 9}) = (-2 *
                                        (2 * mSqq2 + rtSq * (-3 + rCSq + 2 * rC * (-2 + w) + 4 * w) +
                                         mSqq * (2 * rCSq + rtSq + 2 * w - 2 * rC * (1 + w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({7, 6, 8, 1}) = (2 * (1 + rC) * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({7, 6, 9, 1}) = (-2 * (1 + rC) * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({7, 6, 10, 10}) = (-8 * (-1 + rC) * rt * (-1 + w) * wp1Sq) / mSqq;
        result.element({7, 6, 10, 11}) = (-16 * rC * rt * w2m1Sq) / (3. * mSqq);
        result.element({7, 6, 10, 12}) = (-16 * rt * w2m1Sq) / (3. * mSqq);
        result.element({7, 6, 10, 13}) = (16 * rt * (-1 + rC - w) * w2m1) / (3. * mSqq);
        result.element({7, 6, 10, 14}) = (-8 * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({7, 7, 6, 6}) =
            (-2 * wm1Sq * (1 + w) *
             (2 * mSqq * (-1 + rC) * rC * (1 + w) + rtSq * (-3 + 2 * rC * (-2 + w) + rCSq * (1 + 4 * w)))) /
            (9. * mSqq2);
        result.element({7, 7, 6, 7}) =
            (-2 * wm1Sq * (1 + w) * (rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)) + 2 * mSqq * rCSq * w2m1)) /
            (9. * mSqq2);
        result.element({7, 7, 6, 8}) =
            (-2 * wm1Sq * (1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) /
            (9. * mSqq2);
        result.element({7, 7, 6, 9}) =
            (2 * wm1Sq * (1 + w) * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (9. * mSqq2);
        result.element({7, 7, 7, 6}) =
            (2 * wm1Sq * (1 + w) *
             (2 * mSqq * (-1 + rC) * rC * (1 + w) + rtSq * (-3 + 2 * rC * (-2 + w) + rCSq * (1 + 4 * w)))) /
            (9. * mSqq2);
        result.element({7, 7, 7, 7}) =
            (2 * wm1Sq * (1 + w) * (rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)) + 2 * mSqq * rCSq * w2m1)) /
            (9. * mSqq2);
        result.element({7, 7, 7, 8}) =
            (2 * wm1Sq * (1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) /
            (9. * mSqq2);
        result.element({7, 7, 7, 9}) =
            (-2 * wm1Sq * (1 + w) * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (9. * mSqq2);
        result.element({7, 7, 8, 1}) = (2 * rt * wm1Sq * (1 + w) * (-1 + rC * w)) / (3. * mSqq);
        result.element({7, 7, 9, 1}) = (-2 * rt * wm1Sq * (1 + w) * (-1 + rC * w)) / (3. * mSqq);
        result.element({7, 7, 10, 10}) = (-8 * rC * rt * w2m1Sq) / (3. * mSqq);
        result.element({7, 7, 10, 13}) = (8 * rC * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({7, 8, 6, 6}) =
            (-2 * wm1Sq * (1 + w) *
             (rtSq * (-1 + 3 * rCSq - 2 * rC * (-2 + w) - 4 * w) + 2 * mSqq * (-1 + rC) * (1 + w))) /
            (9. * mSqq2);
        result.element({7, 8, 6, 7}) =
            (-2 * wm1Sq * (1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) /
            (9. * mSqq2);
        result.element({7, 8, 6, 8}) =
            (-2 * wm1Sq * (1 + w) * (rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq) + 2 * mSqq * w2m1)) / (9. * mSqq2);
        result.element({7, 8, 6, 9}) = (4 * (mSqq + 2 * rtSq) * (rC - w) * wm1Sq * (1 + w)) / (9. * mSqq2);
        result.element({7, 8, 7, 6}) =
            (2 * wm1Sq * (1 + w) *
             (rtSq * (-1 + 3 * rCSq - 2 * rC * (-2 + w) - 4 * w) + 2 * mSqq * (-1 + rC) * (1 + w))) /
            (9. * mSqq2);
        result.element({7, 8, 7, 7}) =
            (2 * wm1Sq * (1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) /
            (9. * mSqq2);
        result.element({7, 8, 7, 8}) =
            (2 * wm1Sq * (1 + w) * (rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq) + 2 * mSqq * w2m1)) / (9. * mSqq2);
        result.element({7, 8, 7, 9}) = (-4 * (mSqq + 2 * rtSq) * (rC - w) * wm1Sq * (1 + w)) / (9. * mSqq2);
        result.element({7, 8, 8, 1}) = (2 * rt * (rC - w) * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({7, 8, 9, 1}) = (-2 * rt * (rC - w) * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({7, 8, 10, 10}) = (-8 * rt * w2m1Sq) / (3. * mSqq);
        result.element({7, 8, 10, 13}) = (8 * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({7, 9, 6, 6}) = (2 *
                                        (2 * mSqq2 + rtSq * (-3 + rCSq + 2 * rC * (-2 + w) + 4 * w) +
                                         mSqq * (2 * rCSq + rtSq + 2 * w - 2 * rC * (1 + w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({7, 9, 6, 7}) =
            (2 * wm1Sq * (1 + w) * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (9. * mSqq2);
        result.element({7, 9, 6, 8}) = (4 * (mSqq + 2 * rtSq) * (rC - w) * wm1Sq * (1 + w)) / (9. * mSqq2);
        result.element({7, 9, 6, 9}) =
            (-2 * (-1 + w) *
             (4 * mSqq2 + 2 * mSqq * (rCSq + rtSq - 2 * rC * w + wSq) + rtSq * (-3 + rCSq - 2 * rC * w + 4 * wSq))) /
            (9. * mSqq2);
        result.element({7, 9, 7, 6}) = (-2 *
                                        (2 * mSqq2 + rtSq * (-3 + rCSq + 2 * rC * (-2 + w) + 4 * w) +
                                         mSqq * (2 * rCSq + rtSq + 2 * w - 2 * rC * (1 + w))) *
                                        w2m1) /
                                       (9. * mSqq2);
        result.element({7, 9, 7, 7}) =
            (-2 * wm1Sq * (1 + w) * (2 * mSqq * rC * (rC - w) + rtSq * (-3 + rCSq + 2 * rC * w))) / (9. * mSqq2);
        result.element({7, 9, 7, 8}) = (-4 * (mSqq + 2 * rtSq) * (rC - w) * wm1Sq * (1 + w)) / (9. * mSqq2);
        result.element({7, 9, 7, 9}) =
            (2 * (-1 + w) *
             (4 * mSqq2 + 2 * mSqq * (rCSq + rtSq - 2 * rC * w + wSq) + rtSq * (-3 + rCSq - 2 * rC * w + 4 * wSq))) /
            (9. * mSqq2);
        result.element({7, 9, 8, 1}) = (-2 * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({7, 9, 9, 1}) = (2 * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({7, 9, 10, 10}) = (-8 * rt * (1 - 2 * rC + w) * w2m1) / (3. * mSqq);
        result.element({7, 9, 10, 11}) = (8 * rC * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({7, 9, 10, 12}) = (8 * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({7, 9, 10, 13}) = (-8 * rt * (1 + 3 * rC - 2 * w) * (-1 + w)) / (3. * mSqq);
        result.element({7, 9, 10, 14}) = (-8 * rt * w2m1) / (3. * mSqq);
        result.element({8, 0, 6, 2}) = (2 * (-1 + rC) * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({8, 0, 6, 3}) = (2 * rt * (-1 + w) * wp1Sq * (-1 + rC * w)) / (3. * mSqq);
        result.element({8, 0, 6, 4}) = (2 * rt * (rC - w) * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({8, 0, 6, 5}) = (-2 * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({8, 0, 7, 2}) = (2 * (-1 + rC) * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({8, 0, 7, 3}) = (2 * rt * (-1 + w) * wp1Sq * (-1 + rC * w)) / (3. * mSqq);
        result.element({8, 0, 7, 4}) = (2 * rt * (rC - w) * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({8, 0, 7, 5}) = (-2 * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({8, 0, 8, 0}) = (2 * (-1 + w) * wp1Sq) / 3.;
        result.element({8, 0, 9, 0}) = (2 * (-1 + w) * wp1Sq) / 3.;
        result.element({8, 1, 6, 6}) = (-2 * (1 + rC) * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({8, 1, 6, 7}) = (-2 * rt * wm1Sq * (1 + w) * (-1 + rC * w)) / (3. * mSqq);
        result.element({8, 1, 6, 8}) = (-2 * rt * (rC - w) * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({8, 1, 6, 9}) = (2 * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({8, 1, 7, 6}) = (2 * (1 + rC) * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({8, 1, 7, 7}) = (2 * rt * wm1Sq * (1 + w) * (-1 + rC * w)) / (3. * mSqq);
        result.element({8, 1, 7, 8}) = (2 * rt * (rC - w) * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({8, 1, 7, 9}) = (-2 * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({8, 1, 8, 1}) = (2 * wm1Sq * (1 + w)) / 3.;
        result.element({8, 1, 9, 1}) = (-2 * wm1Sq * (1 + w)) / 3.;
        result.element({9, 0, 6, 2}) = (2 * (-1 + rC) * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({9, 0, 6, 3}) = (2 * rt * (-1 + w) * wp1Sq * (-1 + rC * w)) / (3. * mSqq);
        result.element({9, 0, 6, 4}) = (2 * rt * (rC - w) * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({9, 0, 6, 5}) = (-2 * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({9, 0, 7, 2}) = (2 * (-1 + rC) * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({9, 0, 7, 3}) = (2 * rt * (-1 + w) * wp1Sq * (-1 + rC * w)) / (3. * mSqq);
        result.element({9, 0, 7, 4}) = (2 * rt * (rC - w) * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({9, 0, 7, 5}) = (-2 * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({9, 0, 8, 0}) = (2 * (-1 + w) * wp1Sq) / 3.;
        result.element({9, 0, 9, 0}) = (2 * (-1 + w) * wp1Sq) / 3.;
        result.element({9, 1, 6, 6}) = (2 * (1 + rC) * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({9, 1, 6, 7}) = (2 * rt * wm1Sq * (1 + w) * (-1 + rC * w)) / (3. * mSqq);
        result.element({9, 1, 6, 8}) = (2 * rt * (rC - w) * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({9, 1, 6, 9}) = (-2 * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({9, 1, 7, 6}) = (-2 * (1 + rC) * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({9, 1, 7, 7}) = (-2 * rt * wm1Sq * (1 + w) * (-1 + rC * w)) / (3. * mSqq);
        result.element({9, 1, 7, 8}) = (-2 * rt * (rC - w) * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({9, 1, 7, 9}) = (2 * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({9, 1, 8, 1}) = (-2 * wm1Sq * (1 + w)) / 3.;
        result.element({9, 1, 9, 1}) = (2 * wm1Sq * (1 + w)) / 3.;
        result.element({10, 10, 6, 2}) = (-8 * (1 + rC) * rt * wm1Sq * (1 + w)) / mSqq;
        result.element({10, 10, 6, 3}) = (-8 * rC * rt * w2m1Sq) / (3. * mSqq);
        result.element({10, 10, 6, 4}) = (-8 * rt * w2m1Sq) / (3. * mSqq);
        result.element({10, 10, 6, 5}) = (8 * rt * (1 + 2 * rC - w) * w2m1) / (3. * mSqq);
        result.element({10, 10, 6, 6}) = (8 * (-1 + rC) * rt * (-1 + w) * wp1Sq) / mSqq;
        result.element({10, 10, 6, 7}) = (8 * rC * rt * w2m1Sq) / (3. * mSqq);
        result.element({10, 10, 6, 8}) = (8 * rt * w2m1Sq) / (3. * mSqq);
        result.element({10, 10, 6, 9}) = (8 * rt * (1 - 2 * rC + w) * w2m1) / (3. * mSqq);
        result.element({10, 10, 7, 2}) = (-8 * (1 + rC) * rt * wm1Sq * (1 + w)) / mSqq;
        result.element({10, 10, 7, 3}) = (-8 * rC * rt * w2m1Sq) / (3. * mSqq);
        result.element({10, 10, 7, 4}) = (-8 * rt * w2m1Sq) / (3. * mSqq);
        result.element({10, 10, 7, 5}) = (8 * rt * (1 + 2 * rC - w) * w2m1) / (3. * mSqq);
        result.element({10, 10, 7, 6}) = (-8 * (-1 + rC) * rt * (-1 + w) * wp1Sq) / mSqq;
        result.element({10, 10, 7, 7}) = (-8 * rC * rt * w2m1Sq) / (3. * mSqq);
        result.element({10, 10, 7, 8}) = (-8 * rt * w2m1Sq) / (3. * mSqq);
        result.element({10, 10, 7, 9}) = (-8 * rt * (1 - 2 * rC + w) * w2m1) / (3. * mSqq);
        result.element({10, 10, 10, 10}) =
            (64 * (mSqq + 2 * rtSq) * (-4 * rC + (2 + mSqq) * w + 2 * rCSq * w) * w2m1) / (9. * mSqq2);
        result.element({10, 10, 10, 11}) =
            (-32 * (mSqq + 2 * rtSq) * wm1Sq * (1 + w) * (2 + mSqq + 4 * rC - 2 * rCSq * (1 + 2 * w))) / (9. * mSqq2);
        result.element({10, 10, 10, 12}) =
            (32 * (mSqq + 2 * rtSq) * (-2 + mSqq + 4 * rC + 2 * rCSq - 4 * w) * wm1Sq * (1 + w)) / (9. * mSqq2);
        result.element({10, 10, 10, 13}) = (128 * (mSqq + 2 * rtSq) * (rC - w) * w2m1) / (9. * mSqq2);
        result.element({10, 10, 10, 14}) =
            (32 * (mSqq + 2 * rtSq) * (-1 + mSqq + 2 * rC + rCSq - 2 * w) * w2m1) / (9. * mSqq2);
        result.element({10, 10, 10, 15}) = (-32 * (mSqq + 2 * rtSq) * w2m1Sq) / (9. * mSqq);
        result.element({10, 11, 6, 2}) = (-8 * rt * wm1Sq * (1 + w) * (-3 + rC * (-1 + 2 * w))) / (3. * mSqq);
        result.element({10, 11, 6, 3}) = (8 * rC * rt * w2m1Sq) / (3. * mSqq);
        result.element({10, 11, 6, 4}) = (8 * rt * w2m1Sq) / (3. * mSqq);
        result.element({10, 11, 6, 5}) = (8 * (1 + rC) * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({10, 11, 6, 6}) = (16 * rC * rt * w2m1Sq) / (3. * mSqq);
        result.element({10, 11, 6, 9}) = (-8 * rC * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({10, 11, 7, 2}) = (-8 * rt * wm1Sq * (1 + w) * (-3 + rC * (-1 + 2 * w))) / (3. * mSqq);
        result.element({10, 11, 7, 3}) = (8 * rC * rt * w2m1Sq) / (3. * mSqq);
        result.element({10, 11, 7, 4}) = (8 * rt * w2m1Sq) / (3. * mSqq);
        result.element({10, 11, 7, 5}) = (8 * (1 + rC) * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({10, 11, 7, 6}) = (-16 * rC * rt * w2m1Sq) / (3. * mSqq);
        result.element({10, 11, 7, 9}) = (8 * rC * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({10, 11, 10, 10}) =
            (-32 * (mSqq + 2 * rtSq) * wm1Sq * (1 + w) * (2 + mSqq + 4 * rC - 2 * rCSq * (1 + 2 * w))) / (9. * mSqq2);
        result.element({10, 11, 10, 11}) = (32 * wm1Sq * (1 + w) *
                                            (mSqq2 + 4 * rtSq * (1 - 2 * rC * w + rCSq * (-1 + 2 * wSq)) +
                                             2 * mSqq * (1 + rtSq - 2 * rC * w + rCSq * (-1 + 2 * wSq)))) /
                                           (9. * mSqq2);
        result.element({10, 11, 10, 12}) =
            (-32 * (mSqq + 2 * rtSq) * wm1Sq * (1 + w) * (mSqq - 2 * (-2 * rC + w + rCSq * w))) / (9. * mSqq2);
        result.element({10, 11, 10, 13}) =
            (-32 * (mSqq + 2 * rtSq) * (1 + mSqq + 2 * rC - rCSq - 2 * w) * w2m1) / (9. * mSqq2);
        result.element({10, 11, 10, 14}) =
            (-32 * (mSqq + 2 * rtSq) * (mSqq + 2 * rC - w - rCSq * w) * w2m1) / (9. * mSqq2);
        result.element({10, 11, 10, 15}) = (32 * (mSqq + 2 * rtSq) * w2m1Sq) / (9. * mSqq);
        result.element({10, 12, 6, 2}) = (-8 * rt * (1 + 3 * rC - 2 * w) * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({10, 12, 6, 3}) = (-8 * rC * rt * w2m1Sq) / (3. * mSqq);
        result.element({10, 12, 6, 4}) = (-8 * rt * w2m1Sq) / (3. * mSqq);
        result.element({10, 12, 6, 5}) = (16 * rt * (rC - w) * w2m1) / (3. * mSqq);
        result.element({10, 12, 6, 6}) = (16 * rt * w2m1Sq) / (3. * mSqq);
        result.element({10, 12, 6, 9}) = (-8 * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({10, 12, 7, 2}) = (-8 * rt * (1 + 3 * rC - 2 * w) * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({10, 12, 7, 3}) = (-8 * rC * rt * w2m1Sq) / (3. * mSqq);
        result.element({10, 12, 7, 4}) = (-8 * rt * w2m1Sq) / (3. * mSqq);
        result.element({10, 12, 7, 5}) = (16 * rt * (rC - w) * w2m1) / (3. * mSqq);
        result.element({10, 12, 7, 6}) = (-16 * rt * w2m1Sq) / (3. * mSqq);
        result.element({10, 12, 7, 9}) = (8 * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({10, 12, 10, 10}) =
            (32 * (mSqq + 2 * rtSq) * (-2 + mSqq + 4 * rC + 2 * rCSq - 4 * w) * wm1Sq * (1 + w)) / (9. * mSqq2);
        result.element({10, 12, 10, 11}) =
            (-32 * (mSqq + 2 * rtSq) * wm1Sq * (1 + w) * (mSqq - 2 * (-2 * rC + w + rCSq * w))) / (9. * mSqq2);
        result.element({10, 12, 10, 12}) = (32 * wm1Sq * (1 + w) *
                                            (mSqq2 + 4 * rtSq * (-1 + rCSq - 2 * rC * w + 2 * wSq) +
                                             2 * mSqq * (-1 + rCSq + rtSq - 2 * rC * w + 2 * wSq))) /
                                           (9. * mSqq2);
        result.element({10, 12, 10, 13}) =
            (32 * (mSqq + 2 * rtSq) * (-1 + mSqq + rCSq + rC * (2 - 4 * w) - 2 * w + 4 * wSq) * w2m1) / (9. * mSqq2);
        result.element({10, 12, 10, 14}) = (32 *
                                            (mSqq2 + 2 * rtSq * (-1 + rCSq - 2 * rC * w + 2 * wSq) +
                                             mSqq * (-1 + rCSq + 2 * rtSq - 2 * rC * w + 2 * wSq)) *
                                            w2m1) /
                                           (9. * mSqq2);
        result.element({10, 12, 10, 15}) = (-32 * (mSqq + 2 * rtSq) * w2m1Sq) / (9. * mSqq);
        result.element({10, 13, 6, 2}) = (-16 * rt * (1 + rC - w) * w2m1) / (3. * mSqq);
        result.element({10, 13, 6, 3}) = (-8 * rC * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({10, 13, 6, 4}) = (-8 * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({10, 13, 6, 5}) = (8 * rt * (-1 + 3 * rC - 2 * w) * (1 + w)) / (3. * mSqq);
        result.element({10, 13, 6, 6}) = (8 * rt * (2 - 2 * rC + 2 * w) * w2m1) / (3. * mSqq);
        result.element({10, 13, 6, 7}) = (-8 * rC * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({10, 13, 6, 8}) = (-8 * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({10, 13, 6, 9}) = (8 * rt * (1 + 3 * rC - 2 * w) * (-1 + w)) / (3. * mSqq);
        result.element({10, 13, 7, 2}) = (-16 * rt * (1 + rC - w) * w2m1) / (3. * mSqq);
        result.element({10, 13, 7, 3}) = (-8 * rC * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({10, 13, 7, 4}) = (-8 * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({10, 13, 7, 5}) = (8 * rt * (-1 + 3 * rC - 2 * w) * (1 + w)) / (3. * mSqq);
        result.element({10, 13, 7, 6}) = (16 * rt * (-1 + rC - w) * w2m1) / (3. * mSqq);
        result.element({10, 13, 7, 7}) = (8 * rC * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({10, 13, 7, 8}) = (8 * rt * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({10, 13, 7, 9}) = (-8 * rt * (1 + 3 * rC - 2 * w) * (-1 + w)) / (3. * mSqq);
        result.element({10, 13, 10, 10}) = (128 * (mSqq + 2 * rtSq) * (rC - w) * w2m1) / (9. * mSqq2);
        result.element({10, 13, 10, 11}) =
            (-32 * (mSqq + 2 * rtSq) * (1 + mSqq + 2 * rC - rCSq - 2 * w) * w2m1) / (9. * mSqq2);
        result.element({10, 13, 10, 12}) =
            (32 * (mSqq + 2 * rtSq) * (-1 + mSqq + rCSq + rC * (2 - 4 * w) - 2 * w + 4 * wSq) * w2m1) / (9. * mSqq2);
        result.element({10, 13, 10, 13}) = (64 * (mSqq2 * w + 4 * rtSq * (rCSq * w + wCu - rC * (1 + wSq)) +
                                                  2 * mSqq * (rCSq * w - rC * (1 + wSq) + w * (rtSq + wSq)))) /
                                           (9. * mSqq2);
        result.element({10, 13, 10, 14}) =
            (32 * (mSqq + 2 * rtSq) * (1 + w) * (mSqq + 2 * (rCSq + w - rC * (1 + w) + w2m1))) / (9. * mSqq2);
        result.element({10, 13, 10, 15}) = (-32 * (mSqq + 2 * rtSq) * (-1 + w) * wp1Sq) / (9. * mSqq);
        result.element({10, 14, 6, 2}) = (-8 * rt * (1 + 2 * rC - w) * w2m1) / (3. * mSqq);
        result.element({10, 14, 6, 3}) = (-8 * rC * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({10, 14, 6, 4}) = (-8 * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({10, 14, 6, 5}) = (8 * rt * (rC - w) * (1 + w)) / mSqq;
        result.element({10, 14, 6, 6}) = (8 * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({10, 14, 6, 9}) = (8 * rt * w2m1) / (3. * mSqq);
        result.element({10, 14, 7, 2}) = (-8 * rt * (1 + 2 * rC - w) * w2m1) / (3. * mSqq);
        result.element({10, 14, 7, 3}) = (-8 * rC * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({10, 14, 7, 4}) = (-8 * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({10, 14, 7, 5}) = (8 * rt * (rC - w) * (1 + w)) / mSqq;
        result.element({10, 14, 7, 6}) = (-8 * rt * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({10, 14, 7, 9}) = (-8 * rt * w2m1) / (3. * mSqq);
        result.element({10, 14, 10, 10}) =
            (32 * (mSqq + 2 * rtSq) * (-1 + mSqq + 2 * rC + rCSq - 2 * w) * w2m1) / (9. * mSqq2);
        result.element({10, 14, 10, 11}) =
            (-32 * (mSqq + 2 * rtSq) * (mSqq + 2 * rC - w - rCSq * w) * w2m1) / (9. * mSqq2);
        result.element({10, 14, 10, 12}) = (32 *
                                            (mSqq2 + 2 * rtSq * (-1 + rCSq - 2 * rC * w + 2 * wSq) +
                                             mSqq * (-1 + rCSq + 2 * rtSq - 2 * rC * w + 2 * wSq)) *
                                            w2m1) /
                                           (9. * mSqq2);
        result.element({10, 14, 10, 13}) =
            (32 * (mSqq + 2 * rtSq) * (1 + w) * (mSqq + 2 * (rCSq + w - rC * (1 + w) + w2m1))) / (9. * mSqq2);
        result.element({10, 14, 10, 14}) = (32 * (1 + w) *
                                            (mSqq2 + 4 * rtSq * (-1 + rCSq - 2 * rC * w + 2 * wSq) +
                                             2 * mSqq * (-1 + rCSq + rtSq - 2 * rC * w + 2 * wSq))) /
                                           (9. * mSqq2);
        result.element({10, 14, 10, 15}) = (-32 * (mSqq + 2 * rtSq) * (-1 + w) * wp1Sq) / (9. * mSqq);
        result.element({10, 15, 6, 2}) = (8 * (1 + rC) * rt * w2m1Sq) / (3. * mSqq);
        result.element({10, 15, 6, 3}) = (8 * rC * rt * wm1Sq * wp1Cu) / (3. * mSqq);
        result.element({10, 15, 6, 4}) = (8 * rt * wm1Sq * wp1Cu) / (3. * mSqq);
        result.element({10, 15, 6, 5}) = (-8 * rt * (rC - w) * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({10, 15, 7, 2}) = (8 * (1 + rC) * rt * w2m1Sq) / (3. * mSqq);
        result.element({10, 15, 7, 3}) = (8 * rC * rt * wm1Sq * wp1Cu) / (3. * mSqq);
        result.element({10, 15, 7, 4}) = (8 * rt * wm1Sq * wp1Cu) / (3. * mSqq);
        result.element({10, 15, 7, 5}) = (-8 * rt * (rC - w) * (-1 + w) * wp1Sq) / (3. * mSqq);
        result.element({10, 15, 10, 10}) = (-32 * (mSqq + 2 * rtSq) * w2m1Sq) / (9. * mSqq);
        result.element({10, 15, 10, 11}) = (32 * (mSqq + 2 * rtSq) * w2m1Sq) / (9. * mSqq);
        result.element({10, 15, 10, 12}) = (-32 * (mSqq + 2 * rtSq) * w2m1Sq) / (9. * mSqq);
        result.element({10, 15, 10, 13}) = (-32 * (mSqq + 2 * rtSq) * (-1 + w) * wp1Sq) / (9. * mSqq);
        result.element({10, 15, 10, 14}) = (-32 * (mSqq + 2 * rtSq) * (-1 + w) * wp1Sq) / (9. * mSqq);
        result.element({10, 15, 10, 15}) = (32 * (mSqq + 2 * rtSq) * wm1Sq * wp1Cu) / (9. * mSqq);

        result *= RateNorm;

        return result;
    }
    // NOLINTEND(readability-function-size)

} // namespace Hammer
