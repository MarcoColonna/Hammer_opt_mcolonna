///
/// @file  RateLbLcstar12LepNu.cc
/// @brief \f$ \Lambda_b \rightarrow L_c^*(2595) \tau\nu \f$ total rate
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/Rates/RateLbLcstar12LepNu.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    RateLbLcstar12LepNu::RateLbLcstar12LepNu() {
        // Create tensor rank and dimensions
        IndexList dims{{11, 12, 11, 12, _nPoints}};
        string name{"RateLbLcstar12LepNuQ2"};
        auto& pdg = PID::instance();

        addProcessSignature(-PID::LAMBDAB, {PID::LAMBDACSTAR12MINUS, PID::NU_TAU, PID::ANTITAU});
        addIntegrationBoundaries({PS::makeQ2Function(
            pdg.getMass(PID::ANTITAU), pdg.getMass(PID::LAMBDAB) - pdg.getMass(PID::LAMBDACSTAR12MINUS))});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCTAUNU, FF_LBLCSTAR12, WILSON_BCTAUNU_HC,
                                                          FF_LBLCSTAR12_HC, INTEGRATION_INDEX})});

        addProcessSignature(-PID::LAMBDAB, {PID::LAMBDACSTAR12MINUS, PID::NU_MU, PID::ANTIMUON});
        addIntegrationBoundaries({PS::makeQ2Function(
            pdg.getMass(PID::ANTIMUON), pdg.getMass(PID::LAMBDAB) - pdg.getMass(PID::LAMBDACSTAR12MINUS))});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCMUNU, FF_LBLCSTAR12, WILSON_BCMUNU_HC,
                                                          FF_LBLCSTAR12_HC, INTEGRATION_INDEX})});

        addProcessSignature(-PID::LAMBDAB, {PID::LAMBDACSTAR12MINUS, PID::NU_E, PID::POSITRON});
        addIntegrationBoundaries({PS::makeQ2Function(
            pdg.getMass(PID::POSITRON), pdg.getMass(PID::LAMBDAB) - pdg.getMass(PID::LAMBDACSTAR12MINUS))});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCENU, FF_LBLCSTAR12, WILSON_BCENU_HC,
                                                          FF_LBLCSTAR12_HC, INTEGRATION_INDEX})});

        setSignatureIndex();
    }

    Tensor RateLbLcstar12LepNu::evalAtPSPoint(const vector<double>& point) {
        auto labs = getTensor().labels();
        labs.pop_back();
        auto dimensions = getTensor().dims();
        dimensions.pop_back();
        Tensor result{"RateLbLcstar12LepNu", MD::makeEmptySparse(dimensions, labs)};

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
        // const double wp1Sq = (w+1)*(w+1);
        const double rCSq = rC * rC;
        const double rtSq = rt * rt;

        const double RateNorm =
            (pow(GFermi, 2.) * pow(Mb, 3.) * pow(rC, 2.) * pow(mSqq - rtSq, 2.) * sqrt(-1 + wSq)) / (32. * mSqq * pi3);

        // set non-zero tensor elements
        result.element({0, 2, 0, 2}) = (2 * (2 * mSqq2 * (1 + w) + mSqq * (1 - 2 * rC + rCSq + rtSq) * (1 + w) +
                                             rtSq * (-1 + 2 * rC * (-2 + w) + 2 * w + rCSq * (-1 + 2 * w)))) /
                                       (3. * mSqq2);
        result.element({0, 2, 0, 3}) =
            ((-1 + w) * (2 * mSqq * (-1 + rC) * rC * (1 + w) + rtSq * (-3 + 2 * rC * (-2 + w) + rCSq * (1 + 4 * w)))) /
            (3. * mSqq2);
        result.element({0, 2, 0, 4}) =
            ((-1 + w) * (rtSq * (-1 + 3 * rCSq - 2 * rC * (-2 + w) - 4 * w) + 2 * mSqq * (-1 + rC) * (1 + w))) /
            (3. * mSqq2);
        result.element({0, 2, 1, 2}) = (2 * (2 * mSqq2 * (1 + w) + mSqq * (1 - 2 * rC + rCSq + rtSq) * (1 + w) +
                                             rtSq * (-1 + 2 * rC * (-2 + w) + 2 * w + rCSq * (-1 + 2 * w)))) /
                                       (3. * mSqq2);
        result.element({0, 2, 1, 3}) =
            ((-1 + w) * (2 * mSqq * (-1 + rC) * rC * (1 + w) + rtSq * (-3 + 2 * rC * (-2 + w) + rCSq * (1 + 4 * w)))) /
            (3. * mSqq2);
        result.element({0, 2, 1, 4}) =
            ((-1 + w) * (rtSq * (-1 + 3 * rCSq - 2 * rC * (-2 + w) - 4 * w) + 2 * mSqq * (-1 + rC) * (1 + w))) /
            (3. * mSqq2);
        result.element({0, 2, 2, 2}) = (2 * (2 * mSqq2 * (1 + w) + mSqq * (1 - 2 * rC + rCSq + rtSq) * (1 + w) +
                                             rtSq * (-1 + 2 * rC * (-2 + w) + 2 * w + rCSq * (-1 + 2 * w)))) /
                                       (3. * mSqq2);
        result.element({0, 2, 2, 3}) =
            ((-1 + w) * (2 * mSqq * (-1 + rC) * rC * (1 + w) + rtSq * (-3 + 2 * rC * (-2 + w) + rCSq * (1 + 4 * w)))) /
            (3. * mSqq2);
        result.element({0, 2, 2, 4}) =
            ((-1 + w) * (rtSq * (-1 + 3 * rCSq - 2 * rC * (-2 + w) - 4 * w) + 2 * mSqq * (-1 + rC) * (1 + w))) /
            (3. * mSqq2);
        result.element({0, 2, 3, 0}) = -(((1 + rC) * rt * (-1 + w)) / mSqq);
        result.element({0, 2, 4, 0}) = -(((1 + rC) * rt * (-1 + w)) / mSqq);
        result.element({0, 2, 5, 8}) = (12 * (-1 + rC) * rt * (1 + w)) / mSqq;
        result.element({0, 2, 5, 9}) = (4 * rt * (1 + w) * (-3 + rC + 2 * rC * w)) / mSqq;
        result.element({0, 2, 5, 10}) = (4 * rt * (-1 + 3 * rC - 2 * w) * (1 + w)) / mSqq;
        result.element({0, 2, 5, 11}) = (-4 * (-1 + rC) * rt * w2m1) / mSqq;
        result.element({0, 3, 0, 2}) =
            ((-1 + w) * (2 * mSqq * (-1 + rC) * rC * (1 + w) + rtSq * (-3 + 2 * rC * (-2 + w) + rCSq * (1 + 4 * w)))) /
            (3. * mSqq2);
        result.element({0, 3, 0, 3}) =
            ((-1 + w) * (rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)) + 2 * mSqq * rCSq * w2m1)) / (3. * mSqq2);
        result.element({0, 3, 0, 4}) =
            ((-1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) / (3. * mSqq2);
        result.element({0, 3, 1, 2}) =
            ((-1 + w) * (2 * mSqq * (-1 + rC) * rC * (1 + w) + rtSq * (-3 + 2 * rC * (-2 + w) + rCSq * (1 + 4 * w)))) /
            (3. * mSqq2);
        result.element({0, 3, 1, 3}) =
            ((-1 + w) * (rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)) + 2 * mSqq * rCSq * w2m1)) / (3. * mSqq2);
        result.element({0, 3, 1, 4}) =
            ((-1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) / (3. * mSqq2);
        result.element({0, 3, 2, 2}) =
            ((-1 + w) * (2 * mSqq * (-1 + rC) * rC * (1 + w) + rtSq * (-3 + 2 * rC * (-2 + w) + rCSq * (1 + 4 * w)))) /
            (3. * mSqq2);
        result.element({0, 3, 2, 3}) =
            ((-1 + w) * (rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)) + 2 * mSqq * rCSq * w2m1)) / (3. * mSqq2);
        result.element({0, 3, 2, 4}) =
            ((-1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) / (3. * mSqq2);
        result.element({0, 3, 3, 0}) = -((rt * (-1 + w) * (-1 + rC * w)) / mSqq);
        result.element({0, 3, 4, 0}) = -((rt * (-1 + w) * (-1 + rC * w)) / mSqq);
        result.element({0, 3, 5, 8}) = (4 * rC * rt * w2m1) / mSqq;
        result.element({0, 3, 5, 9}) = (4 * rC * rt * w2m1) / mSqq;
        result.element({0, 3, 5, 10}) = (4 * rC * rt * w2m1) / mSqq;
        result.element({0, 3, 5, 11}) = (-4 * rC * rt * wm1Sq * (1 + w)) / mSqq;
        result.element({0, 4, 0, 2}) =
            ((-1 + w) * (rtSq * (-1 + 3 * rCSq - 2 * rC * (-2 + w) - 4 * w) + 2 * mSqq * (-1 + rC) * (1 + w))) /
            (3. * mSqq2);
        result.element({0, 4, 0, 3}) =
            ((-1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) / (3. * mSqq2);
        result.element({0, 4, 0, 4}) =
            ((-1 + w) * (rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq) + 2 * mSqq * w2m1)) / (3. * mSqq2);
        result.element({0, 4, 1, 2}) =
            ((-1 + w) * (rtSq * (-1 + 3 * rCSq - 2 * rC * (-2 + w) - 4 * w) + 2 * mSqq * (-1 + rC) * (1 + w))) /
            (3. * mSqq2);
        result.element({0, 4, 1, 3}) =
            ((-1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) / (3. * mSqq2);
        result.element({0, 4, 1, 4}) =
            ((-1 + w) * (rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq) + 2 * mSqq * w2m1)) / (3. * mSqq2);
        result.element({0, 4, 2, 2}) =
            ((-1 + w) * (rtSq * (-1 + 3 * rCSq - 2 * rC * (-2 + w) - 4 * w) + 2 * mSqq * (-1 + rC) * (1 + w))) /
            (3. * mSqq2);
        result.element({0, 4, 2, 3}) =
            ((-1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) / (3. * mSqq2);
        result.element({0, 4, 2, 4}) =
            ((-1 + w) * (rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq) + 2 * mSqq * w2m1)) / (3. * mSqq2);
        result.element({0, 4, 3, 0}) = -((rt * (rC - w) * (-1 + w)) / mSqq);
        result.element({0, 4, 4, 0}) = -((rt * (rC - w) * (-1 + w)) / mSqq);
        result.element({0, 4, 5, 8}) = (4 * rt * w2m1) / mSqq;
        result.element({0, 4, 5, 9}) = (4 * rt * w2m1) / mSqq;
        result.element({0, 4, 5, 10}) = (4 * rt * w2m1) / mSqq;
        result.element({0, 4, 5, 11}) = (-4 * rt * wm1Sq * (1 + w)) / mSqq;
        result.element({0, 5, 0, 5}) = (2 * (2 * mSqq2 * (-1 + w) + mSqq * (1 + 2 * rC + rCSq + rtSq) * (-1 + w) +
                                             rtSq * (1 + 2 * w - 2 * rC * (2 + w) + rCSq * (1 + 2 * w)))) /
                                       (3. * mSqq2);
        result.element({0, 5, 0, 6}) =
            ((1 + w) * (2 * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2 * rC * (2 + w) + rCSq * (-1 + 4 * w)))) /
            (3. * mSqq2);
        result.element({0, 5, 0, 7}) =
            ((1 + w) * (2 * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3 * rCSq + 4 * w - 2 * rC * (2 + w)))) /
            (3. * mSqq2);
        result.element({0, 5, 1, 5}) = (2 * (2 * mSqq2 * (-1 + w) + mSqq * (1 + 2 * rC + rCSq + rtSq) * (-1 + w) +
                                             rtSq * (1 + 2 * w - 2 * rC * (2 + w) + rCSq * (1 + 2 * w)))) /
                                       (3. * mSqq2);
        result.element({0, 5, 1, 6}) =
            ((1 + w) * (2 * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2 * rC * (2 + w) + rCSq * (-1 + 4 * w)))) /
            (3. * mSqq2);
        result.element({0, 5, 1, 7}) =
            ((1 + w) * (2 * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3 * rCSq + 4 * w - 2 * rC * (2 + w)))) /
            (3. * mSqq2);
        result.element({0, 5, 2, 5}) = (-2 * (2 * mSqq2 * (-1 + w) + mSqq * (1 + 2 * rC + rCSq + rtSq) * (-1 + w) +
                                              rtSq * (1 + 2 * w - 2 * rC * (2 + w) + rCSq * (1 + 2 * w)))) /
                                       (3. * mSqq2);
        result.element({0, 5, 2, 6}) =
            -((1 + w) * (2 * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2 * rC * (2 + w) + rCSq * (-1 + 4 * w)))) /
            (3. * mSqq2);
        result.element({0, 5, 2, 7}) =
            -((1 + w) * (2 * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3 * rCSq + 4 * w - 2 * rC * (2 + w)))) /
            (3. * mSqq2);
        result.element({0, 5, 3, 1}) = ((-1 + rC) * rt * (1 + w)) / mSqq;
        result.element({0, 5, 4, 1}) = -(((-1 + rC) * rt * (1 + w)) / mSqq);
        result.element({0, 5, 5, 8}) = (12 * (1 + rC) * rt * (-1 + w)) / mSqq;
        result.element({0, 5, 5, 9}) = (8 * rC * rt * w2m1) / mSqq;
        result.element({0, 5, 5, 10}) = (8 * rt * w2m1) / mSqq;
        result.element({0, 6, 0, 5}) =
            ((1 + w) * (2 * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2 * rC * (2 + w) + rCSq * (-1 + 4 * w)))) /
            (3. * mSqq2);
        result.element({0, 6, 0, 6}) =
            ((1 + w) * (rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)) + 2 * mSqq * rCSq * w2m1)) / (3. * mSqq2);
        result.element({0, 6, 0, 7}) =
            ((1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) / (3. * mSqq2);
        result.element({0, 6, 1, 5}) =
            ((1 + w) * (2 * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2 * rC * (2 + w) + rCSq * (-1 + 4 * w)))) /
            (3. * mSqq2);
        result.element({0, 6, 1, 6}) =
            ((1 + w) * (rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)) + 2 * mSqq * rCSq * w2m1)) / (3. * mSqq2);
        result.element({0, 6, 1, 7}) =
            ((1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) / (3. * mSqq2);
        result.element({0, 6, 2, 5}) =
            -((1 + w) * (2 * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2 * rC * (2 + w) + rCSq * (-1 + 4 * w)))) /
            (3. * mSqq2);
        result.element({0, 6, 2, 6}) =
            -((1 + w) * (rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)) + 2 * mSqq * rCSq * w2m1)) / (3. * mSqq2);
        result.element({0, 6, 2, 7}) =
            -((1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) / (3. * mSqq2);
        result.element({0, 6, 3, 1}) = (rt * (1 + w) * (-1 + rC * w)) / mSqq;
        result.element({0, 6, 4, 1}) = -((rt * (1 + w) * (-1 + rC * w)) / mSqq);
        result.element({0, 6, 5, 8}) = (4 * rC * rt * w2m1) / mSqq;
        result.element({0, 7, 0, 5}) =
            ((1 + w) * (2 * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3 * rCSq + 4 * w - 2 * rC * (2 + w)))) /
            (3. * mSqq2);
        result.element({0, 7, 0, 6}) =
            ((1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) / (3. * mSqq2);
        result.element({0, 7, 0, 7}) =
            ((1 + w) * (rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq) + 2 * mSqq * w2m1)) / (3. * mSqq2);
        result.element({0, 7, 1, 5}) =
            ((1 + w) * (2 * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3 * rCSq + 4 * w - 2 * rC * (2 + w)))) /
            (3. * mSqq2);
        result.element({0, 7, 1, 6}) =
            ((1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) / (3. * mSqq2);
        result.element({0, 7, 1, 7}) =
            ((1 + w) * (rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq) + 2 * mSqq * w2m1)) / (3. * mSqq2);
        result.element({0, 7, 2, 5}) =
            -((1 + w) * (2 * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3 * rCSq + 4 * w - 2 * rC * (2 + w)))) /
            (3. * mSqq2);
        result.element({0, 7, 2, 6}) =
            -((1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) / (3. * mSqq2);
        result.element({0, 7, 2, 7}) =
            -((1 + w) * (rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq) + 2 * mSqq * w2m1)) / (3. * mSqq2);
        result.element({0, 7, 3, 1}) = (rt * (rC - w) * (1 + w)) / mSqq;
        result.element({0, 7, 4, 1}) = -((rt * (rC - w) * (1 + w)) / mSqq);
        result.element({0, 7, 5, 8}) = (4 * rt * w2m1) / mSqq;
        result.element({1, 2, 0, 2}) = (2 * (2 * mSqq2 * (1 + w) + mSqq * (1 - 2 * rC + rCSq + rtSq) * (1 + w) +
                                             rtSq * (-1 + 2 * rC * (-2 + w) + 2 * w + rCSq * (-1 + 2 * w)))) /
                                       (3. * mSqq2);
        result.element({1, 2, 0, 3}) =
            ((-1 + w) * (2 * mSqq * (-1 + rC) * rC * (1 + w) + rtSq * (-3 + 2 * rC * (-2 + w) + rCSq * (1 + 4 * w)))) /
            (3. * mSqq2);
        result.element({1, 2, 0, 4}) =
            ((-1 + w) * (rtSq * (-1 + 3 * rCSq - 2 * rC * (-2 + w) - 4 * w) + 2 * mSqq * (-1 + rC) * (1 + w))) /
            (3. * mSqq2);
        result.element({1, 2, 1, 2}) = (2 * (2 * mSqq2 * (1 + w) + mSqq * (1 - 2 * rC + rCSq + rtSq) * (1 + w) +
                                             rtSq * (-1 + 2 * rC * (-2 + w) + 2 * w + rCSq * (-1 + 2 * w)))) /
                                       (3. * mSqq2);
        result.element({1, 2, 1, 3}) =
            ((-1 + w) * (2 * mSqq * (-1 + rC) * rC * (1 + w) + rtSq * (-3 + 2 * rC * (-2 + w) + rCSq * (1 + 4 * w)))) /
            (3. * mSqq2);
        result.element({1, 2, 1, 4}) =
            ((-1 + w) * (rtSq * (-1 + 3 * rCSq - 2 * rC * (-2 + w) - 4 * w) + 2 * mSqq * (-1 + rC) * (1 + w))) /
            (3. * mSqq2);
        result.element({1, 2, 2, 2}) = (2 * (2 * mSqq2 * (1 + w) + mSqq * (1 - 2 * rC + rCSq + rtSq) * (1 + w) +
                                             rtSq * (-1 + 2 * rC * (-2 + w) + 2 * w + rCSq * (-1 + 2 * w)))) /
                                       (3. * mSqq2);
        result.element({1, 2, 2, 3}) =
            ((-1 + w) * (2 * mSqq * (-1 + rC) * rC * (1 + w) + rtSq * (-3 + 2 * rC * (-2 + w) + rCSq * (1 + 4 * w)))) /
            (3. * mSqq2);
        result.element({1, 2, 2, 4}) =
            ((-1 + w) * (rtSq * (-1 + 3 * rCSq - 2 * rC * (-2 + w) - 4 * w) + 2 * mSqq * (-1 + rC) * (1 + w))) /
            (3. * mSqq2);
        result.element({1, 2, 3, 0}) = -(((1 + rC) * rt * (-1 + w)) / mSqq);
        result.element({1, 2, 4, 0}) = -(((1 + rC) * rt * (-1 + w)) / mSqq);
        result.element({1, 2, 5, 8}) = (12 * (-1 + rC) * rt * (1 + w)) / mSqq;
        result.element({1, 2, 5, 9}) = (4 * rt * (1 + w) * (-3 + rC + 2 * rC * w)) / mSqq;
        result.element({1, 2, 5, 10}) = (4 * rt * (-1 + 3 * rC - 2 * w) * (1 + w)) / mSqq;
        result.element({1, 2, 5, 11}) = (-4 * (-1 + rC) * rt * w2m1) / mSqq;
        result.element({1, 3, 0, 2}) =
            ((-1 + w) * (2 * mSqq * (-1 + rC) * rC * (1 + w) + rtSq * (-3 + 2 * rC * (-2 + w) + rCSq * (1 + 4 * w)))) /
            (3. * mSqq2);
        result.element({1, 3, 0, 3}) =
            ((-1 + w) * (rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)) + 2 * mSqq * rCSq * w2m1)) / (3. * mSqq2);
        result.element({1, 3, 0, 4}) =
            ((-1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) / (3. * mSqq2);
        result.element({1, 3, 1, 2}) =
            ((-1 + w) * (2 * mSqq * (-1 + rC) * rC * (1 + w) + rtSq * (-3 + 2 * rC * (-2 + w) + rCSq * (1 + 4 * w)))) /
            (3. * mSqq2);
        result.element({1, 3, 1, 3}) =
            ((-1 + w) * (rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)) + 2 * mSqq * rCSq * w2m1)) / (3. * mSqq2);
        result.element({1, 3, 1, 4}) =
            ((-1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) / (3. * mSqq2);
        result.element({1, 3, 2, 2}) =
            ((-1 + w) * (2 * mSqq * (-1 + rC) * rC * (1 + w) + rtSq * (-3 + 2 * rC * (-2 + w) + rCSq * (1 + 4 * w)))) /
            (3. * mSqq2);
        result.element({1, 3, 2, 3}) =
            ((-1 + w) * (rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)) + 2 * mSqq * rCSq * w2m1)) / (3. * mSqq2);
        result.element({1, 3, 2, 4}) =
            ((-1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) / (3. * mSqq2);
        result.element({1, 3, 3, 0}) = -((rt * (-1 + w) * (-1 + rC * w)) / mSqq);
        result.element({1, 3, 4, 0}) = -((rt * (-1 + w) * (-1 + rC * w)) / mSqq);
        result.element({1, 3, 5, 8}) = (4 * rC * rt * w2m1) / mSqq;
        result.element({1, 3, 5, 9}) = (4 * rC * rt * w2m1) / mSqq;
        result.element({1, 3, 5, 10}) = (4 * rC * rt * w2m1) / mSqq;
        result.element({1, 3, 5, 11}) = (-4 * rC * rt * wm1Sq * (1 + w)) / mSqq;
        result.element({1, 4, 0, 2}) =
            ((-1 + w) * (rtSq * (-1 + 3 * rCSq - 2 * rC * (-2 + w) - 4 * w) + 2 * mSqq * (-1 + rC) * (1 + w))) /
            (3. * mSqq2);
        result.element({1, 4, 0, 3}) =
            ((-1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) / (3. * mSqq2);
        result.element({1, 4, 0, 4}) =
            ((-1 + w) * (rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq) + 2 * mSqq * w2m1)) / (3. * mSqq2);
        result.element({1, 4, 1, 2}) =
            ((-1 + w) * (rtSq * (-1 + 3 * rCSq - 2 * rC * (-2 + w) - 4 * w) + 2 * mSqq * (-1 + rC) * (1 + w))) /
            (3. * mSqq2);
        result.element({1, 4, 1, 3}) =
            ((-1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) / (3. * mSqq2);
        result.element({1, 4, 1, 4}) =
            ((-1 + w) * (rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq) + 2 * mSqq * w2m1)) / (3. * mSqq2);
        result.element({1, 4, 2, 2}) =
            ((-1 + w) * (rtSq * (-1 + 3 * rCSq - 2 * rC * (-2 + w) - 4 * w) + 2 * mSqq * (-1 + rC) * (1 + w))) /
            (3. * mSqq2);
        result.element({1, 4, 2, 3}) =
            ((-1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) / (3. * mSqq2);
        result.element({1, 4, 2, 4}) =
            ((-1 + w) * (rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq) + 2 * mSqq * w2m1)) / (3. * mSqq2);
        result.element({1, 4, 3, 0}) = -((rt * (rC - w) * (-1 + w)) / mSqq);
        result.element({1, 4, 4, 0}) = -((rt * (rC - w) * (-1 + w)) / mSqq);
        result.element({1, 4, 5, 8}) = (4 * rt * w2m1) / mSqq;
        result.element({1, 4, 5, 9}) = (4 * rt * w2m1) / mSqq;
        result.element({1, 4, 5, 10}) = (4 * rt * w2m1) / mSqq;
        result.element({1, 4, 5, 11}) = (-4 * rt * wm1Sq * (1 + w)) / mSqq;
        result.element({1, 5, 0, 5}) = (2 * (2 * mSqq2 * (-1 + w) + mSqq * (1 + 2 * rC + rCSq + rtSq) * (-1 + w) +
                                             rtSq * (1 + 2 * w - 2 * rC * (2 + w) + rCSq * (1 + 2 * w)))) /
                                       (3. * mSqq2);
        result.element({1, 5, 0, 6}) =
            ((1 + w) * (2 * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2 * rC * (2 + w) + rCSq * (-1 + 4 * w)))) /
            (3. * mSqq2);
        result.element({1, 5, 0, 7}) =
            ((1 + w) * (2 * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3 * rCSq + 4 * w - 2 * rC * (2 + w)))) /
            (3. * mSqq2);
        result.element({1, 5, 1, 5}) = (2 * (2 * mSqq2 * (-1 + w) + mSqq * (1 + 2 * rC + rCSq + rtSq) * (-1 + w) +
                                             rtSq * (1 + 2 * w - 2 * rC * (2 + w) + rCSq * (1 + 2 * w)))) /
                                       (3. * mSqq2);
        result.element({1, 5, 1, 6}) =
            ((1 + w) * (2 * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2 * rC * (2 + w) + rCSq * (-1 + 4 * w)))) /
            (3. * mSqq2);
        result.element({1, 5, 1, 7}) =
            ((1 + w) * (2 * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3 * rCSq + 4 * w - 2 * rC * (2 + w)))) /
            (3. * mSqq2);
        result.element({1, 5, 2, 5}) = (-2 * (2 * mSqq2 * (-1 + w) + mSqq * (1 + 2 * rC + rCSq + rtSq) * (-1 + w) +
                                              rtSq * (1 + 2 * w - 2 * rC * (2 + w) + rCSq * (1 + 2 * w)))) /
                                       (3. * mSqq2);
        result.element({1, 5, 2, 6}) =
            -((1 + w) * (2 * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2 * rC * (2 + w) + rCSq * (-1 + 4 * w)))) /
            (3. * mSqq2);
        result.element({1, 5, 2, 7}) =
            -((1 + w) * (2 * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3 * rCSq + 4 * w - 2 * rC * (2 + w)))) /
            (3. * mSqq2);
        result.element({1, 5, 3, 1}) = ((-1 + rC) * rt * (1 + w)) / mSqq;
        result.element({1, 5, 4, 1}) = -(((-1 + rC) * rt * (1 + w)) / mSqq);
        result.element({1, 5, 5, 8}) = (12 * (1 + rC) * rt * (-1 + w)) / mSqq;
        result.element({1, 5, 5, 9}) = (8 * rC * rt * w2m1) / mSqq;
        result.element({1, 5, 5, 10}) = (8 * rt * w2m1) / mSqq;
        result.element({1, 6, 0, 5}) =
            ((1 + w) * (2 * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2 * rC * (2 + w) + rCSq * (-1 + 4 * w)))) /
            (3. * mSqq2);
        result.element({1, 6, 0, 6}) =
            ((1 + w) * (rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)) + 2 * mSqq * rCSq * w2m1)) / (3. * mSqq2);
        result.element({1, 6, 0, 7}) =
            ((1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) / (3. * mSqq2);
        result.element({1, 6, 1, 5}) =
            ((1 + w) * (2 * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2 * rC * (2 + w) + rCSq * (-1 + 4 * w)))) /
            (3. * mSqq2);
        result.element({1, 6, 1, 6}) =
            ((1 + w) * (rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)) + 2 * mSqq * rCSq * w2m1)) / (3. * mSqq2);
        result.element({1, 6, 1, 7}) =
            ((1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) / (3. * mSqq2);
        result.element({1, 6, 2, 5}) =
            -((1 + w) * (2 * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2 * rC * (2 + w) + rCSq * (-1 + 4 * w)))) /
            (3. * mSqq2);
        result.element({1, 6, 2, 6}) =
            -((1 + w) * (rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)) + 2 * mSqq * rCSq * w2m1)) / (3. * mSqq2);
        result.element({1, 6, 2, 7}) =
            -((1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) / (3. * mSqq2);
        result.element({1, 6, 3, 1}) = (rt * (1 + w) * (-1 + rC * w)) / mSqq;
        result.element({1, 6, 4, 1}) = -((rt * (1 + w) * (-1 + rC * w)) / mSqq);
        result.element({1, 6, 5, 8}) = (4 * rC * rt * w2m1) / mSqq;
        result.element({1, 7, 0, 5}) =
            ((1 + w) * (2 * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3 * rCSq + 4 * w - 2 * rC * (2 + w)))) /
            (3. * mSqq2);
        result.element({1, 7, 0, 6}) =
            ((1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) / (3. * mSqq2);
        result.element({1, 7, 0, 7}) =
            ((1 + w) * (rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq) + 2 * mSqq * w2m1)) / (3. * mSqq2);
        result.element({1, 7, 1, 5}) =
            ((1 + w) * (2 * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3 * rCSq + 4 * w - 2 * rC * (2 + w)))) /
            (3. * mSqq2);
        result.element({1, 7, 1, 6}) =
            ((1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) / (3. * mSqq2);
        result.element({1, 7, 1, 7}) =
            ((1 + w) * (rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq) + 2 * mSqq * w2m1)) / (3. * mSqq2);
        result.element({1, 7, 2, 5}) =
            -((1 + w) * (2 * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3 * rCSq + 4 * w - 2 * rC * (2 + w)))) /
            (3. * mSqq2);
        result.element({1, 7, 2, 6}) =
            -((1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) / (3. * mSqq2);
        result.element({1, 7, 2, 7}) =
            -((1 + w) * (rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq) + 2 * mSqq * w2m1)) / (3. * mSqq2);
        result.element({1, 7, 3, 1}) = (rt * (rC - w) * (1 + w)) / mSqq;
        result.element({1, 7, 4, 1}) = -((rt * (rC - w) * (1 + w)) / mSqq);
        result.element({1, 7, 5, 8}) = (4 * rt * w2m1) / mSqq;
        result.element({2, 2, 0, 2}) = (2 * (2 * mSqq2 * (1 + w) + mSqq * (1 - 2 * rC + rCSq + rtSq) * (1 + w) +
                                             rtSq * (-1 + 2 * rC * (-2 + w) + 2 * w + rCSq * (-1 + 2 * w)))) /
                                       (3. * mSqq2);
        result.element({2, 2, 0, 3}) =
            ((-1 + w) * (2 * mSqq * (-1 + rC) * rC * (1 + w) + rtSq * (-3 + 2 * rC * (-2 + w) + rCSq * (1 + 4 * w)))) /
            (3. * mSqq2);
        result.element({2, 2, 0, 4}) =
            ((-1 + w) * (rtSq * (-1 + 3 * rCSq - 2 * rC * (-2 + w) - 4 * w) + 2 * mSqq * (-1 + rC) * (1 + w))) /
            (3. * mSqq2);
        result.element({2, 2, 1, 2}) = (2 * (2 * mSqq2 * (1 + w) + mSqq * (1 - 2 * rC + rCSq + rtSq) * (1 + w) +
                                             rtSq * (-1 + 2 * rC * (-2 + w) + 2 * w + rCSq * (-1 + 2 * w)))) /
                                       (3. * mSqq2);
        result.element({2, 2, 1, 3}) =
            ((-1 + w) * (2 * mSqq * (-1 + rC) * rC * (1 + w) + rtSq * (-3 + 2 * rC * (-2 + w) + rCSq * (1 + 4 * w)))) /
            (3. * mSqq2);
        result.element({2, 2, 1, 4}) =
            ((-1 + w) * (rtSq * (-1 + 3 * rCSq - 2 * rC * (-2 + w) - 4 * w) + 2 * mSqq * (-1 + rC) * (1 + w))) /
            (3. * mSqq2);
        result.element({2, 2, 2, 2}) = (2 * (2 * mSqq2 * (1 + w) + mSqq * (1 - 2 * rC + rCSq + rtSq) * (1 + w) +
                                             rtSq * (-1 + 2 * rC * (-2 + w) + 2 * w + rCSq * (-1 + 2 * w)))) /
                                       (3. * mSqq2);
        result.element({2, 2, 2, 3}) =
            ((-1 + w) * (2 * mSqq * (-1 + rC) * rC * (1 + w) + rtSq * (-3 + 2 * rC * (-2 + w) + rCSq * (1 + 4 * w)))) /
            (3. * mSqq2);
        result.element({2, 2, 2, 4}) =
            ((-1 + w) * (rtSq * (-1 + 3 * rCSq - 2 * rC * (-2 + w) - 4 * w) + 2 * mSqq * (-1 + rC) * (1 + w))) /
            (3. * mSqq2);
        result.element({2, 2, 3, 0}) = -(((1 + rC) * rt * (-1 + w)) / mSqq);
        result.element({2, 2, 4, 0}) = -(((1 + rC) * rt * (-1 + w)) / mSqq);
        result.element({2, 2, 5, 8}) = (12 * (-1 + rC) * rt * (1 + w)) / mSqq;
        result.element({2, 2, 5, 9}) = (4 * rt * (1 + w) * (-3 + rC + 2 * rC * w)) / mSqq;
        result.element({2, 2, 5, 10}) = (4 * rt * (-1 + 3 * rC - 2 * w) * (1 + w)) / mSqq;
        result.element({2, 2, 5, 11}) = (-4 * (-1 + rC) * rt * w2m1) / mSqq;
        result.element({2, 3, 0, 2}) =
            ((-1 + w) * (2 * mSqq * (-1 + rC) * rC * (1 + w) + rtSq * (-3 + 2 * rC * (-2 + w) + rCSq * (1 + 4 * w)))) /
            (3. * mSqq2);
        result.element({2, 3, 0, 3}) =
            ((-1 + w) * (rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)) + 2 * mSqq * rCSq * w2m1)) / (3. * mSqq2);
        result.element({2, 3, 0, 4}) =
            ((-1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) / (3. * mSqq2);
        result.element({2, 3, 1, 2}) =
            ((-1 + w) * (2 * mSqq * (-1 + rC) * rC * (1 + w) + rtSq * (-3 + 2 * rC * (-2 + w) + rCSq * (1 + 4 * w)))) /
            (3. * mSqq2);
        result.element({2, 3, 1, 3}) =
            ((-1 + w) * (rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)) + 2 * mSqq * rCSq * w2m1)) / (3. * mSqq2);
        result.element({2, 3, 1, 4}) =
            ((-1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) / (3. * mSqq2);
        result.element({2, 3, 2, 2}) =
            ((-1 + w) * (2 * mSqq * (-1 + rC) * rC * (1 + w) + rtSq * (-3 + 2 * rC * (-2 + w) + rCSq * (1 + 4 * w)))) /
            (3. * mSqq2);
        result.element({2, 3, 2, 3}) =
            ((-1 + w) * (rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)) + 2 * mSqq * rCSq * w2m1)) / (3. * mSqq2);
        result.element({2, 3, 2, 4}) =
            ((-1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) / (3. * mSqq2);
        result.element({2, 3, 3, 0}) = -((rt * (-1 + w) * (-1 + rC * w)) / mSqq);
        result.element({2, 3, 4, 0}) = -((rt * (-1 + w) * (-1 + rC * w)) / mSqq);
        result.element({2, 3, 5, 8}) = (4 * rC * rt * w2m1) / mSqq;
        result.element({2, 3, 5, 9}) = (4 * rC * rt * w2m1) / mSqq;
        result.element({2, 3, 5, 10}) = (4 * rC * rt * w2m1) / mSqq;
        result.element({2, 3, 5, 11}) = (-4 * rC * rt * wm1Sq * (1 + w)) / mSqq;
        result.element({2, 4, 0, 2}) =
            ((-1 + w) * (rtSq * (-1 + 3 * rCSq - 2 * rC * (-2 + w) - 4 * w) + 2 * mSqq * (-1 + rC) * (1 + w))) /
            (3. * mSqq2);
        result.element({2, 4, 0, 3}) =
            ((-1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) / (3. * mSqq2);
        result.element({2, 4, 0, 4}) =
            ((-1 + w) * (rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq) + 2 * mSqq * w2m1)) / (3. * mSqq2);
        result.element({2, 4, 1, 2}) =
            ((-1 + w) * (rtSq * (-1 + 3 * rCSq - 2 * rC * (-2 + w) - 4 * w) + 2 * mSqq * (-1 + rC) * (1 + w))) /
            (3. * mSqq2);
        result.element({2, 4, 1, 3}) =
            ((-1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) / (3. * mSqq2);
        result.element({2, 4, 1, 4}) =
            ((-1 + w) * (rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq) + 2 * mSqq * w2m1)) / (3. * mSqq2);
        result.element({2, 4, 2, 2}) =
            ((-1 + w) * (rtSq * (-1 + 3 * rCSq - 2 * rC * (-2 + w) - 4 * w) + 2 * mSqq * (-1 + rC) * (1 + w))) /
            (3. * mSqq2);
        result.element({2, 4, 2, 3}) =
            ((-1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) / (3. * mSqq2);
        result.element({2, 4, 2, 4}) =
            ((-1 + w) * (rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq) + 2 * mSqq * w2m1)) / (3. * mSqq2);
        result.element({2, 4, 3, 0}) = -((rt * (rC - w) * (-1 + w)) / mSqq);
        result.element({2, 4, 4, 0}) = -((rt * (rC - w) * (-1 + w)) / mSqq);
        result.element({2, 4, 5, 8}) = (4 * rt * w2m1) / mSqq;
        result.element({2, 4, 5, 9}) = (4 * rt * w2m1) / mSqq;
        result.element({2, 4, 5, 10}) = (4 * rt * w2m1) / mSqq;
        result.element({2, 4, 5, 11}) = (-4 * rt * wm1Sq * (1 + w)) / mSqq;
        result.element({2, 5, 0, 5}) = (-2 * (2 * mSqq2 * (-1 + w) + mSqq * (1 + 2 * rC + rCSq + rtSq) * (-1 + w) +
                                              rtSq * (1 + 2 * w - 2 * rC * (2 + w) + rCSq * (1 + 2 * w)))) /
                                       (3. * mSqq2);
        result.element({2, 5, 0, 6}) =
            -((1 + w) * (2 * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2 * rC * (2 + w) + rCSq * (-1 + 4 * w)))) /
            (3. * mSqq2);
        result.element({2, 5, 0, 7}) =
            -((1 + w) * (2 * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3 * rCSq + 4 * w - 2 * rC * (2 + w)))) /
            (3. * mSqq2);
        result.element({2, 5, 1, 5}) = (-2 * (2 * mSqq2 * (-1 + w) + mSqq * (1 + 2 * rC + rCSq + rtSq) * (-1 + w) +
                                              rtSq * (1 + 2 * w - 2 * rC * (2 + w) + rCSq * (1 + 2 * w)))) /
                                       (3. * mSqq2);
        result.element({2, 5, 1, 6}) =
            -((1 + w) * (2 * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2 * rC * (2 + w) + rCSq * (-1 + 4 * w)))) /
            (3. * mSqq2);
        result.element({2, 5, 1, 7}) =
            -((1 + w) * (2 * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3 * rCSq + 4 * w - 2 * rC * (2 + w)))) /
            (3. * mSqq2);
        result.element({2, 5, 2, 5}) = (2 * (2 * mSqq2 * (-1 + w) + mSqq * (1 + 2 * rC + rCSq + rtSq) * (-1 + w) +
                                             rtSq * (1 + 2 * w - 2 * rC * (2 + w) + rCSq * (1 + 2 * w)))) /
                                       (3. * mSqq2);
        result.element({2, 5, 2, 6}) =
            ((1 + w) * (2 * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2 * rC * (2 + w) + rCSq * (-1 + 4 * w)))) /
            (3. * mSqq2);
        result.element({2, 5, 2, 7}) =
            ((1 + w) * (2 * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3 * rCSq + 4 * w - 2 * rC * (2 + w)))) /
            (3. * mSqq2);
        result.element({2, 5, 3, 1}) = -(((-1 + rC) * rt * (1 + w)) / mSqq);
        result.element({2, 5, 4, 1}) = ((-1 + rC) * rt * (1 + w)) / mSqq;
        result.element({2, 5, 5, 8}) = (-12 * (1 + rC) * rt * (-1 + w)) / mSqq;
        result.element({2, 5, 5, 9}) = (-8 * rC * rt * w2m1) / mSqq;
        result.element({2, 5, 5, 10}) = (-8 * rt * w2m1) / mSqq;
        result.element({2, 6, 0, 5}) =
            -((1 + w) * (2 * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2 * rC * (2 + w) + rCSq * (-1 + 4 * w)))) /
            (3. * mSqq2);
        result.element({2, 6, 0, 6}) =
            -((1 + w) * (rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)) + 2 * mSqq * rCSq * w2m1)) / (3. * mSqq2);
        result.element({2, 6, 0, 7}) =
            -((1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) / (3. * mSqq2);
        result.element({2, 6, 1, 5}) =
            -((1 + w) * (2 * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2 * rC * (2 + w) + rCSq * (-1 + 4 * w)))) /
            (3. * mSqq2);
        result.element({2, 6, 1, 6}) =
            -((1 + w) * (rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)) + 2 * mSqq * rCSq * w2m1)) / (3. * mSqq2);
        result.element({2, 6, 1, 7}) =
            -((1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) / (3. * mSqq2);
        result.element({2, 6, 2, 5}) =
            ((1 + w) * (2 * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2 * rC * (2 + w) + rCSq * (-1 + 4 * w)))) /
            (3. * mSqq2);
        result.element({2, 6, 2, 6}) =
            ((1 + w) * (rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)) + 2 * mSqq * rCSq * w2m1)) / (3. * mSqq2);
        result.element({2, 6, 2, 7}) =
            ((1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) / (3. * mSqq2);
        result.element({2, 6, 3, 1}) = -((rt * (1 + w) * (-1 + rC * w)) / mSqq);
        result.element({2, 6, 4, 1}) = (rt * (1 + w) * (-1 + rC * w)) / mSqq;
        result.element({2, 6, 5, 8}) = (-4 * rC * rt * w2m1) / mSqq;
        result.element({2, 7, 0, 5}) =
            -((1 + w) * (2 * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3 * rCSq + 4 * w - 2 * rC * (2 + w)))) /
            (3. * mSqq2);
        result.element({2, 7, 0, 6}) =
            -((1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) / (3. * mSqq2);
        result.element({2, 7, 0, 7}) =
            -((1 + w) * (rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq) + 2 * mSqq * w2m1)) / (3. * mSqq2);
        result.element({2, 7, 1, 5}) =
            -((1 + w) * (2 * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3 * rCSq + 4 * w - 2 * rC * (2 + w)))) /
            (3. * mSqq2);
        result.element({2, 7, 1, 6}) =
            -((1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) / (3. * mSqq2);
        result.element({2, 7, 1, 7}) =
            -((1 + w) * (rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq) + 2 * mSqq * w2m1)) / (3. * mSqq2);
        result.element({2, 7, 2, 5}) =
            ((1 + w) * (2 * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3 * rCSq + 4 * w - 2 * rC * (2 + w)))) /
            (3. * mSqq2);
        result.element({2, 7, 2, 6}) =
            ((1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) / (3. * mSqq2);
        result.element({2, 7, 2, 7}) =
            ((1 + w) * (rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq) + 2 * mSqq * w2m1)) / (3. * mSqq2);
        result.element({2, 7, 3, 1}) = -((rt * (rC - w) * (1 + w)) / mSqq);
        result.element({2, 7, 4, 1}) = (rt * (rC - w) * (1 + w)) / mSqq;
        result.element({2, 7, 5, 8}) = (-4 * rt * w2m1) / mSqq;
        result.element({3, 0, 0, 2}) = -(((1 + rC) * rt * (-1 + w)) / mSqq);
        result.element({3, 0, 0, 3}) = -((rt * (-1 + w) * (-1 + rC * w)) / mSqq);
        result.element({3, 0, 0, 4}) = -((rt * (rC - w) * (-1 + w)) / mSqq);
        result.element({3, 0, 1, 2}) = -(((1 + rC) * rt * (-1 + w)) / mSqq);
        result.element({3, 0, 1, 3}) = -((rt * (-1 + w) * (-1 + rC * w)) / mSqq);
        result.element({3, 0, 1, 4}) = -((rt * (rC - w) * (-1 + w)) / mSqq);
        result.element({3, 0, 2, 2}) = -(((1 + rC) * rt * (-1 + w)) / mSqq);
        result.element({3, 0, 2, 3}) = -((rt * (-1 + w) * (-1 + rC * w)) / mSqq);
        result.element({3, 0, 2, 4}) = -((rt * (rC - w) * (-1 + w)) / mSqq);
        result.element({3, 0, 3, 0}) = -1 + w;
        result.element({3, 0, 4, 0}) = -1 + w;
        result.element({3, 1, 0, 5}) = ((-1 + rC) * rt * (1 + w)) / mSqq;
        result.element({3, 1, 0, 6}) = (rt * (1 + w) * (-1 + rC * w)) / mSqq;
        result.element({3, 1, 0, 7}) = (rt * (rC - w) * (1 + w)) / mSqq;
        result.element({3, 1, 1, 5}) = ((-1 + rC) * rt * (1 + w)) / mSqq;
        result.element({3, 1, 1, 6}) = (rt * (1 + w) * (-1 + rC * w)) / mSqq;
        result.element({3, 1, 1, 7}) = (rt * (rC - w) * (1 + w)) / mSqq;
        result.element({3, 1, 2, 5}) = -(((-1 + rC) * rt * (1 + w)) / mSqq);
        result.element({3, 1, 2, 6}) = -((rt * (1 + w) * (-1 + rC * w)) / mSqq);
        result.element({3, 1, 2, 7}) = -((rt * (rC - w) * (1 + w)) / mSqq);
        result.element({3, 1, 3, 1}) = 1 + w;
        result.element({3, 1, 4, 1}) = -1 - w;
        result.element({4, 0, 0, 2}) = -(((1 + rC) * rt * (-1 + w)) / mSqq);
        result.element({4, 0, 0, 3}) = -((rt * (-1 + w) * (-1 + rC * w)) / mSqq);
        result.element({4, 0, 0, 4}) = -((rt * (rC - w) * (-1 + w)) / mSqq);
        result.element({4, 0, 1, 2}) = -(((1 + rC) * rt * (-1 + w)) / mSqq);
        result.element({4, 0, 1, 3}) = -((rt * (-1 + w) * (-1 + rC * w)) / mSqq);
        result.element({4, 0, 1, 4}) = -((rt * (rC - w) * (-1 + w)) / mSqq);
        result.element({4, 0, 2, 2}) = -(((1 + rC) * rt * (-1 + w)) / mSqq);
        result.element({4, 0, 2, 3}) = -((rt * (-1 + w) * (-1 + rC * w)) / mSqq);
        result.element({4, 0, 2, 4}) = -((rt * (rC - w) * (-1 + w)) / mSqq);
        result.element({4, 0, 3, 0}) = -1 + w;
        result.element({4, 0, 4, 0}) = -1 + w;
        result.element({4, 1, 0, 5}) = -(((-1 + rC) * rt * (1 + w)) / mSqq);
        result.element({4, 1, 0, 6}) = -((rt * (1 + w) * (-1 + rC * w)) / mSqq);
        result.element({4, 1, 0, 7}) = -((rt * (rC - w) * (1 + w)) / mSqq);
        result.element({4, 1, 1, 5}) = -(((-1 + rC) * rt * (1 + w)) / mSqq);
        result.element({4, 1, 1, 6}) = -((rt * (1 + w) * (-1 + rC * w)) / mSqq);
        result.element({4, 1, 1, 7}) = -((rt * (rC - w) * (1 + w)) / mSqq);
        result.element({4, 1, 2, 5}) = ((-1 + rC) * rt * (1 + w)) / mSqq;
        result.element({4, 1, 2, 6}) = (rt * (1 + w) * (-1 + rC * w)) / mSqq;
        result.element({4, 1, 2, 7}) = (rt * (rC - w) * (1 + w)) / mSqq;
        result.element({4, 1, 3, 1}) = -1 - w;
        result.element({4, 1, 4, 1}) = 1 + w;
        result.element({5, 8, 0, 2}) = (12 * (-1 + rC) * rt * (1 + w)) / mSqq;
        result.element({5, 8, 0, 3}) = (4 * rC * rt * w2m1) / mSqq;
        result.element({5, 8, 0, 4}) = (4 * rt * w2m1) / mSqq;
        result.element({5, 8, 0, 5}) = (12 * (1 + rC) * rt * (-1 + w)) / mSqq;
        result.element({5, 8, 0, 6}) = (4 * rC * rt * w2m1) / mSqq;
        result.element({5, 8, 0, 7}) = (4 * rt * w2m1) / mSqq;
        result.element({5, 8, 1, 2}) = (12 * (-1 + rC) * rt * (1 + w)) / mSqq;
        result.element({5, 8, 1, 3}) = (4 * rC * rt * w2m1) / mSqq;
        result.element({5, 8, 1, 4}) = (4 * rt * w2m1) / mSqq;
        result.element({5, 8, 1, 5}) = (12 * (1 + rC) * rt * (-1 + w)) / mSqq;
        result.element({5, 8, 1, 6}) = (4 * rC * rt * w2m1) / mSqq;
        result.element({5, 8, 1, 7}) = (4 * rt * w2m1) / mSqq;
        result.element({5, 8, 2, 2}) = (12 * (-1 + rC) * rt * (1 + w)) / mSqq;
        result.element({5, 8, 2, 3}) = (4 * rC * rt * w2m1) / mSqq;
        result.element({5, 8, 2, 4}) = (4 * rt * w2m1) / mSqq;
        result.element({5, 8, 2, 5}) = (-12 * (1 + rC) * rt * (-1 + w)) / mSqq;
        result.element({5, 8, 2, 6}) = (-4 * rC * rt * w2m1) / mSqq;
        result.element({5, 8, 2, 7}) = (-4 * rt * w2m1) / mSqq;
        result.element({5, 8, 5, 8}) =
            (32 * (mSqq + 2 * rtSq) * (-4 * rC + (2 + mSqq) * w + 2 * rCSq * w)) / (3. * mSqq2);
        result.element({5, 8, 5, 9}) =
            (16 * (mSqq + 2 * rtSq) * (1 + w) * (2 + mSqq - 4 * rC + rCSq * (-2 + 4 * w))) / (3. * mSqq2);
        result.element({5, 8, 5, 10}) =
            (16 * (mSqq + 2 * rtSq) * (1 + w) * (-2 + mSqq - 4 * rC + 2 * rCSq + 4 * w)) / (3. * mSqq2);
        result.element({5, 8, 5, 11}) = (-16 * (mSqq + 2 * rtSq) * w2m1) / (3. * mSqq);
        result.element({5, 9, 0, 2}) = (4 * rt * (1 + w) * (-3 + rC + 2 * rC * w)) / mSqq;
        result.element({5, 9, 0, 3}) = (4 * rC * rt * w2m1) / mSqq;
        result.element({5, 9, 0, 4}) = (4 * rt * w2m1) / mSqq;
        result.element({5, 9, 0, 5}) = (8 * rC * rt * w2m1) / mSqq;
        result.element({5, 9, 1, 2}) = (4 * rt * (1 + w) * (-3 + rC + 2 * rC * w)) / mSqq;
        result.element({5, 9, 1, 3}) = (4 * rC * rt * w2m1) / mSqq;
        result.element({5, 9, 1, 4}) = (4 * rt * w2m1) / mSqq;
        result.element({5, 9, 1, 5}) = (8 * rC * rt * w2m1) / mSqq;
        result.element({5, 9, 2, 2}) = (4 * rt * (1 + w) * (-3 + rC + 2 * rC * w)) / mSqq;
        result.element({5, 9, 2, 3}) = (4 * rC * rt * w2m1) / mSqq;
        result.element({5, 9, 2, 4}) = (4 * rt * w2m1) / mSqq;
        result.element({5, 9, 2, 5}) = (-8 * rC * rt * w2m1) / mSqq;
        result.element({5, 9, 5, 8}) =
            (16 * (mSqq + 2 * rtSq) * (1 + w) * (2 + mSqq - 4 * rC + rCSq * (-2 + 4 * w))) / (3. * mSqq2);
        result.element({5, 9, 5, 9}) = (16 * (1 + w) *
                                        (mSqq2 + 4 * rtSq * (1 - 2 * rC * w + rCSq * (-1 + 2 * wSq)) +
                                         2 * mSqq * (1 + rtSq - 2 * rC * w + rCSq * (-1 + 2 * wSq)))) /
                                       (3. * mSqq2);
        result.element({5, 9, 5, 10}) =
            (16 * (mSqq + 2 * rtSq) * (1 + w) * (mSqq + 2 * (-2 * rC + w + rCSq * w))) / (3. * mSqq2);
        result.element({5, 9, 5, 11}) = (-16 * (mSqq + 2 * rtSq) * w2m1) / (3. * mSqq);
        result.element({5, 10, 0, 2}) = (4 * rt * (-1 + 3 * rC - 2 * w) * (1 + w)) / mSqq;
        result.element({5, 10, 0, 3}) = (4 * rC * rt * w2m1) / mSqq;
        result.element({5, 10, 0, 4}) = (4 * rt * w2m1) / mSqq;
        result.element({5, 10, 0, 5}) = (8 * rt * w2m1) / mSqq;
        result.element({5, 10, 1, 2}) = (4 * rt * (-1 + 3 * rC - 2 * w) * (1 + w)) / mSqq;
        result.element({5, 10, 1, 3}) = (4 * rC * rt * w2m1) / mSqq;
        result.element({5, 10, 1, 4}) = (4 * rt * w2m1) / mSqq;
        result.element({5, 10, 1, 5}) = (8 * rt * w2m1) / mSqq;
        result.element({5, 10, 2, 2}) = (4 * rt * (-1 + 3 * rC - 2 * w) * (1 + w)) / mSqq;
        result.element({5, 10, 2, 3}) = (4 * rC * rt * w2m1) / mSqq;
        result.element({5, 10, 2, 4}) = (4 * rt * w2m1) / mSqq;
        result.element({5, 10, 2, 5}) = (-8 * rt * w2m1) / mSqq;
        result.element({5, 10, 5, 8}) =
            (16 * (mSqq + 2 * rtSq) * (1 + w) * (-2 + mSqq - 4 * rC + 2 * rCSq + 4 * w)) / (3. * mSqq2);
        result.element({5, 10, 5, 9}) =
            (16 * (mSqq + 2 * rtSq) * (1 + w) * (mSqq + 2 * (-2 * rC + w + rCSq * w))) / (3. * mSqq2);
        result.element({5, 10, 5, 10}) = (16 * (1 + w) *
                                          (mSqq2 + 4 * rtSq * (-1 + rCSq - 2 * rC * w + 2 * wSq) +
                                           2 * mSqq * (-1 + rCSq + rtSq - 2 * rC * w + 2 * wSq))) /
                                         (3. * mSqq2);
        result.element({5, 10, 5, 11}) = (-16 * (mSqq + 2 * rtSq) * w2m1) / (3. * mSqq);
        result.element({5, 11, 0, 2}) = (-4 * (-1 + rC) * rt * w2m1) / mSqq;
        result.element({5, 11, 0, 3}) = (-4 * rC * rt * wm1Sq * (1 + w)) / mSqq;
        result.element({5, 11, 0, 4}) = (-4 * rt * wm1Sq * (1 + w)) / mSqq;
        result.element({5, 11, 1, 2}) = (-4 * (-1 + rC) * rt * w2m1) / mSqq;
        result.element({5, 11, 1, 3}) = (-4 * rC * rt * wm1Sq * (1 + w)) / mSqq;
        result.element({5, 11, 1, 4}) = (-4 * rt * wm1Sq * (1 + w)) / mSqq;
        result.element({5, 11, 2, 2}) = (-4 * (-1 + rC) * rt * w2m1) / mSqq;
        result.element({5, 11, 2, 3}) = (-4 * rC * rt * wm1Sq * (1 + w)) / mSqq;
        result.element({5, 11, 2, 4}) = (-4 * rt * wm1Sq * (1 + w)) / mSqq;
        result.element({5, 11, 5, 8}) = (-16 * (mSqq + 2 * rtSq) * w2m1) / (3. * mSqq);
        result.element({5, 11, 5, 9}) = (-16 * (mSqq + 2 * rtSq) * w2m1) / (3. * mSqq);
        result.element({5, 11, 5, 10}) = (-16 * (mSqq + 2 * rtSq) * w2m1) / (3. * mSqq);
        result.element({5, 11, 5, 11}) = (16 * (mSqq + 2 * rtSq) * wm1Sq * (1 + w)) / (3. * mSqq);
        result.element({6, 2, 6, 2}) = (2 * (2 * mSqq2 * (1 + w) + mSqq * (1 - 2 * rC + rCSq + rtSq) * (1 + w) +
                                             rtSq * (-1 + 2 * rC * (-2 + w) + 2 * w + rCSq * (-1 + 2 * w)))) /
                                       (3. * mSqq2);
        result.element({6, 2, 6, 3}) =
            ((-1 + w) * (2 * mSqq * (-1 + rC) * rC * (1 + w) + rtSq * (-3 + 2 * rC * (-2 + w) + rCSq * (1 + 4 * w)))) /
            (3. * mSqq2);
        result.element({6, 2, 6, 4}) =
            ((-1 + w) * (rtSq * (-1 + 3 * rCSq - 2 * rC * (-2 + w) - 4 * w) + 2 * mSqq * (-1 + rC) * (1 + w))) /
            (3. * mSqq2);
        result.element({6, 2, 7, 2}) = (2 * (2 * mSqq2 * (1 + w) + mSqq * (1 - 2 * rC + rCSq + rtSq) * (1 + w) +
                                             rtSq * (-1 + 2 * rC * (-2 + w) + 2 * w + rCSq * (-1 + 2 * w)))) /
                                       (3. * mSqq2);
        result.element({6, 2, 7, 3}) =
            ((-1 + w) * (2 * mSqq * (-1 + rC) * rC * (1 + w) + rtSq * (-3 + 2 * rC * (-2 + w) + rCSq * (1 + 4 * w)))) /
            (3. * mSqq2);
        result.element({6, 2, 7, 4}) =
            ((-1 + w) * (rtSq * (-1 + 3 * rCSq - 2 * rC * (-2 + w) - 4 * w) + 2 * mSqq * (-1 + rC) * (1 + w))) /
            (3. * mSqq2);
        result.element({6, 2, 8, 0}) = -(((1 + rC) * rt * (-1 + w)) / mSqq);
        result.element({6, 2, 9, 0}) = -(((1 + rC) * rt * (-1 + w)) / mSqq);
        result.element({6, 2, 10, 8}) = (12 * (-1 + rC) * rt * (1 + w)) / mSqq;
        result.element({6, 2, 10, 9}) = (4 * rt * (1 + w) * (-3 + rC + 2 * rC * w)) / mSqq;
        result.element({6, 2, 10, 10}) = (4 * rt * (-1 + 3 * rC - 2 * w) * (1 + w)) / mSqq;
        result.element({6, 2, 10, 11}) = (-4 * (-1 + rC) * rt * w2m1) / mSqq;
        result.element({6, 3, 6, 2}) =
            ((-1 + w) * (2 * mSqq * (-1 + rC) * rC * (1 + w) + rtSq * (-3 + 2 * rC * (-2 + w) + rCSq * (1 + 4 * w)))) /
            (3. * mSqq2);
        result.element({6, 3, 6, 3}) =
            ((-1 + w) * (rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)) + 2 * mSqq * rCSq * w2m1)) / (3. * mSqq2);
        result.element({6, 3, 6, 4}) =
            ((-1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) / (3. * mSqq2);
        result.element({6, 3, 7, 2}) =
            ((-1 + w) * (2 * mSqq * (-1 + rC) * rC * (1 + w) + rtSq * (-3 + 2 * rC * (-2 + w) + rCSq * (1 + 4 * w)))) /
            (3. * mSqq2);
        result.element({6, 3, 7, 3}) =
            ((-1 + w) * (rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)) + 2 * mSqq * rCSq * w2m1)) / (3. * mSqq2);
        result.element({6, 3, 7, 4}) =
            ((-1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) / (3. * mSqq2);
        result.element({6, 3, 8, 0}) = -((rt * (-1 + w) * (-1 + rC * w)) / mSqq);
        result.element({6, 3, 9, 0}) = -((rt * (-1 + w) * (-1 + rC * w)) / mSqq);
        result.element({6, 3, 10, 8}) = (4 * rC * rt * w2m1) / mSqq;
        result.element({6, 3, 10, 9}) = (4 * rC * rt * w2m1) / mSqq;
        result.element({6, 3, 10, 10}) = (4 * rC * rt * w2m1) / mSqq;
        result.element({6, 3, 10, 11}) = (-4 * rC * rt * wm1Sq * (1 + w)) / mSqq;
        result.element({6, 4, 6, 2}) =
            ((-1 + w) * (rtSq * (-1 + 3 * rCSq - 2 * rC * (-2 + w) - 4 * w) + 2 * mSqq * (-1 + rC) * (1 + w))) /
            (3. * mSqq2);
        result.element({6, 4, 6, 3}) =
            ((-1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) / (3. * mSqq2);
        result.element({6, 4, 6, 4}) =
            ((-1 + w) * (rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq) + 2 * mSqq * w2m1)) / (3. * mSqq2);
        result.element({6, 4, 7, 2}) =
            ((-1 + w) * (rtSq * (-1 + 3 * rCSq - 2 * rC * (-2 + w) - 4 * w) + 2 * mSqq * (-1 + rC) * (1 + w))) /
            (3. * mSqq2);
        result.element({6, 4, 7, 3}) =
            ((-1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) / (3. * mSqq2);
        result.element({6, 4, 7, 4}) =
            ((-1 + w) * (rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq) + 2 * mSqq * w2m1)) / (3. * mSqq2);
        result.element({6, 4, 8, 0}) = -((rt * (rC - w) * (-1 + w)) / mSqq);
        result.element({6, 4, 9, 0}) = -((rt * (rC - w) * (-1 + w)) / mSqq);
        result.element({6, 4, 10, 8}) = (4 * rt * w2m1) / mSqq;
        result.element({6, 4, 10, 9}) = (4 * rt * w2m1) / mSqq;
        result.element({6, 4, 10, 10}) = (4 * rt * w2m1) / mSqq;
        result.element({6, 4, 10, 11}) = (-4 * rt * wm1Sq * (1 + w)) / mSqq;
        result.element({6, 5, 6, 5}) = (2 * (2 * mSqq2 * (-1 + w) + mSqq * (1 + 2 * rC + rCSq + rtSq) * (-1 + w) +
                                             rtSq * (1 + 2 * w - 2 * rC * (2 + w) + rCSq * (1 + 2 * w)))) /
                                       (3. * mSqq2);
        result.element({6, 5, 6, 6}) =
            ((1 + w) * (2 * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2 * rC * (2 + w) + rCSq * (-1 + 4 * w)))) /
            (3. * mSqq2);
        result.element({6, 5, 6, 7}) =
            ((1 + w) * (2 * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3 * rCSq + 4 * w - 2 * rC * (2 + w)))) /
            (3. * mSqq2);
        result.element({6, 5, 7, 5}) = (-2 * (2 * mSqq2 * (-1 + w) + mSqq * (1 + 2 * rC + rCSq + rtSq) * (-1 + w) +
                                              rtSq * (1 + 2 * w - 2 * rC * (2 + w) + rCSq * (1 + 2 * w)))) /
                                       (3. * mSqq2);
        result.element({6, 5, 7, 6}) =
            -((1 + w) * (2 * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2 * rC * (2 + w) + rCSq * (-1 + 4 * w)))) /
            (3. * mSqq2);
        result.element({6, 5, 7, 7}) =
            -((1 + w) * (2 * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3 * rCSq + 4 * w - 2 * rC * (2 + w)))) /
            (3. * mSqq2);
        result.element({6, 5, 8, 1}) = ((-1 + rC) * rt * (1 + w)) / mSqq;
        result.element({6, 5, 9, 1}) = -(((-1 + rC) * rt * (1 + w)) / mSqq);
        result.element({6, 5, 10, 8}) = (-12 * (1 + rC) * rt * (-1 + w)) / mSqq;
        result.element({6, 5, 10, 9}) = (-8 * rC * rt * w2m1) / mSqq;
        result.element({6, 5, 10, 10}) = (-8 * rt * w2m1) / mSqq;
        result.element({6, 6, 6, 5}) =
            ((1 + w) * (2 * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2 * rC * (2 + w) + rCSq * (-1 + 4 * w)))) /
            (3. * mSqq2);
        result.element({6, 6, 6, 6}) =
            ((1 + w) * (rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)) + 2 * mSqq * rCSq * w2m1)) / (3. * mSqq2);
        result.element({6, 6, 6, 7}) =
            ((1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) / (3. * mSqq2);
        result.element({6, 6, 7, 5}) =
            -((1 + w) * (2 * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2 * rC * (2 + w) + rCSq * (-1 + 4 * w)))) /
            (3. * mSqq2);
        result.element({6, 6, 7, 6}) =
            -((1 + w) * (rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)) + 2 * mSqq * rCSq * w2m1)) / (3. * mSqq2);
        result.element({6, 6, 7, 7}) =
            -((1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) / (3. * mSqq2);
        result.element({6, 6, 8, 1}) = (rt * (1 + w) * (-1 + rC * w)) / mSqq;
        result.element({6, 6, 9, 1}) = -((rt * (1 + w) * (-1 + rC * w)) / mSqq);
        result.element({6, 6, 10, 8}) = (-4 * rC * rt * w2m1) / mSqq;
        result.element({6, 7, 6, 5}) =
            ((1 + w) * (2 * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3 * rCSq + 4 * w - 2 * rC * (2 + w)))) /
            (3. * mSqq2);
        result.element({6, 7, 6, 6}) =
            ((1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) / (3. * mSqq2);
        result.element({6, 7, 6, 7}) =
            ((1 + w) * (rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq) + 2 * mSqq * w2m1)) / (3. * mSqq2);
        result.element({6, 7, 7, 5}) =
            -((1 + w) * (2 * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3 * rCSq + 4 * w - 2 * rC * (2 + w)))) /
            (3. * mSqq2);
        result.element({6, 7, 7, 6}) =
            -((1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) / (3. * mSqq2);
        result.element({6, 7, 7, 7}) =
            -((1 + w) * (rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq) + 2 * mSqq * w2m1)) / (3. * mSqq2);
        result.element({6, 7, 8, 1}) = (rt * (rC - w) * (1 + w)) / mSqq;
        result.element({6, 7, 9, 1}) = -((rt * (rC - w) * (1 + w)) / mSqq);
        result.element({6, 7, 10, 8}) = (-4 * rt * w2m1) / mSqq;
        result.element({7, 2, 6, 2}) = (2 * (2 * mSqq2 * (1 + w) + mSqq * (1 - 2 * rC + rCSq + rtSq) * (1 + w) +
                                             rtSq * (-1 + 2 * rC * (-2 + w) + 2 * w + rCSq * (-1 + 2 * w)))) /
                                       (3. * mSqq2);
        result.element({7, 2, 6, 3}) =
            ((-1 + w) * (2 * mSqq * (-1 + rC) * rC * (1 + w) + rtSq * (-3 + 2 * rC * (-2 + w) + rCSq * (1 + 4 * w)))) /
            (3. * mSqq2);
        result.element({7, 2, 6, 4}) =
            ((-1 + w) * (rtSq * (-1 + 3 * rCSq - 2 * rC * (-2 + w) - 4 * w) + 2 * mSqq * (-1 + rC) * (1 + w))) /
            (3. * mSqq2);
        result.element({7, 2, 7, 2}) = (2 * (2 * mSqq2 * (1 + w) + mSqq * (1 - 2 * rC + rCSq + rtSq) * (1 + w) +
                                             rtSq * (-1 + 2 * rC * (-2 + w) + 2 * w + rCSq * (-1 + 2 * w)))) /
                                       (3. * mSqq2);
        result.element({7, 2, 7, 3}) =
            ((-1 + w) * (2 * mSqq * (-1 + rC) * rC * (1 + w) + rtSq * (-3 + 2 * rC * (-2 + w) + rCSq * (1 + 4 * w)))) /
            (3. * mSqq2);
        result.element({7, 2, 7, 4}) =
            ((-1 + w) * (rtSq * (-1 + 3 * rCSq - 2 * rC * (-2 + w) - 4 * w) + 2 * mSqq * (-1 + rC) * (1 + w))) /
            (3. * mSqq2);
        result.element({7, 2, 8, 0}) = -(((1 + rC) * rt * (-1 + w)) / mSqq);
        result.element({7, 2, 9, 0}) = -(((1 + rC) * rt * (-1 + w)) / mSqq);
        result.element({7, 2, 10, 8}) = (12 * (-1 + rC) * rt * (1 + w)) / mSqq;
        result.element({7, 2, 10, 9}) = (4 * rt * (1 + w) * (-3 + rC + 2 * rC * w)) / mSqq;
        result.element({7, 2, 10, 10}) = (4 * rt * (-1 + 3 * rC - 2 * w) * (1 + w)) / mSqq;
        result.element({7, 2, 10, 11}) = (-4 * (-1 + rC) * rt * w2m1) / mSqq;
        result.element({7, 3, 6, 2}) =
            ((-1 + w) * (2 * mSqq * (-1 + rC) * rC * (1 + w) + rtSq * (-3 + 2 * rC * (-2 + w) + rCSq * (1 + 4 * w)))) /
            (3. * mSqq2);
        result.element({7, 3, 6, 3}) =
            ((-1 + w) * (rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)) + 2 * mSqq * rCSq * w2m1)) / (3. * mSqq2);
        result.element({7, 3, 6, 4}) =
            ((-1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) / (3. * mSqq2);
        result.element({7, 3, 7, 2}) =
            ((-1 + w) * (2 * mSqq * (-1 + rC) * rC * (1 + w) + rtSq * (-3 + 2 * rC * (-2 + w) + rCSq * (1 + 4 * w)))) /
            (3. * mSqq2);
        result.element({7, 3, 7, 3}) =
            ((-1 + w) * (rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)) + 2 * mSqq * rCSq * w2m1)) / (3. * mSqq2);
        result.element({7, 3, 7, 4}) =
            ((-1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) / (3. * mSqq2);
        result.element({7, 3, 8, 0}) = -((rt * (-1 + w) * (-1 + rC * w)) / mSqq);
        result.element({7, 3, 9, 0}) = -((rt * (-1 + w) * (-1 + rC * w)) / mSqq);
        result.element({7, 3, 10, 8}) = (4 * rC * rt * w2m1) / mSqq;
        result.element({7, 3, 10, 9}) = (4 * rC * rt * w2m1) / mSqq;
        result.element({7, 3, 10, 10}) = (4 * rC * rt * w2m1) / mSqq;
        result.element({7, 3, 10, 11}) = (-4 * rC * rt * wm1Sq * (1 + w)) / mSqq;
        result.element({7, 4, 6, 2}) =
            ((-1 + w) * (rtSq * (-1 + 3 * rCSq - 2 * rC * (-2 + w) - 4 * w) + 2 * mSqq * (-1 + rC) * (1 + w))) /
            (3. * mSqq2);
        result.element({7, 4, 6, 3}) =
            ((-1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) / (3. * mSqq2);
        result.element({7, 4, 6, 4}) =
            ((-1 + w) * (rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq) + 2 * mSqq * w2m1)) / (3. * mSqq2);
        result.element({7, 4, 7, 2}) =
            ((-1 + w) * (rtSq * (-1 + 3 * rCSq - 2 * rC * (-2 + w) - 4 * w) + 2 * mSqq * (-1 + rC) * (1 + w))) /
            (3. * mSqq2);
        result.element({7, 4, 7, 3}) =
            ((-1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) / (3. * mSqq2);
        result.element({7, 4, 7, 4}) =
            ((-1 + w) * (rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq) + 2 * mSqq * w2m1)) / (3. * mSqq2);
        result.element({7, 4, 8, 0}) = -((rt * (rC - w) * (-1 + w)) / mSqq);
        result.element({7, 4, 9, 0}) = -((rt * (rC - w) * (-1 + w)) / mSqq);
        result.element({7, 4, 10, 8}) = (4 * rt * w2m1) / mSqq;
        result.element({7, 4, 10, 9}) = (4 * rt * w2m1) / mSqq;
        result.element({7, 4, 10, 10}) = (4 * rt * w2m1) / mSqq;
        result.element({7, 4, 10, 11}) = (-4 * rt * wm1Sq * (1 + w)) / mSqq;
        result.element({7, 5, 6, 5}) = (-2 * (2 * mSqq2 * (-1 + w) + mSqq * (1 + 2 * rC + rCSq + rtSq) * (-1 + w) +
                                              rtSq * (1 + 2 * w - 2 * rC * (2 + w) + rCSq * (1 + 2 * w)))) /
                                       (3. * mSqq2);
        result.element({7, 5, 6, 6}) =
            -((1 + w) * (2 * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2 * rC * (2 + w) + rCSq * (-1 + 4 * w)))) /
            (3. * mSqq2);
        result.element({7, 5, 6, 7}) =
            -((1 + w) * (2 * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3 * rCSq + 4 * w - 2 * rC * (2 + w)))) /
            (3. * mSqq2);
        result.element({7, 5, 7, 5}) = (2 * (2 * mSqq2 * (-1 + w) + mSqq * (1 + 2 * rC + rCSq + rtSq) * (-1 + w) +
                                             rtSq * (1 + 2 * w - 2 * rC * (2 + w) + rCSq * (1 + 2 * w)))) /
                                       (3. * mSqq2);
        result.element({7, 5, 7, 6}) =
            ((1 + w) * (2 * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2 * rC * (2 + w) + rCSq * (-1 + 4 * w)))) /
            (3. * mSqq2);
        result.element({7, 5, 7, 7}) =
            ((1 + w) * (2 * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3 * rCSq + 4 * w - 2 * rC * (2 + w)))) /
            (3. * mSqq2);
        result.element({7, 5, 8, 1}) = -(((-1 + rC) * rt * (1 + w)) / mSqq);
        result.element({7, 5, 9, 1}) = ((-1 + rC) * rt * (1 + w)) / mSqq;
        result.element({7, 5, 10, 8}) = (12 * (1 + rC) * rt * (-1 + w)) / mSqq;
        result.element({7, 5, 10, 9}) = (8 * rC * rt * w2m1) / mSqq;
        result.element({7, 5, 10, 10}) = (8 * rt * w2m1) / mSqq;
        result.element({7, 6, 6, 5}) =
            -((1 + w) * (2 * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2 * rC * (2 + w) + rCSq * (-1 + 4 * w)))) /
            (3. * mSqq2);
        result.element({7, 6, 6, 6}) =
            -((1 + w) * (rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)) + 2 * mSqq * rCSq * w2m1)) / (3. * mSqq2);
        result.element({7, 6, 6, 7}) =
            -((1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) / (3. * mSqq2);
        result.element({7, 6, 7, 5}) =
            ((1 + w) * (2 * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2 * rC * (2 + w) + rCSq * (-1 + 4 * w)))) /
            (3. * mSqq2);
        result.element({7, 6, 7, 6}) =
            ((1 + w) * (rtSq * (3 - 6 * rC * w + rCSq * (-1 + 4 * wSq)) + 2 * mSqq * rCSq * w2m1)) / (3. * mSqq2);
        result.element({7, 6, 7, 7}) =
            ((1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) / (3. * mSqq2);
        result.element({7, 6, 8, 1}) = -((rt * (1 + w) * (-1 + rC * w)) / mSqq);
        result.element({7, 6, 9, 1}) = (rt * (1 + w) * (-1 + rC * w)) / mSqq;
        result.element({7, 6, 10, 8}) = (4 * rC * rt * w2m1) / mSqq;
        result.element({7, 7, 6, 5}) =
            -((1 + w) * (2 * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3 * rCSq + 4 * w - 2 * rC * (2 + w)))) /
            (3. * mSqq2);
        result.element({7, 7, 6, 6}) =
            -((1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) / (3. * mSqq2);
        result.element({7, 7, 6, 7}) =
            -((1 + w) * (rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq) + 2 * mSqq * w2m1)) / (3. * mSqq2);
        result.element({7, 7, 7, 5}) =
            ((1 + w) * (2 * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3 * rCSq + 4 * w - 2 * rC * (2 + w)))) /
            (3. * mSqq2);
        result.element({7, 7, 7, 6}) =
            ((1 + w) * (rtSq * (3 * w + 3 * rCSq * w - 2 * rC * (2 + wSq)) + 2 * mSqq * rC * w2m1)) / (3. * mSqq2);
        result.element({7, 7, 7, 7}) =
            ((1 + w) * (rtSq * (-1 + 3 * rCSq - 6 * rC * w + 4 * wSq) + 2 * mSqq * w2m1)) / (3. * mSqq2);
        result.element({7, 7, 8, 1}) = -((rt * (rC - w) * (1 + w)) / mSqq);
        result.element({7, 7, 9, 1}) = (rt * (rC - w) * (1 + w)) / mSqq;
        result.element({7, 7, 10, 8}) = (4 * rt * w2m1) / mSqq;
        result.element({8, 0, 6, 2}) = -(((1 + rC) * rt * (-1 + w)) / mSqq);
        result.element({8, 0, 6, 3}) = -((rt * (-1 + w) * (-1 + rC * w)) / mSqq);
        result.element({8, 0, 6, 4}) = -((rt * (rC - w) * (-1 + w)) / mSqq);
        result.element({8, 0, 7, 2}) = -(((1 + rC) * rt * (-1 + w)) / mSqq);
        result.element({8, 0, 7, 3}) = -((rt * (-1 + w) * (-1 + rC * w)) / mSqq);
        result.element({8, 0, 7, 4}) = -((rt * (rC - w) * (-1 + w)) / mSqq);
        result.element({8, 0, 8, 0}) = -1 + w;
        result.element({8, 0, 9, 0}) = -1 + w;
        result.element({8, 1, 6, 5}) = ((-1 + rC) * rt * (1 + w)) / mSqq;
        result.element({8, 1, 6, 6}) = (rt * (1 + w) * (-1 + rC * w)) / mSqq;
        result.element({8, 1, 6, 7}) = (rt * (rC - w) * (1 + w)) / mSqq;
        result.element({8, 1, 7, 5}) = -(((-1 + rC) * rt * (1 + w)) / mSqq);
        result.element({8, 1, 7, 6}) = -((rt * (1 + w) * (-1 + rC * w)) / mSqq);
        result.element({8, 1, 7, 7}) = -((rt * (rC - w) * (1 + w)) / mSqq);
        result.element({8, 1, 8, 1}) = 1 + w;
        result.element({8, 1, 9, 1}) = -1 - w;
        result.element({9, 0, 6, 2}) = -(((1 + rC) * rt * (-1 + w)) / mSqq);
        result.element({9, 0, 6, 3}) = -((rt * (-1 + w) * (-1 + rC * w)) / mSqq);
        result.element({9, 0, 6, 4}) = -((rt * (rC - w) * (-1 + w)) / mSqq);
        result.element({9, 0, 7, 2}) = -(((1 + rC) * rt * (-1 + w)) / mSqq);
        result.element({9, 0, 7, 3}) = -((rt * (-1 + w) * (-1 + rC * w)) / mSqq);
        result.element({9, 0, 7, 4}) = -((rt * (rC - w) * (-1 + w)) / mSqq);
        result.element({9, 0, 8, 0}) = -1 + w;
        result.element({9, 0, 9, 0}) = -1 + w;
        result.element({9, 1, 6, 5}) = -(((-1 + rC) * rt * (1 + w)) / mSqq);
        result.element({9, 1, 6, 6}) = -((rt * (1 + w) * (-1 + rC * w)) / mSqq);
        result.element({9, 1, 6, 7}) = -((rt * (rC - w) * (1 + w)) / mSqq);
        result.element({9, 1, 7, 5}) = ((-1 + rC) * rt * (1 + w)) / mSqq;
        result.element({9, 1, 7, 6}) = (rt * (1 + w) * (-1 + rC * w)) / mSqq;
        result.element({9, 1, 7, 7}) = (rt * (rC - w) * (1 + w)) / mSqq;
        result.element({9, 1, 8, 1}) = -1 - w;
        result.element({9, 1, 9, 1}) = 1 + w;
        result.element({10, 8, 6, 2}) = (12 * (-1 + rC) * rt * (1 + w)) / mSqq;
        result.element({10, 8, 6, 3}) = (4 * rC * rt * w2m1) / mSqq;
        result.element({10, 8, 6, 4}) = (4 * rt * w2m1) / mSqq;
        result.element({10, 8, 6, 5}) = (-12 * (1 + rC) * rt * (-1 + w)) / mSqq;
        result.element({10, 8, 6, 6}) = (-4 * rC * rt * w2m1) / mSqq;
        result.element({10, 8, 6, 7}) = (-4 * rt * w2m1) / mSqq;
        result.element({10, 8, 7, 2}) = (12 * (-1 + rC) * rt * (1 + w)) / mSqq;
        result.element({10, 8, 7, 3}) = (4 * rC * rt * w2m1) / mSqq;
        result.element({10, 8, 7, 4}) = (4 * rt * w2m1) / mSqq;
        result.element({10, 8, 7, 5}) = (12 * (1 + rC) * rt * (-1 + w)) / mSqq;
        result.element({10, 8, 7, 6}) = (4 * rC * rt * w2m1) / mSqq;
        result.element({10, 8, 7, 7}) = (4 * rt * w2m1) / mSqq;
        result.element({10, 8, 10, 8}) =
            (32 * (mSqq + 2 * rtSq) * (-4 * rC + (2 + mSqq) * w + 2 * rCSq * w)) / (3. * mSqq2);
        result.element({10, 8, 10, 9}) =
            (16 * (mSqq + 2 * rtSq) * (1 + w) * (2 + mSqq - 4 * rC + rCSq * (-2 + 4 * w))) / (3. * mSqq2);
        result.element({10, 8, 10, 10}) =
            (16 * (mSqq + 2 * rtSq) * (1 + w) * (-2 + mSqq - 4 * rC + 2 * rCSq + 4 * w)) / (3. * mSqq2);
        result.element({10, 8, 10, 11}) = (-16 * (mSqq + 2 * rtSq) * w2m1) / (3. * mSqq);
        result.element({10, 9, 6, 2}) = (4 * rt * (1 + w) * (-3 + rC + 2 * rC * w)) / mSqq;
        result.element({10, 9, 6, 3}) = (4 * rC * rt * w2m1) / mSqq;
        result.element({10, 9, 6, 4}) = (4 * rt * w2m1) / mSqq;
        result.element({10, 9, 6, 5}) = (-8 * rC * rt * w2m1) / mSqq;
        result.element({10, 9, 7, 2}) = (4 * rt * (1 + w) * (-3 + rC + 2 * rC * w)) / mSqq;
        result.element({10, 9, 7, 3}) = (4 * rC * rt * w2m1) / mSqq;
        result.element({10, 9, 7, 4}) = (4 * rt * w2m1) / mSqq;
        result.element({10, 9, 7, 5}) = (8 * rC * rt * w2m1) / mSqq;
        result.element({10, 9, 10, 8}) =
            (16 * (mSqq + 2 * rtSq) * (1 + w) * (2 + mSqq - 4 * rC + rCSq * (-2 + 4 * w))) / (3. * mSqq2);
        result.element({10, 9, 10, 9}) = (16 * (1 + w) *
                                          (mSqq2 + 4 * rtSq * (1 - 2 * rC * w + rCSq * (-1 + 2 * wSq)) +
                                           2 * mSqq * (1 + rtSq - 2 * rC * w + rCSq * (-1 + 2 * wSq)))) /
                                         (3. * mSqq2);
        result.element({10, 9, 10, 10}) =
            (16 * (mSqq + 2 * rtSq) * (1 + w) * (mSqq + 2 * (-2 * rC + w + rCSq * w))) / (3. * mSqq2);
        result.element({10, 9, 10, 11}) = (-16 * (mSqq + 2 * rtSq) * w2m1) / (3. * mSqq);
        result.element({10, 10, 6, 2}) = (4 * rt * (-1 + 3 * rC - 2 * w) * (1 + w)) / mSqq;
        result.element({10, 10, 6, 3}) = (4 * rC * rt * w2m1) / mSqq;
        result.element({10, 10, 6, 4}) = (4 * rt * w2m1) / mSqq;
        result.element({10, 10, 6, 5}) = (-8 * rt * w2m1) / mSqq;
        result.element({10, 10, 7, 2}) = (4 * rt * (-1 + 3 * rC - 2 * w) * (1 + w)) / mSqq;
        result.element({10, 10, 7, 3}) = (4 * rC * rt * w2m1) / mSqq;
        result.element({10, 10, 7, 4}) = (4 * rt * w2m1) / mSqq;
        result.element({10, 10, 7, 5}) = (8 * rt * w2m1) / mSqq;
        result.element({10, 10, 10, 8}) =
            (16 * (mSqq + 2 * rtSq) * (1 + w) * (-2 + mSqq - 4 * rC + 2 * rCSq + 4 * w)) / (3. * mSqq2);
        result.element({10, 10, 10, 9}) =
            (16 * (mSqq + 2 * rtSq) * (1 + w) * (mSqq + 2 * (-2 * rC + w + rCSq * w))) / (3. * mSqq2);
        result.element({10, 10, 10, 10}) = (16 * (1 + w) *
                                            (mSqq2 + 4 * rtSq * (-1 + rCSq - 2 * rC * w + 2 * wSq) +
                                             2 * mSqq * (-1 + rCSq + rtSq - 2 * rC * w + 2 * wSq))) /
                                           (3. * mSqq2);
        result.element({10, 10, 10, 11}) = (-16 * (mSqq + 2 * rtSq) * w2m1) / (3. * mSqq);
        result.element({10, 11, 6, 2}) = (-4 * (-1 + rC) * rt * w2m1) / mSqq;
        result.element({10, 11, 6, 3}) = (-4 * rC * rt * wm1Sq * (1 + w)) / mSqq;
        result.element({10, 11, 6, 4}) = (-4 * rt * wm1Sq * (1 + w)) / mSqq;
        result.element({10, 11, 7, 2}) = (-4 * (-1 + rC) * rt * w2m1) / mSqq;
        result.element({10, 11, 7, 3}) = (-4 * rC * rt * wm1Sq * (1 + w)) / mSqq;
        result.element({10, 11, 7, 4}) = (-4 * rt * wm1Sq * (1 + w)) / mSqq;
        result.element({10, 11, 10, 8}) = (-16 * (mSqq + 2 * rtSq) * w2m1) / (3. * mSqq);
        result.element({10, 11, 10, 9}) = (-16 * (mSqq + 2 * rtSq) * w2m1) / (3. * mSqq);
        result.element({10, 11, 10, 10}) = (-16 * (mSqq + 2 * rtSq) * w2m1) / (3. * mSqq);
        result.element({10, 11, 10, 11}) = (16 * (mSqq + 2 * rtSq) * wm1Sq * (1 + w)) / (3. * mSqq);

        result *= RateNorm;

        return result;
    }

} // namespace Hammer
