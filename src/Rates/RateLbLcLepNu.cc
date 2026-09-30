///
/// @file  RateLbLcLepNu.cc
/// @brief \f$ \Lambda_b \rightarrow \Lambda_c \tau\nu \f$ total rate
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/Rates/RateLbLcLepNu.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    RateLbLcLepNu::RateLbLcLepNu() {
        // Create tensor rank and dimensions
        IndexList dims{{11, 12, 11, 12, _nPoints}};
        string name{"RateLbLcLepNuQ2"};
        auto& pdg = PID::instance();

        addProcessSignature(-PID::LAMBDAB, {PID::LAMBDACMINUS, PID::NU_TAU, PID::ANTITAU});
        addIntegrationBoundaries({PS::makeQ2Function(pdg.getMass(PID::ANTITAU),
                                                     pdg.getMass(PID::LAMBDAB) - pdg.getMass(PID::LAMBDACMINUS))});
        addTensor(Tensor{name, MD::makeEmptySparse(
                                   dims, {WILSON_BCTAUNU, FF_LBLC, WILSON_BCTAUNU_HC, FF_LBLC_HC, INTEGRATION_INDEX})});

        addProcessSignature(-PID::LAMBDAB, {PID::LAMBDACMINUS, PID::NU_MU, PID::ANTIMUON});
        addIntegrationBoundaries({PS::makeQ2Function(pdg.getMass(PID::ANTIMUON),
                                                     pdg.getMass(PID::LAMBDAB) - pdg.getMass(PID::LAMBDACMINUS))});
        addTensor(Tensor{name, MD::makeEmptySparse(
                                   dims, {WILSON_BCMUNU, FF_LBLC, WILSON_BCMUNU_HC, FF_LBLC_HC, INTEGRATION_INDEX})});

        addProcessSignature(-PID::LAMBDAB, {PID::LAMBDACMINUS, PID::NU_E, PID::POSITRON});
        addIntegrationBoundaries({PS::makeQ2Function(pdg.getMass(PID::POSITRON),
                                                     pdg.getMass(PID::LAMBDAB) - pdg.getMass(PID::LAMBDACMINUS))});
        addTensor(Tensor{
            name, MD::makeEmptySparse(dims, {WILSON_BCENU, FF_LBLC, WILSON_BCENU_HC, FF_LBLC_HC, INTEGRATION_INDEX})});

        setSignatureIndex();
    }

    Tensor RateLbLcLepNu::evalAtPSPoint(const vector<double>& point) {
        auto labs = getTensor().labels();
        labs.pop_back();
        auto dimensions = getTensor().dims();
        dimensions.pop_back();
        Tensor result{"RateLbLcLepNu", MD::makeEmptySparse(dimensions, labs)};

        const double Mb = masses()[0];
        const double Mc = masses()[1];
        const double Mt = masses()[3];

        // kinematic objects
        const double Sqq = point[0];
        const double mSqq = Sqq / (Mb * Mb);
        const double w = (Mb * Mb + Mc * Mc - Sqq) / (2. * Mb * Mc);
        const double rC = Mc / Mb;
        const double rt = Mt / Mb;

        const double mSqq2 = pow(mSqq, 2.);
        const double wSq = w * w;
        const double rCSq = rC * rC;
        const double rtSq = pow(rt, 2.);

        const double RateNorm =
            (pow(GFermi, 2.) * pow(Mb, 3.) * pow(rC, 2.) * pow(mSqq - rtSq, 2.) * sqrt(-1 + wSq)) / (32. * mSqq * pi3);

        // set non-zero tensor elements
        result.element({0, 2, 0, 2}) = (2. * (2. * mSqq2 * (-1 + w) + mSqq * (1 + 2. * rC + rCSq + rtSq) * (-1 + w) +
                                              rtSq * (1 + 2. * w - 2. * rC * (2 + w) + rCSq * (1 + 2. * w)))) /
                                       (3. * mSqq2);
        result.element({0, 2, 0, 3}) =
            ((1 + w) * (2. * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2. * rC * (2 + w) + rCSq * (-1 + 4. * w)))) /
            (3. * mSqq2);
        result.element({0, 2, 0, 4}) =
            ((1 + w) * (2. * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3. * rCSq + 4. * w - 2. * rC * (2 + w)))) /
            (3. * mSqq2);
        result.element({0, 2, 1, 0}) = ((-1 + rC) * rt * (1 + w)) / mSqq;
        result.element({0, 2, 3, 0}) = ((-1 + rC) * rt * (1 + w)) / mSqq;
        result.element({0, 2, 5, 2}) = (2. * (2. * mSqq2 * (-1 + w) + mSqq * (1 + 2. * rC + rCSq + rtSq) * (-1 + w) +
                                              rtSq * (1 + 2. * w - 2. * rC * (2 + w) + rCSq * (1 + 2. * w)))) /
                                       (3. * mSqq2);
        result.element({0, 2, 5, 3}) =
            ((1 + w) * (2. * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2. * rC * (2 + w) + rCSq * (-1 + 4. * w)))) /
            (3. * mSqq2);
        result.element({0, 2, 5, 4}) =
            ((1 + w) * (2. * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3. * rCSq + 4. * w - 2. * rC * (2 + w)))) /
            (3. * mSqq2);
        result.element({0, 2, 7, 2}) = (2. * (2. * mSqq2 * (-1 + w) + mSqq * (1 + 2. * rC + rCSq + rtSq) * (-1 + w) +
                                              rtSq * (1 + 2. * w - 2. * rC * (2 + w) + rCSq * (1 + 2. * w)))) /
                                       (3. * mSqq2);
        result.element({0, 2, 7, 3}) =
            ((1 + w) * (2. * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2. * rC * (2 + w) + rCSq * (-1 + 4. * w)))) /
            (3. * mSqq2);
        result.element({0, 2, 7, 4}) =
            ((1 + w) * (2. * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3. * rCSq + 4. * w - 2. * rC * (2 + w)))) /
            (3. * mSqq2);
        result.element({0, 2, 9, 8}) = (-12. * (1 + rC) * rt * (-1 + w)) / mSqq;
        result.element({0, 2, 9, 9}) = (-4. * rt * (-1 + w) * (-3 + rC * (-1 + 2. * w))) / mSqq;
        result.element({0, 2, 9, 10}) = (-4. * rt * (1 + 3. * rC - 2. * w) * (-1 + w)) / mSqq;
        result.element({0, 2, 9, 11}) = (4. * (1 + rC) * rt * (-1 + wSq)) / mSqq;
        result.element({0, 3, 0, 2}) =
            ((1 + w) * (2. * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2. * rC * (2 + w) + rCSq * (-1 + 4. * w)))) /
            (3. * mSqq2);
        result.element({0, 3, 0, 3}) =
            ((1 + w) * (2. * mSqq * rCSq * (-1 + wSq) + rtSq * (3 - 6. * rC * w + rCSq * (-1 + 4. * wSq)))) /
            (3. * mSqq2);
        result.element({0, 3, 0, 4}) =
            ((1 + w) * (2. * mSqq * rC * (-1 + wSq) + rtSq * (3. * w + 3. * rCSq * w - 2. * rC * (2 + wSq)))) /
            (3. * mSqq2);
        result.element({0, 3, 1, 0}) = (rt * (1 + w) * (-1 + rC * w)) / mSqq;
        result.element({0, 3, 3, 0}) = (rt * (1 + w) * (-1 + rC * w)) / mSqq;
        result.element({0, 3, 5, 2}) =
            ((1 + w) * (2. * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2. * rC * (2 + w) + rCSq * (-1 + 4. * w)))) /
            (3. * mSqq2);
        result.element({0, 3, 5, 3}) =
            ((1 + w) * (2. * mSqq * rCSq * (-1 + wSq) + rtSq * (3 - 6. * rC * w + rCSq * (-1 + 4. * wSq)))) /
            (3. * mSqq2);
        result.element({0, 3, 5, 4}) =
            ((1 + w) * (2. * mSqq * rC * (-1 + wSq) + rtSq * (3. * w + 3. * rCSq * w - 2. * rC * (2 + wSq)))) /
            (3. * mSqq2);
        result.element({0, 3, 7, 2}) =
            ((1 + w) * (2. * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2. * rC * (2 + w) + rCSq * (-1 + 4. * w)))) /
            (3. * mSqq2);
        result.element({0, 3, 7, 3}) =
            ((1 + w) * (2. * mSqq * rCSq * (-1 + wSq) + rtSq * (3 - 6. * rC * w + rCSq * (-1 + 4. * wSq)))) /
            (3. * mSqq2);
        result.element({0, 3, 7, 4}) =
            ((1 + w) * (2. * mSqq * rC * (-1 + wSq) + rtSq * (3. * w + 3. * rCSq * w - 2. * rC * (2 + wSq)))) /
            (3. * mSqq2);
        result.element({0, 3, 9, 8}) = (-4. * rC * rt * (-1 + wSq)) / mSqq;
        result.element({0, 3, 9, 9}) = (4. * rC * rt * (-1 + wSq)) / mSqq;
        result.element({0, 3, 9, 10}) = (-4. * rC * rt * (-1 + wSq)) / mSqq;
        result.element({0, 3, 9, 11}) = (4. * rC * rt * (-1 + w) * pow(1 + w, 2.)) / mSqq;
        result.element({0, 4, 0, 2}) =
            ((1 + w) * (2. * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3. * rCSq + 4. * w - 2. * rC * (2 + w)))) /
            (3. * mSqq2);
        result.element({0, 4, 0, 3}) =
            ((1 + w) * (2. * mSqq * rC * (-1 + wSq) + rtSq * (3. * w + 3. * rCSq * w - 2. * rC * (2 + wSq)))) /
            (3. * mSqq2);
        result.element({0, 4, 0, 4}) =
            ((1 + w) * (2. * mSqq * (-1 + wSq) + rtSq * (-1 + 3. * rCSq - 6. * rC * w + 4. * wSq))) / (3. * mSqq2);
        result.element({0, 4, 1, 0}) = (rt * (rC - w) * (1 + w)) / mSqq;
        result.element({0, 4, 3, 0}) = (rt * (rC - w) * (1 + w)) / mSqq;
        result.element({0, 4, 5, 2}) =
            ((1 + w) * (2. * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3. * rCSq + 4. * w - 2. * rC * (2 + w)))) /
            (3. * mSqq2);
        result.element({0, 4, 5, 3}) =
            ((1 + w) * (2. * mSqq * rC * (-1 + wSq) + rtSq * (3. * w + 3. * rCSq * w - 2. * rC * (2 + wSq)))) /
            (3. * mSqq2);
        result.element({0, 4, 5, 4}) =
            ((1 + w) * (2. * mSqq * (-1 + wSq) + rtSq * (-1 + 3. * rCSq - 6. * rC * w + 4. * wSq))) / (3. * mSqq2);
        result.element({0, 4, 7, 2}) =
            ((1 + w) * (2. * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3. * rCSq + 4. * w - 2. * rC * (2 + w)))) /
            (3. * mSqq2);
        result.element({0, 4, 7, 3}) =
            ((1 + w) * (2. * mSqq * rC * (-1 + wSq) + rtSq * (3. * w + 3. * rCSq * w - 2. * rC * (2 + wSq)))) /
            (3. * mSqq2);
        result.element({0, 4, 7, 4}) =
            ((1 + w) * (2. * mSqq * (-1 + wSq) + rtSq * (-1 + 3. * rCSq - 6. * rC * w + 4. * wSq))) / (3. * mSqq2);
        result.element({0, 4, 9, 8}) = (-4. * rt * (-1 + wSq)) / mSqq;
        result.element({0, 4, 9, 9}) = (4. * rt * (-1 + wSq)) / mSqq;
        result.element({0, 4, 9, 10}) = (-4. * rt * (-1 + wSq)) / mSqq;
        result.element({0, 4, 9, 11}) = (4. * rt * (-1 + w) * pow(1 + w, 2.)) / mSqq;
        result.element({0, 5, 0, 5}) = (2. * (2. * mSqq2 * (1 + w) + mSqq * (1 - 2. * rC + rCSq + rtSq) * (1 + w) +
                                              rtSq * (-1 + 2. * rC * (-2 + w) + 2. * w + rCSq * (-1 + 2. * w)))) /
                                       (3. * mSqq2);
        result.element({0, 5, 0, 6}) = ((-1 + w) * (2. * mSqq * (-1 + rC) * rC * (1 + w) +
                                                    rtSq * (-3 + 2. * rC * (-2 + w) + rCSq * (1 + 4. * w)))) /
                                       (3. * mSqq2);
        result.element({0, 5, 0, 7}) =
            ((-1 + w) * (rtSq * (-1 + 3. * rCSq - 2. * rC * (-2 + w) - 4. * w) + 2. * mSqq * (-1 + rC) * (1 + w))) /
            (3. * mSqq2);
        result.element({0, 5, 1, 1}) = ((1 + rC) * rt * (-1 + w)) / mSqq;
        result.element({0, 5, 3, 1}) = -(((1 + rC) * rt * (-1 + w)) / mSqq);
        result.element({0, 5, 5, 5}) = (-2. * (2. * mSqq2 * (1 + w) + mSqq * (1 - 2. * rC + rCSq + rtSq) * (1 + w) +
                                               rtSq * (-1 + 2. * rC * (-2 + w) + 2. * w + rCSq * (-1 + 2. * w)))) /
                                       (3. * mSqq2);
        result.element({0, 5, 5, 6}) = -((-1 + w) * (2. * mSqq * (-1 + rC) * rC * (1 + w) +
                                                     rtSq * (-3 + 2. * rC * (-2 + w) + rCSq * (1 + 4. * w)))) /
                                       (3. * mSqq2);
        result.element({0, 5, 5, 7}) =
            -((-1 + w) * (rtSq * (-1 + 3. * rCSq - 2. * rC * (-2 + w) - 4. * w) + 2. * mSqq * (-1 + rC) * (1 + w))) /
            (3. * mSqq2);
        result.element({0, 5, 7, 5}) = (2. * (2. * mSqq2 * (1 + w) + mSqq * (1 - 2. * rC + rCSq + rtSq) * (1 + w) +
                                              rtSq * (-1 + 2. * rC * (-2 + w) + 2. * w + rCSq * (-1 + 2. * w)))) /
                                       (3. * mSqq2);
        result.element({0, 5, 7, 6}) = ((-1 + w) * (2. * mSqq * (-1 + rC) * rC * (1 + w) +
                                                    rtSq * (-3 + 2. * rC * (-2 + w) + rCSq * (1 + 4. * w)))) /
                                       (3. * mSqq2);
        result.element({0, 5, 7, 7}) =
            ((-1 + w) * (rtSq * (-1 + 3. * rCSq - 2. * rC * (-2 + w) - 4. * w) + 2. * mSqq * (-1 + rC) * (1 + w))) /
            (3. * mSqq2);
        result.element({0, 5, 9, 8}) = (-12. * (-1 + rC) * rt * (1 + w)) / mSqq;
        result.element({0, 5, 9, 9}) = (-8. * rC * rt * (-1 + wSq)) / mSqq;
        result.element({0, 5, 9, 10}) = (-8. * rt * (-1 + wSq)) / mSqq;
        result.element({0, 6, 0, 5}) = ((-1 + w) * (2. * mSqq * (-1 + rC) * rC * (1 + w) +
                                                    rtSq * (-3 + 2. * rC * (-2 + w) + rCSq * (1 + 4. * w)))) /
                                       (3. * mSqq2);
        result.element({0, 6, 0, 6}) =
            ((-1 + w) * (2. * mSqq * rCSq * (-1 + wSq) + rtSq * (3 - 6. * rC * w + rCSq * (-1 + 4. * wSq)))) /
            (3. * mSqq2);
        result.element({0, 6, 0, 7}) =
            ((-1 + w) * (2. * mSqq * rC * (-1 + wSq) + rtSq * (3. * w + 3. * rCSq * w - 2. * rC * (2 + wSq)))) /
            (3. * mSqq2);
        result.element({0, 6, 1, 1}) = (rt * (-1 + w) * (-1 + rC * w)) / mSqq;
        result.element({0, 6, 3, 1}) = -((rt * (-1 + w) * (-1 + rC * w)) / mSqq);
        result.element({0, 6, 5, 5}) = -((-1 + w) * (2. * mSqq * (-1 + rC) * rC * (1 + w) +
                                                     rtSq * (-3 + 2. * rC * (-2 + w) + rCSq * (1 + 4. * w)))) /
                                       (3. * mSqq2);
        result.element({0, 6, 5, 6}) =
            -((-1 + w) * (2. * mSqq * rCSq * (-1 + wSq) + rtSq * (3 - 6. * rC * w + rCSq * (-1 + 4. * wSq)))) /
            (3. * mSqq2);
        result.element({0, 6, 5, 7}) =
            -((-1 + w) * (2. * mSqq * rC * (-1 + wSq) + rtSq * (3. * w + 3. * rCSq * w - 2. * rC * (2 + wSq)))) /
            (3. * mSqq2);
        result.element({0, 6, 7, 5}) = ((-1 + w) * (2. * mSqq * (-1 + rC) * rC * (1 + w) +
                                                    rtSq * (-3 + 2. * rC * (-2 + w) + rCSq * (1 + 4. * w)))) /
                                       (3. * mSqq2);
        result.element({0, 6, 7, 6}) =
            ((-1 + w) * (2. * mSqq * rCSq * (-1 + wSq) + rtSq * (3 - 6. * rC * w + rCSq * (-1 + 4. * wSq)))) /
            (3. * mSqq2);
        result.element({0, 6, 7, 7}) =
            ((-1 + w) * (2. * mSqq * rC * (-1 + wSq) + rtSq * (3. * w + 3. * rCSq * w - 2. * rC * (2 + wSq)))) /
            (3. * mSqq2);
        result.element({0, 6, 9, 8}) = (-4. * rC * rt * (-1 + wSq)) / mSqq;
        result.element({0, 7, 0, 5}) =
            ((-1 + w) * (rtSq * (-1 + 3. * rCSq - 2. * rC * (-2 + w) - 4. * w) + 2. * mSqq * (-1 + rC) * (1 + w))) /
            (3. * mSqq2);
        result.element({0, 7, 0, 6}) =
            ((-1 + w) * (2. * mSqq * rC * (-1 + wSq) + rtSq * (3. * w + 3. * rCSq * w - 2. * rC * (2 + wSq)))) /
            (3. * mSqq2);
        result.element({0, 7, 0, 7}) =
            ((-1 + w) * (2. * mSqq * (-1 + wSq) + rtSq * (-1 + 3. * rCSq - 6. * rC * w + 4. * wSq))) / (3. * mSqq2);
        result.element({0, 7, 1, 1}) = (rt * (rC - w) * (-1 + w)) / mSqq;
        result.element({0, 7, 3, 1}) = -((rt * (rC - w) * (-1 + w)) / mSqq);
        result.element({0, 7, 5, 5}) =
            -((-1 + w) * (rtSq * (-1 + 3. * rCSq - 2. * rC * (-2 + w) - 4. * w) + 2. * mSqq * (-1 + rC) * (1 + w))) /
            (3. * mSqq2);
        result.element({0, 7, 5, 6}) =
            -((-1 + w) * (2. * mSqq * rC * (-1 + wSq) + rtSq * (3. * w + 3. * rCSq * w - 2. * rC * (2 + wSq)))) /
            (3. * mSqq2);
        result.element({0, 7, 5, 7}) =
            -((-1 + w) * (2. * mSqq * (-1 + wSq) + rtSq * (-1 + 3. * rCSq - 6. * rC * w + 4. * wSq))) / (3. * mSqq2);
        result.element({0, 7, 7, 5}) =
            ((-1 + w) * (rtSq * (-1 + 3. * rCSq - 2. * rC * (-2 + w) - 4. * w) + 2. * mSqq * (-1 + rC) * (1 + w))) /
            (3. * mSqq2);
        result.element({0, 7, 7, 6}) =
            ((-1 + w) * (2. * mSqq * rC * (-1 + wSq) + rtSq * (3. * w + 3. * rCSq * w - 2. * rC * (2 + wSq)))) /
            (3. * mSqq2);
        result.element({0, 7, 7, 7}) =
            ((-1 + w) * (2. * mSqq * (-1 + wSq) + rtSq * (-1 + 3. * rCSq - 6. * rC * w + 4. * wSq))) / (3. * mSqq2);
        result.element({0, 7, 9, 8}) = (-4. * rt * (-1 + wSq)) / mSqq;
        result.element({1, 0, 0, 2}) = ((-1 + rC) * rt * (1 + w)) / mSqq;
        result.element({1, 0, 0, 3}) = (rt * (1 + w) * (-1 + rC * w)) / mSqq;
        result.element({1, 0, 0, 4}) = (rt * (rC - w) * (1 + w)) / mSqq;
        result.element({1, 0, 1, 0}) = 1 + w;
        result.element({1, 0, 3, 0}) = 1 + w;
        result.element({1, 0, 5, 2}) = ((-1 + rC) * rt * (1 + w)) / mSqq;
        result.element({1, 0, 5, 3}) = (rt * (1 + w) * (-1 + rC * w)) / mSqq;
        result.element({1, 0, 5, 4}) = (rt * (rC - w) * (1 + w)) / mSqq;
        result.element({1, 0, 7, 2}) = ((-1 + rC) * rt * (1 + w)) / mSqq;
        result.element({1, 0, 7, 3}) = (rt * (1 + w) * (-1 + rC * w)) / mSqq;
        result.element({1, 0, 7, 4}) = (rt * (rC - w) * (1 + w)) / mSqq;
        result.element({1, 1, 0, 5}) = ((1 + rC) * rt * (-1 + w)) / mSqq;
        result.element({1, 1, 0, 6}) = (rt * (-1 + w) * (-1 + rC * w)) / mSqq;
        result.element({1, 1, 0, 7}) = (rt * (rC - w) * (-1 + w)) / mSqq;
        result.element({1, 1, 1, 1}) = -1 + w;
        result.element({1, 1, 3, 1}) = 1 - w;
        result.element({1, 1, 5, 5}) = -(((1 + rC) * rt * (-1 + w)) / mSqq);
        result.element({1, 1, 5, 6}) = -((rt * (-1 + w) * (-1 + rC * w)) / mSqq);
        result.element({1, 1, 5, 7}) = -((rt * (rC - w) * (-1 + w)) / mSqq);
        result.element({1, 1, 7, 5}) = ((1 + rC) * rt * (-1 + w)) / mSqq;
        result.element({1, 1, 7, 6}) = (rt * (-1 + w) * (-1 + rC * w)) / mSqq;
        result.element({1, 1, 7, 7}) = (rt * (rC - w) * (-1 + w)) / mSqq;
        result.element({2, 0, 2, 0}) = 1 + w;
        result.element({2, 0, 4, 0}) = 1 + w;
        result.element({2, 0, 6, 2}) = ((-1 + rC) * rt * (1 + w)) / mSqq;
        result.element({2, 0, 6, 3}) = (rt * (1 + w) * (-1 + rC * w)) / mSqq;
        result.element({2, 0, 6, 4}) = (rt * (rC - w) * (1 + w)) / mSqq;
        result.element({2, 0, 8, 2}) = ((-1 + rC) * rt * (1 + w)) / mSqq;
        result.element({2, 0, 8, 3}) = (rt * (1 + w) * (-1 + rC * w)) / mSqq;
        result.element({2, 0, 8, 4}) = (rt * (rC - w) * (1 + w)) / mSqq;
        result.element({2, 1, 2, 1}) = -1 + w;
        result.element({2, 1, 4, 1}) = 1 - w;
        result.element({2, 1, 6, 5}) = -(((1 + rC) * rt * (-1 + w)) / mSqq);
        result.element({2, 1, 6, 6}) = -((rt * (-1 + w) * (-1 + rC * w)) / mSqq);
        result.element({2, 1, 6, 7}) = -((rt * (rC - w) * (-1 + w)) / mSqq);
        result.element({2, 1, 8, 5}) = ((1 + rC) * rt * (-1 + w)) / mSqq;
        result.element({2, 1, 8, 6}) = (rt * (-1 + w) * (-1 + rC * w)) / mSqq;
        result.element({2, 1, 8, 7}) = (rt * (rC - w) * (-1 + w)) / mSqq;
        result.element({3, 0, 0, 2}) = ((-1 + rC) * rt * (1 + w)) / mSqq;
        result.element({3, 0, 0, 3}) = (rt * (1 + w) * (-1 + rC * w)) / mSqq;
        result.element({3, 0, 0, 4}) = (rt * (rC - w) * (1 + w)) / mSqq;
        result.element({3, 0, 1, 0}) = 1 + w;
        result.element({3, 0, 3, 0}) = 1 + w;
        result.element({3, 0, 5, 2}) = ((-1 + rC) * rt * (1 + w)) / mSqq;
        result.element({3, 0, 5, 3}) = (rt * (1 + w) * (-1 + rC * w)) / mSqq;
        result.element({3, 0, 5, 4}) = (rt * (rC - w) * (1 + w)) / mSqq;
        result.element({3, 0, 7, 2}) = ((-1 + rC) * rt * (1 + w)) / mSqq;
        result.element({3, 0, 7, 3}) = (rt * (1 + w) * (-1 + rC * w)) / mSqq;
        result.element({3, 0, 7, 4}) = (rt * (rC - w) * (1 + w)) / mSqq;
        result.element({3, 1, 0, 5}) = -(((1 + rC) * rt * (-1 + w)) / mSqq);
        result.element({3, 1, 0, 6}) = -((rt * (-1 + w) * (-1 + rC * w)) / mSqq);
        result.element({3, 1, 0, 7}) = -((rt * (rC - w) * (-1 + w)) / mSqq);
        result.element({3, 1, 1, 1}) = 1 - w;
        result.element({3, 1, 3, 1}) = -1 + w;
        result.element({3, 1, 5, 5}) = ((1 + rC) * rt * (-1 + w)) / mSqq;
        result.element({3, 1, 5, 6}) = (rt * (-1 + w) * (-1 + rC * w)) / mSqq;
        result.element({3, 1, 5, 7}) = (rt * (rC - w) * (-1 + w)) / mSqq;
        result.element({3, 1, 7, 5}) = -(((1 + rC) * rt * (-1 + w)) / mSqq);
        result.element({3, 1, 7, 6}) = -((rt * (-1 + w) * (-1 + rC * w)) / mSqq);
        result.element({3, 1, 7, 7}) = -((rt * (rC - w) * (-1 + w)) / mSqq);
        result.element({4, 0, 2, 0}) = 1 + w;
        result.element({4, 0, 4, 0}) = 1 + w;
        result.element({4, 0, 6, 2}) = ((-1 + rC) * rt * (1 + w)) / mSqq;
        result.element({4, 0, 6, 3}) = (rt * (1 + w) * (-1 + rC * w)) / mSqq;
        result.element({4, 0, 6, 4}) = (rt * (rC - w) * (1 + w)) / mSqq;
        result.element({4, 0, 8, 2}) = ((-1 + rC) * rt * (1 + w)) / mSqq;
        result.element({4, 0, 8, 3}) = (rt * (1 + w) * (-1 + rC * w)) / mSqq;
        result.element({4, 0, 8, 4}) = (rt * (rC - w) * (1 + w)) / mSqq;
        result.element({4, 1, 2, 1}) = 1 - w;
        result.element({4, 1, 4, 1}) = -1 + w;
        result.element({4, 1, 6, 5}) = ((1 + rC) * rt * (-1 + w)) / mSqq;
        result.element({4, 1, 6, 6}) = (rt * (-1 + w) * (-1 + rC * w)) / mSqq;
        result.element({4, 1, 6, 7}) = (rt * (rC - w) * (-1 + w)) / mSqq;
        result.element({4, 1, 8, 5}) = -(((1 + rC) * rt * (-1 + w)) / mSqq);
        result.element({4, 1, 8, 6}) = -((rt * (-1 + w) * (-1 + rC * w)) / mSqq);
        result.element({4, 1, 8, 7}) = -((rt * (rC - w) * (-1 + w)) / mSqq);
        result.element({5, 2, 0, 2}) = (2. * (2. * mSqq2 * (-1 + w) + mSqq * (1 + 2. * rC + rCSq + rtSq) * (-1 + w) +
                                              rtSq * (1 + 2. * w - 2. * rC * (2 + w) + rCSq * (1 + 2. * w)))) /
                                       (3. * mSqq2);
        result.element({5, 2, 0, 3}) =
            ((1 + w) * (2. * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2. * rC * (2 + w) + rCSq * (-1 + 4. * w)))) /
            (3. * mSqq2);
        result.element({5, 2, 0, 4}) =
            ((1 + w) * (2. * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3. * rCSq + 4. * w - 2. * rC * (2 + w)))) /
            (3. * mSqq2);
        result.element({5, 2, 1, 0}) = ((-1 + rC) * rt * (1 + w)) / mSqq;
        result.element({5, 2, 3, 0}) = ((-1 + rC) * rt * (1 + w)) / mSqq;
        result.element({5, 2, 5, 2}) = (2. * (2. * mSqq2 * (-1 + w) + mSqq * (1 + 2. * rC + rCSq + rtSq) * (-1 + w) +
                                              rtSq * (1 + 2. * w - 2. * rC * (2 + w) + rCSq * (1 + 2. * w)))) /
                                       (3. * mSqq2);
        result.element({5, 2, 5, 3}) =
            ((1 + w) * (2. * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2. * rC * (2 + w) + rCSq * (-1 + 4. * w)))) /
            (3. * mSqq2);
        result.element({5, 2, 5, 4}) =
            ((1 + w) * (2. * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3. * rCSq + 4. * w - 2. * rC * (2 + w)))) /
            (3. * mSqq2);
        result.element({5, 2, 7, 2}) = (2. * (2. * mSqq2 * (-1 + w) + mSqq * (1 + 2. * rC + rCSq + rtSq) * (-1 + w) +
                                              rtSq * (1 + 2. * w - 2. * rC * (2 + w) + rCSq * (1 + 2. * w)))) /
                                       (3. * mSqq2);
        result.element({5, 2, 7, 3}) =
            ((1 + w) * (2. * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2. * rC * (2 + w) + rCSq * (-1 + 4. * w)))) /
            (3. * mSqq2);
        result.element({5, 2, 7, 4}) =
            ((1 + w) * (2. * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3. * rCSq + 4. * w - 2. * rC * (2 + w)))) /
            (3. * mSqq2);
        result.element({5, 2, 9, 8}) = (-12. * (1 + rC) * rt * (-1 + w)) / mSqq;
        result.element({5, 2, 9, 9}) = (-4. * rt * (-1 + w) * (-3 + rC * (-1 + 2. * w))) / mSqq;
        result.element({5, 2, 9, 10}) = (-4. * rt * (1 + 3. * rC - 2. * w) * (-1 + w)) / mSqq;
        result.element({5, 2, 9, 11}) = (4. * (1 + rC) * rt * (-1 + wSq)) / mSqq;
        result.element({5, 3, 0, 2}) =
            ((1 + w) * (2. * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2. * rC * (2 + w) + rCSq * (-1 + 4. * w)))) /
            (3. * mSqq2);
        result.element({5, 3, 0, 3}) =
            ((1 + w) * (2. * mSqq * rCSq * (-1 + wSq) + rtSq * (3 - 6. * rC * w + rCSq * (-1 + 4. * wSq)))) /
            (3. * mSqq2);
        result.element({5, 3, 0, 4}) =
            ((1 + w) * (2. * mSqq * rC * (-1 + wSq) + rtSq * (3. * w + 3. * rCSq * w - 2. * rC * (2 + wSq)))) /
            (3. * mSqq2);
        result.element({5, 3, 1, 0}) = (rt * (1 + w) * (-1 + rC * w)) / mSqq;
        result.element({5, 3, 3, 0}) = (rt * (1 + w) * (-1 + rC * w)) / mSqq;
        result.element({5, 3, 5, 2}) =
            ((1 + w) * (2. * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2. * rC * (2 + w) + rCSq * (-1 + 4. * w)))) /
            (3. * mSqq2);
        result.element({5, 3, 5, 3}) =
            ((1 + w) * (2. * mSqq * rCSq * (-1 + wSq) + rtSq * (3 - 6. * rC * w + rCSq * (-1 + 4. * wSq)))) /
            (3. * mSqq2);
        result.element({5, 3, 5, 4}) =
            ((1 + w) * (2. * mSqq * rC * (-1 + wSq) + rtSq * (3. * w + 3. * rCSq * w - 2. * rC * (2 + wSq)))) /
            (3. * mSqq2);
        result.element({5, 3, 7, 2}) =
            ((1 + w) * (2. * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2. * rC * (2 + w) + rCSq * (-1 + 4. * w)))) /
            (3. * mSqq2);
        result.element({5, 3, 7, 3}) =
            ((1 + w) * (2. * mSqq * rCSq * (-1 + wSq) + rtSq * (3 - 6. * rC * w + rCSq * (-1 + 4. * wSq)))) /
            (3. * mSqq2);
        result.element({5, 3, 7, 4}) =
            ((1 + w) * (2. * mSqq * rC * (-1 + wSq) + rtSq * (3. * w + 3. * rCSq * w - 2. * rC * (2 + wSq)))) /
            (3. * mSqq2);
        result.element({5, 3, 9, 8}) = (-4. * rC * rt * (-1 + wSq)) / mSqq;
        result.element({5, 3, 9, 9}) = (4. * rC * rt * (-1 + wSq)) / mSqq;
        result.element({5, 3, 9, 10}) = (-4. * rC * rt * (-1 + wSq)) / mSqq;
        result.element({5, 3, 9, 11}) = (4. * rC * rt * (-1 + w) * pow(1 + w, 2.)) / mSqq;
        result.element({5, 4, 0, 2}) =
            ((1 + w) * (2. * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3. * rCSq + 4. * w - 2. * rC * (2 + w)))) /
            (3. * mSqq2);
        result.element({5, 4, 0, 3}) =
            ((1 + w) * (2. * mSqq * rC * (-1 + wSq) + rtSq * (3. * w + 3. * rCSq * w - 2. * rC * (2 + wSq)))) /
            (3. * mSqq2);
        result.element({5, 4, 0, 4}) =
            ((1 + w) * (2. * mSqq * (-1 + wSq) + rtSq * (-1 + 3. * rCSq - 6. * rC * w + 4. * wSq))) / (3. * mSqq2);
        result.element({5, 4, 1, 0}) = (rt * (rC - w) * (1 + w)) / mSqq;
        result.element({5, 4, 3, 0}) = (rt * (rC - w) * (1 + w)) / mSqq;
        result.element({5, 4, 5, 2}) =
            ((1 + w) * (2. * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3. * rCSq + 4. * w - 2. * rC * (2 + w)))) /
            (3. * mSqq2);
        result.element({5, 4, 5, 3}) =
            ((1 + w) * (2. * mSqq * rC * (-1 + wSq) + rtSq * (3. * w + 3. * rCSq * w - 2. * rC * (2 + wSq)))) /
            (3. * mSqq2);
        result.element({5, 4, 5, 4}) =
            ((1 + w) * (2. * mSqq * (-1 + wSq) + rtSq * (-1 + 3. * rCSq - 6. * rC * w + 4. * wSq))) / (3. * mSqq2);
        result.element({5, 4, 7, 2}) =
            ((1 + w) * (2. * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3. * rCSq + 4. * w - 2. * rC * (2 + w)))) /
            (3. * mSqq2);
        result.element({5, 4, 7, 3}) =
            ((1 + w) * (2. * mSqq * rC * (-1 + wSq) + rtSq * (3. * w + 3. * rCSq * w - 2. * rC * (2 + wSq)))) /
            (3. * mSqq2);
        result.element({5, 4, 7, 4}) =
            ((1 + w) * (2. * mSqq * (-1 + wSq) + rtSq * (-1 + 3. * rCSq - 6. * rC * w + 4. * wSq))) / (3. * mSqq2);
        result.element({5, 4, 9, 8}) = (-4. * rt * (-1 + wSq)) / mSqq;
        result.element({5, 4, 9, 9}) = (4. * rt * (-1 + wSq)) / mSqq;
        result.element({5, 4, 9, 10}) = (-4. * rt * (-1 + wSq)) / mSqq;
        result.element({5, 4, 9, 11}) = (4. * rt * (-1 + w) * pow(1 + w, 2.)) / mSqq;
        result.element({5, 5, 0, 5}) = (-2. * (2. * mSqq2 * (1 + w) + mSqq * (1 - 2. * rC + rCSq + rtSq) * (1 + w) +
                                               rtSq * (-1 + 2. * rC * (-2 + w) + 2. * w + rCSq * (-1 + 2. * w)))) /
                                       (3. * mSqq2);
        result.element({5, 5, 0, 6}) = -((-1 + w) * (2. * mSqq * (-1 + rC) * rC * (1 + w) +
                                                     rtSq * (-3 + 2. * rC * (-2 + w) + rCSq * (1 + 4. * w)))) /
                                       (3. * mSqq2);
        result.element({5, 5, 0, 7}) =
            -((-1 + w) * (rtSq * (-1 + 3. * rCSq - 2. * rC * (-2 + w) - 4. * w) + 2. * mSqq * (-1 + rC) * (1 + w))) /
            (3. * mSqq2);
        result.element({5, 5, 1, 1}) = -(((1 + rC) * rt * (-1 + w)) / mSqq);
        result.element({5, 5, 3, 1}) = ((1 + rC) * rt * (-1 + w)) / mSqq;
        result.element({5, 5, 5, 5}) = (2. * (2. * mSqq2 * (1 + w) + mSqq * (1 - 2. * rC + rCSq + rtSq) * (1 + w) +
                                              rtSq * (-1 + 2. * rC * (-2 + w) + 2. * w + rCSq * (-1 + 2. * w)))) /
                                       (3. * mSqq2);
        result.element({5, 5, 5, 6}) = ((-1 + w) * (2. * mSqq * (-1 + rC) * rC * (1 + w) +
                                                    rtSq * (-3 + 2. * rC * (-2 + w) + rCSq * (1 + 4. * w)))) /
                                       (3. * mSqq2);
        result.element({5, 5, 5, 7}) =
            ((-1 + w) * (rtSq * (-1 + 3. * rCSq - 2. * rC * (-2 + w) - 4. * w) + 2. * mSqq * (-1 + rC) * (1 + w))) /
            (3. * mSqq2);
        result.element({5, 5, 7, 5}) = (-2. * (2. * mSqq2 * (1 + w) + mSqq * (1 - 2. * rC + rCSq + rtSq) * (1 + w) +
                                               rtSq * (-1 + 2. * rC * (-2 + w) + 2. * w + rCSq * (-1 + 2. * w)))) /
                                       (3. * mSqq2);
        result.element({5, 5, 7, 6}) = -((-1 + w) * (2. * mSqq * (-1 + rC) * rC * (1 + w) +
                                                     rtSq * (-3 + 2. * rC * (-2 + w) + rCSq * (1 + 4. * w)))) /
                                       (3. * mSqq2);
        result.element({5, 5, 7, 7}) =
            -((-1 + w) * (rtSq * (-1 + 3. * rCSq - 2. * rC * (-2 + w) - 4. * w) + 2. * mSqq * (-1 + rC) * (1 + w))) /
            (3. * mSqq2);
        result.element({5, 5, 9, 8}) = (12. * (-1 + rC) * rt * (1 + w)) / mSqq;
        result.element({5, 5, 9, 9}) = (8. * rC * rt * (-1 + wSq)) / mSqq;
        result.element({5, 5, 9, 10}) = (8. * rt * (-1 + wSq)) / mSqq;
        result.element({5, 6, 0, 5}) = -((-1 + w) * (2. * mSqq * (-1 + rC) * rC * (1 + w) +
                                                     rtSq * (-3 + 2. * rC * (-2 + w) + rCSq * (1 + 4. * w)))) /
                                       (3. * mSqq2);
        result.element({5, 6, 0, 6}) =
            -((-1 + w) * (2. * mSqq * rCSq * (-1 + wSq) + rtSq * (3 - 6. * rC * w + rCSq * (-1 + 4. * wSq)))) /
            (3. * mSqq2);
        result.element({5, 6, 0, 7}) =
            -((-1 + w) * (2. * mSqq * rC * (-1 + wSq) + rtSq * (3. * w + 3. * rCSq * w - 2. * rC * (2 + wSq)))) /
            (3. * mSqq2);
        result.element({5, 6, 1, 1}) = -((rt * (-1 + w) * (-1 + rC * w)) / mSqq);
        result.element({5, 6, 3, 1}) = (rt * (-1 + w) * (-1 + rC * w)) / mSqq;
        result.element({5, 6, 5, 5}) = ((-1 + w) * (2. * mSqq * (-1 + rC) * rC * (1 + w) +
                                                    rtSq * (-3 + 2. * rC * (-2 + w) + rCSq * (1 + 4. * w)))) /
                                       (3. * mSqq2);
        result.element({5, 6, 5, 6}) =
            ((-1 + w) * (2. * mSqq * rCSq * (-1 + wSq) + rtSq * (3 - 6. * rC * w + rCSq * (-1 + 4. * wSq)))) /
            (3. * mSqq2);
        result.element({5, 6, 5, 7}) =
            ((-1 + w) * (2. * mSqq * rC * (-1 + wSq) + rtSq * (3. * w + 3. * rCSq * w - 2. * rC * (2 + wSq)))) /
            (3. * mSqq2);
        result.element({5, 6, 7, 5}) = -((-1 + w) * (2. * mSqq * (-1 + rC) * rC * (1 + w) +
                                                     rtSq * (-3 + 2. * rC * (-2 + w) + rCSq * (1 + 4. * w)))) /
                                       (3. * mSqq2);
        result.element({5, 6, 7, 6}) =
            -((-1 + w) * (2. * mSqq * rCSq * (-1 + wSq) + rtSq * (3 - 6. * rC * w + rCSq * (-1 + 4. * wSq)))) /
            (3. * mSqq2);
        result.element({5, 6, 7, 7}) =
            -((-1 + w) * (2. * mSqq * rC * (-1 + wSq) + rtSq * (3. * w + 3. * rCSq * w - 2. * rC * (2 + wSq)))) /
            (3. * mSqq2);
        result.element({5, 6, 9, 8}) = (4. * rC * rt * (-1 + wSq)) / mSqq;
        result.element({5, 7, 0, 5}) =
            -((-1 + w) * (rtSq * (-1 + 3. * rCSq - 2. * rC * (-2 + w) - 4. * w) + 2. * mSqq * (-1 + rC) * (1 + w))) /
            (3. * mSqq2);
        result.element({5, 7, 0, 6}) =
            -((-1 + w) * (2. * mSqq * rC * (-1 + wSq) + rtSq * (3. * w + 3. * rCSq * w - 2. * rC * (2 + wSq)))) /
            (3. * mSqq2);
        result.element({5, 7, 0, 7}) =
            -((-1 + w) * (2. * mSqq * (-1 + wSq) + rtSq * (-1 + 3. * rCSq - 6. * rC * w + 4. * wSq))) / (3. * mSqq2);
        result.element({5, 7, 1, 1}) = -((rt * (rC - w) * (-1 + w)) / mSqq);
        result.element({5, 7, 3, 1}) = (rt * (rC - w) * (-1 + w)) / mSqq;
        result.element({5, 7, 5, 5}) =
            ((-1 + w) * (rtSq * (-1 + 3. * rCSq - 2. * rC * (-2 + w) - 4. * w) + 2. * mSqq * (-1 + rC) * (1 + w))) /
            (3. * mSqq2);
        result.element({5, 7, 5, 6}) =
            ((-1 + w) * (2. * mSqq * rC * (-1 + wSq) + rtSq * (3. * w + 3. * rCSq * w - 2. * rC * (2 + wSq)))) /
            (3. * mSqq2);
        result.element({5, 7, 5, 7}) =
            ((-1 + w) * (2. * mSqq * (-1 + wSq) + rtSq * (-1 + 3. * rCSq - 6. * rC * w + 4. * wSq))) / (3. * mSqq2);
        result.element({5, 7, 7, 5}) =
            -((-1 + w) * (rtSq * (-1 + 3. * rCSq - 2. * rC * (-2 + w) - 4. * w) + 2. * mSqq * (-1 + rC) * (1 + w))) /
            (3. * mSqq2);
        result.element({5, 7, 7, 6}) =
            -((-1 + w) * (2. * mSqq * rC * (-1 + wSq) + rtSq * (3. * w + 3. * rCSq * w - 2. * rC * (2 + wSq)))) /
            (3. * mSqq2);
        result.element({5, 7, 7, 7}) =
            -((-1 + w) * (2. * mSqq * (-1 + wSq) + rtSq * (-1 + 3. * rCSq - 6. * rC * w + 4. * wSq))) / (3. * mSqq2);
        result.element({5, 7, 9, 8}) = (4. * rt * (-1 + wSq)) / mSqq;
        result.element({6, 2, 2, 0}) = ((-1 + rC) * rt * (1 + w)) / mSqq;
        result.element({6, 2, 4, 0}) = ((-1 + rC) * rt * (1 + w)) / mSqq;
        result.element({6, 2, 6, 2}) = (2. * (2. * mSqq2 * (-1 + w) + mSqq * (1 + 2. * rC + rCSq + rtSq) * (-1 + w) +
                                              rtSq * (1 + 2. * w - 2. * rC * (2 + w) + rCSq * (1 + 2. * w)))) /
                                       (3. * mSqq2);
        result.element({6, 2, 6, 3}) =
            ((1 + w) * (2. * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2. * rC * (2 + w) + rCSq * (-1 + 4. * w)))) /
            (3. * mSqq2);
        result.element({6, 2, 6, 4}) =
            ((1 + w) * (2. * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3. * rCSq + 4. * w - 2. * rC * (2 + w)))) /
            (3. * mSqq2);
        result.element({6, 2, 8, 2}) = (2. * (2. * mSqq2 * (-1 + w) + mSqq * (1 + 2. * rC + rCSq + rtSq) * (-1 + w) +
                                              rtSq * (1 + 2. * w - 2. * rC * (2 + w) + rCSq * (1 + 2. * w)))) /
                                       (3. * mSqq2);
        result.element({6, 2, 8, 3}) =
            ((1 + w) * (2. * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2. * rC * (2 + w) + rCSq * (-1 + 4. * w)))) /
            (3. * mSqq2);
        result.element({6, 2, 8, 4}) =
            ((1 + w) * (2. * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3. * rCSq + 4. * w - 2. * rC * (2 + w)))) /
            (3. * mSqq2);
        result.element({6, 2, 10, 8}) = (-12. * (1 + rC) * rt * (-1 + w)) / mSqq;
        result.element({6, 2, 10, 9}) = (-4. * rt * (-1 + w) * (-3 + rC * (-1 + 2. * w))) / mSqq;
        result.element({6, 2, 10, 10}) = (-4. * rt * (1 + 3. * rC - 2. * w) * (-1 + w)) / mSqq;
        result.element({6, 2, 10, 11}) = (4. * (1 + rC) * rt * (-1 + wSq)) / mSqq;
        result.element({6, 3, 2, 0}) = (rt * (1 + w) * (-1 + rC * w)) / mSqq;
        result.element({6, 3, 4, 0}) = (rt * (1 + w) * (-1 + rC * w)) / mSqq;
        result.element({6, 3, 6, 2}) =
            ((1 + w) * (2. * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2. * rC * (2 + w) + rCSq * (-1 + 4. * w)))) /
            (3. * mSqq2);
        result.element({6, 3, 6, 3}) =
            ((1 + w) * (2. * mSqq * rCSq * (-1 + wSq) + rtSq * (3 - 6. * rC * w + rCSq * (-1 + 4. * wSq)))) /
            (3. * mSqq2);
        result.element({6, 3, 6, 4}) =
            ((1 + w) * (2. * mSqq * rC * (-1 + wSq) + rtSq * (3. * w + 3. * rCSq * w - 2. * rC * (2 + wSq)))) /
            (3. * mSqq2);
        result.element({6, 3, 8, 2}) =
            ((1 + w) * (2. * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2. * rC * (2 + w) + rCSq * (-1 + 4. * w)))) /
            (3. * mSqq2);
        result.element({6, 3, 8, 3}) =
            ((1 + w) * (2. * mSqq * rCSq * (-1 + wSq) + rtSq * (3 - 6. * rC * w + rCSq * (-1 + 4. * wSq)))) /
            (3. * mSqq2);
        result.element({6, 3, 8, 4}) =
            ((1 + w) * (2. * mSqq * rC * (-1 + wSq) + rtSq * (3. * w + 3. * rCSq * w - 2. * rC * (2 + wSq)))) /
            (3. * mSqq2);
        result.element({6, 3, 10, 8}) = (-4. * rC * rt * (-1 + wSq)) / mSqq;
        result.element({6, 3, 10, 9}) = (4. * rC * rt * (-1 + wSq)) / mSqq;
        result.element({6, 3, 10, 10}) = (-4. * rC * rt * (-1 + wSq)) / mSqq;
        result.element({6, 3, 10, 11}) = (4. * rC * rt * (-1 + w) * pow(1 + w, 2.)) / mSqq;
        result.element({6, 4, 2, 0}) = (rt * (rC - w) * (1 + w)) / mSqq;
        result.element({6, 4, 4, 0}) = (rt * (rC - w) * (1 + w)) / mSqq;
        result.element({6, 4, 6, 2}) =
            ((1 + w) * (2. * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3. * rCSq + 4. * w - 2. * rC * (2 + w)))) /
            (3. * mSqq2);
        result.element({6, 4, 6, 3}) =
            ((1 + w) * (2. * mSqq * rC * (-1 + wSq) + rtSq * (3. * w + 3. * rCSq * w - 2. * rC * (2 + wSq)))) /
            (3. * mSqq2);
        result.element({6, 4, 6, 4}) =
            ((1 + w) * (2. * mSqq * (-1 + wSq) + rtSq * (-1 + 3. * rCSq - 6. * rC * w + 4. * wSq))) / (3. * mSqq2);
        result.element({6, 4, 8, 2}) =
            ((1 + w) * (2. * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3. * rCSq + 4. * w - 2. * rC * (2 + w)))) /
            (3. * mSqq2);
        result.element({6, 4, 8, 3}) =
            ((1 + w) * (2. * mSqq * rC * (-1 + wSq) + rtSq * (3. * w + 3. * rCSq * w - 2. * rC * (2 + wSq)))) /
            (3. * mSqq2);
        result.element({6, 4, 8, 4}) =
            ((1 + w) * (2. * mSqq * (-1 + wSq) + rtSq * (-1 + 3. * rCSq - 6. * rC * w + 4. * wSq))) / (3. * mSqq2);
        result.element({6, 4, 10, 8}) = (-4. * rt * (-1 + wSq)) / mSqq;
        result.element({6, 4, 10, 9}) = (4. * rt * (-1 + wSq)) / mSqq;
        result.element({6, 4, 10, 10}) = (-4. * rt * (-1 + wSq)) / mSqq;
        result.element({6, 4, 10, 11}) = (4. * rt * (-1 + w) * pow(1 + w, 2.)) / mSqq;
        result.element({6, 5, 2, 1}) = -(((1 + rC) * rt * (-1 + w)) / mSqq);
        result.element({6, 5, 4, 1}) = ((1 + rC) * rt * (-1 + w)) / mSqq;
        result.element({6, 5, 6, 5}) = (2. * (2. * mSqq2 * (1 + w) + mSqq * (1 - 2. * rC + rCSq + rtSq) * (1 + w) +
                                              rtSq * (-1 + 2. * rC * (-2 + w) + 2. * w + rCSq * (-1 + 2. * w)))) /
                                       (3. * mSqq2);
        result.element({6, 5, 6, 6}) = ((-1 + w) * (2. * mSqq * (-1 + rC) * rC * (1 + w) +
                                                    rtSq * (-3 + 2. * rC * (-2 + w) + rCSq * (1 + 4. * w)))) /
                                       (3. * mSqq2);
        result.element({6, 5, 6, 7}) =
            ((-1 + w) * (rtSq * (-1 + 3. * rCSq - 2. * rC * (-2 + w) - 4. * w) + 2. * mSqq * (-1 + rC) * (1 + w))) /
            (3. * mSqq2);
        result.element({6, 5, 8, 5}) = (-2. * (2. * mSqq2 * (1 + w) + mSqq * (1 - 2. * rC + rCSq + rtSq) * (1 + w) +
                                               rtSq * (-1 + 2. * rC * (-2 + w) + 2. * w + rCSq * (-1 + 2. * w)))) /
                                       (3. * mSqq2);
        result.element({6, 5, 8, 6}) = -((-1 + w) * (2. * mSqq * (-1 + rC) * rC * (1 + w) +
                                                     rtSq * (-3 + 2. * rC * (-2 + w) + rCSq * (1 + 4. * w)))) /
                                       (3. * mSqq2);
        result.element({6, 5, 8, 7}) =
            -((-1 + w) * (rtSq * (-1 + 3. * rCSq - 2. * rC * (-2 + w) - 4. * w) + 2. * mSqq * (-1 + rC) * (1 + w))) /
            (3. * mSqq2);
        result.element({6, 5, 10, 8}) = (-12. * (-1 + rC) * rt * (1 + w)) / mSqq;
        result.element({6, 5, 10, 9}) = (-8. * rC * rt * (-1 + wSq)) / mSqq;
        result.element({6, 5, 10, 10}) = (-8. * rt * (-1 + wSq)) / mSqq;
        result.element({6, 6, 2, 1}) = -((rt * (-1 + w) * (-1 + rC * w)) / mSqq);
        result.element({6, 6, 4, 1}) = (rt * (-1 + w) * (-1 + rC * w)) / mSqq;
        result.element({6, 6, 6, 5}) = ((-1 + w) * (2. * mSqq * (-1 + rC) * rC * (1 + w) +
                                                    rtSq * (-3 + 2. * rC * (-2 + w) + rCSq * (1 + 4. * w)))) /
                                       (3. * mSqq2);
        result.element({6, 6, 6, 6}) =
            ((-1 + w) * (2. * mSqq * rCSq * (-1 + wSq) + rtSq * (3 - 6. * rC * w + rCSq * (-1 + 4. * wSq)))) /
            (3. * mSqq2);
        result.element({6, 6, 6, 7}) =
            ((-1 + w) * (2. * mSqq * rC * (-1 + wSq) + rtSq * (3. * w + 3. * rCSq * w - 2. * rC * (2 + wSq)))) /
            (3. * mSqq2);
        result.element({6, 6, 8, 5}) = -((-1 + w) * (2. * mSqq * (-1 + rC) * rC * (1 + w) +
                                                     rtSq * (-3 + 2. * rC * (-2 + w) + rCSq * (1 + 4. * w)))) /
                                       (3. * mSqq2);
        result.element({6, 6, 8, 6}) =
            -((-1 + w) * (2. * mSqq * rCSq * (-1 + wSq) + rtSq * (3 - 6. * rC * w + rCSq * (-1 + 4. * wSq)))) /
            (3. * mSqq2);
        result.element({6, 6, 8, 7}) =
            -((-1 + w) * (2. * mSqq * rC * (-1 + wSq) + rtSq * (3. * w + 3. * rCSq * w - 2. * rC * (2 + wSq)))) /
            (3. * mSqq2);
        result.element({6, 6, 10, 8}) = (-4. * rC * rt * (-1 + wSq)) / mSqq;
        result.element({6, 7, 2, 1}) = -((rt * (rC - w) * (-1 + w)) / mSqq);
        result.element({6, 7, 4, 1}) = (rt * (rC - w) * (-1 + w)) / mSqq;
        result.element({6, 7, 6, 5}) =
            ((-1 + w) * (rtSq * (-1 + 3. * rCSq - 2. * rC * (-2 + w) - 4. * w) + 2. * mSqq * (-1 + rC) * (1 + w))) /
            (3. * mSqq2);
        result.element({6, 7, 6, 6}) =
            ((-1 + w) * (2. * mSqq * rC * (-1 + wSq) + rtSq * (3. * w + 3. * rCSq * w - 2. * rC * (2 + wSq)))) /
            (3. * mSqq2);
        result.element({6, 7, 6, 7}) =
            ((-1 + w) * (2. * mSqq * (-1 + wSq) + rtSq * (-1 + 3. * rCSq - 6. * rC * w + 4. * wSq))) / (3. * mSqq2);
        result.element({6, 7, 8, 5}) =
            -((-1 + w) * (rtSq * (-1 + 3. * rCSq - 2. * rC * (-2 + w) - 4. * w) + 2. * mSqq * (-1 + rC) * (1 + w))) /
            (3. * mSqq2);
        result.element({6, 7, 8, 6}) =
            -((-1 + w) * (2. * mSqq * rC * (-1 + wSq) + rtSq * (3. * w + 3. * rCSq * w - 2. * rC * (2 + wSq)))) /
            (3. * mSqq2);
        result.element({6, 7, 8, 7}) =
            -((-1 + w) * (2. * mSqq * (-1 + wSq) + rtSq * (-1 + 3. * rCSq - 6. * rC * w + 4. * wSq))) / (3. * mSqq2);
        result.element({6, 7, 10, 8}) = (-4. * rt * (-1 + wSq)) / mSqq;
        result.element({7, 2, 0, 2}) = (2. * (2. * mSqq2 * (-1 + w) + mSqq * (1 + 2. * rC + rCSq + rtSq) * (-1 + w) +
                                              rtSq * (1 + 2. * w - 2. * rC * (2 + w) + rCSq * (1 + 2. * w)))) /
                                       (3. * mSqq2);
        result.element({7, 2, 0, 3}) =
            ((1 + w) * (2. * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2. * rC * (2 + w) + rCSq * (-1 + 4. * w)))) /
            (3. * mSqq2);
        result.element({7, 2, 0, 4}) =
            ((1 + w) * (2. * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3. * rCSq + 4. * w - 2. * rC * (2 + w)))) /
            (3. * mSqq2);
        result.element({7, 2, 1, 0}) = ((-1 + rC) * rt * (1 + w)) / mSqq;
        result.element({7, 2, 3, 0}) = ((-1 + rC) * rt * (1 + w)) / mSqq;
        result.element({7, 2, 5, 2}) = (2. * (2. * mSqq2 * (-1 + w) + mSqq * (1 + 2. * rC + rCSq + rtSq) * (-1 + w) +
                                              rtSq * (1 + 2. * w - 2. * rC * (2 + w) + rCSq * (1 + 2. * w)))) /
                                       (3. * mSqq2);
        result.element({7, 2, 5, 3}) =
            ((1 + w) * (2. * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2. * rC * (2 + w) + rCSq * (-1 + 4. * w)))) /
            (3. * mSqq2);
        result.element({7, 2, 5, 4}) =
            ((1 + w) * (2. * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3. * rCSq + 4. * w - 2. * rC * (2 + w)))) /
            (3. * mSqq2);
        result.element({7, 2, 7, 2}) = (2. * (2. * mSqq2 * (-1 + w) + mSqq * (1 + 2. * rC + rCSq + rtSq) * (-1 + w) +
                                              rtSq * (1 + 2. * w - 2. * rC * (2 + w) + rCSq * (1 + 2. * w)))) /
                                       (3. * mSqq2);
        result.element({7, 2, 7, 3}) =
            ((1 + w) * (2. * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2. * rC * (2 + w) + rCSq * (-1 + 4. * w)))) /
            (3. * mSqq2);
        result.element({7, 2, 7, 4}) =
            ((1 + w) * (2. * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3. * rCSq + 4. * w - 2. * rC * (2 + w)))) /
            (3. * mSqq2);
        result.element({7, 2, 9, 8}) = (-12. * (1 + rC) * rt * (-1 + w)) / mSqq;
        result.element({7, 2, 9, 9}) = (-4. * rt * (-1 + w) * (-3 + rC * (-1 + 2. * w))) / mSqq;
        result.element({7, 2, 9, 10}) = (-4. * rt * (1 + 3. * rC - 2. * w) * (-1 + w)) / mSqq;
        result.element({7, 2, 9, 11}) = (4. * (1 + rC) * rt * (-1 + wSq)) / mSqq;
        result.element({7, 3, 0, 2}) =
            ((1 + w) * (2. * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2. * rC * (2 + w) + rCSq * (-1 + 4. * w)))) /
            (3. * mSqq2);
        result.element({7, 3, 0, 3}) =
            ((1 + w) * (2. * mSqq * rCSq * (-1 + wSq) + rtSq * (3 - 6. * rC * w + rCSq * (-1 + 4. * wSq)))) /
            (3. * mSqq2);
        result.element({7, 3, 0, 4}) =
            ((1 + w) * (2. * mSqq * rC * (-1 + wSq) + rtSq * (3. * w + 3. * rCSq * w - 2. * rC * (2 + wSq)))) /
            (3. * mSqq2);
        result.element({7, 3, 1, 0}) = (rt * (1 + w) * (-1 + rC * w)) / mSqq;
        result.element({7, 3, 3, 0}) = (rt * (1 + w) * (-1 + rC * w)) / mSqq;
        result.element({7, 3, 5, 2}) =
            ((1 + w) * (2. * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2. * rC * (2 + w) + rCSq * (-1 + 4. * w)))) /
            (3. * mSqq2);
        result.element({7, 3, 5, 3}) =
            ((1 + w) * (2. * mSqq * rCSq * (-1 + wSq) + rtSq * (3 - 6. * rC * w + rCSq * (-1 + 4. * wSq)))) /
            (3. * mSqq2);
        result.element({7, 3, 5, 4}) =
            ((1 + w) * (2. * mSqq * rC * (-1 + wSq) + rtSq * (3. * w + 3. * rCSq * w - 2. * rC * (2 + wSq)))) /
            (3. * mSqq2);
        result.element({7, 3, 7, 2}) =
            ((1 + w) * (2. * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2. * rC * (2 + w) + rCSq * (-1 + 4. * w)))) /
            (3. * mSqq2);
        result.element({7, 3, 7, 3}) =
            ((1 + w) * (2. * mSqq * rCSq * (-1 + wSq) + rtSq * (3 - 6. * rC * w + rCSq * (-1 + 4. * wSq)))) /
            (3. * mSqq2);
        result.element({7, 3, 7, 4}) =
            ((1 + w) * (2. * mSqq * rC * (-1 + wSq) + rtSq * (3. * w + 3. * rCSq * w - 2. * rC * (2 + wSq)))) /
            (3. * mSqq2);
        result.element({7, 3, 9, 8}) = (-4. * rC * rt * (-1 + wSq)) / mSqq;
        result.element({7, 3, 9, 9}) = (4. * rC * rt * (-1 + wSq)) / mSqq;
        result.element({7, 3, 9, 10}) = (-4. * rC * rt * (-1 + wSq)) / mSqq;
        result.element({7, 3, 9, 11}) = (4. * rC * rt * (-1 + w) * pow(1 + w, 2.)) / mSqq;
        result.element({7, 4, 0, 2}) =
            ((1 + w) * (2. * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3. * rCSq + 4. * w - 2. * rC * (2 + w)))) /
            (3. * mSqq2);
        result.element({7, 4, 0, 3}) =
            ((1 + w) * (2. * mSqq * rC * (-1 + wSq) + rtSq * (3. * w + 3. * rCSq * w - 2. * rC * (2 + wSq)))) /
            (3. * mSqq2);
        result.element({7, 4, 0, 4}) =
            ((1 + w) * (2. * mSqq * (-1 + wSq) + rtSq * (-1 + 3. * rCSq - 6. * rC * w + 4. * wSq))) / (3. * mSqq2);
        result.element({7, 4, 1, 0}) = (rt * (rC - w) * (1 + w)) / mSqq;
        result.element({7, 4, 3, 0}) = (rt * (rC - w) * (1 + w)) / mSqq;
        result.element({7, 4, 5, 2}) =
            ((1 + w) * (2. * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3. * rCSq + 4. * w - 2. * rC * (2 + w)))) /
            (3. * mSqq2);
        result.element({7, 4, 5, 3}) =
            ((1 + w) * (2. * mSqq * rC * (-1 + wSq) + rtSq * (3. * w + 3. * rCSq * w - 2. * rC * (2 + wSq)))) /
            (3. * mSqq2);
        result.element({7, 4, 5, 4}) =
            ((1 + w) * (2. * mSqq * (-1 + wSq) + rtSq * (-1 + 3. * rCSq - 6. * rC * w + 4. * wSq))) / (3. * mSqq2);
        result.element({7, 4, 7, 2}) =
            ((1 + w) * (2. * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3. * rCSq + 4. * w - 2. * rC * (2 + w)))) /
            (3. * mSqq2);
        result.element({7, 4, 7, 3}) =
            ((1 + w) * (2. * mSqq * rC * (-1 + wSq) + rtSq * (3. * w + 3. * rCSq * w - 2. * rC * (2 + wSq)))) /
            (3. * mSqq2);
        result.element({7, 4, 7, 4}) =
            ((1 + w) * (2. * mSqq * (-1 + wSq) + rtSq * (-1 + 3. * rCSq - 6. * rC * w + 4. * wSq))) / (3. * mSqq2);
        result.element({7, 4, 9, 8}) = (-4. * rt * (-1 + wSq)) / mSqq;
        result.element({7, 4, 9, 9}) = (4. * rt * (-1 + wSq)) / mSqq;
        result.element({7, 4, 9, 10}) = (-4. * rt * (-1 + wSq)) / mSqq;
        result.element({7, 4, 9, 11}) = (4. * rt * (-1 + w) * pow(1 + w, 2.)) / mSqq;
        result.element({7, 5, 0, 5}) = (2. * (2. * mSqq2 * (1 + w) + mSqq * (1 - 2. * rC + rCSq + rtSq) * (1 + w) +
                                              rtSq * (-1 + 2. * rC * (-2 + w) + 2. * w + rCSq * (-1 + 2. * w)))) /
                                       (3. * mSqq2);
        result.element({7, 5, 0, 6}) = ((-1 + w) * (2. * mSqq * (-1 + rC) * rC * (1 + w) +
                                                    rtSq * (-3 + 2. * rC * (-2 + w) + rCSq * (1 + 4. * w)))) /
                                       (3. * mSqq2);
        result.element({7, 5, 0, 7}) =
            ((-1 + w) * (rtSq * (-1 + 3. * rCSq - 2. * rC * (-2 + w) - 4. * w) + 2. * mSqq * (-1 + rC) * (1 + w))) /
            (3. * mSqq2);
        result.element({7, 5, 1, 1}) = ((1 + rC) * rt * (-1 + w)) / mSqq;
        result.element({7, 5, 3, 1}) = -(((1 + rC) * rt * (-1 + w)) / mSqq);
        result.element({7, 5, 5, 5}) = (-2. * (2. * mSqq2 * (1 + w) + mSqq * (1 - 2. * rC + rCSq + rtSq) * (1 + w) +
                                               rtSq * (-1 + 2. * rC * (-2 + w) + 2. * w + rCSq * (-1 + 2. * w)))) /
                                       (3. * mSqq2);
        result.element({7, 5, 5, 6}) = -((-1 + w) * (2. * mSqq * (-1 + rC) * rC * (1 + w) +
                                                     rtSq * (-3 + 2. * rC * (-2 + w) + rCSq * (1 + 4. * w)))) /
                                       (3. * mSqq2);
        result.element({7, 5, 5, 7}) =
            -((-1 + w) * (rtSq * (-1 + 3. * rCSq - 2. * rC * (-2 + w) - 4. * w) + 2. * mSqq * (-1 + rC) * (1 + w))) /
            (3. * mSqq2);
        result.element({7, 5, 7, 5}) = (2. * (2. * mSqq2 * (1 + w) + mSqq * (1 - 2. * rC + rCSq + rtSq) * (1 + w) +
                                              rtSq * (-1 + 2. * rC * (-2 + w) + 2. * w + rCSq * (-1 + 2. * w)))) /
                                       (3. * mSqq2);
        result.element({7, 5, 7, 6}) = ((-1 + w) * (2. * mSqq * (-1 + rC) * rC * (1 + w) +
                                                    rtSq * (-3 + 2. * rC * (-2 + w) + rCSq * (1 + 4. * w)))) /
                                       (3. * mSqq2);
        result.element({7, 5, 7, 7}) =
            ((-1 + w) * (rtSq * (-1 + 3. * rCSq - 2. * rC * (-2 + w) - 4. * w) + 2. * mSqq * (-1 + rC) * (1 + w))) /
            (3. * mSqq2);
        result.element({7, 5, 9, 8}) = (-12. * (-1 + rC) * rt * (1 + w)) / mSqq;
        result.element({7, 5, 9, 9}) = (-8. * rC * rt * (-1 + wSq)) / mSqq;
        result.element({7, 5, 9, 10}) = (-8. * rt * (-1 + wSq)) / mSqq;
        result.element({7, 6, 0, 5}) = ((-1 + w) * (2. * mSqq * (-1 + rC) * rC * (1 + w) +
                                                    rtSq * (-3 + 2. * rC * (-2 + w) + rCSq * (1 + 4. * w)))) /
                                       (3. * mSqq2);
        result.element({7, 6, 0, 6}) =
            ((-1 + w) * (2. * mSqq * rCSq * (-1 + wSq) + rtSq * (3 - 6. * rC * w + rCSq * (-1 + 4. * wSq)))) /
            (3. * mSqq2);
        result.element({7, 6, 0, 7}) =
            ((-1 + w) * (2. * mSqq * rC * (-1 + wSq) + rtSq * (3. * w + 3. * rCSq * w - 2. * rC * (2 + wSq)))) /
            (3. * mSqq2);
        result.element({7, 6, 1, 1}) = (rt * (-1 + w) * (-1 + rC * w)) / mSqq;
        result.element({7, 6, 3, 1}) = -((rt * (-1 + w) * (-1 + rC * w)) / mSqq);
        result.element({7, 6, 5, 5}) = -((-1 + w) * (2. * mSqq * (-1 + rC) * rC * (1 + w) +
                                                     rtSq * (-3 + 2. * rC * (-2 + w) + rCSq * (1 + 4. * w)))) /
                                       (3. * mSqq2);
        result.element({7, 6, 5, 6}) =
            -((-1 + w) * (2. * mSqq * rCSq * (-1 + wSq) + rtSq * (3 - 6. * rC * w + rCSq * (-1 + 4. * wSq)))) /
            (3. * mSqq2);
        result.element({7, 6, 5, 7}) =
            -((-1 + w) * (2. * mSqq * rC * (-1 + wSq) + rtSq * (3. * w + 3. * rCSq * w - 2. * rC * (2 + wSq)))) /
            (3. * mSqq2);
        result.element({7, 6, 7, 5}) = ((-1 + w) * (2. * mSqq * (-1 + rC) * rC * (1 + w) +
                                                    rtSq * (-3 + 2. * rC * (-2 + w) + rCSq * (1 + 4. * w)))) /
                                       (3. * mSqq2);
        result.element({7, 6, 7, 6}) =
            ((-1 + w) * (2. * mSqq * rCSq * (-1 + wSq) + rtSq * (3 - 6. * rC * w + rCSq * (-1 + 4. * wSq)))) /
            (3. * mSqq2);
        result.element({7, 6, 7, 7}) =
            ((-1 + w) * (2. * mSqq * rC * (-1 + wSq) + rtSq * (3. * w + 3. * rCSq * w - 2. * rC * (2 + wSq)))) /
            (3. * mSqq2);
        result.element({7, 6, 9, 8}) = (-4. * rC * rt * (-1 + wSq)) / mSqq;
        result.element({7, 7, 0, 5}) =
            ((-1 + w) * (rtSq * (-1 + 3. * rCSq - 2. * rC * (-2 + w) - 4. * w) + 2. * mSqq * (-1 + rC) * (1 + w))) /
            (3. * mSqq2);
        result.element({7, 7, 0, 6}) =
            ((-1 + w) * (2. * mSqq * rC * (-1 + wSq) + rtSq * (3. * w + 3. * rCSq * w - 2. * rC * (2 + wSq)))) /
            (3. * mSqq2);
        result.element({7, 7, 0, 7}) =
            ((-1 + w) * (2. * mSqq * (-1 + wSq) + rtSq * (-1 + 3. * rCSq - 6. * rC * w + 4. * wSq))) / (3. * mSqq2);
        result.element({7, 7, 1, 1}) = (rt * (rC - w) * (-1 + w)) / mSqq;
        result.element({7, 7, 3, 1}) = -((rt * (rC - w) * (-1 + w)) / mSqq);
        result.element({7, 7, 5, 5}) =
            -((-1 + w) * (rtSq * (-1 + 3. * rCSq - 2. * rC * (-2 + w) - 4. * w) + 2. * mSqq * (-1 + rC) * (1 + w))) /
            (3. * mSqq2);
        result.element({7, 7, 5, 6}) =
            -((-1 + w) * (2. * mSqq * rC * (-1 + wSq) + rtSq * (3. * w + 3. * rCSq * w - 2. * rC * (2 + wSq)))) /
            (3. * mSqq2);
        result.element({7, 7, 5, 7}) =
            -((-1 + w) * (2. * mSqq * (-1 + wSq) + rtSq * (-1 + 3. * rCSq - 6. * rC * w + 4. * wSq))) / (3. * mSqq2);
        result.element({7, 7, 7, 5}) =
            ((-1 + w) * (rtSq * (-1 + 3. * rCSq - 2. * rC * (-2 + w) - 4. * w) + 2. * mSqq * (-1 + rC) * (1 + w))) /
            (3. * mSqq2);
        result.element({7, 7, 7, 6}) =
            ((-1 + w) * (2. * mSqq * rC * (-1 + wSq) + rtSq * (3. * w + 3. * rCSq * w - 2. * rC * (2 + wSq)))) /
            (3. * mSqq2);
        result.element({7, 7, 7, 7}) =
            ((-1 + w) * (2. * mSqq * (-1 + wSq) + rtSq * (-1 + 3. * rCSq - 6. * rC * w + 4. * wSq))) / (3. * mSqq2);
        result.element({7, 7, 9, 8}) = (-4. * rt * (-1 + wSq)) / mSqq;
        result.element({8, 2, 2, 0}) = ((-1 + rC) * rt * (1 + w)) / mSqq;
        result.element({8, 2, 4, 0}) = ((-1 + rC) * rt * (1 + w)) / mSqq;
        result.element({8, 2, 6, 2}) = (2. * (2. * mSqq2 * (-1 + w) + mSqq * (1 + 2. * rC + rCSq + rtSq) * (-1 + w) +
                                              rtSq * (1 + 2. * w - 2. * rC * (2 + w) + rCSq * (1 + 2. * w)))) /
                                       (3. * mSqq2);
        result.element({8, 2, 6, 3}) =
            ((1 + w) * (2. * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2. * rC * (2 + w) + rCSq * (-1 + 4. * w)))) /
            (3. * mSqq2);
        result.element({8, 2, 6, 4}) =
            ((1 + w) * (2. * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3. * rCSq + 4. * w - 2. * rC * (2 + w)))) /
            (3. * mSqq2);
        result.element({8, 2, 8, 2}) = (2. * (2. * mSqq2 * (-1 + w) + mSqq * (1 + 2. * rC + rCSq + rtSq) * (-1 + w) +
                                              rtSq * (1 + 2. * w - 2. * rC * (2 + w) + rCSq * (1 + 2. * w)))) /
                                       (3. * mSqq2);
        result.element({8, 2, 8, 3}) =
            ((1 + w) * (2. * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2. * rC * (2 + w) + rCSq * (-1 + 4. * w)))) /
            (3. * mSqq2);
        result.element({8, 2, 8, 4}) =
            ((1 + w) * (2. * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3. * rCSq + 4. * w - 2. * rC * (2 + w)))) /
            (3. * mSqq2);
        result.element({8, 2, 10, 8}) = (-12. * (1 + rC) * rt * (-1 + w)) / mSqq;
        result.element({8, 2, 10, 9}) = (-4. * rt * (-1 + w) * (-3 + rC * (-1 + 2. * w))) / mSqq;
        result.element({8, 2, 10, 10}) = (-4. * rt * (1 + 3. * rC - 2. * w) * (-1 + w)) / mSqq;
        result.element({8, 2, 10, 11}) = (4. * (1 + rC) * rt * (-1 + wSq)) / mSqq;
        result.element({8, 3, 2, 0}) = (rt * (1 + w) * (-1 + rC * w)) / mSqq;
        result.element({8, 3, 4, 0}) = (rt * (1 + w) * (-1 + rC * w)) / mSqq;
        result.element({8, 3, 6, 2}) =
            ((1 + w) * (2. * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2. * rC * (2 + w) + rCSq * (-1 + 4. * w)))) /
            (3. * mSqq2);
        result.element({8, 3, 6, 3}) =
            ((1 + w) * (2. * mSqq * rCSq * (-1 + wSq) + rtSq * (3 - 6. * rC * w + rCSq * (-1 + 4. * wSq)))) /
            (3. * mSqq2);
        result.element({8, 3, 6, 4}) =
            ((1 + w) * (2. * mSqq * rC * (-1 + wSq) + rtSq * (3. * w + 3. * rCSq * w - 2. * rC * (2 + wSq)))) /
            (3. * mSqq2);
        result.element({8, 3, 8, 2}) =
            ((1 + w) * (2. * mSqq * rC * (1 + rC) * (-1 + w) + rtSq * (3 - 2. * rC * (2 + w) + rCSq * (-1 + 4. * w)))) /
            (3. * mSqq2);
        result.element({8, 3, 8, 3}) =
            ((1 + w) * (2. * mSqq * rCSq * (-1 + wSq) + rtSq * (3 - 6. * rC * w + rCSq * (-1 + 4. * wSq)))) /
            (3. * mSqq2);
        result.element({8, 3, 8, 4}) =
            ((1 + w) * (2. * mSqq * rC * (-1 + wSq) + rtSq * (3. * w + 3. * rCSq * w - 2. * rC * (2 + wSq)))) /
            (3. * mSqq2);
        result.element({8, 3, 10, 8}) = (-4. * rC * rt * (-1 + wSq)) / mSqq;
        result.element({8, 3, 10, 9}) = (4. * rC * rt * (-1 + wSq)) / mSqq;
        result.element({8, 3, 10, 10}) = (-4. * rC * rt * (-1 + wSq)) / mSqq;
        result.element({8, 3, 10, 11}) = (4. * rC * rt * (-1 + w) * pow(1 + w, 2.)) / mSqq;
        result.element({8, 4, 2, 0}) = (rt * (rC - w) * (1 + w)) / mSqq;
        result.element({8, 4, 4, 0}) = (rt * (rC - w) * (1 + w)) / mSqq;
        result.element({8, 4, 6, 2}) =
            ((1 + w) * (2. * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3. * rCSq + 4. * w - 2. * rC * (2 + w)))) /
            (3. * mSqq2);
        result.element({8, 4, 6, 3}) =
            ((1 + w) * (2. * mSqq * rC * (-1 + wSq) + rtSq * (3. * w + 3. * rCSq * w - 2. * rC * (2 + wSq)))) /
            (3. * mSqq2);
        result.element({8, 4, 6, 4}) =
            ((1 + w) * (2. * mSqq * (-1 + wSq) + rtSq * (-1 + 3. * rCSq - 6. * rC * w + 4. * wSq))) / (3. * mSqq2);
        result.element({8, 4, 8, 2}) =
            ((1 + w) * (2. * mSqq * (1 + rC) * (-1 + w) + rtSq * (-1 + 3. * rCSq + 4. * w - 2. * rC * (2 + w)))) /
            (3. * mSqq2);
        result.element({8, 4, 8, 3}) =
            ((1 + w) * (2. * mSqq * rC * (-1 + wSq) + rtSq * (3. * w + 3. * rCSq * w - 2. * rC * (2 + wSq)))) /
            (3. * mSqq2);
        result.element({8, 4, 8, 4}) =
            ((1 + w) * (2. * mSqq * (-1 + wSq) + rtSq * (-1 + 3. * rCSq - 6. * rC * w + 4. * wSq))) / (3. * mSqq2);
        result.element({8, 4, 10, 8}) = (-4. * rt * (-1 + wSq)) / mSqq;
        result.element({8, 4, 10, 9}) = (4. * rt * (-1 + wSq)) / mSqq;
        result.element({8, 4, 10, 10}) = (-4. * rt * (-1 + wSq)) / mSqq;
        result.element({8, 4, 10, 11}) = (4. * rt * (-1 + w) * pow(1 + w, 2.)) / mSqq;
        result.element({8, 5, 2, 1}) = ((1 + rC) * rt * (-1 + w)) / mSqq;
        result.element({8, 5, 4, 1}) = -(((1 + rC) * rt * (-1 + w)) / mSqq);
        result.element({8, 5, 6, 5}) = (-2. * (2. * mSqq2 * (1 + w) + mSqq * (1 - 2. * rC + rCSq + rtSq) * (1 + w) +
                                               rtSq * (-1 + 2. * rC * (-2 + w) + 2. * w + rCSq * (-1 + 2. * w)))) /
                                       (3. * mSqq2);
        result.element({8, 5, 6, 6}) = -((-1 + w) * (2. * mSqq * (-1 + rC) * rC * (1 + w) +
                                                     rtSq * (-3 + 2. * rC * (-2 + w) + rCSq * (1 + 4. * w)))) /
                                       (3. * mSqq2);
        result.element({8, 5, 6, 7}) =
            -((-1 + w) * (rtSq * (-1 + 3. * rCSq - 2. * rC * (-2 + w) - 4. * w) + 2. * mSqq * (-1 + rC) * (1 + w))) /
            (3. * mSqq2);
        result.element({8, 5, 8, 5}) = (2. * (2. * mSqq2 * (1 + w) + mSqq * (1 - 2. * rC + rCSq + rtSq) * (1 + w) +
                                              rtSq * (-1 + 2. * rC * (-2 + w) + 2. * w + rCSq * (-1 + 2. * w)))) /
                                       (3. * mSqq2);
        result.element({8, 5, 8, 6}) = ((-1 + w) * (2. * mSqq * (-1 + rC) * rC * (1 + w) +
                                                    rtSq * (-3 + 2. * rC * (-2 + w) + rCSq * (1 + 4. * w)))) /
                                       (3. * mSqq2);
        result.element({8, 5, 8, 7}) =
            ((-1 + w) * (rtSq * (-1 + 3. * rCSq - 2. * rC * (-2 + w) - 4. * w) + 2. * mSqq * (-1 + rC) * (1 + w))) /
            (3. * mSqq2);
        result.element({8, 5, 10, 8}) = (12. * (-1 + rC) * rt * (1 + w)) / mSqq;
        result.element({8, 5, 10, 9}) = (8. * rC * rt * (-1 + wSq)) / mSqq;
        result.element({8, 5, 10, 10}) = (8. * rt * (-1 + wSq)) / mSqq;
        result.element({8, 6, 2, 1}) = (rt * (-1 + w) * (-1 + rC * w)) / mSqq;
        result.element({8, 6, 4, 1}) = -((rt * (-1 + w) * (-1 + rC * w)) / mSqq);
        result.element({8, 6, 6, 5}) = -((-1 + w) * (2. * mSqq * (-1 + rC) * rC * (1 + w) +
                                                     rtSq * (-3 + 2. * rC * (-2 + w) + rCSq * (1 + 4. * w)))) /
                                       (3. * mSqq2);
        result.element({8, 6, 6, 6}) =
            -((-1 + w) * (2. * mSqq * rCSq * (-1 + wSq) + rtSq * (3 - 6. * rC * w + rCSq * (-1 + 4. * wSq)))) /
            (3. * mSqq2);
        result.element({8, 6, 6, 7}) =
            -((-1 + w) * (2. * mSqq * rC * (-1 + wSq) + rtSq * (3. * w + 3. * rCSq * w - 2. * rC * (2 + wSq)))) /
            (3. * mSqq2);
        result.element({8, 6, 8, 5}) = ((-1 + w) * (2. * mSqq * (-1 + rC) * rC * (1 + w) +
                                                    rtSq * (-3 + 2. * rC * (-2 + w) + rCSq * (1 + 4. * w)))) /
                                       (3. * mSqq2);
        result.element({8, 6, 8, 6}) =
            ((-1 + w) * (2. * mSqq * rCSq * (-1 + wSq) + rtSq * (3 - 6. * rC * w + rCSq * (-1 + 4. * wSq)))) /
            (3. * mSqq2);
        result.element({8, 6, 8, 7}) =
            ((-1 + w) * (2. * mSqq * rC * (-1 + wSq) + rtSq * (3. * w + 3. * rCSq * w - 2. * rC * (2 + wSq)))) /
            (3. * mSqq2);
        result.element({8, 6, 10, 8}) = (4. * rC * rt * (-1 + wSq)) / mSqq;
        result.element({8, 7, 2, 1}) = (rt * (rC - w) * (-1 + w)) / mSqq;
        result.element({8, 7, 4, 1}) = -((rt * (rC - w) * (-1 + w)) / mSqq);
        result.element({8, 7, 6, 5}) =
            -((-1 + w) * (rtSq * (-1 + 3. * rCSq - 2. * rC * (-2 + w) - 4. * w) + 2. * mSqq * (-1 + rC) * (1 + w))) /
            (3. * mSqq2);
        result.element({8, 7, 6, 6}) =
            -((-1 + w) * (2. * mSqq * rC * (-1 + wSq) + rtSq * (3. * w + 3. * rCSq * w - 2. * rC * (2 + wSq)))) /
            (3. * mSqq2);
        result.element({8, 7, 6, 7}) =
            -((-1 + w) * (2. * mSqq * (-1 + wSq) + rtSq * (-1 + 3. * rCSq - 6. * rC * w + 4. * wSq))) / (3. * mSqq2);
        result.element({8, 7, 8, 5}) =
            ((-1 + w) * (rtSq * (-1 + 3. * rCSq - 2. * rC * (-2 + w) - 4. * w) + 2. * mSqq * (-1 + rC) * (1 + w))) /
            (3. * mSqq2);
        result.element({8, 7, 8, 6}) =
            ((-1 + w) * (2. * mSqq * rC * (-1 + wSq) + rtSq * (3. * w + 3. * rCSq * w - 2. * rC * (2 + wSq)))) /
            (3. * mSqq2);
        result.element({8, 7, 8, 7}) =
            ((-1 + w) * (2. * mSqq * (-1 + wSq) + rtSq * (-1 + 3. * rCSq - 6. * rC * w + 4. * wSq))) / (3. * mSqq2);
        result.element({8, 7, 10, 8}) = (4. * rt * (-1 + wSq)) / mSqq;
        result.element({9, 8, 0, 2}) = (-12. * (1 + rC) * rt * (-1 + w)) / mSqq;
        result.element({9, 8, 0, 3}) = (-4. * rC * rt * (-1 + wSq)) / mSqq;
        result.element({9, 8, 0, 4}) = (-4. * rt * (-1 + wSq)) / mSqq;
        result.element({9, 8, 0, 5}) = (-12. * (-1 + rC) * rt * (1 + w)) / mSqq;
        result.element({9, 8, 0, 6}) = (-4. * rC * rt * (-1 + wSq)) / mSqq;
        result.element({9, 8, 0, 7}) = (-4. * rt * (-1 + wSq)) / mSqq;
        result.element({9, 8, 5, 2}) = (-12. * (1 + rC) * rt * (-1 + w)) / mSqq;
        result.element({9, 8, 5, 3}) = (-4. * rC * rt * (-1 + wSq)) / mSqq;
        result.element({9, 8, 5, 4}) = (-4. * rt * (-1 + wSq)) / mSqq;
        result.element({9, 8, 5, 5}) = (12. * (-1 + rC) * rt * (1 + w)) / mSqq;
        result.element({9, 8, 5, 6}) = (4. * rC * rt * (-1 + wSq)) / mSqq;
        result.element({9, 8, 5, 7}) = (4. * rt * (-1 + wSq)) / mSqq;
        result.element({9, 8, 7, 2}) = (-12. * (1 + rC) * rt * (-1 + w)) / mSqq;
        result.element({9, 8, 7, 3}) = (-4. * rC * rt * (-1 + wSq)) / mSqq;
        result.element({9, 8, 7, 4}) = (-4. * rt * (-1 + wSq)) / mSqq;
        result.element({9, 8, 7, 5}) = (-12. * (-1 + rC) * rt * (1 + w)) / mSqq;
        result.element({9, 8, 7, 6}) = (-4. * rC * rt * (-1 + wSq)) / mSqq;
        result.element({9, 8, 7, 7}) = (-4. * rt * (-1 + wSq)) / mSqq;
        result.element({9, 8, 9, 8}) =
            (32. * (mSqq + 2. * rtSq) * (-4. * rC + (2 + mSqq) * w + 2. * rCSq * w)) / (3. * mSqq2);
        result.element({9, 8, 9, 9}) =
            (-16. * (mSqq + 2. * rtSq) * (-1 + w) * (2 + mSqq + 4. * rC - 2. * rCSq * (1 + 2. * w))) / (3. * mSqq2);
        result.element({9, 8, 9, 10}) =
            (16. * (mSqq + 2. * rtSq) * (-2 + mSqq + 4. * rC + 2. * rCSq - 4. * w) * (-1 + w)) / (3. * mSqq2);
        result.element({9, 8, 9, 11}) = (-16. * (mSqq + 2. * rtSq) * (-1 + wSq)) / (3. * mSqq);
        result.element({9, 9, 0, 2}) = (-4. * rt * (-1 + w) * (-3 + rC * (-1 + 2. * w))) / mSqq;
        result.element({9, 9, 0, 3}) = (4. * rC * rt * (-1 + wSq)) / mSqq;
        result.element({9, 9, 0, 4}) = (4. * rt * (-1 + wSq)) / mSqq;
        result.element({9, 9, 0, 5}) = (-8. * rC * rt * (-1 + wSq)) / mSqq;
        result.element({9, 9, 5, 2}) = (-4. * rt * (-1 + w) * (-3 + rC * (-1 + 2. * w))) / mSqq;
        result.element({9, 9, 5, 3}) = (4. * rC * rt * (-1 + wSq)) / mSqq;
        result.element({9, 9, 5, 4}) = (4. * rt * (-1 + wSq)) / mSqq;
        result.element({9, 9, 5, 5}) = (8. * rC * rt * (-1 + wSq)) / mSqq;
        result.element({9, 9, 7, 2}) = (-4. * rt * (-1 + w) * (-3 + rC * (-1 + 2. * w))) / mSqq;
        result.element({9, 9, 7, 3}) = (4. * rC * rt * (-1 + wSq)) / mSqq;
        result.element({9, 9, 7, 4}) = (4. * rt * (-1 + wSq)) / mSqq;
        result.element({9, 9, 7, 5}) = (-8. * rC * rt * (-1 + wSq)) / mSqq;
        result.element({9, 9, 9, 8}) =
            (-16. * (mSqq + 2. * rtSq) * (-1 + w) * (2 + mSqq + 4. * rC - 2. * rCSq * (1 + 2. * w))) / (3. * mSqq2);
        result.element({9, 9, 9, 9}) =
            (16. * (mSqq + 2. * rtSq) * (-1 + w) * (2 + mSqq - 4. * rC * w + rCSq * (-2 + 4. * wSq))) / (3. * mSqq2);
        result.element({9, 9, 9, 10}) =
            (-16. * (mSqq + 2. * rtSq) * (-1 + w) * (mSqq - 2. * (-2. * rC + w + rCSq * w))) / (3. * mSqq2);
        result.element({9, 9, 9, 11}) = (16. * (mSqq + 2. * rtSq) * (-1 + wSq)) / (3. * mSqq);
        result.element({9, 10, 0, 2}) = (-4. * rt * (1 + 3. * rC - 2. * w) * (-1 + w)) / mSqq;
        result.element({9, 10, 0, 3}) = (-4. * rC * rt * (-1 + wSq)) / mSqq;
        result.element({9, 10, 0, 4}) = (-4. * rt * (-1 + wSq)) / mSqq;
        result.element({9, 10, 0, 5}) = (-8. * rt * (-1 + wSq)) / mSqq;
        result.element({9, 10, 5, 2}) = (-4. * rt * (1 + 3. * rC - 2. * w) * (-1 + w)) / mSqq;
        result.element({9, 10, 5, 3}) = (-4. * rC * rt * (-1 + wSq)) / mSqq;
        result.element({9, 10, 5, 4}) = (-4. * rt * (-1 + wSq)) / mSqq;
        result.element({9, 10, 5, 5}) = (8. * rt * (-1 + wSq)) / mSqq;
        result.element({9, 10, 7, 2}) = (-4. * rt * (1 + 3. * rC - 2. * w) * (-1 + w)) / mSqq;
        result.element({9, 10, 7, 3}) = (-4. * rC * rt * (-1 + wSq)) / mSqq;
        result.element({9, 10, 7, 4}) = (-4. * rt * (-1 + wSq)) / mSqq;
        result.element({9, 10, 7, 5}) = (-8. * rt * (-1 + wSq)) / mSqq;
        result.element({9, 10, 9, 8}) =
            (16. * (mSqq + 2. * rtSq) * (-2 + mSqq + 4. * rC + 2. * rCSq - 4. * w) * (-1 + w)) / (3. * mSqq2);
        result.element({9, 10, 9, 9}) =
            (-16. * (mSqq + 2. * rtSq) * (-1 + w) * (mSqq - 2. * (-2. * rC + w + rCSq * w))) / (3. * mSqq2);
        result.element({9, 10, 9, 10}) =
            (16. * (mSqq + 2. * rtSq) * (-1 + w) * (-2 + mSqq + 2. * rCSq - 4. * rC * w + 4. * wSq)) / (3. * mSqq2);
        result.element({9, 10, 9, 11}) = (-16. * (mSqq + 2. * rtSq) * (-1 + wSq)) / (3. * mSqq);
        result.element({9, 11, 0, 2}) = (4. * (1 + rC) * rt * (-1 + wSq)) / mSqq;
        result.element({9, 11, 0, 3}) = (4. * rC * rt * (-1 + w) * pow(1 + w, 2.)) / mSqq;
        result.element({9, 11, 0, 4}) = (4. * rt * (-1 + w) * pow(1 + w, 2.)) / mSqq;
        result.element({9, 11, 5, 2}) = (4. * (1 + rC) * rt * (-1 + wSq)) / mSqq;
        result.element({9, 11, 5, 3}) = (4. * rC * rt * (-1 + w) * pow(1 + w, 2.)) / mSqq;
        result.element({9, 11, 5, 4}) = (4. * rt * (-1 + w) * pow(1 + w, 2.)) / mSqq;
        result.element({9, 11, 7, 2}) = (4. * (1 + rC) * rt * (-1 + wSq)) / mSqq;
        result.element({9, 11, 7, 3}) = (4. * rC * rt * (-1 + w) * pow(1 + w, 2.)) / mSqq;
        result.element({9, 11, 7, 4}) = (4. * rt * (-1 + w) * pow(1 + w, 2.)) / mSqq;
        result.element({9, 11, 9, 8}) = (-16. * (mSqq + 2. * rtSq) * (-1 + wSq)) / (3. * mSqq);
        result.element({9, 11, 9, 9}) = (16. * (mSqq + 2. * rtSq) * (-1 + wSq)) / (3. * mSqq);
        result.element({9, 11, 9, 10}) = (-16. * (mSqq + 2. * rtSq) * (-1 + wSq)) / (3. * mSqq);
        result.element({9, 11, 9, 11}) = (16. * (mSqq + 2. * rtSq) * (-1 + w) * pow(1 + w, 2.)) / (3. * mSqq);
        result.element({10, 8, 6, 2}) = (-12. * (1 + rC) * rt * (-1 + w)) / mSqq;
        result.element({10, 8, 6, 3}) = (-4. * rC * rt * (-1 + wSq)) / mSqq;
        result.element({10, 8, 6, 4}) = (-4. * rt * (-1 + wSq)) / mSqq;
        result.element({10, 8, 6, 5}) = (-12. * (-1 + rC) * rt * (1 + w)) / mSqq;
        result.element({10, 8, 6, 6}) = (-4. * rC * rt * (-1 + wSq)) / mSqq;
        result.element({10, 8, 6, 7}) = (-4. * rt * (-1 + wSq)) / mSqq;
        result.element({10, 8, 8, 2}) = (-12. * (1 + rC) * rt * (-1 + w)) / mSqq;
        result.element({10, 8, 8, 3}) = (-4. * rC * rt * (-1 + wSq)) / mSqq;
        result.element({10, 8, 8, 4}) = (-4. * rt * (-1 + wSq)) / mSqq;
        result.element({10, 8, 8, 5}) = (12. * (-1 + rC) * rt * (1 + w)) / mSqq;
        result.element({10, 8, 8, 6}) = (4. * rC * rt * (-1 + wSq)) / mSqq;
        result.element({10, 8, 8, 7}) = (4. * rt * (-1 + wSq)) / mSqq;
        result.element({10, 8, 10, 8}) =
            (32. * (mSqq + 2. * rtSq) * (-4. * rC + (2 + mSqq) * w + 2. * rCSq * w)) / (3. * mSqq2);
        result.element({10, 8, 10, 9}) =
            (-16. * (mSqq + 2. * rtSq) * (-1 + w) * (2 + mSqq + 4. * rC - 2. * rCSq * (1 + 2. * w))) / (3. * mSqq2);
        result.element({10, 8, 10, 10}) =
            (16. * (mSqq + 2. * rtSq) * (-2 + mSqq + 4. * rC + 2. * rCSq - 4. * w) * (-1 + w)) / (3. * mSqq2);
        result.element({10, 8, 10, 11}) = (-16. * (mSqq + 2. * rtSq) * (-1 + wSq)) / (3. * mSqq);
        result.element({10, 9, 6, 2}) = (-4. * rt * (-1 + w) * (-3 + rC * (-1 + 2. * w))) / mSqq;
        result.element({10, 9, 6, 3}) = (4. * rC * rt * (-1 + wSq)) / mSqq;
        result.element({10, 9, 6, 4}) = (4. * rt * (-1 + wSq)) / mSqq;
        result.element({10, 9, 6, 5}) = (-8. * rC * rt * (-1 + wSq)) / mSqq;
        result.element({10, 9, 8, 2}) = (-4. * rt * (-1 + w) * (-3 + rC * (-1 + 2. * w))) / mSqq;
        result.element({10, 9, 8, 3}) = (4. * rC * rt * (-1 + wSq)) / mSqq;
        result.element({10, 9, 8, 4}) = (4. * rt * (-1 + wSq)) / mSqq;
        result.element({10, 9, 8, 5}) = (8. * rC * rt * (-1 + wSq)) / mSqq;
        result.element({10, 9, 10, 8}) =
            (-16. * (mSqq + 2. * rtSq) * (-1 + w) * (2 + mSqq + 4. * rC - 2. * rCSq * (1 + 2. * w))) / (3. * mSqq2);
        result.element({10, 9, 10, 9}) =
            (16. * (mSqq + 2. * rtSq) * (-1 + w) * (2 + mSqq - 4. * rC * w + rCSq * (-2 + 4. * wSq))) / (3. * mSqq2);
        result.element({10, 9, 10, 10}) =
            (-16. * (mSqq + 2. * rtSq) * (-1 + w) * (mSqq - 2. * (-2. * rC + w + rCSq * w))) / (3. * mSqq2);
        result.element({10, 9, 10, 11}) = (16. * (mSqq + 2. * rtSq) * (-1 + wSq)) / (3. * mSqq);
        result.element({10, 10, 6, 2}) = (-4. * rt * (1 + 3. * rC - 2. * w) * (-1 + w)) / mSqq;
        result.element({10, 10, 6, 3}) = (-4. * rC * rt * (-1 + wSq)) / mSqq;
        result.element({10, 10, 6, 4}) = (-4. * rt * (-1 + wSq)) / mSqq;
        result.element({10, 10, 6, 5}) = (-8. * rt * (-1 + wSq)) / mSqq;
        result.element({10, 10, 8, 2}) = (-4. * rt * (1 + 3. * rC - 2. * w) * (-1 + w)) / mSqq;
        result.element({10, 10, 8, 3}) = (-4. * rC * rt * (-1 + wSq)) / mSqq;
        result.element({10, 10, 8, 4}) = (-4. * rt * (-1 + wSq)) / mSqq;
        result.element({10, 10, 8, 5}) = (8. * rt * (-1 + wSq)) / mSqq;
        result.element({10, 10, 10, 8}) =
            (16. * (mSqq + 2. * rtSq) * (-2 + mSqq + 4. * rC + 2. * rCSq - 4. * w) * (-1 + w)) / (3. * mSqq2);
        result.element({10, 10, 10, 9}) =
            (-16. * (mSqq + 2. * rtSq) * (-1 + w) * (mSqq - 2. * (-2. * rC + w + rCSq * w))) / (3. * mSqq2);
        result.element({10, 10, 10, 10}) =
            (16. * (mSqq + 2. * rtSq) * (-1 + w) * (-2 + mSqq + 2. * rCSq - 4. * rC * w + 4. * wSq)) / (3. * mSqq2);
        result.element({10, 10, 10, 11}) = (-16. * (mSqq + 2. * rtSq) * (-1 + wSq)) / (3. * mSqq);
        result.element({10, 11, 6, 2}) = (4. * (1 + rC) * rt * (-1 + wSq)) / mSqq;
        result.element({10, 11, 6, 3}) = (4. * rC * rt * (-1 + w) * pow(1 + w, 2.)) / mSqq;
        result.element({10, 11, 6, 4}) = (4. * rt * (-1 + w) * pow(1 + w, 2.)) / mSqq;
        result.element({10, 11, 8, 2}) = (4. * (1 + rC) * rt * (-1 + wSq)) / mSqq;
        result.element({10, 11, 8, 3}) = (4. * rC * rt * (-1 + w) * pow(1 + w, 2.)) / mSqq;
        result.element({10, 11, 8, 4}) = (4. * rt * (-1 + w) * pow(1 + w, 2.)) / mSqq;
        result.element({10, 11, 10, 8}) = (-16. * (mSqq + 2. * rtSq) * (-1 + wSq)) / (3. * mSqq);
        result.element({10, 11, 10, 9}) = (16. * (mSqq + 2. * rtSq) * (-1 + wSq)) / (3. * mSqq);
        result.element({10, 11, 10, 10}) = (-16. * (mSqq + 2. * rtSq) * (-1 + wSq)) / (3. * mSqq);
        result.element({10, 11, 10, 11}) = (16. * (mSqq + 2. * rtSq) * (-1 + w) * pow(1 + w, 2.)) / (3. * mSqq);

        result *= RateNorm;

        return result;
    }

} // namespace Hammer
