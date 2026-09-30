///
/// @file  TestIntegrator.cc
/// @brief Tests for Integrator
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#define RATEBDLEPNU_FRIENDS FRIEND_TEST(IntegratorTest, BDLepNu);

#define RATEBDSTARLEPNU_FRIENDS FRIEND_TEST(IntegratorTest, BDstarLepNu);

#define RATETAU3PINU_FRIENDS FRIEND_TEST(IntegratorTest, Tau3PiNu);

#define FFBTODBLPR_FRIENDS FRIEND_TEST(IntegratorTest, BDLepNu);

#define FFBTODSTARBLPR_FRIENDS FRIEND_TEST(IntegratorTest, BDstarLepNu);

#define FFTAUTO3PIRCT_FRIENDS FRIEND_TEST(IntegratorTest, Tau3PiNu);

#include "Hammer/Math/Integrator.hh"
#include "Hammer/Rates/RateBDLepNu.hh"
#include "Hammer/Rates/RateBDstarLepNu.hh"
#include "Hammer/Rates/RateTau3PiNu.hh"
#include "Hammer/FormFactors/BLPR/FFBtoDBLPR.hh"
#include "Hammer/FormFactors/BLPR/FFBtoDstarBLPR.hh"
#include "Hammer/FormFactors/RCT/FFTauto3PiRCT.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Math/PhaseSpace.hh"
#include "Hammer/Math/MultiDim/VectorContainer.hh"
#include "Hammer/Particle.hh"
// May need other stuff here

#include <iostream>

using namespace std;
using namespace complex_literals;

namespace Hammer {

    namespace MD = MultiDimensional;


    TEST(IntegratorTest, intTensor) {


        //        Tensor integrandRate;
        //        integrandRate.setDimensions({2, 2, 19});
        //        integrandRate.setLabels({WILSON_BCTAUNU, WILSON_BCTAUNU_HC, INTEGRATION_INDEX});

        // Obtain integration points for range -1 to 2
        Integrator integT;
        integT.applyRange({[](const vector<double>&) -> pair<double, double> { return make_pair(-1., 2.); }});
        const EvaluationGrid& points = integT.getEvaluationPoints();
        const EvaluationWeights& weights = integT.getEvaluationWeights();

        Tensor integrandFF{"integrandFF",
                           MD::makeEmptyVector({static_cast<IndexType>(points.size())}, {INTEGRATION_INDEX})};
        for (IndexType i = 0; i < static_cast<IndexType>(points.size()); ++i) {
            integrandFF.element({i}) = 1.;
        }

        Tensor integrandRate{"integrandRate",
                             MD::makeEmptyVector({2, 2, static_cast<IndexType>(points.size())},
                                                 {WILSON_BCTAUNU, WILSON_BCTAUNU_HC, INTEGRATION_INDEX})};
        for (IndexType i = 0; i < static_cast<IndexType>(points.size()); ++i) {
            integrandRate.element({0, 0, i}) = sin(points[i][0]) * weights[i];
            integrandRate.element({0, 1, i}) = exp(points[i][0]) * weights[i];
            integrandRate.element({1, 0, i}) = cos(points[i][0]) * weights[i];
            integrandRate.element({1, 1, i}) = weights[i];
        }

        //        Tensor t;
        //        t.setDimensions({2, 2});
        //        t.setLabels({WILSON_BCTAUNU, WILSON_BCTAUNU_HC});
        //        for (size_t i = 0; i < points.size(); ++i) {
        //            t.setValues({sin(points[i]), exp(points[i]), cos(points[i]), 1.});
        //            cout << t.element({0,0}) << '\n';
        //            integrandRate.addAt(t, INTEGRATION_INDEX, i);
        //            cout << integrandRate.element({0,0,i}) << '\n';
        //        }

        // Integrate
        integrandRate.dot(integrandFF, {INTEGRATION_INDEX});

        // Explicit Values
        Tensor integTEval{"integTEval", MD::makeVector({2, 2}, {WILSON_BCTAUNU, WILSON_BCTAUNU_HC},
                                                       {0.9564491424, 7.021176658, 1.750768412, 3.})};

        for (IndexType idx1 = 0; idx1 < 2; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 2; ++idx2) {
                double comRe =
                    compareVals(integrandRate.element({idx1, idx2}).real(), integTEval.element({idx1, idx2}).real());
                EXPECT_NEAR(comRe, 1., 1e-5);
            }
        }
    }

    inline BoundaryFunction makem12Function(const vector<double>& masses) {
        return ([masses](const vector<double>&) -> pair<double, double> { // NOLINT(bugprone-exception-escape)
            return make_pair(masses[1] + masses[2], masses[0] - masses[3]);
        });
    }

    inline BoundaryFunction makem23Function(const vector<double>& masses) {
        return ([masses](const vector<double>& vals) -> pair<double, double> { // NOLINT(bugprone-exception-escape)
            double E2Star = (pow(vals[0], 2.) - pow(masses[1], 2.) + pow(masses[2], 2.)) / (2 * vals[0]);
            double E3Star = (pow(masses[0], 2.) - pow(vals[0], 2.) - pow(masses[3], 2.)) / (2 * vals[0]);
            double p2Star = sqrt(pow(E2Star, 2.) - pow(masses[2], 2.));
            double p3Star = sqrt(pow(E3Star, 2.) - pow(masses[3], 2.));
            double E23Max = sqrt(pow(E2Star + E3Star, 2.) - pow(p2Star - p3Star, 2.));
            double E23Min = sqrt(pow(E2Star + E3Star, 2.) - pow(p2Star + p3Star, 2.));
            return make_pair(E23Min, E23Max);
        });
    }

    inline function<double(const vector<double>&)> pStar(double m1) {
        return ([m1](const vector<double>& vals) -> double {
            double EStar = (pow(vals[0], 2.) - pow(vals[1], 2.) + pow(m1, 2.)) / (2. * vals[0]);
            return sqrt(pow(EStar, 2.) - pow(m1, 2.));
        });
    }

    TEST(IntegratorTest, DalitzSpace) {
        Integrator integT;

        vector<double> masses{1.2, 0.2, 0.3, 0.01};
        IntegrationBoundaries boundaries{makem12Function(masses)};
        integT.applyRange(boundaries);
        const EvaluationGrid& points = integT.getEvaluationPoints();
        const EvaluationWeights& weights = integT.getEvaluationWeights();

        //        auto range = makem12Function(masses)({});
        //        cout << range.first << " " << range.second << '\n';

        double res = 0;
        for (size_t i = 0; i < points.size(); ++i) {
            res += pStar(masses[3])({masses[0], points[i][0]}) * pStar(masses[2])({points[i][0], masses[1]}) /
                   masses[0] * weights[i];
        }
        res *= 1 / (16. * pi3);
        //        cout << "Int: " << res << '\n';

        Integrator integTD;
        IntegrationBoundaries boundariesD{makem12Function(masses), makem23Function(masses)};
        integTD.applyRange(boundariesD);
        const EvaluationGrid& pointsD = integTD.getEvaluationPoints();
        const EvaluationWeights& weightsD = integTD.getEvaluationWeights();

        double resD = 0;
        for (size_t i = 0; i < pointsD.size(); ++i) {
            //        cout << pointsD[i][0] << " " << pointsD[i][1] << "  " << weightsD[i] << '\n';
            resD += 4. * pointsD[i][0] * pointsD[i][1] / pow(masses[0], 2.) * weightsD[i];
        }
        resD *= 1. / (128. * pi3);
        //        cout << "IntD: " << resD << '\n';

        double com = compareVals(res, resD);
        EXPECT_NEAR(com, 1., 1e-4);
    }

    TEST(IntegratorTest, BDLepNu) {
        SettingsHandler sh;
        sh.addSetting<string>("Hammer", "Units", "GeV");

        // Setup Rate and FF integrand tensors
        RateBDLepNu rate;
        rate.setSettingsHandler(sh);
        // Set process to B^0 -> D^- e^+ nu
        rate.setSignatureIndex(5);
        rate.init();
        rate.calcTensor();
        const auto& range = rate.getEvaluationPoints();

        FFBtoDBLPR ffBLPR;
        ffBLPR.setSettingsHandler(sh);
        // Set process to B^0 -> D^*-
        ffBLPR.setSignatureIndex(1);
        ffBLPR.init();
        Tensor ffTen = ffBLPR.getFFPSIntegrand(range);

        // Integrate
        Tensor rateT = rate.getTensor();
        ffTen.dot(rateT, {INTEGRATION_INDEX, FF_BD, FF_BD_HC});

        Tensor intEval{"intEval", MD::makeVector({11, 11}, {WILSON_BCTAUNU, WILSON_BCTAUNU_HC},
                                                 {5.745905878e-12,
                                                  -1.930066489e-15,
                                                  0.,
                                                  -1.930066489e-15,
                                                  0.,
                                                  5.745905878e-12,
                                                  0.,
                                                  5.745905878e-12,
                                                  0.,
                                                  -2.832461217e-15,
                                                  0.,
                                                  -1.930066489e-15,
                                                  5.141348909e-12,
                                                  0.,
                                                  5.141348909e-12,
                                                  0.,
                                                  -1.930066489e-15,
                                                  0.,
                                                  -1.930066489e-15,
                                                  0.,
                                                  0.,
                                                  0.,
                                                  0.,
                                                  0.,
                                                  5.141348909e-12,
                                                  0.,
                                                  5.141348909e-12,
                                                  0.,
                                                  -1.930066489e-15,
                                                  0.,
                                                  -1.930066489e-15,
                                                  0.,
                                                  0.,
                                                  -1.930066489e-15,
                                                  5.141348909e-12,
                                                  0.,
                                                  5.141348909e-12,
                                                  0.,
                                                  -1.930066489e-15,
                                                  0.,
                                                  -1.930066489e-15,
                                                  0.,
                                                  0.,
                                                  0.,
                                                  0.,
                                                  0.,
                                                  5.141348909e-12,
                                                  0.,
                                                  5.141348909e-12,
                                                  0.,
                                                  -1.930066489e-15,
                                                  0.,
                                                  -1.930066489e-15,
                                                  0.,
                                                  0.,
                                                  5.745905878e-12,
                                                  -1.930066489e-15,
                                                  0.,
                                                  -1.930066489e-15,
                                                  0.,
                                                  5.745905878e-12,
                                                  0.,
                                                  5.745905878e-12,
                                                  0.,
                                                  -2.832461217e-15,
                                                  0.,
                                                  0.,
                                                  0.,
                                                  -1.930066489e-15,
                                                  0.,
                                                  -1.930066489e-15,
                                                  0.,
                                                  5.745905878e-12,
                                                  0.,
                                                  5.745905878e-12,
                                                  0.,
                                                  -2.832461217e-15,
                                                  5.745905878e-12,
                                                  -1.930066489e-15,
                                                  0.,
                                                  -1.930066489e-15,
                                                  0.,
                                                  5.745905878e-12,
                                                  0.,
                                                  5.745905878e-12,
                                                  0.,
                                                  -2.832461217e-15,
                                                  0.,
                                                  0.,
                                                  0.,
                                                  -1.930066489e-15,
                                                  0.,
                                                  -1.930066489e-15,
                                                  0.,
                                                  5.745905878e-12,
                                                  0.,
                                                  5.745905878e-12,
                                                  0.,
                                                  -2.832461217e-15,
                                                  -2.832461217e-15,
                                                  0.,
                                                  0.,
                                                  0.,
                                                  0.,
                                                  -2.832461217e-15,
                                                  0.,
                                                  -2.832461217e-15,
                                                  0.,
                                                  4.411459379e-12,
                                                  0.,
                                                  0.,
                                                  0.,
                                                  0.,
                                                  0.,
                                                  0.,
                                                  0.,
                                                  -2.832461217e-15,
                                                  0.,
                                                  -2.832461217e-15,
                                                  0.,
                                                  4.411459379e-12})};

        // Check TEs.
        for (IndexType idx1 = 0; idx1 < 11; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 11; ++idx2) {
                //                    double val = integT.element({idx1, idx2}).real();
                //                    //drop electron mass terms
                //                    if (abs(val) < 1.e-14){
                //                        val = 0.;
                //                        }
                //                    cout << idx1 << " " <<  idx2 << " " <<  integT.element({idx1, idx2}).real() << "
                //                    " << intEval.element({idx1, idx2}).real() << '\n';
                double comRe = compareVals(ffTen.element({idx1, idx2}).real(), intEval.element({idx1, idx2}).real());
                EXPECT_NEAR(comRe, 1., 8e-3);
            }
        }
    }

    //     TEST(IntegratorTest, BDMuVsE) {
    //
    //
    //        SettingsHandler sh;
    //        SettingsConsumer::setSettingsHandler(sh);
    //
    //        //Setup Rate and FF integrand tensors
    //        RateBDLepNu rateE;
    //        RateBDLepNu rateM;
    //        //Set process to B^0 -> D^- e^+ nu
    //        rateE.setSignatureIndex(5);
    //        //Set process to B^0 -> D^- mu^+ nu
    //        rateM.setSignatureIndex(3);
    //        rateE.init();
    //        rateM.init();
    //        auto rangeE = rateE.getQ2Range();
    //        auto rangeM = rateE.getQ2Range();
    //
    //        FFBtoDBLPR ffBLPR;
    //        //Set process to B^0 -> D^*-
    //        ffBLPR.setSignatureIndex(1);
    //        ffBLPR.init();
    //        Tensor ffTenE = ffBLPR.getFFQ2Integrand(rangeE);
    //        Tensor ffTenM = ffBLPR.getFFQ2Integrand(rangeM);
    //
    //        //Integrate
    //        Integrator integTE;
    //        Integrator integTM;
    //        integTE.applyNorm(rangeE.first, rangeE.second);
    //        integTM.applyNorm(rangeM.first, rangeM.second);
    //        integTE.dot(ffTenE);
    //        integTM.dot(ffTenM);
    //        Tensor rateTE = rateE.toTensor();
    //        Tensor rateTM = rateM.toTensor();
    //        integTE.dot(rateTE);
    //        integTM.dot(rateTM);
    //
    //        Tensor wcWpE = makeEmptyTensor("WCWp", {11}, {WILSON_BCENU});
    //        Tensor wcWpM = makeEmptyTensor("WCWp", {11}, {WILSON_BCMUNU});
    //        wcWpE.data() = {1.,0.,0.,0.,0.,0.,0.46,0.,0.,0.,0.};
    //        wcWpM.data() = {1.,0.,0.,0.,0.,0.,0.46,0.,0.,0.,0.};
    //
    //        integTE.dot(outerSquare(wcWpE));
    //        integTM.dot(outerSquare(wcWpM));
    //
    //        cout << "RatioEM: " << integTE.firstEntry().real()/integTM.firstEntry().real() << '\n';
    //
    //     }

    TEST(IntegratorTest, BDstarLepNu) {
        SettingsHandler sh;
        sh.addSetting<string>("Hammer", "Units", "GeV");

        // Setup Rate and FF integrand tensors
        RateBDstarLepNu rate;
        rate.setSettingsHandler(sh);
        // Set process to B^0 -> D^*- e^+ nu
        rate.setSignatureIndex(5);
        rate.init();
        rate.calcTensor();
        const auto& range = rate.getEvaluationPoints();

        FFBtoDstarBLPR ffBLPR;
        ffBLPR.setSettingsHandler(sh);
        // Set process to B^0 -> D^*-
        ffBLPR.setSignatureIndex(1);
        ffBLPR.init();
        Tensor ffTen = ffBLPR.getFFPSIntegrand(range);

        // Integrate
        Tensor rateT = rate.getTensor();
        ffTen.dot(rateT, {INTEGRATION_INDEX, FF_BDSTAR, FF_BDSTAR_HC});

        Tensor intEval{"intEval", MD::makeVector({11, 11}, {WILSON_BCTAUNU, WILSON_BCTAUNU_HC},
                                                 {1.489204009e-11,
                                                  6.614803128e-16,
                                                  0.,
                                                  -6.614803128e-16,
                                                  0.,
                                                  -1.308238231e-11,
                                                  0.,
                                                  1.489204009e-11,
                                                  0.,
                                                  1.683873656e-14,
                                                  0.,
                                                  6.614803128e-16,
                                                  7.704082029e-13,
                                                  0.,
                                                  -7.704082029e-13,
                                                  0.,
                                                  -6.614803128e-16,
                                                  0.,
                                                  6.614803128e-16,
                                                  0.,
                                                  0.,
                                                  0.,
                                                  0.,
                                                  0.,
                                                  7.704082029e-13,
                                                  0.,
                                                  -7.704082029e-13,
                                                  0.,
                                                  -6.614803128e-16,
                                                  0.,
                                                  6.614803128e-16,
                                                  0.,
                                                  0.,
                                                  -6.614803128e-16,
                                                  -7.704082029e-13,
                                                  0.,
                                                  7.704082029e-13,
                                                  0.,
                                                  6.614803128e-16,
                                                  0.,
                                                  -6.614803128e-16,
                                                  0.,
                                                  0.,
                                                  0.,
                                                  0.,
                                                  0.,
                                                  -7.704082029e-13,
                                                  0.,
                                                  7.704082029e-13,
                                                  0.,
                                                  6.614803128e-16,
                                                  0.,
                                                  -6.614803128e-16,
                                                  0.,
                                                  0.,
                                                  -1.308238231e-11,
                                                  -6.614803128e-16,
                                                  0.,
                                                  6.614803128e-16,
                                                  0.,
                                                  1.489204009e-11,
                                                  0.,
                                                  -1.308238231e-11,
                                                  0.,
                                                  -2.700944240e-14,
                                                  0.,
                                                  0.,
                                                  0.,
                                                  -6.614803128e-16,
                                                  0.,
                                                  6.614803128e-16,
                                                  0.,
                                                  1.489204009e-11,
                                                  0.,
                                                  -1.308238231e-11,
                                                  0.,
                                                  1.683873656e-14,
                                                  1.489204009e-11,
                                                  6.614803128e-16,
                                                  0.,
                                                  -6.614803128e-16,
                                                  0.,
                                                  -1.308238231e-11,
                                                  0.,
                                                  1.489204009e-11,
                                                  0.,
                                                  1.683873656e-14,
                                                  0.,
                                                  0.,
                                                  0.,
                                                  6.614803128e-16,
                                                  0.,
                                                  -6.614803128e-16,
                                                  0.,
                                                  -1.308238231e-11,
                                                  0.,
                                                  1.489204009e-11,
                                                  0.,
                                                  -2.700944240e-14,
                                                  1.683873656e-14,
                                                  0.,
                                                  0.,
                                                  0.,
                                                  0.,
                                                  -2.700944240e-14,
                                                  0.,
                                                  1.683873656e-14,
                                                  0.,
                                                  2.471859950e-10,
                                                  0.,
                                                  0.,
                                                  0.,
                                                  0.,
                                                  0.,
                                                  0.,
                                                  0.,
                                                  1.683873656e-14,
                                                  0.,
                                                  -2.700944240e-14,
                                                  0.,
                                                  2.471859950e-10})};

        // Check TEs.
        for (IndexType idx1 = 0; idx1 < 11; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 11; ++idx2) {
                //                    double val = integT.element({idx1, idx2}).real();
                //                    if (abs(val) < 1.e-13){
                //                        val = 0.;
                //                        }
                //                    cout << idx1 << " " <<  idx2 << " " <<integT.element({idx1, idx2}).real() << "  "
                //                    << intEval.element({idx1, idx2}).real() << '\n';
                double comRe = compareVals(ffTen.element({idx1, idx2}).real(), intEval.element({idx1, idx2}).real());
                EXPECT_NEAR(comRe, 1., 5.1e-3);
            }
        }
    }

    //    TEST(IntegratorTest, Tau3PiNu) {
    //
    //        SettingsHandler sh;
    //        SettingsConsumer::setSettingsHandler(sh);
    //
    //        //Setup Rate and FF integrand tensors
    //        RateTau3PiNu rate;
    //        rate.setSignatureIndex(0);
    //        rate.init();
    //        rate.calcTensor();
    //        const auto& range = rate.getEvaluationPoints();
    //
    //        FFTauto3PiRCT ffRCT;
    //        ffRCT.setSignatureIndex(0);
    //        ffRCT.init();
    //
    //        Tensor ffTen = ffRCT.getFFPSIntegrand(range);
    //
    //        //Integrate
    //        Tensor rateT = rate.toTensor();
    //
    //        auto ratedata = rateT.data();
    //        auto ffdata = ffTen.data();
    //        complex<double> res = 0.;
    //        for(size_t idx = 0; idx < ratedata.size(); ++idx){
    //            res += ffdata[idx]*ratedata[idx];
    //        }
    //        cout << "Int: " << res.real() << '\n';
    //
    ////        ffTen.dot(rateT);
    ////        cout << "Int: " << ffTen.firstEntry().real() << '\n';
    //
    //    }

} // namespace Hammer
