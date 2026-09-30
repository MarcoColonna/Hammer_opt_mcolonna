///
/// @file  RateBDstarLepNu.cc
/// @brief \f$ B \rightarrow D^* \tau\nu \f$ total rate
/// @brief \f$ B_s \rightarrow D_s^* \tau\nu \f$ total rate
/// @brief \f$ B_c \rightarrow J/\psi \tau\nu \f$ total rate
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/Rates/RateBDstarLepNu.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"
#include "Hammer/Tools/Pdg.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    RateBDstarLepNu::RateBDstarLepNu() {
        // Create tensor rank and dimensions
        IndexList dims{{11, 8, 11, 8, _nPoints}};
        string name{"RateBDstarLepNuQ2"};
        auto& pdg = PID::instance();
        addProcessSignature(PID::BPLUS, {-PID::DSTAR, PID::NU_TAU, PID::ANTITAU});
        addIntegrationBoundaries(
            {PS::makeQ2Function(pdg.getMass(PID::ANTITAU), pdg.getMass(PID::BPLUS) - pdg.getMass(PID::DSTAR))});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCTAUNU, FF_BDSTAR, WILSON_BCTAUNU_HC, FF_BDSTAR_HC,
                                                          INTEGRATION_INDEX})});

        addProcessSignature(PID::BZERO, {PID::DSTARMINUS, PID::NU_TAU, PID::ANTITAU});
        addIntegrationBoundaries(
            {PS::makeQ2Function(pdg.getMass(PID::ANTITAU), pdg.getMass(PID::BZERO) - pdg.getMass(PID::DSTARMINUS))});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCTAUNU, FF_BDSTAR, WILSON_BCTAUNU_HC, FF_BDSTAR_HC,
                                                          INTEGRATION_INDEX})});

        addProcessSignature(PID::BPLUS, {-PID::DSTAR, PID::NU_MU, PID::ANTIMUON});
        addIntegrationBoundaries(
            {PS::makeQ2Function(pdg.getMass(PID::ANTIMUON), pdg.getMass(PID::BPLUS) - pdg.getMass(PID::DSTAR))});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCMUNU, FF_BDSTAR, WILSON_BCMUNU_HC, FF_BDSTAR_HC,
                                                          INTEGRATION_INDEX})});

        addProcessSignature(PID::BZERO, {PID::DSTARMINUS, PID::NU_MU, PID::ANTIMUON});
        addIntegrationBoundaries(
            {PS::makeQ2Function(pdg.getMass(PID::ANTIMUON), pdg.getMass(PID::BZERO) - pdg.getMass(PID::DSTARMINUS))});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCMUNU, FF_BDSTAR, WILSON_BCMUNU_HC, FF_BDSTAR_HC,
                                                          INTEGRATION_INDEX})});

        addProcessSignature(PID::BPLUS, {-PID::DSTAR, PID::NU_E, PID::POSITRON});
        addIntegrationBoundaries(
            {PS::makeQ2Function(pdg.getMass(PID::POSITRON), pdg.getMass(PID::BPLUS) - pdg.getMass(PID::DSTAR))});
        addTensor(Tensor{name, MD::makeEmptySparse(
                                   dims, {WILSON_BCENU, FF_BDSTAR, WILSON_BCENU_HC, FF_BDSTAR_HC, INTEGRATION_INDEX})});

        addProcessSignature(PID::BZERO, {PID::DSTARMINUS, PID::NU_E, PID::POSITRON});
        addIntegrationBoundaries(
            {PS::makeQ2Function(pdg.getMass(PID::POSITRON), pdg.getMass(PID::BZERO) - pdg.getMass(PID::DSTARMINUS))});
        addTensor(Tensor{name, MD::makeEmptySparse(
                                   dims, {WILSON_BCENU, FF_BDSTAR, WILSON_BCENU_HC, FF_BDSTAR_HC, INTEGRATION_INDEX})});

        // bs -> cs
        addProcessSignature(PID::BS, {PID::DSSTARMINUS, PID::NU_TAU, PID::ANTITAU});
        addIntegrationBoundaries(
            {PS::makeQ2Function(pdg.getMass(PID::ANTITAU), pdg.getMass(PID::BS) - pdg.getMass(PID::DSSTARMINUS))});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCTAUNU, FF_BSDSSTAR, WILSON_BCTAUNU_HC,
                                                          FF_BSDSSTAR_HC, INTEGRATION_INDEX})});

        addProcessSignature(PID::BS, {PID::DSSTARMINUS, PID::NU_MU, PID::ANTIMUON});
        addIntegrationBoundaries(
            {PS::makeQ2Function(pdg.getMass(PID::ANTIMUON), pdg.getMass(PID::BS) - pdg.getMass(PID::DSSTARMINUS))});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCMUNU, FF_BSDSSTAR, WILSON_BCMUNU_HC, FF_BSDSSTAR_HC,
                                                          INTEGRATION_INDEX})});

        addProcessSignature(PID::BS, {PID::DSSTARMINUS, PID::NU_E, PID::POSITRON});
        addIntegrationBoundaries(
            {PS::makeQ2Function(pdg.getMass(PID::POSITRON), pdg.getMass(PID::BS) - pdg.getMass(PID::DSSTARMINUS))});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCENU, FF_BSDSSTAR, WILSON_BCENU_HC, FF_BSDSSTAR_HC,
                                                          INTEGRATION_INDEX})});

        // bc -> cc
        addProcessSignature(PID::BCPLUS, {PID::JPSI, PID::NU_TAU, PID::ANTITAU});
        addIntegrationBoundaries(
            {PS::makeQ2Function(pdg.getMass(PID::ANTITAU), pdg.getMass(PID::BCPLUS) - pdg.getMass(PID::PID::JPSI))});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCTAUNU, FF_BCJPSI, WILSON_BCTAUNU_HC, FF_BCJPSI_HC,
                                                          INTEGRATION_INDEX})});

        addProcessSignature(PID::BCPLUS, {PID::JPSI, PID::NU_MU, PID::ANTIMUON});
        addIntegrationBoundaries(
            {PS::makeQ2Function(pdg.getMass(PID::ANTIMUON), pdg.getMass(PID::BCPLUS) - pdg.getMass(PID::PID::JPSI))});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BCMUNU, FF_BCJPSI, WILSON_BCMUNU_HC, FF_BCJPSI_HC,
                                                          INTEGRATION_INDEX})});

        addProcessSignature(PID::BCPLUS, {PID::JPSI, PID::NU_E, PID::POSITRON});
        addIntegrationBoundaries(
            {PS::makeQ2Function(pdg.getMass(PID::POSITRON), pdg.getMass(PID::BCPLUS) - pdg.getMass(PID::PID::JPSI))});
        addTensor(Tensor{name, MD::makeEmptySparse(
                                   dims, {WILSON_BCENU, FF_BCJPSI, WILSON_BCENU_HC, FF_BCJPSI_HC, INTEGRATION_INDEX})});

        setSignatureIndex();
    }

    Tensor RateBDstarLepNu::evalAtPSPoint(const vector<double>& point) {
        auto labs = getTensor().labels();
        labs.pop_back();
        auto dimensions = getTensor().dims();
        dimensions.pop_back();
        Tensor result{"RateBDstarLepNu", MD::makeEmptySparse(dimensions, labs)};

        const double Mb = masses()[0];
        const double Mc = masses()[1];
        const double Mt = masses()[3];

        const double Mb2 = pow(Mb, 2.);
        const double Mb4 = pow(Mb, 4.);
        const double Mc2 = pow(Mc, 2.);
        const double Mc4 = pow(Mc, 4.);
        const double Mt2 = pow(Mt, 2.);

        // kinematic objects
        const double Sqq = point[0];
        const double Sqq2 = pow(Sqq, 2.);
        const double Ew = (Mb2 - Mc2 + Sqq) / (2 * Mb);
        const double Ew2 = pow(Ew, 2.);
        const double Pw = sqrt(Ew2 - Sqq);
        const double Pw2 = pow(Pw, 2.);
        const double Pw4 = pow(Pw, 4.);

        double RateNorm = (pow(GFermi, 2.) * Pw * pow(-Mt2 + Sqq, 2.)) / (32. * Mb2 * pi3 * Sqq);

        // set non-zero tensor elements
        result.element({0, 1, 0, 1}) =
            (2. * Mb2 * Pw2 * (2. * Mt2 + Sqq) + 3. * Mc2 * Sqq * (Mt2 + 2. * Sqq)) / (6. * Mc2 * Sqq2);
        result.element({0, 1, 0, 3}) = (Mb2 * Mt2 * Pw2) / (2. * Mc2 * Sqq);
        result.element({0, 1, 0, 4}) =
            (Mb2 * Pw2 * (2. * Mb2 * (2. * Mt2 + Sqq) - 2. * Mc2 * (2. * Mt2 + Sqq) - Sqq * (Mt2 + 2. * Sqq))) /
            (6. * Mc2 * Sqq2);
        result.element({0, 1, 1, 0}) = -(Mb2 * Mt * Pw2) / (2. * Mc2 * Sqq);
        result.element({0, 1, 3, 0}) = (Mb2 * Mt * Pw2) / (2. * Mc2 * Sqq);
        result.element({0, 1, 5, 1}) =
            -(2. * Mb2 * Pw2 * (2. * Mt2 + Sqq) + 3. * Mc2 * Sqq * (Mt2 + 2. * Sqq)) / (6. * Mc2 * Sqq2);
        result.element({0, 1, 5, 3}) = -(Mb2 * Mt2 * Pw2) / (2. * Mc2 * Sqq);
        result.element({0, 1, 5, 4}) =
            -(Mb2 * Pw2 * (2. * Mb2 * (2. * Mt2 + Sqq) - 2. * Mc2 * (2. * Mt2 + Sqq) - Sqq * (Mt2 + 2. * Sqq))) /
            (6. * Mc2 * Sqq2);
        result.element({0, 1, 7, 1}) =
            (2. * Mb2 * Pw2 * (2. * Mt2 + Sqq) + 3. * Mc2 * Sqq * (Mt2 + 2. * Sqq)) / (6. * Mc2 * Sqq2);
        result.element({0, 1, 7, 3}) = (Mb2 * Mt2 * Pw2) / (2. * Mc2 * Sqq);
        result.element({0, 1, 7, 4}) =
            (Mb2 * Pw2 * (2. * Mb2 * (2. * Mt2 + Sqq) - 2. * Mc2 * (2. * Mt2 + Sqq) - Sqq * (Mt2 + 2. * Sqq))) /
            (6. * Mc2 * Sqq2);
        result.element({0, 1, 9, 5}) = (2. * Mb2 * Mt * Pw2 * (-Mb2 + Mc2 + Sqq)) / (Mc2 * Sqq);
        result.element({0, 1, 9, 6}) = -2. * Mt * (3 + (Mb2 * Pw2) / (Mc2 * Sqq));
        result.element({0, 1, 9, 7}) = (-2. * Mt * (-3. * Mc4 + Mb2 * (3. * Mc2 + Pw2))) / (Mc2 * Sqq);
        result.element({0, 2, 0, 2}) = (4. * Mb2 * Pw2 * (Mt2 + 2. * Sqq)) / (3. * Sqq);
        result.element({0, 2, 5, 2}) = (4. * Mb2 * Pw2 * (Mt2 + 2. * Sqq)) / (3. * Sqq);
        result.element({0, 2, 7, 2}) = (4. * Mb2 * Pw2 * (Mt2 + 2. * Sqq)) / (3. * Sqq);
        result.element({0, 2, 9, 7}) = (16. * Mb2 * Mt * Pw2) / Sqq;
        result.element({0, 3, 0, 1}) = (Mb2 * Mt2 * Pw2) / (2. * Mc2 * Sqq);
        result.element({0, 3, 0, 3}) = (Mb2 * Mt2 * Pw2) / (2. * Mc2);
        result.element({0, 3, 0, 4}) = ((Mb4 - Mb2 * Mc2) * Mt2 * Pw2) / (2. * Mc2 * Sqq);
        result.element({0, 3, 1, 0}) = -(Mb2 * Mt * Pw2) / (2. * Mc2);
        result.element({0, 3, 3, 0}) = (Mb2 * Mt * Pw2) / (2. * Mc2);
        result.element({0, 3, 5, 1}) = -(Mb2 * Mt2 * Pw2) / (2. * Mc2 * Sqq);
        result.element({0, 3, 5, 3}) = -(Mb2 * Mt2 * Pw2) / (2. * Mc2);
        result.element({0, 3, 5, 4}) = ((-Mb4 + Mb2 * Mc2) * Mt2 * Pw2) / (2. * Mc2 * Sqq);
        result.element({0, 3, 7, 1}) = (Mb2 * Mt2 * Pw2) / (2. * Mc2 * Sqq);
        result.element({0, 3, 7, 3}) = (Mb2 * Mt2 * Pw2) / (2. * Mc2);
        result.element({0, 3, 7, 4}) = ((Mb4 - Mb2 * Mc2) * Mt2 * Pw2) / (2. * Mc2 * Sqq);
        result.element({0, 4, 0, 1}) =
            (Mb2 * Pw2 * (2. * Mb2 * (2. * Mt2 + Sqq) - 2. * Mc2 * (2. * Mt2 + Sqq) - Sqq * (Mt2 + 2. * Sqq))) /
            (6. * Mc2 * Sqq2);
        result.element({0, 4, 0, 3}) = ((Mb4 - Mb2 * Mc2) * Mt2 * Pw2) / (2. * Mc2 * Sqq);
        result.element({0, 4, 0, 4}) =
            (Mb2 * Mt2 * Pw2 * (3. * pow(Mb2 - Mc2, 2.) + 4. * Mb2 * Pw2) + 8. * Mb4 * Pw4 * Sqq) / (6. * Mc2 * Sqq2);
        result.element({0, 4, 1, 0}) = ((-Mb4 + Mb2 * Mc2) * Mt * Pw2) / (2. * Mc2 * Sqq);
        result.element({0, 4, 3, 0}) = ((Mb4 - Mb2 * Mc2) * Mt * Pw2) / (2. * Mc2 * Sqq);
        result.element({0, 4, 5, 1}) =
            -(Mb2 * Pw2 * (2. * Mb2 * (2. * Mt2 + Sqq) - 2. * Mc2 * (2. * Mt2 + Sqq) - Sqq * (Mt2 + 2. * Sqq))) /
            (6. * Mc2 * Sqq2);
        result.element({0, 4, 5, 3}) = ((-Mb4 + Mb2 * Mc2) * Mt2 * Pw2) / (2. * Mc2 * Sqq);
        result.element({0, 4, 5, 4}) =
            (Mb2 * Mt2 * Pw2 * (-3. * pow(Mb2 - Mc2, 2.) - 4. * Mb2 * Pw2) - 8. * Mb4 * Pw4 * Sqq) / (6. * Mc2 * Sqq2);
        result.element({0, 4, 7, 1}) =
            (Mb2 * Pw2 * (2. * Mb2 * (2. * Mt2 + Sqq) - 2. * Mc2 * (2. * Mt2 + Sqq) - Sqq * (Mt2 + 2. * Sqq))) /
            (6. * Mc2 * Sqq2);
        result.element({0, 4, 7, 3}) = ((Mb4 - Mb2 * Mc2) * Mt2 * Pw2) / (2. * Mc2 * Sqq);
        result.element({0, 4, 7, 4}) =
            (Mb2 * Mt2 * Pw2 * (3. * pow(Mb2 - Mc2, 2.) + 4. * Mb2 * Pw2) + 8. * Mb4 * Pw4 * Sqq) / (6. * Mc2 * Sqq2);
        result.element({0, 4, 9, 5}) = (-8. * Mb4 * Mt * Pw4) / (Mc2 * Sqq);
        result.element({0, 4, 9, 6}) = (-2. * Mb2 * Mt * Pw2 * (Mb2 - Mc2 - Sqq)) / (Mc2 * Sqq);
        result.element({0, 4, 9, 7}) = (-2. * Mb2 * Mt * Pw2 * (Mb2 + 3. * Mc2 - Sqq)) / (Mc2 * Sqq);
        result.element({1, 0, 0, 1}) = -(Mb2 * Mt * Pw2) / (2. * Mc2 * Sqq);
        result.element({1, 0, 0, 3}) = -(Mb2 * Mt * Pw2) / (2. * Mc2);
        result.element({1, 0, 0, 4}) = ((-Mb4 + Mb2 * Mc2) * Mt * Pw2) / (2. * Mc2 * Sqq);
        result.element({1, 0, 1, 0}) = (Mb2 * Pw2) / (2. * Mc2);
        result.element({1, 0, 3, 0}) = -(Mb2 * Pw2) / (2. * Mc2);
        result.element({1, 0, 5, 1}) = (Mb2 * Mt * Pw2) / (2. * Mc2 * Sqq);
        result.element({1, 0, 5, 3}) = (Mb2 * Mt * Pw2) / (2. * Mc2);
        result.element({1, 0, 5, 4}) = ((Mb4 - Mb2 * Mc2) * Mt * Pw2) / (2. * Mc2 * Sqq);
        result.element({1, 0, 7, 1}) = -(Mb2 * Mt * Pw2) / (2. * Mc2 * Sqq);
        result.element({1, 0, 7, 3}) = -(Mb2 * Mt * Pw2) / (2. * Mc2);
        result.element({1, 0, 7, 4}) = ((-Mb4 + Mb2 * Mc2) * Mt * Pw2) / (2. * Mc2 * Sqq);
        result.element({2, 0, 2, 0}) = (Mb2 * Pw2) / (2. * Mc2);
        result.element({2, 0, 4, 0}) = -(Mb2 * Pw2) / (2. * Mc2);
        result.element({2, 0, 6, 1}) = (Mb2 * Mt * Pw2) / (2. * Mc2 * Sqq);
        result.element({2, 0, 6, 3}) = (Mb2 * Mt * Pw2) / (2. * Mc2);
        result.element({2, 0, 6, 4}) = ((Mb4 - Mb2 * Mc2) * Mt * Pw2) / (2. * Mc2 * Sqq);
        result.element({2, 0, 8, 1}) = -(Mb2 * Mt * Pw2) / (2. * Mc2 * Sqq);
        result.element({2, 0, 8, 3}) = -(Mb2 * Mt * Pw2) / (2. * Mc2);
        result.element({2, 0, 8, 4}) = ((-Mb4 + Mb2 * Mc2) * Mt * Pw2) / (2. * Mc2 * Sqq);
        result.element({3, 0, 0, 1}) = (Mb2 * Mt * Pw2) / (2. * Mc2 * Sqq);
        result.element({3, 0, 0, 3}) = (Mb2 * Mt * Pw2) / (2. * Mc2);
        result.element({3, 0, 0, 4}) = ((Mb4 - Mb2 * Mc2) * Mt * Pw2) / (2. * Mc2 * Sqq);
        result.element({3, 0, 1, 0}) = -(Mb2 * Pw2) / (2. * Mc2);
        result.element({3, 0, 3, 0}) = (Mb2 * Pw2) / (2. * Mc2);
        result.element({3, 0, 5, 1}) = -(Mb2 * Mt * Pw2) / (2. * Mc2 * Sqq);
        result.element({3, 0, 5, 3}) = -(Mb2 * Mt * Pw2) / (2. * Mc2);
        result.element({3, 0, 5, 4}) = ((-Mb4 + Mb2 * Mc2) * Mt * Pw2) / (2. * Mc2 * Sqq);
        result.element({3, 0, 7, 1}) = (Mb2 * Mt * Pw2) / (2. * Mc2 * Sqq);
        result.element({3, 0, 7, 3}) = (Mb2 * Mt * Pw2) / (2. * Mc2);
        result.element({3, 0, 7, 4}) = ((Mb4 - Mb2 * Mc2) * Mt * Pw2) / (2. * Mc2 * Sqq);
        result.element({4, 0, 2, 0}) = -(Mb2 * Pw2) / (2. * Mc2);
        result.element({4, 0, 4, 0}) = (Mb2 * Pw2) / (2. * Mc2);
        result.element({4, 0, 6, 1}) = -(Mb2 * Mt * Pw2) / (2. * Mc2 * Sqq);
        result.element({4, 0, 6, 3}) = -(Mb2 * Mt * Pw2) / (2. * Mc2);
        result.element({4, 0, 6, 4}) = ((-Mb4 + Mb2 * Mc2) * Mt * Pw2) / (2. * Mc2 * Sqq);
        result.element({4, 0, 8, 1}) = (Mb2 * Mt * Pw2) / (2. * Mc2 * Sqq);
        result.element({4, 0, 8, 3}) = (Mb2 * Mt * Pw2) / (2. * Mc2);
        result.element({4, 0, 8, 4}) = ((Mb4 - Mb2 * Mc2) * Mt * Pw2) / (2. * Mc2 * Sqq);
        result.element({5, 1, 0, 1}) =
            -(2. * Mb2 * Pw2 * (2. * Mt2 + Sqq) + 3. * Mc2 * Sqq * (Mt2 + 2. * Sqq)) / (6. * Mc2 * Sqq2);
        result.element({5, 1, 0, 3}) = -(Mb2 * Mt2 * Pw2) / (2. * Mc2 * Sqq);
        result.element({5, 1, 0, 4}) =
            -(Mb2 * Pw2 * (2. * Mb2 * (2. * Mt2 + Sqq) - 2. * Mc2 * (2. * Mt2 + Sqq) - Sqq * (Mt2 + 2. * Sqq))) /
            (6. * Mc2 * Sqq2);
        result.element({5, 1, 1, 0}) = (Mb2 * Mt * Pw2) / (2. * Mc2 * Sqq);
        result.element({5, 1, 3, 0}) = -(Mb2 * Mt * Pw2) / (2. * Mc2 * Sqq);
        result.element({5, 1, 5, 1}) =
            (2. * Mb2 * Pw2 * (2. * Mt2 + Sqq) + 3. * Mc2 * Sqq * (Mt2 + 2. * Sqq)) / (6. * Mc2 * Sqq2);
        result.element({5, 1, 5, 3}) = (Mb2 * Mt2 * Pw2) / (2. * Mc2 * Sqq);
        result.element({5, 1, 5, 4}) =
            (Mb2 * Pw2 * (2. * Mb2 * (2. * Mt2 + Sqq) - 2. * Mc2 * (2. * Mt2 + Sqq) - Sqq * (Mt2 + 2. * Sqq))) /
            (6. * Mc2 * Sqq2);
        result.element({5, 1, 7, 1}) =
            -(2. * Mb2 * Pw2 * (2. * Mt2 + Sqq) + 3. * Mc2 * Sqq * (Mt2 + 2. * Sqq)) / (6. * Mc2 * Sqq2);
        result.element({5, 1, 7, 3}) = -(Mb2 * Mt2 * Pw2) / (2. * Mc2 * Sqq);
        result.element({5, 1, 7, 4}) =
            -(Mb2 * Pw2 * (2. * Mb2 * (2. * Mt2 + Sqq) - 2. * Mc2 * (2. * Mt2 + Sqq) - Sqq * (Mt2 + 2. * Sqq))) /
            (6. * Mc2 * Sqq2);
        result.element({5, 1, 9, 5}) = (-2. * Mb2 * Mt * Pw2 * (-Mb2 + Mc2 + Sqq)) / (Mc2 * Sqq);
        result.element({5, 1, 9, 6}) = 2. * Mt * (3 + (Mb2 * Pw2) / (Mc2 * Sqq));
        result.element({5, 1, 9, 7}) = (2. * Mt * (-3. * Mc4 + Mb2 * (3. * Mc2 + Pw2))) / (Mc2 * Sqq);
        result.element({5, 2, 0, 2}) = (4. * Mb2 * Pw2 * (Mt2 + 2. * Sqq)) / (3. * Sqq);
        result.element({5, 2, 5, 2}) = (4. * Mb2 * Pw2 * (Mt2 + 2. * Sqq)) / (3. * Sqq);
        result.element({5, 2, 7, 2}) = (4. * Mb2 * Pw2 * (Mt2 + 2. * Sqq)) / (3. * Sqq);
        result.element({5, 2, 9, 7}) = (16. * Mb2 * Mt * Pw2) / Sqq;
        result.element({5, 3, 0, 1}) = -(Mb2 * Mt2 * Pw2) / (2. * Mc2 * Sqq);
        result.element({5, 3, 0, 3}) = -(Mb2 * Mt2 * Pw2) / (2. * Mc2);
        result.element({5, 3, 0, 4}) = ((-Mb4 + Mb2 * Mc2) * Mt2 * Pw2) / (2. * Mc2 * Sqq);
        result.element({5, 3, 1, 0}) = (Mb2 * Mt * Pw2) / (2. * Mc2);
        result.element({5, 3, 3, 0}) = -(Mb2 * Mt * Pw2) / (2. * Mc2);
        result.element({5, 3, 5, 1}) = (Mb2 * Mt2 * Pw2) / (2. * Mc2 * Sqq);
        result.element({5, 3, 5, 3}) = (Mb2 * Mt2 * Pw2) / (2. * Mc2);
        result.element({5, 3, 5, 4}) = ((Mb4 - Mb2 * Mc2) * Mt2 * Pw2) / (2. * Mc2 * Sqq);
        result.element({5, 3, 7, 1}) = -(Mb2 * Mt2 * Pw2) / (2. * Mc2 * Sqq);
        result.element({5, 3, 7, 3}) = -(Mb2 * Mt2 * Pw2) / (2. * Mc2);
        result.element({5, 3, 7, 4}) = ((-Mb4 + Mb2 * Mc2) * Mt2 * Pw2) / (2. * Mc2 * Sqq);
        result.element({5, 4, 0, 1}) =
            -(Mb2 * Pw2 * (2. * Mb2 * (2. * Mt2 + Sqq) - 2. * Mc2 * (2. * Mt2 + Sqq) - Sqq * (Mt2 + 2. * Sqq))) /
            (6. * Mc2 * Sqq2);
        result.element({5, 4, 0, 3}) = ((-Mb4 + Mb2 * Mc2) * Mt2 * Pw2) / (2. * Mc2 * Sqq);
        result.element({5, 4, 0, 4}) =
            (Mb2 * Mt2 * Pw2 * (-3. * pow(Mb2 - Mc2, 2.) - 4. * Mb2 * Pw2) - 8. * Mb4 * Pw4 * Sqq) / (6. * Mc2 * Sqq2);
        result.element({5, 4, 1, 0}) = ((Mb4 - Mb2 * Mc2) * Mt * Pw2) / (2. * Mc2 * Sqq);
        result.element({5, 4, 3, 0}) = ((-Mb4 + Mb2 * Mc2) * Mt * Pw2) / (2. * Mc2 * Sqq);
        result.element({5, 4, 5, 1}) =
            (Mb2 * Pw2 * (2. * Mb2 * (2. * Mt2 + Sqq) - 2. * Mc2 * (2. * Mt2 + Sqq) - Sqq * (Mt2 + 2. * Sqq))) /
            (6. * Mc2 * Sqq2);
        result.element({5, 4, 5, 3}) = ((Mb4 - Mb2 * Mc2) * Mt2 * Pw2) / (2. * Mc2 * Sqq);
        result.element({5, 4, 5, 4}) =
            (Mb2 * Mt2 * Pw2 * (3. * pow(Mb2 - Mc2, 2.) + 4. * Mb2 * Pw2) + 8. * Mb4 * Pw4 * Sqq) / (6. * Mc2 * Sqq2);
        result.element({5, 4, 7, 1}) =
            -(Mb2 * Pw2 * (2. * Mb2 * (2. * Mt2 + Sqq) - 2. * Mc2 * (2. * Mt2 + Sqq) - Sqq * (Mt2 + 2. * Sqq))) /
            (6. * Mc2 * Sqq2);
        result.element({5, 4, 7, 3}) = ((-Mb4 + Mb2 * Mc2) * Mt2 * Pw2) / (2. * Mc2 * Sqq);
        result.element({5, 4, 7, 4}) =
            (Mb2 * Mt2 * Pw2 * (-3. * pow(Mb2 - Mc2, 2.) - 4. * Mb2 * Pw2) - 8. * Mb4 * Pw4 * Sqq) / (6. * Mc2 * Sqq2);
        result.element({5, 4, 9, 5}) = (8. * Mb4 * Mt * Pw4) / (Mc2 * Sqq);
        result.element({5, 4, 9, 6}) = (2. * Mb2 * Mt * Pw2 * (Mb2 - Mc2 - Sqq)) / (Mc2 * Sqq);
        result.element({5, 4, 9, 7}) = (2. * Mb2 * Mt * Pw2 * (Mb2 + 3. * Mc2 - Sqq)) / (Mc2 * Sqq);
        result.element({6, 1, 2, 0}) = (Mb2 * Mt * Pw2) / (2. * Mc2 * Sqq);
        result.element({6, 1, 4, 0}) = -(Mb2 * Mt * Pw2) / (2. * Mc2 * Sqq);
        result.element({6, 1, 6, 1}) =
            (2. * Mb2 * Pw2 * (2. * Mt2 + Sqq) + 3. * Mc2 * Sqq * (Mt2 + 2. * Sqq)) / (6. * Mc2 * Sqq2);
        result.element({6, 1, 6, 3}) = (Mb2 * Mt2 * Pw2) / (2. * Mc2 * Sqq);
        result.element({6, 1, 6, 4}) =
            (Mb2 * Pw2 * (2. * Mb2 * (2. * Mt2 + Sqq) - 2. * Mc2 * (2. * Mt2 + Sqq) - Sqq * (Mt2 + 2. * Sqq))) /
            (6. * Mc2 * Sqq2);
        result.element({6, 1, 8, 1}) =
            -(2. * Mb2 * Pw2 * (2. * Mt2 + Sqq) + 3. * Mc2 * Sqq * (Mt2 + 2. * Sqq)) / (6. * Mc2 * Sqq2);
        result.element({6, 1, 8, 3}) = -(Mb2 * Mt2 * Pw2) / (2. * Mc2 * Sqq);
        result.element({6, 1, 8, 4}) =
            -(Mb2 * Pw2 * (2. * Mb2 * (2. * Mt2 + Sqq) - 2. * Mc2 * (2. * Mt2 + Sqq) - Sqq * (Mt2 + 2. * Sqq))) /
            (6. * Mc2 * Sqq2);
        result.element({6, 1, 10, 5}) = (2. * Mb2 * Mt * Pw2 * (-Mb2 + Mc2 + Sqq)) / (Mc2 * Sqq);
        result.element({6, 1, 10, 6}) = -2. * Mt * (3 + (Mb2 * Pw2) / (Mc2 * Sqq));
        result.element({6, 1, 10, 7}) = (-2. * Mt * (-3. * Mc4 + Mb2 * (3. * Mc2 + Pw2))) / (Mc2 * Sqq);
        result.element({6, 2, 6, 2}) = (4. * Mb2 * Pw2 * (Mt2 + 2. * Sqq)) / (3. * Sqq);
        result.element({6, 2, 8, 2}) = (4. * Mb2 * Pw2 * (Mt2 + 2. * Sqq)) / (3. * Sqq);
        result.element({6, 2, 10, 7}) = (16. * Mb2 * Mt * Pw2) / Sqq;
        result.element({6, 3, 2, 0}) = (Mb2 * Mt * Pw2) / (2. * Mc2);
        result.element({6, 3, 4, 0}) = -(Mb2 * Mt * Pw2) / (2. * Mc2);
        result.element({6, 3, 6, 1}) = (Mb2 * Mt2 * Pw2) / (2. * Mc2 * Sqq);
        result.element({6, 3, 6, 3}) = (Mb2 * Mt2 * Pw2) / (2. * Mc2);
        result.element({6, 3, 6, 4}) = ((Mb4 - Mb2 * Mc2) * Mt2 * Pw2) / (2. * Mc2 * Sqq);
        result.element({6, 3, 8, 1}) = -(Mb2 * Mt2 * Pw2) / (2. * Mc2 * Sqq);
        result.element({6, 3, 8, 3}) = -(Mb2 * Mt2 * Pw2) / (2. * Mc2);
        result.element({6, 3, 8, 4}) = ((-Mb4 + Mb2 * Mc2) * Mt2 * Pw2) / (2. * Mc2 * Sqq);
        result.element({6, 4, 2, 0}) = ((Mb4 - Mb2 * Mc2) * Mt * Pw2) / (2. * Mc2 * Sqq);
        result.element({6, 4, 4, 0}) = ((-Mb4 + Mb2 * Mc2) * Mt * Pw2) / (2. * Mc2 * Sqq);
        result.element({6, 4, 6, 1}) =
            (Mb2 * Pw2 * (2. * Mb2 * (2. * Mt2 + Sqq) - 2. * Mc2 * (2. * Mt2 + Sqq) - Sqq * (Mt2 + 2. * Sqq))) /
            (6. * Mc2 * Sqq2);
        result.element({6, 4, 6, 3}) = ((Mb4 - Mb2 * Mc2) * Mt2 * Pw2) / (2. * Mc2 * Sqq);
        result.element({6, 4, 6, 4}) =
            (Mb2 * Mt2 * Pw2 * (3. * pow(Mb2 - Mc2, 2.) + 4. * Mb2 * Pw2) + 8. * Mb4 * Pw4 * Sqq) / (6. * Mc2 * Sqq2);
        result.element({6, 4, 8, 1}) =
            -(Mb2 * Pw2 * (2. * Mb2 * (2. * Mt2 + Sqq) - 2. * Mc2 * (2. * Mt2 + Sqq) - Sqq * (Mt2 + 2. * Sqq))) /
            (6. * Mc2 * Sqq2);
        result.element({6, 4, 8, 3}) = ((-Mb4 + Mb2 * Mc2) * Mt2 * Pw2) / (2. * Mc2 * Sqq);
        result.element({6, 4, 8, 4}) =
            (Mb2 * Mt2 * Pw2 * (-3. * pow(Mb2 - Mc2, 2.) - 4. * Mb2 * Pw2) - 8. * Mb4 * Pw4 * Sqq) / (6. * Mc2 * Sqq2);
        result.element({6, 4, 10, 5}) = (-8. * Mb4 * Mt * Pw4) / (Mc2 * Sqq);
        result.element({6, 4, 10, 6}) = (-2. * Mb2 * Mt * Pw2 * (Mb2 - Mc2 - Sqq)) / (Mc2 * Sqq);
        result.element({6, 4, 10, 7}) = (-2. * Mb2 * Mt * Pw2 * (Mb2 + 3. * Mc2 - Sqq)) / (Mc2 * Sqq);
        result.element({7, 1, 0, 1}) =
            (2. * Mb2 * Pw2 * (2. * Mt2 + Sqq) + 3. * Mc2 * Sqq * (Mt2 + 2. * Sqq)) / (6. * Mc2 * Sqq2);
        result.element({7, 1, 0, 3}) = (Mb2 * Mt2 * Pw2) / (2. * Mc2 * Sqq);
        result.element({7, 1, 0, 4}) =
            (Mb2 * Pw2 * (2. * Mb2 * (2. * Mt2 + Sqq) - 2. * Mc2 * (2. * Mt2 + Sqq) - Sqq * (Mt2 + 2. * Sqq))) /
            (6. * Mc2 * Sqq2);
        result.element({7, 1, 1, 0}) = -(Mb2 * Mt * Pw2) / (2. * Mc2 * Sqq);
        result.element({7, 1, 3, 0}) = (Mb2 * Mt * Pw2) / (2. * Mc2 * Sqq);
        result.element({7, 1, 5, 1}) =
            -(2. * Mb2 * Pw2 * (2. * Mt2 + Sqq) + 3. * Mc2 * Sqq * (Mt2 + 2. * Sqq)) / (6. * Mc2 * Sqq2);
        result.element({7, 1, 5, 3}) = -(Mb2 * Mt2 * Pw2) / (2. * Mc2 * Sqq);
        result.element({7, 1, 5, 4}) =
            -(Mb2 * Pw2 * (2. * Mb2 * (2. * Mt2 + Sqq) - 2. * Mc2 * (2. * Mt2 + Sqq) - Sqq * (Mt2 + 2. * Sqq))) /
            (6. * Mc2 * Sqq2);
        result.element({7, 1, 7, 1}) =
            (2. * Mb2 * Pw2 * (2. * Mt2 + Sqq) + 3. * Mc2 * Sqq * (Mt2 + 2. * Sqq)) / (6. * Mc2 * Sqq2);
        result.element({7, 1, 7, 3}) = (Mb2 * Mt2 * Pw2) / (2. * Mc2 * Sqq);
        result.element({7, 1, 7, 4}) =
            (Mb2 * Pw2 * (2. * Mb2 * (2. * Mt2 + Sqq) - 2. * Mc2 * (2. * Mt2 + Sqq) - Sqq * (Mt2 + 2. * Sqq))) /
            (6. * Mc2 * Sqq2);
        result.element({7, 1, 9, 5}) = (2. * Mb2 * Mt * Pw2 * (-Mb2 + Mc2 + Sqq)) / (Mc2 * Sqq);
        result.element({7, 1, 9, 6}) = -2. * Mt * (3 + (Mb2 * Pw2) / (Mc2 * Sqq));
        result.element({7, 1, 9, 7}) = (-2. * Mt * (-3. * Mc4 + Mb2 * (3. * Mc2 + Pw2))) / (Mc2 * Sqq);
        result.element({7, 2, 0, 2}) = (4. * Mb2 * Pw2 * (Mt2 + 2. * Sqq)) / (3. * Sqq);
        result.element({7, 2, 5, 2}) = (4. * Mb2 * Pw2 * (Mt2 + 2. * Sqq)) / (3. * Sqq);
        result.element({7, 2, 7, 2}) = (4. * Mb2 * Pw2 * (Mt2 + 2. * Sqq)) / (3. * Sqq);
        result.element({7, 2, 9, 7}) = (16. * Mb2 * Mt * Pw2) / Sqq;
        result.element({7, 3, 0, 1}) = (Mb2 * Mt2 * Pw2) / (2. * Mc2 * Sqq);
        result.element({7, 3, 0, 3}) = (Mb2 * Mt2 * Pw2) / (2. * Mc2);
        result.element({7, 3, 0, 4}) = ((Mb4 - Mb2 * Mc2) * Mt2 * Pw2) / (2. * Mc2 * Sqq);
        result.element({7, 3, 1, 0}) = -(Mb2 * Mt * Pw2) / (2. * Mc2);
        result.element({7, 3, 3, 0}) = (Mb2 * Mt * Pw2) / (2. * Mc2);
        result.element({7, 3, 5, 1}) = -(Mb2 * Mt2 * Pw2) / (2. * Mc2 * Sqq);
        result.element({7, 3, 5, 3}) = -(Mb2 * Mt2 * Pw2) / (2. * Mc2);
        result.element({7, 3, 5, 4}) = ((-Mb4 + Mb2 * Mc2) * Mt2 * Pw2) / (2. * Mc2 * Sqq);
        result.element({7, 3, 7, 1}) = (Mb2 * Mt2 * Pw2) / (2. * Mc2 * Sqq);
        result.element({7, 3, 7, 3}) = (Mb2 * Mt2 * Pw2) / (2. * Mc2);
        result.element({7, 3, 7, 4}) = ((Mb4 - Mb2 * Mc2) * Mt2 * Pw2) / (2. * Mc2 * Sqq);
        result.element({7, 4, 0, 1}) =
            (Mb2 * Pw2 * (2. * Mb2 * (2. * Mt2 + Sqq) - 2. * Mc2 * (2. * Mt2 + Sqq) - Sqq * (Mt2 + 2. * Sqq))) /
            (6. * Mc2 * Sqq2);
        result.element({7, 4, 0, 3}) = ((Mb4 - Mb2 * Mc2) * Mt2 * Pw2) / (2. * Mc2 * Sqq);
        result.element({7, 4, 0, 4}) =
            (Mb2 * Mt2 * Pw2 * (3. * pow(Mb2 - Mc2, 2.) + 4. * Mb2 * Pw2) + 8. * Mb4 * Pw4 * Sqq) / (6. * Mc2 * Sqq2);
        result.element({7, 4, 1, 0}) = ((-Mb4 + Mb2 * Mc2) * Mt * Pw2) / (2. * Mc2 * Sqq);
        result.element({7, 4, 3, 0}) = ((Mb4 - Mb2 * Mc2) * Mt * Pw2) / (2. * Mc2 * Sqq);
        result.element({7, 4, 5, 1}) =
            -(Mb2 * Pw2 * (2. * Mb2 * (2. * Mt2 + Sqq) - 2. * Mc2 * (2. * Mt2 + Sqq) - Sqq * (Mt2 + 2. * Sqq))) /
            (6. * Mc2 * Sqq2);
        result.element({7, 4, 5, 3}) = ((-Mb4 + Mb2 * Mc2) * Mt2 * Pw2) / (2. * Mc2 * Sqq);
        result.element({7, 4, 5, 4}) =
            (Mb2 * Mt2 * Pw2 * (-3. * pow(Mb2 - Mc2, 2.) - 4. * Mb2 * Pw2) - 8. * Mb4 * Pw4 * Sqq) / (6. * Mc2 * Sqq2);
        result.element({7, 4, 7, 1}) =
            (Mb2 * Pw2 * (2. * Mb2 * (2. * Mt2 + Sqq) - 2. * Mc2 * (2. * Mt2 + Sqq) - Sqq * (Mt2 + 2. * Sqq))) /
            (6. * Mc2 * Sqq2);
        result.element({7, 4, 7, 3}) = ((Mb4 - Mb2 * Mc2) * Mt2 * Pw2) / (2. * Mc2 * Sqq);
        result.element({7, 4, 7, 4}) =
            (Mb2 * Mt2 * Pw2 * (3. * pow(Mb2 - Mc2, 2.) + 4. * Mb2 * Pw2) + 8. * Mb4 * Pw4 * Sqq) / (6. * Mc2 * Sqq2);
        result.element({7, 4, 9, 5}) = (-8. * Mb4 * Mt * Pw4) / (Mc2 * Sqq);
        result.element({7, 4, 9, 6}) = (-2. * Mb2 * Mt * Pw2 * (Mb2 - Mc2 - Sqq)) / (Mc2 * Sqq);
        result.element({7, 4, 9, 7}) = (-2. * Mb2 * Mt * Pw2 * (Mb2 + 3. * Mc2 - Sqq)) / (Mc2 * Sqq);
        result.element({8, 1, 2, 0}) = -(Mb2 * Mt * Pw2) / (2. * Mc2 * Sqq);
        result.element({8, 1, 4, 0}) = (Mb2 * Mt * Pw2) / (2. * Mc2 * Sqq);
        result.element({8, 1, 6, 1}) =
            -(2. * Mb2 * Pw2 * (2. * Mt2 + Sqq) + 3. * Mc2 * Sqq * (Mt2 + 2. * Sqq)) / (6. * Mc2 * Sqq2);
        result.element({8, 1, 6, 3}) = -(Mb2 * Mt2 * Pw2) / (2. * Mc2 * Sqq);
        result.element({8, 1, 6, 4}) =
            -(Mb2 * Pw2 * (2. * Mb2 * (2. * Mt2 + Sqq) - 2. * Mc2 * (2. * Mt2 + Sqq) - Sqq * (Mt2 + 2. * Sqq))) /
            (6. * Mc2 * Sqq2);
        result.element({8, 1, 8, 1}) =
            (2. * Mb2 * Pw2 * (2. * Mt2 + Sqq) + 3. * Mc2 * Sqq * (Mt2 + 2. * Sqq)) / (6. * Mc2 * Sqq2);
        result.element({8, 1, 8, 3}) = (Mb2 * Mt2 * Pw2) / (2. * Mc2 * Sqq);
        result.element({8, 1, 8, 4}) =
            (Mb2 * Pw2 * (2. * Mb2 * (2. * Mt2 + Sqq) - 2. * Mc2 * (2. * Mt2 + Sqq) - Sqq * (Mt2 + 2. * Sqq))) /
            (6. * Mc2 * Sqq2);
        result.element({8, 1, 10, 5}) = (-2. * Mb2 * Mt * Pw2 * (-Mb2 + Mc2 + Sqq)) / (Mc2 * Sqq);
        result.element({8, 1, 10, 6}) = 2. * Mt * (3 + (Mb2 * Pw2) / (Mc2 * Sqq));
        result.element({8, 1, 10, 7}) = (2. * Mt * (-3. * Mc4 + Mb2 * (3. * Mc2 + Pw2))) / (Mc2 * Sqq);
        result.element({8, 2, 6, 2}) = (4. * Mb2 * Pw2 * (Mt2 + 2. * Sqq)) / (3. * Sqq);
        result.element({8, 2, 8, 2}) = (4. * Mb2 * Pw2 * (Mt2 + 2. * Sqq)) / (3. * Sqq);
        result.element({8, 2, 10, 7}) = (16. * Mb2 * Mt * Pw2) / Sqq;
        result.element({8, 3, 2, 0}) = -(Mb2 * Mt * Pw2) / (2. * Mc2);
        result.element({8, 3, 4, 0}) = (Mb2 * Mt * Pw2) / (2. * Mc2);
        result.element({8, 3, 6, 1}) = -(Mb2 * Mt2 * Pw2) / (2. * Mc2 * Sqq);
        result.element({8, 3, 6, 3}) = -(Mb2 * Mt2 * Pw2) / (2. * Mc2);
        result.element({8, 3, 6, 4}) = ((-Mb4 + Mb2 * Mc2) * Mt2 * Pw2) / (2. * Mc2 * Sqq);
        result.element({8, 3, 8, 1}) = (Mb2 * Mt2 * Pw2) / (2. * Mc2 * Sqq);
        result.element({8, 3, 8, 3}) = (Mb2 * Mt2 * Pw2) / (2. * Mc2);
        result.element({8, 3, 8, 4}) = ((Mb4 - Mb2 * Mc2) * Mt2 * Pw2) / (2. * Mc2 * Sqq);
        result.element({8, 4, 2, 0}) = ((-Mb4 + Mb2 * Mc2) * Mt * Pw2) / (2. * Mc2 * Sqq);
        result.element({8, 4, 4, 0}) = ((Mb4 - Mb2 * Mc2) * Mt * Pw2) / (2. * Mc2 * Sqq);
        result.element({8, 4, 6, 1}) =
            -(Mb2 * Pw2 * (2. * Mb2 * (2. * Mt2 + Sqq) - 2. * Mc2 * (2. * Mt2 + Sqq) - Sqq * (Mt2 + 2. * Sqq))) /
            (6. * Mc2 * Sqq2);
        result.element({8, 4, 6, 3}) = ((-Mb4 + Mb2 * Mc2) * Mt2 * Pw2) / (2. * Mc2 * Sqq);
        result.element({8, 4, 6, 4}) =
            (Mb2 * Mt2 * Pw2 * (-3. * pow(Mb2 - Mc2, 2.) - 4. * Mb2 * Pw2) - 8. * Mb4 * Pw4 * Sqq) / (6. * Mc2 * Sqq2);
        result.element({8, 4, 8, 1}) =
            (Mb2 * Pw2 * (2. * Mb2 * (2. * Mt2 + Sqq) - 2. * Mc2 * (2. * Mt2 + Sqq) - Sqq * (Mt2 + 2. * Sqq))) /
            (6. * Mc2 * Sqq2);
        result.element({8, 4, 8, 3}) = ((Mb4 - Mb2 * Mc2) * Mt2 * Pw2) / (2. * Mc2 * Sqq);
        result.element({8, 4, 8, 4}) =
            (Mb2 * Mt2 * Pw2 * (3. * pow(Mb2 - Mc2, 2.) + 4. * Mb2 * Pw2) + 8. * Mb4 * Pw4 * Sqq) / (6. * Mc2 * Sqq2);
        result.element({8, 4, 10, 5}) = (8. * Mb4 * Mt * Pw4) / (Mc2 * Sqq);
        result.element({8, 4, 10, 6}) = (2. * Mb2 * Mt * Pw2 * (Mb2 - Mc2 - Sqq)) / (Mc2 * Sqq);
        result.element({8, 4, 10, 7}) = (2. * Mb2 * Mt * Pw2 * (Mb2 + 3. * Mc2 - Sqq)) / (Mc2 * Sqq);
        result.element({9, 5, 0, 1}) = (2. * Mb2 * Mt * Pw2 * (-Mb2 + Mc2 + Sqq)) / (Mc2 * Sqq);
        result.element({9, 5, 0, 4}) = (-8. * Mb4 * Mt * Pw4) / (Mc2 * Sqq);
        result.element({9, 5, 5, 1}) = (-2. * Mb2 * Mt * Pw2 * (-Mb2 + Mc2 + Sqq)) / (Mc2 * Sqq);
        result.element({9, 5, 5, 4}) = (8. * Mb4 * Mt * Pw4) / (Mc2 * Sqq);
        result.element({9, 5, 7, 1}) = (2. * Mb2 * Mt * Pw2 * (-Mb2 + Mc2 + Sqq)) / (Mc2 * Sqq);
        result.element({9, 5, 7, 4}) = (-8. * Mb4 * Mt * Pw4) / (Mc2 * Sqq);
        result.element({9, 5, 9, 5}) = (32. * Mb4 * Pw4 * (2. * Mt2 + Sqq)) / (3. * Mc2 * Sqq);
        result.element({9, 5, 9, 6}) = (8. * Mb2 * Pw2 * (Mb2 - Mc2 - Sqq) * (2. * Mt2 + Sqq)) / (3. * Mc2 * Sqq);
        result.element({9, 5, 9, 7}) = (8. * Mb2 * Pw2 * (Mb2 + 3. * Mc2 - Sqq) * (2. * Mt2 + Sqq)) / (3. * Mc2 * Sqq);
        result.element({9, 6, 0, 1}) = -2. * Mt * (3 + (Mb2 * Pw2) / (Mc2 * Sqq));
        result.element({9, 6, 0, 4}) = (-2. * Mb2 * Mt * Pw2 * (Mb2 - Mc2 - Sqq)) / (Mc2 * Sqq);
        result.element({9, 6, 5, 1}) = 2. * Mt * (3 + (Mb2 * Pw2) / (Mc2 * Sqq));
        result.element({9, 6, 5, 4}) = (2. * Mb2 * Mt * Pw2 * (Mb2 - Mc2 - Sqq)) / (Mc2 * Sqq);
        result.element({9, 6, 7, 1}) = -2. * Mt * (3 + (Mb2 * Pw2) / (Mc2 * Sqq));
        result.element({9, 6, 7, 4}) = (-2. * Mb2 * Mt * Pw2 * (Mb2 - Mc2 - Sqq)) / (Mc2 * Sqq);
        result.element({9, 6, 9, 5}) = (8. * Mb2 * Pw2 * (Mb2 - Mc2 - Sqq) * (2. * Mt2 + Sqq)) / (3. * Mc2 * Sqq);
        result.element({9, 6, 9, 6}) = (8. * (2. * Mt2 + Sqq) * (Mb2 * Pw2 + 3. * Mc2 * Sqq)) / (3. * Mc2 * Sqq);
        result.element({9, 6, 9, 7}) =
            (8. * (-3. * Mc4 + Mb2 * (3. * Mc2 + Pw2)) * (2. * Mt2 + Sqq)) / (3. * Mc2 * Sqq);
        result.element({9, 7, 0, 1}) = (-2. * Mt * (-3. * Mc4 + Mb2 * (3. * Mc2 + Pw2))) / (Mc2 * Sqq);
        result.element({9, 7, 0, 2}) = (16. * Mb2 * Mt * Pw2) / Sqq;
        result.element({9, 7, 0, 4}) = (-2. * Mb2 * Mt * Pw2 * (Mb2 + 3. * Mc2 - Sqq)) / (Mc2 * Sqq);
        result.element({9, 7, 5, 1}) = (2. * Mt * (-3. * Mc4 + Mb2 * (3. * Mc2 + Pw2))) / (Mc2 * Sqq);
        result.element({9, 7, 5, 2}) = (16. * Mb2 * Mt * Pw2) / Sqq;
        result.element({9, 7, 5, 4}) = (2. * Mb2 * Mt * Pw2 * (Mb2 + 3. * Mc2 - Sqq)) / (Mc2 * Sqq);
        result.element({9, 7, 7, 1}) = (-2. * Mt * (-3. * Mc4 + Mb2 * (3. * Mc2 + Pw2))) / (Mc2 * Sqq);
        result.element({9, 7, 7, 2}) = (16. * Mb2 * Mt * Pw2) / Sqq;
        result.element({9, 7, 7, 4}) = (-2. * Mb2 * Mt * Pw2 * (Mb2 + 3. * Mc2 - Sqq)) / (Mc2 * Sqq);
        result.element({9, 7, 9, 5}) = (8. * Mb2 * Pw2 * (Mb2 + 3. * Mc2 - Sqq) * (2. * Mt2 + Sqq)) / (3. * Mc2 * Sqq);
        result.element({9, 7, 9, 6}) =
            (8. * (-3. * Mc4 + Mb2 * (3. * Mc2 + Pw2)) * (2. * Mt2 + Sqq)) / (3. * Mc2 * Sqq);
        result.element({9, 7, 9, 7}) =
            (8. * (2. * Mt2 + Sqq) *
             (3. * Mc2 * (2. * Mc2 - Sqq) * Sqq + Mb2 * (Pw2 * Sqq + 2. * Mc2 * (8. * Pw2 + 3. * Sqq)))) /
            (3. * Mc2 * Sqq2);
        result.element({10, 5, 6, 1}) = (2. * Mb2 * Mt * Pw2 * (-Mb2 + Mc2 + Sqq)) / (Mc2 * Sqq);
        result.element({10, 5, 6, 4}) = (-8. * Mb4 * Mt * Pw4) / (Mc2 * Sqq);
        result.element({10, 5, 8, 1}) = (-2. * Mb2 * Mt * Pw2 * (-Mb2 + Mc2 + Sqq)) / (Mc2 * Sqq);
        result.element({10, 5, 8, 4}) = (8. * Mb4 * Mt * Pw4) / (Mc2 * Sqq);
        result.element({10, 5, 10, 5}) = (32. * Mb4 * Pw4 * (2. * Mt2 + Sqq)) / (3. * Mc2 * Sqq);
        result.element({10, 5, 10, 6}) = (8. * Mb2 * Pw2 * (Mb2 - Mc2 - Sqq) * (2. * Mt2 + Sqq)) / (3. * Mc2 * Sqq);
        result.element({10, 5, 10, 7}) =
            (8. * Mb2 * Pw2 * (Mb2 + 3. * Mc2 - Sqq) * (2. * Mt2 + Sqq)) / (3. * Mc2 * Sqq);
        result.element({10, 6, 6, 1}) = -2. * Mt * (3 + (Mb2 * Pw2) / (Mc2 * Sqq));
        result.element({10, 6, 6, 4}) = (-2. * Mb2 * Mt * Pw2 * (Mb2 - Mc2 - Sqq)) / (Mc2 * Sqq);
        result.element({10, 6, 8, 1}) = 2. * Mt * (3 + (Mb2 * Pw2) / (Mc2 * Sqq));
        result.element({10, 6, 8, 4}) = (2. * Mb2 * Mt * Pw2 * (Mb2 - Mc2 - Sqq)) / (Mc2 * Sqq);
        result.element({10, 6, 10, 5}) = (8. * Mb2 * Pw2 * (Mb2 - Mc2 - Sqq) * (2. * Mt2 + Sqq)) / (3. * Mc2 * Sqq);
        result.element({10, 6, 10, 6}) = (8. * (2. * Mt2 + Sqq) * (Mb2 * Pw2 + 3. * Mc2 * Sqq)) / (3. * Mc2 * Sqq);
        result.element({10, 6, 10, 7}) =
            (8. * (-3. * Mc4 + Mb2 * (3. * Mc2 + Pw2)) * (2. * Mt2 + Sqq)) / (3. * Mc2 * Sqq);
        result.element({10, 7, 6, 1}) = (-2. * Mt * (-3. * Mc4 + Mb2 * (3. * Mc2 + Pw2))) / (Mc2 * Sqq);
        result.element({10, 7, 6, 2}) = (16. * Mb2 * Mt * Pw2) / Sqq;
        result.element({10, 7, 6, 4}) = (-2. * Mb2 * Mt * Pw2 * (Mb2 + 3. * Mc2 - Sqq)) / (Mc2 * Sqq);
        result.element({10, 7, 8, 1}) = (2. * Mt * (-3. * Mc4 + Mb2 * (3. * Mc2 + Pw2))) / (Mc2 * Sqq);
        result.element({10, 7, 8, 2}) = (16. * Mb2 * Mt * Pw2) / Sqq;
        result.element({10, 7, 8, 4}) = (2. * Mb2 * Mt * Pw2 * (Mb2 + 3. * Mc2 - Sqq)) / (Mc2 * Sqq);
        result.element({10, 7, 10, 5}) =
            (8. * Mb2 * Pw2 * (Mb2 + 3. * Mc2 - Sqq) * (2. * Mt2 + Sqq)) / (3. * Mc2 * Sqq);
        result.element({10, 7, 10, 6}) =
            (8. * (-3. * Mc4 + Mb2 * (3. * Mc2 + Pw2)) * (2. * Mt2 + Sqq)) / (3. * Mc2 * Sqq);
        result.element({10, 7, 10, 7}) =
            (8. * (2. * Mt2 + Sqq) *
             (3. * Mc2 * (2. * Mc2 - Sqq) * Sqq + Mb2 * (Pw2 * Sqq + 2. * Mc2 * (8. * Pw2 + 3. * Sqq)))) /
            (3. * Mc2 * Sqq2);


        result *= RateNorm;

        return result;
    }

} // namespace Hammer
