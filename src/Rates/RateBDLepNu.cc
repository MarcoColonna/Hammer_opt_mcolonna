///
/// @file  RateBDLepNu.cc
/// @brief \f$ B \rightarrow D \tau\nu \f$ total rate
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/Rates/RateBDLepNu.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    RateBDLepNu::RateBDLepNu() {
        // Create tensor rank and dimensions
        IndexList dims{{11, 4, 11, 4, _nPoints}};
        string name{"RateBDLepNuQ2"};
        auto& pdg = PID::instance();
        addProcessSignature(PID::BPLUS, {-PID::D0, PID::NU_TAU, PID::ANTITAU});
        addIntegrationBoundaries(
            {PS::makeQ2Function(pdg.getMass(PID::ANTITAU), pdg.getMass(PID::BPLUS) - pdg.getMass(PID::D0))});
        addTensor(Tensor{
            name, MD::makeEmptySparse(dims, {WILSON_BCTAUNU, FF_BD, WILSON_BCTAUNU_HC, FF_BD_HC, INTEGRATION_INDEX})});

        addProcessSignature(PID::BZERO, {PID::DMINUS, PID::NU_TAU, PID::ANTITAU});
        addIntegrationBoundaries(
            {PS::makeQ2Function(pdg.getMass(PID::ANTITAU), pdg.getMass(PID::BZERO) - pdg.getMass(PID::DMINUS))});
        addTensor(Tensor{
            name, MD::makeEmptySparse(dims, {WILSON_BCTAUNU, FF_BD, WILSON_BCTAUNU_HC, FF_BD_HC, INTEGRATION_INDEX})});

        addProcessSignature(PID::BPLUS, {-PID::D0, PID::NU_MU, PID::ANTIMUON});
        addIntegrationBoundaries(
            {PS::makeQ2Function(pdg.getMass(PID::ANTIMUON), pdg.getMass(PID::BPLUS) - pdg.getMass(PID::D0))});
        addTensor(Tensor{
            name, MD::makeEmptySparse(dims, {WILSON_BCMUNU, FF_BD, WILSON_BCMUNU_HC, FF_BD_HC, INTEGRATION_INDEX})});

        addProcessSignature(PID::BZERO, {PID::DMINUS, PID::NU_MU, PID::ANTIMUON});
        addIntegrationBoundaries(
            {PS::makeQ2Function(pdg.getMass(PID::ANTIMUON), pdg.getMass(PID::BZERO) - pdg.getMass(PID::DMINUS))});
        addTensor(Tensor{
            name, MD::makeEmptySparse(dims, {WILSON_BCMUNU, FF_BD, WILSON_BCMUNU_HC, FF_BD_HC, INTEGRATION_INDEX})});

        addProcessSignature(PID::BPLUS, {-PID::D0, PID::NU_E, PID::POSITRON});
        addIntegrationBoundaries(
            {PS::makeQ2Function(pdg.getMass(PID::POSITRON), pdg.getMass(PID::BPLUS) - pdg.getMass(PID::D0))});
        addTensor(Tensor{
            name, MD::makeEmptySparse(dims, {WILSON_BCENU, FF_BD, WILSON_BCENU_HC, FF_BD_HC, INTEGRATION_INDEX})});

        addProcessSignature(PID::BZERO, {PID::DMINUS, PID::NU_E, PID::POSITRON});
        addIntegrationBoundaries(
            {PS::makeQ2Function(pdg.getMass(PID::POSITRON), pdg.getMass(PID::BZERO) - pdg.getMass(PID::DMINUS))});
        addTensor(Tensor{
            name, MD::makeEmptySparse(dims, {WILSON_BCENU, FF_BD, WILSON_BCENU_HC, FF_BD_HC, INTEGRATION_INDEX})});

        // bs -> cs
        addProcessSignature(PID::BS, {PID::DSMINUS, PID::NU_TAU, PID::ANTITAU});
        addIntegrationBoundaries(
            {PS::makeQ2Function(pdg.getMass(PID::ANTITAU), pdg.getMass(PID::BS) - pdg.getMass(PID::DSMINUS))});
        addTensor(Tensor{name, MD::makeEmptySparse(
                                   dims, {WILSON_BCTAUNU, FF_BSDS, WILSON_BCTAUNU_HC, FF_BSDS_HC, INTEGRATION_INDEX})});

        addProcessSignature(PID::BS, {PID::DSMINUS, PID::NU_MU, PID::ANTIMUON});
        addIntegrationBoundaries(
            {PS::makeQ2Function(pdg.getMass(PID::ANTIMUON), pdg.getMass(PID::BS) - pdg.getMass(PID::DSMINUS))});
        addTensor(Tensor{name, MD::makeEmptySparse(
                                   dims, {WILSON_BCMUNU, FF_BSDS, WILSON_BCMUNU_HC, FF_BSDS_HC, INTEGRATION_INDEX})});

        addProcessSignature(PID::BS, {PID::DSMINUS, PID::NU_E, PID::POSITRON});
        addIntegrationBoundaries(
            {PS::makeQ2Function(pdg.getMass(PID::POSITRON), pdg.getMass(PID::BS) - pdg.getMass(PID::DSMINUS))});
        addTensor(Tensor{
            name, MD::makeEmptySparse(dims, {WILSON_BCENU, FF_BSDS, WILSON_BCENU_HC, FF_BSDS_HC, INTEGRATION_INDEX})});

        // b -> u
        addProcessSignature(PID::BPLUS, {PID::PI0, PID::NU_TAU, PID::ANTITAU});
        addIntegrationBoundaries(
            {PS::makeQ2Function(pdg.getMass(PID::ANTITAU), pdg.getMass(PID::BPLUS) - pdg.getMass(PID::PI0))});
        addTensor(Tensor{name, MD::makeEmptySparse(
                                   dims, {WILSON_BUTAUNU, FF_BPI, WILSON_BUTAUNU_HC, FF_BPI_HC, INTEGRATION_INDEX})});

        addProcessSignature(PID::BZERO, {PID::PIMINUS, PID::NU_TAU, PID::ANTITAU});
        addIntegrationBoundaries(
            {PS::makeQ2Function(pdg.getMass(PID::ANTITAU), pdg.getMass(PID::BZERO) - pdg.getMass(PID::PIMINUS))});
        addTensor(Tensor{name, MD::makeEmptySparse(
                                   dims, {WILSON_BUTAUNU, FF_BPI, WILSON_BUTAUNU_HC, FF_BPI_HC, INTEGRATION_INDEX})});

        addProcessSignature(PID::BPLUS, {PID::PI0, PID::NU_MU, PID::ANTIMUON});
        addIntegrationBoundaries(
            {PS::makeQ2Function(pdg.getMass(PID::ANTIMUON), pdg.getMass(PID::BPLUS) - pdg.getMass(PID::PI0))});
        addTensor(Tensor{
            name, MD::makeEmptySparse(dims, {WILSON_BUMUNU, FF_BPI, WILSON_BUMUNU_HC, FF_BPI_HC, INTEGRATION_INDEX})});

        addProcessSignature(PID::BZERO, {PID::PIMINUS, PID::NU_MU, PID::ANTIMUON});
        addIntegrationBoundaries(
            {PS::makeQ2Function(pdg.getMass(PID::ANTIMUON), pdg.getMass(PID::BZERO) - pdg.getMass(PID::PIMINUS))});
        addTensor(Tensor{
            name, MD::makeEmptySparse(dims, {WILSON_BUMUNU, FF_BPI, WILSON_BUMUNU_HC, FF_BPI_HC, INTEGRATION_INDEX})});

        addProcessSignature(PID::BPLUS, {PID::PI0, PID::NU_E, PID::POSITRON});
        addIntegrationBoundaries(
            {PS::makeQ2Function(pdg.getMass(PID::POSITRON), pdg.getMass(PID::BPLUS) - pdg.getMass(PID::PI0))});
        addTensor(Tensor{
            name, MD::makeEmptySparse(dims, {WILSON_BUENU, FF_BPI, WILSON_BUENU_HC, FF_BPI_HC, INTEGRATION_INDEX})});

        addProcessSignature(PID::BZERO, {PID::PIMINUS, PID::NU_E, PID::POSITRON});
        addIntegrationBoundaries(
            {PS::makeQ2Function(pdg.getMass(PID::POSITRON), pdg.getMass(PID::BZERO) - pdg.getMass(PID::PIMINUS))});
        addTensor(Tensor{
            name, MD::makeEmptySparse(dims, {WILSON_BUENU, FF_BPI, WILSON_BUENU_HC, FF_BPI_HC, INTEGRATION_INDEX})});

        // bs -> us
        addProcessSignature(PID::BS, {PID::KMINUS, PID::NU_TAU, PID::ANTITAU});
        addIntegrationBoundaries(
            {PS::makeQ2Function(pdg.getMass(PID::ANTITAU), pdg.getMass(PID::BS) - pdg.getMass(PID::KMINUS))});
        addTensor(Tensor{name, MD::makeEmptySparse(
                                   dims, {WILSON_BUTAUNU, FF_BSK, WILSON_BUTAUNU_HC, FF_BSK_HC, INTEGRATION_INDEX})});

        addProcessSignature(PID::BS, {PID::KMINUS, PID::NU_MU, PID::ANTIMUON});
        addIntegrationBoundaries(
            {PS::makeQ2Function(pdg.getMass(PID::ANTIMUON), pdg.getMass(PID::BS) - pdg.getMass(PID::KMINUS))});
        addTensor(Tensor{
            name, MD::makeEmptySparse(dims, {WILSON_BUMUNU, FF_BSK, WILSON_BUMUNU_HC, FF_BSK_HC, INTEGRATION_INDEX})});

        addProcessSignature(PID::BS, {PID::KMINUS, PID::NU_E, PID::POSITRON});
        addIntegrationBoundaries(
            {PS::makeQ2Function(pdg.getMass(PID::POSITRON), pdg.getMass(PID::BS) - pdg.getMass(PID::KMINUS))});
        addTensor(Tensor{
            name, MD::makeEmptySparse(dims, {WILSON_BUENU, FF_BSK, WILSON_BUENU_HC, FF_BSK_HC, INTEGRATION_INDEX})});

        setSignatureIndex();
    }

    Tensor RateBDLepNu::evalAtPSPoint(const vector<double>& point) {
        auto labs = getTensor().labels();
        labs.pop_back();
        auto dimensions = getTensor().dims();
        dimensions.pop_back();
        Tensor result{"RateBDLepNu", MD::makeEmptySparse(dimensions, labs)};

        const double Mb = masses()[0];
        const double Mc = masses()[1];
        const double Mt = masses()[3];

        const double Mb2 = Mb * Mb;
        const double Mc2 = Mc * Mc;
        const double Mt2 = Mt * Mt;

        // kinematic objects
        const double Sqq = point[0];
        const double Sqq2 = pow(Sqq, 2.);
        const double Ew = (Mb2 - Mc2 + Sqq) / (2 * Mb);
        const double Ew2 = Ew * Ew;
        const double Pw = sqrt(Ew2 - Sqq);
        const double Pw2 = Pw * Pw;
        const double Mb2Mc2Sq = pow(Mb2 - Mc2, 2.);

        // combinatoric for neutral unflavored mesons
        const double Cv2 = (getPdgIds().second[0] == PID::PI0) ? 0.5 : 1.;

        double RateNorm = Cv2 * (pow(GFermi, 2.) * Pw * pow(-Mt2 + Sqq, 2.)) / (32. * Mb2 * pi3 * Sqq);

        // set non-zero tensor elements
        result.element({0, 1, 0, 1}) = (Mb2Mc2Sq * Mt2) / (2. * Sqq2);
        result.element({0, 1, 1, 0}) = ((-Mb2 + Mc2) * Mt) / (2. * Sqq);
        result.element({0, 1, 3, 0}) = ((-Mb2 + Mc2) * Mt) / (2. * Sqq);
        result.element({0, 1, 5, 1}) = (Mb2Mc2Sq * Mt2) / (2. * Sqq2);
        result.element({0, 1, 7, 1}) = (Mb2Mc2Sq * Mt2) / (2. * Sqq2);
        result.element({0, 2, 0, 2}) = (2. * Mb2 * Pw2 * (Mt2 + 2. * Sqq)) / (3. * Sqq2);
        result.element({0, 2, 5, 2}) = (2. * Mb2 * Pw2 * (Mt2 + 2. * Sqq)) / (3. * Sqq2);
        result.element({0, 2, 7, 2}) = (2. * Mb2 * Pw2 * (Mt2 + 2. * Sqq)) / (3. * Sqq2);
        result.element({0, 2, 9, 3}) = (8. * Mb2 * Mt * (-Ew2 + Sqq)) / Sqq;
        result.element({1, 0, 0, 1}) = ((-Mb2 + Mc2) * Mt) / (2. * Sqq);
        result.element({1, 0, 1, 0}) = 0.5;
        result.element({1, 0, 3, 0}) = 0.5;
        result.element({1, 0, 5, 1}) = ((-Mb2 + Mc2) * Mt) / (2. * Sqq);
        result.element({1, 0, 7, 1}) = ((-Mb2 + Mc2) * Mt) / (2. * Sqq);
        result.element({2, 0, 2, 0}) = 0.5;
        result.element({2, 0, 4, 0}) = 0.5;
        result.element({2, 0, 6, 1}) = ((-Mb2 + Mc2) * Mt) / (2. * Sqq);
        result.element({2, 0, 8, 1}) = ((-Mb2 + Mc2) * Mt) / (2. * Sqq);
        result.element({3, 0, 0, 1}) = ((-Mb2 + Mc2) * Mt) / (2. * Sqq);
        result.element({3, 0, 1, 0}) = 0.5;
        result.element({3, 0, 3, 0}) = 0.5;
        result.element({3, 0, 5, 1}) = ((-Mb2 + Mc2) * Mt) / (2. * Sqq);
        result.element({3, 0, 7, 1}) = ((-Mb2 + Mc2) * Mt) / (2. * Sqq);
        result.element({4, 0, 2, 0}) = 0.5;
        result.element({4, 0, 4, 0}) = 0.5;
        result.element({4, 0, 6, 1}) = ((-Mb2 + Mc2) * Mt) / (2. * Sqq);
        result.element({4, 0, 8, 1}) = ((-Mb2 + Mc2) * Mt) / (2. * Sqq);
        result.element({5, 1, 0, 1}) = (Mb2Mc2Sq * Mt2) / (2. * Sqq2);
        result.element({5, 1, 1, 0}) = ((-Mb2 + Mc2) * Mt) / (2. * Sqq);
        result.element({5, 1, 3, 0}) = ((-Mb2 + Mc2) * Mt) / (2. * Sqq);
        result.element({5, 1, 5, 1}) = (Mb2Mc2Sq * Mt2) / (2. * Sqq2);
        result.element({5, 1, 7, 1}) = (Mb2Mc2Sq * Mt2) / (2. * Sqq2);
        result.element({5, 2, 0, 2}) = (2. * Mb2 * Pw2 * (Mt2 + 2. * Sqq)) / (3. * Sqq2);
        result.element({5, 2, 5, 2}) = (2. * Mb2 * Pw2 * (Mt2 + 2. * Sqq)) / (3. * Sqq2);
        result.element({5, 2, 7, 2}) = (2. * Mb2 * Pw2 * (Mt2 + 2. * Sqq)) / (3. * Sqq2);
        result.element({5, 2, 9, 3}) = (8. * Mb2 * Mt * (-Ew2 + Sqq)) / Sqq;
        result.element({6, 1, 2, 0}) = ((-Mb2 + Mc2) * Mt) / (2. * Sqq);
        result.element({6, 1, 4, 0}) = ((-Mb2 + Mc2) * Mt) / (2. * Sqq);
        result.element({6, 1, 6, 1}) = (Mb2Mc2Sq * Mt2) / (2. * Sqq2);
        result.element({6, 1, 8, 1}) = (Mb2Mc2Sq * Mt2) / (2. * Sqq2);
        result.element({6, 2, 6, 2}) = (2. * Mb2 * Pw2 * (Mt2 + 2. * Sqq)) / (3. * Sqq2);
        result.element({6, 2, 8, 2}) = (2. * Mb2 * Pw2 * (Mt2 + 2. * Sqq)) / (3. * Sqq2);
        result.element({6, 2, 10, 3}) = (8. * Mb2 * Mt * (-Ew2 + Sqq)) / Sqq;
        result.element({7, 1, 0, 1}) = (Mb2Mc2Sq * Mt2) / (2. * Sqq2);
        result.element({7, 1, 1, 0}) = ((-Mb2 + Mc2) * Mt) / (2. * Sqq);
        result.element({7, 1, 3, 0}) = ((-Mb2 + Mc2) * Mt) / (2. * Sqq);
        result.element({7, 1, 5, 1}) = (Mb2Mc2Sq * Mt2) / (2. * Sqq2);
        result.element({7, 1, 7, 1}) = (Mb2Mc2Sq * Mt2) / (2. * Sqq2);
        result.element({7, 2, 0, 2}) = (2. * Mb2 * Pw2 * (Mt2 + 2. * Sqq)) / (3. * Sqq2);
        result.element({7, 2, 5, 2}) = (2. * Mb2 * Pw2 * (Mt2 + 2. * Sqq)) / (3. * Sqq2);
        result.element({7, 2, 7, 2}) = (2. * Mb2 * Pw2 * (Mt2 + 2. * Sqq)) / (3. * Sqq2);
        result.element({7, 2, 9, 3}) = (8. * Mb2 * Mt * (-Ew2 + Sqq)) / Sqq;
        result.element({8, 1, 2, 0}) = ((-Mb2 + Mc2) * Mt) / (2. * Sqq);
        result.element({8, 1, 4, 0}) = ((-Mb2 + Mc2) * Mt) / (2. * Sqq);
        result.element({8, 1, 6, 1}) = (Mb2Mc2Sq * Mt2) / (2. * Sqq2);
        result.element({8, 1, 8, 1}) = (Mb2Mc2Sq * Mt2) / (2. * Sqq2);
        result.element({8, 2, 6, 2}) = (2. * Mb2 * Pw2 * (Mt2 + 2. * Sqq)) / (3. * Sqq2);
        result.element({8, 2, 8, 2}) = (2. * Mb2 * Pw2 * (Mt2 + 2. * Sqq)) / (3. * Sqq2);
        result.element({8, 2, 10, 3}) = (8. * Mb2 * Mt * (-Ew2 + Sqq)) / Sqq;
        result.element({9, 3, 0, 2}) = (8. * Mb2 * Mt * (-Ew2 + Sqq)) / Sqq;
        result.element({9, 3, 5, 2}) = (8. * Mb2 * Mt * (-Ew2 + Sqq)) / Sqq;
        result.element({9, 3, 7, 2}) = (8. * Mb2 * Mt * (-Ew2 + Sqq)) / Sqq;
        result.element({9, 3, 9, 3}) = (32. * Mb2 * Pw2 * (2. * Mt2 + Sqq)) / (3. * Sqq);
        result.element({10, 3, 6, 2}) = (8. * Mb2 * Mt * (-Ew2 + Sqq)) / Sqq;
        result.element({10, 3, 8, 2}) = (8. * Mb2 * Mt * (-Ew2 + Sqq)) / Sqq;
        result.element({10, 3, 10, 3}) = (32. * Mb2 * Pw2 * (2. * Mt2 + Sqq)) / (3. * Sqq);

        result *= RateNorm;

        return result;
    }

} // namespace Hammer
