///
/// @file  RateBRhoOmegaLepNu.cc
/// @brief \f$ B \rightarrow \rho \tau\nu \f$ total rate
/// @brief \f$ B \rightarrow \omega \tau\nu \f$ total rate
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/Rates/RateBRhoOmegaLepNu.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"
#include "Hammer/Tools/Pdg.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    RateBRhoOmegaLepNu::RateBRhoOmegaLepNu() {
        // Create tensor rank and dimensions
        IndexList dims{{11, 8, 11, 8, _nPoints}};
        string name{"RateBRhoOmegaLepNuQ2"};
        auto& pdg = PID::instance();
        // rho
        addProcessSignature(PID::BPLUS, {PID::RHO0, PID::NU_TAU, PID::ANTITAU});
        addIntegrationBoundaries(
            {PS::makeQ2Function(pdg.getMass(PID::ANTITAU), pdg.getMass(PID::BPLUS) - pdg.getMass(PID::RHO0))});
        addTensor(Tensor{name, MD::makeEmptySparse(
                                   dims, {WILSON_BUTAUNU, FF_BRHO, WILSON_BUTAUNU_HC, FF_BRHO_HC, INTEGRATION_INDEX})});

        addProcessSignature(PID::BZERO, {PID::RHOMINUS, PID::NU_TAU, PID::ANTITAU});
        addIntegrationBoundaries(
            {PS::makeQ2Function(pdg.getMass(PID::ANTITAU), pdg.getMass(PID::BZERO) - pdg.getMass(PID::RHOMINUS))});
        addTensor(Tensor{name, MD::makeEmptySparse(
                                   dims, {WILSON_BUTAUNU, FF_BRHO, WILSON_BUTAUNU_HC, FF_BRHO_HC, INTEGRATION_INDEX})});

        addProcessSignature(PID::BPLUS, {PID::RHO0, PID::NU_MU, PID::ANTIMUON});
        addIntegrationBoundaries(
            {PS::makeQ2Function(pdg.getMass(PID::ANTIMUON), pdg.getMass(PID::BPLUS) - pdg.getMass(PID::RHO0))});
        addTensor(Tensor{name, MD::makeEmptySparse(
                                   dims, {WILSON_BUMUNU, FF_BRHO, WILSON_BUMUNU_HC, FF_BRHO_HC, INTEGRATION_INDEX})});

        addProcessSignature(PID::BZERO, {PID::RHOMINUS, PID::NU_MU, PID::ANTIMUON});
        addIntegrationBoundaries(
            {PS::makeQ2Function(pdg.getMass(PID::ANTIMUON), pdg.getMass(PID::BZERO) - pdg.getMass(PID::RHOMINUS))});
        addTensor(Tensor{name, MD::makeEmptySparse(
                                   dims, {WILSON_BUMUNU, FF_BRHO, WILSON_BUMUNU_HC, FF_BRHO_HC, INTEGRATION_INDEX})});

        addProcessSignature(PID::BPLUS, {PID::RHO0, PID::NU_E, PID::POSITRON});
        addIntegrationBoundaries(
            {PS::makeQ2Function(pdg.getMass(PID::POSITRON), pdg.getMass(PID::BPLUS) - pdg.getMass(PID::RHO0))});
        addTensor(Tensor{
            name, MD::makeEmptySparse(dims, {WILSON_BUENU, FF_BRHO, WILSON_BUENU_HC, FF_BRHO_HC, INTEGRATION_INDEX})});

        addProcessSignature(PID::BZERO, {PID::RHOMINUS, PID::NU_E, PID::POSITRON});
        addIntegrationBoundaries(
            {PS::makeQ2Function(pdg.getMass(PID::POSITRON), pdg.getMass(PID::BZERO) - pdg.getMass(PID::RHOMINUS))});
        addTensor(Tensor{
            name, MD::makeEmptySparse(dims, {WILSON_BUENU, FF_BRHO, WILSON_BUENU_HC, FF_BRHO_HC, INTEGRATION_INDEX})});

        // omega
        addProcessSignature(PID::BPLUS, {PID::OMEGA, PID::NU_TAU, PID::ANTITAU});
        addIntegrationBoundaries(
            {PS::makeQ2Function(pdg.getMass(PID::ANTITAU), pdg.getMass(PID::BPLUS) - pdg.getMass(PID::OMEGA))});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BUTAUNU, FF_BOMEGA, WILSON_BUTAUNU_HC, FF_BOMEGA_HC,
                                                          INTEGRATION_INDEX})});

        addProcessSignature(PID::BPLUS, {PID::OMEGA, PID::NU_MU, PID::ANTIMUON});
        addIntegrationBoundaries(
            {PS::makeQ2Function(pdg.getMass(PID::ANTIMUON), pdg.getMass(PID::BPLUS) - pdg.getMass(PID::OMEGA))});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {WILSON_BUMUNU, FF_BOMEGA, WILSON_BUMUNU_HC, FF_BOMEGA_HC,
                                                          INTEGRATION_INDEX})});

        addProcessSignature(PID::BPLUS, {PID::OMEGA, PID::NU_E, PID::POSITRON});
        addIntegrationBoundaries(
            {PS::makeQ2Function(pdg.getMass(PID::POSITRON), pdg.getMass(PID::BPLUS) - pdg.getMass(PID::OMEGA))});
        addTensor(Tensor{name, MD::makeEmptySparse(
                                   dims, {WILSON_BUENU, FF_BOMEGA, WILSON_BUENU_HC, FF_BOMEGA_HC, INTEGRATION_INDEX})});

        setSignatureIndex();
    }

    Tensor RateBRhoOmegaLepNu::evalAtPSPoint(const vector<double>& point) {
        auto labs = getTensor().labels();
        labs.pop_back();
        auto dimensions = getTensor().dims();
        dimensions.pop_back();
        Tensor result{"RateBRhoOmegaLepNu", MD::makeEmptySparse(dimensions, labs)};

        const double Mb = masses()[0];
        const double Mu = masses()[1];
        const double Mt = masses()[3];

        const double Mb2 = Mb * Mb;
        const double Mu2 = Mu * Mu;
        const double Mt2 = Mt * Mt;
        const double MbMu2 = (Mb + Mu) * (Mb + Mu);

        // kinematic objects
        const double Sqq = point[0];
        const double Sqq2 = Sqq * Sqq;
        const double Ew = (Mb2 - Mu2 + Sqq) / (2 * Mb);
        const double Pw = sqrt(Ew * Ew - Sqq);
        const double Pw2 = Pw * Pw;

        // combinatoric for neutral unflavored mesons
        const double Cv2 = ((getPdgIds().second[0] == PID::RHO0) || (getPdgIds().second[0] == PID::OMEGA)) ? 0.5 : 1.;

        double RateNorm = Cv2 * pow(GFermi, 2.) * Pw * pow(-Mt2 + Sqq, 2.) / (32. * Mb2 * pi3 * Sqq);

        //Fix sign convention in T2 and T23 definitions, FF indices 6 and 7.
        const double signfix = -1.;
        // set non-zero tensor elements
        result.element({0, 1, 0, 1}) = (4. * Mb2 * Pw2 * (Mt2 + 2. * Sqq)) / (3. * MbMu2 * Sqq);
        result.element({0, 1, 5, 1}) = (4. * Mb2 * Pw2 * (Mt2 + 2. * Sqq)) / (3. * MbMu2 * Sqq);
        result.element({0, 1, 7, 1}) = (4. * Mb2 * Pw2 * (Mt2 + 2. * Sqq)) / (3. * MbMu2 * Sqq);
        result.element({0, 1, 9, 5}) = (-16. * Mb2 * Mt * Pw2) / ((Mb + Mu) * Sqq);
        result.element({0, 2, 0, 2}) = (2. * Mb2 * Mt2 * Pw2) / Sqq2;
        result.element({0, 2, 1, 0}) = -((Mb2 * Mt * Pw2) / (Mu * Sqq));
        result.element({0, 2, 3, 0}) = (Mb2 * Mt * Pw2) / (Mu * Sqq);
        result.element({0, 2, 5, 2}) = (-2. * Mb2 * Mt2 * Pw2) / Sqq2;
        result.element({0, 2, 7, 2}) = (2. * Mb2 * Mt2 * Pw2) / Sqq2;
        result.element({0, 3, 0, 3}) = (MbMu2 * (Mt2 + 2. * Sqq)) / (3. * Sqq);
        result.element({0, 3, 5, 3}) = -(MbMu2 * (Mt2 + 2. * Sqq)) / (3. * Sqq);
        result.element({0, 3, 7, 3}) = (MbMu2 * (Mt2 + 2. * Sqq)) / (3. * Sqq);
        result.element({0, 3, 9, 6}) = signfix * (-4. * Mt * (Mb - Mu) * MbMu2) / Sqq;
        result.element({0, 4, 0, 4}) = (32. * Mb2 * Mu2 * (Mt2 + 2. * Sqq)) / (3. * Sqq2);
        result.element({0, 4, 5, 4}) = (-32. * Mb2 * Mu2 * (Mt2 + 2. * Sqq)) / (3. * Sqq2);
        result.element({0, 4, 7, 4}) = (32. * Mb2 * Mu2 * (Mt2 + 2. * Sqq)) / (3. * Sqq2);
        result.element({0, 4, 9, 7}) = signfix * (-64. * Mb2 * Mt * Mu2) / ((Mb + Mu) * Sqq);
        result.element({1, 0, 0, 2}) = -((Mb2 * Mt * Pw2) / (Mu * Sqq));
        result.element({1, 0, 1, 0}) = (Mb2 * Pw2) / (2. * Mu2);
        result.element({1, 0, 3, 0}) = -(Mb2 * Pw2) / (2. * Mu2);
        result.element({1, 0, 5, 2}) = (Mb2 * Mt * Pw2) / (Mu * Sqq);
        result.element({1, 0, 7, 2}) = -((Mb2 * Mt * Pw2) / (Mu * Sqq));
        result.element({2, 0, 2, 0}) = (Mb2 * Pw2) / (2. * Mu2);
        result.element({2, 0, 4, 0}) = -(Mb2 * Pw2) / (2. * Mu2);
        result.element({2, 0, 6, 2}) = (Mb2 * Mt * Pw2) / (Mu * Sqq);
        result.element({2, 0, 8, 2}) = -((Mb2 * Mt * Pw2) / (Mu * Sqq));
        result.element({3, 0, 0, 2}) = (Mb2 * Mt * Pw2) / (Mu * Sqq);
        result.element({3, 0, 1, 0}) = -(Mb2 * Pw2) / (2. * Mu2);
        result.element({3, 0, 3, 0}) = (Mb2 * Pw2) / (2. * Mu2);
        result.element({3, 0, 5, 2}) = -((Mb2 * Mt * Pw2) / (Mu * Sqq));
        result.element({3, 0, 7, 2}) = (Mb2 * Mt * Pw2) / (Mu * Sqq);
        result.element({4, 0, 2, 0}) = -(Mb2 * Pw2) / (2. * Mu2);
        result.element({4, 0, 4, 0}) = (Mb2 * Pw2) / (2. * Mu2);
        result.element({4, 0, 6, 2}) = -((Mb2 * Mt * Pw2) / (Mu * Sqq));
        result.element({4, 0, 8, 2}) = (Mb2 * Mt * Pw2) / (Mu * Sqq);
        result.element({5, 1, 0, 1}) = (4. * Mb2 * Pw2 * (Mt2 + 2. * Sqq)) / (3. * MbMu2 * Sqq);
        result.element({5, 1, 5, 1}) = (4. * Mb2 * Pw2 * (Mt2 + 2. * Sqq)) / (3. * MbMu2 * Sqq);
        result.element({5, 1, 7, 1}) = (4. * Mb2 * Pw2 * (Mt2 + 2. * Sqq)) / (3. * MbMu2 * Sqq);
        result.element({5, 1, 9, 5}) = (-16. * Mb2 * Mt * Pw2) / ((Mb + Mu) * Sqq);
        result.element({5, 2, 0, 2}) = (-2. * Mb2 * Mt2 * Pw2) / Sqq2;
        result.element({5, 2, 1, 0}) = (Mb2 * Mt * Pw2) / (Mu * Sqq);
        result.element({5, 2, 3, 0}) = -((Mb2 * Mt * Pw2) / (Mu * Sqq));
        result.element({5, 2, 5, 2}) = (2. * Mb2 * Mt2 * Pw2) / Sqq2;
        result.element({5, 2, 7, 2}) = (-2. * Mb2 * Mt2 * Pw2) / Sqq2;
        result.element({5, 3, 0, 3}) = -(MbMu2 * (Mt2 + 2. * Sqq)) / (3. * Sqq);
        result.element({5, 3, 5, 3}) = (MbMu2 * (Mt2 + 2. * Sqq)) / (3. * Sqq);
        result.element({5, 3, 7, 3}) = -(MbMu2 * (Mt2 + 2. * Sqq)) / (3. * Sqq);
        result.element({5, 3, 9, 6}) = signfix * (4. * Mt * (Mb - Mu) * MbMu2) / Sqq;
        result.element({5, 4, 0, 4}) = (-32. * Mb2 * Mu2 * (Mt2 + 2. * Sqq)) / (3. * Sqq2);
        result.element({5, 4, 5, 4}) = (32. * Mb2 * Mu2 * (Mt2 + 2. * Sqq)) / (3. * Sqq2);
        result.element({5, 4, 7, 4}) = (-32. * Mb2 * Mu2 * (Mt2 + 2. * Sqq)) / (3. * Sqq2);
        result.element({5, 4, 9, 7}) = signfix * (64. * Mb2 * Mt * Mu2) / ((Mb + Mu) * Sqq);
        result.element({6, 1, 6, 1}) = (4. * Mb2 * Pw2 * (Mt2 + 2. * Sqq)) / (3. * MbMu2 * Sqq);
        result.element({6, 1, 8, 1}) = (4. * Mb2 * Pw2 * (Mt2 + 2. * Sqq)) / (3. * MbMu2 * Sqq);
        result.element({6, 1, 10, 5}) = (-16. * Mb2 * Mt * Pw2) / ((Mb + Mu) * Sqq);
        result.element({6, 2, 2, 0}) = (Mb2 * Mt * Pw2) / (Mu * Sqq);
        result.element({6, 2, 4, 0}) = -((Mb2 * Mt * Pw2) / (Mu * Sqq));
        result.element({6, 2, 6, 2}) = (2. * Mb2 * Mt2 * Pw2) / Sqq2;
        result.element({6, 2, 8, 2}) = (-2. * Mb2 * Mt2 * Pw2) / Sqq2;
        result.element({6, 3, 6, 3}) = (MbMu2 * (Mt2 + 2. * Sqq)) / (3. * Sqq);
        result.element({6, 3, 8, 3}) = -(MbMu2 * (Mt2 + 2. * Sqq)) / (3. * Sqq);
        result.element({6, 3, 10, 6}) = signfix * (-4. * Mt * (Mb - Mu) * MbMu2) / Sqq;
        result.element({6, 4, 6, 4}) = (32. * Mb2 * Mu2 * (Mt2 + 2. * Sqq)) / (3. * Sqq2);
        result.element({6, 4, 8, 4}) = (-32. * Mb2 * Mu2 * (Mt2 + 2. * Sqq)) / (3. * Sqq2);
        result.element({6, 4, 10, 7}) = signfix * (-64. * Mb2 * Mt * Mu2) / ((Mb + Mu) * Sqq);
        result.element({7, 1, 0, 1}) = (4. * Mb2 * Pw2 * (Mt2 + 2. * Sqq)) / (3. * MbMu2 * Sqq);
        result.element({7, 1, 5, 1}) = (4. * Mb2 * Pw2 * (Mt2 + 2. * Sqq)) / (3. * MbMu2 * Sqq);
        result.element({7, 1, 7, 1}) = (4. * Mb2 * Pw2 * (Mt2 + 2. * Sqq)) / (3. * MbMu2 * Sqq);
        result.element({7, 1, 9, 5}) = (-16. * Mb2 * Mt * Pw2) / ((Mb + Mu) * Sqq);
        result.element({7, 2, 0, 2}) = (2. * Mb2 * Mt2 * Pw2) / Sqq2;
        result.element({7, 2, 1, 0}) = -((Mb2 * Mt * Pw2) / (Mu * Sqq));
        result.element({7, 2, 3, 0}) = (Mb2 * Mt * Pw2) / (Mu * Sqq);
        result.element({7, 2, 5, 2}) = (-2. * Mb2 * Mt2 * Pw2) / Sqq2;
        result.element({7, 2, 7, 2}) = (2. * Mb2 * Mt2 * Pw2) / Sqq2;
        result.element({7, 3, 0, 3}) = (MbMu2 * (Mt2 + 2. * Sqq)) / (3. * Sqq);
        result.element({7, 3, 5, 3}) = -(MbMu2 * (Mt2 + 2. * Sqq)) / (3. * Sqq);
        result.element({7, 3, 7, 3}) = (MbMu2 * (Mt2 + 2. * Sqq)) / (3. * Sqq);
        result.element({7, 3, 9, 6}) = signfix * (-4. * Mt * (Mb - Mu) * MbMu2) / Sqq;
        result.element({7, 4, 0, 4}) = (32. * Mb2 * Mu2 * (Mt2 + 2. * Sqq)) / (3. * Sqq2);
        result.element({7, 4, 5, 4}) = (-32. * Mb2 * Mu2 * (Mt2 + 2. * Sqq)) / (3. * Sqq2);
        result.element({7, 4, 7, 4}) = (32. * Mb2 * Mu2 * (Mt2 + 2. * Sqq)) / (3. * Sqq2);
        result.element({7, 4, 9, 7}) = signfix * (-64. * Mb2 * Mt * Mu2) / ((Mb + Mu) * Sqq);
        result.element({8, 1, 6, 1}) = (4. * Mb2 * Pw2 * (Mt2 + 2. * Sqq)) / (3. * MbMu2 * Sqq);
        result.element({8, 1, 8, 1}) = (4. * Mb2 * Pw2 * (Mt2 + 2. * Sqq)) / (3. * MbMu2 * Sqq);
        result.element({8, 1, 10, 5}) = (-16. * Mb2 * Mt * Pw2) / ((Mb + Mu) * Sqq);
        result.element({8, 2, 2, 0}) = -((Mb2 * Mt * Pw2) / (Mu * Sqq));
        result.element({8, 2, 4, 0}) = (Mb2 * Mt * Pw2) / (Mu * Sqq);
        result.element({8, 2, 6, 2}) = (-2. * Mb2 * Mt2 * Pw2) / Sqq2;
        result.element({8, 2, 8, 2}) = (2. * Mb2 * Mt2 * Pw2) / Sqq2;
        result.element({8, 3, 6, 3}) = -(MbMu2 * (Mt2 + 2. * Sqq)) / (3. * Sqq);
        result.element({8, 3, 8, 3}) = (MbMu2 * (Mt2 + 2. * Sqq)) / (3. * Sqq);
        result.element({8, 3, 10, 6}) = signfix * (4. * Mt * (Mb - Mu) * MbMu2) / Sqq;
        result.element({8, 4, 6, 4}) = (-32. * Mb2 * Mu2 * (Mt2 + 2. * Sqq)) / (3. * Sqq2);
        result.element({8, 4, 8, 4}) = (32. * Mb2 * Mu2 * (Mt2 + 2. * Sqq)) / (3. * Sqq2);
        result.element({8, 4, 10, 7}) = signfix * (64. * Mb2 * Mt * Mu2) / ((Mb + Mu) * Sqq);
        result.element({9, 5, 0, 1}) = (-16. * Mb2 * Mt * Pw2) / ((Mb + Mu) * Sqq);
        result.element({9, 5, 5, 1}) = (-16. * Mb2 * Mt * Pw2) / ((Mb + Mu) * Sqq);
        result.element({9, 5, 7, 1}) = (-16. * Mb2 * Mt * Pw2) / ((Mb + Mu) * Sqq);
        result.element({9, 5, 9, 5}) = (64. * Mb2 * Pw2 * (2. * Mt2 + Sqq)) / (3. * Sqq2);
        result.element({9, 6, 0, 3}) = signfix * (-4. * Mt * (Mb - Mu) * MbMu2) / Sqq;
        result.element({9, 6, 5, 3}) = signfix * (4. * Mt * (Mb - Mu) * MbMu2) / Sqq;
        result.element({9, 6, 7, 3}) = signfix * (-4. * Mt * (Mb - Mu) * MbMu2) / Sqq;
        result.element({9, 6, 9, 6}) =
            signfix * signfix * (16. * (2. * Mt2 + Sqq) * ((2. * Mu2 - Sqq) * Sqq + 2. * Mb2 * (2. * Pw2 + Sqq))) / (3. * Sqq2);
        result.element({9, 7, 0, 4}) = signfix * (-64. * Mb2 * Mt * Mu2) / ((Mb + Mu) * Sqq);
        result.element({9, 7, 5, 4}) = signfix * (64. * Mb2 * Mt * Mu2) / ((Mb + Mu) * Sqq);
        result.element({9, 7, 7, 4}) = signfix * (-64. * Mb2 * Mt * Mu2) / ((Mb + Mu) * Sqq);
        result.element({9, 7, 9, 7}) = signfix * signfix * (128. * Mb2 * Mu2 * (2. * Mt2 + Sqq)) / (3. * MbMu2 * Sqq);
        result.element({10, 5, 6, 1}) = (-16. * Mb2 * Mt * Pw2) / ((Mb + Mu) * Sqq);
        result.element({10, 5, 8, 1}) = (-16. * Mb2 * Mt * Pw2) / ((Mb + Mu) * Sqq);
        result.element({10, 5, 10, 5}) = (64. * Mb2 * Pw2 * (2. * Mt2 + Sqq)) / (3. * Sqq2);
        result.element({10, 6, 6, 3}) = signfix * (-4. * Mt * (Mb - Mu) * MbMu2) / Sqq;
        result.element({10, 6, 8, 3}) = signfix * (4. * Mt * (Mb - Mu) * MbMu2) / Sqq;
        result.element({10, 6, 10, 6}) =
            signfix * signfix * (16. * (2. * Mt2 + Sqq) * ((2. * Mu2 - Sqq) * Sqq + 2. * Mb2 * (2. * Pw2 + Sqq))) / (3. * Sqq2);
        result.element({10, 7, 6, 4}) = signfix * (-64. * Mb2 * Mt * Mu2) / ((Mb + Mu) * Sqq);
        result.element({10, 7, 8, 4}) = signfix * (64. * Mb2 * Mt * Mu2) / ((Mb + Mu) * Sqq);
        result.element({10, 7, 10, 7}) = signfix * signfix * (128. * Mb2 * Mu2 * (2. * Mt2 + Sqq)) / (3. * MbMu2 * Sqq);

        result *= RateNorm;

        return result;
    }

} // namespace Hammer
