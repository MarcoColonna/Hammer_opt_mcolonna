///
/// @file  RateTau3PiNu.cc
/// @brief \f$ \tau \rightarrow 3\pi \nu \f$ total rate
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/Rates/RateTau3PiNu.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace PS {

        inline BoundaryFunction makeS1Function(double Mp) {
            return ([Mp](const vector<double>& vals) -> pair<double, double> {
                return make_pair(4. * Mp * Mp, pow(sqrt(vals[0]) - Mp, 2.));
            });
        }

        inline BoundaryFunction makeS2Function(double Mp) {
            return ([Mp](const vector<double>& vals) -> pair<double, double> {
                double Mp2 = Mp * Mp;
                double E12 = (vals[0] - Mp2 - vals[1]) / (2 * sqrt(vals[1]));
                double sqTerm = (vals[0] - vals[1]) / 2. + 3. * Mp2 / 2.;
                double cosTerm = 2. * sqrt(E12 * E12 - Mp2) * sqrt(vals[1] / 4. - Mp2);
                return make_pair(sqTerm - cosTerm, sqTerm + cosTerm);
            });
        }
    } // namespace PS

    namespace MD = MultiDimensional;

    RateTau3PiNu::RateTau3PiNu() {
        // Create tensor rank and dimensions
        string name{"RateTau3PiInt"};
        auto totNumPoints = static_cast<IndexType>(_nPoints * _nPoints * _nPoints);
        IndexList dims{{3, 3, totNumPoints}};
        auto& pdg = PID::instance();
        // Ordering is piminus then piplus as anti-tau parent has pdg < 0
        addProcessSignature(PID::ANTITAU, {PID::PIMINUS, PID::PIPLUS, PID::PIPLUS, PID::NU_TAUBAR});
        addIntegrationBoundaries({PS::makeQ2Function(3. * pdg.getMass(PID::PIPLUS), pdg.getMass(PID::TAU)),
                                  PS::makeS1Function(pdg.getMass(PID::PIPLUS)),
                                  PS::makeS2Function(pdg.getMass(PID::PIPLUS))});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_TAU3PI, FF_TAU3PI_HC, INTEGRATION_INDEX})});

        setSignatureIndex();
    }

    Tensor RateTau3PiNu::evalAtPSPoint(const vector<double>& point) {
        auto labs = getTensor().labels();
        labs.pop_back();
        auto dimensions = getTensor().dims();
        dimensions.pop_back();
        Tensor result{"RateTau3Pi", MD::makeEmptySparse(dimensions, labs)};

        const double Mt = masses()[0];
        const double Mp = masses()[1];
        const double Mt2 = Mt * Mt;
        const double Mp2 = Mp * Mp;

        const double Sqp = point[0];
        const double s1 = point[1];
        const double s2 = point[2];

        const double Ep1 = (Sqp - s1 + Mp2) / (2 * sqrt(Sqp));
        const double Ep2 = (Sqp - s2 + Mp2) / (2 * sqrt(Sqp));
        const double Kp12 = Ep1 * Ep1 - Mp2;
        // const double Kp1 = sqrt(Kp12);
        const double Kp22 = Ep2 * Ep2 - Mp2;
        // const double Kp2 = sqrt(Kp22);

        const double Kp1Kp2CosTp12 = (2 * Ep1 * Ep2 + Mp2 - 2 * (Ep1 + Ep2) * sqrt(Sqp) + Sqp) / 2.;

        double RateNorm = pow(GFermi, 2.) * pow(Mt2 - Sqp, 2.) / (4096. * pow(Mt, 3.) * pow(pi, 5.) * Sqp);

        // set non-zero tensor elements
        result.element({0, 0}) = ((Mt2 + 2. * Sqp) * (4. * Kp12 + Kp22 + 4. * Kp1Kp2CosTp12)) / (3. * Sqp);
        result.element({0, 1}) = ((Mt2 + 2. * Sqp) * (2. * (Kp12 + Kp22) + 5. * Kp1Kp2CosTp12)) / (3. * Sqp);
        result.element({1, 0}) = ((Mt2 + 2. * Sqp) * (2. * (Kp12 + Kp22) + 5. * Kp1Kp2CosTp12)) / (3. * Sqp);
        result.element({1, 1}) = ((Mt2 + 2. * Sqp) * (Kp12 + 4. * Kp22 + 4. * Kp1Kp2CosTp12)) / (3. * Sqp);
        result.element({2, 2}) = Mt2;

        // Include extra 1/2 factor following tauola code (origin to be understood)
        result *= 0.5 * RateNorm;

        return result;
    }

} // namespace Hammer


//        inline BoundaryFunction makeS1Function(double Mpi) {
//            return ([Mpi](const vector<double>& vals) -> pair<double, double> {
//                return make_pair(pow(2. * Mpi, 2.), pow(sqrt(vals[0]) - Mpi, 2.));
//            });
//        }
//
//        inline BoundaryFunction makeS2Function(double Mpi) {
//            double Mpi2 = pow(Mpi, 2.);
//            return ([Mpi2](const vector<double>& vals) -> pair<double, double> {
//                double lambda1 = sqrt(vals[0] * vals[0] + vals[1] * vals[1] + Mpi2 * Mpi2 -
//                                      2. * (vals[0] * vals[1] + vals[1] * Mpi2 + vals[0] * Mpi2));
//                double lambda2 = sqrt(vals[1] * vals[1] - 4. * vals[1] * Mpi2);
//                double coeff = 1. / (4. * vals[1]);
//                double a = pow(vals[0] - Mpi2, 2.);
//                return make_pair(coeff * (a - pow(lambda1 + lambda2, 2.)), coeff * (a - pow(lambda1 - lambda2, 2.)));
//            });
//        }
