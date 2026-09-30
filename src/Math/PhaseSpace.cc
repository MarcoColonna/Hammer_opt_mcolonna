///
/// @file  PhaseSpace.cc
/// @brief Phase space integrals
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include <boost/math/special_functions/ellint_3.hpp>

#include "Hammer/Math/PhaseSpace.hh"
#include "Hammer/Exceptions.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/Integrator.hh"
#include "Hammer/Math/Utils.hh"

using namespace std;

namespace Hammer {

    static BoundaryFunction makeMjFunction(double parentMass, vector<double> masses, size_t j) {
        return ([parentMass, masses = std::move(masses), j](const vector<double>& vals) -> pair<double, double> {
            size_t N = masses.size();
            double MjMin = 0;
            double MjMax = 0;
            if (j > 1 && j < N) { // j is number of particles in invariant mass
                for (size_t idx = 0; idx < j; ++idx) {
                    MjMin += masses[idx];
                }
                if (j < N - 1) {
                    // j = N-2 takes val[0] - masses[N-2] as upper range, etc
                    MjMax = vals[N - 2 - j] - masses[j];
                } else {
                    // j = N-1 takes Parent - masses[N-1]
                    MjMax = parentMass - masses[j];
                }
            }
            return make_pair(MjMin, MjMax);
        });
    }

    static double kStar(double mParent, double mk, double mSibling) {
        double EStar = (mParent * mParent - mSibling * mSibling + mk * mk) / (2. * mParent);
        return sqrt(EStar * EStar - mk * mk);
    }

    double phaseSpaceNBody(const double parentMass, const vector<double>& daughterMasses) {
        size_t N = daughterMasses.size();
        double res = 0.;
        switch (N) {
        case 1: {
            throw Error("Implementing a 1->1 decay is perhaps a poor idea.");
        }
        case 2: {
            res += kStar(parentMass, daughterMasses[1], daughterMasses[0]) / (4. * pi * parentMass);
            break;
        }
        default:
            Integrator integ;
            IntegrationBoundaries boundaries{};
            for (size_t j = 2; j < N; ++j) {
                boundaries.push_back(makeMjFunction(parentMass, daughterMasses, N + 1 - j));
            }
            integ.applyRange(boundaries);
            const EvaluationGrid& points = integ.getEvaluationPoints();
            const EvaluationWeights& weights = integ.getEvaluationWeights();

            for (size_t i = 0; i < points.size(); ++i) {
                double respart =
                    kStar(parentMass, daughterMasses[N - 1], points[i][0]) * weights[i] / (4. * pi * parentMass);
                for (size_t j = 0; j < N - 3; ++j) {
                    respart *= kStar(points[i][j], daughterMasses[N - 2 - j], points[i][j + 1]) / (4. * pi2);
                }
                respart *= kStar(points[i][N - 3], daughterMasses[1], daughterMasses[0]) / (4. * pi2);
                res += respart;
            }
            break;
        }
        return res;
    }

    static double phaseSpace2(const double m0, const double m1, const double m2) {
        const double coeff = 1. / (8. * pi);
        const double m0sq = m0 * m0;
        const double m1sq = m1 * m1;
        const double m2sq = m2 * m2;
        const double kallen = m0sq * m0sq + m1sq * m1sq + m2sq * m2sq - 2. * (m0sq * m1sq + m1sq * m2sq + m0sq * m2sq);
        const double result = coeff * sqrt(kallen) / m0sq;
        return result;
    }


    // taken from http://cds.cern.ch/record/682356/files/0311075.pdf
    static double phaseSpace3(const double m0, const double m1, const double m2, const double m3) {
        const double coeff = 1. / (128. * pi3);

        const double m0sq = m0 * m0;
        const double m1sq = m1 * m1;
        const double m2sq = m2 * m2;
        const double m3sq = m3 * m3;

        const double sq0p3 = pow(m0 + m3, 2.);
        const double sq0m3 = pow(m0 - m3, 2.);
        const double sq1p2 = pow(m1 + m2, 2.);
        const double sq1m2 = pow(m1 - m2, 2.);

        const double qpp = sq0p3 - sq1p2;
        const double qpm = sq0p3 - sq1m2;
        const double qmp = sq0m3 - sq1p2;
        const double qmm = sq0m3 - sq1m2;

        const double Qplus = qpp * qmm;
        const double Qminus = qmp * qpm;

        double k = regularize(sqrt(Qminus / Qplus), 1., 1.e-10, -1);

        const double alpha1sq = isZero(m1 + m2) ? 1. : regularize(qmp / qmm, 1., 1.e-10, -1);
        const double alpha2sq = isZero(m1 + m2) ? 1. : sq1m2 / sq1p2 * alpha1sq;


        // Boost elliptic ell_int1/2(k) and ellint_3(k,n) functions require -1 <= k <=1 and n < 1
        const double ellipticE = boost::math::ellint_2(k);
        const double ellipticK = boost::math::ellint_1(k);
        const double ellipticPi1 = boost::math::ellint_3(k, alpha1sq);
        const double ellipticPi2 = boost::math::ellint_3(k, alpha2sq);

        const double expr = 0.5 * Qplus * (m1sq + m2sq + m3sq + m0sq) * ellipticE +
                            4. * m1 * m2 *
                                ((sq0m3 - sq1m2) * (sq0p3 - m3 * m0 + m1 * m2) * ellipticK +
                                 2. * ((m1sq + m2sq) * (m3sq + m0sq) - 2. * (m1sq * m2sq + m3sq * m0sq)) * ellipticPi1 -
                                 2. * pow(m0sq - m3sq, 2.) * ellipticPi2);

        const double result = coeff * expr / m0sq / sqrt(Qplus);

        return result;
    }


    double phaseSpaceN(const double parentMass, const vector<double>& daughterMasses) {
        switch (daughterMasses.size()) {
        case 2: {
            return phaseSpace2(parentMass, daughterMasses[0], daughterMasses[1]);
        }
        case 3: {
            return phaseSpace3(parentMass, daughterMasses[0], daughterMasses[1], daughterMasses[2]);
        }
        case 4: {
            function<double(double)> f = [parentMass, &daughterMasses](double m) -> double {
                const double f1 = phaseSpace3(m, daughterMasses[0], daughterMasses[1], daughterMasses[2]);
                const double f2 = phaseSpace2(parentMass, m, daughterMasses[3]);
                return 2 * m * f1 * f2 / (2 * pi);
            };
            return integrate(f, daughterMasses[0] + daughterMasses[1] + daughterMasses[2],
                             parentMass - daughterMasses[3]);
        }
        case 5: {
            function<double(double)> f = [parentMass, &daughterMasses](double m) -> double {
                const double f1 = phaseSpace3(m, daughterMasses[0], daughterMasses[1], daughterMasses[2]);
                const double f2 = phaseSpace3(parentMass, m, daughterMasses[3], daughterMasses[4]);
                return 2 * m * f1 * f2 / (2 * pi);
            };
            return integrate(f, daughterMasses[0] + daughterMasses[1] + daughterMasses[2],
                             parentMass - daughterMasses[3] - daughterMasses[4]);
        }
        default:
            throw Error("Unimplemented phase space integral: n>5");
        }
    }

} // namespace Hammer
