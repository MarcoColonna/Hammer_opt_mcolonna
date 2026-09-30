///
/// @file  AsCoeffs.cc
/// @brief Hammer class for HQET \f$ \alpha_s \f$ corrections
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-

#include <cmath>

#include "Hammer/FormFactors/AsCoeffs.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/Differentiator.hh"

using namespace std;

namespace Hammer {

    void AsCoeffs::updateVars(double w, double z) const {
        double adjW = max(w, 1.);
        bool upW = !isZero(adjW - _w);
        bool upZ = !isZero(z - _z);
        if (upW) {
            _w = adjW;
            _wSq = _w * _w;
            _sqrt1wSq = sqrt(_wSq - 1);
            _wP = _w + _sqrt1wSq;
            _wM = _w - _sqrt1wSq;
            _lnwP = log(_wP);
            _rW = fabs(_w - 1.) > numeric_limits<double>::min() ? _lnwP / _sqrt1wSq : 1.;
        }
        if (upZ) {
            _z = z;
            _zSq = _z * _z;
            _zCu = _z * _z * _z;
            _zp1Sq = (_z + 1.) * (_z + 1.);
            _zm1Sq = (_z - 1.) * (_z - 1.);
            ;
            _lnz = log(_z);
            _wZ = (_z + 1. / _z) / 2.;
        }
        if (upW || upZ) {
            _wmwZSq = (_w - _wZ) * (_w - _wZ);
            _mOmega = fabs(_w - 1.) > numeric_limits<double>::min()
                          ? (_w / (2. * _sqrt1wSq) *
                                 (2 * DiLog(1 - _wM * _z) - 2 * DiLog(1 - _wP * _z) + DiLog(1 - _wP * _wP) -
                                  DiLog(1 - _wM * _wM)) -
                             _w * _rW * _lnz + 1.)
                          : (-1. + (1 + _z) / (_z - 1.) * _lnz);
        }
    }

    double AsCoeffs::CS(double w, double z) const {
        updateVars(w, z);
        return (2. * _mOmega * (_w - _wZ) * _z + (-1. + _zSq) * _lnz - (-1. + _w) * _zp1Sq * _rW) /
               (3. * (_w - _wZ) * _z);
    }

    double AsCoeffs::CP(double w, double z) const {
        updateVars(w, z);
        return (2. * _mOmega * (_w - _wZ) * _z + (-1. + _zSq) * _lnz - (1. + _w) * _zm1Sq * _rW) /
               (3. * (_w - _wZ) * _z);
    }

    double AsCoeffs::CV1(double w, double z) const {
        updateVars(w, z);
        return (4. * _mOmega * (_w - _wZ) * _z + 12 * (-_w + _wZ) * _z + _lnz - _zSq * _lnz +
                2. * (1. + _w) * (-1. + (-1. + 3. * _w - _z) * _z) * _rW) /
               (6. * (_w - _wZ) * _z);
    }

    double AsCoeffs::CV2(double w, double z) const {
        updateVars(w, z);
        return (-(_z * (-2. * (_w - _wZ) * (-1. + _z) + (3. - 2. * _w + (2. - 4. * _w) * _z + _zSq) * _lnz)) +
                (-2. + (-1. + 5. * _w + 2. * _wSq) * _z - 2. * _w * (1. + 2. * _w) * _zSq + (1. + _w) * _zCu) * _rW) /
               (6. * _wmwZSq * _zSq);
    }

    double AsCoeffs::CV3(double w, double z) const {
        updateVars(w, z);
        return (-2. * (_w - _wZ) * (-1. + _z) * _z + (1. + (2. - 4. * _w) * _z + (3. - 2. * _w) * _zSq) * _lnz +
                (1. + 2. * _wSq * (-2. + _z) * _z - _zSq - 2. * _zCu + _w * (1. - 2. * _z + 5. * _zSq)) * _rW) /
               (6. * _wmwZSq * _z);
    }

    double AsCoeffs::CA1(double w, double z) const {
        updateVars(w, z);
        return (-12 * _w * _z + 4. * _mOmega * (_w - _wZ) * _z + 12 * _wZ * _z + _lnz - _zSq * _lnz +
                2. * (-1. + _w) * (-1. + _z + 3. * _w * _z - _zSq) * _rW) /
               (6. * (_w - _wZ) * _z);
    }

    double AsCoeffs::CA2(double w, double z) const {
        updateVars(w, z);
        return (-(_z * (-2. * (_w - _wZ) * (1. + _z) + (3. + 2. * _w + (-2. - 4. * _w) * _z + _zSq) * _lnz)) +
                (-2. + (1. + 5. * _w - 2. * _wSq) * _z + 2. * (1. - 2. * _w) * _w * _zSq + (-1. + _w) * _zCu) * _rW) /
               (6. * _wmwZSq * _zSq);
    }

    double AsCoeffs::CA3(double w, double z) const {
        updateVars(w, z);
        return (-2. * (_w - _wZ) * _z * (1. + _z) + (-1. + (2. + 4. * _w) * _z + (-3. - 2. * _w) * _zSq) * _lnz +
                (1. - _zSq + 2. * _zCu + 2. * _wSq * _z * (2. + _z) - _w * (1. + 2. * _z + 5. * _zSq)) * _rW) /
               (6. * _wmwZSq * _z);
    }

    double AsCoeffs::CT1(double w, double z) const {
        updateVars(w, z);
        return (-6. * _w * _z + 2. * _mOmega * (_w - _wZ) * _z + 6. * _wZ * _z + _lnz - _zSq * _lnz +
                (-1. + _w) * (-1. + (2. + 4. * _w) * _z - _zSq) * _rW) /
               (3. * (_w - _wZ) * _z);
    }

    double AsCoeffs::CT2(double w, double z) const {
        updateVars(w, z);
        return (-2. * ((1. - _w * _z) * _rW + _z * _lnz)) / (3. * (-_w + _wZ) * _z);
    }

    double AsCoeffs::CT3(double w, double z) const {
        updateVars(w, z);
        return (-2. * ((_w - _z) * _rW + _lnz)) / (3. * (-_w + _wZ));
    }

    double AsCoeffs::derCS(double w, double z) const {
        auto f = [this, z](double x) -> double { return CS(x, z); };
        return der<1>(f, w, {1., std::numeric_limits<double>::max()});
    }

    double AsCoeffs::derCP(double w, double z) const {
        auto f = [this, z](double x) -> double { return CP(x, z); };
        return der<1>(f, w, {1., std::numeric_limits<double>::max()});
    }

    double AsCoeffs::derCV1(double w, double z) const {
        auto f = [this, z](double x) -> double { return CV1(x, z); };
        return der<1>(f, w, {1., std::numeric_limits<double>::max()});
    }

    double AsCoeffs::derCV2(double w, double z) const {
        auto f = [this, z](double x) -> double { return CV2(x, z); };
        return der<1>(f, w, {1., std::numeric_limits<double>::max()});
    }

    double AsCoeffs::derCV3(double w, double z) const {
        auto f = [this, z](double x) -> double { return CV3(x, z); };
        return der<1>(f, w, {1., std::numeric_limits<double>::max()});
    }

    double AsCoeffs::derCA1(double w, double z) const {
        auto f = [this, z](double x) -> double { return CA1(x, z); };
        return der<1>(f, w, {1., std::numeric_limits<double>::max()});
    }

    double AsCoeffs::derCA2(double w, double z) const {
        auto f = [this, z](double x) -> double { return CA2(x, z); };
        return der<1>(f, w, {1., std::numeric_limits<double>::max()});
    }

    double AsCoeffs::derCA3(double w, double z) const {
        auto f = [this, z](double x) -> double { return CA3(x, z); };
        return der<1>(f, w, {1., std::numeric_limits<double>::max()});
    }

    double AsCoeffs::derCT1(double w, double z) const {
        auto f = [this, z](double x) -> double { return CT1(x, z); };
        return der<1>(f, w, {1., std::numeric_limits<double>::max()});
    }

    double AsCoeffs::derCT2(double w, double z) const {
        auto f = [this, z](double x) -> double { return CT2(x, z); };
        return der<1>(f, w, {1., std::numeric_limits<double>::max()});
    }

    double AsCoeffs::derCT3(double w, double z) const {
        auto f = [this, z](double x) -> double { return CT3(x, z); };
        return der<1>(f, w, {1., std::numeric_limits<double>::max()});
    }

    double AsCoeffs::DiLog(double z) const {
        // ====================================================================
        // This file is part of Dilogarithm.
        //
        // Dilogarithm is licenced under the GNU Lesser General Public License
        // (GNU LGPL) version 3.
        // Code translated by R.Brun from CERNLIB DILOG function C332
        // ====================================================================
        constexpr double HF = 0.5;
        constexpr double PI3 = pi2 / 3.;
        constexpr double PI6 = pi2 / 6.;
        constexpr double PI12 = pi2 / 12.;
        constexpr array<double, 20> C{
            0.42996693560813697,  0.40975987533077105, -0.01858843665014592, 0.00145751084062268,
            -0.00014304184442340, 0.00001588415541880, -0.00000190784959387, 0.00000024195180854,
            -0.00000003193341274, 0.00000000434545063, -0.00000000060578480, 0.00000000008612098,
            -0.00000000001244332, 0.00000000000182256, -0.00000000000027007, 0.00000000000004042,
            -0.00000000000000610, 0.00000000000000093, -0.00000000000000014, 0.00000000000000002};

        double H;

        if (isZero(_z - 1.)) {
            H = PI6;
        } else if (isZero(_z + 1.)) {
            H = -PI12;
        } else {
            double T;
            double Y;
            double S;
            double A;
            double ALFA;
            double B1;
            double B2;
            double B0;
            T = -z;
            if (T <= -2) {
                Y = -1 / (1 + T);
                S = 1;
                B1 = log(-T);
                B2 = log1p(1 / T);
                A = -PI3 + HF * (B1 * B1 - B2 * B2);
            } else if (T < -1) {
                Y = -1 - T;
                S = -1;
                A = log(-T);
                A = -PI6 + A * (A + log1p(1 / T));
            } else if (T <= -0.5) {
                Y = -(1 + T) / T;
                S = 1;
                A = log(-T);
                A = -PI6 + A * (-HF * A + log1p(T));
            } else if (T < 0) {
                Y = -T / (1 + T);
                S = -1;
                B1 = log1p(T);
                A = HF * B1 * B1;
            } else if (T <= 1) {
                Y = T;
                S = 1;
                A = 0;
            } else {
                Y = 1 / T;
                S = -1;
                B1 = log(T);
                A = PI6 + HF * B1 * B1;
            }
            H = Y + Y - 1;
            ALFA = H + H;
            B1 = 0;
            B2 = 0;
            B0 = 0;
            for (int i = 19; i >= 0; i--) {
                B0 = C[static_cast<size_t>(i)] + ALFA * B1 - B2;
                B2 = B1;
                B1 = B0;
            }
            H = -(S * (B0 - H * B2) + A);
        }
        return H;
    }

} // namespace Hammer
