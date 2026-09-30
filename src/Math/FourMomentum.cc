///
/// @file  FourMomentum.cc
/// @brief Hammer four momentum class
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>
#include <string>

#include "Hammer/Math/FourMomentum.hh"
#include "Hammer/Exceptions.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Tools/HammerRoot.hh"

using namespace std;

namespace Hammer {

    FourMomentum::FourMomentum() : _vec{{0., 0., 0., 0.}} {
    }

    // set the four momentum in format FourMomentum({E,px,py,pz})
    FourMomentum::FourMomentum(const array<double, 4>& other) : _vec(other) {
    }

    // set the four momentum coordinates
    FourMomentum::FourMomentum(const double E, const double px, const double py, const double pz) noexcept
        : _vec{{E, px, py, pz}} {
    }

    FourMomentum& FourMomentum::setE(double E) {
        _vec[0] = E;
        return *this;
    }
    FourMomentum& FourMomentum::setPx(double px) {
        _vec[1] = px;
        return *this;
    }

    FourMomentum& FourMomentum::setPy(double py) {
        _vec[2] = py;
        return *this;
    }

    FourMomentum& FourMomentum::setPz(double pz) {
        _vec[3] = pz;
        return *this;
    }

    FourMomentum FourMomentum::fromPM(double px, double py, double pz, double mass) {
        return FourMomentum{sqrt(px * px + py * py + pz * pz + mass * mass), px, py, pz};
    }

    FourMomentum FourMomentum::fromPtEtaPhiM(double pt, double eta, double phi, double mass) {
        return FourMomentum{sqrt(mass * mass + pow(pt * cosh(eta), 2.)), pt * cos(phi), pt * sin(phi), pt * sinh(eta)};
    }

    FourMomentum FourMomentum::fromEtaPhiME(double eta, double phi, double mass, double E) {
        if (E <= mass) {
            throw RangeError("Invalid mass and energy relation");
        }
        double pt = sqrt((E * E - mass * mass) / (1 + pow(sinh(eta), 2.)));
        return FourMomentum{E, pt * cos(phi), pt * sin(phi), pt * sinh(eta)};
    }

#ifdef HAVE_ROOT

    FourMomentum FourMomentum::fromRoot(const TLorentzVector& v) {
        return FourMomentum{v.E(), v.Px(), v.Py(), v.Pz()};
    }

#endif

    // Values of energy and momentum components.
    double FourMomentum::E() const {
        return _vec[0];
    }

    double FourMomentum::px() const {
        return _vec[1];
    }

    double FourMomentum::py() const {
        return _vec[2];
    }

    double FourMomentum::pz() const {
        return _vec[3];
    }


    // extract mass
    double FourMomentum::mass2() const {
        return (E() * E() - p2());
    }

    double FourMomentum::mass() const {
        if (fuzzyLess(mass2(), 0.0)) { // less than 'a little bit negative'
            throw Error("Rosebud! Negative mass^2: " + to_string(mass2()));
        }
        return sqrt(fabs(mass2()));
    }

    // absolute value of spatial momentum, p_T, rapidity, pseudorapidity, azimuthal and polar angles
    double FourMomentum::p2() const {
        return px() * px() + py() * py() + pz() * pz();
    }

    double FourMomentum::p() const {
        return sqrt(p2());
    }

    double FourMomentum::pt() const {
        return sqrt(px() * px() + py() * py());
    }

    double FourMomentum::rapidity() const {
        return acosh(E() / mass());
    }

    double FourMomentum::theta() const {
        return atan2(pt(), pz());
    }

    double FourMomentum::eta() const {
        return -log(tan(theta() / 2.));
    }

    double FourMomentum::phi() const {
        return atan2(py(), px());
    }

    array<double, 3> FourMomentum::pVec() const {
        return array<double, 3>{{px(), py(), pz()}};
    }

    // dot operator, trace -2 metric
    double FourMomentum::dot(const FourMomentum& v) const {
        return E() * v.E() - px() * v.px() - py() * v.py() - pz() * v.pz();
    }

    double operator*(const FourMomentum& v, const FourMomentum& w) {
        return dot(v, w);
    }

    double dot(const FourMomentum& v, const FourMomentum& w) {
        return v.dot(w);
    }

    FourMomentum& FourMomentum::operator+=(const FourMomentum& v) {
        for (size_t i = 0; i <= 3; ++i) {
            _vec[i] += v._vec[i];
        }
        return *this;
    }

    FourMomentum& FourMomentum::operator-=(const FourMomentum& v) {
        for (size_t i = 0; i <= 3; ++i) {
            _vec[i] -= v._vec[i];
        }
        return *this;
    }

    FourMomentum FourMomentum::operator-() const {
        FourMomentum res{-_vec[0], -_vec[1], -_vec[2], -_vec[3]};
        return res;
    }

    FourMomentum& FourMomentum::operator*=(double a) {
        for (size_t i = 0; i <= 3; ++i) {
            _vec[i] *= a;
        }
        return *this;
    }

    FourMomentum& FourMomentum::operator/=(double a) {
        for (size_t i = 0; i <= 3; ++i) {
            _vec[i] /= a;
        }
        return *this;
    }

    FourMomentum FourMomentum::PFlip() const {
        FourMomentum res{_vec[0], -_vec[1], -_vec[2], -_vec[3]};
        return res;
    }

    // boost variables
    double FourMomentum::gamma() const {
        double m = mass();
        if (isZero(m)) {
            throw Error("Rosebud! Zero mass: " + to_string(m));
        }
        return E() / m;
    }

    double FourMomentum::beta() const {
        double energy = E();
        if (isZero(energy)) {
            throw Error("Rosebud! Zero energy: " + to_string(energy));
        }
        return p() / energy;
    }

    array<double, 3> FourMomentum::boostVector() const {
        double energy = E();
        if (isZero(energy)) {
            throw Error("Rosebud! Zero energy: " + to_string(energy));
        }
        return {{_vec[1] / energy, _vec[2] / energy, _vec[3] / energy}};
    }

    void FourMomentum::boostBy(const std::array<double, 3>& v) {
        const double vSq = ::Hammer::dot(v, v);
        ;
        const double gamma = 1.0 / sqrt(1.0 - vSq);
        const double vp = ::Hammer::dot(v, pVec());
        const double gamma2 = vSq > 0 ? (gamma - 1.0) / vSq : 0.0;

        _vec[1] = (px() + (gamma2 * vp - gamma * E()) * v[0]);
        _vec[2] = (py() + (gamma2 * vp - gamma * E()) * v[1]);
        _vec[3] = (pz() + (gamma2 * vp - gamma * E()) * v[2]);
        _vec[0] = gamma * (E() - vp);
    }

    FourMomentum& FourMomentum::boostToRestFrameOf(const FourMomentum& v) {
        boostBy(v.boostVector());
        return *this;
    }

    FourMomentum& FourMomentum::boostToRestFrameOf(const std::array<double, 3>& v) {
        boostBy(v);
        return *this;
    }

    FourMomentum& FourMomentum::boostFromRestFrameOf(const FourMomentum& v) {
        boostBy(v.PFlip().boostVector());
        return *this;
    }

    // Delta collider coordinates
    double angle(const FourMomentum& v, const FourMomentum& w) {
        return acos(costheta(v.pVec(), w.pVec()));
    }

    double deltaPhi(const FourMomentum& v, const FourMomentum& w) {
        return abs(v.phi() - w.phi());
    }

    double deltaEta(const FourMomentum& v, const FourMomentum& w) {
        return abs(v.eta() - w.eta());
    }

    double deltaR(const FourMomentum& v, const FourMomentum& w) {
        return sqrt(pow(deltaPhi(v, w), 2) + pow(deltaEta(v, w), 2));
    }

    // epsilon contraction.
    double epsilon(const FourMomentum& a, const FourMomentum& b, const FourMomentum& c, const FourMomentum& d) {
        array<double, 4> p = {{a.E(), a.px(), a.py(), a.pz()}};
        array<double, 4> q = {{b.E(), b.px(), b.py(), b.pz()}};
        array<double, 4> r = {{c.E(), c.px(), c.py(), c.pz()}};
        array<double, 4> s = {{d.E(), d.px(), d.py(), d.pz()}};
        double temp = +(-p[3] * q[2] * r[1] * s[0] + p[2] * q[3] * r[1] * s[0] + p[3] * q[1] * r[2] * s[0] -
                        p[1] * q[3] * r[2] * s[0] - p[2] * q[1] * r[3] * s[0] + p[1] * q[2] * r[3] * s[0] +
                        p[3] * q[2] * r[0] * s[1] - p[2] * q[3] * r[0] * s[1] - p[3] * q[0] * r[2] * s[1] +
                        p[0] * q[3] * r[2] * s[1] + p[2] * q[0] * r[3] * s[1] - p[0] * q[2] * r[3] * s[1] -
                        p[3] * q[1] * r[0] * s[2] + p[1] * q[3] * r[0] * s[2] + p[3] * q[0] * r[1] * s[2] -
                        p[0] * q[3] * r[1] * s[2] - p[1] * q[0] * r[3] * s[2] + p[0] * q[1] * r[3] * s[2] +
                        p[2] * q[1] * r[0] * s[3] - p[1] * q[2] * r[0] * s[3] - p[2] * q[0] * r[1] * s[3] +
                        p[0] * q[2] * r[1] * s[3] + p[1] * q[0] * r[2] * s[3] - p[0] * q[1] * r[2] * s[3]);
        return temp;
    }

    FourMomentum boostToRestFrameOf(const FourMomentum& mom, double vx, double vy, double vz) {
        FourMomentum result = mom;
        result.boostToRestFrameOf(array<double, 3>{{vx, vy, vz}});
        return result;
    }

    FourMomentum boostToRestFrameOf(const FourMomentum& mom, const array<double, 3>& v) {
        FourMomentum result = mom;
        result.boostToRestFrameOf(v);
        return result;
    }

    FourMomentum boostToRestFrameOf(const FourMomentum& mom, const FourMomentum& v) {
        FourMomentum result = mom;
        result.boostToRestFrameOf(v);
        return result;
    }

    double dot(const array<double, 3>& a, const array<double, 3>& b) {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    }

    double costheta(const array<double, 3>& a, const array<double, 3>& b) {
        double res = dot(a, b) / sqrt(dot(a, a) * dot(b, b));
        if (is_ieee754_nan(res)) {
            throw Error("Nan in Cosine. Cock. In Hungarian.");
        }
        return res;
    }

} // namespace Hammer
