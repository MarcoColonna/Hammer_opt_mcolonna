///
/// @file  Pdg.cc
/// @brief Hammer particle data class
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>
#include <regex>

#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Tools/ParticleUtils.hh"

using namespace std;

namespace Hammer {

    PID* PID::_thePID = nullptr;

    PID::~PID() {
        _masses.clear();
        _widths.clear();
        _names.clear();
        _brs.clear();
        if (_thePID != nullptr) {
            PID* pTemp = _thePID;
            _thePID = nullptr;
            delete pTemp;
        }
    }

    PdgId PID::toPdgCode(const string& name) const {
        auto it = _names.find(name);
        if (it == _names.end() || it->second.size() != 1) {
            return 0;
        }
        return it->second.front();
    }

    vector<PdgId> PID::toPdgList(const string& name) const {
        auto it = _names.find(name);
        if (it == _names.end()) {
            return {};
        }
        return it->second;
    }

    double PID::getMass(PdgId id) const {
        auto it = _masses.find(abspid(id));
        if (it != _masses.end()) {
            return it->second;
        }
        return -1.;
    }

    double PID::getMass(const string& name) const {
        return getMass(toPdgCode(name));
    }


    void PID::setMass(PdgId id, double value) {
        auto it = _masses.find(abspid(id));
        if (it != _masses.end()) {
            it->second = value;
        } else {
            _masses.insert(make_pair(abspid(id), value));
        }
    }

    void PID::setMass(const string& name, double value) {
        setMass(toPdgCode(name), value);
    }

    double PID::getWidth(PdgId id) const {
        auto it = _widths.find(abspid(id));
        if (it != _widths.end()) {
            return it->second;
        }
        return -1.;
    }

    double PID::getWidth(const string& name) const {
        return getWidth(toPdgCode(name));
    }

    void PID::setWidth(PdgId id, double value) {
        auto it = _widths.find(abspid(id));
        if (it != _widths.end()) {
            it->second = value;
        } else {
            _widths.insert(make_pair(abspid(id), value));
        }
    }

    void PID::setWidth(const string& name, double value) {
        setWidth(toPdgCode(name), value);
    }

    size_t PID::getSpinMultiplicity(PdgId id) {
        if (abspid(id) <= 18 || id == PHOTON) {
            return 2.;
        }
        if (isMeson(id) || isBaryon(id)) {
            return digit(nj, id);
        }
        return 1.;
    }

    size_t PID::getSpinMultiplicity(const string& name) const {
        return getSpinMultiplicity(toPdgCode(name));
    }

    size_t PID::getSpinMultiplicities(const vector<PdgId>& ids) {
        size_t result = 1ul;
        for (auto elem : ids) {
            result *= getSpinMultiplicity(elem);
        }
        return result;
    }

    size_t PID::getSpinMultiplicities(const vector<string>& names) const {
        size_t result = 1ul;
        for (const auto& elem : names) {
            result *= getSpinMultiplicity(elem);
        }
        return result;
    }

    int PID::getThreeCharge(PdgId pid) {
        int charge = 0;
        unsigned short q1 = digit(nq1, pid);
        unsigned short q2 = digit(nq2, pid);
        unsigned short q3 = digit(nq3, pid);
        int ida = abspid(pid);
        int sid = fundamentalID(pid);
        static constexpr std::array<int, 100> ch100{
            -1, 2, -1, 2, -1, 2, -1, 2, 0, 0, -3, 0, -3, 0, -3, 0, -3, 0, 0, 0, 0, 0, 0, 3, 0,
            0,  0, 0,  0, 0,  0, 0,  0, 3, 0, 0,  3, 0,  0, 0,  0, -1, 0, 0, 0, 0, 0, 0, 0, 0,
            0,  6, 3,  6, 0,  0, 0,  0, 0, 0, 0,  0, 0,  0, 0,  0, 0,  0, 0, 0, 0, 0, 0, 0, 0,
            0,  0, 0,  0, 0,  0, 0,  0, 0, 0, 0,  0, 0,  0, 0,  0, 0,  0, 0, 0, 0, 0, 0, 0, 0};

        if (ida == 0 || extraBits(pid) > 0) { // ion or illegal
            return 0;
        }
        if (digit(nj, pid) == 0) { // KL, Ks, or undefined
            return 0;
        }
        if (sid > 0 && sid <= 100) { // use table
            charge = ch100[static_cast<size_t>(sid - 1)];
        } else if (isMeson(pid)) { // mesons
            if (q2 == 3 || q2 == 5) {
                charge = ch100[q3 - 1] - ch100[q2 - 1];
            } else {
                charge = ch100[q2 - 1] - ch100[q3 - 1];
            }
        } else if (isDiQuark(pid)) { // diquarks
            charge = ch100[q2 - 1] + ch100[q1 - 1];
        } else if (isBaryon(pid)) { // baryons
            charge = ch100[q3 - 1] + ch100[q2 - 1] + ch100[q1 - 1];
        } else { // unknown
            return 0;
        }

        if (charge == 0) {
            return 0;
        }
        if (pid < 0) {
            charge = -charge;
        }

        return charge;
    }

    int PID::getThreeCharge(const vector<PdgId>& ids) {
        int result = 0;
        for (auto elem : ids) {
            result += getThreeCharge(elem);
        }
        return result;
    }

    int PID::getLeptonNumber(PdgId id) {
        if (abspid(id) < 11 || abspid(id) > 18) {
            return 0;
        }
        return (id > 0) ? 1 : -1;
    }

    int PID::getLeptonNumber(const vector<PdgId>& ids) {
        int result = 0;
        for (auto elem : ids) {
            result += getLeptonNumber(elem);
        }
        return result;
    }

    int PID::getBaryonNumber(PdgId id) {
        if (isBaryon(id)) {
            return (id > 0) ? 1 : -1;
        }
        return 0;
    }

    int PID::getBaryonNumber(const vector<PdgId>& ids) {
        int result = 0;
        for (auto elem : ids) {
            result += getBaryonNumber(elem);
        }
        return result;
    }

    tuple<int, int, int> PID::getLeptonFlavorNumber(PdgId id) {
        if (!isLepton(id)) {
            return make_tuple(0, 0, 0);
        }
        if (abspid(id) - 11 < 2) {
            return make_tuple((id > 0) ? 1 : -1, 0, 0);
        }
        if (abspid(id) - 13 < 2) {
            return make_tuple(0, (id > 0) ? 1 : -1, 0);
        }
        return make_tuple(0, 0, (id > 0) ? 1 : -1);
    }

    tuple<int, int, int> PID::getLeptonFlavorNumber(const vector<PdgId>& ids) {
        tuple<int, int, int> result = make_tuple(0, 0, 0);
        for (auto elem : ids) {
            tuple<int, int, int> tmp = getLeptonFlavorNumber(elem);
            get<0>(result) += get<0>(tmp);
            get<1>(result) += get<1>(tmp);
            get<2>(result) += get<2>(tmp);
        }
        return result;
    }

    PID* PID::getPIDInstance() {
        if (_thePID == nullptr) {
            _thePID = new PID();
            _thePID->init();
        }
        return _thePID;
    }

    PID& PID::instance() {
        return *PID::getPIDInstance();
    }

    void PID::addBR(PdgId parent, const vector<PdgId>& daughters, const double br) {
        auto tmp = combineDaughters(daughters, {});
        _brs[processID(parent, tmp)] = make_pair(br, parent);
    }

    map<HashId, double> PID::getPartialWidths() {
        map<HashId, double> res;
        for (auto& elem : _brs) {
            auto it = _widths.find(abspid(elem.second.second));
            if (it != _widths.end()) {
                res.insert({elem.first, (elem.second.first) * (it->second)});
            }
        }
        return res;
    }

    void PID::init() {
        _masses[ELECTRON] = 5.11e-4;
        _masses[MUON] = 0.1056583745;
        _masses[TAU] = 1.77686;
        _masses[NU_E] = 0.;
        _masses[NU_MU] = 0.;
        _masses[NU_TAU] = 0.;
        _masses[PHOTON] = 0.;
        _masses[PROTON] = 0.938272081;
        _masses[NEUTRON] = 0.939565413;
        _masses[PI0] = 0.1349770;
        _masses[PIPLUS] = 0.13957061;
        _masses[RHO0] = 0.77526;
        _masses[RHOPLUS] = 0.77511;
        _masses[OMEGA] = 0.78265;
        _masses[K0L] = 0.497611;
        _masses[K0S] = 0.497611;
        _masses[K0] = 0.497611;
        _masses[KPLUS] = 0.493677;
        _masses[D0] = 1.86483;
        _masses[DPLUS] = 1.86959;
        _masses[DSTAR] = 2.00685;
        _masses[DSTARPLUS] = 2.01026;
        _masses[DSSD0STAR] = 2.300;
        _masses[DSSD0STARPLUS] = 2.349;
        _masses[DSSD1STAR] = 2.427;
        _masses[DSSD1STARPLUS] = 2.427;
        _masses[DSSD1] = 2.421;
        _masses[DSSD1PLUS] = 2.423;
        _masses[DSSD2STAR] = 2.461;
        _masses[DSSD2STARPLUS] = 2.465;
        _masses[DSPLUS] = 1.96828;
        _masses[DSSTARPLUS] = 2.1122;
        _masses[DSSDS0STARPLUS] = 2.3178;
        _masses[DSSDS1STARPLUS] = 2.4595;
        _masses[DSSDS1PLUS] = 2.5351;
        _masses[DSSDS2STARPLUS] = 2.5691;
        _masses[BZERO] = 5.27963;
        _masses[BPLUS] = 5.27932;
        _masses[BS] = 5.36689;
        _masses[LAMBDAB] = 5.620;
        _masses[LAMBDACPLUS] = 2.286;
        _masses[LAMBDACSTAR12PLUS] = 2.592;
        _masses[LAMBDACSTAR32PLUS] = 2.628;
        _masses[BCPLUS] = 6.275;
        _masses[JPSI] = 3.09690;

        const double hbarGevs = 6.582e-25;
        _widths[TAU] = hbarGevs / 2.903e-13;
        _widths[RHO0] = 0.147;
        _widths[RHOPLUS] = 0.149;
        _widths[OMEGA] = 8.49e-3;
        _widths[D0] = hbarGevs / 4.101e-13;
        _widths[DPLUS] = hbarGevs / 1.040e-12;
        _widths[DSTAR] = 2.1e-3;
        _widths[DSTARPLUS] = 83.4e-6;
        _widths[DSSD0STAR] = 0.274;
        _widths[DSSD0STARPLUS] = 0.221;
        _widths[DSSD1STAR] = 0.384;
        _widths[DSSD1STARPLUS] = 0.384;
        _widths[DSSD1] = 31.7e-3;
        _widths[DSSD1PLUS] = 25e-3;
        _widths[DSSD2STAR] = 47.5e-3;
        _widths[DSSD2STARPLUS] = 46.7e-3;
        _widths[DSPLUS] = hbarGevs / 5.04e-13;
        _widths[DSSTARPLUS] = 1.9e-3;
        _widths[DSSDS0STARPLUS] = 3.8e-3;
        _widths[DSSDS1STARPLUS] = 3.5e-3;
        _widths[DSSDS1PLUS] = 0.92e-3;
        _widths[DSSDS2STARPLUS] = 16.9e-3;
        _widths[BZERO] = hbarGevs / 1.52e-12;
        _widths[BPLUS] = hbarGevs / 1.638e-12;
        _widths[BS] = hbarGevs / 1.527e-12;
        _widths[LAMBDAB] = hbarGevs / 1.470e-12;
        _widths[LAMBDACPLUS] = hbarGevs / 2.e-13;
        _widths[LAMBDACSTAR12PLUS] = 0.0026;
        _widths[LAMBDACSTAR32PLUS] = 0.001;
        _widths[BCPLUS] = hbarGevs / 0.510e-12;
        _widths[JPSI] = 92.9e-6;

        addBR(ANTITAU, {NU_TAUBAR, NU_MU, ANTIMUON}, 0.1739);
        addBR(ANTITAU, {NU_TAUBAR, NU_E, POSITRON}, 0.1782);
        addBR(ANTITAU, {PIPLUS, NU_TAUBAR}, 0.1082);
        addBR(RHO0, {PIPLUS, PIMINUS}, 1.0);
        addBR(RHOPLUS, {PIPLUS, PI0}, 1.0);
        addBR(OMEGA, {PIPLUS, PIMINUS, PI0}, 0.893);
        addBR(DSTAR, {D0, PI0}, 0.647);
        addBR(DSTAR, {D0, GAMMA}, 0.353);
        addBR(DSTARMINUS, {DMINUS, PI0}, 0.307);
        addBR(DSTARMINUS, {-D0, PIMINUS}, 0.677);
        addBR(DSTARMINUS, {DMINUS, GAMMA}, 0.016);
        addBR(DSSTARMINUS, {DSMINUS, PI0}, 0.058);
        addBR(DSSTARMINUS, {DSMINUS, GAMMA}, 0.935);

        addBR(DSSD0STAR, {DPLUS, PIMINUS}, 0.67);       // placeholder
        addBR(DSSD0STAR, {D0, PI0}, 0.33);              // placeholder
        addBR(DSSD0STARMINUS, {-D0, PIMINUS}, 0.67);    // placeholder
        addBR(DSSD0STARMINUS, {DMINUS, PI0}, 0.33);     // placeholder
        addBR(DSSD1STAR, {DSTARPLUS, PIMINUS}, 0.67);   // placeholder
        addBR(DSSD1STAR, {DSTAR, PI0}, 0.33);           // placeholder
        addBR(DSSD1STARMINUS, {-DSTAR, PIMINUS}, 0.67); // placeholder
        addBR(DSSD1STARMINUS, {DSTARMINUS, PI0}, 0.33); // placeholder
        addBR(DSSD1, {DSTARPLUS, PIMINUS}, 0.42);       // placeholder
        addBR(DSSD1, {DSTAR, PI0}, 0.21);               // placeholder
        addBR(DSSD1MINUS, {-DSTAR, PIMINUS}, 0.42);     // placeholder
        addBR(DSSD1MINUS, {DSTARMINUS, PI0}, 0.21);     // placeholder
        addBR(DSSD2STAR, {DSTARPLUS, PIMINUS}, 0.26);   // placeholder
        addBR(DSSD2STAR, {DSTAR, PI0}, 0.13);           // placeholder
        addBR(DSSD2STAR, {DPLUS, PIMINUS}, 0.40);       // placeholder
        addBR(DSSD2STAR, {D0, PI0}, 0.20);              // placeholder
        addBR(DSSD2STARMINUS, {-DSTAR, PIMINUS}, 0.26); // placeholder
        addBR(DSSD2STARMINUS, {DSTARMINUS, PI0}, 0.13); // placeholder
        addBR(DSSD2STARMINUS, {-D0, PIMINUS}, 0.40);    // placeholder
        addBR(DSSD2STARMINUS, {DMINUS, PI0}, 0.20);     // placeholder

        addBR(JPSI, {ELECTRON, POSITRON}, 0.00597);
        addBR(JPSI, {MUON, ANTIMUON}, 0.00596);

        _names["Tau"] = {TAU, ANTITAU};
        _names["Nu"] = {NU_E, NU_EBAR, NU_MU, NU_MUBAR, NU_TAU, NU_TAUBAR};
        _names["Lep"] = {TAU, ANTITAU, MUON, ANTIMUON, ELECTRON, POSITRON};
        _names["Ell"] = {MUON, ANTIMUON, ELECTRON, POSITRON};
        _names["E"] = {POSITRON, ELECTRON};
        _names["Mu"] = {MUON, ANTIMUON};
        _names["Gamma"] = {GAMMA};

        _names["D"] = {DPLUS, DMINUS, D0, -D0};
        _names["D*"] = {DSTAR, DSTARMINUS, DSTARPLUS, -DSTAR};
        _names["D**0*"] = {DSSD0STAR, DSSD0STARMINUS, DSSD0STARPLUS, -DSSD0STAR};
        _names["D**1*"] = {DSSD1STAR, DSSD1STARMINUS, DSSD1STARPLUS, -DSSD1STAR};
        _names["D**1"] = {DSSD1, DSSD1MINUS, DSSD1PLUS, -DSSD1};
        _names["D**2*"] = {DSSD2STAR, DSSD2STARMINUS, DSSD2STARPLUS, -DSSD2STAR};
        _names["D**"] = {DSSD0STAR,     DSSD0STARMINUS, DSSD0STARPLUS, -DSSD0STAR, DSSD1STAR, DSSD1STARMINUS,
                         DSSD1STARPLUS, -DSSD1STAR,     DSSD1,         DSSD1MINUS, DSSD1PLUS, -DSSD1,
                         DSSD2STAR,     DSSD2STARMINUS, DSSD2STARPLUS, -DSSD2STAR};

        _names["Ds"] = {DSPLUS, DSMINUS};
        _names["Ds*"] = {DSSTARMINUS, DSSTARPLUS};
        _names["Ds**0*"] = {DSSDS0STARMINUS, DSSDS0STARPLUS};
        _names["Ds**1*"] = {DSSDS1STARMINUS, DSSDS1STARPLUS};
        _names["Ds**1"] = {DSSDS1MINUS, DSSDS1PLUS};
        _names["Ds**2*"] = {DSSDS2STARMINUS, DSSDS2STARPLUS};
        _names["Ds**"] = {DSSDS0STARMINUS, DSSDS0STARPLUS, DSSDS1STARMINUS, DSSDS1STARPLUS,
                          DSSDS1MINUS,     DSSDS1PLUS,     DSSDS2STARMINUS, DSSDS2STARPLUS};

        _names["B"] = {BZERO, BMINUS, BPLUS, -BZERO};
        _names["Bs"] = {BS, -BS};
        _names["Bc"] = {BCPLUS, BCMINUS};
        _names["Lb"] = {LAMBDAB, -LAMBDAB};
        _names["Lc"] = {LAMBDACPLUS, LAMBDACMINUS};
        _names["Lc*2595"] = {LAMBDACSTAR12PLUS, LAMBDACSTAR12MINUS};
        _names["Lc*2625"] = {LAMBDACSTAR32PLUS, LAMBDACSTAR32MINUS};
        _names["Jpsi"] = {JPSI};
        _names["K"] = {KPLUS, KMINUS, K0L, K0S, K0, -K0};
        _names["Pi"] = {PI0, PIPLUS, PIMINUS};
        _names["Rho"] = {RHO0, RHOPLUS, RHOMINUS};

        _names["Tau+"] = {ANTITAU};
        _names["Tau-"] = {TAU};
        _names["E+"] = {POSITRON};
        _names["E-"] = {ELECTRON};
        _names["Mu+"] = {ANTIMUON};
        _names["Mu-"] = {MUON};

        _names["Rho0"] = {RHO0};
        _names["Rho+"] = {RHOPLUS};
        _names["Rho-"] = {RHOMINUS};
        _names["Rho+-"] = {RHOPLUS, RHOMINUS};
        _names["Omega"] = {OMEGA};

        _names["K+"] = {KPLUS};
        _names["K-"] = {KMINUS};
        _names["K+-"] = {KPLUS, KMINUS};
        _names["K0S"] = {K0S};
        _names["K0L"] = {K0L};
        _names["K0"] = {K0};
        _names["K0bar"] = {-K0};
        _names["Kabar"] = {K0, -K0};

        _names["D0"] = {D0};
        _names["D+"] = {DPLUS};
        _names["D-"] = {DMINUS};
        _names["D0bar"] = {-D0};
        _names["D+-"] = {DPLUS, DMINUS};
        _names["Dabar"] = {D0, -D0};
        _names["D*0"] = {DSTAR};
        _names["D*+"] = {DSTARPLUS};
        _names["D*-"] = {DSTARMINUS};
        _names["D*0bar"] = {-DSTAR};
        _names["D*+-"] = {DSTARPLUS, DSTARMINUS};
        _names["D*abar"] = {DSTAR, -DSTAR};

        _names["D**0*0"] = {DSSD0STAR};
        _names["D**0*+"] = {DSSD0STARPLUS};
        _names["D**0*-"] = {DSSD0STARMINUS};
        _names["D**0*0bar"] = {-DSSD0STAR};
        _names["D**0*+-"] = {DSSD0STARPLUS, DSSD0STARMINUS};
        _names["D**0*abar"] = {DSSD0STAR, -DSSD0STAR};
        _names["D**1*0"] = {DSSD1STAR};
        _names["D**1*+"] = {DSSD1STARPLUS};
        _names["D**1*-"] = {DSSD1STARMINUS};
        _names["D**1*0bar"] = {-DSSD1STAR};
        _names["D**1*+-"] = {DSSD1STARPLUS, DSSD1STARMINUS};
        _names["D**1*abar"] = {DSSD1STAR, -DSSD1STAR};
        _names["D**10"] = {DSSD1};
        _names["D**1+"] = {DSSD1PLUS};
        _names["D**1-"] = {DSSD1MINUS};
        _names["D**10bar"] = {-DSSD1};
        _names["D**1+-"] = {DSSD1PLUS, DSSD1MINUS};
        _names["D**1abar"] = {DSSD1, -DSSD1};
        _names["D**2*0"] = {DSSD2STAR};
        _names["D**2*+"] = {DSSD2STARPLUS};
        _names["D**2*-"] = {DSSD2STARMINUS};
        _names["D**2*0bar"] = {-DSSD2STAR};
        _names["D**2*+-"] = {DSSD2STARPLUS, DSSD2STARMINUS};
        _names["D**2*abar"] = {DSSD2STAR, -DSSD2STAR};

        _names["Lc+"] = {LAMBDACPLUS};
        _names["Lc-"] = {LAMBDACMINUS};
        _names["Lc*2595+"] = {LAMBDACSTAR12PLUS};
        _names["Lc*2595-"] = {LAMBDACSTAR12MINUS};
        _names["Lc*2625+"] = {LAMBDACSTAR32PLUS};
        _names["Lc*2625-"] = {LAMBDACSTAR32MINUS};

        _names["Ds+"] = {DSPLUS};
        _names["Ds-"] = {DSMINUS};
        _names["Ds+-"] = {DSPLUS, DSMINUS};
        _names["Ds*+"] = {DSSTARPLUS};
        _names["Ds*-"] = {DSSTARMINUS};
        _names["Ds*+-"] = {DSSTARPLUS, DSSTARMINUS};

        _names["Ds**0*+"] = {DSSDS0STARPLUS};
        _names["Ds**0*-"] = {DSSDS0STARMINUS};
        _names["Ds**0*+-"] = {DSSDS0STARPLUS, DSSDS0STARMINUS};
        _names["Ds**1*+"] = {DSSDS1STARPLUS};
        _names["Ds**1*-"] = {DSSDS1STARMINUS};
        _names["Ds**1*+-"] = {DSSDS1STARPLUS, DSSDS1STARMINUS};
        _names["Ds**1+"] = {DSSDS1PLUS};
        _names["Ds**1-"] = {DSSDS1MINUS};
        _names["Ds**1+-"] = {DSSDS1PLUS, DSSDS1MINUS};
        _names["Ds**2*+"] = {DSSDS2STARPLUS};
        _names["Ds**2*-"] = {DSSDS2STARMINUS};
        _names["Ds**2*+-"] = {DSSDS2STARPLUS, DSSDS2STARMINUS};

        _names["B0"] = {BZERO};
        _names["B+"] = {BPLUS};
        _names["B-"] = {BMINUS};
        _names["B0bar"] = {-BZERO};
        _names["B+-"] = {BPLUS, BMINUS};
        _names["Babar"] = {BZERO, -BZERO};
        _names["Bs0"] = {BS};
        _names["Bs0bar"] = {-BS};
        _names["Bc+"] = {BCPLUS};
        _names["Bc-"] = {BCMINUS};
        _names["Lb0"] = {LAMBDAB};
        _names["Lb0bar"] = {-LAMBDAB};
        _names["Pi0"] = {PI0};
        _names["Pi+"] = {PIPLUS};
        _names["Pi-"] = {PIMINUS};
        _names["Pi+-"] = {PIPLUS, PIMINUS};
        _names["Nut"] = {NU_TAU};
        _names["Nutbar"] = {NU_TAUBAR};
        _names["Num"] = {NU_MU};
        _names["Numbar"] = {NU_MUBAR};
        _names["Nue"] = {NU_E};
        _names["Nuebar"] = {NU_EBAR};
        _names["W+"] = {WBOSON};
        _names["W-"] = {-WBOSON};
        _names["W"] = {WBOSON, -WBOSON};
    }

    unsigned short PID::digit(location loc, PdgId pid) {
        //  PID digits (base 10) are: n nr nl nq1 nq2 nq3 nj
        //  the location enum provides a convenient index into the PID
        int numerator = static_cast<int>(pow(10.0, (loc - 1)));
        return static_cast<unsigned short>((abspid(pid) / numerator) % 10);
    }

    PdgId PID::abspid(PdgId pid) {
        return (pid < 0) ? -pid : pid;
    }

    PdgId PID::extraBits(PdgId pid) {
        return abspid(pid) / 10000000;
    }

    PdgId PID::fundamentalID(PdgId pid) {
        if (extraBits(pid) > 0) {
            return 0;
        }

        if (digit(nq2, pid) == 0 && digit(nq1, pid) == 0) {
            return abspid(pid) % 10000;
        }
        if (abspid(pid) <= 100) {
            return abspid(pid);
        }
        return 0;
    }

    bool PID::isMeson(PdgId pid) {
        if (extraBits(pid) > 0) {
            return false;
        }

        if (abspid(pid) <= 100) {
            return false;
        }

        if (fundamentalID(pid) <= 100 && fundamentalID(pid) > 0) {
            return false;
        }

        int aid = abspid(pid);

        if (aid == 130 || aid == 310) { //|| aid == 111 || aid == 113 ) {
            return (pid >= 0);
        }

        // EvtGen uses some odd numbers
        if (aid == 150 || aid == 350 || aid == 510 || aid == 530) {
            return true;
        }

        // pomeron, etc.
        if (pid == 110 || pid == 990 || pid == 9990) {
            return true;
        }

        if (digit(nj, pid) > 0 && digit(nq3, pid) > 0 && digit(nq2, pid) > 0 && digit(nq1, pid) == 0) {
            // check for illegal antiparticles
            return (digit(nq3, pid) != digit(nq2, pid) || pid >= 0);
        }

        return false;
    }

    //  check to see if this is a valid baryon
    bool PID::isBaryon(PdgId pid) {
        if (extraBits(pid) > 0) {
            return false;
        }

        if (abspid(pid) <= 100) {
            return false;
        }

        if (fundamentalID(pid) <= 100 && fundamentalID(pid) > 0) {
            return false;
        }

        if (abspid(pid) == 2110 || abspid(pid) == 2210) {
            return true;
        }

        if (digit(nj, pid) > 0 && digit(nq3, pid) > 0 && digit(nq2, pid) > 0 && digit(nq1, pid) > 0) {
            return true;
        }

        return false;
    }

    //  check to see if this is a valid diquark
    bool PID::isDiQuark(PdgId pid) {
        if (extraBits(pid) > 0) {
            return false;
        }

        if (abspid(pid) <= 100) {
            return false;
        }

        if (fundamentalID(pid) <= 100 && fundamentalID(pid) > 0) {
            return false;
        }

        if (digit(nj, pid) > 0 && digit(nq3, pid) == 0 && digit(nq2, pid) > 0 &&
            digit(nq1, pid) > 0) { // diquark signature
            // EvtGen uses the diquarks for quark pairs, so, for instance,
            //   5501 is a valid "diquark" for EvtGen
            // if( digit(nj) == 1 && digit(nq2) == digit(nq1) ) {     // illegal
            //   return false;
            //} else {
            return true;
            //}
        }

        return false;
    }

    // is this a valid hadron ID?
    bool PID::isHadron(PdgId pid) {
        if (extraBits(pid) > 0) {
            return false;
        }

        if (isMeson(pid)) {
            return true;
        }

        if (isBaryon(pid)) {
            return true;
        }

        if (isPentaquark(pid)) {
            return true;
        }

        return false;
    }
    // is this a valid lepton ID?
    bool PID::isLepton(PdgId pid) {
        if (extraBits(pid) > 0) {
            return false;
        }

        if (fundamentalID(pid) >= 11 && fundamentalID(pid) <= 18) {
            return true;
        }

        return false;
    }
    // is this a neutrino ID?
    bool PID::isNeutrino(PdgId pid) {
        if (extraBits(pid) > 0) {
            return false;
        }

        if (fundamentalID(pid) == 12 || fundamentalID(pid) == 14 || fundamentalID(pid) == 16) {
            return true;
        }

        return false;
    }

    //  check to see if this is a valid pentaquark
    bool PID::isPentaquark(PdgId pid) {
        // a pentaquark is of the form 9abcdej,
        // where j is the spin and a, b, c, d, and e are quarks
        if (extraBits(pid) > 0) {
            return false;
        }

        if (digit(n, pid) != 9) {
            return false;
        }

        if (digit(nr, pid) == 9 || digit(nr, pid) == 0) {
            return false;
        }

        if (digit(nj, pid) == 9 || digit(nl, pid) == 0) {
            return false;
        }

        if (digit(nq1, pid) == 0) {
            return false;
        }

        if (digit(nq2, pid) == 0) {
            return false;
        }

        if (digit(nq3, pid) == 0) {
            return false;
        }

        if (digit(nj, pid) == 0) {
            return false;
        }

        // check ordering
        if (digit(nq2, pid) > digit(nq1, pid)) {
            return false;
        }

        if (digit(nq1, pid) > digit(nl, pid)) {
            return false;
        }

        if (digit(nl, pid) > digit(nr, pid)) {
            return false;
        }

        return true;
    }

    vector<pair<PdgId, vector<PdgId>>> PID::expandToValidVertices(const string& name) const {
        regex regexpression("[A-Z][^A-Z]*");
        regex_iterator<string::const_iterator> rit(name.begin(), name.end(), regexpression);
        regex_iterator<string::const_iterator> rend;
        vector<string> names{};
        while (rit != rend) {
            names.push_back(rit->str());
            ++rit;
        }
        vector<pair<PdgId, vector<PdgId>>> results;
        vector<vector<PdgId>> tempRes;
        for (const auto& elem : names) {
            auto res = toPdgList(elem);
            if (res.empty()) {
                return {};
            }
            if (tempRes.empty()) {
                transform(res.begin(), res.end(), back_inserter(tempRes),
                          [](PdgId elem2) -> vector<PdgId> { return {elem2}; });
            } else {
                vector<vector<PdgId>> tmpMat(tempRes.size() * res.size());
                for (size_t i = 0; i < tempRes.size(); ++i) {
                    for (size_t j = 0; j < res.size(); ++j) {
                        tmpMat[i * res.size() + j] = tempRes[i];
                        tmpMat[i * res.size() + j].push_back(res[j]);
                    }
                }
                tempRes.swap(tmpMat);
            }
        }
        for (auto& elemres : tempRes) {
            PdgId parent = elemres[0];
            if (getLeptonNumber(elemres) - 2 * getLeptonNumber(parent) != 0) {
                continue;
            }
            if (getBaryonNumber(elemres) - 2 * getBaryonNumber(parent) != 0) {
                continue;
            }
            if (getThreeCharge(elemres) - 2 * getThreeCharge(parent) != 0) {
                continue;
            }
            auto left = getLeptonFlavorNumber(elemres);
            auto right = getLeptonFlavorNumber(parent);
            if (get<0>(left) - 2 * get<0>(right) != 0 || get<1>(left) - 2 * get<1>(right) != 0 ||
                get<2>(left) - 2 * get<2>(right) != 0) {
                continue;
            }
            elemres.erase(elemres.begin());
            auto daughters = combineDaughters(elemres);
            bool found = false;
            for (auto elem : results) {
                if (elem.first != parent || elem.second.size() != daughters.size()) {
                    continue;
                }
                bool equal = true;
                for (size_t j = 0; j < daughters.size(); ++j) {
                    if (elem.second[j] != daughters[j]) {
                        equal = false;
                        break;
                    }
                }
                if (equal) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                results.emplace_back(parent, daughters);
                //                results.push_back(make_pair(flipSign(parent), flipSigns(daughters)));
            }
        }
        return results;
    }

    vector<HashId> PID::expandToValidVertexUIDs(const string& name, const bool& hadonly, const bool& flipped) const {
        auto tmp = expandToValidVertices(name);
        vector<HashId> result;
        result.reserve(tmp.size());
        for (const auto& elem : tmp) {
            if(!flipped){
                result.push_back(processID(elem.first, elem.second));
            } else {
                result.push_back(processID(flipSign(elem.first), flipSigns(elem.second)));
            }
        }
        if (hadonly) {
            auto tmpW = expandToValidVertices(name + "W");
            auto tmpN = expandToValidVertices(name + "Nu");
            tmpW.insert(tmpW.end(), tmpN.begin(), tmpN.end());
            for (const auto& elem : tmpW) {
                vector<PdgId> hadDaughters = elem.second;
                hadDaughters.erase(
                    remove_if(hadDaughters.begin(), hadDaughters.end(), [](PdgId idx) { return (abs(idx) < 100); }),
                    hadDaughters.end());
                //                elem.second.erase(remove(elem.second.begin(), elem.second.end(), 24),
                //                elem.second.end()); elem.second.erase(remove(elem.second.begin(), elem.second.end(),
                //                -24), elem.second.end());
                if(!flipped){
                    result.push_back(processID(elem.first, hadDaughters));
                } else {
                    result.push_back(processID(flipSign(elem.first), flipSigns(hadDaughters)));
                }
            }
        }
        return result;
    }


} // namespace Hammer
