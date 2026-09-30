///
/// @file  FFBLPRBase.cc
/// @brief Hammer base class for BLPR form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/FormFactors/BLPR/FFBLPRBase.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Math/Constants.hh"

using namespace std;

namespace Hammer {

    FFBLPRBase::FFBLPRBase() {
        setGroup("BLPR");
    }

    void FFBLPRBase::addRefs() const {
        if (!getSettingsHandler()->checkReference("Bernlochner:2017jka")) {
            string ref =
                "@article{Bernlochner:2017jka,\n"
                "      author         = \"Bernlochner, Florian U. and Ligeti, Zoltan and Papucci, Michele and "
                "Robinson, Dean J.\",\n"
                "      title          = \"{Combined analysis of semileptonic $B$ decays to $D$ and $D^*$: "
                "$R(D^{(*)})$, $|V_{cb}|$, and new physics}\",\n"
                "      journal        = \"Phys. Rev.\",\n"
                "      volume         = \"D95\",\n"
                "      year           = \"2017\",\n"
                "      number         = \"11\",\n"
                "      pages          = \"115008\",\n"
                "      doi            = \"10.1103/PhysRevD.95.115008, 10.1103/PhysRevD.97.059902\",\n"
                "      note           = \"[erratum: Phys. Rev.D97,no.5,059902(2018)]\",\n"
                "      eprint         = \"1703.05330\",\n"
                "      archivePrefix  = \"arXiv\",\n"
                "      primaryClass   = \"hep-ph\",\n"
                "      SLACcitation   = \"%%CITATION = ARXIV:1703.05330;%%\"\n"
                "}\n";
            getSettingsHandler()->addReference("Bernlochner:2017jka", ref);
        }
    }

    void FFBLPRBase::defineSettings() {
        addSetting<double>("RhoSq", 1.24);
        addSetting<double>("a", 1.509 / sqrt2);
        addSetting<double>("dV20", 0.);

        // 1S scheme: mb1S = 4710, mbBar = 5313, lambda1 = -3 10^5 MeV^2, delta mb-mc = 3400, alpha_s = 26/100
        addSetting<double>("as", 0.26);
        addSetting<double>("mb", 4.710);
        addSetting<double>("mc", 4.710 - 3.400);
        addSetting<double>("la", 0.57115);

        // Renormalon improvement (1S scheme)
        addSetting<double>("ebR/eb", 0.861);
        addSetting<double>("ecR/ec", 0.822);

        addSetting<double>("chi21", -0.06);
        addSetting<double>("chi2p", -0.00);
        addSetting<double>("chi3p", 0.05);
        addSetting<double>("eta1", 0.30);
        addSetting<double>("etap", -0.05);

        fillIWs();
        fillLis();
    }

    void FFBLPRBase::fillIWs() {
        _mIWs[ETA] = [this](double w) -> double {
            return (*getSetting<double>("eta1")) + (*getSetting<double>("etap")) * (w - 1.);
        };
        _mIWs[CHI2] = [this](double w) -> double {
            return (*getSetting<double>("chi21")) + (w - 1.) * (*getSetting<double>("chi2p"));
        };
        _mIWs[CHI3] = [this](double w) -> double { return (*getSetting<double>("chi3p")) * (w - 1.); };
    }

    void FFBLPRBase::fillLis() {
        auto unused = [](double) -> double { return 0.; };
        _mLi[0] = unused;
        _mLi[1] = [this](double w) -> double { return 4. * (3. * _mIWs[CHI3](w) - (w - 1.) * _mIWs[CHI2](w)); };
        _mLi[2] = [this](double w) -> double { return -4. * _mIWs[CHI3](w); };
        _mLi[3] = [this](double w) -> double { return 4. * _mIWs[CHI2](w); };
        _mLi[4] = [this](double w) -> double { return 2. * _mIWs[ETA](w) - 1.; };
        _mLi[5] = [](double) -> double { return -1.; };
        _mLi[6] = [this](double w) -> double { return -2. * (_mIWs[ETA](w) + 1.) / (w + 1.); };
    }

    double FFBLPRBase::ZofW(double w) const {
        return (sqrt(w + 1) - sqrt2 * _a) / (sqrt(w + 1) + sqrt2 * _a);
    }

    pair<double, double> FFBLPRBase::Xi(double w) const {
        // Variables for leading IW function (derived from G(1))
        const double V21 = 57.0;
        const double V20 = 7.5 + (*getSetting<double>("dV20"));
        const double RhoSq = *getSetting<double>("RhoSq");
        const double rD = 1867. / 5280.;
        const double eb = (*getSetting<double>("la")) / (*getSetting<double>("mb") * 2.);
        const double ec = (*getSetting<double>("la")) / (*getSetting<double>("mc") * 2.);
        const double as = (*getSetting<double>("as")) / pi;
        const double chi21 = (*getSetting<double>("chi21"));
        const double chi2p = (*getSetting<double>("chi2p"));
        const double chi3p = (*getSetting<double>("chi3p"));
        const double etap = (*getSetting<double>("etap"));
        const double zBC = (*getSetting<double>("mc")) / (*getSetting<double>("mb"));
        const double a = *getSetting<double>("a");
        const double a2 = a * a;
        const double a4 = a2 * a2;
        const double a6 = a4 * a2;
        // Derivatives (approx) at w_0 = 2a^2-1. Calculated once and stored.
        if (!_initCs || !isZero(a - _a) || !isZero(zBC - _zBC)) {
            _w0 = 2 * a2 - 1.;
            _a = a;
            _zBC = zBC;

            _mCv1z = _asCorrections.CV1(_w0, zBC);
            _mCv2z = _asCorrections.CV2(_w0, zBC);
            _mCv3z = _asCorrections.CV3(_w0, zBC);

            _mCv1zp = (_asCorrections.CV1(_w0 + 1e-5, zBC) - _mCv1z) / 1e-5;
            _mCv2zp = (_asCorrections.CV2(_w0 + 1e-5, zBC) - _mCv2z) / 1e-5;
            _mCv3zp = (_asCorrections.CV3(_w0 + 1e-5, zBC) - _mCv3z) / 1e-5;

            _mCv1zpp = (_mCv1zp - (_mCv1z - _asCorrections.CV1(_w0 - 1e-5, zBC)) / 1e-5) / 1e-5;
            _mCv2zpp = (_mCv2zp - (_mCv2z - _asCorrections.CV2(_w0 - 1e-5, zBC)) / 1e-5) / 1e-5;
            _mCv3zpp = (_mCv3zp - (_mCv3z - _asCorrections.CV3(_w0 - 1e-5, zBC)) / 1e-5) / 1e-5;

            _initCs = true;
        }
        const double zCon = ZofW(w);
        const double Xi = 64 * a4 * RhoSq - 16 * a2 - V21;
        double xiIW = 1 - 8 * a2 * RhoSq * zCon +
                      zCon * zCon *
                          (V21 * RhoSq - V20 + (eb - ec) * (2 * Xi * etap * (1 - rD) / (1 + rD)) +
                           (eb + ec) * (Xi * (12 * chi3p - 4 * chi21) - 16 * ((a2 - 1) * Xi - 16 * a4) * chi2p) +
                           as * (Xi * (_mCv1zp + (_mCv3z + rD * _mCv2z) / (1 + rD)) +
                                 2 * a2 * (Xi - 32 * a2) * (_mCv3zp + rD * _mCv2zp) / (1 + rD) -
                                 64 * a6 * (_mCv3zpp + rD * _mCv2zpp) / (1 + rD) - 32 * a4 * _mCv1zpp));
        double xiIWN = 1 - 8 * a2 * RhoSq * (1 - a) / (1 + a) +
                       pow((1 - a) / (1 + a), 2.) *
                           (V21 * RhoSq - V20 + (eb - ec) * (2 * Xi * etap * (1 - rD) / (1 + rD)) +
                            (eb + ec) * (Xi * (12 * chi3p - 4 * chi21) - 16 * ((a2 - 1) * Xi - 16 * a4) * chi2p) +
                            as * (Xi * (_mCv1zp + (_mCv3z + rD * _mCv2z) / (1 + rD)) +
                                  2 * a2 * (Xi - 32 * a2) * (_mCv3zp + rD * _mCv2zp) / (1 + rD) -
                                  64 * a6 * (_mCv3zpp + rD * _mCv2zpp) / (1 + rD) - 32 * a4 * _mCv1zpp));
        return {xiIW, xiIWN};
    }


} // namespace Hammer
