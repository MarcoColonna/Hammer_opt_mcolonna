///
/// @file  FFBLPRXPBase.cc
/// @brief Hammer base class for BLPR form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/FormFactors/BLPRXP/FFBLPRXPBase.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Math/Constants.hh"

using namespace std;

namespace Hammer {

    FFBLPRXPBase::FFBLPRXPBase() {
        setGroup("BLPRXP");
    }

    void FFBLPRXPBase::addRefs() const {
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
        if (!getSettingsHandler()->checkReference("Bernlochner:2022ywh")) {
            string ref =
                "@article{Bernlochner:2022ywh,\n"
                "    author = \"Bernlochner, Florian U. and Ligeti, Zoltan and Papucci, Michele and Prim, Markus T. "
                "and Robinson, Dean J. and Xiong, Chenglu\",\n"
                "    title = \"{Constrained second-order power corrections in HQET: $R(D^{(*)})$, $|V_{cb}|$, and new "
                "physics}\",\n"
                "    eprint = \"2206.11281\",\n"
                "    archivePrefix = \"arXiv\",\n"
                "    primaryClass = \"hep-ph\",\n"
                "    reportNumber = \"CALT-TH-2022-022\",\n"
                "    month = \"6\",\n"
                "    year = \"2022\"\n"
                "}\n";
            getSettingsHandler()->addReference("Bernlochner:2022ywh", ref);
        }
    }

    void FFBLPRXPBase::defineSettings() {
        addSetting<double>("RhoStSq", 1.104);
        addSetting<double>("a", 1.509 / sqrt2);
        addSetting<double>("cSt", 2.392);

        // 1S scheme: mb1S = 4710, mbBar = 5313, lambda1 = -3 10^5 MeV^2, delta mb-mc = 3400, alpha_s = 26/100
        addSetting<double>("as", 0.27);
        addSetting<double>("mb", 4.707);
        addSetting<double>("mc", 4.707 - 3.407);
        addSetting<double>("mBBar", 5.313);
        addSetting<double>("mDBar", 1.973);
        addSetting<double>("rho1", -0.363);
        addSetting<double>("la2", 0.12);

        addSetting<double>("chi21", -0.116);
        addSetting<double>("chi2p", 0.0);
        addSetting<double>("chi3p", 0.0);
        addSetting<double>("eta1", 0.336);
        addSetting<double>("etap", 0.0);
        addSetting<double>("phi1p", 0.252);
        addSetting<double>("beta21", 0.0);
        addSetting<double>("beta3p", 0.0);

        // Switches
        addSetting<bool>("With1OverMb2", true);
        addSetting<bool>("With1OverMbMc", true);
        addSetting<bool>("WithAsOverM", true);

        fillIWs();
        fillLisMis();
    }

    void FFBLPRXPBase::fillIWs() {
        _mIWs[CHI2] = [this](double w) -> double {
            return (*getSetting<double>("chi21")) + (w - 1.) * (*getSetting<double>("chi2p"));
        };
        _mIWs[CHI3] = [this](double w) -> double { return (*getSetting<double>("chi3p")) * (w - 1.); };
        _mIWs[ETA] = [this](double w) -> double {
            return (*getSetting<double>("eta1")) + (*getSetting<double>("etap")) * (w - 1.);
        };
        _mIWs[BETA1] = [this](double) -> double { return _la1OverlaB2 / 4.; };
        _mIWs[BETA2] = [this](double) -> double { return (*getSetting<double>("beta21")); };
        _mIWs[BETA3] = [this](double w) -> double {
            return _la2OverlaB2 / 8. + (*getSetting<double>("beta3p")) * (w - 1.);
        };
        _mIWs[PHI1] = [this](double w) -> double {
            return (_la1OverlaB2 / 3. - _la2OverlaB2 / 2.) / 2. + (*getSetting<double>("phi1p")) * (w - 1.);
        };
        _mIWs[PHI1Q] = [this](double) -> double { return (*getSetting<double>("phi1p")); };
    }

    void FFBLPRXPBase::fillLisMis() {
        auto unused = [](double) -> double { return 0.; };

        // 1/M
        _mLi1[0] = unused;
        _mLi1[1] = [this](double w) -> double { return 4. * (3. * _mIWs[CHI3](w) - (w - 1.) * _mIWs[CHI2](w)); };
        _mLi1[2] = [this](double w) -> double { return -4. * _mIWs[CHI3](w); };
        _mLi1[3] = [this](double w) -> double { return 4. * _mIWs[CHI2](w); };
        _mLi1[4] = [this](double w) -> double { return 2. * _mIWs[ETA](w) - 1.; };
        _mLi1[5] = [](double) -> double { return -1.; };
        _mLi1[6] = [this](double w) -> double { return -2. * (_mIWs[ETA](w) + 1.) / (w + 1.); };

        // 1/M^2
        _mLi2[0] = unused;
        for (size_t i = 0; i < 8; ++i) {
            _mMi[i] = unused;
        }
        for (size_t i = 11; i < 25; ++i) {
            _mMi[i] = unused;
        }

        _mLi2[1] = [this](double w) -> double {
            return 2. * _mIWs[BETA1](w) + 4. * (3. * _mIWs[BETA3](w) - (w - 1.) * _mIWs[BETA2](w));
        };
        _mLi2[2] = [this](double w) -> double { return 2. * _mIWs[BETA1](w) - 4. * _mIWs[BETA3](w); };
        _mLi2[3] = [this](double w) -> double { return 4. * _mIWs[BETA2](w); };
        _mLi2[4] = [this](double w) -> double { return 3 * _la2OverlaB2 + 2. * (w + 1.) * _mIWs[PHI1](w); };
        _mLi2[5] = [this](double w) -> double { return _la2OverlaB2 + 2. * (w + 1.) * _mIWs[PHI1](w); };
        _mLi2[6] = [this](double w) -> double { return 4. * _mIWs[PHI1](w); };
        _mMi[8] = [this](double w) -> double {
            return (_la1OverlaB2 + 6. * _la2OverlaB2 / (w + 1.) - 2. * (w - 1.) * _mIWs[PHI1](w) -
                    2. * (2. * _mIWs[ETA](w) - 1.) * (w - 1.) / (w + 1.));
        };
        _mMi[9] = [this](double w) -> double {
            return 3. * _la2OverlaB2 / (w + 1.) + 2. * _mIWs[PHI1](w) - (2. * _mIWs[ETA](w) - 1.) * (w - 1.) / (w + 1.);
        };
        _mMi[10] = [this](double w) -> double {
            return (_la1OverlaB2 / 3. - _la2OverlaB2 * (w + 4.) / (2. * (w + 1.)) + 2. * (w + 2.) * _mIWs[PHI1Q](w) -
                    (2. * _mIWs[ETA](w) - 1.) / (w + 1.));
        };
    }

    double FFBLPRXPBase::ZofW(double w) const {
        return (sqrt(w + 1) - sqrt2 * _a) / (sqrt(w + 1) + sqrt2 * _a);
    }

    double FFBLPRXPBase::Xi(double w) const {
        // Variables for leading IW function (derived from G(1))
        const double cSt = *getSetting<double>("cSt");
        const double RhoSq = *getSetting<double>("RhoStSq");
        const double a2 = _a * _a;
        const double a4 = a2 * a2;
        const double zCon = ZofW(w);
        const double zCon1 = (1 - _a) / (1 + _a);
        const double XiVal = (1. - 8 * a2 * RhoSq * zCon + 16. * (2 * cSt * a4 - RhoSq * a2) * zCon * zCon) /
                             (1. - 8 * a2 * RhoSq * zCon1 + 16. * (2 * cSt * a4 - RhoSq * a2) * zCon1 * zCon1);
        return XiVal;
    }

    void FFBLPRXPBase::calcConstants() {
        const double a = *getSetting<double>("a");
        const double mb = *getSetting<double>("mb");
        const double mc = *getSetting<double>("mc");
        const double la2 = *getSetting<double>("la2");
        const double rho1 = *getSetting<double>("rho1");
        const double mBBar = *getSetting<double>("mBBar");
        const double mDBar = *getSetting<double>("mDBar");
        const bool changed = !(isZero(_a - a) && isZero(_mb - mb) && isZero(_mc - mc) && isZero(_la2 - la2) &&
                               isZero(_rho1 - rho1) && isZero(_mBBar - mBBar) && isZero(_mDBar - mDBar));
        if (changed) {
            _a = a;
            _mb = mb;
            _mc = mc;
            _mBBar = mBBar;
            _mDBar = mDBar;
            _la2 = la2;
            _rho1 = rho1;
            _zBC = _mc / _mb;
            _cmagb = 3. / 2. * (0.5 * log(_zBC) + 13. / 9.);
            _cmagc = -3. / 2. * (0.5 * log(_zBC) - 13. / 9.);
            const double laB = (_mb * _mBBar - _mc * _mDBar) / (_mb - _mc) - (_mb + _mc) + _rho1 / (4. * _mb * _mc);
            const double la1 = (2. * _mb * _mc) / (_mb - _mc) * (_mBBar - _mDBar - (_mb - _mc)) +
                               _rho1 * (_mb + _mc) / (2. * _mb * _mc);
            _eb = laB / (2. * _mb);
            _ec = laB / (2. * _mc);
            _la2OverlaB2 = _la2 / (laB * laB);
            _la1OverlaB2 = la1 / (laB * laB);
            // upsilon corrections
            const double corr1S =
                2. * pow((*getSetting<double>("as")) / 3. * 0.796, 2.); // runned up from sqrt(mb mc) to mb1S
            const double dmbc = _mb - _mc;
            const double mb_0 = _mb;
            const double mb_1 = _mb * corr1S;
            const double mc_0 = _mc;
            const double mc_1 = _mb * corr1S;
            const double LambdaBar_0 =
                (mb_0 * _mBBar - mc_0 * _mDBar) / dmbc - (2. * mb_0 - dmbc) + _rho1 / (4 * mb_0 * mc_0);
            const double LambdaBar_1 = (mb_1 * _mBBar - mc_1 * _mDBar) / dmbc - 2. * mb_1;
            _upsilonb = LambdaBar_1 / LambdaBar_0 - mb_1 / mb_0;
            _upsilonc = LambdaBar_1 / LambdaBar_0 - mc_1 / mc_0;
        }
    }

} // namespace Hammer
