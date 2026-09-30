///
/// @file  FFBLRSXPBase.cc
/// @brief Hammer base class for BLRSXP form factors
//

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/FormFactors/BLRSXP/FFBLRSXPBase.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Math/Constants.hh"

using namespace std;

namespace Hammer {

    FFBLRSXPBase::FFBLRSXPBase() {
        setGroup("BLRSXP");
    }

    void FFBLRSXPBase::addRefs() const {
        if (!getSettingsHandler()->checkReference("Bernlochner:2022ywh")) {
            string ref1 =
                "@article{Bernlochner:2022ywh,\n"
                "    author = \"Bernlochner, Florian U. and Ligeti, Zoltan and Papucci, Michele and Prim, Markus T. "
                "and Robinson, Dean J. and Xiong, Chenglu\",\n"
                "    title = \"{Constrained second-order power corrections in HQET: R(D(*)), |Vcb|, and new "
                "physics}\",\n"
                "    eprint = \"2206.11281\",\n"
                "    archivePrefix = \"arXiv\",\n"
                "    primaryClass = \"hep-ph\",\n"
                "    reportNumber = \"CALT-TH-2022-022\",\n"
                "    doi = \"10.1103/PhysRevD.106.096015\",\n"
                "    journal = \"Phys. Rev. D\",\n"
                "    volume = \"106\",\n"
                "    number = \"9\",\n"
                "    pages = \"096015\",\n"
                "    year = \"2022\"\n"
                "}\n";
            getSettingsHandler()->addReference("Bernlochner:2022ywh", ref1);
        }
        if (!getSettingsHandler()->checkReference("Bernlochner:2023jkp")) {
            string ref2 =
                "@article{Bernlochner:2023jkp,\n"
                "    author = \"Bernlochner, Florian U. and Papucci, Michele and Robinson, Dean J.\",\n"
                "    title = \"{Appraising constrained second-order power corrections in HQET with $\\Lambda_b \\to "
                "\\Lambda_c l \\nu$}\",\n"
                "    eprint = \"2312.07758\",\n"
                "    archivePrefix = \"arXiv\",\n"
                "    primaryClass = \"hep-ph\",\n"
                "    reportNumber = \"CALT-TH-2023-049\",\n"
                "    month = \"12\",\n"
                "    year = \"2023\"\n"
                "}\n";
            getSettingsHandler()->addReference("Bernlochner:2023jkp", ref2);
        }
    }

    void FFBLRSXPBase::defineSettings() {
        // 1S scheme
        addSetting<double>("as", 0.27);
        addSetting<double>("asmb", 0.215);   // as(mb)
        addSetting<double>("mb1S", 4.7241);  // GeV
        addSetting<double>("dmbmc", 3.4056); // GeV
        addSetting<double>("mLb", 5.620);    // GeV
        addSetting<double>("mLc", 2.286);    // GeV

        addSetting<double>("z1", -2.1638);
        addSetting<double>("z2", 3.3683);
        addSetting<double>("ph1p", 0.092666);
        addSetting<double>("ph1pp", -4.8788);
        addSetting<double>("rho1", -0.17446); // GeV^3
    }

    void FFBLRSXPBase::fillKisMis(const double& /*unused*/) {
        auto unused = [](double) -> double { return 0.; };

        _ph1p = (*getSetting<double>("ph1p"));
        _ph1pp = (*getSetting<double>("ph1pp"));

        _mPHI1N = [this](double w) -> double { return _ph1p + 0.5 * (w - 1) * _ph1pp; };

        // NLO
        _mKi1[0] = unused;
        _mKi1[1] = unused;
        _mKi1[2] = [](double w) -> double { return -1. / (1 + w); };

        // NNLO
        _lam1onla1S2 = _lam1 / (_la1S * _la1S);
        _mKi2[0] = unused;
        _mKi2[1] = [this](double) -> double { return _lam1onla1S2 / 2; };
        _mKi2[2] = [this](double w) -> double { return _lam1onla1S2 / 3 + 2 * (w - 1) * _mPHI1N(w); };
        for (size_t i = 0; i < 3; ++i) {
            _mMi[i] = unused;
        }
        _mMi[3] = [this](double w) -> double {
            return ((w - 2) * _lam1onla1S2 / 6 + (w * w + 2) * _mPHI1N(w) - (w - 2) / (2. * (w + 1))) / (w + 1);
        };
        _mMi[4] = [this](double w) -> double {
            return (3 - w) * _lam1onla1S2 / 6 - w * (w - 1) * _mPHI1N(w) + (w - 1) / (2. * (w + 1));
        };
    }

    void FFBLRSXPBase::calcConstants(const double& unitres) {
        // Renormalon improvement scheme
        const double aS = (*getSetting<double>("as")) / pi;
        const double aSmb = (*getSetting<double>("asmb"));
        const double mb1S = unitres * (*getSetting<double>("mb1S"));
        const double dmbmc = unitres * (*getSetting<double>("dmbmc"));
        const double mLb = unitres * (*getSetting<double>("mLb"));
        const double mLc = unitres * (*getSetting<double>("mLc"));
        const double rho1 = unitres * unitres * unitres * (*getSetting<double>("rho1"));
        const bool changed =
            !(isZero(_aS - aS) && isZero(_aSmb - aSmb) && isZero(_mb1S - mb1S) && isZero(_dmbmc - dmbmc) &&
              isZero(_mLb - mLb) && isZero(_mLc - mLc) && isZero(_rho1 - rho1));
        if (changed) {
            _aS = aS;
            _aSmb = aSmb;
            _mb1S = mb1S;
            _dmbmc = dmbmc;
            _mLb = mLb;
            _mLc = mLc;
            _rho1 = rho1;
            const double mc1S = mb1S - dmbmc;
            _mc1S = mc1S;
            const double la1S = (mb1S * mLb - mc1S * mLc) / dmbmc - (mb1S + mc1S) + rho1 / (4. * mc1S * mb1S);
            const double lam1 =
                2. * mb1S * mc1S * (mLb - mLc - dmbmc) / dmbmc + rho1 * (mb1S + mc1S) / (2. * mc1S * mb1S);
            _la1S = la1S;
            _lam1 = lam1;
            _zBC = mc1S / mb1S;

            _eB0 = la1S / (2. * mb1S);
            _eC0 = la1S / (2. * mc1S);
            const double eps = 2. * pow(aSmb, 2.) / 9.;
            // Explicit expansion to O(eps) N.B. different to BLRS bc of lam1/2mq(mb1S) term vs lam1/2mq1S
            _eB = _eB0 - eps * (mLc + dmbmc + rho1 * (mb1S + 2 * mc1S) / (4. * mb1S * mc1S * mc1S)) / (2. * mb1S);
            _eC = _eC0 - eps * ((mLb - dmbmc) * (mb1S / mc1S) + rho1 * (2 * mb1S + mc1S) / (4. * mb1S * mc1S * mc1S)) /
                             (2. * mc1S);
        }
    }

} // namespace Hammer
