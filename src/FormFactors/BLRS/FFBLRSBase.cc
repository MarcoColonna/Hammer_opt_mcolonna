///
/// @file  FFBLRSBase.cc
/// @brief Hammer base class for BLRS form factors
//

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/FormFactors/BLRS/FFBLRSBase.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Math/Constants.hh"

using namespace std;

namespace Hammer {

    FFBLRSBase::FFBLRSBase() {
        setGroup("BLRS");
    }

    void FFBLRSBase::addRefs() const {
        if (!getSettingsHandler()->checkReference("Bernlochner:2018kxh")) {
            string ref1 =
                "   @article{Bernlochner:2018kxh,\n"
                "   author         = \"Bernlochner, Florian U. and Ligeti, Zoltan and Robinson, Dean J. and Sutcliffe, "
                "William L.\",\n"
                "   title          = \"{New predictions for $\\Lambda_b\\to\\Lambda_c$ semileptonic decays and tests "
                "of heavy quark symmetry}\",\n"
                "   journal        = \"Phys. Rev. Lett.\",\n"
                "   volume         = \"121\",\n"
                "   year           = \"2018\",\n"
                "   number         = \"20\",\n"
                "   pages          = \"202001\",\n"
                "   doi            = \"10.1103/PhysRevLett.121.202001\",\n"
                "   eprint         = \"1808.09464\",\n"
                "   archivePrefix  = \"arXiv\",\n"
                "   primaryClass   = \"hep-ph\",\n"
                "   SLACcitation   = \"%%CITATION = ARXIV:1808.09464;%%\"\n"
                "}\n";
            getSettingsHandler()->addReference("Bernlochner:2018kxh", ref1);
        }
        if (!getSettingsHandler()->checkReference("Bernlochner:2018bfn")) {
            string ref2 =
                "@article{Bernlochner:2018bfn,\n"
                "   author         = \"Bernlochner, Florian U. and Ligeti, Zoltan and Robinson, Dean J. and Sutcliffe, "
                "William L.\",\n"
                "   title          = \"{Precise predictions for $\\Lambda_b \\to \\Lambda_c$ semileptonic decays}\",\n"
                "   journal        = \"Phys. Rev.\",\n"
                "   volume         = \"D99\",\n"
                "   year           = \"2019\",\n"
                "   number         = \"5\",\n"
                "   pages          = \"055008\",\n"
                "   doi            = \"10.1103/PhysRevD.99.055008\",\n"
                "   eprint         = \"1812.07593\",\n"
                "   archivePrefix  = \"arXiv\",\n"
                "   primaryClass   = \"hep-ph\",\n"
                "   SLACcitation   = \"%%CITATION = ARXIV:1812.07593;%%\"\n"
                "}\n";
            getSettingsHandler()->addReference("Bernlochner:2018bfn", ref2);
        }
    }

    void FFBLRSBase::defineSettings() {
        // 1S scheme: mb1S = 4721.5, delta mb-mc = 3400, alpha_s = 26/100
        addSetting<double>("as", 0.26);
        addSetting<double>("mb1S", 4.7215); // GeV
        addSetting<double>("dmbmc", 3.400); // GeV
        addSetting<double>("mLb", 5.620);   // GeV
        addSetting<double>("mLc", 2.286);   // GeV

        addSetting<double>("z1", -2.037);
        addSetting<double>("z2", 3.159);
        addSetting<double>("b1", -0.459); // GeV^2
        addSetting<double>("b2", -0.390); // GeV^2
    }

    void FFBLRSBase::fillKisMis(const double& unitres) {
        auto unused = [](double) -> double { return 0.; };

        _b1 = unitres * unitres * (*getSetting<double>("b1"));
        _b2 = unitres * unitres * (*getSetting<double>("b2"));

        // NLO
        _mKi1[0] = unused;
        _mKi1[1] = unused;
        _mKi1[2] = [](double w) -> double { return -1. / (1 + w); };

        // NNLO
        _mKi2[0] = unused;
        _mKi2[1] = [this](double w) -> double { return _b1 / (_la1S * _la1S) + (w - 1.) * _b2 / (2. * _la1S * _la1S); };
        _mKi2[2] = [this](double) -> double { return _b2 / (2. * _la1S * _la1S); };
        //        for(size_t i = 0; i < 5; ++i) {
        //            _mMi[i] = unused;
        //        }
    }

    void FFBLRSBase::calcConstants(const double& unitres) {
        // Renormalon improvement scheme
        const double aS = (*getSetting<double>("as")) / pi;
        const double mb1S = unitres * (*getSetting<double>("mb1S"));
        const double dmbmc = unitres * (*getSetting<double>("dmbmc"));
        const double mLb = unitres * (*getSetting<double>("mLb"));
        const double mLc = unitres * (*getSetting<double>("mLc"));
        const bool changed = !(isZero(_aS - aS) && isZero(_mb1S - mb1S) && isZero(_dmbmc - dmbmc) &&
                               isZero(_mLb - mLb) && isZero(_mLc - mLc));
        if (changed) {
            _aS = aS;
            _mb1S = mb1S;
            _dmbmc = dmbmc;
            _mLb = mLb;
            _mLc = mLc;
            const double mc1S = mb1S - dmbmc;
            _mc1S = mc1S;
            const double la1S = (mb1S * (mLb - mb1S) - mc1S * (mLc - mc1S)) / (mb1S - mc1S);
            _la1S = la1S;
            _zBC = mc1S / mb1S;

            _eB0 = la1S / (2. * mb1S);
            _eC0 = la1S / (2. * mc1S);
            const double eps = 2. * pow(aS * pi, 2.) / 9.;
            // Explicit expansion to O(eps)
            _eB = _eB0 + eps * (dmbmc * mc1S - mb1S * mLb + mc1S * mLc) / (mb1S - mc1S) / (2. * mb1S);
            _eC =
                _eC0 + eps * ((mb1S * (dmbmc * mb1S - mb1S * mLb + mc1S * mLc)) / (mc1S * (mb1S - mc1S))) / (2. * mc1S);
        }
    }

} // namespace Hammer
