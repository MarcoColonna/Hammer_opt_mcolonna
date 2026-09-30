///
/// @file  FFBGLXBase.hh
/// @brief Hammer base class for BGL Extended form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/FormFactors/BGL/FFBGLXBase.hh"
#include "Hammer/Math/Constants.hh"

using namespace std;

namespace Hammer {

    FFBGLXBase::FFBGLXBase() {
        setGroup("BGLX");
    }

    void FFBGLXBase::addRefs() const {
        if (!getSettingsHandler()->checkReference("Boyd:1995sq")) {
            string ref1 =
                "@article{Boyd:1995sq,\n"
                "      author         = \"Boyd, C. Glenn and Grinstein, Benjamin and Lebed, Richard F.\",\n"
                "      title          = \"{Model independent determinations of anti-B ---> D (lepton), D* (lepton) "
                "anti-neutrino form-factors}\",\n"
                "      journal        = \"Nucl. Phys.\",\n"
                "      volume         = \"B461\",\n"
                "      year           = \"1996\",\n"
                "      pages          = \"493-511\",\n"
                "      doi            = \"10.1016/0550-3213(95)00653-2\",\n"
                "      eprint         = \"hep-ph/9508211\",\n"
                "      archivePrefix  = \"arXiv\",\n"
                "      primaryClass   = \"hep-ph\",\n"
                "      reportNumber   = \"UCSD-PTH-95-11\",\n"
                "      SLACcitation   = \"%%CITATION = HEP-PH/9508211;%%\"\n"
                "}\n";
            getSettingsHandler()->addReference("Boyd:1995sq", ref1);
        }
        if (!getSettingsHandler()->checkReference("Boyd:1997kz")) {
            string ref2 =
                "@article{Boyd:1997kz,\n"
                "      author         = \"Boyd, C. Glenn and Grinstein, Benjamin and Lebed, Richard F.\",\n"
                "      title          = \"{Precision corrections to dispersive bounds on form-factors}\",\n"
                "      journal        = \"Phys. Rev.\",\n"
                "      volume         = \"D56\",\n"
                "      year           = \"1997\",\n"
                "      pages          = \"6895-6911\",\n"
                "      doi            = \"10.1103/PhysRevD.56.6895\",\n"
                "      eprint         = \"hep-ph/9705252\",\n"
                "      archivePrefix  = \"arXiv\",\n"
                "      primaryClass   = \"hep-ph\",\n"
                "      reportNumber   = \"CMU-HEP-97-07A, UCSD-PTH-97-12\",\n"
                "      SLACcitation   = \"%%CITATION = HEP-PH/9705252;%%\"\n"
                "}\n";
            getSettingsHandler()->addReference("Boyd:1997kz", ref2);
        }
        if (!getSettingsHandler()->checkReference("Bordone:2025jur")) {
            string ref3 =
                "@article{Bordone:2025jur,\n"
                "    author = \"Bordone, Marzia and Gubernari, Nico and Jung, Martin and van Dyk, Danny\",\n"
                "    title = \"{Challenging $ {\\overline{B}}_{(s)}\\to {D}_{(s)}^{\\left(\\ast \\right)} $ form factors with the heavy quark expansion}\",\n"
                "    eprint = \"2507.03569\",\n"
                "    archivePrefix = \"arXiv\",\n"
                "    primaryClass = \"hep-ph\",\n"
                "    reportNumber = \"CERN-TH-2025-092, EOS-2025-03, IPPP/25/25, P3H-25-032, SI-HEP-2025-10, ZU-TH 35/25\",\n"
                "    doi = \"10.1007/JHEP11(2025)051\",\n"
                "    journal = \"JHEP\",\n"
                "    volume = \"11\",\n"
                "    pages = \"051\",\n"
                "    year = \"2025\"\n"
                "}\n";
            getSettingsHandler()->addReference("Bordone:2025jur", ref3);
        }
    }
    
//    void FFBGLXBase::defineSettings() {
//        // values from hep-ph/9705252, u = 0.33
//        addSetting<double>("Vcb", 41.5e-3);
//        
//        addSetting<double>("ChimT", 3.068e-4); // GeV^-2  f, F1
//        addSetting<double>("ChipT", 5.280e-4); // GeV^-2  g, fp
//        addSetting<double>("ChimL", 2.466e-3); // P1|F2
//        addSetting<double>("ChipL", 4.832e-3); // f0
//        
//        addSetting<double>("Chi1T", 4.98e-4 * 2.6/2.0); // GeV^-2  from 2507.03569 NB different matching scale and isospin factor
//        addSetting<double>("Chi1AT", 2.77e-4 * 2.6/2.0); // GeV^-2 from 2507.03569 NB different matching scale and isospin factor
//        
//        // pole mass from hep-ph/9705252
//        addSetting<double>("mb", 4.9);       // GeV
//        addSetting<double>("mc", 4.9 - 3.4); // GeV (using mb-mc = 3.4 GeV)
//        
//        vector<double> BcStates1P{6.730, 6.736, 7.135, 7.142};// GeV
//        vector<double> BcStates1M{6.337, 6.899, 7.012, 7.280};// GeV
//        vector<double> BcStates0M{6.275, 6.842, 7.250};       // GeV
//        vector<double> BcStates0P{6.716, 7.121};       // GeV
//        addSetting<vector<double>>("BcStates1P", BcStates1P);   //f, F1, T2, T23
//        addSetting<vector<double>>("BcStates1M", BcStates1M);   //g, fp, fT, T1
//        addSetting<vector<double>>("BcStates0M", BcStates0M);   //P1|F2
//        addSetting<vector<double>>("BcStates0P", BcStates0P);   //f0
//    
//    }
    
    void FFBGLXBase::setTZ(const double& t,  const double& tp, const double& tm, const double& t0) {
        _tp = tp;
        _tm = tm;
        _t0 = t0;
        _t = t;
        
        const double sqrttpt = sqrt(tp - t);
        const double sqrttptm = sqrt(tp - tm);
        const double sqrttpt0 = sqrt(tp - t0);
        _z = (sqrttpt - sqrttpt0)/(sqrttpt + sqrttpt0);
        _zm = (sqrttptm - sqrttpt0)/(sqrttptm + sqrttpt0);
    }
    
    void FFBGLXBase::fillOuterFactors(Tq2Val tq2) {
        const double t = (tq2 == TQ2 ? _t : _tm);
        const double sqrttpt = sqrt(_tp - t);
        const double sqrttpt0 = sqrt(_tp - _t0);
        const double sqrttptm = sqrt(_tp - _tm);
        const double sqrttp = sqrt(_tp);
        
        _mPhiA1 = pow(sqrttpt, 0.5);
        _mPhiA3 = pow(sqrttpt, 1.5);
        _mPhiB1 = pow(sqrttpt + sqrttptm, 0.5);
        _mPhiB3 = pow(sqrttpt + sqrttptm, 1.5);
        _mPhiC1 = pow(sqrttpt + sqrttp, -(1. + 3.));
        _mPhiC2 = pow(sqrttpt + sqrttp, -(2. + 3.));
        
        _mPhiPre = pow(sqrttpt/sqrttpt0, 0.5) * (sqrttpt + sqrttpt0);
    }
    
    void FFBGLXBase::fillOutersJ1P() {  
        vector<Tq2Val> tq2s{TQ2, TQ2M}; 
        for(auto tq2 : tq2s){
            fillOuterFactors(tq2);
        
            //Defined up to sqrt(nI/chi) prefactor
            _mPhi[F][tq2] = 1./sqrt(24. * pi) * _mPhiPre * _mPhiA1 * _mPhiB1 * _mPhiC1;
            _mPhi[F1][tq2] = 1./sqrt(48. * pi) * _mPhiPre * _mPhiA1 * _mPhiB1 * _mPhiC2;
            _mPhi[T2][tq2] = 1./sqrt(24./(_tp*_tm) * pi) * _mPhiPre * _mPhiA1 * _mPhiB1 * _mPhiC2;
            _mPhi[T23][tq2] = 1./sqrt(48.*_tp/((_tp-_tm)*(_tp-_tm)) * pi) * _mPhiPre * _mPhiA1 * _mPhiB1 * _mPhiC1; 
        }
    
    }
    
    void FFBGLXBase::fillOutersJ1M() {
        vector<Tq2Val> tq2s{TQ2, TQ2M}; 
        for(auto tq2 : tq2s){
            fillOuterFactors(tq2);
        
            //Defined up to sqrt(nI/chi) prefactor
            _mPhi[G][tq2] = 1./sqrt(96. * pi) * _mPhiPre * _mPhiA3 * _mPhiB3 * _mPhiC1;
            _mPhi[FP][tq2] = 1./sqrt(48. * pi) * _mPhiPre * _mPhiA3 * _mPhiB3 * _mPhiC2;
            _mPhi[T1][tq2] = 1./sqrt(24. * pi) * _mPhiPre * _mPhiA3 * _mPhiB3 * _mPhiC2;
            _mPhi[FT][tq2] = 1./sqrt(48. * _tp * pi) * _mPhiPre * _mPhiA3 * _mPhiB3 * _mPhiC1;
        }
    }
    
    void FFBGLXBase::fillOutersJ0P() {
        fillOuterFactors(TQ2);
    
        //Defined up to sqrt(nI/chi) prefactor
        _mPhi[F0][TQ2] = 1./sqrt(16. * pi) * _mPhiPre * _mPhiA1 * _mPhiB1 * _mPhiC1;
    }
    
    void FFBGLXBase::fillOutersJ0M() {
        fillOuterFactors(TQ2);
        
        //Defined up to sqrt(nI/chi) prefactor
        _mPhi[F2][TQ2] = 1./sqrt(64. * pi) * _mPhiPre * _mPhiA3 * _mPhiB3 * _mPhiC1;
    }
    
    void FFBGLXBase::fillBlaschkes(const vector<double>& BcStates, const PjName J) {
        const double sqrttpt0 = sqrt(_tp - _t0);
        
        _mPJ[J] = 1.;
        for (double mBc : BcStates) {
            const double mBc2 = mBc*mBc*_units*_units;
            if (mBc2 > _tp) {
                MSG_WARNING("Subthreshold pole at " + to_string(sqrt(mBc2)) + " exceeds branch point at "  + to_string(sqrt(_tp)) + ". Pole discarded!");
                continue;
            }
            const double sqtpbc = sqrt(_tp - mBc*mBc*_units*_units);
            const double zBc = (sqtpbc - sqrttpt0)/(sqtpbc + sqrttpt0);
            _mPJ[J] *= ((_z - zBc) / (1 - _z * zBc));
        } 
    }
    
} // namespace Hammer
