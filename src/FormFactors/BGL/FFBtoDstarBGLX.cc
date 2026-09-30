///
/// @file  FFBtoDstarBGLX.cc
/// @brief \f$ B \rightarrow D \f$ BGL Extended form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/FormFactors/BGL/FFBtoDstarBGLX.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFBtoDstarBGLX::FFBtoDstarBGLX() {
        // Create tensor rank and dimensions
        IndexList dims = {8};
        string name{"FFBtoDstarBGLX"};

        setPrefix("BtoD*");
        addProcessSignature(PID::BPLUS, {-PID::DSTAR});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BDSTAR})});

        addProcessSignature(PID::BZERO, {PID::DSTARMINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BDSTAR})});

        setPrefix("BstoDs*");
        addProcessSignature(PID::BS, {PID::DSSTARMINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BSDSSTAR})});

        setSignatureIndex();
    }

    void FFBtoDstarBGLX::defineSettings() {
        //_mFFErrNames = ;
        setPath(getFFErrPrefixGroup().get());
        setUnits("GeV");

        // Optimized z-expansion
        addSetting<bool>("OptZ", false);
        //Use externally specified branch points rather than even kinematics
        addSetting<bool>("ExtBranches", false);
        
        // values from hep-ph/9705252, u = 0.33
        addSetting<double>("Vcb", 41.5e-3);        
        addSetting<double>("ChimT", 3.068e-4); // GeV^-2  f, F1
        addSetting<double>("ChipT", 5.280e-4); // GeV^-2  g
        addSetting<double>("ChimL", 2.466e-3); // P1|F2
        
        addSetting<double>("Chi1T", 4.98e-4 * 2.6/2.0); // GeV^-2  from 2507.03569 NB different matching scale and sym factor
        addSetting<double>("Chi1AT", 2.77e-4 * 2.6/2.0); // GeV^-2 from 2507.03569 NB different matching scale and sym factor
        
        // pole mass from hep-ph/9705252
        addSetting<double>("mb", 4.9);       // GeV
        addSetting<double>("mc", 4.9 - 3.4); // GeV (using mb-mc = 3.4 GeV)

        // Using 222 fit from 1902.09553
        vector<double> agvec = {0.00038, 0.026905, 0.};
        vector<double> afvec = {0.00055, -0.0020370, 0.};
        vector<double> aF1vec = {-0.000433, 0.005353};
        // Approximate values from 1707.09509 (adapted by rescaling wrt ChimL value)
        vector<double> aP1vec = {0.007, -0.036};
        // From a fit to synthetic data of 2304.03137. (NB aP1vec data is consistent with this)
        // T1(q^2 = 0) = T2(q^2 = 0) is not (yet) directly enforced, but obeyed by these settings
        vector<double> aT1vec = {0.0004695, -0.002800, -0.04079};
        vector<double> aT2vec = {0.0001230, 0.0006294, -0.03221};
        vector<double> aT23vec = {0.0003665, 0.001322, -0.033100};
        
        addSetting<vector<double>>("agvec", agvec);
        addSetting<vector<double>>("afvec", afvec);
        addSetting<vector<double>>("aF1vec", aF1vec);
        addSetting<vector<double>>("aP1vec", aP1vec);
        addSetting<vector<double>>("aT1vec", aT1vec);
        addSetting<vector<double>>("aT2vec", aT2vec);
        addSetting<vector<double>>("aT23vec", aT23vec);

        vector<double> BcStatesJ1P{6.730, 6.736, 7.135, 7.142}; // GeV
        vector<double> BcStatesJ1M{6.337, 6.899, 7.012, 7.280}; // GeV
        vector<double> BcStatesJ0M{6.275, 6.842, 7.250};       // GeV
        addSetting<vector<double>>("BcStatesJ1P", BcStatesJ1P);   // f, F1, T2, T23
        addSetting<vector<double>>("BcStatesJ1M", BcStatesJ1M);   // g, T1
        addSetting<vector<double>>("BcStatesJ0M", BcStatesJ0M);   // P1|F2

        //Branch points for each current from lowest mass hadron pair
        addSetting<double>("ExtBranchJ1P", 5.324 + 1.867); // GeV mB* + mD
        addSetting<double>("ExtBranchJ1M", 5.280 + 1.867); // GeV mB + mD
        addSetting<double>("ExtBranchJ0M", 5.324 + 1.867); // GeV mB* + mD
        
        addSetting<double>("nc", 2.6); // Effective SU(3) symmetry factor
        
        setInitialized();
    }

    void FFBtoDstarBGLX::evalAtPSPoint(const vector<double>& point, const vector<double>& masses) {
        Tensor& result = getTensor();
        result.clearData();

        if (!isInitialized()) {
            MSG_WARNING("Warning, Settings have not been defined!");
        }

        double Mb = 0.;
        double Mc = 0.;
        double unitres = 1.;
        tie(Mb, Mc, unitres) = getParentDaughterHadMasses(masses);
        const double Sqq = point[0]; // q2
        const double Mb2 = Mb * Mb;
        const double Mc2 = Mc * Mc;
        const double Mb3 = Mb2 * Mb;
        const double rC = Mc / Mb;
        const double rC2 = rC * rC;
        const double sqrC = sqrt(rC);

        double w = getW(Sqq, Mb, Mc);
        // safety measure if w==1.0
        if (isZero(w - 1.0)) {
            w += 1e-6;
        }
        const double w2 = w * w;


        const vector<double>& ag = (*getSetting<vector<double>>("agvec"));
        const vector<double>& af = (*getSetting<vector<double>>("afvec"));
        const vector<double>& aF1 = (*getSetting<vector<double>>("aF1vec"));
        const vector<double>& aP1 = (*getSetting<vector<double>>("aP1vec")); //F2
        const vector<double>& aT1 = (*getSetting<vector<double>>("aT1vec"));
        const vector<double>& aT2 = (*getSetting<vector<double>>("aT2vec"));
        const vector<double>& aT23 = (*getSetting<vector<double>>("aT23vec"));
        
        const double nc = (*getSetting<double>("nc"));
        const double etaEWVcb = 1.0066 * (*getSetting<double>("Vcb"));
        const double chimT = (*getSetting<double>("ChimT")) / (unitres * unitres);
        const double chipT = (*getSetting<double>("ChipT")) / (unitres * unitres);
        const double chimL = (*getSetting<double>("ChimL"));
        const double chi1T = (*getSetting<double>("Chi1T")) / (unitres * unitres);
        const double chi1AT = (*getSetting<double>("Chi1AT")) / (unitres * unitres);
        
        const bool OptZ = (*getSetting<bool>("OptZ"));
        const bool useextbranch = (*getSetting<bool>("ExtBranches")); 
        //Default: Branch point coincides with event mass pair threshold
        double tp = (Mb + Mc) * (Mb + Mc);
        const double tm = (Mb - Mc) * (Mb - Mc);
        double t0 = tm; // Choose z(q2max) = 0
        if(OptZ) { t0 = tp * (1. - sqrt(1. - tm/tp)); } 

        setTZ(Sqq, tp, tm, t0);  //Set internal tp, tm, t0 and z, z(q2max) values. Can be called for each current separately 
        
        vector<double>& BcStatesJ1P = (*getSetting<vector<double>>("BcStatesJ1P")); //should remain in GeV per fillBlaschkes
        vector<double>& BcStatesJ1M = (*getSetting<vector<double>>("BcStatesJ1M"));
        vector<double>& BcStatesJ0M = (*getSetting<vector<double>>("BcStatesJ0M"));

        const double mb = (*getSetting<double>("mb")) * unitres;
        const double mc = (*getSetting<double>("mc")) * unitres;
        
        // *** J = 1P current FFs: f, F1, T2, T23
        if(useextbranch){
            tp = pow((*getSetting<double>("ExtBranchJ1P")) * unitres, 2);
            if(OptZ) { t0 = tp * (1. - sqrt(1. - tm/tp)); } 
            setTZ(Sqq, tp, tm, t0);
        }
        fillOutersJ1P();
        fillBlaschkes(BcStatesJ1P, J1P);
        
        double f = 0;
        double fmax = 0;
        double zpow = 1.;
        double zmpow = 1.;
        for (size_t n = 0; n < af.size(); ++n) {
            f += af[n] * zpow;
            fmax += af[n] * zmpow;
            zpow *= _z;
            zmpow *= _zm;
        }
        f /= (_mPJ[J1P] * _mPhi[F].at(TQ2) * sqrt(nc/chimT));
        
        // Leading F1 coefficient c_0 set by constraint F1(q^2=q^2_max) = (Mb - Mc) f(q^2=q^2_max), i.e at t = tm.
        // Implemented to permit general choice of t0.
        double f1 = fmax * (Mb - Mc) * _mPhi[F1].at(TQ2M) / _mPhi[F].at(TQ2M); 
        zpow = _z;
        zmpow = _zm;
        for (size_t n = 0; n < aF1.size(); ++n) {
            f1 += aF1[n] * (zpow - zmpow);
            zpow *= _z;
            zmpow *= _zm;
        }
        f1 /= (_mPJ[J1P] * _mPhi[F1].at(TQ2) * sqrt(nc/chimT));
        
        double t2 = 0;
        zpow = 1.;
        for (size_t n = 0; n < aT2.size(); ++n) {
            t2 += aT2[n] * zpow;
            zpow *= _z;
        }
        t2 /= (_mPJ[J1P] * _mPhi[T2].at(TQ2) * sqrt(nc/chi1AT));
        
        double t23 = 0;
        zpow = 1.;
        for (size_t n = 0; n < aT23.size(); ++n) {
            t23 += aT23[n] * zpow;
            zpow *= _z;
        }
        t23 /= (_mPJ[J1P] * _mPhi[T23].at(TQ2) * sqrt(nc/chi1AT));
        
        // *** J = 1M current FFs: g, T1
        if(useextbranch){
            tp = pow((*getSetting<double>("ExtBranchJ1M")) * unitres, 2);
            if(OptZ) { t0 = tp * (1. - sqrt(1. - tm/tp)); } 
            setTZ(Sqq, tp, tm, t0);
        }
        fillOutersJ1M();
        fillBlaschkes(BcStatesJ1M, J1M);
        
        double g = 0;
        zpow = 1.;
        for (size_t n = 0; n < ag.size(); ++n) {
            g += ag[n] * zpow;
            zpow *= _z;
        }
        g /= (_mPJ[J1M] * _mPhi[G].at(TQ2) * sqrt(nc/chipT));
        
        double t1 = 0;
        zpow = 1.;
        for (size_t n = 0; n < aT1.size(); ++n) {
            t1 += aT1[n] * zpow;
            zpow *= _z;
        }
        t1 /= (_mPJ[J1M] * _mPhi[T1].at(TQ2) * sqrt(nc/chi1T));
        
        //J = 0M current FFs: P1|F2 (conversion from F2 to P1 by scaling factor sqrC/(1+rC). See 1707.09509)
        if(useextbranch){
            tp = pow((*getSetting<double>("ExtBranchJ0M")) * unitres, 2);
            if(OptZ) { t0 = tp * (1. - sqrt(1. - tm/tp)); } 
            setTZ(Sqq, tp, tm, t0);
        }
        fillOutersJ0M();
        fillBlaschkes(BcStatesJ0M, J0M);
        
        double f2 = 0;
        zpow = 1.;
        for (size_t n = 0; n < aP1.size(); ++n) {
            f2 += aP1[n] * zpow;
            zpow *= _z;
        }
        f2 /= (_mPJ[J0M] * _mPhi[F2].at(TQ2) * sqrt(nc/chimL));
        double p1 = sqrC / (1. + rC) * f2;
         
        // Mapping to amplitude FF basis
        const double Fpf = (rC - w) / (2. * rC * Mb2 * (w2 - 1));
        const double FpF1 = 1. / (2. * rC * Mb3 * (w2 - 1));
        const double Fmf = (rC + w) / (2. * rC * Mb2 * (w2 - 1));
        const double FmF1 = 1. / (2. * rC * Mb3 * (w2 - 1)) * (rC2 - 1) / (1 + rC2 - 2 * rC * w);
        const double FmP1 = sqrC * (rC + 1) / (Mb * (1 + rC2 - 2 * rC * w));

        const double fm = (Fmf * f + FmF1 * f1 + FmP1 * p1);
        const double fp = (Fpf * f + FpF1 * f1);

        // Pseudoscalar from eqn of motion
        // const double fps = -(f + fp * (Mb2 - Mc2) + fm * Sqq) / (mb + mc);
        const double fps = -sqrC * (1. + rC) * Mb * p1 / (mb + mc);
        // Transformation to t3
        const double Ft2 = (Mb2 - Mc2)*(Mb2 + 3. * Mc2 - Sqq)/ ( (Mb2 - Mc2 - Sqq)*(Mb2 - Mc2 - Sqq) - 4*Mc2*Sqq );
        const double Ft23 = 8.*Mb*Mc2*(Mb - Mc)/( (Mb2 - Mc2 - Sqq)*(Mb2 - Mc2 - Sqq) - 4*Mc2*Sqq );
        
        const double t3 = Ft2*t2 - Ft23*t23;
        
        // set elements, mapping to amplitude FF basis
        // Fps (dim 0)
        result.element({0}) = fps;
        // Ff (dim 1)
        result.element({1}) = f;
        // Fg (dim -1) Definition in BGL (hep-ph/9705252) differs from Manohar, Wise or 1610.02045 by factor of 2
        result.element({2}) = g / 2.;
        // Fm (dim -1)
        result.element({3}) = fm;
        // Fp (dim -1)
        result.element({4}) = fp;
        // Fzt (dim -2)
        result.element({5}) = -(t1 - t2)/Sqq + t3/(Mb2 - Mc2);
        // Fmt (dim 0)
        result.element({6}) = (Mb2 - Mc2)/Sqq * (t1 - t2);
        // Fpt (dim 0)
        result.element({7}) = -t1;
        
        result *= (1. / etaEWVcb);
        
    }

    unique_ptr<FormFactorBase> FFBtoDstarBGLX::clone(const string& label) {
        MAKE_CLONE(FFBtoDstarBGLX, label);
    }

} // namespace Hammer
