///
/// @file  FFLbtoLcstar32PCR.cc
/// @brief \f$ \Lambda_b \rightarrow Lambda_c^*(2625) \f$ PCR form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/FormFactors/PCR/FFLbtoLcstar32PCR.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFLbtoLcstar32PCR::FFLbtoLcstar32PCR() {
        // Create tensor rank and dimensions
        IndexList dims = {16};
        string name{"FFLbtoLcstar32PCR"};

        setPrefix("LbtoLc*2625");
        addProcessSignature(-PID::LAMBDAB, {PID::LAMBDACSTAR32MINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_LBLCSTAR32})});

        setSignatureIndex();
    }

    void FFLbtoLcstar32PCR::defineSettings() {
        //_mFFErrNames = ;
        setPath(getFFErrPrefixGroup().get());
        setUnits("GeV");

        addSetting<double>("mQ", 5.28);
        addSetting<double>("mq", 1.89);
        addSetting<double>("md", 0.40);
        addSetting<double>("aL", 0.59);
        addSetting<double>("aLp", 0.47);
        addSetting<double>("k2", 1.0);

        setInitialized();
    }

    void FFLbtoLcstar32PCR::eval(const Particle& parent, const ParticleList& daughters,
                                 const ParticleList& /*references*/) {
        // Momenta
        const FourMomentum& pLbmes = parent.momentum();
        const FourMomentum& pLcmes = daughters[0].momentum();
        // const FourMomentum& pTau = daughters[2].momentum();

        // kinematic objects
        const double Mb = pLbmes.mass();
        const double Mc = pLcmes.mass();
        const double Sqq = Mb * Mb + Mc * Mc - 2. * (pLbmes * pLcmes);

        evalAtPSPoint({Sqq}, {Mb, Mc});
    }

    void FFLbtoLcstar32PCR::evalAtPSPoint(const vector<double>& point, const vector<double>& masses) {
        Tensor& result = getTensor();
        result.clearData();

        if (!isInitialized()) {
            MSG_WARNING("Warning, Settings have not been defined!");
        }

        double Mb = 0.;
        double Mc = 0.;
        double unitres = 1.;
        tie(Mb, Mc, unitres) = getParentDaughterHadMasses(masses);


        double Sqq = point[0];

        // const double sqSqq = sqrt(Sqq);
        double w = getW(Sqq, Mb, Mc);
        // safety measure if w==1.0
        if (isZero(w - 1.0)) {
            w += 1e-6;
        }

        // Parameters
        const double mQ = (*getSetting<double>("mQ"));
        const double mq = (*getSetting<double>("mq"));
        const double md = (*getSetting<double>("md"));
        const double aL = (*getSetting<double>("aL"));
        const double aLp = (*getSetting<double>("aLp"));
        const double k2 = (*getSetting<double>("k2"));

        const double aL2 = aL * aL;
        const double aLp2 = aLp * aLp;
        const double aLLp2 = 0.5 * (aL2 + aLp2);
        const double rho2 = 3. * md * md / (2. * k2 * aLLp2);

        const double IH = -1. / sqrt3 * pow(aL * aLp / aLLp2, 2.5) * exp(-rho2 * (w * w - 1.));

        // set elements F_1234 and G_1234
        result.element({2}) = 3. * md / aL * (1. + (md / aLLp2) * ((aLp2 / mq) + (aL2 / mQ)));
        result.element({3}) =
            -((3. * md * md / mq) * (aLp2 / (aLLp2 * aL2)) - 5. * aL * aLp2 * md / (4. * aLLp2 * mq * mQ));
        result.element({4}) = -(3. * md * md * aL / (mQ * aLLp2) + aL / (2. * mQ));
        result.element({5}) = aL / mQ;
        result.element({6}) = 3. * md / aL - (aL / (2. * mQ)) * (1. + 3. * md * aLp2 / (2. * aLLp2 * mq));
        result.element({7}) = -((3. * md * md / mq) * (aLp2 / (aLLp2 * aL)) +
                                aL * aLp2 * md / (4. * aLLp2 * aLLp2 * mq * mQ) * (aLLp2 + 12. * md * md));
        result.element({8}) =
            aL / (mQ * aLLp2) * (aLLp2 / 2. + 3. * md * md + aLp2 * md / (mq * aLLp2) * (aLLp2 + 6. * md * md));
        result.element({9}) = -(aL / mQ + md / (mq * mQ) * aLp2 * aL / aLLp2);

        result *= IH;

        //        for(size_t idx = 0; idx < 8; ++idx){
        //            const IndexType idxp = static_cast<IndexType>(idx + 2u);
        //            cout << result.element({idxp}).real() << " ";
        //        }
        //        cout << '\n';
    }

    std::unique_ptr<FormFactorBase> FFLbtoLcstar32PCR::clone(const std::string& label) {
        MAKE_CLONE(FFLbtoLcstar32PCR, label);
    }

} // namespace Hammer
