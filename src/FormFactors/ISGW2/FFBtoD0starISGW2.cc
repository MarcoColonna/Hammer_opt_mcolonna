///
/// @file  FFBtoD0starISGW2.cc
/// @brief \f$ B \rightarrow D_0^* \f$ ISGW2 form factors
/// @brief Ported directly from EvtGen
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/FormFactors/ISGW2/FFBtoD0starISGW2.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFBtoD0starISGW2::FFBtoD0starISGW2() {
        // Create tensor rank and dimensions
        IndexList dims = {4};
        string name{"FFBtoD0starISGW2"};

        setPrefix("BtoD**0*");
        addProcessSignature(PID::BPLUS, {-PID::DSSD0STAR});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BDSSD0STAR})});

        addProcessSignature(PID::BZERO, {PID::DSSD0STARMINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BDSSD0STAR})});

        setPrefix("BstoDs**0*");
        addProcessSignature(PID::BS, {PID::DSSDS0STARMINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BSDSSDS0STAR})});

        setSignatureIndex();
    }

    void FFBtoD0starISGW2::defineSettings() {
        //_mFFErrNames = ;
        setPath(getFFErrPrefixGroup().get());
        setUnits("GeV");

        setInitialized();
    }

    void FFBtoD0starISGW2::evalAtPSPoint(const vector<double>& point, const vector<double>& masses) {
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
        double t = Sqq / (unitres * unitres);

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wswitch-enum"
        // Ported from EvtGen EvtISGW2FF3P0
        switch (result.labels()[0]) {
        case FF_BDSSD0STAR:
            msb = 5.2;
            msd = 0.33;
            bb2 = 0.431 * 0.431;
            mbb = 5.31;

            msq = 1.82;
            bx2 = 0.33 * 0.33;
            mbx = (3.0 * 2.49 + 2.40) / 4.0;
            nfp = 3.0;
            break;
        case FF_BSDSSDS0STAR:
            msb = 5.2;
            msd = 0.55;
            bb2 = 0.54 * 0.54;
            mbb = 5.38;

            msq = 1.82;
            bx2 = 0.41 * 0.41;
            mbx = (3.0 * 2.54 + 2.46) / 4.0;
            nfp = 3.0;
            break;
        default:
            MSG_ERROR("Unknown assignments for parametrization " + getFFErrPrefixGroup().get() + ".");
        }
#pragma clang diagnostic pop

        const double mtb = msb + msd;
        const double mtx = msq + msd;

        const double mb = Mb / unitres;
        const double mx = Mc / unitres;

        const double bbx2 = 0.5 * (bb2 + bx2);
        const double tm = (mb - mx) * (mb - mx);
        if (t > tm) {
            t = 0.99 * tm;
        }

        const double mqm = 0.1;
        const double r2 = 3.0 / (4.0 * msb * msq) + 3 * msd * msd / (2 * mbb * mbx * bbx2) +
                          (16.0 / (mbb * mbx * (33.0 - 2.0 * nfp))) * log(Getas(mqm) / Getas(msq));

        const double f5 =
            sqrt(mtx / mtb) * pow(sqrt(bx2 * bb2) / bbx2, 5.0 / 2.0) / (pow((1.0 + r2 * (tm - t) / 18.0), 3.0));

        const double f5uppum = f5 * pow((mbb / mtb), -0.5) * pow((mbx / mtx), 0.5);
        const double f5upmum = f5 * pow((mbb / mtb), 0.5) * pow((mbx / mtx), -0.5);

        const double uppum = -1.0 * f5uppum * sqrt(2.0 / (3.0 * bb2)) * msd;
        const double upmum = 1.0 * f5upmum * sqrt(2.0 / (3.0 * bb2)) * msd * mtb / mtx;

        // FF Basis map
        const double gppgm = mb / sqrt(mb * mx) * uppum;
        const double gpmgm = mx / sqrt(mb * mx) * upmum;

        const double gp = (gppgm + gpmgm) / 2.0;
        const double gm = (gppgm - gpmgm) / 2.0;

        // Set elements
        result.element({1}) = gp;
        result.element({2}) = gm;
    }

    std::unique_ptr<FormFactorBase> FFBtoD0starISGW2::clone(const std::string& label) {
        MAKE_CLONE(FFBtoD0starISGW2, label);
    }

} // namespace Hammer
