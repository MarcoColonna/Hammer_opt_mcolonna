///
/// @file  FFBtoD2starISGW2.cc
/// @brief \f$ B \rightarrow D_2^* \f$ ISGW2 form factors
/// @brief Ported directly from EvtGen
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/FormFactors/ISGW2/FFBtoD2starISGW2.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFBtoD2starISGW2::FFBtoD2starISGW2() {
        // Create tensor rank and dimensions
        IndexList dims = {8};
        string name{"FFBtoD2starISGW2"};

        setPrefix("BtoD**2*");
        addProcessSignature(PID::BPLUS, {-PID::DSSD2STAR});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BDSSD2STAR})});

        addProcessSignature(PID::BZERO, {PID::DSSD2STARMINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BDSSD2STAR})});

        setPrefix("BstoDs**2*");
        addProcessSignature(PID::BS, {PID::DSSDS2STARMINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BSDSSDS2STAR})});

        setSignatureIndex();
    }

    void FFBtoD2starISGW2::defineSettings() {
        //_mFFErrNames = ;
        setPath(getFFErrPrefixGroup().get());
        setUnits("GeV");

        setInitialized();
    }

    void FFBtoD2starISGW2::evalAtPSPoint(const vector<double>& point, const vector<double>& masses) {
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
        // Ported from EvtGen EvtISGW2FF3P2
        switch (result.labels()[0]) {
        case FF_BDSSD2STAR:
            msb = 5.2;
            msd = 0.33;
            bb2 = 0.431 * 0.431;
            mbb = 5.31;

            msq = 1.82;
            bx2 = 0.33 * 0.33;
            mbx = (5.0 * 2.46 + 3.0 * 2.42) / 8.0;
            nfp = 3.0;
            break;
        case FF_BSDSSDS2STAR:
            msb = 5.2;
            msd = 0.55;
            bb2 = 0.54 * 0.54;
            mbb = 5.38;

            msq = 1.82;
            bx2 = 0.41 * 0.41;
            mbx = (5.0 * 2.61 + 3.0 * 2.54) / 8.0;
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

        const double mup = 1.0 / (1.0 / msq + 1.0 / msb);
        const double mum = 1.0 / (1.0 / msq - 1.0 / msb);
        const double bbx2 = 0.5 * (bb2 + bx2);
        const double tm = (mb - mx) * (mb - mx);
        if (t > tm) {
            t = 0.99 * tm;
        }
        const double wt = 1.0 + (tm - t) / (2.0 * mbb * mbx);

        const double mqm = 0.1;
        const double r2 = 3.0 / (4.0 * msb * msq) + 3 * msd * msd / (2 * mbb * mbx * bbx2) +
                          (16.0 / (mbb * mbx * (33.0 - 2.0 * nfp))) * log(Getas(mqm) / Getas(msq));

        const double f5 =
            sqrt(mtx / mtb) * pow(sqrt(bx2 * bb2) / bbx2, 5.0 / 2.0) / (pow((1.0 + r2 * (tm - t) / 18.0), 3.0));

        const double f5h = f5 * pow((mbb / mtb), -1.5) * pow((mbx / mtx), -0.5);
        const double f5k = f5 * pow((mbb / mtb), -0.5) * pow((mbx / mtx), 0.5);
        const double f5bppbm = f5 * pow((mbb / mtb), -2.5) * pow((mbx / mtx), 0.5);
        const double f5bpmbm = f5 * pow((mbb / mtb), -1.5) * pow((mbx / mtx), -0.5);

        const double hf =
            f5h * (msd / (sqrt(8.0 * bb2) * mtb)) * ((1.0 / msq) - (msd * bb2 / (2.0 * mum * mtx * bbx2)));
        const double kf = f5k * (msd / (sqrt(2.0 * bb2))) * (1.0 + wt);
        const double bppbm = ((msd * msd * f5bppbm * bx2) / (sqrt(32.0 * bb2) * msq * msb * mtb * bbx2)) *
                             (1.0 - (msd * bx2 / (2.0 * mtb * bbx2)));
        const double bpmbm = -1.0 * (msd * f5bpmbm / (sqrt(2.0 * bb2) * msb * mtx)) *
                             (1.0 - ((msd * msb * bx2) / (2.0 * mup * mtb * bbx2)) +
                              ((msd * bx2 * (1.0 - ((msd * bx2) / (2.0 * mtb * bbx2)))) / (4.0 * msq * bbx2)));

        // FF Basis map
        const double sqmbmx = sqrt(mb * mx);
        const double Ka1 = kf * mb / sqmbmx;
        const double Ka2 = mb * mb * mb * bppbm / sqmbmx;
        const double Ka3 = bpmbm * mb * sqmbmx;
        const double Kv = 2 * hf * mb * sqmbmx;

        // Set elements
        result.element({1}) = Ka1;
        result.element({2}) = Ka2;
        result.element({3}) = Ka3;
        result.element({4}) = Kv;
    }

    std::unique_ptr<FormFactorBase> FFBtoD2starISGW2::clone(const std::string& label) {
        MAKE_CLONE(FFBtoD2starISGW2, label);
    }

} // namespace Hammer
