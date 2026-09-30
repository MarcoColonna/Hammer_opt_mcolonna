///
/// @file  FFBtoDISGW2.cc
/// @brief \f$ B \rightarrow D \f$ ISGW2 form factors
/// @brief Ported directly from EvtGen
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/FormFactors/ISGW2/FFBtoDISGW2.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFBtoDISGW2::FFBtoDISGW2() {
        // Create tensor rank and dimensions
        IndexList dims = {4};
        string name{"FFBtoDISGW2"};

        setPrefix("BtoD");
        addProcessSignature(PID::BPLUS, {-PID::D0});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BD})});

        addProcessSignature(PID::BZERO, {PID::DMINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BD})});

        setPrefix("BstoDs");
        addProcessSignature(PID::BS, {PID::DSMINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BSDS})});

        setSignatureIndex();
    }

    void FFBtoDISGW2::defineSettings() {
        //_mFFErrNames = ;
        setPath(getFFErrPrefixGroup().get());
        setUnits("GeV");

        setInitialized();
    }

    void FFBtoDISGW2::evalAtPSPoint(const vector<double>& point, const vector<double>& masses) {
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
        // Ported from EvtGen EvtISGW2FF1S0
        switch (result.labels()[0]) {
        case FF_BD:
            msb = 5.2;
            msd = 0.33;
            bb2 = 0.431 * 0.431;
            mbb = 5.31;
            nf = 4.0;

            msq = 1.82;
            bx2 = 0.45 * 0.45;
            mbx = 0.75 * 2.01 + 0.25 * 1.87;
            nfp = 3.0;
            break;
        case FF_BSDS:
            msb = 5.2;
            msd = 0.55;
            bb2 = 0.54 * 0.54;
            mbb = 5.38;
            nf = 4.0;

            msq = 1.82;
            bx2 = 0.56 * 0.56;
            mbx = 0.75 * 2.11 + 0.25 * 1.97;
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
        const double bbx2 = 0.5 * (bb2 + bx2);
        const double tm = (mb - mx) * (mb - mx);
        if (t > tm) {
            t = 0.99 * tm;
        }

        const double mqm = 0.1;
        const double r2 = 3.0 / (4.0 * msb * msq) + 3 * msd * msd / (2 * mbb * mbx * bbx2) +
                          (16.0 / (mbb * mbx * (33.0 - 2.0 * nfp))) * log(Getas(mqm, mqm) / Getas(msq, msq));

        const double f3 = sqrt(mtx / mtb) * pow(sqrt(bx2 * bb2) / bbx2, 1.5) /
                          ((1.0 + r2 * (tm - t) / 12.0) * (1.0 + r2 * (tm - t) / 12.0));

        const double ai = -1.0 * (6.0 / (33.0 - 2.0 * nf));
        const double cji = pow((Getas(msb, msb) / Getas(msq, msq)), ai);

        const double zji = msq / msb;

        const double gammaji = GetGammaji(zji);
        const double chiji = -1.0 - (gammaji / (1 - zji));
        const double betaji_fppfm = gammaji - (2.0 / 3.0) * chiji;
        const double betaji_fpmfm = gammaji + (2.0 / 3.0) * chiji;
        const double rfppfm = cji * (1.0 + betaji_fppfm * Getas(msq, sqrt(msb * msq)) / pi);
        const double rfpmfm = cji * (1.0 + betaji_fpmfm * Getas(msq, sqrt(msb * msq)) / pi);
        const double f3fppfm = f3 * pow((mbb / mtb), -0.5) * pow((mbx / mtx), 0.5);
        const double f3fpmfm = f3 * pow((mbb / mtb), 0.5) * pow((mbx / mtx), -0.5);
        const double fppfm =
            f3fppfm * rfppfm * (2.0 - ((mtx / msq) * (1 - ((msd * msq * bb2) / (2.0 * mup * mtx * bbx2)))));
        const double fpmfm = f3fpmfm * rfpmfm * (mtb / msq) * (1 - ((msd * msq * bb2) / (2.0 * mup * mtx * bbx2)));

        const double fp = (fppfm + fpmfm) / 2.0;
        const double fm = (fppfm - fpmfm) / 2.0;

        const double f0 = (fm / ((mb * mb - mx * mx) / t)) + fp;

        // Temporary scalar NP form factors
        const double MbcSqq = pow(mb + mx, 2.) - t;
        const double fs = fp * MbcSqq / (mb + mx);
        //        double ft = fp / (Mb + Mc);

        // set elements, mapping to amplitude FF basis
        // Fs
        result.element({0}) = unitres * fs;
        // Fz
        result.element({1}) = f0;
        // Fp
        result.element({2}) = fp;
        // Ft
        //  result.element({3}) = 0;
    }

    std::unique_ptr<FormFactorBase> FFBtoDISGW2::clone(const std::string& label) {
        MAKE_CLONE(FFBtoDISGW2, label);
    }

} // namespace Hammer
