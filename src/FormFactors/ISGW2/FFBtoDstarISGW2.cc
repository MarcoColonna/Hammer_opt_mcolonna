///
/// @file  FFBtoDstarISGW2.cc
/// @brief \f$ B \rightarrow D^* \f$ ISGW2 form factors
/// @brief Ported directly from EvtGen
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/FormFactors/ISGW2/FFBtoDstarISGW2.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFBtoDstarISGW2::FFBtoDstarISGW2() {
        // Create tensor rank and dimensions
        IndexList dims = {8};
        string name{"FFBtoDstarISGW2"};

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

    void FFBtoDstarISGW2::defineSettings() {
        //_mFFErrNames = ;
        setPath(getFFErrPrefixGroup().get());
        setUnits("GeV");

        setInitialized();
    }

    void FFBtoDstarISGW2::evalAtPSPoint(const vector<double>& point, const vector<double>& masses) {
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
        // Ported from EvtGen EvtISGW2FF3S1
        switch (result.labels()[0]) {
        case FF_BDSTAR:
            msb = 5.2;
            msd = 0.33;
            bb2 = 0.431 * 0.431;
            mbb = 5.31;
            nf = 4.0;
            cf = 0.989;
            msq = 1.82;
            bx2 = 0.38 * 0.38;
            mbx = 0.75 * 2.01 + 0.25 * 1.87;
            nfp = 3.0;
            break;
        case FF_BSDSSTAR:
            msb = 5.2;
            msd = 0.55;
            bb2 = 0.54 * 0.54;
            mbb = 5.38;
            nf = 4.0;
            cf = 0.984;
            msq = 1.82;
            bx2 = 0.49 * 0.49;
            mbx = 0.75 * 2.11 + 0.25 * 1.97;
            nfp = 3.0;
            break;
        default:
            MSG_ERROR("Unknown assignments for parametrization " + getFFErrPrefixGroup().get() + ".");
        }
#pragma clang diagnostic pop

        const double mtb = msb + msd;
        const double mtx = msq + msd;

        const double mup = 1.0 / (1.0 / msq + 1.0 / msb);
        const double mum = 1.0 / (1.0 / msq - 1.0 / msb);
        const double bbx2 = 0.5 * (bb2 + bx2);
        const double mb = Mb / unitres;
        const double mx = Mc / unitres;
        const double tm = (mb - mx) * (mb - mx);
        if (t > tm) {
            t = 0.99 * tm;
        }

        const double wt = 1.0 + (tm - t) / (2.0 * mbb * mbx);
        const double mqm = 0.1;

        const double r2 = 3.0 / (4.0 * msb * msq) + 3 * msd * msd / (2 * mbb * mbx * bbx2) +
                          (16.0 / (mbb * mbx * (33.0 - 2.0 * nfp))) * log(Getas(mqm, mqm) / Getas(msq, msq));
        const double ai = -1.0 * (6.0 / (33.0 - 2.0 * nf));

        const double cji = pow((Getas(msb, msb) / Getas(msq, msq)), ai);
        const double zji = msq / msb;

        const double gammaji = GetGammaji(zji);
        const double chiji = -1.0 - (gammaji / (1 - zji));

        const double betaji_g = (2.0 / 3.0) + gammaji;
        const double betaji_f = (-2.0 / 3.0) + gammaji;
        const double betaji_appam = -1.0 - chiji + (4.0 / (3.0 * (1.0 - zji))) +
                                    (2.0 * (1 + zji) * gammaji / (3.0 * (1.0 - zji) * (1.0 - zji)));

        const double betaji_apmam = (1.0 / 3.0) - chiji - (4.0 / (3.0 * (1.0 - zji))) -
                                    (2.0 * (1 + zji) * gammaji / (3.0 * (1.0 - zji) * (1.0 - zji))) + gammaji;

        const double r_g = cji * (1 + (betaji_g * Getas(msq, sqrt(mb * msq)) / pi));
        const double r_f = cji * (1 + (betaji_f * Getas(msq, sqrt(mb * msq)) / pi));
        const double r_apmam = cji * (1 + (betaji_apmam * Getas(msq, sqrt(mb * msq)) / pi));


        const double f3 = sqrt(mtx / mtb) * pow(sqrt(bx2 * bb2) / bbx2, 1.5) /
                          ((1.0 + r2 * (tm - t) / 12.0) * (1.0 + r2 * (tm - t) / 12.0));

        const double f3f = sqrt(mbx * mbb / (mtx * mtb)) * f3;
        const double f3g = sqrt(mtx * mtb / (mbx * mbb)) * f3;
        const double f3appam = sqrt(mtb * mtb * mtb * mbx / (mbb * mbb * mbb * mtx)) * f3;
        const double f3apmam = sqrt(mtx * mtb / (mbx * mbb)) * f3;

        const double f = cf * mtb * (1 + wt + msd * (wt - 1) / (2 * mup)) * f3f * r_f;
        const double g = 0.5 * (1 / msq - msd * bb2 / (2 * mum * mtx * bbx2)) * f3g * r_g;

        const double appam = cji *
                             (msd * bx2 * (1 - msd * bx2 / (2 * mtb * bbx2)) / ((1 + wt) * msq * msb * bbx2) -
                              betaji_appam * Getas(msq, sqrt(msq * mb)) / (mtb * pi)) *
                             f3appam;

        const double apmam =
            -1.0 *
            (mtb / msb - msd * bx2 / (2 * mup * bbx2) +
             wt * msd * mtb * bx2 * (1 - msd * bx2 / (2 * mtb * bbx2)) / ((wt + 1) * msq * msb * bbx2)) *
            f3apmam * r_apmam / mtx;

        const double ap = 0.5 * (appam + apmam);
        const double am = 0.5 * (appam - apmam);

        // const double vf = (g)*(Mb+Mc);
        // const double a1f = (f)/(Mb+Mc);
        // const double a2f = -1.0*(ap)*(Mb+Mc);
        // const double a3f = ((Mb+Mc)/(2.0*Mc))*a1f-((Mb-Mc)/(2.0*Mc))*a2f;
        // const double a0f = a3f + ( (Sqq*am)/(2.0*Mc));

        // resulting FF
        //  const double fv1=(Mb+Mc)/sqrt(Mb*Mc)*a1f;
        //  const double fv2=-(Mb/Mc)*(sqrt(Mb+Mc))*((a2f/(Mb+Mc))+((2.0*Mc/Sqq)*(a3f-a0f)));
        //  const double fv3=-sqrt(Mb*Mc)*((a2f/(Mb+Mc))-((2.0*Mc/Sqq)*(a3f-a0f)));
        //  const double fa=-2.*sqrt(Mb*Mc)*vf/(Mb+Mc);

        // Temporary NP form factors
        const double MbcSqq = pow(mb + mx, 2.) - t;
        const double fs = -2 * mx * f / MbcSqq;
        //        double Fpt = -(Mb + Mc) * f / MbcSqq;
        //        double Fmt = (Mb - Mc) * f / MbcSqq;
        //        double Fzt = 0.;

        // set elements, mapping to amplitude FF basis
        // Fs
        result.element({0}) = fs;
        // Ff
        result.element({1}) = f * unitres;
        // Fg
        result.element({2}) = g / unitres;
        // Fm
        result.element({3}) = am / unitres;
        // Fp
        result.element({4}) = ap / unitres;
        // Fzt
        //  result.element({5}) = 0;
        // Fmt
        //  result.element({6}) = 0;
        // Fpt
        //  result.element({7}) = 0;
    }

    std::unique_ptr<FormFactorBase> FFBtoDstarISGW2::clone(const std::string& label) {
        MAKE_CLONE(FFBtoDstarISGW2, label);
    }

} // namespace Hammer
