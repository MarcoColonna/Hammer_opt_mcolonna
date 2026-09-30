///
/// @file  FFBtoD1starISGW2.cc
/// @brief \f$ B \rightarrow D_1^* \f$ ISGW2 form factors
/// @brief Ported directly from EvtGen
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/FormFactors/ISGW2/FFBtoD1starISGW2.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFBtoD1starISGW2::FFBtoD1starISGW2() {
        // Create tensor rank and dimensions
        IndexList dims = {8};
        string name{"FFBtoD1starISGW2"};

        setPrefix("BtoD**1*");
        addProcessSignature(PID::BPLUS, {-PID::DSSD1STAR});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BDSSD1STAR})});

        addProcessSignature(PID::BZERO, {PID::DSSD1STARMINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BDSSD1STAR})});

        setPrefix("BstoDs**1*");
        addProcessSignature(PID::BS, {PID::DSSDS1STARMINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BSDSSDS1STAR})});

        setSignatureIndex();
    }

    void FFBtoD1starISGW2::defineSettings() {
        //_mFFErrNames = ;
        setPath(getFFErrPrefixGroup().get());
        setUnits("GeV");

        addSetting<bool>("SmearQ2", true); // EvtGen default smearing

        setInitialized();
    }

    void FFBtoD1starISGW2::evalAtPSPoint(const vector<double>& point, const vector<double>& masses) {
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
        // Ported from EvtGen EvtISGW2FF3P1
        switch (result.labels()[0]) {
        case FF_BDSSD1STAR:
            msb = 5.2;
            msd = 0.33;
            bb2 = 0.431 * 0.431;
            mbb = 5.31;

            msq = 1.82;
            bx2 = 0.33 * 0.33;
            mbx = (3.0 * 2.49 + 2.40) / 4.0;
            nfp = 3.0;
            break;
        case FF_BSDSSDS1STAR:
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

        const double f5q = f5 * pow((mbb / mtb), -0.5) * pow((mbx / mtx), -0.5);
        const double f5l = f5 * pow((mbb / mtb), 0.5) * pow((mbx / mtx), 0.5);
        const double f5cppcm = f5 * pow((mbb / mtb), -1.5) * pow((mbx / mtx), 0.5);
        const double f5cpmcm = f5 * pow((mbb / mtb), -0.5) * pow((mbx / mtx), -0.5);

        const double ql = f5q * sqrt(1.0 / 6.0) * msd / (sqrt(bb2) * mtx) * (1.0 - bb2 * mtb / (4.0 * msd * msq * msb));
        const double ll = f5l * sqrt(2.0 / 3.0) * mtb * sqrt(bb2) *
                          (1.0 / (2.0 * msq) - 3.0 / (2.0 * msb) +
                           msd * mtx * (wt - 1) / bb2 * (1.0 / msq - msd * bb2 / (2.0 * mum * mtx * bbx2)));
        const double cppcm = msd * msd * bx2 * f5cppcm / (sqrt(6.0) * mtb * msq * sqrt(bb2) * bbx2);
        const double cpmcm =
            -sqrt(2.0 / 3.0) * msd * f5cpmcm / (sqrt(bb2) * mtx) * (1 + msd * bx2 / (2.0 * msq * bbx2));

        // FF Basis map
        const double sqmbmx = sqrt(mb * mx);
        const double Gv1 = ll / sqmbmx;
        const double Gv2 = mb * mb * cppcm / sqmbmx;
        const double Gv3 = cpmcm * sqmbmx;
        const double Ga = 2 * ql * sqmbmx;

        const bool smear = (*getSetting<bool>("SmearQ2"));
        double smearQ2 = 1.;
        if (smear) { // Ad hoc smearing, following EvtGen EvtISGW2FF3P1
            const double q2max = (mb - mx) * (mb - mx);
            const double mbmean = this->masses()[0];
            const double mxmean = this->masses()[1];
            const double q2maxmean = (mbmean - mxmean) * (mbmean - mxmean);
            smearQ2 = sqrt(q2maxmean / q2max);
        }

        // Set elements
        result.element({1}) = Gv1 * smearQ2;
        result.element({2}) = Gv2 * smearQ2;
        result.element({3}) = Gv3 * smearQ2;
        result.element({4}) = Ga * smearQ2;
    }

    std::unique_ptr<FormFactorBase> FFBtoD1starISGW2::clone(const std::string& label) {
        MAKE_CLONE(FFBtoD1starISGW2, label);
    }

} // namespace Hammer
