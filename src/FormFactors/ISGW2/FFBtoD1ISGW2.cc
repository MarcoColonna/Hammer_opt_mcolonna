///
/// @file  FFBtoD1ISGW2.cc
/// @brief \f$ B \rightarrow D_1 \f$ ISGW2 form factors
/// @brief Ported directly from EvtGen
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <cmath>

#include "Hammer/FormFactors/ISGW2/FFBtoD1ISGW2.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    FFBtoD1ISGW2::FFBtoD1ISGW2() {
        // Create tensor rank and dimensions
        IndexList dims = {8};
        string name{"FFBtoD1ISGW2"};

        setPrefix("BtoD**1");
        addProcessSignature(PID::BPLUS, {-PID::DSSD1});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BDSSD1})});

        addProcessSignature(PID::BZERO, {PID::DSSD1MINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BDSSD1})});

        setPrefix("BstoDs**1");
        addProcessSignature(PID::BS, {PID::DSSDS1MINUS});
        addTensor(Tensor{name, MD::makeEmptySparse(dims, {FF_BSDSSDS1})});

        setSignatureIndex();
    }

    void FFBtoD1ISGW2::defineSettings() {
        //_mFFErrNames = ;
        setPath(getFFErrPrefixGroup().get());
        setUnits("GeV");

        addSetting<bool>("SmearQ2", true); // EvtGen default smearing

        setInitialized();
    }

    void FFBtoD1ISGW2::evalAtPSPoint(const vector<double>& point, const vector<double>& masses) {
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
        // Ported from EvtGen EvtISGW2FF1P1
        switch (result.labels()[0]) {
        case FF_BDSSD1:
            msb = 5.2;
            msd = 0.33;
            bb2 = 0.431 * 0.431;
            mbb = 5.31;

            msq = 1.82;
            bx2 = 0.33 * 0.33;
            mbx = (5.0 * 2.46 + 3.0 * 2.42) / 8.0;
            nfp = 3.0;
            break;
        case FF_BSDSSDS1:
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
                          (16.0 / (mbb * mbx * (33.0 - 2.0 * nfp))) * log(Getas(mqm, mqm) / Getas(msq, msq));

        const double f5 =
            sqrt(mtx / mtb) * pow(sqrt(bx2 * bb2) / bbx2, 5.0 / 2.0) / (pow((1.0 + r2 * (tm - t) / 18.0), 3.0));

        const double f5v = f5 * pow((mbb / mtb), -0.5) * pow((mbx / mtx), -0.5);
        const double f5r = f5 * pow((mbb / mtb), 0.5) * pow((mbx / mtx), 0.5);
        const double f5sppsm = f5 * pow((mbb / mtb), -1.5) * pow((mbx / mtx), 0.5);
        const double f5spmsm = f5 * pow((mbb / mtb), -0.5) * pow((mbx / mtx), -0.5);

        const double vv =
            -msd * f5v / (2.0 * sqrt(3.0 * bb2) * mtx) * ((wt + 1) / 2.0 + bb2 * mtb / (2.0 * msd * msq * msb));
        const double rr = -2. * mtb * sqrt(bb2 / 3.0) * f5r *
                          (1.0 / msq + mtx * msd * (wt - 1) / (2.0 * bb2) *
                                           ((wt + 1) / (2.0 * msq) - msd * bb2 / (2.0 * mum * mtx * bbx2)));
        const double sppsm = -sqrt(3.0) * msd * f5sppsm / (2.0 * sqrt(bb2) * mtb) *
                             (1 - msd / (3.0 * msq) - msd * bb2 / (3.0 * bbx2) * (1.0 / (2.0 * mum) - 1.0 / mup));
        const double spmsm = -msd * f5spmsm / (2.0 * sqrt(3.0 * bb2) * mtx) *
                             ((2 - wt) * mtx / msq + msd * bb2 / bbx2 * (1.0 / (2.0 * mum) - 1.0 / mup));

        // FF Basis map
        const double sqmbmx = sqrt(mb * mx);
        const double Fv1 = rr / sqmbmx;
        const double Fv2 = mb * mb * sppsm / sqmbmx;
        const double Fv3 = spmsm * sqmbmx;
        const double Fa = 2 * vv * sqmbmx;

        const bool smear = (*getSetting<bool>("SmearQ2"));
        double smearQ2 = 1.;
        if (smear) { // Ad hoc smearing, following EvtGen EvtISGW2FF1P1
            const double q2max = (mb - mx) * (mb - mx);
            const double mbmean = this->masses()[0];
            const double mxmean = this->masses()[1];
            const double q2maxmean = (mbmean - mxmean) * (mbmean - mxmean);
            smearQ2 = sqrt(q2maxmean / q2max);
        }

        // Set elements
        result.element({1}) = Fv1 * smearQ2;
        result.element({2}) = Fv2 * smearQ2;
        result.element({3}) = Fv3 * smearQ2;
        result.element({4}) = Fa * smearQ2;
    }

    std::unique_ptr<FormFactorBase> FFBtoD1ISGW2::clone(const std::string& label) {
        MAKE_CLONE(FFBtoD1ISGW2, label);
    }

} // namespace Hammer
