///
/// @file  TestFourMomentum.cc
/// @brief Tests for FourMomentum
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "Hammer/Math/FourMomentum.hh"
#include "gtest/gtest.h"

#include <cmath>
#include <array>


namespace Hammer {

    TEST(FourMomentumTest, AssignmentAccess) {
        FourMomentum v = FourMomentum::fromPtEtaPhiM(10., 1.5, 0., 91.2);
        EXPECT_DOUBLE_EQ(v.pt(), 10.);
        EXPECT_DOUBLE_EQ(v.eta(), 1.5);
        EXPECT_DOUBLE_EQ(v.phi(), 0.);
        EXPECT_DOUBLE_EQ(v.mass(), 91.2);
        FourMomentum v2 = FourMomentum::fromEtaPhiME(1.5, 0., 91.2, 100.);
        EXPECT_DOUBLE_EQ(v2.E(), 100.);
        EXPECT_DOUBLE_EQ(v2.eta(), 1.5);
        EXPECT_DOUBLE_EQ(v2.phi(), 0.);
        EXPECT_DOUBLE_EQ(v2.mass(), 91.2);
        FourMomentum v3{{{20., 2., 3., 5.}}};
        EXPECT_DOUBLE_EQ(v3.E(), 20.);
        EXPECT_DOUBLE_EQ(v3.px(), 2.);
        EXPECT_DOUBLE_EQ(v3.py(), 3.);
        EXPECT_DOUBLE_EQ(v3.pz(), 5.);
        FourMomentum v4 = FourMomentum::fromPM(2., 3., 5., 91.2);
        EXPECT_DOUBLE_EQ(v4.px(), 2.);
        EXPECT_DOUBLE_EQ(v4.py(), 3.);
        EXPECT_DOUBLE_EQ(v4.pz(), 5.);
        EXPECT_DOUBLE_EQ(v4.mass(), 91.2);
        v4.setE(20.);
        EXPECT_EQ(v4.E(), 20.);
    }

    TEST(FourMomentumTest, Access) {
        // mass2, mass, p2, p, pt for a null-mass vector (E^2 = px^2 + py^2)
        FourMomentum v{5.0, 3.0, 4.0, 0.0};
        EXPECT_NEAR(v.mass2(), 0., 1e-10);
        EXPECT_NEAR(v.p2(), 25., 1e-10);
        EXPECT_NEAR(v.p(), 5., 1e-10);
        EXPECT_NEAR(v.pt(), 5., 1e-10);

        // pVec returns the spatial 3-vector
        auto pvec = v.pVec();
        EXPECT_NEAR(pvec[0], 3., 1e-10);
        EXPECT_NEAR(pvec[1], 4., 1e-10);
        EXPECT_NEAR(pvec[2], 0., 1e-10);

        // dot product: E*E - px*px - py*py - pz*pz
        FourMomentum v3{5.0, 3.0, 4.0, 0.0};
        FourMomentum v4{5.0, -3.0, -4.0, 0.0};
        EXPECT_NEAR(v3.dot(v4), 5. * 5. - 3. * (-3.) - 4. * (-4.) - 0., 1e-10); // = 25+9+16 = 50

        // setPx / setPy / setPz chaining
        FourMomentum v5{10., 1., 2., 3.};
        v5.setPx(5.).setPy(6.).setPz(7.);
        EXPECT_NEAR(v5.px(), 5., 1e-10);
        EXPECT_NEAR(v5.py(), 6., 1e-10);
        EXPECT_NEAR(v5.pz(), 7., 1e-10);

        // Arithmetic operators (via boost::operators)
        FourMomentum va{5., 1., 2., 3.};
        FourMomentum vb{3., 0., 1., 1.};

        auto vc = va + vb;
        EXPECT_NEAR(vc.E(), 8., 1e-10);
        EXPECT_NEAR(vc.px(), 1., 1e-10);
        EXPECT_NEAR(vc.py(), 3., 1e-10);
        EXPECT_NEAR(vc.pz(), 4., 1e-10);

        auto vd = va - vb;
        EXPECT_NEAR(vd.E(), 2., 1e-10);
        EXPECT_NEAR(vd.px(), 1., 1e-10);

        auto ve = va * 2.0;
        EXPECT_NEAR(ve.E(), 10., 1e-10);
        EXPECT_NEAR(ve.px(), 2., 1e-10);

        auto vf = va / 2.0;
        EXPECT_NEAR(vf.E(), 2.5, 1e-10);
        EXPECT_NEAR(vf.px(), 0.5, 1e-10);

        auto vneg = -va;
        EXPECT_NEAR(vneg.E(), -5., 1e-10);
        EXPECT_NEAR(vneg.px(), -1., 1e-10);
        EXPECT_NEAR(vneg.py(), -2., 1e-10);
        EXPECT_NEAR(vneg.pz(), -3., 1e-10);

        // PFlip: flips only the 3-momentum, not the energy
        auto vflip = va.PFlip();
        EXPECT_NEAR(vflip.E(), va.E(), 1e-10);
        EXPECT_NEAR(vflip.px(), -va.px(), 1e-10);
        EXPECT_NEAR(vflip.py(), -va.py(), 1e-10);
        EXPECT_NEAR(vflip.pz(), -va.pz(), 1e-10);
    }

    TEST(FourMomentumTest, Angles) {
        FourMomentum v1{5., 0., 0., 5.}; // pointing in z direction -> theta = 0
        EXPECT_NEAR(v1.theta(), 0., 1e-10);

        FourMomentum v2{5., 5., 0., 0.}; // pointing in x direction -> theta = pi/2
        EXPECT_NEAR(v2.theta(), M_PI / 2., 1e-10);

        EXPECT_NEAR(v1.phi(), 0., 1e-10); // pz-only -> phi undefined but atan2(0,0)=0
        EXPECT_NEAR(v2.phi(), 0., 1e-10); // px > 0, py = 0 -> phi = 0

        FourMomentum v3{5., 0., 5., 0.}; // pointing in y direction -> phi = pi/2
        EXPECT_NEAR(v3.phi(), M_PI / 2., 1e-10);

        EXPECT_NEAR(v2.eta(), 0., 1e-10);

        FourMomentum v4{100., 0., 0., 99.};
        EXPECT_GT(v4.eta(), 3.); // large positive pseudorapidity

        FourMomentum v6 = FourMomentum::fromPM(0., 0., 0., 5.); // at rest, E = 5 = mass
        EXPECT_NEAR(v6.rapidity(), 0., 1e-10);

        FourMomentum v7 = FourMomentum::fromPM(0., 0., 4., 3.); // E=5, pz=4, mass=3
        EXPECT_NEAR(v7.rapidity(), std::acosh(5. / 3.), 1e-10);
    }

    TEST(FourMomentumTest, Boosts) {
        FourMomentum v1{4., 2., 0., 0.};
        auto bv = v1.boostVector();
        EXPECT_NEAR(bv[0], 0.5, 1e-10);
        EXPECT_NEAR(bv[1], 0., 1e-10);
        EXPECT_NEAR(bv[2], 0., 1e-10);

        EXPECT_NEAR(v1.beta(), 0.5, 1e-10);

        FourMomentum massive{5., 3., 0., 0.}; // mass = sqrt(25-9) = 4
        double expected_gamma = 5. / 4.;
        EXPECT_NEAR(massive.gamma(), expected_gamma, 1e-10);

        FourMomentum pB{5.280, 0., 0., 0.};
        FourMomentum pDs{2.66, 1.39, 0.40, 0.97};
        FourMomentum pDs_copy = pDs;
        pDs_copy.boostToRestFrameOf(pB);
        EXPECT_NEAR(pDs_copy.E(), pDs.E(), 1e-8);
        EXPECT_NEAR(pDs_copy.px(), pDs.px(), 1e-8);
        EXPECT_NEAR(pDs_copy.py(), pDs.py(), 1e-8);
        EXPECT_NEAR(pDs_copy.pz(), pDs.pz(), 1e-8);

        FourMomentum parent{5., 2., 0., 0.}; // mass = sqrt(25-4) = sqrt(21)
        FourMomentum daughter{2., 1., 0., 0.};
        FourMomentum roundtrip = daughter;
        roundtrip.boostToRestFrameOf(parent);
        roundtrip.boostFromRestFrameOf(parent);
        EXPECT_NEAR(roundtrip.E(), daughter.E(), 1e-10);
        EXPECT_NEAR(roundtrip.px(), daughter.px(), 1e-10);
        EXPECT_NEAR(roundtrip.py(), daughter.py(), 1e-10);
        EXPECT_NEAR(roundtrip.pz(), daughter.pz(), 1e-10);

        FourMomentum parent2{5., 2., 0., 0.};
        std::array<double, 3> bvec = parent2.boostVector();
        FourMomentum daughter2{3., 1., 0., 0.};
        FourMomentum daughter2_member = daughter2;
        daughter2_member.boostToRestFrameOf(bvec);
        FourMomentum daughter2_free = boostToRestFrameOf(daughter2, bvec);
        EXPECT_NEAR(daughter2_member.E(), daughter2_free.E(), 1e-10);
        EXPECT_NEAR(daughter2_member.px(), daughter2_free.px(), 1e-10);

        FourMomentum daughter3{3., 1., 0., 0.};
        FourMomentum daughter3_free = boostToRestFrameOf(daughter3, bvec[0], bvec[1], bvec[2]);
        EXPECT_NEAR(daughter3_free.E(), daughter2_free.E(), 1e-10);
        EXPECT_NEAR(daughter3_free.px(), daughter2_free.px(), 1e-10);

        FourMomentum daughter4 = boostToRestFrameOf(daughter2, parent2);
        EXPECT_NEAR(daughter4.E(), daughter2_member.E(), 1e-10);
        EXPECT_NEAR(daughter4.px(), daughter2_member.px(), 1e-10);
    }

} // namespace Hammer
