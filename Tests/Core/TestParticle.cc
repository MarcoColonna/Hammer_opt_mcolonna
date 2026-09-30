///
/// @file  TestParticle.cc
/// @brief Tests for Particle
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "Hammer/Particle.hh"
#include "gtest/gtest.h"


namespace Hammer {

    TEST(ParticleTest, CreationAssignment) {
        Particle p;
        const Particle* pp = &p;
        EXPECT_TRUE(pp == &(p.setPdgId(22)));
        Particle p2(FourMomentum(.139, 0., 0., 0.), 131);
        const Particle* pp2 = &p2;
        EXPECT_TRUE(pp2 == &(p2.setMomentum(FourMomentum(25., 0., 0., -25.))));
    }

    TEST(ParticleTest, Access) {
        Particle p(FourMomentum(.139, 0., 0., 0.), 131);
        EXPECT_NEAR(p.momentum().E(), .139, .001);
        EXPECT_EQ(p.pdgId(), 131);
        p.setMomentum(FourMomentum(25., 0., 0., -25.));
        EXPECT_NEAR(p.p().pz(), -25., 0.001);
    }

} // namespace Hammer
