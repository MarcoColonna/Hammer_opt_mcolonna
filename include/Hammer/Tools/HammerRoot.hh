///
/// @file  HammerRoot.hh
/// @brief ROOT includes
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_HammerRoot_HH
#define HAMMER_HammerRoot_HH

#include "Hammer/Tools/HammerRoot.fhh"
#include "Hammer/Config/HammerConfig.hh"

#ifdef HAVE_ROOT


#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wnonportable-include-path"
#include INCLRCONFIG


#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TLorentzVector.h"
#include "Math/Vector4D.h"
#pragma clang diagnostic pop

#endif

#endif
