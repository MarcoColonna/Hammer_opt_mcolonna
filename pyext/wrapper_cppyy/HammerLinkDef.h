//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#ifdef __CLING__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;
#pragma link C++ nestedtypedef;

#pragma link C++ class Hammer::FourMomentum;
#pragma link C++ class Hammer::PID;
#pragma link C++ class Hammer::Particle;
#pragma link C++ class Hammer::Process;
#pragma link C++ class Hammer::Log;
#pragma link C++ class Hammer::Hammer;
#pragma link C++ class Hammer::IOBuffer;
#pragma link C++ class Hammer::IOBuffers;

#pragma link C++ struct Hammer::BinContents;
#pragma link C++ struct Hammer::HistoInfo;
#pragma link C++ enum Hammer::WTerm;
#pragma link C++ enum Hammer::PAction;
#pragma link C++ enum Hammer::RecordType;

#ifdef HAVE_ROOT
#pragma link C++ class Hammer::RootIOBuffer;
#endif

#endif