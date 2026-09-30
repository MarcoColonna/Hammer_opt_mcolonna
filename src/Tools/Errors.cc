///
/// @file  Errors.cc
/// @brief Hammer errors reporting helper functions
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <iostream>

#include "Hammer/Exceptions.hh"

using namespace std;

namespace Assert {

    void HandleAssert(const char* message, const char* condition, const char* fileName, long lineNumber) {
        cerr << "Assert Failed: \"" << message << "\"" << '\n';
        cerr << "Condition: " << condition << '\n';
        cerr << "File: " << fileName << '\n';
        cerr << "Line: " << lineNumber << '\n';
        cerr << "Application now terminating";
    }

} // namespace Assert
