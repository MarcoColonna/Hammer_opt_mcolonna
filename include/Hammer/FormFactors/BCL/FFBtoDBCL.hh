///
/// @file  FFBtoDBCL.hh
/// @brief \f$ B \rightarrow D \f$ BCL form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BDBCL
#define HAMMER_FF_BDBCL

#include "Hammer/FormFactors/BCL/FFBCLBase.hh"

namespace Hammer {

    class FFBtoDBCL final : public FFBCLBase {

    public:

        FFBtoDBCL();

        FFBtoDBCL(const FFBtoDBCL& other) = default;
        FFBtoDBCL& operator=(const FFBtoDBCL& other) = delete;
        FFBtoDBCL(FFBtoDBCL&& other) = delete;
        FFBtoDBCL& operator=(FFBtoDBCL&& other) = delete;
        ~FFBtoDBCL() final = default;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;

        void addRefs() const final;
    };

} // namespace Hammer

#endif
