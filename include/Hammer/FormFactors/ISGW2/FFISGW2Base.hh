///
/// @file  FFISGW2Base.hh
/// @brief Hammer base class for ISGW2 form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_ISGW2_BASE
#define HAMMER_FF_ISGW2_BASE

#include "Hammer/FormFactorBase.hh"

namespace Hammer {

    /// @brief Base class for ISGW2 form factors
    /// implementation matched to EvtGen
    ///
    /// Provides the rate contents as tensor...
    ///
    /// @ingroup FormFactors
    class FFISGW2Base : public FF1to1Base {

    public:

        FFISGW2Base();

        FFISGW2Base(const FFISGW2Base& other) = default;
        FFISGW2Base& operator=(const FFISGW2Base& other) = delete;
        FFISGW2Base(FFISGW2Base&& other) = delete;
        FFISGW2Base& operator=(FFISGW2Base&& other) = delete;
        ~FFISGW2Base() override = default;

    protected:

        /// @brief
        /// @param[in] z
        /// @return
        static double GetGammaji(double z);

        /// @brief
        /// @param[in] mq1
        /// @param[in] mq2
        /// @return
        static double Getas(double mq1, double mq2);

        /// @brief
        /// @param[in] m
        /// @return
        static double Getas(double m);

        /// @brief
        /// @param[in] point
        /// @param[in] masses
        /// @return
        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) override = 0;

        /// @brief
        void defineSettings() override = 0;

        void addRefs() const override;

        double msb{0.};
        double msd{0.};
        double bb2{0.};
        double mbb{0.};
        double nf{0.};
        double cf{0.};
        double msq{0.};
        double bx2{0.};
        double mbx{0.};
        double nfp{0.};
    };

} // namespace Hammer

#endif
