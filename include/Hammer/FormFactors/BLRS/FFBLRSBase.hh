///
/// @file  FFBLRSBase.hh
/// @brief Hammer base class for BLRS form factors
//

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BLRS_BASE
#define HAMMER_FF_BLRS_BASE

#include "Hammer/FormFactorBase.hh"
#include "Hammer/FormFactors/AsCoeffs.hh"

namespace Hammer {

    /// @brief Base class for BLR form factors
    ///
    /// @ingroup FormFactors
    class FFBLRSBase : public FF1to1Base {

    public:

        FFBLRSBase();

        FFBLRSBase(const FFBLRSBase& other) = default;
        FFBLRSBase& operator=(const FFBLRSBase& other) = delete;
        FFBLRSBase(FFBLRSBase&& other) = delete;
        FFBLRSBase& operator=(FFBLRSBase&& other) = delete;
        ~FFBLRSBase() override = default;

    protected:

        /// @brief
        /// @param[in] point
        /// @param[in] masses
        /// @return
        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) override = 0;

        /// @brief
        void defineSettings() override;

        void addRefs() const override;

        void calcConstants(const double& unitres);

        void fillKisMis(const double& unitres);

        AsCoeffs _asCorrections;

        std::array<std::function<double(double)>, 3> _mKi1;
        std::array<std::function<double(double)>, 3> _mKi2;
        // std::array<std::function<double(double)>, 5> _mMi;

        double _aS{0.};
        double _zBC{0.};
        double _mb1S{0.};
        double _mc1S{0.};

        double _eB0{0.};
        double _eC0{0.};
        double _eB{0.};
        double _eC{0.};

    private:

        double _dmbmc{0.};
        double _mLb{0.};
        double _mLc{0.};
        double _la1S{0.};
        double _b1{0.};
        double _b2{0.};
    };

} // namespace Hammer

#endif
