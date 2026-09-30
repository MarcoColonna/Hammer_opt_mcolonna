///
/// @file  FFBLPRXPBase.hh
/// @brief Hammer base class for BLPR-XP form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BLPRXP_BASE
#define HAMMER_FF_BLPRXP_BASE

#ifndef FFBLPRXPBASE_FRIENDS
#define FFBLPRXPBASE_FRIENDS typedef bool DummyFriend
#endif

#include <cstdint>
#include <vector>
#include <array>
#include <functional>

#include "Hammer/FormFactorBase.hh"
#include "Hammer/FormFactors/AsCoeffs.hh"

namespace Hammer {

    /// @brief Base class for BLPR-XP form factors
    ///
    /// @ingroup FormFactors
    class FFBLPRXPBase : public FF1to1Base {

        FFBLPRXPBASE_FRIENDS;

    public:

        FFBLPRXPBase();

        FFBLPRXPBase(const FFBLPRXPBase& other) = default;
        FFBLPRXPBase& operator=(const FFBLPRXPBase& other) = delete;
        FFBLPRXPBase(FFBLPRXPBase&& other) = delete;
        FFBLPRXPBase& operator=(FFBLPRXPBase&& other) = delete;
        ~FFBLPRXPBase() override = default;

    protected:

        /// @brief
        /// @param[in] point
        /// @param[in] masses
        /// @return
        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) override = 0;

        /// @brief
        void defineSettings() override;

        void addRefs() const override;

        double Xi(double w) const; // numerator and denominator separately

        double ZofW(double w) const;

        void calcConstants();

    private:

        void fillIWs();
        void fillLisMis();

    protected:

        AsCoeffs _asCorrections;

        enum IwLabel : uint8_t { CHI2 = 0, CHI3, ETA, BETA1, BETA2, BETA3, PHI1, PHI1Q };

        std::array<std::function<double(double)>, 8> _mIWs;
        std::array<std::function<double(double)>, 7> _mLi1;
        std::array<std::function<double(double)>, 7> _mLi2;
        std::array<std::function<double(double)>, 25> _mMi;

        double _a{0.};
        double _zBC{0.};

        double _la2OverlaB2{0.};
        double _la1OverlaB2{0.};
        double _eb{0.};
        double _ec{0.};
        double _upsilonc{0.};
        double _upsilonb{0.};
        double _cmagc{0.};
        double _cmagb{0.};

    private:

        double _mb{0.};
        double _mc{0.};
        double _la2{0.};
        double _rho1{0.};
        double _mBBar{0.};
        double _mDBar{0.};
    };

} // namespace Hammer

#endif
