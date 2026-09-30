///
/// @file  FFBGLXBase.hh
/// @brief Hammer base class for BGL Extended form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BGLX_BASE
#define HAMMER_FF_BGLX_BASE

#include "Hammer/FormFactorBase.hh"

namespace Hammer {

    /// @brief Base class for BGL form factors
    ///
    /// @ingroup FormFactors
    class FFBGLXBase : public FF1to1Base {

    public:

        FFBGLXBase();

        FFBGLXBase(const FFBGLXBase& other) = default;
        FFBGLXBase& operator=(const FFBGLXBase& other) = delete;
        FFBGLXBase(FFBGLXBase&& other) = delete;
        FFBGLXBase& operator=(FFBGLXBase&& other) = delete;
        ~FFBGLXBase() override = default;

    protected:

        /// @brief
        /// @param[in] point
        /// @param[in] masses
        /// @return
        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) override = 0;

        /// @brief
        void defineSettings() override = 0;

        void addRefs() const override;
        
        enum FfName : uint8_t {FP = 0, F0, FT, G, F, F1, F2, T1, T2, T23};
        
        enum PjName : uint8_t {J1P = 0, J1M, J0P, J0M};
        
        enum Tq2Val : uint8_t {TQ2 = 0, TQ2M};
        
        void setTZ(const double& t, const double& tp, const double& tm, const double& t0);
        
        void fillOutersJ1P();
        
        void fillOutersJ1M();
        
        void fillOutersJ0M();
        
        void fillOutersJ0P();
        
        void fillBlaschkes(const std::vector<double>& BcStates, const PjName J);
        
        std::array<std::map<Tq2Val, double>, 10> _mPhi{}; 
        std::array<double, 4> _mPJ{};
        
        double _z{0.};
        double _zm{0.}; //z(q^2 = q^2 max)
        
    private:
         
        void fillOuterFactors(Tq2Val tq2 = TQ2);
        
        double _tp{0.};
        double _tm{0.};
        double _t0{0.};
        double _t{0.};

        double _mPhiA1{0.};
        double _mPhiA3{0.};
        double _mPhiB1{0.};
        double _mPhiB3{0.};
        double _mPhiC1{0.};
        double _mPhiC2{0.};
        double _mPhiPre{0.};
        
    };

} // namespace Hammer

#endif
