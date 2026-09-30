///
/// @file  FFBLPRBase.hh
/// @brief Hammer base class for BLPR form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BLPR_BASE
#define HAMMER_FF_BLPR_BASE

#include <cstdint>
#include <vector>
#include <array>
#include <functional>

#include "Hammer/FormFactorBase.hh"
#include "Hammer/FormFactors/AsCoeffs.hh"

namespace Hammer {

    /// @brief Base class for BLPR form factors
    ///
    /// @ingroup FormFactors
    class FFBLPRBase : public FF1to1Base {

    public:

        FFBLPRBase();

        FFBLPRBase(const FFBLPRBase& other) = default;
        FFBLPRBase& operator=(const FFBLPRBase& other) = delete;
        FFBLPRBase(FFBLPRBase&& other) = delete;
        FFBLPRBase& operator=(FFBLPRBase&& other) = delete;
        ~FFBLPRBase() override = default;

    protected:

        /// @brief
        /// @param[in] point
        /// @param[in] masses
        /// @return
        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) override = 0;

        /// @brief
        void defineSettings() override;

        void addRefs() const override;

        std::pair<double, double> Xi(double w) const; // numerator and denominator separately

        double ZofW(double w) const;

    private:

        void fillIWs();
        void fillLis();

    protected:

        AsCoeffs _asCorrections;

        enum IwLabel : uint8_t { ETA = 0, CHI2, CHI3 };

        std::array<std::function<double(double)>, 3> _mIWs;
        std::array<std::function<double(double)>, 7> _mLi;

        mutable double _a{0.};
        mutable double _zBC{0.};

        mutable double _mCv1z{0.};
        mutable double _mCv2z{0.};
        mutable double _mCv3z{0.};
        mutable double _mCv1zp{0.};
        mutable double _mCv2zp{0.};
        mutable double _mCv3zp{0.};
        mutable double _mCv1zpp{0.};
        mutable double _mCv2zpp{0.};
        mutable double _mCv3zpp{0.};

    private:

        mutable bool _initCs = false;
        mutable double _w0{0.};
    };

} // namespace Hammer

#endif
