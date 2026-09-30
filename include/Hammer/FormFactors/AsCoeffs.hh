///
/// @file  AsCoeffs.hh
/// @brief Hammer class for HQET \f$ \alpha_s \f$ corrections
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_ASCOEFFS
#define HAMMER_FF_ASCOEFFS

namespace Hammer {

    /// @brief Base class for HQET \f$ \alpha_s \f$ corrections
    ///
    /// @ingroup FormFactors
    class AsCoeffs {

    public:

        AsCoeffs() = default;

        AsCoeffs(const AsCoeffs& other) = default;
        AsCoeffs& operator=(const AsCoeffs& other) = delete;
        AsCoeffs(AsCoeffs&& other) = delete;
        AsCoeffs& operator=(AsCoeffs&& other) = delete;
        virtual ~AsCoeffs() = default;

        /// @brief
        /// @param[in] w
        /// @param[in] z
        /// @return
        double CS(double w, double z) const;

        /// @brief
        /// @param[in] w
        /// @param[in] z
        /// @return
        double CP(double w, double z) const;

        /// @brief
        /// @param[in] w
        /// @param[in] z
        /// @return
        double CV1(double w, double z) const;

        /// @brief
        /// @param[in] w
        /// @param[in] z
        /// @return
        double CV2(double w, double z) const;

        /// @brief
        /// @param[in] w
        /// @param[in] z
        /// @return
        double CV3(double w, double z) const;

        /// @brief
        /// @param[in] w
        /// @param[in] z
        /// @return
        double CA1(double w, double z) const;

        /// @brief
        /// @param[in] w
        /// @param[in] z
        /// @return
        double CA2(double w, double z) const;

        /// @brief
        /// @param[in] w
        /// @param[in] z
        /// @return
        double CA3(double w, double z) const;

        /// @brief
        /// @param[in] w
        /// @param[in] z
        /// @return
        double CT1(double w, double z) const;

        /// @brief
        /// @param[in] w
        /// @param[in] z
        /// @return
        double CT2(double w, double z) const;

        /// @brief
        /// @param[in] w
        /// @param[in] z
        /// @return
        double CT3(double w, double z) const;

        /// @brief
        /// @param[in] w
        /// @param[in] z
        /// @return
        double derCS(double w, double z) const;

        /// @brief
        /// @param[in] w
        /// @param[in] z
        /// @return
        double derCP(double w, double z) const;

        /// @brief
        /// @param[in] w
        /// @param[in] z
        /// @return
        double derCV1(double w, double z) const;

        /// @brief
        /// @param[in] w
        /// @param[in] z
        /// @return
        double derCV2(double w, double z) const;

        /// @brief
        /// @param[in] w
        /// @param[in] z
        /// @return
        double derCV3(double w, double z) const;

        /// @brief
        /// @param[in] w
        /// @param[in] z
        /// @return
        double derCA1(double w, double z) const;

        /// @brief
        /// @param[in] w
        /// @param[in] z
        /// @return
        double derCA2(double w, double z) const;

        /// @brief
        /// @param[in] w
        /// @param[in] z
        /// @return
        double derCA3(double w, double z) const;

        /// @brief
        /// @param[in] w
        /// @param[in] z
        /// @return
        double derCT1(double w, double z) const;

        /// @brief
        /// @param[in] w
        /// @param[in] z
        /// @return
        double derCT2(double w, double z) const;

        /// @brief
        /// @param[in] w
        /// @param[in] z
        /// @return
        double derCT3(double w, double z) const;

    private:

        // Initialization of common variables & functions
        void updateVars(double w, double z) const;

        // DiLog function
        double DiLog(double z) const;

        mutable double _wSq{1.};
        mutable double _sqrt1wSq{0.};
        mutable double _zSq{1.};
        mutable double _zCu{1.};
        mutable double _zm1Sq{0.};
        mutable double _zp1Sq{4.};
        mutable double _lnz{0.};
        mutable double _wZ{1.};
        mutable double _wP{1.};
        mutable double _wM{1.};
        mutable double _wmwZSq{0.};
        mutable double _lnwP{0.};
        mutable double _rW{1.};
        mutable double _mOmega{-1.};
        mutable double _w{1.};
        mutable double _z{1.};
    };

} // namespace Hammer

#endif
