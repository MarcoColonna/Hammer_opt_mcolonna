///
/// @file  Units.hh
/// @brief Unit conversion factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_MATH_UNITS
#define HAMMER_MATH_UNITS

#include <map>
#include <string>

namespace Hammer {

    static constexpr double eV = 1.e-9;
    static constexpr double eV2 = eV * eV;
    static constexpr double keV = 1.e-6;
    static constexpr double keV2 = keV * keV;
    static constexpr double MeV = 1.e-3;
    static constexpr double MeV2 = MeV * MeV;
    static constexpr double GeV = 1.;
    static constexpr double GeV2 = GeV * GeV;
    static constexpr double TeV = 1.e3;
    static constexpr double TeV2 = TeV * TeV;

    // static const double cm =

    /// @brief Hammer class for dealing with units
    ///
    /// Organized as a singleton accessed by the `instance()` method
    ///
    /// @ingroup Math

    class Units {

    public:

        Units(const Units&) = delete;
        Units& operator=(const Units&) = delete;
        Units(Units&&) = delete;
        Units& operator=(Units&&) = delete;

    protected:

        Units() = default;
        ~Units() = default;

        /// @brief
        /// @return
        static Units* getUnitsInstance();

        /// @brief initializes the units map
        void init();

    public:

        /// @brief
        /// @return
        static Units& instance();

        /// @brief
        /// @return
        [[nodiscard]] double getUnitsRescalingToMC(const std::string& mcunits, const std::string& localunits) const;

    private:

        static Units* _theUnits;

        std::map<std::string, double> _unitfactors;
    };

} // namespace Hammer

#endif
