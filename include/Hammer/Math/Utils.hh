///
/// @file  Utils.hh
/// @brief Hammer math utilities class
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_MATH_UTILS
#define HAMMER_MATH_UTILS

#include <algorithm>
#include <complex>
#include <limits>
#include <climits>
#include <cstdint>
#include <cstring>
#include <array>

namespace Hammer {

    inline constexpr double precision = 0.001;

    /// @brief
    /// @param[in] val
    /// @return
    [[nodiscard]] inline bool isZero(const std::complex<double>& val) {
        return ((fabs(val.real()) < std::numeric_limits<double>::min()) &&
                (fabs(val.imag()) < std::numeric_limits<double>::min()));
    }

    /// @brief
    /// @param[in] val
    /// @return
    [[nodiscard]] inline bool isZero(double val) {
        return (fabs(val) < std::numeric_limits<double>::min());
    }

    /// @brief
    /// @param[in] val1
    /// @param[in] val2
    /// @return
    [[nodiscard]] inline bool fuzzyLess(double val1, double val2) {
        return (val1 - val2 < -1. * std::max(precision, std::numeric_limits<double>::min()));
    }

    /// @brief
    /// @param[in] val1
    /// @param[in] val2
    /// @return
    [[nodiscard]] double compareVals(double val1, double val2);

    /// @brief
    /// @param[in] val1
    /// @param[in] val2
    /// @return
    [[nodiscard]] std::complex<double> compareVals(const std::complex<double>& val1, const std::complex<double>& val2);


    [[nodiscard]] inline double regularize(double regularVal, double problematicValue,
                                           double delta = std::numeric_limits<double>::min(), int direction = 0) {
        if (fabs(problematicValue - regularVal) > delta) {
            return regularVal;
        }
        double signedSide = (problematicValue > regularVal) ? -1.0 : 1.0;
        return problematicValue + (delta * (direction == 0 ? signedSide : 1. * direction));
    }

    template <typename T>
    typename std::enable_if_t<std::is_arithmetic_v<typename std::remove_reference_t<T>>, double> toDouble(T value) {
        return static_cast<double>(value);
    }

    template <typename T>
    typename std::enable_if_t<
        std::is_same_v<typename std::remove_reference_t<typename std::remove_cv_t<T>>, std::string>, double>
    toDouble(T value) {
        return stod(value);
    }

    template <typename T>
    typename std::enable_if_t<
        std::is_same_v<typename std::remove_reference_t<typename std::remove_cv_t<T>>, std::complex<double>>, double>
    toDouble(T value) {
        return value.real();
    }

    template <typename T>
    typename std::enable_if_t<std::is_unsigned_v<T> && std::is_integral_v<T>, T> minPadding(T value) {
        if (value == 0u) {
            return value;
        }
        if (value == 1u) {
            return 0;
        }
        T pad = 0u;
        --value;
        while (value != 0u) {
            value = static_cast<T>(value >> 1);
            ++pad;
        }
        return pad;
    }

    // taken from https://stackoverflow.com/a/42138465
    using Fp_info = std::numeric_limits<double>;

    // taken from https://stackoverflow.com/a/42138465
    [[nodiscard]] inline auto is_ieee754_nan(double const x) -> bool {
        static constexpr bool is_claimed_ieee754 = Fp_info::is_iec559;
        static constexpr int n_bits_per_byte = CHAR_BIT;
        using Byte = unsigned char;

        static_assert(is_claimed_ieee754, "!");
        static_assert(n_bits_per_byte == 8, "!");
        static_assert(sizeof(x) == sizeof(uint64_t), "!");

#ifdef _MSC_VER
        uint64_t const bits = reinterpret_cast<uint64_t const&>(x);
#else
        std::array<Byte, sizeof(x)> bytes{};
        memcpy(bytes.data(), &x, sizeof(x));
        uint64_t int_value;
        memcpy(&int_value, bytes.data(), sizeof(x));
        uint64_t const& bits = int_value;
#endif

        static constexpr uint64_t sign_mask = 0x8000000000000000;
        static constexpr uint64_t exp_mask = 0x7FF0000000000000;
        static constexpr uint64_t mantissa_mask = 0x000FFFFFFFFFFFFF;

        (void) sign_mask;
        return (bits & exp_mask) == exp_mask and (bits & mantissa_mask) != 0;
    }

    [[nodiscard]] inline auto is_ieee754_nan(std::complex<double> x) -> bool {
        return is_ieee754_nan(x.real()) || is_ieee754_nan(x.imag());
    }

} // namespace Hammer

#endif
