///
/// @file  Exceptions.hh
/// @brief Hammer exception definitions
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_EXCEPTIONS
#define HAMMER_EXCEPTIONS

#include <stdexcept>
#include <string>

namespace Hammer {

    /// @brief Generic error class
    /// @ingroup Base
    class Error : public std::runtime_error {

    public:

        explicit Error(const std::string& message) noexcept : runtime_error(message) {
        }

        Error(const Error&) = default;
        Error(Error&&) = default;
        Error& operator=(const Error&) = default;
        Error& operator=(Error&&) = default;
        ~Error() noexcept override = default;
    };

    /// @brief Invalid index label error class
    /// @ingroup Base
    class IndexLabelError : public Error {

    public:

        explicit IndexLabelError(const std::string& message) noexcept : Error(message) {
        }
    };

    /// @brief Out-of-range error class
    /// @ingroup Base
    class RangeError : public Error {

    public:

        explicit RangeError(const std::string& message) noexcept : Error(message) {
        }
    };

    /// @brief Numerical issues error class
    /// @ingroup Base
    class NumericalError : public Error {

    public:

        explicit NumericalError(const std::string& message) noexcept : Error(message) {
        }
    };

    /// @brief Invalid phase space point error class
    /// @ingroup Base
    class PhaseSpaceError : public Error {

    public:

        explicit PhaseSpaceError(const std::string& message) noexcept : Error(message) {
        }
    };

    /// @brief Initialization error class
    /// @ingroup Base
    class InitializationError : public Error {

    public:

        explicit InitializationError(const std::string& message) noexcept : Error(message) {
        }
    };

} // namespace Hammer

namespace Assert {
    void HandleAssert(const char* message, const char* condition, const char* fileName, long lineNumber);
}

#ifndef NDEBUG

#define ASSERT_MSG(x, msg)                                     \
    do {                                                       \
        if (!(x)) {                                            \
            Assert::HandleAssert(msg, #x, __FILE__, __LINE__); \
            abort();                                           \
        }                                                      \
    } while (0)

#define ASSERT(x)                                                       \
    do {                                                                \
        if (!(x)) {                                                     \
            Assert::HandleAssert("Assertion!", #x, __FILE__, __LINE__); \
            abort();                                                    \
        }                                                               \
    } while (0)

#else // Disabled Asserts

#define ASSERT(x)            \
    do {                     \
        (void) sizeof(true); \
    } while (0)

#define ASSERT_MSG(x, msg)                       \
    do {                                         \
        (void) sizeof(true), (void) sizeof(msg); \
    } while (0)

#endif

#endif
