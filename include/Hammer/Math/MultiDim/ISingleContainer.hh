///
/// @file  ISingleContainer.hh
/// @brief Interface class for single container data structure
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_MATH_ISINGLECONTAINER
#define HAMMER_MATH_ISINGLECONTAINER


#include "Hammer/Math/MultiDim/IContainer.hh"

namespace Hammer::MultiDimensional {

    class ISingleContainer : public IContainer {

    public:

        ~ISingleContainer() override = default;

        class ItBase {
        public:

            virtual ~ItBase() = default;
            [[nodiscard]] virtual IContainer::ElementType value() const = 0;
            [[nodiscard]] virtual PositionType position() const = 0;
            virtual void next(int n = 1) = 0;
            [[nodiscard]] virtual bool isSame(const ItBase& other) const = 0;
            [[nodiscard]] virtual bool isAligned() const = 0;
            [[nodiscard]] virtual ptrdiff_t distanceFrom(const ItBase& other) const = 0;
        };

        using NonZeroIt = std::unique_ptr<ItBase>;

        [[nodiscard]] virtual NonZeroIt firstNonZero() const = 0;
        [[nodiscard]] virtual NonZeroIt endNonZero() const = 0;
    };

    inline void next(ISingleContainer::NonZeroIt it) {
        it->next();
    }

    inline bool operator==(const ISingleContainer::ItBase& lhs, const ISingleContainer::ItBase& rhs) {
        return lhs.isSame(rhs);
    }

    inline bool operator!=(const ISingleContainer::ItBase& lhs, const ISingleContainer::ItBase& rhs) {
        return !(lhs == rhs);
    }

} // namespace Hammer::MultiDimensional

#endif
