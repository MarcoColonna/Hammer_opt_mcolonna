///
/// @file  Loki.hh
/// @brief Double dispatcher code based on examples from A.Alexandrescu's book
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_TOOLS_LOKI
#define HAMMER_TOOLS_LOKI

#include <cstdint>
#include <type_traits>

namespace Loki {

    template <typename Test, template <typename...> class Ref>
    struct isSpecialization : std::false_type {};

    template <template <typename...> class Ref, typename... Args>
    struct isSpecialization<Ref<Args...>, Ref> : std::true_type {};


    template <typename... Ts>
    struct TypeList {
        static constexpr std::size_t size{sizeof...(Ts)};
    };

    using NullType = TypeList<>;

    namespace TL {

        template <typename, typename>
        struct IndexOf {};

        template <typename T, typename... Ts>
        struct IndexOf<T, TypeList<T, Ts...>> : std::integral_constant<std::size_t, 0> {};

        template <typename T, typename TOther, typename... Ts>
        struct IndexOf<T, TypeList<TOther, Ts...>>
            : std::integral_constant<std::size_t, 1 + IndexOf<T, TypeList<Ts...>>::value> {};

    } // namespace TL

    template <class Executor, class BaseLhs, class TypesLhs, bool symmetric = true, class BaseRhs = BaseLhs,
              class TypesRhs = TypesLhs, typename ResultType = void>
    class StaticDoubleDispatcher {

    private:

        template <bool swapArgs, class SomeLhs, class SomeRhs>
        struct InvocationTraits {
            static ResultType DoDispatch(SomeLhs& lhs, SomeRhs& rhs, Executor& exec) {
                return exec(lhs, rhs);
            }
        };
        template <class SomeLhs, class SomeRhs>
        struct InvocationTraits<true, SomeLhs, SomeRhs> {
            static ResultType DoDispatch(SomeLhs& lhs, SomeRhs& rhs, Executor& exec) {
                return exec(rhs, lhs);
            }
        };

    public:

        template <class SomeLhs>
        static ResultType DispatchRhs(SomeLhs& lhs, BaseRhs& rhs, Executor& exec, NullType /*unused*/) {
            return exec.error(lhs, rhs);
        }

        template <class SomeLhs, class Head, class... Ts>
        static ResultType DispatchRhs(SomeLhs& lhs, BaseRhs& rhs, Executor& exec, TypeList<Head, Ts...> /*unused*/) {
            if (Head* p2 = dynamic_cast<Head*>(&rhs)) {
                enum : std::uint8_t {
                    swapArgs = symmetric &&
                               int(TL::IndexOf<Head, TypesRhs>::value) < int(TL::IndexOf<SomeLhs, TypesLhs>::value)
                };

                using CallTraits = InvocationTraits<swapArgs, SomeLhs, Head>;

                return CallTraits::DoDispatch(lhs, *p2, exec);
            }
            return DispatchRhs(lhs, rhs, exec, TypeList<Ts...>());
        }

        static ResultType DispatchLhs(BaseLhs& lhs, BaseRhs& rhs, Executor& exec, NullType /*unused*/) {
            return exec.error(lhs, rhs);
        }

        template <class Head, class... Ts>
        static ResultType DispatchLhs(BaseLhs& lhs, BaseRhs& rhs, Executor& exec, TypeList<Head, Ts...> /*unused*/) {
            if (Head* p1 = dynamic_cast<Head*>(&lhs)) {
                return DispatchRhs(*p1, rhs, exec, TypesRhs());
            }
            return DispatchLhs(lhs, rhs, exec, TypeList<Ts...>());
        }

        static ResultType Go(BaseLhs& lhs, BaseRhs& rhs, Executor& exec) {
            return DispatchLhs(lhs, rhs, exec, TypesLhs());
        }
    };

    template <class Executor, class BaseLhs, class TypesLhs, typename ResultType = void>
    class StaticSingleDispatcher {

    private:

        template <class SomeLhs>
        struct InvocationTraits {
            static ResultType DoDispatch(SomeLhs& lhs, Executor& exec) {
                return exec(lhs);
            }
        };

    public:

        static ResultType Dispatch(BaseLhs& lhs, Executor& exec, NullType /*unused*/) {
            return exec.error(lhs);
        }

        template <class Head, class... Ts>
        static ResultType Dispatch(BaseLhs& lhs, Executor& exec, TypeList<Head, Ts...> /*unused*/) {
            if (Head* p1 = dynamic_cast<Head*>(&lhs)) {
                using CallTraits = InvocationTraits<Head>;

                return CallTraits::DoDispatch(*p1, exec);
            }
            return Dispatch(lhs, exec, TypeList<Ts...>());
        }

        static ResultType Go(BaseLhs& lhs, Executor& exec) {
            return Dispatch(lhs, exec, TypesLhs());
        }
    };
} // namespace Loki

#endif
