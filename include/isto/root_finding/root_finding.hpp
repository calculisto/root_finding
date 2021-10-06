#include <utility>
#include <type_traits>
#include <vector>
#include <tuple>
#pragma once
    namespace isto::root_finding
{
    struct
options_t
{
        int
    max_iter = 100;
};

    struct
no_convergence_e
{};

    struct
zero_derivative_e
{};
/*
    namespace
info
{
        struct
    none_t
    {};
        constexpr auto
    none = none_t {};

        struct
    iterations_t
    {};
        constexpr auto
    iterations = iterations_t {};
        struct
    iterations_data_t
    {
            int
        iteration_count;
    };

        struct
    convergence_t
    {};
        constexpr auto
    convergence = convergence_t {};
        template <class Value, class FunctionResult, class DerivativeResult>
        struct
    convergence_data_newton_t
    {
            std::vector <std::tuple <
                  Value
                , FunctionResult
                , DerivativeResult
            >>
        convergence;
    };
        template <class Value, class FunctionResult>
        struct
    convergence_data_zhang_t
    {
            std::vector <std::tuple <
                  Value
                , Value
                , FunctionResult
                , FunctionResult
            >>
        convergence;
    };
} // namespace info
    namespace
detail_newton
{
        template <class...>
        struct
    info_data;

        template <class... Ts>
        struct
    info_data <info::none_t, Ts...>
    {
            using
        type = int const;
    };
        template <class... Ts>
        struct
    info_data <info::iterations_t, Ts...>
    {
            using
        type = info::iterations_data_t;
    };
        template <
              class Function
            , class Derivative
            , class Value
            , class Predicate
        >
        struct
    info_data <info::convergence_t, Function, Derivative, Value, Predicate>
    {
            using
        type = info::convergence_data_newton_t <
              Value
            , std::invoke_result_t <Function, Value>
            , std::invoke_result_t <Derivative, Value>
        >;
    };
        template <class... Ts>
        using
    info_data_t = info_data <Ts...>::type;
} // namespace detail_newton

*/
    struct
info_tag_t
{
        int
    code;
        friend bool
    operator == (info_tag_t, info_tag_t) = default;
};

    template <info_tag_t>
    struct
info_t
{};

    namespace
info
{
        namespace
    tag
    {
            constexpr auto
        none = info_tag_t { 0 };
            constexpr auto
        iterations = info_tag_t { 1 };
            constexpr auto
        convergence = info_tag_t { 2 };
    }
        constexpr auto
    none = info_t <tag::none> {};
        constexpr auto
    iterations = info_t <tag::iterations> {};
        constexpr auto
    convergence = info_t <tag::convergence> {};
        namespace
    data
    {
            struct
        iterations_t
        {
                int
            iteration_count;
        };
            template <
                  class Value
                , class FunctionResult
                , class DerivativeResult
            >
            struct
        convergence_newton_t
        {
                std::vector <std::tuple <
                      Value
                    , FunctionResult
                    , DerivativeResult
                >>
            convergence;
        };
            template <
                  class Value
                , class FunctionResult
            >
            struct
        convergence_zhang_t
        {
                std::vector <std::tuple <
                      Value
                    , Value
                    , FunctionResult
                    , FunctionResult
                >>
            convergence;
        };

        // Select the right data type
            template <info_tag_t, class...>
            struct
        select_newton;

            template <class... Ts>
            struct
        select_newton <tag::none, Ts...>
        {
                using
            type = int const;
        };
            template <class... Ts>
            struct
        select_newton <tag::iterations, Ts...>
        {
                using
            type = iterations_t;
        };
            template <
                  class Function
                , class Derivative
                , class Value
                , class Predicate
            >
            struct
        select_newton <
              tag::convergence
            , Function
            , Derivative
            , Value
            , Predicate
        >{
                using
            type = convergence_newton_t <
                  Value
                , std::invoke_result_t <Function, Value>
                , std::invoke_result_t <Derivative, Value>
            >;
        };
            template <info_tag_t Tag, class... Ts>
            using
        select_newton_t = select_newton <Tag, Ts...>::type;

            template <info_tag_t, class...>
            struct
        select_zhang;

            template <class... Ts>
            struct
        select_zhang <tag::none, Ts...>
        {
                using
            type = int const;
        };
            template <class... Ts>
            struct
        select_zhang <tag::iterations, Ts...>
        {
                using
            type = iterations_t;
        };
            template <
                  class Function
                , class Value
                , class Predicate
            >
            struct
        select_zhang <
              tag::convergence
            , Function
            , Value
            , Predicate
        >{
                using
            type = convergence_zhang_t <
                  Value
                , std::invoke_result_t <Function, Value>
            >;
        };
            template <info_tag_t Tag, class... Ts>
            using
        select_zhang_t = select_zhang <Tag, Ts...>::type;
    } // namespace data

} // namespace info
/*
    namespace
detail_newton
{
        template <info_tag_t, class...>
        struct
    info_data;

        template <class... Ts>
        struct
    info_data <info::tag::none, Ts...>
    {
            using
        type = int const;
    };
        template <class... Ts>
        struct
    info_data <info::tag::iterations, Ts...>
    {
            using
        type = info::data::iterations_t;
    };
        template <
              class Function
            , class Derivative
            , class Value
            , class Predicate
        >
        struct
    info_data <
          info::tag::convergence
        , Function
        , Derivative
        , Value
        , Predicate
    >{
            using
        type = info::data::convergence_newton_t <
              Value
            , std::invoke_result_t <Function, Value>
            , std::invoke_result_t <Derivative, Value>
        >;
    };
        template <info_tag_t Tag, class... Ts>
        using
    info_data_t = info_data <Tag, Ts...>::type;
} // namespace detail_newton
*/

    template <
          class Function
        , class Derivative
        , class Value
        , class Predicate
        , info_tag_t InfoTag = info::tag::none
    >
    auto
newton (
      Function&&       function
    , Derivative&&     derivative
    , Value const&     initial_guess
    , Predicate&&      converged
    , options_t const& options = options_t {}
    , [[ maybe_unused ]] info_t <InfoTag> info = info::none
){
        constexpr static auto
    need_info_iterations = InfoTag == info::tag::iterations;
        constexpr static auto
    need_info_convergence = InfoTag == info::tag::convergence;;
        constexpr static auto
    need_info = need_info_iterations || need_info_convergence;

        [[maybe_unused]]
        auto
    info_data = info::data::select_newton_t <
          InfoTag
        , Function
        , Derivative
        , Value
        , Predicate
    > {};

        Value
    current = initial_guess;
    for (int i = 0; i < options.max_iter; ++i)
    {
        // TODO: check functions evaluations?
            auto
        f = std::forward <Function> (function) (current);
        if (std::forward <Predicate> (converged) (f))
        {
            if constexpr (need_info_iterations)
            {
                info_data.iteration_count = i;
            }
            if constexpr (need_info)
            {
                return std::pair { current, info_data };
            }
            else
            {
                return current;
            }
        }
        // TODO: check functions evaluations?
            auto
        df = std::forward <Derivative> (derivative) (current);
        if (df == 0.)
        {
            throw zero_derivative_e {};
        }
        current -= f / df;
        if constexpr (need_info_convergence)
        {
            info_data.convergence.push_back ({current, f, df});
        }
    }
    throw no_convergence_e {};
}

    struct
no_single_root_between_brackets_e
{};
/*
    namespace
detail_zhang
{
        template <class...>
        struct
    info_data;

        template <class... Ts>
        struct
    info_data <info::none_t, Ts...>
    {
            using
        type = int const;
    };
        template <class... Ts>
        struct
    info_data <info::iterations_t, Ts...>
    {
            using
        type = info::iterations_data_t;
    };
        template <
              class Function
            , class Value
            , class Predicate
        >
        struct
    info_data <info::convergence_t, Function, Value, Predicate>
    {
            using
        type = info::convergence_data_zhang_t <
              Value
            , std::invoke_result_t <Function, Value>
        >;
    };
        template <class... Ts>
        using
    info_data_t = info_data <Ts...>::type;
} // namespace detail_zhang
*/
    template <
          class Function
        , class Value
        , class Predicate
        , info_tag_t InfoTag = info::tag::none
    >
    auto
zhang (
      Function&&       function
    , Value            a // bracket 1
    , Value            b // bracket 2
    , Predicate&&      converged
    , options_t const& options = options_t {}
    , [[ maybe_unused ]] info_t <InfoTag> info = info::none
){
        constexpr static auto
    need_info_iterations = InfoTag == info::tag::iterations;
        constexpr static auto
    need_info_convergence = InfoTag == info::tag::convergence;;
        constexpr static auto
    need_info = need_info_iterations || need_info_convergence;

        [[maybe_unused]]
        auto
    info_data = info::data::select_zhang_t <
          InfoTag
        , Function
        , Value
        , Predicate
    > {};

        using std::swap;
    if (b < a)
    {
            using std::swap;
        swap (a, b);
    }
    // TODO: check functions evaluations?
        auto
    fa = std::forward <Function> (function) (a);
        auto
    fb = std::forward <Function> (function) (b);
    if (fa * fb > 0)
    {
        throw no_single_root_between_brackets_e {};
    }
    for (int i = 0; i < options.max_iter; ++i)
    {
            auto
        c = (a + b) / 2;
        // here,
            auto
        fc = std::forward <Function> (function) (c);
            auto
        s = (fa != fc && fb != fc) ?
            b - fb * (b - a) / (fb - fa)
        :
              a * fb * fc / ((fa - fb) * (fa - fc)) 
            + b * fa * fc / ((fb - fa) * (fb - fc)) 
            + c * fa * fb / ((fc - fa) * (fc - fb))
        ;
        // and here
            auto
        fs = std::forward <Function> (function) (s);
        if (c > s)
        {
            swap (s , c);
            swap (fs, fc);
        }
        if (fs * fc < 0)
        {
            a  = s;
            b  = c;
            fa = fs;
            fb = fc;
        }
        else
        {
            if (fs * fb < 0)
            {
                a  = c;
                fa = fc;
            }
            else
            {
                b  = s;
                fb = fs;
            }
        }
        if constexpr (need_info_convergence)
        {
            info_data.convergence.push_back ({ a, b, fa, fb });
        }
        if (std::forward <Predicate> (converged) (a, b, fa, fb))
        {
            if constexpr (need_info_iterations)
            {
                info_data.iteration_count = i;
            }
            if constexpr (need_info)
            {
                return std::pair { (a + b) / 2,  info_data };
            }
            else
            {
                return (a + b) / 2;
            }
        }
    }
    throw no_convergence_e {};
};

} // namespace isto::root_finding


#if 0
// TODO: use std::expected.
#include <string>
#include "tl/expected.hpp"
    using tl::expected;
    using tl::make_unexpected;

    namespace
isto::misc
{
    template <class Value>
    class
result_t
    : public expected <Value, std::string>
{
public:
        using
    expected_type = expected <Value, std::string>;
    // TODO: use `explicit(something)`
        template <class Other= Value>
        /*explicit*/ constexpr
    result_t (Other&& value, int iteration_count)
        : expected_type { std::forward <Other> (value) }
        , iteration_count_m { iteration_count }
    {}
        using expected_type::expected_type;
        constexpr int
    iteration_count () const noexcept
    {
        return iteration_count_m;
    }
private:
        int
    iteration_count_m = -1;
};
    /* {{{
    template <class Value>
    class
result_t
{
public:
        template <class T = Value>
        explicit constexpr 
    result_t (T&& value, int iteration_count)
        : has_value_m { true }
        , value_m { std::forward <T> (value) }
        , iteration_count_m { iteration_count }
    {}
        explicit constexpr 
    result_t (const char* error)
        : has_value_m { false }
        , error_m { error }
    {}
        [[nodiscard]]
        constexpr bool 
    has_value () const noexcept
    {
        return has_value_m;
    }
        constexpr explicit 
    operator bool() const noexcept
    {
        return has_value_m;
    }
        constexpr int
    iteration_count () const noexcept
    {
        return iteration_count_m;
    }
        constexpr const Value&
    value () const& 
    {
        if (!has_value ())
        {
            throw std::runtime_error { error_m };
        }
        return value_m;
    }

    operator Value ()
    {
        return value ();
    }
        constexpr const char*
    error () const noexcept
    {
        return error_m;
    }
private:
        bool
    has_value_m = false;
        Value
    value_m;
        int
    iteration_count_m = -1;
        const char*
    error_m;
};
}}} */

// TODO: if `function` returns an `expected`, forward the errors.
    template <
          class Function
        , class Derivative
        , class Value
        , class Predicate
        , class Result = result_t <Value>
        , class Return = std::invoke_result_t <Derivative, Value>
    >
    Result
newton (
      Function&& function
    , Derivative&& derivative
    , Value initial_guess
    , Predicate&& converged
){
        Value
    current = initial_guess;
        constexpr int
    max_iter = 100;
    for (int i = 0; i < max_iter; ++i)
    {
        // TODO: What do if this fails?
            auto
        f = std::forward <Function> (function) (current);
        if (std::forward <Predicate> (converged) (f))
        {
            return { current, i };
        }
        // TODO: Or this?
            auto
        df = std::forward <Derivative> (derivative) (current);
        if (df == static_cast <Return> (0))
        {
            return make_unexpected ("Zero derivative encountered");
        }
        current -= f / df;
    }
    return make_unexpected ("No convergence");
}


    template <
          class Function
        , class Value
        , class Predicate
        , class Result = result_t <Value>
        , class Return = std::invoke_result_t <Function, Value>
    >
    Result
brent (
      Function&& function
    , Value a // bracket 1
    , Value b // bracket 2
    , Predicate&& converged
){
        using std::swap;
    // This is in fact Zhang's improved method.
    if (b < a)
    {
        swap (a, b);
    }
        auto
    fa = std::forward <Function> (function) (a);
        auto
    fb = std::forward <Function> (function) (b);
    if (fa * fb > static_cast <Return> (0))
    {
        return make_unexpected ("No single root between brackets");
    }
        constexpr int
    max_iter = 100;
    for (int i = 0; i < max_iter; ++i)
    {
            auto
        c = (a + b) / 2;
            auto
        fc = std::forward <Function> (function) (c);
            auto
        s = (fa != fc && fb != fc) ?
            b - fb * (b - a) / (fb - fa)
        :
              a * fb * fc / ((fa - fb) * (fa - fc)) 
            + b * fa * fc / ((fb - fa) * (fb - fc)) 
            + c * fa * fb / ((fc - fa) * (fc - fb))
        ;
            auto
        fs = std::forward <Function> (function) (s);
        if (c > s)
        {
            swap (s , c);
            swap (fs, fc);
        }
        if (fs * fc < static_cast <Return> (0))
        {
            a  = s;
            b  = c;
            fa = fs;
            fb = fc;
        }
        else
        {
            if (fs * fb < static_cast <Return> (0))
            {
                a  = c;
                fa = fc;
            }
            else
            {
                b  = s;
                fb = fs;
            }
        }
        if (std::forward <Predicate> (converged) (a, b, fa, fb))
        {
            return { (a + b) / 2, i };
        }
    }
    return make_unexpected ("No convergence");
};

} // namespace isto::misc

// vim: foldmethod=marker
#endif
