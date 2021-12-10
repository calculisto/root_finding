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
info //{{{
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
        iterations_newton_t
        {
                bool
            converged = true;
                bool
            zero_derivative = false;
                int
            iteration_count;
        };
            struct
        iterations_zhang_t
        {
                bool
            converged = true;
                bool
            no_single_root_between_bracket = false;
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
                bool
            converged = true;
                bool
            zero_derivative = false;
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
                bool
            converged = true;
                bool
            no_single_root_between_bracket = false;
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
            type = iterations_newton_t;
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
            type = iterations_zhang_t;
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

} // namespace info }}}

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
            if constexpr (need_info)
            {
                info_data.converged = false;
                info_data.zero_derivative = true;
                return std::pair { current, info_data };
            }
            else
            {
                throw zero_derivative_e {};
            }
        }
        current -= f / df;
        if constexpr (need_info_convergence)
        {
            info_data.convergence.push_back ({current, f, df});
        }
    }
    if constexpr (need_info)
    {
        info_data.converged = false;
        return std::pair { current, info_data };
    }
    else
    {
        throw no_convergence_e {};
    }
}

    struct
no_single_root_between_brackets_e
{};

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
        if constexpr (need_info)
        {
            info_data.converged = false;
            info_data.no_single_root_between_bracket = true;
            return std::pair { (a + b) / 2,  info_data };
        }
        else
        {
            throw no_single_root_between_brackets_e {};
        }
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
    if constexpr (need_info)
    {
        info_data.converged = false;
        return std::pair { (a + b) / 2, info_data };
    }
    else
    {
        throw no_convergence_e {};
    }
};

} // namespace isto::root_finding
// vim: foldmethod=marker
