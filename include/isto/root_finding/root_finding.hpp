#pragma once
#include <utility>
#include <type_traits>
#include <vector>
#include <tuple>

    namespace 
isto::root_finding
{
        namespace
    defaults
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
    } // namespace defaults

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
        base_t
        {
                bool
            converged = true;
                bool
            function_threw = false;
        };
            struct
        iterations_base_t
            : base_t
        {
                int
            iteration_count;
        };

        // Select the right data type
            template <class T, info_tag_t, class... Ts>
            struct
        select
        {};

            template <class T, class... Ts>
            struct
        select <T, tag::none, Ts...>
        {
                using
            type = int const;
        };

            template <class T, info_tag_t Tag, class... Ts>
            using
        select_t = typename select <T, Tag, Ts...>::type;
    } // namespace data

} // namespace info

// Newton method
    struct
NewtonTag
{};

// What options it can takes
    using
newton_options_t = defaults::options_t;

// What it might throw
    using
newton_no_convergence_e = defaults::no_convergence_e;

    struct
newton_zero_derivative_e
{};

// What info it might return
    namespace
info::data
{
        struct
    iterations_newton_t
        : iterations_base_t
    {
            bool
        zero_derivative = false;
            bool
        derivative_threw = false;
    };

        template <
              class Value
            , class FunctionResult
            , class DerivativeResult
        >
        struct
    convergence_newton_t
        : base_t
    {
            bool
        zero_derivative = false;
            bool
        derivative_threw = false;
            std::vector <std::tuple <
                  Value
                , FunctionResult
                , DerivativeResult
            >>
        convergence;
    };

        template <class... Ts>
        struct
    select <NewtonTag, tag::iterations, Ts...>
    {
            using
        type = iterations_newton_t;
    };

        template <
              class Function
            , class Derivative
            , class Value
        >
        struct
    select <NewtonTag, tag::convergence, Function, Derivative, Value>
    {
            using
        type = convergence_newton_t <
              Value
            , std::invoke_result_t <Function, Value>
            , std::invoke_result_t <Derivative, Value>
        >;
    };
} // namespace info::data

// The function itself
    template <
          class Function
        , class Derivative
        , class Value
        , class Predicate
        , info_tag_t InfoTag = info::tag::none
        , class FunctionReturn = std::invoke_result_t <Function, Value>
        , class DerivativeReturn = std::invoke_result_t <Derivative, Value>
    >
    requires 
           std::invocable <Function, Value> 
        && std::invocable <Derivative, Value>
    auto
newton (
      Function&&       function
    , Derivative&&     derivative
    , Value const&     initial_guess
    , Predicate&&      converged
    , newton_options_t const& options = newton_options_t {}
    ,   [[maybe_unused]] 
      info_t <InfoTag> info = info::none
){
        constexpr static auto
    need_info_iterations = InfoTag == info::tag::iterations;
        constexpr static auto
    need_info_convergence = InfoTag == info::tag::convergence;;
        constexpr static auto
    need_info = need_info_iterations || need_info_convergence;

        [[maybe_unused]]
        auto
    info_data = info::data::select_t <
          NewtonTag
        , InfoTag
        , Function
        , Derivative
        , Value
    > {};

        Value
    current = initial_guess;
    for (int i = 0; i < options.max_iter; ++i)
    {
            auto
        f = FunctionReturn {};
        try 
        {
            f = std::forward <Function> (function) (current);
        }
        catch (...)
        {
            if constexpr (need_info)
            {
                info_data.converged = false;
                info_data.function_threw = true;
                return std::pair { current, info_data };
            }
            else
            {
                throw;
            }
        }
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
            auto
        df = DerivativeReturn {};
        try
        {
            df = std::forward <Derivative> (derivative) (current);
        }
        catch (...)
        {
            if constexpr (need_info)
            {
                info_data.converged = false;
                info_data.derivative_threw = true;
                return std::pair { current, info_data };
            }
            else
            {
                throw;
            }
        }
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
                throw newton_zero_derivative_e {};
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
        throw newton_no_convergence_e {};
    }
}

// Zhang method (i.e. better Brent).
    struct
ZhangTag
{};

    using
zhang_options_t = defaults::options_t;

    using
zhang_no_convergence_e = defaults::no_convergence_e;

    struct
zhang_no_single_root_between_brackets_e
{};

    namespace
info::data
{
        struct
    iterations_zhang_t
        : iterations_base_t
    {
            bool
        no_single_root_between_bracket = false;
    };

        template <
              class Value
            , class FunctionResult
        >
        struct
    convergence_zhang_t
        : base_t
    {
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
        template <class... Ts>
        struct
    select <ZhangTag, tag::iterations, Ts...>
    {
            using
        type = iterations_zhang_t;
    };

        template <
              class Function
            , class Value
        >
        struct
    select <ZhangTag, tag::convergence, Function, Value>
    {
            using
        type = convergence_zhang_t <
              Value
            , std::invoke_result_t <Function, Value>
        >;
    };
} // namespace info::data

    template <
          class Function
        , class Value
        , class Predicate
        , info_tag_t InfoTag = info::tag::none
        , class FunctionReturn = std::invoke_result_t <Function, Value>
    >
    auto
zhang (
      Function&&  function
    , Value       a // bracket 1
    , Value       b // bracket 2
    , Predicate&& converged
    , zhang_options_t const& options = zhang_options_t {}
    ,   [[maybe_unused]] 
      info_t <InfoTag> info = info::none
){
        constexpr static auto
    need_info_iterations = InfoTag == info::tag::iterations;
        constexpr static auto
    need_info_convergence = InfoTag == info::tag::convergence;;
        constexpr static auto
    need_info = need_info_iterations || need_info_convergence;

        [[maybe_unused]]
        auto
    info_data = info::data::select_t <
          ZhangTag
        , InfoTag
        , Function
        , Value
    > {};

        using std::swap;
    if (b < a)
    {
            using std::swap;
        swap (a, b);
    }
        FunctionReturn
      fa
    , fb
    ;
    try
    {
        fa = std::forward <Function> (function) (a);
        fb = std::forward <Function> (function) (b);
    }
    catch (...)
    {
        if constexpr (need_info)
        {
            info_data.converged = false;
            info_data.function_threw = true;
            return std::pair { (a + b) / 2, info_data };
        }
        else
        {
            throw;
        }
    }
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
            throw zhang_no_single_root_between_brackets_e {};
        }
    }
    for (int i = 0; i < options.max_iter; ++i)
    {
            auto
        c = (a + b) / 2;
            auto
        fc = FunctionReturn {};
        try
        {
            fc = std::forward <Function> (function) (c);
        }
        catch (...)
        {
            if constexpr (need_info)
            {
                info_data.converged = false;
                info_data.function_threw = true;
                return std::pair { (a + b) / 2, info_data };
            }
            else
            {
                throw;
            }
        }
            auto
        s = (fa != fc && fb != fc) ?
            b - fb * (b - a) / (fb - fa)
        :
              a * fb * fc / ((fa - fb) * (fa - fc)) 
            + b * fa * fc / ((fb - fa) * (fb - fc)) 
            + c * fa * fb / ((fc - fa) * (fc - fb))
        ;
            auto
        fs = FunctionReturn {};
        try
        {
            fs = std::forward <Function> (function) (s);
        }
        catch (...)
        {
            if constexpr (need_info)
            {
                info_data.converged = false;
                info_data.function_threw = true;
                return std::pair { (a + b) / 2, info_data };
            }
            else
            {
                throw;
            }
        }
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
        throw zhang_no_convergence_e {};
    }
};

} // namespace isto::root_finding
