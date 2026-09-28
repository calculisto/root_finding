#pragma once
#include <cmath>
#include <utility>
#include <type_traits>
#include <vector>
#include <tuple>
#include <numbers>
#include <numbers>
#include <array>
#include <ranges>
#include <valarray>
#include <functional>
// For automatic differenciation
#include <calculisto/auto_diff/dual.hpp>
    using calculisto::auto_diff::dual_t;
// For multidimensions
#include <eigen3/Eigen/Dense>

    namespace
calculisto::root_finding
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
        base_iterations_t
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

//------------------------------------------------------------------------------
// Newton method
    struct
NewtonTag
{};

// How it stops, by default
    template <
          class Value
        , class FunctionResult
    >
    bool
newton_default_converged (
      Value            const& current
    , Value            const& past
    , FunctionResult   const& result
){
        using std::fabs;
    return
           fabs ((past - current) / current)
           <
           std::numeric_limits <Value>::epsilon ()
        || result == 0;
    ;
}

// Build an alternative simple convergence predicate
    template <class Value>
    auto
make_newton_simple_converged (Value const& tolerance)
{
    return [=](auto const& current, auto const& past, auto const& result)
    {
            using std::fabs;
        return fabs ((past - current) / current) < tolerance || result == 0;
    };
}

// What options it can take
    template <
          class Value
        , class FunctionResult
        , class DerivativeResult
    >
    struct
newton_options_t
{
        int
    max_iter = 100;
        std::function <bool (
          Value const&
        , Value const&
        , FunctionResult const&
    )>
    converged = &newton_default_converged <Value, FunctionResult>;
};

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
    newton_iterations_t
        : base_iterations_t
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
        type = newton_iterations_t;
    };

        template <
              class FunctionResult
            , class DerivativeResult
            , class Value
        >
        struct
    select <NewtonTag, tag::convergence, FunctionResult, DerivativeResult, Value>
    {
            using
        type = convergence_newton_t <
              Value
            , FunctionResult
            , DerivativeResult
        >;
    };
} // namespace info::data

// The function itself
    template <
          class Function
        , class Derivative
        , class Value
        , info_tag_t InfoTag = info::tag::none
        , class FunctionResult = std::invoke_result_t <Function, Value>
        , class DerivativeResult = std::invoke_result_t <Derivative, Value>
    >
    requires
           std::invocable <Function, Value>
        && std::invocable <Derivative, Value>
    auto
newton (
      Function&&       function
    , Derivative&&     derivative
    , Value const&     initial_guess
    , newton_options_t <Value, FunctionResult, DerivativeResult> const& options = {}
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
        , FunctionResult
        , DerivativeResult
        , Value
    > {};

        Value
    past;
        Value
    current = initial_guess;
    for (int i = 0; i < options.max_iter; ++i)
    {
            auto
        f = FunctionResult {};
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
            auto
        df = DerivativeResult {};
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
        past = current;
        current -= f / df;
        if constexpr (need_info_convergence)
        {
            info_data.convergence.push_back ({current, f, df});
        }
        if (options.converged (current, past, f))
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

//------------------------------------------------------------------------------
// Newton method with automatic differenciation
    template <
          class Function
        , class Value
        , info_tag_t InfoTag = info::tag::none
        , class FunctionResult = std::invoke_result_t <Function, Value>
        , class DerivativeResult = FunctionResult
        , class Derivative = Function
    >
    requires std::invocable <Function, Value>
    auto
newton (
      Function&&       function
    , Value const&     initial_guess
    , newton_options_t <Value, FunctionResult, DerivativeResult> const& options = {}
    ,   [[maybe_unused]]
      info_t <InfoTag> info = info::none
){
    static_assert (
          std::same_as <FunctionResult, Value>
        , "automatic differentiation requires this"
    );
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
        , FunctionResult
        , DerivativeResult
        , Value
    > {};

        using
    dual_type = dual_t <1, Value>;

        auto 
    past = dual_type {};
        auto
    current = dual_type { initial_guess, 0 };
    for (int i = 0; i < options.max_iter; ++i)
    {
            auto
        f_df = dual_type {};
        try
        {
            f_df = std::forward <Function> (function) (current);
        }
        catch (...)
        {
            if constexpr (need_info)
            {
                info_data.converged = false;
                info_data.function_threw = true;
                return std::pair { current.value, info_data };
            }
            else
            {
                throw;
            }
        }
            const auto
        f = f_df.value;
            const auto
        df = f_df.differentials[0];
        if (df == 0.)
        {
            if constexpr (need_info)
            {
                info_data.converged = false;
                info_data.zero_derivative = true;
                return std::pair { current.value, info_data };
            }
            else
            {
                throw newton_zero_derivative_e {};
            }
        }
        past = current;
        current -= f / df;
        // FIXME: needed?
        current.differentials[0] = 1;

        if constexpr (need_info_convergence)
        {
            info_data.convergence.push_back ({current.value, f, df});
        }
        if (options.converged (current.value, past.value, f))
        {
            if constexpr (need_info_iterations)
            {
                info_data.iteration_count = i;
            }
            if constexpr (need_info)
            {
                return std::pair { current.value, info_data };
            }
            else
            {
                return current.value;
            }
        }
    }
    if constexpr (need_info)
    {
        info_data.converged = false;
        return std::pair { current.value, info_data };
    }
    else
    {
        throw newton_no_convergence_e {};
    }
}

//------------------------------------------------------------------------------
// Zhang method (i.e. better Brent).
    struct
ZhangTag
{};

    template <class Value, class FunctionResult>
    bool
zhang_default_converged (
      Value const& a
    , Value const& b
    , FunctionResult const& fa
    , FunctionResult const& fb
){
        using std::fabs;
    return
           fa == 0
        || fb == 0
        || fabs (b - a) < std::numeric_limits <Value>::epsilon ()
    ;
}

    template <class Value>
    auto
make_zhang_simple_converged (Value const& tolerance)
{
    return [=] <class FunctionResult> (
          Value const& a
        , Value const& b
        , FunctionResult const& fa
        , FunctionResult const& fb
    ) {
            using std::fabs;
    return
           fa == 0
        || fb == 0
        || fabs (b - a) < tolerance
    ;
    };
}

    template <class Value, class FunctionResult>
    struct
zhang_options_t
{
        int
    max_iter = 100;
        std::function <bool (
              Value const&
            , Value const&
            , FunctionResult const&
            , FunctionResult const&
        )>
    converged = zhang_default_converged <Value, FunctionResult>;
};

    using
zhang_no_convergence_e = defaults::no_convergence_e;

    struct
zhang_no_single_root_between_brackets_e
{};

    namespace
info::data
{
        struct
    zhang_iterations_t
        : base_iterations_t
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
        type = zhang_iterations_t;
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
        , info_tag_t InfoTag = info::tag::none
        , class FunctionResult = std::invoke_result_t <Function, Value>
    >
    auto
zhang (
      Function&&  function
    , Value       a // bracket 1
    , Value       b // bracket 2
    , zhang_options_t<Value, FunctionResult> const& options = {}
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
        FunctionResult
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
        fc = FunctionResult {};
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
        fs = FunctionResult {};
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
        if (options.converged (a, b, fa, fb))
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

//------------------------------------------------------------------------------
// Bracket an extremum
    struct
BracketExtremaTag
{};
// https://stackoverflow.com/a/58876657/1622545
    struct
bracket_minimum_options_t
{
        int
    max_iter = 100;
        double
    gold = std::numbers::phi;
};

    struct
bracket_minimum_no_convergence_e
{};

    namespace
info::data
{
        struct
    bracket_minimum_iterations_t
        : base_iterations_t
    {};

        template <
              class Value
            , class FunctionResult
        >
        struct
    bracket_minimum_convergence_t
        : base_t
    {
            std::vector <std::array <std::pair <
                  Value
                , FunctionResult
            >, 3>>
        convergence;
    };
        template <class... Ts>
        struct
    select <BracketExtremaTag, tag::iterations, Ts...>
    {
            using
        type = bracket_minimum_iterations_t;
    };

        template <
              class Function
            , class Value
        >
        struct
    select <BracketExtremaTag, tag::convergence, Function, Value>
    {
            using
        type = bracket_minimum_convergence_t <
              Value
            , std::invoke_result_t <Function, Value>
        >;
    };
} // namespace info::data

    template <
          class Function
        , class Value
        , info_tag_t InfoTag = info::tag::none
        , class FunctionResult = std::invoke_result_t <Function, Value>
        , class Return = std::tuple <Value, Value, FunctionResult, FunctionResult>
    >
    requires std::invocable <Function, Value>
    auto
bracket_minimum (
      Function&& function
    , Value      a
    , Value      b
    , bracket_minimum_options_t const& options = {}
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
          BracketExtremaTag
        , InfoTag
        , Function
        , Value
    > {};

        FunctionResult
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
            return std::pair { Return {}, info_data };
        }
        else
        {
            throw;
        }
    }
        auto
    h = b - a;
    if (fa < fb)
    {
        h = -h;
            using std::swap;
        swap (a, b);
        swap (fa, fb);
    }
    for (auto i = 0; i < options.max_iter; ++i)
    {
            const auto
        c = b + h;
            FunctionResult
        fc;
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
                return std::pair { Return { a, c, fa, fc }, info_data };
            }
            else
            {
                throw;
            }
        }
        if constexpr (need_info_convergence)
        {
            info_data.convergence.push_back ({{{ a, fa }, { b, fb }, { c, fc }}});
        }
        if (fc > fb)
        {
            if constexpr (need_info_iterations)
            {
                info_data.iteration_count = i;
            }
            if constexpr (need_info)
            {
                return std::pair { Return { a, c, fa, fc }, info_data };
            }
            else
            {
                return Return { a, c, fa, fc };
            }
        }
        a  = b;
        fa = fb;
        b  = c;
        fb = fc;
        h *= options.gold;
    }
    if constexpr (need_info)
    {
        info_data.converged = false;
        return std::pair { Return {}, info_data };
    }
    else
    {
        throw bracket_minimum_no_convergence_e {};
    }
}

//------------------------------------------------------------------------------
// Golden section search
    struct
GoldenSectionTag
{};

    template <class Value>
    struct
golden_section_options_t
{
        Value
    tolerance = std::numeric_limits <Value>::epsilon ();
        bracket_minimum_options_t
    bracket_minimum_options = bracket_minimum_options_t {};
};

    using
golden_section_no_convergence_e = defaults::no_convergence_e;

    namespace
info::data
{
        struct
    golden_section_iterations_t
        : base_iterations_t
    {
            bracket_minimum_iterations_t
        bracket_minimum_info;
    };

        template <
              class Value
            , class FunctionResult
        >
        struct
    golden_section_convergence_t
        : base_t
    {
            bracket_minimum_convergence_t <Value, FunctionResult>
        bracket_minimum_info;
            std::vector <std::tuple <
                  std::pair <Value, FunctionResult>
                , std::pair <Value, FunctionResult>
                , std::pair <Value, FunctionResult>
                , std::pair <Value, FunctionResult>
            >>
        convergence;
    };

        template <class... Ts>
        struct
    select <GoldenSectionTag, tag::iterations, Ts...>
    {
            using
        type = golden_section_iterations_t;
    };

        template <
              class Function
            , class Value
        >
        struct
    select <GoldenSectionTag, tag::convergence, Function, Value>
    {
            using
        type = golden_section_convergence_t <
              Value
            , std::invoke_result_t <Function, Value>
        >;
    };
} // namespace info::data

    template <
          class Function
        , class Value
        , info_tag_t InfoTag = info::tag::none
        , class FunctionResult = std::invoke_result_t <Function, Value>
    >
    requires std::invocable <Function, Value>
    auto
golden_section (
      Function&&        function
    , Value             a
    , Value             b
    , golden_section_options_t <Value> const& options = {}
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
          GoldenSectionTag
        , InfoTag
        , Function
        , Value
    > {};

        FunctionResult
      fa
    , fb
    ;
    if constexpr (need_info)
    {
            auto
        [r, in] = bracket_minimum (
              std::forward <Function> (function)
            , a
            , b
            , options.bracket_minimum_options
            , info
        );
            std::tie
        (a, b, fa, fb) = std::move (r);
        info_data.bracket_minimum_info = std::move (in);
        if (!info_data.bracket_minimum_info.converged)
        {
            info_data.converged = false;
            return std::pair { 0., info_data};
        }
    }
    else
    {
            std::tie
        (a, b, fa, fb) = bracket_minimum (
              std::forward <Function> (function)
            , a
            , b
            , options.bracket_minimum_options
        );
    }
    if (a > b)
    {
            using std::swap;
        swap (a, b);
    }
        using std::numbers::phi;
        auto
    h = b - a;
        auto
    c = a + h / phi / phi;
        auto
    d = a + h / phi;
        FunctionResult
      fc
    , fd
    ;
    try
    {
        fc = std::forward <Function> (function) (c);
        fd = std::forward <Function> (function) (d);
    }
    catch (...)
    {
        if constexpr (need_info)
        {
            info_data.converged = false;
            info_data.function_threw = true;
            return std::pair { 0., info_data };
        }
        else
        {
            throw;
        }
    }
    // FIXME: check for overflow
        using std::ceil;
        const auto
    n = static_cast <int> (std::round (
        ceil (log (options.tolerance / h) / log (1. / phi))
    ));
    if constexpr (need_info_convergence)
    {
        info_data.convergence.push_back ({
              { a, fa }
            , { c, fc }
            , { d, fd }
            , { b, fb }
        });
    }
    for (auto i = 0; i < n; ++i)
    {
        if (fc < fd)
        {
            b  = d;
            d  = c;
            fd = fc;
            h  = h / phi;
            c  = a + h / phi / phi;
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
                    return std::pair { 0., info_data };
                }
                else
                {
                    throw;
                }
            }
        }
        else
        {
            a  = c;
            c  = d;
            fc = fd;
            h  = h / phi;
            d  = a + h / phi;
            try
            {
                fd = std::forward <Function> (function) (d);
            }
            catch (...)
            {
                if constexpr (need_info)
                {
                    info_data.converged = false;
                    info_data.function_threw = true;
                    return std::pair { 0., info_data };
                }
                else
                {
                    throw;
                }
            }
        }
        if constexpr (need_info_convergence)
        {
            info_data.convergence.push_back ({ { a, fa }, { c, fc }, { d, fd }, { b, fb } });
        }
    }
    if constexpr (need_info_iterations)
    {
        info_data.iteration_count = n;
    }
    if (fc < fd)
    {
        if constexpr (need_info)
        {
            return std::pair { (a + d) / 2., info_data };
        }
        else
        {
            return (a + d) / 2.;
        }
    }
    if constexpr (need_info)
    {
        return std::pair { (c + b) / 2., info_data };
    }
    else
    {
        return (c + b) / 2.;
    }
}

//------------------------------------------------------------------------------
// Powell
    struct
PowellTag
{};

    template <class Value, class FunctionResult>
    struct
powell_options_t
{
        int
    max_iter = 100;
        FunctionResult
    tolerance = std::numeric_limits <FunctionResult>::epsilon ();
        golden_section_options_t <Value>
    golden_section_options = {};
};

    using
powell_no_convergence_e = defaults::no_convergence_e;

    namespace
info::data
{
        struct
    powell_iterations_t
        : base_iterations_t
    {
            golden_section_iterations_t
        golden_section_info;
    };

        template <
              class Function
            , class Value
            , class FunctionResult = std::invoke_result_t <Function, std::valarray <Value>>
        >
        struct
    powell_convergence_t
        : base_t
    {
            std::vector <golden_section_convergence_t <Value, FunctionResult>>
        golden_section_info;
            std::vector <std::tuple <int, int, FunctionResult, std::valarray <Value>>>
        convergence;
    };

        template <class... Ts>
        struct
    select <PowellTag, tag::iterations, Ts...>
    {
            using
        type = powell_iterations_t;
    };

        template <
              class Function
            , class Value
        >
        struct
    select <PowellTag, tag::convergence, Function, Value>
    {
            using
        type = powell_convergence_t <
              Function
            , Value
        >;
    };
} // namespace info::data

/* Don't know how to properly mix lambda capture and perfect forwarding, so we
 * copy the functor. See
 * https://stackoverflow.com/q/54418941/1622545
 */
    template <
          class Function
        , class Value
        , info_tag_t InfoTag = info::tag::none
        , class FunctionResult = std::invoke_result_t <Function, std::valarray <Value>>
    >
    requires std::invocable <Function, std::valarray <Value>>
    auto
powell (
      Function                function
    , std::valarray <Value>&& init
    , powell_options_t <Value, FunctionResult> const& options = {}
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
          PowellTag
        , InfoTag
        , Function
        , Value
    > {};

        const auto
    n = init.size ();
        auto
    xi = std::valarray (std::valarray (0., n), n);
    for (auto i = 0u; i < n; ++i)
    {
        xi[i][i] = 1.;
    }
        auto
    p = std::valarray (init);
        FunctionResult
    f;
    if constexpr (need_info)
    {
        try
        {
            f = function (p);
        }
        catch (...)
        {
            info_data.converged = false;
            info_data.function_threw = true;
            return std::pair { p, info_data };
        }
    }
    for (int j = 1; j <= options.max_iter; ++j)
    {
            auto
        p0 = p;
            auto
        f0 = f;
            auto
        delta = std::numeric_limits <Value>::min ();
            auto
        max_index = 0;
        for (auto i = 0u; i < n; ++i)
        {
                auto
            f_ = [&](auto lambda){ return function (p + lambda * xi[i]); };
                Value
            lambda;
            if constexpr (need_info)
            {
                    auto
                [ lam, inf ] = golden_section (f_, 0., 0.1, options.golden_section_options, info);
                info_data.golden_section_info.push_back (std::move (inf));
                if (!info_data.golden_section_info.back ().converged)
                {
                    info_data.converged = false;
                    return std::pair { p, info_data };
                }
                lambda = std::move (lam);
            }
            else
            {
                lambda = golden_section (f_, 0., 0.1, options.golden_section_options);
            }
            p += lambda * xi[i];
                const auto
            f_prev = f;
            if constexpr (need_info)
            {
                try
                {
                    f = function (p);
                }
                catch (...)
                {
                    info_data.converged = false;
                    info_data.function_threw = true;
                    return std::pair { p, info_data };
                }
            }
            if (f_prev - f > delta)
            {
                delta = f_prev - f;
                max_index = i;
            }
            if constexpr (need_info_convergence)
            {
                info_data.convergence.push_back ({j, i, f, p});
            }
        }
            const auto
        f3 = function (2. * p - p0);
        if (f3 < f0 && (f0 - 2. * f + f3)
            * pow (f0 - f - delta, 2.) < 0.5 * pow (f0 - f3, 2.)
        ){
                const auto
            xi_ = p - p0;
                auto
            f_ = [&](auto lambda){ return function (p + lambda * xi_); };
                Value
            lambda;
            if constexpr (need_info)
            {
                    auto
                [ lam, inf ] = golden_section (f_, 0., 0.1, options.golden_section_options, info);
                info_data.golden_section_info.push_back (std::move (inf));
                if (!info_data.golden_section_info.back ().converged)
                {
                    info_data.converged = false;
                    return std::pair { p, info_data };
                }
                lambda = std::move (lam);
            }
            else
            {
                lambda = golden_section (f_, 0., 0.1, options.golden_section_options);
            }
            xi[max_index] = xi_;
            p += lambda * xi_;
            if constexpr (need_info)
            {
                try
                {
                    f = function (p);
                }
                catch (...)
                {
                    info_data.converged = false;
                    info_data.function_threw = true;
                    return std::pair { p, info_data };
                }
            }
            if constexpr (need_info_convergence)
            {
                info_data.convergence.push_back ({j, n, f, p});
            }
        }
        if (fabs (f - f0) < options.tolerance)
        {
            if constexpr (need_info)
            {
                return std::pair { p, info_data };
            }
            else
            {
                return p;
            }
        }
    }
    info_data.converged = false;
    return std::pair { p, info_data };
}

//------------------------------------------------------------------------------
// Multidimensional Newton method with automatic differentiation
    struct
MultidimensionalNewtonTag
{};

    namespace
detail
{
    template <class T, class U>
    struct
same_container_as;

    template <class T, class U>
    struct
same_container_as <std::vector <T>, U>
{
        using
    type = std::vector <U>;
};
    template <class T, std::size_t N, class U>
    struct
same_container_as <std::array <T, N>, U>
{
        using
    type = std::array <U, N>;
};
    template <class T, std::size_t N, class U>
    struct
same_container_as <T[N], U>
{
        using
    type = U[N];
};
    template <class T, int Rows, int Cols, class U>
    struct
same_container_as <Eigen::Matrix <T, Rows, Cols>, U>
{
        using
    type = Eigen::Matrix <U, Rows, Cols>;
};
} // namespace detail

    template <class T, class U>
    using
same_container_as_t = detail::same_container_as <std::remove_cvref_t <T>, U>::type;

    namespace
detail::tests
{
    static_assert (std::same_as <same_container_as_t <std::vector <double>, int>, std::vector <int>>);
    static_assert (std::same_as <same_container_as_t <std::array <double, 3>, int>, std::array <int, 3>>);
    static_assert (std::same_as <same_container_as_t <double[3], int>, int[3]>);
    static_assert (std::same_as <same_container_as_t <Eigen::Matrix <double, 3, 3>, int>, Eigen::Matrix <int, 3, 3>>);
} // namespace detail::tests

    namespace
detail
{
        template <
              std::ranges::random_access_range Range
            , class DualValue = std::ranges::range_value_t <Range>
            , class Value = DualValue::value_type
        >
        auto
    values (Range const& x)
    {
            auto
        r = same_container_as_t <Range, Value> {};
        for (auto const& [index, dual]: std::views::enumerate (x))
        {
            r[index] = dual.value;
        }
        return r;
    }
} // namespace detail

    namespace
info::data
{
        template <class... Ts>
        struct
    select <MultidimensionalNewtonTag, tag::iterations, Ts...>
    {
            using
        type = newton_iterations_t;
    };

        template <
              class FunctionResult
            , class DerivativeResult
            , class Range
        >
        struct
    select <MultidimensionalNewtonTag, tag::convergence, FunctionResult, DerivativeResult, Range>
    {
            using
        type = convergence_newton_t <
              Range
            , FunctionResult
            , DerivativeResult
        >;
    };
} // namespace info::data

    template <std::size_t Size, class Value>
    bool
default_multidimensional_newton_convergence_test (
      Eigen::Matrix <Value, Size, 1> const& current
    , Eigen::Matrix <Value, Size, 1> const& past
    , Eigen::Matrix <Value, Size, 1> const& result
){
    return
            (past - current).norm () / current.norm ()
            <
            std::numeric_limits <Value>::epsilon ()
        || result.squaredNorm () == 0
    ;
};
    template <
          std::size_t Size
        , class Range
        , class FunctionResult
        , class DerivativeResult
        , class Value = std::ranges::range_value_t <Range>
    >
    struct
multidimensional_newton_options_t
{
        int
    max_iter = 100;
        std::function <bool (
          Range const&
        , Range const&
        , FunctionResult const&
    )>
    converged = &default_multidimensional_newton_convergence_test <Size, Value>;
};

    template <
          std::size_t Size
        , class Function
        , std::ranges::random_access_range Range
        , info_tag_t InfoTag = info::tag::none
        , class FunctionResult = std::invoke_result_t <Function, Range>
        , class Value = std::ranges::range_value_t <Range>
    >
    requires (
           std::invocable <Function, Range>
        && Size < std::numeric_limits <int>::max ()
    )
    auto
newton (
      Function&& function
    , Range&&    initial_guess
    , multidimensional_newton_options_t <
          Size
        , std::remove_cvref_t <Range>
        , FunctionResult
        , FunctionResult
      > const& options = {}
    ,   [[maybe_unused]]
      info_t <InfoTag> info = info::none
){
    if (std::ranges::size (initial_guess) != Size)
    {
        // FIXME: throw something
    }

        constexpr static auto
    need_info_iterations = InfoTag == info::tag::iterations;
        constexpr static auto
    need_info_convergence = InfoTag == info::tag::convergence;;
        constexpr static auto
    need_info = need_info_iterations || need_info_convergence;

        [[maybe_unused]]
        auto
    info_data = info::data::select_t <
          MultidimensionalNewtonTag
        , InfoTag
        , std::remove_cvref_t <Range>
        , Eigen::Matrix <Value, Size, Size>
        , std::remove_cvref_t <Range>
    > {};

        using
    dual_type = dual_t <Size, Value>;
        using
    container_of_dual = same_container_as_t <Range, dual_type>;

        auto
    past = container_of_dual {};
        auto
    current = container_of_dual {};

    for (std::size_t i = 0; i < Size; ++i)
    {
        current[i] = dual_type { initial_guess[i], i };
    }
        auto
    jacobian = Eigen::Matrix <Value, Size, Size> {};
        auto
    b = Eigen::Matrix <Value, Size, 1> {};
        auto
    f_df = container_of_dual {};
    for (int iter = 0; iter < options.max_iter; ++iter)
    {
        try
        {
            f_df = std::forward <Function> (function) (current);
        }
        catch (...)
        {
            if constexpr (need_info)
            {
                info_data.converged = false;
                info_data.function_threw = true;
                return std::pair { detail::values (current), info_data };
            }
            else
            {
                throw;
            }
        }
        for (std::size_t i = 0; i < Size; ++i)
        {
            b[i] = -f_df[i].value;
            for (std::size_t j = 0; j < Size; ++j)
            {
                jacobian(i, j) = f_df[i].differentials[j];
            }
        }
        if constexpr (need_info_convergence)
        {
            info_data.convergence.push_back ({ detail::values (current), detail::values (f_df), jacobian});
        }
        past = current;
            const auto
        delta = jacobian.fullPivLu ().solve(b);
        for (std::size_t i = 0; i < Size; ++i)
        {
            current[i].value += delta[i];
        }
            auto const
        current_value = detail::values (current);
            auto const
        current_f = detail::values (f_df);
        if (options.converged (current_value, detail::values (past), current_f))
        {
            if constexpr (need_info_iterations)
            {
                info_data.iteration_count = iter;
            }
            if constexpr (need_info)
            {
                return std::pair { current_value, info_data };
            }
            else
            {
                return current_value;
            }
        }
    }
    if constexpr (need_info)
    {
        info_data.converged = false;
        return std::pair { detail::values (current), info_data };
    }
    else
    {
        throw newton_no_convergence_e {};
    }
}
//------------------------------------------------------------------------------
// Multidimensional Newton method with a jacobian
    template <
          std::size_t Size
        , class Function
        , class Jacobian
        , std::ranges::random_access_range Range
        , info_tag_t InfoTag = info::tag::none
        , class FunctionResult = std::invoke_result_t <Function, Range>
        , class JacobianResult = std::invoke_result_t <Jacobian, Range>
        , class Value = std::ranges::range_value_t <Range>
    >
    requires (
           std::invocable <Function, Range>
        && std::invocable <Jacobian, Range>
        && Size < std::numeric_limits <int>::max ()
    )
    auto
newton (
      Function&& function
    , Jacobian&& jacobian
    , Range&&    initial_guess
    , multidimensional_newton_options_t <
          Size
        , std::remove_cvref_t <Range>
        , FunctionResult
        , FunctionResult
      > const& options = {}
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
        , FunctionResult
        , JacobianResult
        , Range
    > {};

        using 
    ActualRange = std::remove_cvref_t <Range>;
        ActualRange
    past;
        ActualRange
    current = initial_guess;
    for (int i = 0; i < options.max_iter; ++i)
    {
            auto
        f = FunctionResult {};
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
            auto
        j = JacobianResult {};
        try
        {
            j = std::forward <Jacobian> (jacobian) (current);
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
        /*
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
        */
        past = current;
            const auto
        delta = j.colPivHouseholderQr().solve(-f);
        current += delta;
        if constexpr (need_info_convergence)
        {
            info_data.convergence.push_back ({current, f, j});
        }
        if (options.converged (current, past, f))
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
} // namespace calculisto::root_finding
