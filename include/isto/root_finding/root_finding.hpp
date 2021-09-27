#include <utility>
#include <type_traits>
#pragma once
    namespace isto::root_finding
{
    struct
info_none_t
{};

    struct
info_iterations_t
{
        int
    iteration_count;
};

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

    template <
          class Function
        , class Derivative
        , class Value
        , class Predicate
        , class Info = info_none_t
        , class Return = std::invoke_result_t <Function, Value>
    >
    auto
newton (
      Function&&       function
    , Derivative&&     derivative
    , Value&&          initial_guess
    , Predicate&&      converged
    , options_t const& options = options_t {}
    , Info&&           info = info_none_t {}
){
    // TODO: check functions evaluations?
        Value
    current = initial_guess;
    for (int i = 0; i < options.max_iter; ++i)
    {
        // Here,
            auto
        f = std::forward <Function> (function) (current);
        if (std::forward <Predicate> (converged) (f))
        {
            if constexpr (std::is_same_v <Info, info_iterations_t&>)
            {
                info.iteration_count = i;
            }
            return current;
        }
        // and here.
            auto
        df = std::forward <Derivative> (derivative) (current);
        if (df == 0.)
        {
            // TODO: add some info?
            throw zero_derivative_e {};
        }
        current -= f / df;
    }
    throw no_convergence_e {};
}

    struct
no_single_root_between_brackets_e
{};

    template <
          class Function
        , class Value
        , class Predicate
        , class Info = info_none_t
        , class Return = std::invoke_result_t <Function, Value>
    >
    auto
zhang (
      Function&&       function
    , Value            a // bracket 1
    , Value            b // bracket 2
    , Predicate&&      converged
    , options_t const& options = options_t {}
    , Info&&           info = info_none_t {}
){
    // TODO: check functions evaluations?
        using std::swap;
    if (b < a)
    {
            using std::swap;
        swap (a, b);
    }
    // Here (x2),
        auto
    fa = std::forward <Function> (function) (a);
        auto
    fb = std::forward <Function> (function) (b);
    if (fa * fb > static_cast <Return> (0))
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
            if constexpr (std::is_same_v <Info, info_iterations_t&>)
            {
                info.iteration_count = i;
            }
            return (a + b) / 2;
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
