#include <doctest/doctest.h>
#include <fmt/ranges.h>
#include "../include/calculisto/root_finding/root_finding.hpp"
    using namespace calculisto::root_finding;
#include <cmath>
    using std::cos, std::sin, std::pow;

    namespace
doctest
{
    template <std::ranges::range Range>
    struct
StringMaker <Range>
{
        static String
    convert(Range const& r)
    {
        return fmt::format ("{}", r).c_str ();
    }
};
}

    auto
f1 = [](auto x){ return cos (x) - pow (x, 3.0); };
    auto
df1 = [](double x){ return -sin (x) - 3 * pow (x, 2.0); };
    auto const
target1 = 0.8654740331016144466206859011862287477929;

    auto
f2 = [](auto x){ return cos (x) - x; };
    auto
df2 = [](double x){ return -sin (x) - 1; };
    auto const
target2 = 0.7390851332151606416553120876738734040134 ;

    auto
f3 = [](auto x){ throw int {}; return static_cast <decltype (x)> (1); };
    auto
df3 = [](double) { return 0; };

    auto
f4 = [](double x){ return static_cast <decltype (x)> (1); };
    auto
df4 = [](double) { throw int {}; return 1.; };

    bool
cvg1 (double curr, double old, double f)
{
    return f == 0 || fabs (curr - old) < 1e-12;
}
// -----------------------------------------------------------------------------
TEST_CASE("Newton")
{
    SUBCASE("newton")
    {
            auto
        r = newton (f1, df1, 1.0);
        CHECK(r == doctest::Approx { target1 });

            auto
        s = newton (f2, df2, 1.0);
        CHECK(s == doctest::Approx { target2 });
    }
    SUBCASE("newton, with custom convergence criterion")
    {
            auto
        r = newton (f1, df1, 1.0, { .converged = cvg1 });
        CHECK(r == doctest::Approx { target1 });
    }
    SUBCASE("newton, throws if zero derivative")
    {
        CHECK_THROWS_AS(
              newton (f1, [](auto){ return 0.0; }, 1.0)
            , newton_zero_derivative_e
        );
    }
    SUBCASE("newton, user function throws")
    {
        CHECK_THROWS_AS(
              newton (f3, df3, 1.0)
            , int
        );
        CHECK_THROWS_AS(
              newton (f4, df4, 1.0)
            , int
        );
    }
    SUBCASE("newton, with options")
    {
        CHECK_THROWS_AS(
              newton (f1, df1, 1.0, { .max_iter = 1 })
            , newton_no_convergence_e
        );
    }
    SUBCASE("newton, with info (iteration count)")
    {
            auto const
        [ result, info ] = newton (f1, df1, 1.0, { /*default options*/ }, info::iterations);
        CHECK(info.iteration_count > 1);
    }
    SUBCASE("newton, with info (convergence)")
    {
            auto const
        [ result, info ] = newton (f1, df1, 1.0, { /*default options*/ }, info::convergence);
        CHECK(info.convergence.size () > 1);
        for (auto&& [v, f ,df]: info.convergence)
        {
            MESSAGE (v, ", ", f, ", ", df);
        }
    }
    SUBCASE("newton, with info convergence does not throw no_convergence_e!")
    {
            auto const
        [ result, info ] = newton (f1, df1, 1.0, { .max_iter = 3 }, info::convergence);
        CHECK(info.convergence.size () == 3);
        CHECK(info.converged == false);
        for (auto&& [v, f ,df]: info.convergence)
        {
            MESSAGE (v, ", ", f, ", ", df);
        }
    }
    SUBCASE("In fact, with any info, never throw")
    {
            auto const
        [ result, info ] = newton (f1, [](auto){ return 0.0; }, 1.0, {}, info::iterations);
        CHECK(!info.converged);
        CHECK(info.zero_derivative);
    }
    SUBCASE("newton, user function throws, with info")
    {
            auto const
        [ result, info ] = newton (f3, df3, 1.0, {}, info::iterations);
        CHECK(!info.converged);
        CHECK(info.function_threw);
    }
    SUBCASE("newton, user derivative throws, with info")
    {
            auto const
        [ result, info ] = newton (f4, df4, 1.0, {}, info::iterations);
        CHECK(!info.converged);
        CHECK(info.derivative_threw);
    }
    SUBCASE("Use alternative convergence predicate")
    {
            auto
        r = newton (f1, df1, 1.0, { .converged = make_newton_simple_converged (1e-8) });
        CHECK(r == doctest::Approx { target1 });
    }
}
// -----------------------------------------------------------------------------
TEST_CASE("Newton with automatic differenciation")
{
    SUBCASE("newton")
    {
            auto
        r = newton (f1, 1.0);
        CHECK(r == doctest::Approx { target1 });

            auto
        s = newton (f2, 1.0);
        CHECK(s == doctest::Approx { target2 });
    }
    SUBCASE("newton, with custom convergence criterion")
    {
            auto
        r = newton (f1, 1.0, { .converged = cvg1 });
        CHECK(r == doctest::Approx { target1 });
    }
    SUBCASE("newton, throws if zero derivative")
    {
        CHECK_THROWS_AS(
              newton (f1, [](auto){ return 0.0; }, 1.0)
            , newton_zero_derivative_e
        );
    }
    SUBCASE("newton, user function throws")
    {
        CHECK_THROWS_AS(
              newton (f3, 1.0)
            , int
        );
    }
    SUBCASE("newton, with options")
    {
        CHECK_THROWS_AS(
              newton (f1, 1.0, { .max_iter = 1 })
            , newton_no_convergence_e
        );
    }
    SUBCASE("newton, with info (iteration count)")
    {
            auto const
        [ result, info ] = newton (f1, 1.0, { /*default options*/ }, info::iterations);
        CHECK(info.iteration_count > 1);
    }
    SUBCASE("newton, with info (convergence)")
    {
            auto const
        [ result, info ] = newton (f1, 1.0, { /*default options*/ }, info::convergence);
        CHECK(info.convergence.size () > 1);
        for (auto&& [v, f ,df]: info.convergence)
        {
            MESSAGE (v, ", ", f, ", ", df);
        }
    }
    SUBCASE("newton, with info convergence does not throw no_convergence_e!")
    {
            auto const
        [ result, info ] = newton (f1, 1.0, { .max_iter = 3 }, info::convergence);
        CHECK(info.convergence.size () == 3);
        CHECK(info.converged == false);
        for (auto&& [v, f ,df]: info.convergence)
        {
            MESSAGE (v, ", ", f, ", ", df);
        }
    }
    SUBCASE("newton, user function throws, with info")
    {
            auto const
        [ result, info ] = newton (f3, 1.0, {}, info::iterations);
        CHECK(!info.converged);
        CHECK(info.function_threw);
    }
    SUBCASE("Use alternative convergence predicate")
    {
            auto
        r = newton (f1, 1.0, { .converged = make_newton_simple_converged (1e-8) });
        CHECK(r == doctest::Approx { target1 });
    }
}
// -----------------------------------------------------------------------------
TEST_CASE("Zhang")
{
        auto
    cvg2 = [](
          double a
        , double b
        , double fa
        , double fb
    ){
            constexpr double
        tol = 1e-12;
        return
               fa == 0.0
            || fb == 0.0
            || fabs (b - a) < tol
        ;
    };

    SUBCASE("zhang")
    {
            auto
        r = zhang (f1, 0.0, 10.0);
        CHECK(r == doctest::Approx { target1 });

            auto
        s = zhang (f2, 0.0, 10.0);
        CHECK(s == doctest::Approx { target2 });
    }
    SUBCASE("zhang, with custom stopping criterion")
    {
            auto
        r = zhang (f1, 0.0, 10.0, { .converged = cvg2 });
        CHECK(r == doctest::Approx { target1 });
    }
    SUBCASE("zhang, throws if no single root between brackets")
    {
        CHECK_THROWS_AS(
              zhang (f1, 0.0, 0.1);
            , zhang_no_single_root_between_brackets_e
        );
    }
    SUBCASE("zhang, user function throws")
    {
        CHECK_THROWS_AS(
              zhang (f3, 0.0, 0.1);
            , int
        );
    }
    SUBCASE("zhang, with options")
    {
        CHECK_THROWS_AS(
              zhang (f1, 0.0, 10.0, { .max_iter = 1 })
            , zhang_no_convergence_e
        );
    }
    SUBCASE("zhang, with info (iteration count)")
    {
            auto const
        [ result, info ] = zhang (f1, 0.0, 10.0, { /*default options*/ }, info::iterations);
        CHECK(info.iteration_count > 1);
    }
    SUBCASE("zhang, with info (convergence)")
    {
            auto const
        [ result, info ] = zhang (f1, 0.0, 10.0, { /*default options*/ }, info::convergence);
        CHECK(info.convergence.size () > 1);
        for (auto&& [ a, b, fa, fb ]: info.convergence)
        {
            MESSAGE (a, ", ", b, ", ", fa, ", ", fb);
        }
    }
    SUBCASE("zhang, with info convergence does not throw no_convergence_e!")
    {
            auto const
        [ result, info ] = zhang (f1, 0.0, 10.0, { .max_iter = 3 }, info::convergence);
        CHECK(info.convergence.size () == 3);
        CHECK(info.converged == false);
        for (auto&& [ a, b, fa, fb ]: info.convergence)
        {
            MESSAGE (a, ", ", b, ", ", fa, ", ", fb);
        }
    }
    SUBCASE("In fact, with any info, never throw")
    {
            auto const
        [ result, info ] = zhang (f1, 0.0, 0.1, {}, info::iterations);
        CHECK(!info.converged);
        CHECK(info.no_single_root_between_bracket);
    }
    SUBCASE("zhang, user function throws, with info")
    {
            auto const
        [ result, info ] = zhang (f3, 0.0, 0.1, {}, info::iterations);
        CHECK(!info.converged);
        CHECK(info.function_threw);
    }
    SUBCASE("zhang, with generated stopping criterion")
    {
            auto
        r = zhang (f1, 0.0, 10.0, { .converged = make_zhang_simple_converged (1e-8) });
        CHECK(r == doctest::Approx { target1 });
    }
}
// -----------------------------------------------------------------------------
    auto
f5 = [](auto x) { return x * x; };
    auto
f6 = [](auto x) { return x; };
TEST_CASE("Bracket minimum")
{
    SUBCASE("bracket_minimum")
    {
            auto
        [a, b, fa, fb] = bracket_minimum (f5, -10., -11.);
        CHECK(a < 0.);
        CHECK(b > 0.);
    }
    SUBCASE("bracket_minimum, throws if no single root between brackets")
    {
        CHECK_THROWS_AS(
              bracket_minimum (f6, 0.0, 0.1);
            , bracket_minimum_no_convergence_e
        );
    }
    SUBCASE("bracket_minimum, user function throws")
    {
        CHECK_THROWS_AS(
              bracket_minimum (f3, 0.0, 0.1);
            , int
        );
    }
    SUBCASE("bracket_minimum, with options")
    {
        CHECK_THROWS_AS(
              bracket_minimum (f5, -10., -9., { .max_iter = 1 })
            , bracket_minimum_no_convergence_e
        );
    }
    SUBCASE("bracket_minimum, with info (iteration count)")
    {
            auto const
        [ result, info ] = bracket_minimum (f5, 10.0, 11.0, { /*default options*/ }, info::iterations);
        CHECK(info.iteration_count > 1);
    }
    SUBCASE("bracket_minimum, with info (convergence)")
    {
            auto const
        [ result, info ] = bracket_minimum (f5, 10.0, 11.0, { /*default options*/ }, info::convergence);
        CHECK(info.convergence.size () > 1);
        for (auto&& [ a, b, c ]: info.convergence)
        {
            MESSAGE (fmt::format (
                  "[{}, {}], [{}, {}], [{}, {}]"
                , a.first, a.second
                , b.first, b.second
                , c.first, c.second
            ));
        }
        MESSAGE(fmt::format(
              "Result: [{}, {}], [{}, {}]"
            , std::get <0> (result)
            , std::get <2> (result)
            , std::get <1> (result)
            , std::get <3> (result)
        ));
    }
    SUBCASE("bracket_minimum, with info convergence does not throw no_convergence_e!")
    {
            auto const
        [ result, info ] = bracket_minimum (f5, 20e4, 20e4 + 1e-4, { .max_iter = 2 }, info::convergence);
        CHECK(info.convergence.size () == 2);
        CHECK(info.converged == false);
    }
    SUBCASE("bracket_minimum, user function throws, with info")
    {
            auto const
        [ result, info ] = bracket_minimum (f3, 0.0, 0.1, {}, info::iterations);
        CHECK(!info.converged);
        CHECK(info.function_threw);
    }
}
// -----------------------------------------------------------------------------
TEST_CASE("Golden section")
{
    SUBCASE("golden_section")
    {
            auto
        r = golden_section (f5, -10., -11.);
        CHECK(r == doctest::Approx { 0. });
    }
    SUBCASE("golden_section, throws if no single root between brackets")
    {
        CHECK_THROWS_AS(
              golden_section (f6, 0.0, 0.1);
            , bracket_minimum_no_convergence_e
        );
    }
    SUBCASE("golden_section, user function throws")
    {
        CHECK_THROWS_AS(
              golden_section (f3, 0.0, 0.1);
            , int
        );
    }
    SUBCASE("golden_section, with options")
    {
            auto
        r = golden_section (f5, -10., -9., { .tolerance = 1e-6 });
        CHECK(r == doctest::Approx { 0. });
    }
    SUBCASE("golden_section, with info (iteration count)")
    {
            auto const
        [ result, info ] = golden_section (f5, 10., 11., { /*default options*/ }, info::iterations);
        CHECK(info.iteration_count > 1);
        CHECK(info.bracket_minimum_info.iteration_count > 1);
    }
    SUBCASE("golden_section, with info (convergence)")
    {
            auto const
        [ result, info ] = golden_section (f5, 10., 11., { /*default options*/ }, info::convergence);
        CHECK(info.convergence.size () > 1);
        MESSAGE("Bracket minimum convergence:");
        for (auto&& [ a, c, b ]: info.bracket_minimum_info.convergence)
        {
            MESSAGE (fmt::format (
                  "  [{}, {}], [{}, {}], [{}, {}]"
                , a.first, a.second
                , c.first, c.second
                , b.first, b.second
            ));
        }
        MESSAGE("Golden section convergence:");
        for (auto&& [ a, c, d, b ]: info.convergence)
        {
            MESSAGE (fmt::format (
                  "  [{}, {}], [{}, {}], [{}, {}], [{}, {}]"
                  , a.first, a.second
                  , c.first, c.second
                  , d.first, d.second
                  , b.first, b.second
            ));
        }
        MESSAGE("Result: ", result);
    }
    SUBCASE("golden_section, user function throws, with info")
    {
            auto const
        [ result, info ] = golden_section (f3, 0.0, 0.1, {}, info::iterations);
        CHECK(!info.converged);
        CHECK(info.bracket_minimum_info.function_threw);
    }
    SUBCASE("golden_section, passing options to bracket_minimum")
    {
            auto const
        [ result, info ] = golden_section (f5, 10., 11., { .bracket_minimum_options = { .max_iter = 1 } }, info::convergence);
        CHECK(!info.converged);
        CHECK(!info.bracket_minimum_info.converged);
    }
}
// -----------------------------------------------------------------------------
TEST_CASE("Powel")
{
        auto
    rosenbrock = [](std::valarray <double> const& x)
    {
            using std::pow;
        return pow (1. - x[0], 2.) + 100. * pow (x[1] - pow (x[0], 2.), 2.);
    };
    SUBCASE("powell")
    {
            auto
        [ r, info ] = powell (rosenbrock, std::valarray { 0.1, 0.1 }, {}, info::convergence);
        CHECK(r[0] == doctest::Approx { 1. });
        CHECK(r[1] == doctest::Approx { 1. });
        for (auto&& [ i, j, f, p ]: info.convergence)
        {
            MESSAGE("iter: ", i, ", direction: ", j, ", f= ", f, ", p= ", p);
            MESSAGE("Bracket minimum convergence:");
            for (auto&& [a, b, c]: info.golden_section_info.at (i).bracket_minimum_info.convergence)
            {
                MESSAGE (fmt::format (
                      "  [{}, {}], [{}, {}], [{}, {}]"
                    , a.first, a.second
                    , b.first, b.second
                    , c.first, c.second
                ));
            }
            MESSAGE("Golden section search convergence:");
            for (auto&& [a, b, c, d]: info.golden_section_info.at (i).convergence)
            {
                MESSAGE (fmt::format (
                      "  [{}, {}], [{}, {}], [{}, {}], [{}, {}]"
                      , a.first, a.second
                      , c.first, c.second
                      , d.first, d.second
                      , b.first, b.second
                ));
            }
        }

    }
}
// -----------------------------------------------------------------------------
TEST_CASE("Multidimensional Newton with automatic differenciation")
{
        auto
    mf = [] <class T> (Eigen::Matrix <T, 2, 1> const& x)
    {
            auto
        r = Eigen::Matrix <T, 2, 1> {};
            const auto&
        b = x (0);
            const auto&
        h = x (1);
        r(0) = 2272000 / h / b + 2 * b + 5 + (h + b + 10) * (-2272000 / h / b / b + 2);
        r(1) = 2272000 / h / b + 2 * b + 5 - (h + b + 10) * (2272000 / h / h / b);
        return r;
    };
        auto
    mg = [] <class T> (Eigen::Matrix <T, 2, 1> const& x)
    {
            auto
        r = Eigen::Matrix <T, 2, 1> {};
            const auto&
        X = x (0);
            const auto&
        Y = x (1);
            constexpr auto
        a = 1.;
            constexpr auto
        b = 10.;
        r (0) = a * (1 - X);
        r (1) = b * (Y - X * X);
        return r;
    };
    SUBCASE("newton")
    {
            const auto
        x0 = Eigen::Vector2d { 100., 100. };
            auto
        r = newton <2> (mf, x0);
        CHECK(r[0] == doctest::Approx { 65.67176514794734 });
        CHECK(r[1] == doctest::Approx { 138.56848974267868 });
            const auto
        y0 = Eigen::Vector2d { -10., -5. };
            auto
        s = newton <2> (mg, y0);
        CHECK(s[0] == doctest::Approx { 1 });
        CHECK(s[1] == doctest::Approx { 1 });
    }
    SUBCASE("newton, with custom convergence criterion")
    {
            auto
        r = newton (f1, 1.0, { .converged = cvg1 });
        CHECK(r == doctest::Approx { target1 });
    }
    SUBCASE("newton, with info (iteration count)")
    {
            auto const
        [ result, info ] = newton (f1, 1.0, { /*default options*/ }, info::iterations);
        CHECK(info.iteration_count > 1);
    }
    #if 0
    SUBCASE("newton, throws if zero derivative")
    {
        CHECK_THROWS_AS(
              newton (f1, [](auto){ return 0.0; }, 1.0)
            , newton_zero_derivative_e
        );
    }
    SUBCASE("newton, user function throws")
    {
        CHECK_THROWS_AS(
              newton (f3, 1.0)
            , int
        );
    }
    SUBCASE("newton, with options")
    {
        CHECK_THROWS_AS(
              newton (f1, 1.0, { .max_iter = 1 })
            , newton_no_convergence_e
        );
    }
    SUBCASE("newton, with info (iteration count)")
    {
            auto const
        [ result, info ] = newton (f1, 1.0, { /*default options*/ }, info::iterations);
        CHECK(info.iteration_count > 1);
    }
    SUBCASE("newton, with info (convergence)")
    {
            auto const
        [ result, info ] = newton (f1, 1.0, { /*default options*/ }, info::convergence);
        CHECK(info.convergence.size () > 1);
        for (auto&& [v, f ,df]: info.convergence)
        {
            MESSAGE (v, ", ", f, ", ", df);
        }
    }
    SUBCASE("newton, with info convergence does not throw no_convergence_e!")
    {
            auto const
        [ result, info ] = newton (f1, 1.0, { .max_iter = 3 }, info::convergence);
        CHECK(info.convergence.size () == 3);
        CHECK(info.converged == false);
        for (auto&& [v, f ,df]: info.convergence)
        {
            MESSAGE (v, ", ", f, ", ", df);
        }
    }
    SUBCASE("newton, user function throws, with info")
    {
            auto const
        [ result, info ] = newton (f3, 1.0, {}, info::iterations);
        CHECK(!info.converged);
        CHECK(info.function_threw);
    }
    SUBCASE("Use alternative convergence predicate")
    {
            auto
        r = newton (f1, 1.0, { .converged = make_newton_simple_converged (1e-8) });
        CHECK(r == doctest::Approx { target1 });
    }
    #endif

}
// -----------------------------------------------------------------------------
TEST_CASE("Multidimensional Newton with jacobian")
{
        auto
    mf1 = [] <class T> (Eigen::Matrix <T, 2, 1> const& x)
    {
            auto
        r = Eigen::Matrix <T, 2, 1> {};
            const auto&
        b = x (0);
            const auto&
        h = x (1);
        r (0) = 2 * h + 4 * b + 25 - 2272000 / b / b - 10 * 2272000 / h / b / b;
        r (1) = 2 * b + 5 - 2272000 / h / h - 10 * 2272000 / h / h / b;
        return r;
    };
        auto
    mj1 = [] <class T> (Eigen::Matrix <T, 2, 1> const& x)
    {
            auto
        r = Eigen::Matrix <T, 2, 2> {};
            const auto&
        b = x (0);
            const auto&
        h = x (1);
        r (0, 0) = 4 + 2 * 2272000 / b / b / b + 20 * 2272000  / (h * h * h * b);
        r (1, 1) = 2 * 2272000 * (b + 10) / (h * h * h * b);
        r (0, 1) = 2 + 10 * 2272000 / (h * h * h * b);
        r (1, 0) = r (0, 1);
        return r;
    };
        auto
    mf2 = [] <class T> (Eigen::Matrix <T, 2, 1> const& x)
    {
            auto
        r = Eigen::Matrix <T, 2, 1> {};
            const auto&
        X = x (0);
            const auto&
        Y = x (1);
            constexpr auto
        a = 1.;
            constexpr auto
        b = 10.;
        r (0) = a * (1 - X);
        r (1) = b * (Y - X * X);
        return r;
    };
        auto
    mj2 = [] <class T> (Eigen::Matrix <T, 2, 1> const& x)
    {
            auto
        r = Eigen::Matrix <T, 2, 2> {};
            const auto&
        X = x (0);
            constexpr auto
        a = 1.;
            constexpr auto
        b = 10.;
        r (0, 0) = -a;
        r (1, 1) = b;
        r (0, 1) = 0;
        r (1, 0) = -2 * b * X;
        return r;
    };
    SUBCASE("newton")
    {
            const auto
        x0 = Eigen::Vector2d { 100., 100. };
            auto
        r = newton <2> (mf1, mj1, x0);
        CHECK(r[0] == doctest::Approx { 65.67176514794734 });
        CHECK(r[1] == doctest::Approx { 138.56848974267868 });
            const auto
        y0 = Eigen::Vector2d { -10., -5. };
            auto
        s = newton <2> (mf2, mj2, y0);
        CHECK(s[0] == doctest::Approx { 1 });
        CHECK(s[1] == doctest::Approx { 1 });
    }
    #if 0
    SUBCASE("newton, with info (convergence)")
    {
            const auto
        x0 = Eigen::Vector2d { 100., 100. };
            auto const
        [ result, info ] = newton <2> (mf1, mj1, x0, { /*default options*/ }, info::convergence);
        for (auto const& [x, f, j]: info.convergence)
        {
            MESSAGE(x, f);
        }
    }
    SUBCASE("newton, with custom convergence criterion")
    {
            auto
        r = newton (mf1, 1.0, { .converged = cvg1 });
        CHECK(r == doctest::Approx { target1 });
    }
    SUBCASE("newton, throws if zero derivative")
    {
        CHECK_THROWS_AS(
              newton (f1, [](auto){ return 0.0; }, 1.0)
            , newton_zero_derivative_e
        );
    }
    SUBCASE("newton, user function throws")
    {
        CHECK_THROWS_AS(
              newton (f3, 1.0)
            , int
        );
    }
    SUBCASE("newton, with options")
    {
        CHECK_THROWS_AS(
              newton (f1, 1.0, { .max_iter = 1 })
            , newton_no_convergence_e
        );
    }
    SUBCASE("newton, with info (iteration count)")
    {
            auto const
        [ result, info ] = newton (f1, 1.0, { /*default options*/ }, info::iterations);
        CHECK(info.iteration_count > 1);
    }
    SUBCASE("newton, with info (convergence)")
    {
            auto const
        [ result, info ] = newton (f1, 1.0, { /*default options*/ }, info::convergence);
        CHECK(info.convergence.size () > 1);
        for (auto&& [v, f ,df]: info.convergence)
        {
            MESSAGE (v, ", ", f, ", ", df);
        }
    }
    SUBCASE("newton, with info convergence does not throw no_convergence_e!")
    {
            auto const
        [ result, info ] = newton (f1, 1.0, { .max_iter = 3 }, info::convergence);
        CHECK(info.convergence.size () == 3);
        CHECK(info.converged == false);
        for (auto&& [v, f ,df]: info.convergence)
        {
            MESSAGE (v, ", ", f, ", ", df);
        }
    }
    SUBCASE("newton, user function throws, with info")
    {
            auto const
        [ result, info ] = newton (f3, 1.0, {}, info::iterations);
        CHECK(!info.converged);
        CHECK(info.function_threw);
    }
    SUBCASE("Use alternative convergence predicate")
    {
            auto
        r = newton (f1, 1.0, { .converged = make_newton_simple_converged (1e-8) });
        CHECK(r == doctest::Approx { target1 });
    }
    #endif

}
