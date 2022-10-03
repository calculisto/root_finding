#include <doctest/doctest.h>
#include "../include/isto/root_finding/root_finding.hpp"
    using namespace isto::root_finding;
#include <cmath>
    using std::cos, std::sin, std::pow;

TEST_CASE("root_finding.hpp")
{
        auto
    f1 = [](double x){ return cos (x) - pow (x, 3.0); };
        auto
    df1 = [](double x){ return -sin (x) - 3 * pow (x, 2.0); };
        auto const
    target1 = 0.8654740331016144466206859011862287477929;
     
        auto
    f2 = [](double x){ return cos (x) - x; };
        auto
    df2 = [](double x){ return -sin (x) - 1; };
        auto const
    target2 = 0.7390851332151606416553120876738734040134 ;

        auto
    f3 = [](double){ throw int {}; return 1.; };
        auto
    df3 = [](double) { return 0; };

        auto
    f4 = [](double){ return 1.; };
        auto
    df4 = [](double) { throw int {}; return 1.; };

        auto
    cvg1 = [](double x){ return fabs (x) < 1e-12; };

    SUBCASE("newton")
    {
            auto
        r = newton (f1, df1, 1.0, cvg1);
        CHECK(r == doctest::Approx { target1 });

            auto
        s = newton (f2, df2, 1.0, cvg1); 
        CHECK(s == doctest::Approx { target2 });
    }
    SUBCASE("newton, throws if zero derivative")
    {
        CHECK_THROWS_AS(
              newton (f1, [](auto){ return 0.0; }, 1.0, cvg1)
            , newton_zero_derivative_e
        );
    }
    SUBCASE("newton, user function throws")
    {
        CHECK_THROWS_AS(
              newton (f3, df3, 1.0, cvg1)
            , int
        );
        CHECK_THROWS_AS(
              newton (f4, df4, 1.0, cvg1)
            , int
        );
    }
    SUBCASE("newton, with options")
    {
        CHECK_THROWS_AS(
              newton (f1, df1, 1.0, cvg1, { .max_iter = 1 })
            , newton_no_convergence_e
        );
    }
    SUBCASE("newton, with info (iteration count)")
    {
            auto const
        [ result, info ] = newton (f1, df1, 1.0, cvg1, { /*default options*/ }, info::iterations);
        CHECK(info.iteration_count > 1);
    }
    SUBCASE("newton, with info (convergence)")
    {
            auto const
        [ result, info ] = newton (f1, df1, 1.0, cvg1, { /*default options*/ }, info::convergence);
        CHECK(info.convergence.size () > 1);
        for (auto&& [v, f ,df]: info.convergence)
        {
            MESSAGE (v, ", ", f, ", ", df);
        }
    }
    SUBCASE("newton, with info convergence does not throw no_convergence_e!")
    {
            auto const
        [ result, info ] = newton (f1, df1, 1.0, cvg1, { .max_iter = 3 }, info::convergence);
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
        [ result, info ] = newton (f1, [](auto){ return 0.0; }, 1.0, cvg1, {}, info::iterations);
        CHECK(!info.converged);
        CHECK(info.zero_derivative);
    }
    SUBCASE("newton, user function throws, with info")
    {
            auto const
        [ result, info ] = newton (f3, df3, 1.0, cvg1, {}, info::iterations);
        CHECK(!info.converged);
        CHECK(info.function_threw);
    }
    SUBCASE("newton, user derivative throws, with info")
    {
            auto const
        [ result, info ] = newton (f4, df4, 1.0, cvg1, {}, info::iterations);
        CHECK(!info.converged);
        CHECK(info.derivative_threw);
    }

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
        r = zhang (f1, 0.0, 10.0, cvg2);
        CHECK(r == doctest::Approx { target1 });

            auto
        s = zhang (f2, 0.0, 10.0, cvg2); 
        CHECK(s == doctest::Approx { target2 });
    }
    SUBCASE("zhang, throws if no single root between brackets")
    {
        CHECK_THROWS_AS(
              zhang (f1, 0.0, 0.1, cvg2);
            , zhang_no_single_root_between_brackets_e
        );
    }
    SUBCASE("zhang, user function throws")
    {
        CHECK_THROWS_AS(
              zhang (f3, 0.0, 0.1, cvg2);
            , int
        );
    }
    SUBCASE("zhang, with options")
    {
        CHECK_THROWS_AS(
              zhang (f1, 0.0, 10.0, cvg2, { .max_iter = 1 })
            , zhang_no_convergence_e
        );
    }
    SUBCASE("zhang, with info (iteration count)")
    {
            auto const
        [ result, info ] = zhang (f1, 0.0, 10.0, cvg2, { /*default options*/ }, info::iterations);
        CHECK(info.iteration_count > 1);
    }
    SUBCASE("zhang, with info (convergence)")
    {
            auto const
        [ result, info ] = zhang (f1, 0.0, 10.0, cvg2, { /*default options*/ }, info::convergence);
        CHECK(info.convergence.size () > 1);
        for (auto&& [ a, b, fa, fb ]: info.convergence)
        {
            MESSAGE (a, ", ", b, ", ", fa, ", ", fb);
        }
    }
    SUBCASE("zhang, with info convergence does not throw no_convergence_e!")
    {
            auto const
        [ result, info ] = zhang (f1, 0.0, 10.0, cvg2, { .max_iter = 3 }, info::convergence);
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
        [ result, info ] = zhang (f1, 0.0, 0.1, cvg2, {}, info::iterations);
        CHECK(!info.converged);
        CHECK(info.no_single_root_between_bracket);
    }
    SUBCASE("zhang, user function throws, with info")
    {
            auto const
        [ result, info ] = zhang (f3, 0.0, 0.1, cvg2, {}, info::iterations);
        CHECK(!info.converged);
        CHECK(info.function_threw);
    }

        auto
    f5 = [](auto x) { return x * x; };
        auto
    f6 = [](auto x) { return x; };
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
        for (auto&& [ a, fa ]: info.convergence)
        {
            MESSAGE (a, ", ", fa);
        }
    }
    SUBCASE("bracket_minimum, with info convergence does not throw no_convergence_e!")
    {
            auto const
        [ result, info ] = bracket_minimum (f5, 20e4, 20e4 + 1e-4, { .max_iter = 2 }, info::convergence);
        CHECK(info.convergence.size () == 4);
        CHECK(info.converged == false);
        for (auto&& [ a, fa ]: info.convergence)
        {
            MESSAGE (a, ", ", fa);
        }
    }
    SUBCASE("bracket_minimum, user function throws, with info")
    {
            auto const
        [ result, info ] = bracket_minimum (f3, 0.0, 0.1, {}, info::iterations);
        CHECK(!info.converged);
        CHECK(info.function_threw);
    }
// -----------------------------------------------------------------------------
    SUBCASE("golden_section")
    {
            auto
        r = golden_section (f5, -10., -11., 1e-10);
        CHECK(r == doctest::Approx { 0. });
    }
    SUBCASE("golden_section, throws if no single root between brackets")
    {
        CHECK_THROWS_AS(
              golden_section (f6, 0.0, 0.1, 1e-10);
            , bracket_minimum_no_convergence_e
        );
    }
    SUBCASE("golden_section, user function throws")
    {
        CHECK_THROWS_AS(
              golden_section (f3, 0.0, 0.1, 1e-10);
            , int
        );
    }
    // gss has no max_iter option.
    /*
    SUBCASE("golden_section, with options")
    {
        CHECK_THROWS_AS(
              golden_section (f5, -10., -9., 1e-10, { .max_iter = 1 })
            , golden_section_no_convergence_e
        );
    }
    */
    SUBCASE("golden_section, with info (iteration count)")
    {
            auto const
        [ result, info ] = golden_section (f5, 10., 11., 1e-10, { /*default options*/ }, info::iterations);
        CHECK(info.iteration_count > 1);
        CHECK(info.bracket_minimum_info.iteration_count > 1);
    }
    SUBCASE("golden_section, with info (convergence)")
    {
            auto const
        [ result, info ] = golden_section (f5, 10., 11., 1e-10, { /*default options*/ }, info::convergence);
        CHECK(info.convergence.size () > 1);
        for (auto&& [ a, c, d, b ]: info.convergence)
        {
            MESSAGE (a, ", ", c, ", ", d, ", ", b);
        }
        for (auto&& [ a, fa ]: info.bracket_minimum_info.convergence)
        {
            MESSAGE (a, ", ", fa);
        }
    }
    /*
    SUBCASE("golden_section, with info convergence does not throw no_convergence_e!")
    {
            auto const
        [ result, info ] = golden_section (f5, 20e4, 20e4 + 1e-4, 1e-10, { .max_iter = 2 }, info::convergence);
        CHECK(info.convergence.size () == 4);
        CHECK(info.converged == false);
        for (auto&& [ a, fa ]: info.convergence)
        {
            MESSAGE (a, ", ", fa);
        }
    }
    */
    SUBCASE("golden_section, user function throws, with info")
    {
            auto const
        [ result, info ] = golden_section (f3, 0.0, 0.1, 1e-10, {}, info::iterations);
        CHECK(!info.converged);
        CHECK(info.bracket_minimum_failed);
    }
    SUBCASE("golden_section, passing options to bracket_minimum")
    {
            auto const
        [ result, info ] = golden_section (f5, 10., 11., 1e-10, { .bracket_minimum_options = { .max_iter = 1 } }, info::convergence);
        CHECK(!info.converged);
        CHECK(info.bracket_minimum_failed);
    }
// -----------------------------------------------------------------------------
}
