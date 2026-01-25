#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <ricatti.hpp>

TEST_CASE("ricatti") 
{
    class Lhotka_Voltera : public verlet::ricatti_core<double> 
    {
        public:
            Lhotka_Voltera() : verlet::ricatti_core<double>(2) {};

            void set_b_coef() override
            {
                B_COEF(0, 0) = 2.0/3.0;
                B_COEF(1, 1) = -1.0;
            }

            void set_c_coef() override
            {
                C_COEF(0, 0, 1) = -4.0/3.0;
                C_COEF(1, 0, 1) = 1.0;
            }   
    };

    Lhotka_Voltera core;
    core.u = { 1.0, 1.0 };
   
    auto invariant = [](const std::vector<double>& u) -> double 
    {
        return u[0] - std::log(u[0]) + 4.0/3.0 * u[1] - 2.0/3.0 * std::log(u[1]);
    };
 
    double i1 = invariant(core.u);

    for(double h = 0.1; h > 1e-7; h *= 0.1) {
        core.step(h);
        double i2 = invariant(core.u);
        REQUIRE(i2 == Catch::Approx(i1).epsilon(0.01 * h * h));
        i1 = i2;
    }
}