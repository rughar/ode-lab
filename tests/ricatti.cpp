#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <ricatti.hpp>

TEST_CASE("ricatti") 
{
class LotkaVolterra_core : public verlet::ricatti_core<double> 
{
    public:
        double invariant() {
            return (coef->C[1][0][1] + coef->C[1][1][0]) * x[0] - (coef->C[0][0][1] + coef->C[0][1][0]) * x[1] + coef->B[1][1] * std::log(x[0]) - coef->B[0][0] * std::log(x[1]);
        }
};      
    verlet::ricatti_coef<double>  coefficients;
    coefficients.init(2);

    coefficients.B[0][0] = 2.0/3.0;      // Growth rate of prey
    coefficients.B[1][1] = -1.0;         // Death rate of predators        
    coefficients.C[0][0][1] = -4.0/3.0;  // Non-linear effect on prey
    coefficients.C[1][0][1] = 1.0;       // Non-linear effect on predators

    coefficients.symmetrize_C();

    LotkaVolterra_core core;
    core.init(&coefficients);
    core.x = {1.0, 1.0};
    
    double i1 = core.invariant();

    for(double h = 1.0; h > 1e-7; h *= 0.1) {
        core.step(h);
        double i2 = core.invariant();
        REQUIRE(i2 == Catch::Approx(i1).epsilon(0.01 * h * h));
        i1 = i2;
    }
}