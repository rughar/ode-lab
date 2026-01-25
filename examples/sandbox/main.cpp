#include <iostream>
#include <ricatti.hpp>
//#include <verlet.hpp>
#include <fstream>

class LotkaVolterra_core : public verlet::ricatti_core<double> 
{
    public:
        double invariant() {
            return (coef->C[1][0][1] + coef->C[1][1][0]) * x[0] - (coef->C[0][0][1] + coef->C[0][1][0]) * x[1] + coef->B[1][1] * std::log(x[0]) - coef->B[0][0] * std::log(x[1]);
        }
};      

int main() {

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

    std::ofstream out("trajectory.txt");
    
    for (size_t i = 0; i < 100; ++i) {
        core.step(0.1);
        out << core.invariant() << std::endl;
    }

    out.close();

    return 0;

}