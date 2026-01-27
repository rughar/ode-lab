#include <iostream>
#include <ricatti.hpp>
#include <fstream>

int main()
{

    class Lhotka_Voltera : public verlet::ricatti_core<double>
    {
    public:
        Lhotka_Voltera() : verlet::ricatti_core<double>(2) {};

        void set_coef() override
        {
            a_coef(0) = -0.0;

            b_coef(0, 0) = 2.0 / 3.0;
            b_coef(1, 1) = -1.0;

            c_coef(0, 0, 1) = -4.0 / 3.0;
            c_coef(1, 0, 1) = 1.0;
        }
    };

    Lhotka_Voltera core;
    core.u = {1.0, 1.0};

    std::ofstream out("trajectory.txt");

    double h = core.suggest_first_stepsize(1.0, 0.1);

    for (size_t i = 0; i < 500; ++i)
    {
        // core.step(0.1);
        core.step_adaptive(h, 0.1, 0.3, 2.0);
        // out << core.u[0] << " " << core.u[1] << std::endl;
        out << h << std::endl;
    }
    out.close();

    return 0;
}