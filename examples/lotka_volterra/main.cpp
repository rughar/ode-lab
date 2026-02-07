#include <iostream>
#include <qode1.hpp>
#include <fstream>

int main()
{

    class Lotka_Voltera : public qode::qode1_core<double>
    {
    public:
        Lotka_Voltera() : qode::qode1_core<double>(2) {};

        void set_coef() override
        {
            b_coef(0, 0) = 2.0 / 3.0;
            b_coef(1, 1) = -1.0;

            c_coef(0, 0, 1) = -4.0 / 3.0;
            c_coef(1, 0, 1) = 1.0;
        }
    };

    Lotka_Voltera core;
    core.x = {1.0, 1.0};

    std::ofstream out("trajectory.txt");

    double h = core.suggest_first_stepsize(1.0, 0.03);

    for (size_t i = 0; i < 500; ++i)
    {
        // core.step(0.1);
        core.step_adaptive(h, 0.03);
        // out << core.u[0] << " " << core.u[1] << std::endl;
        out << h << std::endl;
    }
    out.close();

    return 0;
}