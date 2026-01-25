#include <iostream>
#include <ricatti.hpp>
#include <fstream>

int main() {

class Lhotka_Voltera : public verlet::ricatti_core<double> 
    {
        public:
            Lhotka_Voltera() : verlet::ricatti_core<double>(2) {};

            void set_coef() override
            {
                A_COEF(0) = -0.01;
  
                B_COEF(0, 0) = 2.0/3.0;
                B_COEF(1, 1) = -1.0;
 
                C_COEF(0, 0, 1) = -4.0/3.0;
                C_COEF(1, 0, 1) = 1.0;
            }   
    };

    Lhotka_Voltera core;
    core.u = {1.0, 1.0};
   
    std::ofstream out("trajectory.txt");

    for(size_t i = 0; i < 1000; ++i) {
        core.step(0.1);
        out << core.u[0] << " " << core.u[1] << std::endl;
    }
    out.close();

    return 0;
}