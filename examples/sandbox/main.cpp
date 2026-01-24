#include <iostream>
#include <ricatti.hpp>
#include <fstream>



class LotkaVolterra : public verlet::ricatti_core<double> {
    public:
        LotkaVolterra() 
        {
            init(2);   
            set_B();
            set_C(); 
        }

        void set_B() override 
        {
            B[0][0] = 2.0/3.0;      // Growth rate of prey
            B[1][1] = -1.0;         // Death rate of predators
        }

        void set_C() override 
        {
            C[0][0][1] = -4.0/3.0;  // Non-linear effect on prey
            C[1][0][1] = 1.0;       // Non-linear effect on predators

            symmetrize_C();
        }

        template<class V>
        double invariant(const V& x) {
            return (C[1][0][1] + C[1][1][0]) * x[0] - (C[0][0][1] + C[0][1][0]) * x[1] + B[1][1] * std::log(x[0]) - B[0][0] * std::log(x[1]);
        }
};      



int main() {
   
    std::ofstream out("trajectory.txt");

    verlet::ricatti_data<double> rd;
    rd.init(2);
   
    LotkaVolterra rc;
   
    std::vector<double> x = {1.0, 1.0};

    for (size_t i = 0; i < 1000; ++i) {
        rc.step(0.1, x, rd);
        out << rc.invariant(x) << std::endl;
    }

    out.close();

    return 0;

}