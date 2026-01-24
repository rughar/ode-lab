#include <iostream>
#include <ricatti.hpp>



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
};      



int main() {
   
    verlet::ricatti_data<double> rd;
    rd.init(2);
   
    LotkaVolterra rc;
   
    std::vector<double> x = {1.0, 1.0};

    for (size_t i = 0; i < 100; ++i) 
        rc.step(0.1, x, rd);

    std::cout << "x after step: " << x[0] << std::endl;

    
    return 0;

}