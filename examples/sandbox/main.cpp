#include <iostream>
#include <verlet.hpp>
#include <fstream>


//class Oscillator : public verlet::verlet_core<double> 
//{
//    private:
//        
//        verlet::ricatti_coef<double> coef;
//
//    public:
        
//        Oscillator() 
//        {
//            init(1);
//            coef.init(1);
//            attach_coef(&coef);
//        }
 
//        void update_coefs() override
//        {
//            coef.A[0] = -x[0];
//        }
//};

int main() {


//    Oscillator osc;
//    osc.x[0] = 1.0;

    std::ofstream out("trajectory.txt");

//    for(size_t i = 0; i< 1000; ++i) {
//        out << osc.x[0] * osc.x[0] + osc.u[0] * osc.u[0] << std::endl;
//        osc.step(0.1);
//    }   
    out.close();

    return 0;

}