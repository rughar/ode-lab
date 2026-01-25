#include <iostream>
#include <verlet.hpp>
#include <fstream>


int main() {

    verlet::verlet_core<double> core;   
    core.init(2);
    core.x[0] = 1.0;
    core.step(0.1);

    std::ofstream out("trajectory.txt");

    //out << core.invariant() << std::endl;
    out.close();

    return 0;

}