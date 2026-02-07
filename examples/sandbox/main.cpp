#include <iostream>
#include <verlet.hpp>
#include <fstream>
#include <minijacobian.hpp>
#include <rkgl.hpp>

int main()
{

  //    Oscillator osc;
  //    osc.x[0] = 1.0;

  //    std::ofstream out("trajectory.txt");

  //    for(size_t i = 0; i< 1000; ++i) {
  //        out << osc.x[0] * osc.x[0] + osc.u[0] * osc.u[0] << std::endl;
  //        osc.step(0.1);
  //    }
  //    out.close();
  rkgl::rkgl<double, 3> core;
  core.set(2);

  rkgl::mini_jacobian<double> mj;
  mj.set(2);

  auto f = [](const double *x, double *y)
  {
    y[0] = x[1];
    y[1] = -x[0];
  };

  std::vector<double> x = {1.0, 0.0};
  mj.evaluate(f, x.data(), 0.1);

  core.step(f, mj, x.data(), 1.0, 1e-14);



  return 0;
}