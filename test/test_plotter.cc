#include "CepGen/Utils/Plotter.h"
#include <iostream>
#include <cmath>

int main() {
  cepgen::utils::Graph1D graph1d;
  for (double x = -M_PI; x <= M_PI; x += 0.25)
    graph1d.addPoint(x, sin(x));
  graph1d.draw(std::cout);

  return 0;
}
