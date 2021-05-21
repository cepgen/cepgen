#include "CepGen/Utils/Plotter.h"
#include <iostream>
#include <cmath>

int main() {
  cepgen::utils::Graph1D graph;
  for (double x = -M_PI; x <= M_PI; x += 0.2)
    graph.addPoint(x, sin(x));
  graph.draw(std::cout);
  return 0;
}
