#include "CepGen/Utils/Plotter.h"
#include <iostream>
#include <cmath>

int main() {
  cepgen::utils::Graph1D graph1d;
  for (double x = -M_PI; x <= M_PI; x += 0.25)
    graph1d.addPoint(x, sin(x));
  graph1d.draw(std::cout);

  cepgen::utils::Graph2D graph2d;
  for (double x = -5.; x < 5.; x += 0.5)
    for (double y = -5.; y < 5.; y += 0.5)
      graph2d.addPoint(x, y, (sin(x) / x) * (sin(y) / y));
  graph2d.draw(std::cout);
  return 0;
}
