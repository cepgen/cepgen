#include "CepGen/Utils/Plotter.h"
#include <iostream>
#include <cmath>
#include <random>

int main() {
  // test 1D graph
  cepgen::utils::Graph1D graph1d;
  for (double x = -M_PI; x <= M_PI; x += 0.25)
    graph1d.addPoint(x, sin(x));
  graph1d.draw(std::cout);

  // test 2D graph
  cepgen::utils::Graph2D graph2d;
  for (double x = -5.; x < 5.; x += 0.5)
    for (double y = -5.; y < 5.; y += 0.2)
      graph2d.addPoint(x, y, (sin(x) / x) * (sin(y) / y));
  graph2d.draw(std::cout);

  std::default_random_engine gen;

  // test 1D histogram
  cepgen::utils::Hist1D hist1d(20, {-5., 5.});
  std::normal_distribution<double> gaus(0., 1.);
  for (size_t i = 0; i < 5000; ++i)
    hist1d.fill(gaus(gen));
  hist1d.draw(std::cout);

  // test 2d histogram
  cepgen::utils::Hist2D hist2d(20, {-5., 5.}, 50, {-5., 5.});
  std::normal_distribution<double> gaus1(0., 1.), gaus2(0., 1.);
  for (size_t i = 0; i < 1000; ++i)
    for (size_t j = 0; j < 1000; ++j)
      hist2d.fill(gaus1(gen), gaus2(gen));
  hist2d.draw(std::cout);

  return 0;
}
