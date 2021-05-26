#include "CepGen/Utils/Plotter.h"
#include <iostream>
#include <cmath>
#include <random>

using namespace std;

int main() {
  cout << "---------- 1D graph ----------" << endl;

  // test 1D graph
  cepgen::utils::Graph1D graph1d;
  for (double x = -M_PI; x <= M_PI; x += 0.25)
    graph1d.addPoint(x, sin(x));
  graph1d.draw(cout);

  cout << "---------- 2D graph ----------" << endl;

  // test 2D graph
  cepgen::utils::Graph2D graph2d;
  for (double x = -5.; x < 5.; x += 0.5)
    for (double y = -5.; y < 5.; y += 0.2)
      graph2d.addPoint(x, y, (sin(x) / x) * (sin(y) / y));
  graph2d.setLog(true);
  graph2d.draw(cout);

  default_random_engine gen;

  cout << "-------- 1D histogram --------" << endl;

  // test 1D histogram
  cepgen::utils::Hist1D hist1d(20, {-5., 5.});
  normal_distribution<double> gaus(0., 1.);
  for (size_t i = 0; i < 1000; ++i)
    hist1d.fill(gaus(gen));
  hist1d.setXlabel("Random variable");
  hist1d.setYlabel("Occurrences");
  hist1d.draw(cout);

  cout << "-------- 2D histogram --------" << endl;

  // test 2d histogram
  cepgen::utils::Hist2D hist2d(20, {-5., 5.}, 50, {-5., 5.});
  normal_distribution<double> gaus1(0., 1.), gaus2(0., 1.);
  for (size_t i = 0; i < 1000; ++i)
    for (size_t j = 0; j < 1000; ++j)
      hist2d.fill(gaus1(gen), gaus2(gen));
  hist2d.draw(cout);

  return 0;
}
