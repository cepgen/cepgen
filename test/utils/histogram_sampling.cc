/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2024-2025  Laurent Forthomme
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <cmath>
#include <random>

#include "CepGen/Generator.h"
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Modules/RandomGeneratorFactory.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Histogram.h"
#include "CepGen/Utils/RandomGenerator.h"
#include "CepGen/Utils/Test.h"

using namespace std;

int main(int argc, char* argv[]) {
  int num_samples, num_samples_ini;
  string rng_name, plotter;
  double precision;

  cepgen::initialise();

  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("num-sample,n", "number of events to sample", &num_samples, 1000)
      .addOptionalArgument("random-generator,r", "type of random number generator to use", &rng_name, "stl")
      .addOptionalArgument("plotter,p", "type of plotter to use", &plotter, "")
      .addOptionalArgument(
          "num-sample-ini", "number of events to sample for initial histogram", &num_samples_ini, 1'000'000)
      .addOptionalArgument("precision", "magnitude of precision to expect from hist compatibilities", &precision, 0.5)
      .parse();

  CG_TEST_SET_PRECISION(precision);

  std::unique_ptr<cepgen::utils::Drawer> plt;
  if (!plotter.empty())
    plt = cepgen::DrawerFactory::get().build(plotter);

  // define the CepGen random number generator
  auto rng = cepgen::RandomGeneratorFactory::get().build(rng_name);

  {  // 1D histogram testing
    // first define the histograms (original, and one to be resampled)
    auto hist = cepgen::utils::Hist1D(100, {-10., 10.}, "base", "Base"),
         hist_resampled = cepgen::utils::Hist1D(100, {-10., 10.}, "resampled", "Resampled");

    for (int i = 0; i < num_samples_ini; ++i)
      hist.fill(rng->breitWigner(0., 1.));

    // sample the original histogram
    for (int i = 0; i < num_samples; ++i)
      hist_resampled.fill(hist.sample(*rng));

    CG_TEST_EQUIV(hist_resampled.mean(), hist.mean(), "histograms mean");
    CG_TEST_EQUIV(hist_resampled.rms(), hist.rms(), "histograms rms");

    if (plt) {
      (void)plt->draw(hist);
      (void)plt->draw(hist_resampled);
      hist_resampled.scale(hist.integral() / hist_resampled.integral());
      (void)plt->draw({&hist, &hist_resampled}, "histograms_generated_resampled");
    }
  }

  {  // 2D histogram testing
    // first define the histograms (original, and one to be resampled)
    auto hist = cepgen::utils::Hist2D(100, {-10., 10.}, 100, {-10., 10.}, "base2d", "Base"),
         hist_resampled = cepgen::utils::Hist2D(100, {-10., 10.}, 100, {-10., 10.}, "resampled2d", "Resampled");

    for (int i = 0; i < num_samples_ini; ++i)
      hist.fill(rng->breitWigner(0., 1.), rng->breitWigner(0., 1.));

    // sample the original histogram
    for (int i = 0; i < num_samples; ++i)
      hist_resampled.fill(hist.sample(*rng));

    CG_TEST_EQUIV(hist_resampled.meanX(), hist.meanX(), "histograms mean X");
    CG_TEST_EQUIV(hist_resampled.rmsX(), hist.rmsX(), "histograms rms X");
    CG_TEST_EQUIV(hist_resampled.meanY(), hist.meanY(), "histograms mean Y");
    CG_TEST_EQUIV(hist_resampled.rmsY(), hist.rmsY(), "histograms rms Y");

    if (plt) {
      (void)plt->draw(hist);
      (void)plt->draw(hist_resampled);
    }
  }

  CG_TEST_SUMMARY;
}
