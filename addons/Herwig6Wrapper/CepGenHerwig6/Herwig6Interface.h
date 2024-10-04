/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2024  Laurent Forthomme
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

#ifndef CepGenHerwig6_Herwig6Interface_h
#define CepGenHerwig6_Herwig6Interface_h

#include <array>

extern "C" {
/// Basic parameters (and quantities derived from them)
extern struct {
  std::array<double, 16> afch[2];
  double alphem, b1lim, betaf, btclm, cafac, cffac, clmax, clpow;
  std::array<double, 2> clsmr;
  double cspeed, ensof, etamix, f0mix, f1mix, f2mix, gamh, gamw, gamz, gamzp, gev2nb, h1mix, pdiqk, psgmx;
  std::array<double, 4> pgspl;
  double phimix, pifac, prsof;
  std::array<double, 2> psplt;
  double ptrms, pxrms, qcdl3, qcdl5, qcdlam, qdiqk;
  std::array<double, 16> qfch;
  double qg, qspac, qv, scab1, swein, tmtop;
  std::array<double, 16> vfch[2];
  std::array<double, 3> vckm[3];
  double vgcut, vqcut, vpcut, zbinm, effmin, omhmix, et2mix, ph3mix, gcutme;
  int ioprem, iprint, ispac, lrsud, lwsud;
  std::array<int, 2> modpdf;
  int nbtry, ncolo, nctry, ndtry, netry, nflav, ngspl, nstru, nstry, nzbin;
  std::array<int, 2> iop4jt;
  int nprfmt, azsoft, azspin;
  std::array<int, 2> cldir;
  int hardme, nospac, prndec, prvtx, softme, zprime, prndef, prntex, prnweb;
} hwpram_;
}

/// Herwig 6 utilities namespace
namespace cepgen::herwig6 {
  void initialise();
  double hwuaem(double q2);
  double hwualf(int mode, double q2);
  double hwsfun(double xbj, double q2, int idhad, int nset, int ibeam);
}  // namespace cepgen::herwig6

#endif
