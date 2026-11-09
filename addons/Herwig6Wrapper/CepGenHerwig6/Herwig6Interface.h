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

namespace cepgen::herwig6 {
  static constexpr int NMXHEP = 4000;  ///< Size of the /HEPEVT/ arrays (Herwig 6.5)
  static constexpr int NMXRES = 500;   ///< Size of the particle properties arrays (Herwig 6.5)
}  // namespace cepgen::herwig6

extern "C" {
/// Standard event record
extern struct {
  int nevhep, nhep;
  int isthep[cepgen::herwig6::NMXHEP], idhep[cepgen::herwig6::NMXHEP];
  int jmohep[cepgen::herwig6::NMXHEP][2], jdahep[cepgen::herwig6::NMXHEP][2];
  double phep[cepgen::herwig6::NMXHEP][5], vhep[cepgen::herwig6::NMXHEP][4];
} hepevt_;
/// Beam particle names
extern struct {
  char part1[8], part2[8];
} hwbmch_;
/// Beams, process and number of events
extern struct {
  double ebeam1, ebeam2, pbeam1, pbeam2;
  int iproc, maxev;
} hwproc_;
/// Event-level bookkeeping
extern struct {
  double avwgt, evwgt, gamwt, tlout, wbigst, wgtmax, wgtsum, wsqsum;
  int idhw[cepgen::herwig6::NMXHEP];
  int ierror, istat, lwevt, maxer, maxpr, nowgt, nrn[2], numer, numeru, nwgts, gensof;
} hwevnt_;
/// Particle properties
extern struct {
  double rltim[cepgen::herwig6::NMXRES + 1], rmass[cepgen::herwig6::NMXRES + 1], rspin[cepgen::herwig6::NMXRES + 1];
  int ichrg[cepgen::herwig6::NMXRES + 1], idpdg[cepgen::herwig6::NMXRES + 1], iflav[cepgen::herwig6::NMXRES + 1], nres;
  int vtocdk[cepgen::herwig6::NMXRES + 1], vtordk[cepgen::herwig6::NMXRES + 1];
  int qorqqb[cepgen::herwig6::NMXRES + 1], qborqq[cepgen::herwig6::NMXRES + 1];
} hwprop_;
/// Particle names
extern struct {
  char rname[cepgen::herwig6::NMXRES + 1][8];
  char txname[cepgen::herwig6::NMXRES + 1][2][37];
} hwunam_;
/// Spin correlations steering
extern struct {
  int ndecsy, nsearch, lrdec, lwdec, syspin, threeb, fourb;
  char taudec[6];
} hwdspn_;
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
  void hwuinc();             ///< Compute the parameter-dependent constants
  void hweini();             ///< Initialise the elementary process
  void hwuine();             ///< Initialise an event
  void hwcfor();             ///< Form clusters from colour-connected partons
  void hwcdec();             ///< Decay clusters into primary hadrons
  void hwdhad();             ///< Decay unstable hadrons
  void hwdhvy();             ///< Decay heavy-flavour hadrons
  void hwufne();             ///< Finalise an event
  bool compatibleVersion();  ///< Is the Herwig 6 linked library version compatible with tested common block layouts?
}  // namespace cepgen::herwig6

#endif
