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

#ifndef CepGenAddOns_BasesWrapper_BasesCommonBlocks_h
#define CepGenAddOns_BasesWrapper_BasesCommonBlocks_h

extern "C" {
void bsinit_();
void bases_(double (*fxn)(double[]), double& s, double& sigma, double& ctime, int& it1, int& it2);
void spinfo_(int&);
void spring_(double (*func)(double*), int& mxtry);

static constexpr size_t mxdim = 50;
extern struct {
  std::array<double, mxdim> xl, xu;
  int ndim, nwild;
  std::array<int, mxdim> ig;
  int ncall;
} bparm1_;
extern struct {
  double acc1, acc2;
  int itmx1, itmx2;
} bparm2_;
extern struct {
  int intv, ipnt, nloop, mloop;
} bscntl_;
extern struct {
  int mxtryp, nevent, ntrial, miss;
} sprng2_;
}

#endif
