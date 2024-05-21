/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2019-2022  Laurent Forthomme
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

#ifndef CepGenAddOns_PAWWrapper_PAWCommons_h
#define CepGenAddOns_PAWWrapper_PAWCommons_h

extern "C" {
#include <cfortran/cfortran.h>

#define PAWC_SIZE 50000
#define PAWC_ALIGNMENT 64

typedef struct {
  float PAW[PAWC_SIZE];
} PAWC_DEF;
#define PAWC COMMON_BLOCK(PAWC, pawc)

typedef struct {
  int iquest[100] __attribute__((aligned(PAWC_ALIGNMENT)));
} QUEST_DEF;
#define QUEST COMMON_BLOCK(QUEST, quest)
}

#endif
