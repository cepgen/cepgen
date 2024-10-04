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

#include "CepGen/Core/Exception.h"
#include "CepGen/Generator.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Test.h"
#include "CepGenMadGraph/Utils.h"

using namespace std;
using namespace std::string_literals;

int main(int argc, char* argv[]) {
  cepgen::ArgumentsParser(argc, argv).parse();
  cepgen::initialise();

  for (const auto& part : vector<pair<cepgen::spdgid_t, string> >{
           {11, "e"},    {12, "ve"},   {13, "mu"},   {14, "vm"},   {15, "ta"}, {16, "vt"}, {-11, "e+"}, {-12, "ve~"},
           {-13, "mu+"}, {-14, "vm~"}, {-15, "ta+"}, {-16, "vt~"}, {22, "a"},  {23, "z"},  {-24, "w-"}, {24, "w+"},
           {25, "h"},    {1, "d"},     {2, "u"},     {3, "s"},     {4, "c"},   {5, "b"},   {6, "t"},    {-1, "d~"},
           {-2, "u~"},   {-3, "s~"},   {-4, "c~"},   {-5, "b~"},   {-6, "t~"}}) {
    const auto mg_prop = cepgen::mg5amc::describeParticle(part.second, "sm");
    const auto cg_prop = cepgen::PDG::get()(part.first);
    const auto name = part.second + "/" + cepgen::PDG::get().name(part.first);
    CG_TEST_EQUAL(mg_prop.pdgid, std::labs(part.first), name + " PDG");
    CG_TEST_EQUAL(mg_prop.fermion, cg_prop.fermion, name + " fermion/boson");
    CG_TEST_EQUIV(mg_prop.mass, cg_prop.mass, name + " mass");
    CG_TEST_EQUIV(mg_prop.width, cg_prop.width, name + " width");
    CG_TEST_EQUAL(mg_prop.charges, cg_prop.charges, name + " charges");
  }

  CG_TEST_SUMMARY;
}
