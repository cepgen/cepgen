/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2021  Laurent Forthomme
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

#include <fstream>
#include <string>

#include "CepGen/Core/Exception.h"
#include "CepGen/Physics/MCDFileParser.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/String.h"

namespace pdg {
  const std::unordered_map<std::string, short> MCDFileParser::MAP_CHARGE_STR = {
      {"-", -3}, {"--", -6}, {"+", +3}, {"++", +6}, {"0", 0}, {"-1/3", -1}, {"-2/3", -2}, {"+1/3", +1}, {"+2/3", +2}};

  void MCDFileParser::parse(const std::string& path) {
    std::ifstream ifile(path);
    if (!ifile.is_open())
      throw CG_FATAL("MCDFileParser") << "Failed to parse MCD file \"" << path << "\"!";
    std::string line;
    while (std::getline(ifile, line)) {
      if (line[0] == '*')  // skip comments
        continue;
      std::vector<int> pdg_ids;
      std::vector<short> charges;
      double mass, width;
      std::string part_name, part_charge_int;
      {  // pdg ids
        std::istringstream ss(line.substr(PDG_BEG, PDG_END));
        std::string buf;
        // split for each PDG id
        while (ss >> buf)
          pdg_ids.emplace_back(std::stoi(buf));
      }
      {                                      // mass + error(s)
        double mass_err_low, mass_err_high;  // unused
        std::istringstream oss(line.substr(MASS_BEG, MASS_END));
        oss >> mass >> mass_err_low >> mass_err_high;
      }
      {                                        // width + error(s)
        double width_err_low, width_err_high;  // unused
        std::istringstream oss(line.substr(WIDTH_BEG, WIDTH_END));
        oss >> width >> width_err_low >> width_err_high;
      }
      {  // name + charge
        std::istringstream oss(line.substr(AUX_BEG));
        oss >> part_name >> part_charge_int;
        std::istringstream oss_ch(part_charge_int);
        std::string charge_int;
        // split by ','
        while (std::getline(oss_ch, charge_int, ',')) {
          if (MAP_CHARGE_STR.count(charge_int) == 0)
            throw CG_FATAL("MCDFileParser") << "Failed to retrieve an integer charge "
                                            << "for string \"" << charge_int << "\"!";
          charges.emplace_back(MAP_CHARGE_STR.at(charge_int));
        }
      }
      if (pdg_ids.size() != charges.size())
        throw CG_FATAL("MCDFileParser") << "Error while parsing the MCD file \"" << path << "\".\n\t"
                                        << "Invalid PDG ids / charges vectors sizes: " << pdg_ids.size()
                                        << " != " << charges.size() << ".";
      cepgen::ParticleProperties prop;
      prop.name = part_name;
      prop.descr = part_name;
      prop.colours = 1;
      prop.mass = mass;
      prop.width = width;
      prop.fermion = false;
      for (size_t i = 0; i < pdg_ids.size(); ++i) {
        prop.pdgid = (cepgen::pdgid_t)pdg_ids.at(i);
        prop.charge = charges.at(i);
        switch (pdg_ids.at(i)) {
          // start with quarks
          case 1:
          case 2:
          case 3:
          case 4:
          case 5:
          case 6:
            prop.colours = 3;
            prop.fermion = true;
            break;
          // then move to leptons/neutrinos
          case 11:
          case 12:
          case 13:
          case 14:
          case 15:
          case 16:
            prop.colours = 1;
            prop.fermion = true;
            break;
          // then gluons
          case 21:
            prop.colours = 9;
            prop.fermion = false;
            break;
          // and finally the rest
          default:
            break;
        }
        cepgen::PDG::get().define(prop);
      }
    }
    CG_INFO("MCDFileParser") << cepgen::utils::s("particle", cepgen::PDG::get().size()) << " defined from \"" << path
                             << "\". ";
  }
}  // namespace pdg
