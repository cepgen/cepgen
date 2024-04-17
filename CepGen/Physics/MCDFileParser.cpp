/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2019-2024  Laurent Forthomme
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
#include "CepGen/Utils/Filesystem.h"
#include "CepGen/Utils/String.h"

namespace pdg {
  const std::unordered_map<std::string, short> MCDFileParser::MAP_CHARGE_STR = {
      {"-", -3}, {"--", -6}, {"+", +3}, {"++", +6}, {"0", 0}, {"-1/3", -1}, {"-2/3", -2}, {"+1/3", +1}, {"+2/3", +2}};

  void MCDFileParser::parse(const std::string& path) {
    for (const auto& line : cepgen::utils::split(cepgen::utils::readFile(path), '\n')) {
      if (line[0] == '*')  // skip comments
        continue;
      std::vector<int> pdg_ids, charges;
      double mass{0.}, width{0.};
      std::string part_name;
      {  // pdg ids
        std::istringstream ss(line.substr(PDG_BEG, PDG_END));
        std::string buf;
        // split for each PDG id
        while (ss >> buf)
          pdg_ids.emplace_back(std::stoi(buf));
      }
      {                                              // mass + error(s)
        double mass_err_low{0.}, mass_err_high{0.};  // unused
        const auto mass_substr = cepgen::utils::trim(line.substr(MASS_BEG, MASS_END));
        if (!mass_substr.empty()) {
          std::istringstream oss(mass_substr);
          oss >> mass >> mass_err_low >> mass_err_high;
        }
      }
      {                                                // width + error(s)
        double width_err_low{0.}, width_err_high{0.};  // unused
        const auto width_substr = cepgen::utils::trim(line.substr(WIDTH_BEG, WIDTH_END));
        if (!width_substr.empty()) {
          std::istringstream oss(width_substr);
          oss >> width >> width_err_low >> width_err_high;
        }
      }
      {  // name + charge
        std::istringstream oss(line.substr(AUX_BEG));
        std::string part_charge_int;
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
        if (const auto ch = charges.at(i); ch != 0)
          prop.charges = {ch, -ch};
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
    CG_DEBUG("MCDFileParser") << cepgen::utils::s("particle", cepgen::PDG::get().size()) << " defined from \"" << path
                              << "\". ";
  }
}  // namespace pdg
