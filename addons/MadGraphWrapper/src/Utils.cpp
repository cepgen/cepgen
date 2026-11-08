/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2020-2025  Laurent Forthomme
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

#ifndef MADGRAPH_BIN
#error "*** MADGRAPH_BIN variable not set! ***"
#endif

#include <array>
#include <fstream>

#include "CepGen/Core/Exception.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/Caller.h"
#include "CepGen/Utils/Filesystem.h"
#include "CepGen/Utils/String.h"
#include "CepGenMadGraph/Utils.h"
#include "CepGenPython/Environment.h"
#include "CepGenPython/ObjectPtr.h"
#include "CepGenPython/Utils.h"

using namespace std::string_literals;

namespace cepgen::mg5amc {
  ProcessParticles unpackProcessParticles(const std::string& process_name) {
    ProcessParticles out;
    const auto process_no_removals = utils::split(utils::trim(process_name), '/').at(0);
    // dirty fix to specify incoming- and outgoing states as extracted from the mg5_aMC process string
    const auto primary_process = utils::split(process_no_removals, ',', true).at(0);
    const auto parts = utils::split(primary_process, '>', true);
    if (parts.size() != 2)
      throw CG_FATAL("MadGraphInterface:unpackProcessParticles")
          << "Unable to unpack particles from process name: \"" << process_name << "\" -> " << parts << "!";
    out.first = utils::split(parts.at(0), ' ', true);  // incoming parton-like particles
    CG_DEBUG("MadGraphInterface:unpackProcessParticles") << "Primary particles: " << out.first << ".";
    if (out.first.size() != 2)
      throw CG_FATAL("MadGraphInterface:unpackProcessParticles")
          << "Unable to unpack particles from primary particles list: \"" << parts.at(0) << "\" -> " << out.first
          << "!";
    out.second = utils::split(parts.at(1), ' ', true);  // outgoing, central system
    CG_DEBUG("MadGraphInterface:unpackProcessParticles") << "Outgoing system: " << out.second << ".";
    return out;
  }

  ParticleProperties describeParticle(const std::string& particle_name, const std::string& physics_model) {
    ParametersList plist_part;
    {  // this part retrieves the list of parameters for a given particle name, using a python call to MadGraph
      python::Environment env(ParametersList().set("name", "MadGraph5_aMC__describeParticles"s));
      const std::string name_part_dict = "part_dict";
      std::vector<std::string> cmds;
      if (!physics_model.empty()) {
        cmds.emplace_back("set auto_convert_model T");
        cmds.emplace_back("import model " + physics_model);
      }
      try {
        cmds.emplace_back("display particles " + particle_name);
        std::string py_output;
        bool found_properties{false};
        const auto tmp_path = fs::temp_directory_path() / "mg5_aMC_part_query.dat";
        if (!utils::isWriteable(tmp_path))
          throw CG_ERROR("MadGraphInterface:describeParticle")
              << "Temporary path '" << tmp_path << "' is not writeable.";
        for (auto line : runCommand(cmds, tmp_path, true)) {
          if (!found_properties) {
            if (line.find("has the following properties") != std::string::npos)
              found_properties = true;
            continue;
          }
          if (utils::startsWith(utils::trim(line), "'spin(2s+1 format)'"))  // SUPER hacky...
            line = utils::replaceAll(line,
                                     {{"(2s+1 format)"s, ""s},
                                      {/*1*/ " (scalar)"s, ""s},
                                      {/*2*/ " (fermion)"s, ""s},
                                      {/*3*/ " (vector)"s, ""s}});
          if (utils::startsWith(line, "exit"))
            break;
          py_output += line;
        }
        CG_DEBUG("MadGraphInterface:describeParticle") << "Will unpack the following attributes:\n" << py_output;
        if (py_output.empty())
          throw CG_ERROR("MadGraphInterface:describeParticle")
              << "No output retrieved from MadGraph command '" << cmds << "'. See the possible message output above.";
        if (const auto mod = python::ObjectPtr::defineModule("part", name_part_dict + "=" + py_output); mod) {
          if (const auto part_prop = mod.attribute(name_part_dict); part_prop)
            plist_part = part_prop.value<ParametersList>();
        } else
          throw CG_ERROR("MadGraphInterface:describeParticle")
              << "Error while parsing the MadGraph python output for particle '" << particle_name << "' of model '"
              << physics_model << ". Python output:\n"
              << py_output;
      } catch (const Exception&) {
        switch (particle_name[particle_name.size() - 1]) {
          case '+':
          case '-':
            throw;
          default:
            return describeParticle(particle_name + "+", physics_model);
        }
      }
    }
    // recast all the properties retrieved from the MG output into CepGen-specific particle properties
    if (const auto pdg_id = plist_part.get<int>("pdg_code", 0); pdg_id == 0)
      throw CG_FATAL("MadGraphInterface:describeParticle")
          << "Failed to retrieve a 'pdg_code' key to the unpacked particle properties: " << plist_part << ".";
    CG_DEBUG("MadGraphInterface:describeParticle") << "List of parameters retrieved from MadGraph on particle '"
                                                   << particle_name << "' from model '" << physics_model << "':\n"
                                                   << plist_part << ".";
    ParticleProperties props;
    if (auto name = plist_part.get<std::string>("name"); !name.empty()) {
      if (const auto last_chr = name[name.size() - 1]; last_chr == '-' || last_chr == '+')
        name.pop_back();
      props.name = name;
      props.human_name = name;
    }
    props.pdgid = plist_part.get<int>("pdg_code");
    plist_part.fill<int>("color"s, props.colours);  //FIXME might require some additional massaging
    props.mass = plist_part.has<double>("mass"s) ? plist_part.get<double>("mass"s) : PDG::get().mass(props.pdgid);
    props.width = plist_part.has<double>("width"s) ? plist_part.get<double>("width"s) : PDG::get().width(props.pdgid);
    if (plist_part.has<double>("charge"s)) {
      const auto ch = std::floor(plist_part.get<double>("charge"s) * 3.);
      if (ch != 0) {
        props.charges.emplace_back(ch);
        if (!plist_part.get<bool>("self_antipart"s))
          props.charges.emplace_back(-ch);
      }
    }
    props.fermion = plist_part.get<int>("spin", 0) % 2 == 0;
    CG_DEBUG("MadGraphInterface:describeParticle")
        << "Particle '" << particle_name << "' of model '" << physics_model
        << "' was successfully described from MG5 with properties: " << props << ".";
    return props;
  }

  std::vector<std::string> runCommand(const std::vector<std::string>& commands_list,
                                      const std::string& card_path,
                                      bool keep_output) {
    CG_DEBUG("MadGraphInterface:runCommand")
        << "Will run the following commands: " << commands_list << " with the following card path: " << card_path
        << ". Will keep output? " << std::boolalpha << keep_output << ".";
    std::ofstream tmp_card(card_path);
    for (const auto& command : commands_list)
      tmp_card << command << "\n";
    tmp_card << "exit\n";
    tmp_card.close();
    std::vector<std::string> output;
    const auto commands = std::vector<std::string>{MADGRAPH_BIN, "-f", card_path};
    CG_DEBUG("MadGraphInterface:runCommand")
        << "Calling mg5_aMC with the following command(s):\n\t'" << commands << "'.";
    {
      utils::Caller caller;
      for (const auto& line : utils::split(caller.call(commands), '\n'))
        if (!utils::startsWith(line, "MG5_aMC>"))  // skip the prompt lines
          output.emplace_back(line);
    }
    CG_DEBUG("MadGraphInterface:runCommand") << "\nCommands:\n"
                                             << commands_list << "\nOutput:\n"
                                             << utils::merge(output, "\n");
    if (!keep_output) {  // drop the steering card after usage
      fs::remove(card_path);
      CG_DEBUG("MadGraphInterface:runCommand") << "Steering card file '" << card_path << "' was removed.";
    }
    return output;
  }

  std::string normalise(const std::string& process_name, const std::string& physics_model) {
    return (!physics_model.empty() ? physics_model + "__" : "") +
           utils::replaceAll(process_name,
                             {{" ", "_"}, {">", "_to_"}, {"+", "p"}, {"-", "m"}, {"~", "bar"}, {"/", "_exc_"}});
  }
}  // namespace cepgen::mg5amc
