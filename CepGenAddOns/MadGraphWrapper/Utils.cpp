/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2020-2024  Laurent Forthomme
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
#include "CepGenAddOns/MadGraphWrapper/Utils.h"
#include "CepGenAddOns/PythonWrapper/Environment.h"
#include "CepGenAddOns/PythonWrapper/ObjectPtr.h"
#include "CepGenAddOns/PythonWrapper/PythonUtils.h"

namespace cepgen {
  namespace mg5amc {
    ProcessParticles unpackProcessParticles(const std::string& proc) {
      ProcessParticles out;
      auto trim_all = [](std::vector<std::string> coll) -> std::vector<std::string> {
        std::for_each(coll.begin(), coll.end(), [](std::string& it) { it = utils::trim(it); });
        return coll;
      };
      //--- dirty fix to specify incoming- and outgoing states
      //    as extracted from the mg5_aMC process string
      const auto prim_proc = utils::split(utils::trim(proc), ',')[0];
      auto parts = trim_all(utils::split(prim_proc, '>'));
      if (parts.size() != 2)
        throw CG_FATAL("MadGraphInterface:unpackProcessParticles")
            << "Unable to unpack particles from process name: \"" << proc << "\" -> " << parts << "!";
      //--- incoming parton-like particles
      auto prim_parts = trim_all(utils::split(parts[0], ' '));
      CG_DEBUG("MadGraphInterface:unpackProcessParticles") << "Primary particles: " << prim_parts;
      if (prim_parts.size() != 2)
        throw CG_FATAL("MadGraphInterface:unpackProcessParticles")
            << "Unable to unpack particles from primary particles list: \"" << parts[0] << "\" -> " << prim_parts
            << "!";
      for (const auto& p : prim_parts)
        out.first.emplace_back(p);
      //---- outgoing system
      auto dec_parts = trim_all(utils::split(parts[1], ' '));
      CG_DEBUG("MadGraphInterface:unpackProcessParticles") << "Outgoing system: " << dec_parts;
      for (auto& p : dec_parts)
        out.second.emplace_back(p);
      return out;
    }

    ParticleProperties describeParticle(const std::string& part_name, const std::string& model) {
      ParametersList plist_part;
      {  // this part retrieves the list of parameters for a given particle name, using a python call to MadGraph
        python::Environment env({});
        const std::string name_part_dict = "part_dict";
        std::vector<std::string> cmds;
        if (!model.empty()) {
          cmds.emplace_back("set auto_convert_model T");
          cmds.emplace_back("import model " + model);
        }
        try {
          cmds.emplace_back("display particles " + part_name);
          std::string py_output;
          bool found_properties{false};
          for (const auto& line : runCommand(cmds, "/tmp/mg5_aMC_part_query.dat", true)) {
            if (!found_properties) {
              if (line.find("has the following properties") != std::string::npos)
                found_properties = true;
              continue;
            }
            if (utils::startsWith(line, "exit"))
              break;
            py_output += line;
          }
          if (py_output.empty())
            throw CG_ERROR("MadGraphInterface:describeParticle")
                << "No output retrieved from MadGraph command '" << cmds << "'. See the possible message output above.";
          if (auto mod = python::ObjectPtr::defineModule("part", name_part_dict + "=" + py_output); mod) {
            if (auto part_prop = mod.attribute(name_part_dict); part_prop)
              plist_part = part_prop.value<ParametersList>();
          } else
            throw CG_ERROR("MadGraphInterface:describeParticle")
                << "Error while parsing the MadGraph python output for particle '" << part_name << "' of model '"
                << model << ". Python output:\n"
                << py_output;
        } catch (const Exception& exc) {
          switch (part_name[part_name.size() - 1]) {
            case '+':
            case '-':
              throw;
            default:
              return describeParticle(part_name + "+", model);
          }
        }
      }
      // recast all the properties retrieved from the MG output into CepGen-specific particle properties
      const auto pdg_id = plist_part.get<int>("pdg_code", 0);
      if (pdg_id == 0)
        throw CG_FATAL("MadGraphInterface:describeParticle")
            << "Failed to retrieve a 'pdg_code' key to the unpacked particle properties: " << plist_part << ".";
      CG_DEBUG("MadGraphInterface:describeParticle") << "List of parameters retrieved from MadGraph on particle '"
                                                     << part_name << "' from model '" << model << "':\n"
                                                     << plist_part << ".";
      ParticleProperties props;
      if (auto name = plist_part.get<std::string>("name"); !name.empty()) {
        if (const auto last_chr = name[name.size() - 1]; last_chr == '-' || last_chr == '+')
          name.pop_back();
        props.name = name;
        props.descr = name;
      }
      props.pdgid = plist_part.get<int>("pdg_code");
      plist_part.fill<int>("color", props.colours);  //FIXME might not be correct
      props.mass = plist_part.has<double>("mass") ? plist_part.get<double>("mass") : PDG::get().mass(props.pdgid);
      props.width = plist_part.has<double>("width") ? plist_part.get<double>("width") : PDG::get().width(props.pdgid);
      if (plist_part.has<double>("charge")) {
        const auto ch = std::floor(plist_part.get<double>("charge") * 3.);
        if (ch != 0) {
          props.charges.emplace_back(ch);
          if (!plist_part.get<bool>("self_antipart"))
            props.charges.emplace_back(-ch);
        }
      }
      props.fermion = plist_part.get<int>("spin", 0) % 2 == 0;
      CG_DEBUG("MadGraphInterface:describeParticle")
          << "Particle '" << part_name << "' of model '" << model
          << "' was successfully described from MG5 with properties: " << props << ".";
      return props;
    }

    std::vector<std::string> runCommand(const std::vector<std::string>& cmds,
                                        const std::string& card_path,
                                        bool keep_output) {
      std::ofstream tmp_card(card_path);
      for (const auto& cmd : cmds)
        tmp_card << cmd << "\n";
      tmp_card << "exit\n";
      tmp_card.close();
      std::vector<std::string> output;
      {
        utils::Caller caller;
        for (const auto& line : utils::split(caller.call({MADGRAPH_BIN, "-f", card_path}), '\n'))
          if (!utils::startsWith(line, "MG5_aMC>"))
            output.emplace_back(line);
      }
      CG_DEBUG("MadGraphInterface:runCommand") << "\nCommands:\n"
                                               << cmds << "\nOutput:\n"
                                               << utils::merge(output, "\n");
      if (!keep_output) {
        fs::remove(card_path);
        CG_DEBUG("MadGraphInterface:runCommand") << "Steering card file '" << card_path << "' was removed.";
      }
      return output;
    }

    std::string normalise(const std::string& proc_name, const std::string& model) {
      return (!model.empty() ? model + "__" : "") +
             utils::replaceAll(proc_name, {{" ", "_"}, {">", "_to_"}, {"+", "p"}, {"-", "m"}, {"~", "bar"}});
    }
  }  // namespace mg5amc
}  // namespace cepgen
