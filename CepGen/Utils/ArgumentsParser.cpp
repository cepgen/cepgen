/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2019-2023  Laurent Forthomme
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

#include <strings.h>

#include <algorithm>
#include <fstream>
#include <sstream>

#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/String.h"
#include "CepGen/Version.h"

namespace cepgen {
  ArgumentsParser::ArgumentsParser(int argc, char* argv[])
      : command_name_(argc > 0 ? argv[0] : ""),
        help_str_({{"help,h"}}),
        version_str_({{"version,v"}}),
        config_str_({{"cmd,c"}}),
        debug_str_({{"debug,d"}}) {
    //--- first remove the program name
    std::vector<std::string> args_tmp;
    if (argc > 1) {
      args_tmp.resize(argc - 1);
      std::copy(argv + 1, argv + argc, args_tmp.begin());
    }
    //--- then loop on user arguments to identify word -> value pairs
    for (auto it_arg = args_tmp.begin(); it_arg != args_tmp.end(); ++it_arg) {
      auto arg_val = utils::split(*it_arg, '=');  // particular case for --arg=value
      //--- check if help message is requested
      for (const auto& str : help_str_)
        if (arg_val.at(0) == "--" + str.name.at(0) || (str.name.size() > 1 && arg_val.at(0) == "-" + str.name.at(1)))
          help_req_ = true;
      //--- check if version message is requested
      for (const auto& str : version_str_)
        if (arg_val.at(0) == "--" + str.name.at(0) || (str.name.size() > 1 && arg_val.at(0) == "-" + str.name.at(1)))
          version_req_ = true;
      //--- check if debugging is requested
      for (const auto& str : debug_str_)
        if (arg_val.at(0) == "--" + str.name.at(0) || (str.name.size() > 1 && arg_val.at(0) == "-" + str.name.at(1))) {
          CG_LOG_LEVEL(debug);
          if (arg_val.size() > 1)
            utils::Logger::get().output().reset(new std::ofstream(arg_val.at(1)));
          debug_req_ = true;
        }
      //--- check if configuration word is requested
      auto it = std::find_if(config_str_.begin(), config_str_.end(), [&arg_val](const auto& str) {
        return (arg_val.at(0) == "--" + str.name.at(0) ||
                (str.name.size() > 1 && arg_val.at(0) == "-" + str.name.at(1)));
      });
      if (it != config_str_.end()) {
        // if a configuration word is found, all the remaining flags are parsed as such
        extra_config_ = std::vector<std::string>(it_arg + 1, args_tmp.end());
        break;
      }
      //--- parse arguments if word found after
      if (arg_val.size() == 1 && arg_val.at(0)[0] == '-' && it_arg != std::prev(args_tmp.end())) {
        const auto& word = *std::next(it_arg);
        if (word[0] != '-') {
          arg_val.emplace_back(*std::next(it_arg));
          ++it_arg;
        }
      }
      args_.emplace_back(std::make_pair(arg_val.at(0), arg_val.size() > 1 ? arg_val.at(1) : ""));
    }
  }

  void ArgumentsParser::print_help() const { CG_LOG << help_message(); }

  void ArgumentsParser::print_version() const { CG_LOG << cepgen::version::banner; }

  void ArgumentsParser::dump() const {
    CG_INFO("ArgumentsParser").log([&](auto& info) {
      info << "List of parameters retrieved from command-line:";
      for (const auto& par : params_)
        info << "\n\t[--" << par.name.at(0) << (par.name.size() > 1 ? "|-" + par.name.at(1) : "")
             << (par.optional ? ", optional" : "") << "] = " << par.value;
    });
  }

  ArgumentsParser& ArgumentsParser::parse() {
    if (help_req_) {
      print_help();
      exit(0);
    }
    //--- dump the generator version
    if (version_req_) {
      print_version();
      exit(0);
    }
    if (debug_req_)
      CG_DEBUG("ArgumentsParser") << "Debugging mode enabled.";
    //--- loop over all parameters
    size_t i = 0;
    for (auto& par : params_) {
      if (par.name.empty()) {
        //--- no argument name ; fetching by index
        if (i >= args_.size())
          throw CG_FATAL("ArgumentsParser") << help_message() << " Failed to retrieve required <arg" << i << ">.";
        par.value = !par.boolean() ? args_.at(i).second : "1";
      } else {
        // for each parameter, loop over arguments to find correspondence
        auto it = std::find_if(args_.begin(), args_.end(), [&i, &par](const auto& arg) {
          if (arg.first != "--" + par.name.at(0) && (par.name.size() < 2 || arg.first != "-" + par.name.at(1)))
            return false;
          par.value = arg.second;
          if (par.boolean()) {  // all particular cases for boolean arguments
            const auto word = utils::tolower(arg.second);
            if (word.empty() || word == "1" || word == "on" || word == "yes" || word == "true")
              par.value = "1";  // if the flag is set, enabled by default
            else
              par.value = "0";
          }
          ++i;
          return true;
        });
        if (it == args_.end()) {
          if (args_.size() > i && args_.at(i).first[0] != '-')
            par.value = args_.at(i).first;
          else if (!par.optional)  // no match
            throw CG_FATAL("ArgumentsParser")
                << help_message() << " The following parameter was not set: '" << par.name.at(0) << "'.";
        }
      }
      par.parse();
      CG_DEBUG("ArgumentsParser") << "Parameter '" << i << "|--" << par.name.at(0)
                                  << (par.name.size() > 1 ? "|-" + par.name.at(1) : "") << "'"
                                  << " has value '" << par.value << "'.";
      ++i;
    }
    return *this;
  }

  std::string ArgumentsParser::operator[](std::string name) const {
    for (const auto& par : params_) {
      if ("--" + par.name.at(0) == name)
        return par.value;
      if (par.name.size() > 1 && "-" + par.name.at(1) == name)
        return par.value;
    }
    throw CG_FATAL("ArgumentsParser") << "The parameter \"" << name << "\" was not declared "
                                      << "in the arguments parser constructor!";
  }

  std::string ArgumentsParser::help_message() const {
    std::ostringstream oss;
    std::vector<std::pair<Parameter, size_t> > req_params, opt_params;
    oss << "Usage: " << command_name_;
    size_t i = 0;
    for (const auto& par : params_) {
      if (par.optional) {
        opt_params.emplace_back(std::make_pair(par, i));
        oss << " [";
      } else {
        req_params.emplace_back(std::make_pair(par, i));
        oss << " ";
      }
      oss << (!par.name.at(0).empty() ? "--" : " <arg" + std::to_string(i) + ">") << par.name.at(0);
      if (par.name.size() > 1)
        oss << (!par.name.at(0).empty() ? "|" : "") << "-" << par.name.at(1);
      if (par.optional)
        oss << "]";
      ++i;
    }
    if (req_params.size() > 0) {
      oss << "\n    " << utils::s("required argument", req_params.size(), false) << ":";
      for (const auto& par : req_params)
        oss << utils::format(
            "\n\t%s%-18s\t%-30s",
            (par.first.name.size() > 1 ? "-" + par.first.name.at(1) + "|" : "").c_str(),
            (!par.first.name.at(0).empty() ? "--" + par.first.name.at(0) : "<arg" + std::to_string(par.second) + ">")
                .c_str(),
            par.first.description.c_str());
    }
    if (opt_params.size() > 0) {
      oss << "\n    " << utils::s("optional argument", opt_params.size(), false) << ":";
      for (const auto& par : opt_params)
        oss << utils::format(
            "\n\t%s%-18s\t%-30s\tdef: '%s'",
            (par.first.name.size() > 1 ? "-" + par.first.name.at(1) + "|" : "").c_str(),
            (!par.first.name.at(0).empty() ? "--" + par.first.name.at(0) : "<arg" + std::to_string(par.second) + ">")
                .c_str(),
            par.first.description.c_str(),
            par.first.value.c_str());
    }
    oss << std::endl;
    return oss.str();
  }

  //----- simple parameters

  ArgumentsParser::Parameter::Parameter(std::string pname, std::string pdesc, std::string* var, std::string def)
      : name(utils::split(pname, ',')), description(pdesc), value(def), str_variable_(var) {}

  ArgumentsParser::Parameter::Parameter(std::string pname, std::string pdesc, double* var, double def)
      : name(utils::split(pname, ',')), description(pdesc), value(utils::format("%g", def)), float_variable_(var) {}

  ArgumentsParser::Parameter::Parameter(std::string pname, std::string pdesc, int* var, int def)
      : name(utils::split(pname, ',')), description(pdesc), value(utils::format("%+i", def)), int_variable_(var) {}

  ArgumentsParser::Parameter::Parameter(std::string pname, std::string pdesc, unsigned int* var, unsigned int def)
      : name(utils::split(pname, ',')), description(pdesc), value(std::to_string(def)), uint_variable_(var) {}

  ArgumentsParser::Parameter::Parameter(std::string pname, std::string pdesc, bool* var, bool def)
      : name(utils::split(pname, ',')), description(pdesc), value(utils::format("%d", def)), bool_variable_(var) {}

  ArgumentsParser::Parameter::Parameter(std::string pname, std::string pdesc, Limits* var, Limits def)
      : name(utils::split(pname, ',')),
        description(pdesc),
        value(utils::format("%g,%g", def.min(), def.max())),
        lim_variable_(var) {}

  //----- vector of parameters

  ArgumentsParser::Parameter::Parameter(std::string pname,
                                        std::string pdesc,
                                        std::vector<std::string>* var,
                                        std::vector<std::string> def)
      : name(utils::split(pname, ',')), description(pdesc), value(utils::merge(def, ",")), vec_str_variable_(var) {}

  ArgumentsParser::Parameter::Parameter(std::string pname,
                                        std::string pdesc,
                                        std::vector<int>* var,
                                        std::vector<int> def)
      : name(utils::split(pname, ',')), description(pdesc), value(utils::merge(def, ",")), vec_int_variable_(var) {}

  ArgumentsParser::Parameter::Parameter(std::string pname,
                                        std::string pdesc,
                                        std::vector<double>* var,
                                        std::vector<double> def)
      : name(utils::split(pname, ',')), description(pdesc), value(utils::merge(def, ",")), vec_float_variable_(var) {}

  bool ArgumentsParser::Parameter::matches(const std::string& key) const {
    if (key == "--" + name.at(0))
      return true;
    if (name.size() > 1 && key == "-" + name.at(1))
      return true;
    return false;
  }

  ArgumentsParser::Parameter& ArgumentsParser::Parameter::parse() {
    CG_DEBUG("ArgumentsParser:Parameter:parse") << "Parsing argument " << name << ".";
    if (str_variable_) {
      *str_variable_ = value;
      return *this;
    }
    if (float_variable_)
      try {
        *float_variable_ = std::stod(value);
        return *this;
      } catch (const std::invalid_argument&) {
        throw CG_FATAL("ArgumentsParser:Parameter:parse") << "Failed to parse variable '" << name << "' as float!";
      }
    if (int_variable_)
      try {
        *int_variable_ = std::stoi(value);
        return *this;
      } catch (const std::invalid_argument&) {
        throw CG_FATAL("ArgumentsParser:Parameter:parse") << "Failed to parse variable '" << name << "' as integer!";
      }
    if (uint_variable_)
      try {
        *uint_variable_ = std::stoi(value);
        return *this;
      } catch (const std::invalid_argument&) {
        throw CG_FATAL("ArgumentsParser:Parameter:parse")
            << "Failed to parse variable '" << name << "' as unsigned integer!";
      }
    if (bool_variable_) {
      try {
        *bool_variable_ = (std::stoi(value) != 0);
        return *this;
      } catch (const std::invalid_argument&) {
        *bool_variable_ = (strcasecmp("true", value.c_str()) == 0 || strcasecmp("yes", value.c_str()) == 0 ||
                           strcasecmp("on", value.c_str()) == 0 || strcasecmp("1", value.c_str()) == 0) &&
                          strcasecmp("false", value.c_str()) != 0 && strcasecmp("no", value.c_str()) != 0 &&
                          strcasecmp("off", value.c_str()) != 0 && strcasecmp("0", value.c_str()) != 0;
      }
    }
    if (vec_str_variable_) {
      *vec_str_variable_ = utils::split(value, ',');
      return *this;
    }
    if (vec_int_variable_) {
      vec_int_variable_->clear();
      const auto buf = utils::split(value, ',');
      std::transform(buf.begin(), buf.end(), std::back_inserter(*vec_int_variable_), [](const std::string& str) {
        return std::stoi(str);
      });
      return *this;
    }
    if (vec_float_variable_ || lim_variable_) {
      std::vector<double> vec_flt;
      const auto buf = utils::split(value, ',');
      std::transform(buf.begin(), buf.end(), std::back_inserter(vec_flt), [](const std::string& str) {
        try {
          return std::stod(str);
        } catch (const std::invalid_argument&) {
          return Limits::INVALID;
        }
      });
      if (vec_float_variable_)
        *vec_float_variable_ = vec_flt;
      else if (vec_flt.size() == 2) {
        if (vec_flt.at(0) != Limits::INVALID)
          lim_variable_->min() = vec_flt.at(0);
        if (vec_flt.at(1) != Limits::INVALID)
          lim_variable_->max() = vec_flt.at(1);
      }
      return *this;
    }
    throw CG_FATAL("Parameter") << "Failed to parse parameter \"" << name.at(0) << "\"!";
  }
}  // namespace cepgen
