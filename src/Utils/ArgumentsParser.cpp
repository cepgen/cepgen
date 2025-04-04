/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2019-2025  Laurent Forthomme
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
#include "CepGen/Generator.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/String.h"
#include "CepGen/Version.h"

using namespace cepgen;

ArgumentsParser::ArgumentsParser(int argc, char* argv[]) : command_name_(argc > 0 ? argv[0] : "") {
  std::vector<std::string> args;
  if (argc > 1) {  // drop remove the program name
    args.resize(argc - 1);
    std::copy(argv + 1, argv + argc, args.begin());
  }
  for (auto it_arg = args.begin(); it_arg != args.end(); ++it_arg) {  // identify "arg->value" pairs in arguments
    auto arg_val = utils::split(*it_arg, '=');
    const auto argument = arg_val.at(0);
    std::string value;
    if (arg_val.size() > 1)  // "--arg=value" case
      value = arg_val.at(1);
    else if (it_arg != args.end() - 1 && (*(it_arg + 1))[0] != '-') {  // "--arg value" case
      value = *(it_arg + 1);
      ++it_arg;  // skip the next word, as it is already parsed as value for this argument one
    }
    // at this point, "arg->value" pairs are identified ;
    // parse specific commands that do not require further parsing
    if (Parameter{"help,h"}.matches(argument))  // help message is requested
      help_req_ = true;
    else if (Parameter{"version,v"}.matches(argument))  // version message is requested
      version_req_ = true;
    else if (Parameter{"debug,d"}.matches(argument)) {  // debugging is enabled
      CG_LOG_LEVEL(debug);
      if (!value.empty()) {
        if (value.find(":") == std::string::npos)
          utils::Logger::get().setOutput(new std::ofstream(value));
        else {
          const auto tokens = utils::split(value, ':');
          utils::Logger::get().setOutput(new std::ofstream(tokens.at(1)));
          for (const auto& tok : utils::split(tokens.at(0), ';'))
            utils::Logger::get().addExceptionRule(tok);
        }
      }
      debug_req_ = true;
    } else if (Parameter{"add-ons,a"}.matches(argument) && !value.empty())  // extra add-ons are requested
      add_ons_ = utils::split(value, ';');
    else if (Parameter{"cmd,c"}.matches(argument) && !value.empty())  // configuration word is added
      extra_config_.emplace_back(value);
    else  // script/executable-specific argument
      args_.emplace_back(std::make_pair(argument, value));
  }
  CG_DEBUG("ArgumentsParser") << "List of arguments retrieved: " << args_ << ". "
                              << "List of configuration words retrieved: " << extra_config_ << ".";
}

void ArgumentsParser::print_help() const { CG_LOG << help_message(); }

void ArgumentsParser::print_version() { CG_LOG << cepgen::version::banner; }

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
    std::exit(EXIT_SUCCESS);
  }
  //--- dump the generator version
  if (version_req_) {
    print_version();
    std::exit(EXIT_SUCCESS);
  }
  if (debug_req_)
    CG_DEBUG("ArgumentsParser") << "Debugging mode enabled.";
  CG_DEBUG("ArgumentsParser") << "List of add-ons to be loaded along: " << add_ons_ << ".";
  for (const auto& addon : add_ons_)
    try {
      cepgen::loadLibrary(addon);  // loading of an additional plugin into the runtime environment manager
    } catch (const cepgen::Exception& e) {
      e.dump();
    }
  size_t i = 0;
  for (auto& par : params_) {  // loop over all parameters
    if (par.name.empty()) {    // no argument name ; fetching by index
      if (i >= args_.size())
        throw CG_FATAL("ArgumentsParser") << help_message() << " Failed to retrieve required <arg" << i << ">.";
      par.value = par.boolean() ? "1" : args_.at(i).second;
    } else if (std::find_if(args_.begin(), args_.end(), [&i, &par](const auto& arg) {
                 if (!par.matches(arg.first))
                   return false;
                 if (par.boolean()) {  // all particular cases for boolean arguments
                   if (const auto word = utils::toLower(arg.second);
                       word.empty() || word == "1" || word == "on" || word == "yes" || word == "true")
                     par.value = "1";  // if the flag is set, enabled by default
                   else
                     par.value = "0";
                 } else
                   par.value = arg.second;
                 ++i;
                 return true;
               }) == args_.end()) {  // for each parameter, loop over arguments to find correspondence
      if (args_.size() > i && args_.at(i).first[0] != '-')
        par.value = args_.at(i).first;
      else if (!par.optional)  // no match
        throw CG_FATAL("ArgumentsParser")
            << help_message() << " The following parameter was not set: '" << par.name.at(0) << "'.";
    }
    par.parse();
    CG_DEBUG("ArgumentsParser") << "Parameter '" << i << "|--" << par.name.at(0)
                                << (par.name.size() > 1 ? "|-" + par.name.at(1) : "") << "'"
                                << " has value '" << par.value << "'.";
    ++i;
  }
  return *this;
}

std::string ArgumentsParser::operator[](const std::string& name) const {
  for (const auto& par : params_)
    if (par.matches(name))
      return par.value;
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
  oss << "\n   "
      << "debugging:\n\t"
      << "-d|--debug<=output file|=mod1,mod2,...:output file>\t"
      << "redirect to output file, enable module(s)";
  oss << "\n   "
      << "add-ons:\n\t"
      << "-a|--add-ons mod1<;mod2<;...>>\t"
      << "load additional module(s) into the RTE";
  oss << std::endl;
  return oss.str();
}

//----- simple parameters

ArgumentsParser::Parameter::Parameter(const std::string& param_name,
                                      const std::string& param_description,
                                      std::string* var,
                                      const std::string& def)
    : name(utils::split(param_name, ',')), description(param_description), value(def), str_variable_(var) {}

ArgumentsParser::Parameter::Parameter(const std::string& param_name,
                                      const std::string& param_description,
                                      double* var,
                                      const double& def)
    : name(utils::split(param_name, ',')),
      description(param_description),
      value(utils::format("%g", def)),
      float_variable_(var) {}

ArgumentsParser::Parameter::Parameter(const std::string& param_name,
                                      const std::string& param_description,
                                      int* var,
                                      const int& def)
    : name(utils::split(param_name, ',')),
      description(param_description),
      value(utils::format("%+i", def)),
      int_variable_(var) {}

ArgumentsParser::Parameter::Parameter(const std::string& param_name,
                                      const std::string& param_description,
                                      unsigned int* var,
                                      const unsigned int& def)
    : name(utils::split(param_name, ',')),
      description(param_description),
      value(std::to_string(def)),
      uint_variable_(var) {}

ArgumentsParser::Parameter::Parameter(const std::string& param_name,
                                      const std::string& param_description,
                                      bool* var,
                                      const bool& def)
    : name(utils::split(param_name, ',')),
      description(param_description),
      value(utils::format("%d", def)),
      bool_variable_(var) {}

ArgumentsParser::Parameter::Parameter(const std::string& param_name,
                                      const std::string& param_description,
                                      Limits* var,
                                      const Limits& def)
    : name(utils::split(param_name, ',')),
      description(param_description),
      value(utils::format("%g,%g", def.min(), def.max())),
      lim_variable_(var) {}

//----- vector of parameters

ArgumentsParser::Parameter::Parameter(const std::string& param_name,
                                      const std::string& param_description,
                                      std::vector<std::string>* var,
                                      const std::vector<std::string>& def)
    : name(utils::split(param_name, ',')),
      description(param_description),
      value(utils::merge(def, ";")),
      vec_str_variable_(var) {}

ArgumentsParser::Parameter::Parameter(const std::string& param_name,
                                      const std::string& param_description,
                                      std::vector<int>* var,
                                      const std::vector<int>& def)
    : name(utils::split(param_name, ',')),
      description(param_description),
      value(utils::merge(def, ";")),
      vec_int_variable_(var) {}

ArgumentsParser::Parameter::Parameter(const std::string& param_name,
                                      const std::string& param_description,
                                      std::vector<double>* var,
                                      const std::vector<double>& def)
    : name(utils::split(param_name, ',')),
      description(param_description),
      value(utils::merge(def, ";")),
      vec_float_variable_(var) {}

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
      *bool_variable_ = std::stoi(value) != 0;
      return *this;
    } catch (const std::invalid_argument&) {
      *bool_variable_ = (strcasecmp("true", value.c_str()) == 0 || strcasecmp("yes", value.c_str()) == 0 ||
                         strcasecmp("on", value.c_str()) == 0 || strcasecmp("1", value.c_str()) == 0) &&
                        strcasecmp("false", value.c_str()) != 0 && strcasecmp("no", value.c_str()) != 0 &&
                        strcasecmp("off", value.c_str()) != 0 && strcasecmp("0", value.c_str()) != 0;
    }
  }
  if (vec_str_variable_) {
    *vec_str_variable_ = utils::split(value, ';', true);
    return *this;
  }
  if (vec_int_variable_) {
    vec_int_variable_->clear();
    const auto buf = utils::split(value, ';');
    std::transform(buf.begin(), buf.end(), std::back_inserter(*vec_int_variable_), [](const std::string& str) {
      return std::stoi(str);
    });
    return *this;
  }
  auto unpack_floats = [](const std::string& value, char delim) {
    std::vector<double> vec_flt;
    const auto buf = utils::split(value, delim);
    std::transform(buf.begin(), buf.end(), std::back_inserter(vec_flt), [](const std::string& str) {
      try {
        return std::stod(str);
      } catch (const std::invalid_argument&) {
        return Limits::INVALID;
      }
    });
    return vec_flt;
  };
  if (vec_float_variable_) {
    *vec_float_variable_ = unpack_floats(value, ';');
    return *this;
  }
  if (lim_variable_) {
    const auto vec_flt = unpack_floats(value, ',');
    if (vec_flt.size() == 2) {
      if (vec_flt.at(0) != Limits::INVALID)
        lim_variable_->min() = vec_flt.at(0);
      if (vec_flt.at(1) != Limits::INVALID)
        lim_variable_->max() = vec_flt.at(1);
    }
    return *this;
  }
  throw CG_FATAL("Parameter") << "Failed to parse parameter \"" << name.at(0) << "\"!";
}
