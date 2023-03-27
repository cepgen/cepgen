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

#include "CepGen/Cards/Handler.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Modules/CardsHandlerFactory.h"
#include "CepGen/Parameters.h"
#include "CepGen/Utils/Filesystem.h"

namespace cepgen {
  namespace card {
    Handler::Handler(const ParametersList& params)
        : NamedModule(params), filename_(steer<std::string>("filename")), rt_params_(new Parameters) {
      if (!filename_.empty())
        parseFile(filename_, rt_params_);
    }

    Parameters* Handler::parseString(const std::string& filename) {
      try {
        auto parser = CardsHandlerFactory::get().build(utils::fileExtension(filename));
        return parser->parseString(filename, new Parameters);
      } catch (const std::invalid_argument& err) {
        throw CG_FATAL("Cards:handler") << "Failed to parse the steering card at \"" << filename << "\"! "
                                        << err.what();
      }
    }

    Parameters* Handler::parseFile(const std::string& filename) {
      try {
        auto parser = CardsHandlerFactory::get().build(utils::fileExtension(filename));
        return parser->parseFile(filename, new Parameters);
      } catch (const std::invalid_argument& err) {
        throw CG_FATAL("Cards:handler") << "Failed to parse the steering card at \"" << filename << "\"! "
                                        << err.what();
      }
    }

    void Handler::write(const Parameters* params, const std::string& filename) {
      try {
        auto writer = CardsHandlerFactory::get().build(utils::fileExtension(filename));
        writer->pack(params);
        return writer->write(filename);
      } catch (const std::invalid_argument& err) {
        throw CG_FATAL("Cards:handler") << "Failed to write the configuration to \"" << filename << "\"! "
                                        << err.what();
      }
    }

    ParametersDescription Handler::description() {
      auto desc = ParametersDescription();
      desc.add<std::string>("filename", "").setDescription("Steering card to parse");
      return desc;
    }
  }  // namespace card
}  // namespace cepgen
