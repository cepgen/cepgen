/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2022  Laurent Forthomme
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

#include "CepGen/Modules/EventExporterFactory.h"
#include "CepGen/OutputModules/IntegratedEventVariablesHandler.h"

namespace cepgen {
  struct TopdrawerHistsHandler final : IntegratedEventVariablesHandler {
    using IntegratedEventVariablesHandler::IntegratedEventVariablesHandler;
    static ParametersDescription description() {
      auto desc = IntegratedEventVariablesHandler::description();
      desc.setDescription("Topdrawer event histogramming tool");
      desc.add<std::string>("plotter", "topdrawer");
      return desc;
    }
  };
}  // namespace cepgen

REGISTER_IO_MODULE("topdrawer", TopdrawerHistsHandler)
