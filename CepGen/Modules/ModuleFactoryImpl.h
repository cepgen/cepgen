/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2021-2024  Laurent Forthomme
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
#include "CepGen/CollinearFluxes/CollinearFlux.h"
#include "CepGen/Core/GeneratorWorker.h"
#include "CepGen/Event/Event.h"
#include "CepGen/EventFilter/EventExporter.h"
#include "CepGen/EventFilter/EventImporter.h"
#include "CepGen/EventFilter/EventModifier.h"
#include "CepGen/FormFactors/Parameterisation.h"
#include "CepGen/Integration/AnalyticIntegrator.h"
#include "CepGen/Integration/Integrator.h"
#include "CepGen/KTFluxes/KTFlux.h"
#include "CepGen/Physics/Coupling.h"
#include "CepGen/Process/PhaseSpaceGenerator.h"
#include "CepGen/Process/Process.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/StructureFunctions/SigmaRatio.h"
#include "CepGen/Utils/Derivator.h"
#include "CepGen/Utils/DocumentationGenerator.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Functional.h"
#include "CepGen/Utils/RandomGenerator.h"

namespace cepgen {
  template class ModuleFactory<card::Handler, std::string>;
  template class ModuleFactory<CollinearFlux, std::string>;
  template class ModuleFactory<Coupling, std::string>;
  template class ModuleFactory<utils::Derivator, std::string>;
  template class ModuleFactory<utils::DocumentationGenerator, std::string>;
  template class ModuleFactory<utils::Drawer, std::string>;
  template class ModuleFactory<EventImporter, std::string>;
  template class ModuleFactory<EventModifier, std::string>;
  template class ModuleFactory<EventExporter, std::string>;
  template class ModuleFactory<formfac::Parameterisation, std::string>;
  template class ModuleFactory<GeneratorWorker, std::string>;
  template class ModuleFactory<Integrator, std::string>;
  template class ModuleFactory<AnalyticIntegrator, std::string>;
  template class ModuleFactory<KTFlux, std::string>;
  template class ModuleFactory<PhaseSpaceGenerator, std::string>;
  template class ModuleFactory<proc::Process, std::string>;
  template class ModuleFactory<sigrat::Parameterisation, std::string>;
  template class ModuleFactory<strfun::Parameterisation, std::string>;
  template class ModuleFactory<utils::Functional, std::string>;
  template class ModuleFactory<utils::RandomGenerator, std::string>;
}  // namespace cepgen
