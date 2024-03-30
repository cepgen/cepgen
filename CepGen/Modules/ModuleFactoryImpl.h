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
  template class ModuleFactory<card::Handler>;
  template class ModuleFactory<CollinearFlux>;
  template class ModuleFactory<Coupling>;
  template class ModuleFactory<utils::Derivator>;
  template class ModuleFactory<utils::DocumentationGenerator>;
  template class ModuleFactory<utils::Drawer>;
  template class ModuleFactory<EventImporter>;
  template class ModuleFactory<EventModifier>;
  template class ModuleFactory<EventExporter>;
  template class ModuleFactory<formfac::Parameterisation>;
  template class ModuleFactory<GeneratorWorker>;
  template class ModuleFactory<Integrator>;
  template class ModuleFactory<AnalyticIntegrator>;
  template class ModuleFactory<KTFlux>;
  template class ModuleFactory<PhaseSpaceGenerator>;
  template class ModuleFactory<proc::Process>;
  template class ModuleFactory<sigrat::Parameterisation>;
  template class ModuleFactory<strfun::Parameterisation>;
  template class ModuleFactory<utils::Functional>;
  template class ModuleFactory<utils::RandomGenerator>;
}  // namespace cepgen
