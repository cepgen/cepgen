/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2025  Laurent Forthomme
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
#include "CepGen/Core/RunParameters.h"
#include "CepGen/Modules/CardsHandlerFactory.h"
#include "CepGen/Utils/Filesystem.h"

using namespace cepgen;

CardsHandlerFactory& CardsHandlerFactory::get() {
  static CardsHandlerFactory instance;
  return instance;
}

std::unique_ptr<card::Handler> CardsHandlerFactory::buildFromFilename(const std::string& filename) const {
  return build(utils::fileExtension(filename));
}
