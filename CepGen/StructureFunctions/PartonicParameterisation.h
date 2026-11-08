/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023-2025  Laurent Forthomme
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

#ifndef CepGen_StructureFunctions_PartonicParameterisation_h
#define CepGen_StructureFunctions_PartonicParameterisation_h

#include <array>

#include "CepGen/StructureFunctions/Parameterisation.h"

namespace cepgen::strfun {
  /// Generic partonic level perturbative structure functions built from an external PDFs grid
  class PartonicParameterisation : public Parameterisation {
  public:
    explicit PartonicParameterisation(const ParametersList&);
    ~PartonicParameterisation() override = default;

    /// Quarks types
    enum class Mode { full = 0, valence = 1, sea = 2 };
    friend std::ostream& operator<<(std::ostream&, const Mode& mode);

    static ParametersDescription description();
    void eval() final;

  protected:
    virtual double evalxQ2(int flavour, double xbj, double q2) = 0;
    const unsigned short num_flavours_;  ///< Number of quark flavours considered in the SF building
    const Mode mode_;                    ///< Quarks types considered in the SF building

    static constexpr std::array<short, 6> QUARK_PDG_IDS = {{1, 2, 3, 4, 5, 6}};
    static constexpr std::array<short, 6> Q_TIMES_3 = {{
        -1 /*d*/, 2 /*u*/, -1 /*s*/, 2 /*c*/, -1 /*b*/, 2 /*t*/
    }};
  };
}  // namespace cepgen::strfun

#endif
