/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2019-2024  Laurent Forthomme
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

#include <LHAPDF/LHAPDF.h>

#include "CepGen/Modules/CouplingFactory.h"
#include "CepGen/Physics/Coupling.h"

#if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION == 6
#define LHAPDF_GE_6 1
#endif

namespace cepgen {
  namespace lhapdf {
    /// A perturbative PDF-oriented \f$\alpha_S(Q^2)\f$ evaluator
    class AlphaS final : public Coupling {
    public:
      explicit AlphaS(const ParametersList& params)
          : Coupling(params)
#ifdef LHAPDF_GE_6
            ,
            lhapdf_(LHAPDF::mkPDF(steer<std::string>("pdfSet"), steer<int>("pdfMember"))) {
      }
#else
      {
        LHAPDF::initPDFSet(steer<std::string>("pdfSet"), LHAPDF::LHGRID, steer<int>("pdfMember"));
      }
#endif
      static ParametersDescription description() {
        auto desc = Coupling::description();
        desc.setDescription("LHAPDF pert.PDF-orient.evol.algo.");
        desc.add<std::string>("pdfSet", "cteq66");
        desc.add<int>("pdfMember", 0);
        return desc;
      }

      double operator()(double q) const override {
#ifdef LHAPDF_GE_6
        return lhapdf_->alphasQ(q);
#else
        return LHAPDF::alphasPDF(q);
#endif
      }

    private:
#ifdef LHAPDF_GE_6
      const std::unique_ptr<LHAPDF::PDF> lhapdf_;
#endif
    };
  }  // namespace lhapdf
}  // namespace cepgen
using AlphaSLHAPDF = cepgen::lhapdf::AlphaS;
REGISTER_ALPHAS_MODULE("lhapdf", AlphaSLHAPDF);
