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

#include <LHAPDF/LHAPDF.h>

#include "CepGen/Physics/AlphaS.h"

#if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION == 6
#define LHAPDF_GE_6 1
#endif

namespace cepgen {
  /// A perturbative PDF-oriented \f$\alpha_S(Q^2)\f$ evaluator
  class AlphaSLHAPDF : public AlphaS {
  public:
    explicit AlphaSLHAPDF(const ParametersList& params)
        : AlphaS(params)
#ifdef LHAPDF_GE_6
          ,
          lhapdf_(LHAPDF::mkPDF(params.get<std::string>("pdfSet", "cteq6"), params.get<int>("pdfMember", 0))) {
    }
#else
    {
      LHAPDF::initPDFSet(params.get<std::string>("pdfSet"), LHAPDF::LHGRID, params.get<int>("pdfMember", 0));
    }
#endif
    static std::string description() { return "Perturbative PDF-oriented evolution algorithm"; }

    double operator()(double q) const override {
#ifdef LHAPDF_GE_6
      return lhapdf_->alphasQ(q);
#else
      return LHAPDF::alphasPDF(q);
#endif
    }

  private:
#ifdef LHAPDF_GE_6
    std::unique_ptr<LHAPDF::PDF> lhapdf_;
#endif
  };
}  // namespace cepgen

REGISTER_ALPHAS_MODULE("lhapdf", AlphaSLHAPDF)
