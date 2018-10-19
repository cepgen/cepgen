#ifndef CepGen_StructureFunctions_Schaefer_h
#define CepGen_StructureFunctions_Schaefer_h

#include "CepGen/StructureFunctions/StructureFunctions.h"

#include <vector>
#include <memory>

namespace cepgen
{
  namespace strfun
  {
    /// LUX-like hybrid modelling of \f$F_{2,L}\f$ structure functions
    class Schaefer : public Parameterisation
    {
      public:
        Schaefer();
        explicit Schaefer( const ParametersList& params );
        Schaefer& operator()( double xbj, double q2 ) override;

      private:
        std::string description() const override;
        double rho( double w2 ) const;
        void initialise();
        /// Transition \f$Q^2\f$ before reaching the continuum/perturbative regions
        double q2_cut_;
        /// Transition \f$W^2\f$ between:
        /// - resonances and hybrid continuum/resonances low-\f$Q^2\f$ regions,
        /// - hybrid continuum/resonances and continuum low-\f$Q^2\f$ regions, or
        /// - continuum and perturbative high-\f$Q^2\f$ regions
        std::vector<double> w2_lim_;
        /// Enable/disable the HT correction
        bool higher_twist_;
        /// Resonances-dominated region (low-\f$Q^2/W^2\f$) modelling
        std::shared_ptr<Parameterisation> resonances_model_;
        /// Perturbative region (high-\f$Q^2/W^2\f$) modelling
        std::shared_ptr<Parameterisation> perturbative_model_;
        /// Continuum regions modelling
        std::shared_ptr<Parameterisation> continuum_model_;
        bool initialised_;
        double inv_omega_range_;
    };
  }
}

#endif
