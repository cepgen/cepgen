#ifndef CepGen_StructureFunctions_Schaefer_h
#define CepGen_StructureFunctions_Schaefer_h

#include "CepGen/StructureFunctions/StructureFunctions.h"
#include <memory>

namespace CepGen
{
  namespace SF
  {
    /// LUX-like hybrid modelling of \f$F_{2/L}\f$ structure functions
    class Schaefer : public StructureFunctions
    {
      public:
        /// Standard parameterisation for this SF set
        struct Parameterisation
        {
          static Parameterisation mstwGrid(); ///< "Standard" parameterisation with MSTW grid NNLO perturbative model
          static Parameterisation mstwParton(); ///< "Standard" parameterisation with partonic MSTW perturbative model
          static Parameterisation cteq(); ///< "Standard" parameterisation with partonic CTEQ perturbative model
          double q2_cut; ///< Transition \f$Q^2\f$ before reaching the continuum/perturbative regions
          double w2_lo; ///< Transition \f$W^2\f$ between resonances and hybrid continuum/resonances low-\f$Q^2\f$ regions
          double w2_hi; ///< Transition \f$W^2\f$ between hybrid continuum/resonances and continuum low-\f$Q^2\f$ regions, or continuum and perturbative high-\f$Q^2\f$ regions
          std::shared_ptr<StructureFunctions> resonances_model; ///< Resonances-dominated region (low-\f$Q^2/W^2\f$) modelling
          std::shared_ptr<StructureFunctions> perturbative_model; ///< Perturbative region (high-\f$Q^2/W^2\f$) modelling
          std::shared_ptr<StructureFunctions> continuum_model; ///< Continuum regions modelling
          bool higher_twist; ///< Enable/disable the HT correction
        };
        Schaefer( const Parameterisation& param = Parameterisation::mstwGrid() );
        Schaefer& operator()( double xbj, double q2 ) override;

        Parameterisation params;

      private:
        std::string description() const override;
        double rho( double w2 ) const;
        void initialise();
        bool initialised_;
        double inv_omega_range_;
    };
  }
}

#endif
