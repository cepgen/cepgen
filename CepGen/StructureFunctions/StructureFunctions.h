#ifndef CepGen_StructureFunctions_StructureFunctions_h
#define CepGen_StructureFunctions_StructureFunctions_h

#include "CepGen/StructureFunctions/SigmaRatio.h"

#include <iostream>
#include <map>

namespace CepGen
{
  namespace SF
  {
    /// Proton structure function to be used in the outgoing state description
    /// \note Values correspond to the LPAIR legacy steering card values
    enum struct Type {
      Invalid             = 0,
      Electron            = 1,
      ElasticProton       = 2,
      SuriYennie          = 11,
      SzczurekUleshchenko = 12,
      BlockDurandHa       = 13,
      FioreBrasse         = 101,
      ChristyBosted       = 102,
      CLAS                = 103,
      ALLM91              = 201,
      ALLM97              = 202,
      GD07p               = 203,
      GD11p               = 204,
      MSTWgrid            = 205,
      Schaefer            = 301,
      LHAPDF              = 401,
    };
  }
  std::ostream& operator<<( std::ostream&, const SF::Type& );

  class StructureFunctionsFactory;
  /// Generic placeholder for the parameterisation of nucleon structure functions
  class StructureFunctions
  {
    public:
      /// Copy constructor
      StructureFunctions( const StructureFunctions& sf ) :
        type( sf.type ), F2( sf.F2 ), FL( sf.FL ), old_vals_( sf.old_vals_ ) {}
      /// Standard SF parameterisation constructor
      StructureFunctions( const SF::Type& type = SF::Type::Invalid, double f2 = 0., double fl = 0. ) :
        type( type ), F2( f2 ), FL( fl ), old_vals_({ 0., 0. }) {}
      ~StructureFunctions() {}

      /// Human-readable description of this SF parameterisation
      friend std::ostream& operator<<( std::ostream&, const StructureFunctions& );
      /// Assign from another SF parameterisation object
      StructureFunctions& operator=( const StructureFunctions& sf ) {
        type = sf.type, F2 = sf.F2, FL = sf.FL, old_vals_ = sf.old_vals_;
        return *this;
      }

      /// Build a SF parameterisation for a given type
      static StructureFunctions builder( const SF::Type& );

      /// Compute all relevant structure functions for a given (xbj,Q²) couple
      virtual StructureFunctions& operator()( double xbj, double q2 ) { return *this; }
      /// Compute the longitudinal structure function for a given point
      virtual void computeFL( double xbj, double q2, const SF::SigmaRatio& ratio = SF::E143Ratio() );
      /// Compute the longitudinal structure function for a given point
      virtual void computeFL( double xbj, double q2, double r );
      /// Compute the F₁ structure function for a given point
      double F1( double xbj, double q2 ) const;

      /// Interpolation type of structure functions
      SF::Type type;
      double F2; ///< Last computed transverse structure function value
      double FL; ///< Last computed longitudinal structure function value

    protected:
      virtual std::string description() const; ///< Human-readable description of this SF set
      static const double mp_; ///< Proton mass, in GeV/c²
      static const double mp2_; ///< Squared proton mass, in GeV²/c⁴
      std::pair<double,double> old_vals_; ///< Last (xbj,Q²) couple computed

    private:
      std::string name_;
  };
}

#endif
