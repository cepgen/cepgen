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
  class StructureFunctions
  {
    public:
      StructureFunctions( const StructureFunctions& sf ) :
        type( sf.type ), F2( sf.F2 ), FL( sf.FL ), old_vals_( sf.old_vals_ ) {}
      StructureFunctions( const SF::Type& type = SF::Type::Invalid, double f2 = 0., double fl = 0. ) :
        type( type ), F2( f2 ), FL( fl ), old_vals_({ 0., 0. }) {}
      ~StructureFunctions() {}

      friend std::ostream& operator<<( std::ostream&, const StructureFunctions& );
      StructureFunctions& operator=( const StructureFunctions& sf ) {
        type = sf.type, F2 = sf.F2, FL = sf.FL, old_vals_ = sf.old_vals_;
        return *this;
      }

      static StructureFunctions builder( const SF::Type& );

      virtual StructureFunctions& operator()( double xbj, double q2 ) { return *this; }
      virtual void computeFL( double xbj, double q2, const SF::SigmaRatio& ratio = SF::E143Ratio() );
      virtual void computeFL( double xbj, double q2, double r );
      double F1( double xbj, double q2 ) const;

      SF::Type type;
      double F2, FL;

    protected:
      virtual std::string description() const;
      static const double mp_, mp2_;
      std::pair<double,double> old_vals_;

    private:
      std::string name_;
  };
}

#endif
