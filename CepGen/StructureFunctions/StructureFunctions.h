#ifndef CepGen_StructureFunctions_StructureFunctions_h
#define CepGen_StructureFunctions_StructureFunctions_h

#include <iostream>
#include "SigmaRatio.h"

namespace CepGen
{
  class StructureFunctions
  {
    public:
      StructureFunctions( double f2 = 0., double fl = 0. ) :
        F2( f2 ), FL( fl ) {}
      /// Proton structure function to be used in the outgoing state description
      /// \note Values correspond to the LPAIR legacy steering card values
      enum Type {
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
        Schaefer            = 301
      };

      double F2, FL;
      void computeFL( double q2, double xbj, const SF::SigmaRatio& ratio = SF::E143Ratio() );
      void computeFL( double q2, double xbj, double r );
      double F1( double q2, double xbj ) const;

    private:
      std::string name_;
  };
  std::ostream& operator<<( std::ostream&, const StructureFunctions& );
  std::ostream& operator<<( std::ostream&, const StructureFunctions::Type& );
}

#endif
