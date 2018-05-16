#ifndef CepGen_Processes_TestProcess_h
#define CepGen_Processes_TestProcess_h

#include "GenericProcess.h"
#include "CepGen/Core/Functional.h"

namespace CepGen
{
  namespace Process
  {
    /// Generic process to test the Vegas instance
    template<size_t N>
    class TestProcess : public GenericProcess
    {
      public:
        TestProcess() :
          GenericProcess( "test", ".oO TEST PROCESS Oo.", false ),
          funct_( "1./(1.-cos(x*_pi)*cos(y*_pi)*cos(z*_pi))", { { "x", "y", "z" } } ) {}
        TestProcess( const char* formula, std::array<std::string,N> args ) :
          GenericProcess( "test", Form( ".oO TEST PROCESS (%s) Oo.", formula ), false ),
          funct_( formula, args ) {}

        ProcessPtr clone() const override { return ProcessPtr( new TestProcess<N>( *this ) ); }

        void addEventContent() override {}
        /// Number of dimensions on which to perform the integration
        unsigned int numDimensions( const Kinematics::Mode& ) const override { return N; }
        /// Generic formula to compute a weight out of a point in the phase space
        double computeWeight() override {
          std::array<double,N> args;
          std::copy_n( x_.begin(), N, args.begin() );
          return funct_.eval( args );
        }
        /// Dummy function to be called on events generation
        void fillKinematics( bool ) override { return; }

      private:
        Functional<N> funct_;
    };
  }
}

#endif

