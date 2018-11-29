#ifndef CepGen_Processes_TestProcess_h
#define CepGen_Processes_TestProcess_h

#include "CepGen/Processes/GenericProcess.h"
#include "CepGen/Processes/ProcessesHandler.h"
#include "CepGen/Core/Functional.h"

namespace cepgen
{
  class ParametersList;
  namespace proc
  {
    /// Generic process to test the Vegas instance
    template<size_t N=3>
    class TestProcess : public GenericProcess
    {
      public:
        TestProcess( const ParametersList& params = ParametersList() ) :
          GenericProcess( params, "test", ".oO TEST PROCESS Oo.", false ),
          funct_( "1./(1.-cos(x*_pi)*cos(y*_pi)*cos(z*_pi))", { "x", "y", "z" } ) {}
        TestProcess( const char* formula, const std::vector<std::string>& args ) :
          GenericProcess( ParametersList(), "test", Form( ".oO TEST PROCESS (%s) Oo.", formula ), false ),
          funct_( formula, args ) {}

        ProcessPtr clone( const ParametersList& params ) const override { return ProcessPtr( new TestProcess<N>( *this ) ); }

        void addEventContent() override {}
        /// Number of dimensions on which to perform the integration
        unsigned int numDimensions() const override { return N; }
        /// Generic formula to compute a weight out of a point in the phase space
        double computeWeight() override {
          std::array<double,N> args;
          std::copy_n( x_.begin(), N, args.begin() );
          return funct_.eval( args );
        }
        /// Dummy function to be called on events generation
        void fillKinematics( bool ) override { return; }

      private:
        utils::Functional<N> funct_;
    };
    // register process and define aliases
    typedef TestProcess<1> TestProcess1D;
    typedef TestProcess<2> TestProcess2D;
    typedef TestProcess<3> TestProcess3D;
    REGISTER_PROCESS( test_1d_process, TestProcess1D )
    REGISTER_PROCESS( test_2d_process, TestProcess2D )
    REGISTER_PROCESS( test_3d_process, TestProcess3D )
  }
}

#endif
