#ifndef CepGen_Processes_TestProcess_h
#define CepGen_Processes_TestProcess_h

#include "GenericProcess.h"
#include "CepGen/Core/Functional.h"

namespace CepGen
{
  namespace Process
  {
    /// Generic process to test the Vegas instance
    class TestProcess : public GenericProcess
    {
      public:
        TestProcess();
        TestProcess( const char* formula, std::vector<std::string> args );
        ~TestProcess() {}

        void addEventContent() {}
        /// Number of dimensions on which to perform the integration
        unsigned int numDimensions( const Kinematics::ProcessMode& ) const { return num_args_; }
        /// Generic formula to compute a weight out of a point in the phase space
        double computeWeight();
        /// Dummy function to be called on events generation
        void fillKinematics( bool ) { return; }

      private:
        bool use_default_formula_;
        size_t num_args_;
        Functional<2> funct2_;
        Functional<3> funct3_;
    };
  }
}

#endif

