#ifndef CepGen_Processes_Test_h
#define CepGen_Processes_Test_h

#include "GenericProcess.h"

namespace CepGen
{
  namespace Process
  {
    /// Generic process to test the Vegas instance
    class TestProcess : public GenericProcess
    {
      public:
        TestProcess();
        ~TestProcess() {}

        void addEventContent() {}
        /// Number of dimensions on which to perform the integration
        unsigned int numDimensions( const Kinematics::ProcessMode& ) const { return 3; }
        /// Generic formula to compute a weight out of a point in the phase space
        double computeWeight();
        /// Dummy function to be called on events generation
        void fillKinematics( bool ) { return; }

      private:
    };
  }
}

#endif

