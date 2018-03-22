#ifndef CepGen_Tests_abort_h
#define CepGen_Tests_abort_h

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"
#include <csignal>

namespace CepGen
{
  extern volatile int gSignal;
}

class AbortHandler
{
  public:
    AbortHandler( int flags = SA_SIGINFO ) {
      action_.sa_sigaction = handle_ctrl_c;
      sigemptyset( &action_.sa_mask );
      action_.sa_flags = flags;
      if ( sigaction( SIGINT, &action_, nullptr ) != 0
        || sigaction( SIGTERM, &action_, nullptr ) != 0 )
        throw CepGen::Exception( __PRETTY_FUNCTION__, "Failed to initialise the C-c handler!", CepGen::FatalError );
    }

  private:
    static void handle_ctrl_c( int signal, siginfo_t*, void* ) {
      CepGen::gSignal = signal;
    }
    struct sigaction action_;
};

#endif
