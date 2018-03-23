#ifndef CepGen_Tests_abort_h
#define CepGen_Tests_abort_h

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"
#include <csignal>

namespace CepGen
{
  extern volatile int gSignal;
  struct RunAbortedException : Exception { using Exception::Exception; };
}

class AbortHandler
{
  public:
    AbortHandler( int flags = SA_SIGINFO ) {
      action_.sa_sigaction = handle_ctrl_c;
      sigemptyset( &action_.sa_mask );
      action_.sa_flags = flags;
      init();
    }
    void setMT( bool mt_on = true ) {
      if ( mt_on )
        action_.sa_sigaction = handle_ctrl_c_mt;
      else
        action_.sa_sigaction = handle_ctrl_c;
      init();
    }

  private:
    static void handle_ctrl_c_mt( int signal, siginfo_t*, void* ) {
      CepGen::gSignal = signal;
    }
    static void handle_ctrl_c( int signal, siginfo_t*, void* ) {
      CepGen::gSignal = signal;
      throw CepGen::RunAbortedException( __PRETTY_FUNCTION__, "Run aborted...", CepGen::JustWarning );
    }
    void init() {
      if ( sigaction( SIGINT, &action_, nullptr ) != 0
        || sigaction( SIGTERM, &action_, nullptr ) != 0 )
        throw CepGen::Exception( __PRETTY_FUNCTION__, "Failed to initialise the C-c handler!", CepGen::FatalError );
    }
    struct sigaction action_;
};

#endif
