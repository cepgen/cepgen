#ifndef CepGen_Tests_abort_h
#define CepGen_Tests_abort_h

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"
#include <csignal>

namespace CepGen
{
  extern volatile int gSignal;
}

struct AbortHandler
{
  AbortHandler() {
    memset( &action, 0, sizeof( struct sigaction ) );
    action.sa_sigaction = handle_ctrl_c;
    action.sa_flags = SA_SIGINFO;
    sigaction( SIGINT, &action, nullptr );
    sigaction( SIGTERM, &action, nullptr );
  }
  static void handle_ctrl_c( int signal, siginfo_t*, void* ) {
    CepGen::gSignal = signal;
  }
  struct sigaction action;
};

#endif
