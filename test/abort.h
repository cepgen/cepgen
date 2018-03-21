#ifndef CepGen_Tests_abort_h
#define CepGen_Tests_abort_h

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"
#include <csignal>

namespace CepGen
{
  struct AbortException : Exception { using Exception::Exception; };
}

void handle_ctrl_c( int signal ) {
  throw CepGen::AbortException( __PRETTY_FUNCTION__, Form( "Process aborted with signal %d", signal ), CepGen::JustWarning );
}

struct AbortHandler
{
  AbortHandler() {
    std::signal( SIGINT, handle_ctrl_c );
    std::signal( SIGTERM, handle_ctrl_c );
  }
};

#endif
