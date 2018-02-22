#ifndef CepGen_Tests_abort_h
#define CepGen_Tests_abort_h

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"
#include <signal.h>

void handle_ctrl_c( int signal ) { throw CepGen::Exception( __PRETTY_FUNCTION__, Form( "Process aborted with signal %d", signal ), CepGen::JustWarning ); }

class AbortHandler
{
  public:
    AbortHandler( int flags = 0 ) {
      handler_.sa_handler = handle_ctrl_c;
      sigemptyset( &handler_.sa_mask );
      handler_.sa_flags = flags;
      sigaction( SIGINT, &handler_, NULL );
    }

  private:
    struct sigaction handler_;
};

#endif
