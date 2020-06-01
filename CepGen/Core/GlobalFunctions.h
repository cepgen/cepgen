#ifndef CepGen_Core_GlobalFunctions_h
#define CepGen_Core_GlobalFunctions_h

#include <string>

namespace cepgen
{
  /// Import a shared library in the runtime environment
  bool loadLibrary( const std::string& );
  /// Launch the initialisation procedure
  void initialise();
  /// Dump this program's header into the standard output stream
  void printHeader();
  /// List the modules registered in the runtime database
  void dumpModules();
}

#endif

