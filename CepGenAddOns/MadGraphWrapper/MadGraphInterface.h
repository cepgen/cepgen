#ifndef CepGenAddOns_MadGraphWrapper_MadGraphInterface_h
#define CepGenAddOns_MadGraphWrapper_MadGraphInterface_h

#include "CepGen/Core/ParametersList.h"

#include <string>

namespace cepgen
{
  class MadGraphInterface
  {
    public:
      MadGraphInterface( const ParametersList& );

      std::string run() const;

    private:
      static std::string runCommand( const std::string& );
      static std::string generateLibrary( const std::string&, const std::string&, const std::string& );
      static std::string generateProcess( const std::string& );

      void prepareCard() const;
      std::string prepareMadGraphProcess() const;

      const std::string proc_;
      const std::string model_;
      const std::string card_path_;
      const std::string standalone_cpp_path_;
      const std::string tmp_dir_;
      const std::string log_filename_;
  };
}

#endif
