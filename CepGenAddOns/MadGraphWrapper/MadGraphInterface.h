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
      std::string runCommand( const std::string& ) const;
      void prepareCard() const;
      std::string prepareMadGraphProcess() const;

      std::string generateProcess() const;
      std::string generateLibrary( const std::string& ) const;

      const std::string proc_;
      const std::string model_;
      const std::string card_path_;
      const std::string tmp_dir_;
      const std::string log_filename_;
  };
}

#endif
