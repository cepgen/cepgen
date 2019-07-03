#include "CepGen/IO/ExportHandler.h"

#include "CepGen/Parameters.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Event/Event.h"

#include "modules/Delphes.h"

#include <sstream>

namespace cepgen
{
  namespace output
  {
    /**
     * \brief Export handler for Delphes
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date Jul 2019
     */
    class DelphesHandler : public GenericExportHandler
    {
      public:
        /// Class constructor
        explicit DelphesHandler( const ParametersList& );

        void initialise( const Parameters& ) override;
        /// Writer operator
        void operator<<( const Event& ) override;
        void setCrossSection( double, double ) override {}

      private:
    };

    DelphesHandler::DelphesHandler( const ParametersList& params )
    {}

    void
    DelphesHandler::initialise( const Parameters& params )
    {}

    void
    DelphesHandler::operator<<( const Event& ev )
    {}
  }
}

REGISTER_IO_MODULE( delphes, DelphesHandler )
