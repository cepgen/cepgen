#ifndef CepGen_Modules_ExportModule_h
#define CepGen_Modules_ExportModule_h

#include "CepGen/Modules/NamedModule.h"

#include <iosfwd>
#include <string>

namespace cepgen
{
  class Event;
  class Parameters;
  /// Location for all output generators
  namespace io
  {
    /**
     * \brief Output format handler for events export
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date Sep 2016
     */
    class ExportModule : public NamedModule<std::string>
    {
      public:
        /// Class constructor
        /// \param[in] params User-controlled steering parameters for this module
        explicit ExportModule( const ParametersList& params );
        virtual ~ExportModule();

        /// Unique name of the output module
        const std::string& name() const { return name_; }
        /// Set the process cross section and its associated error
        virtual void setCrossSection( double /*xsec*/, double /*err_xsec*/ ) {}
        /// Set the event number
        void setEventNumber( const unsigned int& ev_id ) { event_num_ = ev_id; }

        /// Initialise the handler and its inner parameterisation
        virtual void initialise( const Parameters& ) = 0;
        /// Writer operator
        virtual void operator<<( const Event& ) = 0;

      protected:
        /// Print a banner containing all runtime parameters information
        static std::string banner( const Parameters&, const std::string& prep = "" );
        /// Event index
        unsigned long long event_num_;
    };
  }
}

#endif

