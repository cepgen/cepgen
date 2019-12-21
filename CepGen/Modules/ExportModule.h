#ifndef CepGen_Modules_ExportModule_h
#define CepGen_Modules_ExportModule_h

#include "CepGen/Core/ParametersList.h"

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
    class ExportModule
    {
      public:
        /// Class constructor
        explicit ExportModule( const ParametersList& );
        virtual ~ExportModule() = default;
        const std::string& name() const { return name_; }
        /// Initialise the handler and its inner parameterisation
        virtual void initialise( const Parameters& ) = 0;
        /// Global list of steering parameters
        const ParametersList& parameters() { return params_; }
        /// Set the process cross section and its associated error
        virtual void setCrossSection( double /*xsec*/, double /*err_xsec*/ ) {}
        /// Set the event number
        void setEventNumber( const unsigned int& ev_id ) { event_num_ = ev_id; }
        /// Writer operator
        virtual void operator<<( const Event& ) = 0;

      protected:
        static std::string banner( const Parameters&, const std::string& prep = "" );
        /// Set of parameters to steer this output module
        const ParametersList params_;
        /// Module unique name
        const std::string name_;
        /// Event index
        unsigned long long event_num_;
    };
  }
}

#endif

