#ifndef CepGen_IO_GenericExportHandler_h
#define CepGen_IO_GenericExportHandler_h

#include <iosfwd>

namespace cepgen
{
  class Event;
  class Parameters;
  /// Location for all output generators
  namespace output
  {
    /**
     * \brief Output format handler for events export
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date Sep 2016
     */
    class GenericExportHandler
    {
      public:
        /// Class constructor
        explicit GenericExportHandler();
        /// Initialise the handler and its inner parameterisation
        virtual void initialise( const Parameters& ) = 0;
        /// Set the process cross section and its associated error
        virtual void setCrossSection( double /*xsec*/, double /*err_xsec*/ ) {}
        /// Set the event number
        void setEventNumber( const unsigned int& ev_id ) { event_num_ = ev_id; }
        /// Writer operator
        virtual void operator<<( const Event& ) = 0;

      protected:
        static std::string banner( const Parameters& );
        /// Event index
        unsigned int event_num_;
    };
  }
}

#endif
