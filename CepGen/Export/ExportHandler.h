#ifndef CepGen_Export_ExportHandler_h
#define CepGen_Export_ExportHandler_h

#include "CepGen/Parameters.h"
#include "CepGen/Physics/Event.h"

namespace CepGen
{
  /// Location for all output generators
  namespace OutputHandler
  {
    /**
     * \brief Output format handler for events export
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date Sep 2016
     */
    class ExportHandler
    {
      public:
        /// All types of output available for export
        enum OutputType {
          HepMC, LHE
        };
        friend std::ostream& operator<<( std::ostream& os, const OutputType& type ) {
          switch ( type ) {
            case HepMC: return os << "HepMC ASCII";
            case LHE: return os << "LHEF";
          }
          return os;
        }

      public:
        /// Class constructor
        /// \param[in] type Requested output type
        ExportHandler( const OutputType& type ) :
          type_( type ), cross_sect_( 0. ), cross_sect_err_( 0. ), event_num_( 0. ) {}
        virtual ~ExportHandler() {}
        virtual void initialise( const Parameters& ) = 0;
        /// Set the process cross section and its associated error
        void setCrossSection( const float& xsec, const float& err_xsec ) {
          cross_sect_ = xsec;
          cross_sect_err_ = err_xsec;
        }
        /// Set the event number
        void setEventNumber( const unsigned int& ev_id ) { event_num_ = ev_id; }
        /// Writer operator
        virtual void operator<<( const Event* ) = 0;

      protected:
        /// Type of output requested
        OutputType type_;
        /// Process cross section
        float cross_sect_;
        /// Error on process cross section
        float cross_sect_err_;
        /// Event number in generation
        unsigned int event_num_;
    };
  }
}

#endif
