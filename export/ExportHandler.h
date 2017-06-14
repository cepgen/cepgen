#ifndef OutputHandler_ExportHandler_h
#define OutputHandler_ExportHandler_h

#include "physics/Event.h"

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

   public:
    /// Class constructor
    /// \param[in] type Requested output type
    ExportHandler( const OutputType& type ) : type_( type ) {}
    virtual ~ExportHandler() {}
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

#endif
