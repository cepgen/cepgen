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
    ExportHandler( const OutputType& type ) : fType( type ) {;}
    virtual ~ExportHandler() {;}
    /// Set the process cross section and its associated error
    void SetCrossSection( const float& xsec, const float& err_xsec ) {
      fCrossSect = xsec;
      fCrossSectErr = err_xsec;
    }
    /// Set the event number
    void SetEventNumber( const unsigned int& ev_id ) { fEventNum = ev_id; }
    /// Writer operator
    virtual void operator<<( const Event* ) {;}

   protected:
    /// Type of output requested
    OutputType fType;
    /// Process cross section
    float fCrossSect;
    /// Error on process cross section
    float fCrossSectErr;
    /// Event number in generation
    unsigned int fEventNum;
    
  };
}

#endif
