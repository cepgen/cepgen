#ifndef CepGen_IO_PythiaEventInterface_h
#define CepGen_IO_PythiaEventInterface_h

#include <Pythia8/Pythia.h>
#include <unordered_map>

namespace cepgen {
  class Parameters;
  class Event;
  class Particle;
}

namespace Pythia8
{
  /// Interfacing between CepGen and Pythia8 event definitions
  class CepGenEvent : public LHAup
  {
    public:
      /// List of particles to be included to the event content
      enum struct Type
      {
        centralAndPartons, ///< only include initiators and central system
        centralAndBeamRemnants, ///< include undissociated beam remnants and central system
        centralAndFullBeamRemnants ///< include dissociated beam remnants and central system
      };
      explicit CepGenEvent();
      /// Initialise this conversion object with CepGen parameters
      void initialise( const cepgen::Parameters& );
      /// Feed a new CepGen event to this conversion object
      /// \param[in] ev CepGen event to be fed
      /// \param[in] type Type of storage
      void feedEvent( const cepgen::Event& ev, const Type& type );
      /// Set the cross section for a given process
      /// \param[in] id Process identifier
      /// \param[in] xsec Process cross section, in pb
      /// \param[in] xsec_err Uncertainty on process cross section, in pb
      void setCrossSection( int id, double xsec, double xsec_err );
      /// Specify new process attributes
      /// \param[in] id Process identifier
      /// \param[in] xsec Process cross section, in pb
      /// \param[in] q2_scale Hard event scale \f$Q^2\f$, in GeV\f$^2\f$
      /// \param[in] alpha_qed \f$\alpha_{\rm em}\f$ for this process
      /// \param[in] alpha_qcd \f$\alpha_{\rm s}\f$ for this process
      void setProcess( int id, double xsec, double q2_scale, double alpha_qed, double alpha_qcd );

      /// Feed comments to the LHEF block
      void addComments( const std::string& comments ) { osLHEF << comments; }

      /// Retrieve the CepGen particle index given its Pythia8 event id
      /// \param[in] py_id Pythia8 particle id
      /// \return CepGen particle id
      unsigned short cepgenId( unsigned short py_id ) const;
      /// Retrieve the Pythia8 particle index given its CepGen event id
      /// \param[in] cg_id CepGen particle id
      /// \return Pythia8 particle id
      unsigned short pythiaId( unsigned short cg_id ) const;
      /// Add a CepGen particle to the event content
      void addCepGenParticle( const cepgen::Particle& part, int status = INVALID_ID,
                              const std::pair<int,int>& mothers = { 0, 0 },
                              const std::pair<int,int>& colours = { 0, 0 } );
      /// Register a new Pythia8 / CepGen particle mapping
      /// \param[in] py_id Pythia8 particle id
      /// \param[in] cg_id CepGen particle id
      void addCorresp( unsigned short py_id, unsigned short cg_id );
      /// Print all Pythia8 / CepGen Particles correspondences
      void dumpCorresp() const;

      static constexpr unsigned short INVALID_ID = 999; ///< Invalid id association
      static constexpr unsigned short MIN_COLOUR_INDEX = 501; ///< Minimal colour indexing number

      inline bool setInit() override { return true; }
      inline bool setEvent( int ) override { return true; }

    private:
      std::pair<int,int> findMothers( const cepgen::Event& ev, const cepgen::Particle& p ) const;
      static const double mp_, mp2_;
      bool inel1_, inel2_;
      std::unordered_map<unsigned short, unsigned short> py_cg_corresp_;
      const cepgen::Parameters* params_; // borrowed
  };
}
#endif

