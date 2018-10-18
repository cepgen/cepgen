#include <Pythia8/Pythia.h>

#ifdef PYTHIA8
namespace cepgen {
  class Parameters;
  class Event;
}

namespace Pythia8
{
  class CepGenEvent : public Pythia8::LHAup
  {
    public:
      explicit CepGenEvent();
      void initialise( const cepgen::Parameters& );
      void feedEvent( const cepgen::Event&, bool full );
      inline bool setInit() override { return true; }
      inline bool setEvent( int ) override { return true; }
      void setCrossSection( int id, double xsec, double xsec_err );
      void setProcess( int id, double xsec, double q2_scale, double alpha_qed, double alpha_qcd );

      void addComments( const std::string& );

      unsigned short cepgenId( unsigned short py_id ) const;
      unsigned short pythiaId( unsigned short cg_id ) const;
      void addCorresp( unsigned short py_id, unsigned short cg_id );
      void dumpCorresp() const;

      static constexpr unsigned short INVALID_ID = 999;
    private:
      static const double mp_, mp2_;
      std::vector<std::pair<unsigned short, unsigned short> > py_cg_corresp_;
      const cepgen::Parameters* params_;
  };
}
#endif
