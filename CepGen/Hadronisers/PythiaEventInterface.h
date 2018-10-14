#include <Pythia8/Pythia.h>

#ifdef PYTHIA8
namespace cepgen {
  class Parameters;
  class Event;
  enum class KinematicsMode;
}

namespace Pythia8
{
  class CepGenEvent : public Pythia8::LHAup
  {
    public:
      explicit CepGenEvent( const cepgen::Parameters* );
      void feedEvent( const cepgen::Event& ev, bool full, const cepgen::KinematicsMode& );
      bool setInit() override;
      bool setEvent( int ) override;
      void setCrossSection( int id, double xsec, double xsec_err );
      void setProcess( int id, double xsec, double q2_scale, double alpha_qed, double alpha_qcd );

      unsigned short cepgenId( unsigned short py_id ) const;
      unsigned short pythiaId( unsigned short cg_id ) const;
      void addCorresp( unsigned short py_id, unsigned short cg_id );
      void dumpCorresp() const;

      static constexpr unsigned short invalid_id = 999;
    private:
      static const double mp_, mp2_;
      std::vector<std::pair<unsigned short, unsigned short> > py_cg_corresp_;
      const cepgen::Parameters* params_;
  };
}
#endif
