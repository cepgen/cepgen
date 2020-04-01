#include "CepGen/Physics/PDG.h"
#include "CepGen/Core/Exception.h"

#include <iomanip>

namespace cepgen
{
  PDG::PDG()
  {
    define( { invalid, "[...]", "", 0, -1, -1., 0, false } );
    define( { diffractiveProton, "diff_proton", "p\u002A", 0, 0., 0., 3, false } );
    define( { pomeron, "pomeron", "\u2119", 0, 0., 0., 0, false } );
    define( { reggeon, "reggeon", "\u211D", 0, 0., 0., 0, false } );
  }

  //--------------------------------------------------------------------

  std::ostream&
  operator<<( std::ostream& os, const ParticleProperties& prop )
  {
    return os << prop.name << "{"
      << "id=" << prop.pdgid
      << ",desc=" << prop.description
      << ",colours=" << prop.colours
      << ",mass=" << prop.mass
      << ",width=" << prop.width
      << ",charge=" << prop.charge
      << ( prop.fermion ? ",fermion" : "" )
      << "}";
  }

  //--------------------------------------------------------------------

  PDG&
  PDG::get()
  {
    static PDG instance;
    return instance;
  }

  bool
  PDG::has( pdgid_t id ) const
  {
    return particles_.count( id ) > 0;
  }

  const ParticleProperties&
  PDG::operator()( pdgid_t id ) const
  {
    if ( particles_.count( id ) > 0 )
      return particles_.at( id );
    throw CG_ERROR( "PDG" )
      << "No particle with PDG id " << id << " in the catalogue.";
  }

  void
  PDG::define( const ParticleProperties& props )
  {
    CG_DEBUG( "PDG:define" ) << "Adding a new particle with "
      << "PDG id=" << std::setw( 8 ) << props.pdgid << ", " << props;
    particles_[props.pdgid] = props;
  }

  const std::vector<pdgid_t>
  PDG::particles() const
  {
    std::vector<pdgid_t> out;
    for ( const auto& pt : particles_ )
      out.emplace_back( pt.first );
    return out;
  }

  const std::string&
  PDG::name( pdgid_t id ) const
  {
    return operator()( id ).description;
  }

  double
  PDG::colours( pdgid_t id ) const
  {
    return operator()( id ).colours;
  }

  double
  PDG::mass( pdgid_t id ) const
  {
    return operator()( id ).mass;
  }

  double
  PDG::width( pdgid_t id ) const
  {
    return operator()( id ).width;
  }

  double
  PDG::charge( pdgid_t id ) const
  {
    return operator()( id ).charge/3.;
  }

  size_t
  PDG::size() const
  {
    return particles_.size();
  }

  void
  PDG::dump() const
  {
    //--- first build a sorted vector out of the (unsorted) map
    std::vector<std::pair<pdgid_t,ParticleProperties> > tmp;
    for ( const auto& prt : particles_ )
      tmp.emplace_back( prt.first, prt.second );
    std::sort( tmp.begin(), tmp.end(),
      []( const std::pair<pdgid_t,ParticleProperties>& a,
          const std::pair<pdgid_t,ParticleProperties>& b ) {
        return a.first < b.first;
      } );
    //--- then the proper dump begins
    std::ostringstream oss;
    for ( const auto& prt : tmp )
      if ( prt.first != PDG::invalid )
        oss << "\n" << prt.second;
    CG_INFO( "PDG" ) << "List of particles registered:" << oss.str();
  }
}
