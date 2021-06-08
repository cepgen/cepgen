#ifndef CepGen_Event_EventBrowser_h
#define CepGen_Event_EventBrowser_h

#include "CepGen/Event/Particle.h"
#include <regex>

namespace cepgen {
  class Event;
  namespace utils {
    /**
     * \brief A user-friendly browser for the Event content
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date Jul 2019
     */
    class EventBrowser {
    public:
      EventBrowser() = default;
      /// Get/compute a variable value
      double get(const Event& ev, const std::string& var) const;

    private:
      /// Retrieve a named variable from a particle
      double variable(const Event&, const Particle&, const std::string&) const;
      /// Retrieve a named variable from the whole event
      static double variable(const Event&, const std::string&);

      static const std::regex rgx_select_id_, rgx_select_role_;
      static constexpr double INVALID_OUTPUT = -999.;

      //--- auxiliary helper maps
      const std::unordered_map<std::string, Particle::Role> role_str_ = {{"ib1", Particle::Role::IncomingBeam1},
                                                                         {"ib2", Particle::Role::IncomingBeam2},
                                                                         {"ob1", Particle::Role::OutgoingBeam1},
                                                                         {"ob2", Particle::Role::OutgoingBeam2},
                                                                         {"pa1", Particle::Role::Parton1},
                                                                         {"pa2", Particle::Role::Parton2},
                                                                         {"cs", Particle::Role::CentralSystem},
                                                                         {"int", Particle::Role::Intermediate}};
      typedef double (Momentum::*pMethod)() const;
      /// Mapping of string variables to momentum getter methods
      const std::unordered_map<std::string, pMethod> m_mom_str_ = {{"px", &Momentum::px},
                                                                   {"py", &Momentum::py},
                                                                   {"pz", &Momentum::pz},
                                                                   {"pt", &Momentum::pt},
                                                                   {"eta", &Momentum::eta},
                                                                   {"phi", &Momentum::phi},
                                                                   {"m", &Momentum::mass},
                                                                   {"e", &Momentum::energy},
                                                                   {"p", &Momentum::p},
                                                                   {"pt2", &Momentum::pt2},
                                                                   {"th", &Momentum::theta},
                                                                   {"y", &Momentum::rapidity}};
    };
  }  // namespace utils
}  // namespace cepgen

#endif
