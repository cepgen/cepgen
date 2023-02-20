/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2023  Laurent Forthomme
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CepGen_Process_Process_h
#define CepGen_Process_Process_h

#include <cstddef>  // size_t
#include <map>
#include <memory>
#include <vector>

#include "CepGen/Event/Particle.h"
#include "CepGen/Modules/NamedModule.h"
#include "CepGen/Physics/Coupling.h"
#include "CepGen/Physics/Kinematics.h"

namespace cepgen {
  class Event;
  /// Location for all physics processes to be generated
  namespace proc {
    /// \brief Class template to define any process to compute using this MC integrator/events generator
    /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
    /// \date Jan 2014
    class Process : public NamedModule<std::string> {
    public:
      /// Default constructor for an undefined process
      /// \param[in] params Process-level parameters
      /// \param[in] has_event Do we generate the associated event structure?
      explicit Process(const ParametersList&);
      /// Copy constructor for a user process
      Process(const Process&);
      virtual ~Process() = default;

      static ParametersDescription description();

      /// Reset process prior to the phase space and variables definition
      void clear();

      /// Assignment operator
      Process& operator=(const Process&);

      /// Copy all process attributes into a new object
      virtual std::unique_ptr<Process> clone() const;
      /// Set the incoming and outgoing state to be expected in the process
      inline virtual void addEventContent() = 0;
      /// Compute the phase space point weight
      virtual double computeWeight() = 0;
      /// Compute the incoming state kinematics
      virtual void prepareKinematics() {}
      /// Fill the Event object with the particles' kinematics
      /// \param[in] symmetrise Symmetrise the event? (randomise the production of positively- and negatively-charged outgoing central particles)
      virtual void fillKinematics(bool symmetrise = false) = 0;

      /// Restore the event object to its initial state
      void clearEvent();
      /// Return a constant reference to the process kinematics
      const Kinematics& kinematics() const { return kin_; }
      /// Return a reference to the process kinematics
      Kinematics& kinematics() { return kin_; }
      /// Initialise the process once the kinematics has been set
      void initialise();
      /// Compute the weight for a phase-space point
      double weight(const std::vector<double>&);
      /// Dump the evaluated point's coordinates in the standard output stream
      void dumpPoint() const;
      /// List all variables handled by this generic process
      void dumpVariables(std::ostream* = nullptr) const;

      ///Get the number of dimensions on which the integration is performed
      inline size_t ndim() const { return mapped_variables_.size(); }

      /// Does the process contain (and hold) an event?
      bool hasEvent() const { return (bool)event_; }
      /// Event object containing all the generated Particle objects and their relationships
      const Event& event() const;
      /// Non-const event retrieval method
      Event& event();
      /// Event pointer retrieval method
      Event* eventPtr();

    protected:
      /// Map of all incoming state particles in the process
      typedef std::map<Particle::Role, pdgid_t> IncomingState;
      /// Map of all outgoing particles in the process
      typedef std::map<Particle::Role, std::vector<pdgid_t> > OutgoingState;

      Momentum& pA();                                 ///< Positive-z incoming beam particle's 4-momentum
      const Momentum& pA() const;                     ///< Positive-z incoming beam particle's 4-momentum
      double mA() const { return std::sqrt(mA2()); }  ///< Positive-z incoming beam particle's mass
      double mA2() const { return mA2_; }             ///< Positive-z incoming beam particle's squared mass

      Momentum& pB();                                 ///< Negative-z incoming beam particle's 4-momentum
      const Momentum& pB() const;                     ///< Negative-z incoming beam particle's 4-momentum
      double mB() const { return std::sqrt(mB2()); }  ///< Negative-z incoming beam particle's mass
      double mB2() const { return mB2_; }             ///< Negative-z incoming beam particle's squared mass

      Momentum& pX();                                 ///< Positive-z outgoing beam particle's 4-momentum
      const Momentum& pX() const;                     ///< Positive-z outgoing beam particle's 4-momentum
      double mX() const { return std::sqrt(mX2()); }  ///< Positive-z outgoing beam particle's mass
      double& mX2() { return mX2_; }                  ///< Positive-z outgoing beam particle's squared mass
      double mX2() const { return mX2_; }             ///< Positive-z outgoing beam particle's squared mass

      Momentum& pY();                                 ///< Negative-z outgoing beam particle's 4-momentum
      const Momentum& pY() const;                     ///< Negative-z outgoing beam particle's 4-momentum
      double mY() const { return std::sqrt(mY2()); }  ///< Negative-z outgoinging beam particle's mass
      double& mY2() { return mY2_; }                  ///< Negative-z outgoing beam particle's squared mass
      double mY2() const { return mY2_; }             ///< Negative-z outgoing beam particle's squared mass

      Momentum& q1();                    ///< Positive-z incoming parton's 4-momentum
      const Momentum& q1() const;        ///< Positive-z incoming parton's 4-momentum
      double& t1() { return t1_; }       ///< Positive-z incoming parton's squared mass
      double t1() const { return t1_; }  ///< Positive-z incoming parton's squared mass

      Momentum& q2();                    ///< Negative-z incoming parton's 4-momentum
      const Momentum& q2() const;        ///< Negative-z incoming parton's 4-momentum
      double& t2() { return t2_; }       ///< Negative-z incoming parton's squared mass
      double t2() const { return t2_; }  ///< Negative-z incoming parton's squared mass

      Momentum& pc(size_t);              ///< Central particle's 4-momentum
      const Momentum& pc(size_t) const;  ///< Central particle's 4-momentum

      double s() const { return s_; }        ///< Two-beam squared centre of mass energy
      double sqrtS() const { return sqs_; }  ///< Two-beam centre of mass energy

      //--- Mandelstam variables
      double shat() const;  ///< \f$\hat s=(p_1+p_2)^2=(p_3+...)^2\f$

      double mp_;   ///< Proton mass, in GeV/c\f$^2\f$
      double mp2_;  ///< Squared proton mass, in GeV\f$^2\f$/c\f$^4\f$
      /**
         * \brief Type of mapping to apply on the variable
         */
      enum class Mapping {
        linear = 0,   ///< a linear \f${\rm d}x\f$ mapping
        exponential,  ///< an exponential \f$\frac{\dot{x}}{x} = \dot{\log x}\f$ mapping
        square,       ///< a square \f${\rm d}x^2=2x\cdot\dot{x}\f$ mapping
        power_law     ///< a power-law mapping inherited from LPAIR
        /**
           * Define modified variables of integration to avoid peaks integrations (see \cite Vermaseren:1982cz for details):
           * - \f$y_{\rm out} = x_{\rm min}\left(\frac{x_{\rm max}}{x_{\rm min}}\right)^{\rm exp}\f$ the new variable
           * - \f${\rm d}y_{\rm out} = x_{\rm min}\left(\frac{x_{\rm max}}{x_{\rm min}}\right)^{\rm exp}\log\frac{x_{\rm min}}{x_{\rm max}}\f$, the new variable's differential form
           * \note This method overrides the set of `mapxx` subroutines in ILPAIR, with a slight difference according to the sign of the
           *  \f${\rm d}y_{\rm out}\f$ parameter :
           *  - left unchanged :
           * > `mapw2`, `mapxq`, `mapwx`, `maps2`
           *  - opposite sign :
           * > `mapt1`, `mapt2`
           */
      };
      /// Human-friendly printout of the mapping type
      friend std::ostream& operator<<(std::ostream&, const Mapping&);
      /// Register a variable to be handled and populated whenever
      ///  a new phase space point weight is to be calculated.
      /// \note To be run once per generation (before any point computation)
      /// \param[out] out Reference to the variable to be mapped
      /// \param[in] type Type of mapping to apply
      /// \param[in] in Integration limits
      /// \param[in] default_limits Limits to apply if none retrieved from the user configuration
      /// \param[in] description Human-readable description of the variable
      Process& defineVariable(double& out,
                              const Mapping& type,
                              Limits in = {0., 1.},
                              const Limits& default_limits = {0., 1.},
                              const std::string& description = "");
      /// Generate and initialise all variables handled by this process
      /// \note To be run at each point computation (therefore, to be optimised!)
      void generateVariables() const;
      /// Phase space point-dependent component of the Jacobian weight of the point in the phase space for integration
      double jacobian() const;

      /// Set the incoming and outgoing states to be defined in this process (and prepare the Event object accordingly)
      void setEventContent(const IncomingState& ini, const OutgoingState& fin);

      /// Compute the electromagnetic running coupling algorithm at a given scale
      double alphaEM(double q) const;
      /// Compute the strong coupling algorithm at a given scale
      double alphaS(double q) const;

      /// Numerical limits for sanity comparisons
      static constexpr double NUM_LIMITS = 1.e-3;  // MeV/mm-level

    private:
      /// \f$s\f$, squared centre of mass energy of the incoming particles' system, in \f$\mathrm{GeV}^2\f$
      double s_{-1.};
      /// \f$\sqrt s\f$, centre of mass energy of the incoming particles' system (in GeV)
      double sqs_{-1.};
      double mA2_{-1.};  ///< first incoming beam particle squared mass
      double mB2_{-1.};  ///< second incoming beam particle squared mass
      double mX2_{-1.};  ///< First diffractive state squared mass
      double mY2_{-1.};  ///< Second diffractive state squared mass
      double t1_{-1.};   ///< First parton virtuality
      double t2_{-1.};   ///< Second parton virtuality
      /// Electromagnetic running coupling algorithm
      std::unique_ptr<Coupling> alphaem_;
      /// Strong running coupling algorithm
      std::unique_ptr<Coupling> alphas_;
      /// Handler to a variable mapped by this process
      struct MappingVariable {
        std::string description;  ///< Human-readable description of the variable
        Limits limits;            ///< Kinematic limits to apply on the variable
        double& value;            ///< Reference to the process variable to generate/map
        Mapping type;             ///< Interpolation type
        size_t index;             ///< Corresponding integration variable
      };
      /// Collection of variables to be mapped at the weight generation stage
      std::vector<MappingVariable> mapped_variables_;
      /// Point coordinate for matrix element computation
      std::vector<double> point_coord_;
      /// Phase space point-independent component of the Jacobian weight of the point in the phase space for integration
      double base_jacobian_{1.};
      /// Set of cuts to apply on the final phase space
      Kinematics kin_{ParametersList()};
      /// Event object containing all the information on all particles in the system
      std::unique_ptr<Event> event_;
    };
    /// Helper typedef for a Process unique pointer
    typedef std::unique_ptr<Process> ProcessPtr;
  }  // namespace proc
}  // namespace cepgen

#endif
