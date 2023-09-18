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
#include <random>
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
      Process(const Process&);  ///< Copy constructor for a user process
      virtual ~Process() = default;
      Process& operator=(const Process&);  ///< Assignment operator

      static ParametersDescription description();

      virtual std::unique_ptr<Process> clone() const;  ///< Copy all process attributes into a new object
      /// Compute the phase space point weight
      virtual double computeWeight() = 0;
      /// Fill the Event object with the particles' kinematics
      /// \param[in] symmetrise Symmetrise the event? (randomise the production of positively- and negatively-charged outgoing central particles)
      virtual void fillKinematics(bool symmetrise = false) = 0;

      void clear();       ///< Reset process prior to the phase space and variables definition
      void clearEvent();  ///< Restore the event object to its initial state
      void initialise();  ///< Initialise the process once the kinematics has been set

      const Kinematics& kinematics() const { return kin_; }  ///< Constant reference to the process kinematics
      Kinematics& kinematics() { return kin_; }              ///< Reference to the process kinematics

      // debugging utilities
      double weight(const std::vector<double>&);      ///< Compute the weight for a phase-space point
      void dumpPoint(std::ostream* = nullptr) const;  ///< Dump the coordinate of the phase-space point being evaluated
      void dumpVariables(std::ostream* = nullptr) const;  ///< List all variables handled by this generic process

      /// Number of dimensions on which the integration is performed
      inline size_t ndim() const { return mapped_variables_.size(); }

      bool hasEvent() const { return (bool)event_; }  ///< Does the process contain (and hold) an event?
      const Event& event() const;                     ///< Handled particles objects and their relationships
      Event& event();                                 ///< Event object read/write accessor
      Event* eventPtr();                              ///< Event pointer read/write accessor

      const Momentum& pA() const;                     ///< Positive-z incoming beam particle's 4-momentum
      double mA2() const { return mA2_; }             ///< Positive-z incoming beam particle's squared mass
      double mA() const { return std::sqrt(mA2()); }  ///< Positive-z incoming beam particle's mass
      const Momentum& pB() const;                     ///< Negative-z incoming beam particle's 4-momentum
      double mB2() const { return mB2_; }             ///< Negative-z incoming beam particle's squared mass
      double mB() const { return std::sqrt(mB2()); }  ///< Negative-z incoming beam particle's mass
      const Momentum& pX() const;                     ///< Positive-z outgoing beam particle's 4-momentum
      double mX2() const { return mX2_; }             ///< Positive-z outgoing beam particle's squared mass
      double mX() const { return std::sqrt(mX2()); }  ///< Positive-z outgoing beam particle's mass
      const Momentum& pY() const;                     ///< Negative-z outgoing beam particle's 4-momentum
      double mY2() const { return mY2_; }             ///< Negative-z outgoing beam particle's squared mass
      double mY() const { return std::sqrt(mY2()); }  ///< Negative-z outgoing beam particle's mass

      double x1() const { return x1_; }  ///< Positive-z incoming parton's fractional momentum
      double t1() const { return t1_; }  ///< Positive-z incoming parton's squared mass
      const Momentum& q1() const;        ///< Positive-z incoming parton's 4-momentum
      double x2() const { return x2_; }  ///< Negative-z incoming parton's fractional momentum
      double t2() const { return t2_; }  ///< Negative-z incoming parton's squared mass
      const Momentum& q2() const;        ///< Negative-z incoming parton's 4-momentum

      Momentum& q1();  ///< Positive-z incoming parton's 4-momentum
      Momentum& q2();  ///< Negative-z incoming parton's 4-momentum

    protected:
      static constexpr double NUM_LIMITS = 1.e-3;  ///< Numerical limits for sanity comparisons (MeV/mm-level)

      virtual void addEventContent() = 0;  ///< Set the incoming and outgoing state to be expected in the process
      virtual void prepareKinematics() {}  ///< Compute the incoming state kinematics

      typedef std::map<Particle::Role, pdgid_t> IncomingState;   ///< Map of all incoming state particles in the process
      typedef std::map<Particle::Role, pdgids_t> OutgoingState;  ///< Map of all outgoing particles in the process

      Momentum& pA();  ///< Positive-z incoming beam particle's 4-momentum
      Momentum& pB();  ///< Negative-z incoming beam particle's 4-momentum
      Momentum& pX();  ///< Positive-z outgoing beam particle's 4-momentum
      Momentum& pY();  ///< Negative-z outgoing beam particle's 4-momentum

      double& mX2() { return mX2_; }  ///< Positive-z outgoing beam particle's squared mass
      double& mY2() { return mY2_; }  ///< Negative-z outgoing beam particle's squared mass
      double& t1() { return t1_; }    ///< Positive-z incoming parton's squared mass
      double& x1() { return x1_; }    ///< Positive-z incoming parton's fractional momentum
      double& t2() { return t2_; }    ///< Negative-z incoming parton's squared mass
      double& x2() { return x2_; }    ///< Negative-z incoming parton's fractional momentum

      Momentum& pc(size_t);              ///< Central particle's 4-momentum
      const Momentum& pc(size_t) const;  ///< Central particle's 4-momentum

      double s() const { return s_; }        ///< Two-beam squared centre of mass energy
      double sqrtS() const { return sqs_; }  ///< Two-beam centre of mass energy

      //--- Mandelstam variables
      double shat() const;  ///< \f$\hat s=(p_1+p_2)^2=(p_3+...)^2\f$

      double mp_;   ///< Proton mass, in GeV/c\f$^2\f$
      double mp2_;  ///< Squared proton mass, in GeV\f$^2\f$/c\f$^4\f$

      /// Type of mapping to apply on the variable
      enum class Mapping {
        linear = 0,   ///< a linear \f${\rm d}x\f$ mapping
        exponential,  ///< an exponential \f$\frac{\dot{x}}{x} = \dot{\log x}\f$ mapping
        square,       ///< a square \f${\rm d}x^2=2x\cdot\dot{x}\f$ mapping
        power_law     ///< a power-law mapping inherited from LPAIR
      };
      /// Human-friendly printout of the mapping type
      friend std::ostream& operator<<(std::ostream&, const Mapping&);
      /// Register a variable to be handled and populated whenever a new phase space point weight is to be calculated.
      /// \note To be run once per generation (before any point computation)
      /// \param[out] out Reference to the variable to be mapped
      /// \param[in] type Type of mapping to apply
      /// \param[in] lim Integration limits
      /// \param[in] description Human-readable description of the variable
      Process& defineVariable(double& out, const Mapping& type, const Limits& lim, const std::string& description = "");
      /// Generate and initialise all variables handled by this process
      /// \return Phase space point-dependent component of the Jacobian weight of the point in the phase space for integration
      /// \note To be run at each point computation (therefore, to be optimised!)
      double generateVariables() const;

      /// Set the incoming and outgoing states to be defined in this process (and prepare the Event object accordingly)
      void setEventContent(const IncomingState& ini, const OutgoingState& fin);

      double alphaEM(double q) const;  ///< Compute the electromagnetic running coupling algorithm at a given scale
      double alphaS(double q) const;   ///< Compute the strong coupling algorithm at a given scale

      std::default_random_engine rnd_gen_;  ///< Random number generator engine

    private:
      double s_{-1.};    ///< \f$s\f$, squared centre of mass energy of the two-beam system, in \f$\mathrm{GeV}^2\f$
      double sqs_{-1.};  ///< \f$\sqrt s\f$, centre of mass energy of the two-beam system (in GeV)
      double mA2_{-1.};  ///< first incoming beam particle squared mass
      double mB2_{-1.};  ///< second incoming beam particle squared mass
      double mX2_{-1.};  ///< First diffractive state squared mass
      double mY2_{-1.};  ///< Second diffractive state squared mass
      double t1_{-1.};   ///< First parton virtuality
      double t2_{-1.};   ///< Second parton virtuality
      double x1_{0.};    ///< First parton fractional momentum
      double x2_{0.};    ///< Second parton fractional momentum
      std::unique_ptr<Coupling> alphaem_;  ///< Electromagnetic running coupling algorithm
      std::unique_ptr<Coupling> alphas_;   ///< Strong running coupling algorithm
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
      std::vector<double> point_coord_;  ///< Point coordinate for matrix element computation
      /// Phase space point-independent component of the Jacobian weight of the point in the phase space for integration
      double base_jacobian_{1.};
      Kinematics kin_{ParametersList()};  ///< Set of cuts to apply on the final phase space
      std::unique_ptr<Event> event_;      ///< Event object tracking all information on all particles in the system
    };
    /// Helper typedef for a Process unique pointer
    typedef std::unique_ptr<Process> ProcessPtr;
  }  // namespace proc
}  // namespace cepgen

#endif
