/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2024  Laurent Forthomme
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

#include "CepGen/Event/Event.h"
#include "CepGen/Physics/Coupling.h"
#include "CepGen/Physics/Kinematics.h"
#include "CepGen/Utils/ProcessVariablesAnalyser.h"
#include "CepGen/Utils/RandomGenerator.h"

/// Location for all physics processes to be generated
namespace cepgen::proc {
  /// Class template to define any process to compute using this MC integrator/events generator
  /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
  /// \date Jan 2014
  class Process : public NamedModule<Process> {
  public:
    explicit Process(const ParametersList&);
    Process(const Process&);  ///< Copy constructor for a user process
    virtual ~Process() = default;
    Process& operator=(const Process&);  ///< Assignment operator

    static ParametersDescription description();

    virtual std::unique_ptr<Process> clone() const;  ///< Copy all process attributes into a new object
    virtual double computeWeight() = 0;              ///< Compute the phase space point weight

    void clear();       ///< Reset process prior to the phase space and variables definition
    void clearEvent();  ///< Restore the event object to its initial state
    void initialise();  ///< Initialise the process once the kinematics has been set

    inline const Kinematics& kinematics() const { return kin_; }  ///< Constant reference to the process kinematics
    inline Kinematics& kinematics() { return kin_; }              ///< Reference to the process kinematics
    void setKinematics();

    // debugging utilities
    double weight(const std::vector<double>&);      ///< Compute the weight for a phase-space point
    void dumpPoint(std::ostream* = nullptr) const;  ///< Dump the coordinate of the phase-space point being evaluated
    void dumpVariables(std::ostream* = nullptr) const;  ///< List all variables handled by this generic process

    inline size_t ndim() const { return mapped_variables_.size(); }  ///< Number of dimensions to perform integration

    /// Does the process contain (and hold) an event?
    inline bool hasEvent() const { return static_cast<bool>(event_); }
    const Event& event() const;  ///< Handled particles objects and their relationships
    Event& event();              ///< Event object read/write accessor
    Event* eventPtr();           ///< Event pointer read/write accessor

    const Momentum& pA() const;                            ///< Positive-z incoming beam particle's 4-momentum
    inline double mA2() const { return mA2_; }             ///< Positive-z incoming beam particle's squared mass
    inline double mA() const { return std::sqrt(mA2()); }  ///< Positive-z incoming beam particle's mass
    const Momentum& pB() const;                            ///< Negative-z incoming beam particle's 4-momentum
    inline double mB2() const { return mB2_; }             ///< Negative-z incoming beam particle's squared mass
    inline double mB() const { return std::sqrt(mB2()); }  ///< Negative-z incoming beam particle's mass
    const Momentum& pX() const;                            ///< Positive-z outgoing beam particle's 4-momentum
    inline double mX2() const { return mX2_; }             ///< Positive-z outgoing beam particle's squared mass
    inline double mX() const { return std::sqrt(mX2()); }  ///< Positive-z outgoing beam particle's mass
    const Momentum& pY() const;                            ///< Negative-z outgoing beam particle's 4-momentum
    inline double mY2() const { return mY2_; }             ///< Negative-z outgoing beam particle's squared mass
    inline double mY() const { return std::sqrt(mY2()); }  ///< Negative-z outgoing beam particle's mass

    inline double x1() const { return x1_; }  ///< Positive-z incoming parton's fractional momentum
    inline double t1() const { return t1_; }  ///< Positive-z incoming parton's squared mass
    const Momentum& q1() const;               ///< Positive-z incoming parton's 4-momentum
    inline double x2() const { return x2_; }  ///< Negative-z incoming parton's fractional momentum
    inline double t2() const { return t2_; }  ///< Negative-z incoming parton's squared mass
    const Momentum& q2() const;               ///< Negative-z incoming parton's 4-momentum

    Momentum& q1();  ///< Positive-z incoming parton's 4-momentum
    Momentum& q2();  ///< Negative-z incoming parton's 4-momentum

    inline double wCM() const { return wcm_; }  ///< Two-parton centre of mass energy

    Momentum& pA();  ///< Positive-z incoming beam particle's 4-momentum
    Momentum& pB();  ///< Negative-z incoming beam particle's 4-momentum
    Momentum& pX();  ///< Positive-z outgoing beam particle's 4-momentum
    Momentum& pY();  ///< Negative-z outgoing beam particle's 4-momentum

    inline double& mX2() { return mX2_; }  ///< Positive-z outgoing beam particle's squared mass
    inline double& mY2() { return mY2_; }  ///< Negative-z outgoing beam particle's squared mass
    inline double& t1() { return t1_; }    ///< Positive-z incoming parton's squared mass
    inline double& x1() { return x1_; }    ///< Positive-z incoming parton's fractional momentum
    inline double& t2() { return t2_; }    ///< Negative-z incoming parton's squared mass
    inline double& x2() { return x2_; }    ///< Negative-z incoming parton's fractional momentum

    Momentum& pc(size_t);              ///< Central particle's 4-momentum
    const Momentum& pc(size_t) const;  ///< Central particle's 4-momentum

    inline double s() const { return s_; }                   ///< Two-beam squared centre of mass energy
    inline double sqrtS() const { return sqs_; }             ///< Two-beam centre of mass energy
    inline double inverseSqrtS() const { return inv_sqs_; }  ///< Inverse two-beam centre of mass energy

    utils::RandomGenerator& randomGenerator() const;  ///< Accessor for this process' random number generator

  protected:
    virtual void addEventContent() = 0;         ///< Set the incoming and outgoing state to be expected in the process
    inline virtual void prepareKinematics() {}  ///< Compute the incoming state kinematics
    virtual void fillKinematics() = 0;          ///< Fill the Event object with the particles' kinematics

    //--- Mandelstam variables
    double shat() const;  ///< \f$\hat s=(p_1+p_2)^2=(p_3+...)^2\f$

    double mp_;   ///< Proton mass, in GeV/c\f$^2\f$
    double mp2_;  ///< Squared proton mass, in GeV\f$^2\f$/c\f$^4\f$

  public:
    /// Type of mapping to apply on the variable
    enum class Mapping {
      linear = 0,   ///< a linear \f${\rm d}x\f$ mapping
      exponential,  ///< an exponential \f$\frac{\dot{x}}{x} = \dot{\log x}\f$ mapping
      square,       ///< a square \f${\rm d}x^2=2x\cdot\dot{x}\f$ mapping
      power_law     ///< a power-law mapping inherited from LPAIR
    };
    friend std::ostream& operator<<(std::ostream&, const Mapping&);  ///< Human-friendly printout of the mapping type
    /// Register a variable to be handled and populated whenever a new phase space point weight is to be calculated.
    /// \note To be run once per generation (before any point computation)
    /// \param[out] out Reference to the variable to be mapped
    /// \param[in] type Type of mapping to apply
    /// \param[in] lim Integration limits
    /// \param[in] name Computer-readable variable name
    /// \param[in] description Human-readable description of the variable
    Process& defineVariable(double& out,
                            const Mapping& type,
                            const Limits& lim,
                            const std::string& name,
                            const std::string& description = "");
    /// Retrieve the physical value for one variable
    double variableValue(size_t i, double x) const;

  protected:
    /// Generate and initialise all variables handled by this process
    /// \return Phase space point-dependent component of the Jacobian weight of the point in the phase space for integration
    /// \note To be run at each point computation (therefore, to be optimised!)
    double generateVariables() const;

    /// Set the incoming and outgoing states to be defined in this process (and prepare the Event object accordingly)
    void setEventContent(const std::unordered_map<Particle::Role, spdgids_t>&);

    double alphaEM(double q) const;  ///< Compute the electromagnetic running coupling algorithm at a given scale
    double alphaS(double q) const;   ///< Compute the strong coupling algorithm at a given scale

    std::unique_ptr<utils::RandomGenerator> rnd_gen_;  ///< Process-local random number generator engine

    inline const std::vector<double>& lastCoordinates() const { return point_coord_; }  ///< Last coordinates fed

  private:
    double s_{-1.};        ///< \f$s\f$, squared centre of mass energy of the two-beam system, in GeV\f${}^2\f$
    double sqs_{-1.};      ///< \f$\sqrt s\f$, centre of mass energy of the two-beam system (in GeV)
    double inv_sqs_{-1.};  ///< \f$s^{-1/2}\f$, inverse CM energy of the two-beam system (in GeV\f${}^{-1}\f$)
    double wcm_{-1.};      ///< two-parton centre of mass energy
    double mA2_{-1.};      ///< first incoming beam particle squared mass
    double mB2_{-1.};      ///< second incoming beam particle squared mass
    double mX2_{-1.};      ///< First diffractive state squared mass
    double mY2_{-1.};      ///< Second diffractive state squared mass
    double t1_{-1.};       ///< First parton virtuality
    double t2_{-1.};       ///< Second parton virtuality
    double x1_{0.};        ///< First parton fractional momentum
    double x2_{0.};        ///< Second parton fractional momentum
    std::unique_ptr<Coupling> alphaem_;  ///< Electromagnetic running coupling algorithm
    std::unique_ptr<Coupling> alphas_;   ///< Strong running coupling algorithm
    /// Handler to a variable mapped by this process
    struct MappingVariable {
      std::string name;         ///< Variable name for debugging
      std::string description;  ///< Human-readable description of the variable
      Limits limits;            ///< Kinematic limits to apply on the variable
      double& value;            ///< Reference to the process variable to generate/map
      Mapping type;             ///< Interpolation type
      size_t index;             ///< Corresponding integration variable
    };
    std::vector<MappingVariable> mapped_variables_;  ///< Collection of variables mapped at the weight generation stage
    std::vector<double> point_coord_;                ///< Point coordinate for matrix element computation
    double base_jacobian_{1.};  ///< Phase space point-independent component of the Jacobian weight for integration
    Kinematics kin_{ParametersList()};  ///< Set of cuts to apply on the final phase space
    std::unique_ptr<Event> event_;      ///< Event object tracking all information on all particles in the system
    friend class utils::ProcessVariablesAnalyser;
  };
  typedef std::unique_ptr<Process> ProcessPtr;  ///< Helper typedef for a Process unique pointer
}  // namespace cepgen::proc

#endif
