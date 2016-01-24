#ifndef GenericProcess_h
#define GenericProcess_h

#include "Event.h"
#include "Physics.h"
#include "Kinematics.h"

/**
 * Class template to define any process to compute using this MC integrator/events generator
 * \author Laurent Forthomme <laurent.forthomme@uclouvain.be>
 * \date Jan 2014
 */
class GenericProcess
{
 public:
  /// Type of outgoing process kinematics to be considered (elastic/dissociative final states)
  enum ProcessMode {
    ElasticElastic = 1,
    ElasticInelastic = 2,
    InelasticElastic = 3,
    InelasticInelastic = 4
  };
  /// Human-readable format of a process mode (elastic/dissociative parts)
  friend std::ostream& operator<<(std::ostream& os, const GenericProcess::ProcessMode& pm);
  
  /// Proton structure function to be used in the outgoing state description
  enum StructureFunctions {
    Electron = 1,
    ElasticProton = 2,
    SuriYennie = 11,
    SuriYennieLowQ2 = 12,
    SzczurekUleshchenko = 15,
    FioreVal = 101,
    FioreSea = 102,
    Fiore = 103
  };
  /// Human-readable format of a structure function object
  friend std::ostream& operator<<(std::ostream& os, const GenericProcess::StructureFunctions& sf);

  /// Generic map of particles with their role in the process
  typedef std::map<Particle::Role,Particle::ParticleCode> ParticlesRoleMap;
  /// Pair of particle with their associated role in the process
  typedef std::pair<Particle::Role,Particle::ParticleCode> ParticleWithRole;
  /// Map of all incoming state particles in the process
  typedef ParticlesRoleMap IncomingState;
  /// Map of all outgoing particles in the process
  typedef ParticlesRoleMap OutgoingState;
 
  /// Default constructor for an undefined process
  /// \param[in] name_ Human-readable format of the process name
  GenericProcess(std::string name_="<invalid process>");
  virtual ~GenericProcess();

  /// Restore the Event object to its initial state
  inline void ClearEvent() { fEvent->Restore(); }
  /// Set the kinematics of the incoming state particles
  void SetIncomingKinematics(Particle::Momentum p1, Particle::Momentum p2);
  /// Compute the incoming state kinematics
  void PrepareKinematics();
  
  // --- virtual (process-defined) methods

  /// Set the incoming and outgoing state to be expected in the process
  inline virtual void AddEventContent() {;}
  /// Prepare the process for its integration over the whole phase space
  inline virtual void BeforeComputeWeight() {;}
  /// Compute the weight for this point in the phase-space
  inline virtual double ComputeWeight() { throw Exception(__PRETTY_FUNCTION__, "Calling ComputeWeight on an invalid process!", Fatal); }
  /// Fill the Event object with the particles' kinematics
  /// \param[in] symmetrise_ Symmetrise the event? (randomise the production of positively-
  /// and negatively-charged outgoing central particles)
  inline virtual void FillKinematics(bool symmetrise_=false) {
    Information("Virtual method called");
    if (symmetrise_) Information("The kinematics is symmetrised");
  }
  /// Return the number of dimensions on which the integration has to be performed
  /// \return Number of dimensions on which to integrate
  inline virtual int GetNdim(ProcessMode) const { return 10; }
  /**
   * Sets the phase space point to compute the weight associated to it.
   * @brief Sets the phase space point to compute
   * @param[in] ndim_ The number of dimensions of the point in the phase space
   * @param[in] x_[] The (@a ndim_)-dimensional point in the phase space on
   * which the kinematics and the cross-section are computed
   */
  void SetPoint(const unsigned int ndim_,double x_[]);
  /// Dump the evaluated point's coordinates in the standard output stream
  void DumpPoint(const ExceptionType& et);
  /// Set the list of kinematic cuts to apply on the outgoing particles' final state
  /// \param[in] cuts_ The Cuts object containing the kinematic parameters
  inline virtual void SetKinematics(Kinematics cuts_) { fCuts=cuts_; }
  /// \return Event object containing all the generated Particle objects
  /// Complete list of Particle with their role in the process for the point considered
  /// in the phase space, returned as an Event object.
  inline Event* GetEvent() { return fEvent; }
  ///Get the number of dimensions on which the integration is performed
  inline unsigned int ndim() const { return fNumDimensions; }
  /// Get the value of a component of the @a fNumDimensions -dimensional point considered
  inline double x(const unsigned int idx_) { return (idx_>=fNumDimensions)?-1.:fX[idx_]; }
  /// Get a human-readable name of the process considered
  inline std::string GetName() { return fName; }
  
 protected:
  /// Set the incoming and outgoing states to be defined in this process (and prepare the Event object accordingly)
  void SetEventContent(IncomingState is, OutgoingState os);
 
  /// Get a list of pointers to the particles with a given role in the process
  /// \param[in] role role in the process for the particle to retrieve
  /// \return A vector of pointers to Particle objects associated to the role
  inline ParticlesRef GetParticles(const Particle::Role& role) { return fEvent->GetByRole(role); }
  /// Get the pointer to one particle in the event (using its role)
  inline Particle* GetParticle(const Particle::Role& role, unsigned int id=0) {
    if (id==0) return fEvent->GetOneByRole(role);
    ParticlesRef pp = fEvent->GetByRole(role);
    if (!pp.size() or id>pp.size()) return 0;
    return pp.at(id);
  }
  /// Get the pointer to one particle in the event (using its identifier)
  inline Particle* GetParticle(unsigned int id) { return fEvent->GetById(id); }
  
  // --- 
  
  /// Array of @a fNumDimensions components representing the point on which the weight in the cross-section is computed
  double* fX;
  /// List of incoming state particles (including intermediate partons)
  IncomingState fIncomingState;
  /// List of outgoing state particles
  OutgoingState fOutgoingState;
  /// \f$s\f$, squared centre of mass energy of the incoming particles' system, in \f$\mathrm{GeV}^2\f$
  double fS;
  /// \f$\sqrt s\f$, centre of mass energy of the incoming particles' system (in GeV)
  double fSqS;
  /// Number of dimensions on which the integration has to be performed.
  unsigned int fNumDimensions;
  /// Set of cuts to apply on the final phase space
  Kinematics fCuts;
  /// Event object containing all the information on the in- and outgoing particles
  Event* fEvent;
  /// Is the phase space point set?
  bool fIsPointSet;
  /// Are the event's incoming particles set?
  bool fIsInStateSet;
  /// Are the event's outgoing particles set?
  bool fIsOutStateSet;
  /// Is the full event's kinematic set?
  bool fIsKinematicSet;
  /// Name of the process (useful for logging and debugging)
  std::string fName;
  
 private:
  /**
   * Is the system's kinematics well defined and compatible with the process ?
   * This check is mandatory to perform the (@a fNumDimensions)-dimensional point's
   * cross-section computation.
   * @brief Is the system's kinematics well defined?
   * @return A boolean stating if the input kinematics and the final states are
   * well defined
   */
  inline bool IsKinematicsDefined() {
    if (GetParticles(Particle::IncomingBeam1).size()!=0 and GetParticles(Particle::IncomingBeam2).size()!=0) fIsInStateSet = true;
    if  ((GetParticles(Particle::OutgoingBeam1).size()!=0   and GetParticles(Particle::OutgoingBeam2).size()!=0)
     and (GetParticles(Particle::CentralParticle1).size()!=0 or GetParticles(Particle::CentralParticle2).size()!=0)) fIsOutStateSet = true;
    fIsKinematicSet = fIsInStateSet and fIsOutStateSet;
    return fIsKinematicSet;
  }

};

#endif
