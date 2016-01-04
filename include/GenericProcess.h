#ifndef GenericProcess_h
#define GenericProcess_h

#include "Kinematics.h"
#include "Event.h"
#include "Physics.h"

/**
 * Class template to define any process to compute using this MC
 * integrator/events generator
 * @author Laurent Forthomme <laurent.forthomme@uclouvain.be>
 * @date January 2014
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
  friend std::ostream& operator<<(std::ostream& os, const GenericProcess::StructureFunctions& sf);

  typedef std::map<Particle::Role,Particle::ParticleCode> ParticlesRoleMap;
  typedef std::pair<Particle::Role,Particle::ParticleCode> ParticleWithRole;
  typedef ParticlesRoleMap IncomingState;
  typedef ParticlesRoleMap OutgoingState;
  
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
  /**
   * Fills the private Event object with all the Particle object contained
   * in this event.
   * @param[in] symmetrise_ Do we have to symmetrise the event (randomise the
   * production of the positively- and negatively-charged lepton) ?
   * @brief Fills the Event object with the particles' kinematics
   */
  inline virtual void FillKinematics(bool symmetrise_=false) {
    Information("Virtual method called");
    if (symmetrise_) Information("The kinematics is symmetrised");
  }
  /**
   * @brief Returns the number of dimensions on which the integration has to be performed
   * @param[in] process_mode_ Type of subprocess to consider :
   *  - 1: Elastic-elastic
   *  - 2: Elastic-inelastic
   *  - 3: Inelastic-elastic
   *  - 4: Inelastic-inelastic
   * @return Number of dimensions on which to integrate
   */
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
  /**
   * @brief Sets the list of kinematic cuts to apply on the outgoing particles'
   * final state
   * @param[in] cuts_ The Cuts object containing the kinematic parameters
   */
  inline virtual void SetKinematics(Kinematics cuts_) { fCuts=cuts_; }
  /**
   * Returns the complete list of Particle with their role in the process for
   * the point considered in the phase space as an Event object.
   * @brief Get the event content (list of particles with an assigned role)
   * @return The Event object containing all the generated Particle objects
   */
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
