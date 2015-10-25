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
  enum ProcessMode {
    ElasticElastic = 1,
    ElasticInelastic = 2,
    InelasticElastic = 3,
    InelasticInelastic = 4
  };
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
  enum ParticleRole {
    IncomingBeam1 = 1,
    IncomingBeam2 = 2,
    Parton1 = 41,
    Parton2 = 42,
    CentralSystem = 4,
    OutgoingBeam1 = 3,
    OutgoingBeam2 = 5,
    CentralParticle1 = 6,
    CentralParticle2 = 7
  };
  friend std::ostream& operator<<(std::ostream& os, const GenericProcess::ProcessMode& pm);
  friend std::ostream& operator<<(std::ostream& os, const GenericProcess::StructureFunctions& sf);

  typedef std::map<ParticleRole,Particle::ParticleCode> ParticlesRoleMap;
  typedef std::pair<ParticleRole,Particle::ParticleCode> ParticleWithRole;
  typedef ParticlesRoleMap IncomingState;
  typedef ParticlesRoleMap OutgoingState;
  
  GenericProcess(std::string name_="<invalid process>");
  virtual ~GenericProcess();

  inline virtual void AddEventContent() {;}
  void ClearEvent();
  void PrepareKinematics();
  inline virtual void BeforeComputeWeight() {;}
  /**
   * @brief Returns the weight for this point in the phase-space
   */
  inline virtual double ComputeWeight() { throw Exception(__PRETTY_FUNCTION__, "Calling ComputeWeight on an invalid process!", Fatal); }
  /**
   * Fills the private Event object with all the Particle object contained
   * in this event.
   * @param[in] symmetrise_ Do we have to symmetrise the event (randomise the
   * production of the positively- and negatively-charged lepton) ?
   * @brief Fills the Event object with the particles' kinematics
   */
  inline virtual void FillKinematics(bool symmetrise_=false) { if (symmetrise_) std::cout << "symmetrised" << std::endl; }
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
  /**
   * @brief Dumps the evaluated point's coordinates in the standard output stream
   */
  void DumpPoint(const ExceptionType& et);
  /**
   * @brief Sets the list of kinematic cuts to apply on the outgoing particles'
   * final state
   * @param[in] cuts_ The Cuts object containing the kinematic parameters
   */
  inline virtual void SetKinematics(Kinematics cuts_) { fCuts=cuts_; }
  /**
   * Is the system's kinematics well defined and compatible with the process ?
   * This check is mandatory to perform the (@a fNumDimensions)-dimensional point's
   * cross-section computation.
   * @brief Is the system's kinematics well defined?
   * @return A boolean stating if the input kinematics and the final states are
   * well defined
   */
  inline bool IsKinematicsDefined() {
    if   (fEvent->GetByRole(1).size()!=0 and fEvent->GetByRole(2).size()!=0) fIsInStateSet = true;
    if  ((fEvent->GetByRole(3).size()!=0 and fEvent->GetByRole(5).size()!=0)
     and (fEvent->GetByRole(6).size()!=0 or  fEvent->GetByRole(7).size()!=0)) fIsOutStateSet = true;
    fIsKinematicSet = fIsInStateSet and fIsOutStateSet;
    return fIsKinematicSet;
  }
  /**
   * Returns the complete list of Particle with their role in the process for
   * the point considered in the phase space as an Event object.
   * @brief Returns the event content (list of particles with an assigned role)
   * @return The Event object containing all the generated Particle objects
   */
  inline Event* GetEvent() { return fEvent; }
  /**
   * @brief Returns the number of dimensions on which the integration is performed
   */
  inline unsigned int ndim() const { return fNumDimensions; }
  /**
   * @brief Returns the value of a component of the @a fNumDimensions -dimensional point considered
   */
  inline double x(const unsigned int idx_) { return (idx_>=fNumDimensions)?-1.:fX[idx_]; }
  /**
   * @brief Returns the human-readable name of the process considered
   */
  inline std::string GetName() { return fName; }
 protected:
  inline virtual void AddEventKinematics() {;}
  void SetEventContent(IncomingState is, OutgoingState os);
  /**
   * Specifies the incoming particles' kinematics as well as their properties
   * using two Particle objects.
   * @brief Sets the momentum and PDG id for the incoming particles
   * @param[in] ip1_ Information on the first incoming particle
   * @param[in] ip2_ Information on the second incoming particle
   * @return A boolean stating whether or not the incoming kinematics is
   * properly set for this event
   */
  inline void SetIncomingParticles(Particle ip1_,Particle ip2_) { 
    double k = 0., *p1 = ip1_.P4(), *p2 = ip2_.P4();
    ip1_.role=(ip1_.Pz()>0.)?1:2; fEvent->AddParticle(ip1_);
    ip2_.role=(ip2_.Pz()>0.)?1:2; fEvent->AddParticle(ip2_);
    for (int i=0; i<3; i++) k += p1[i]*p2[i];
    fS = ip1_.M2()+ip2_.M2()+2.*(ip1_.E()*ip2_.E()-k);
    fSqS = sqrt(fS);
  }
  /**
   * @brief Sets the PDG id for the outgoing particles
   * @param[in] part_ Role of the particle in the process
   * @param[in] pdgId_ Particle ID according to the PDG convention
   * @param[in] mothRole_ Integer role of the outgoing particle's mother
   * @return A boolean stating whether or not the outgoing kinematics is
   * properly set for this event
   */
  inline virtual void SetOutgoingParticles(int part_, Particle::ParticleCode pdgId_, int mothRole_=-1) {
    fEvent->AddParticle(Particle(part_, pdgId_));
    if (mothRole_!=-1) fEvent->GetOneByRole(part_)->SetMother(fEvent->GetOneByRole(mothRole_));
  };
  /**
   * @brief Array of @a fNumDimensions components representing the point on which the
   * weight in the cross-section is computed
   */
  double* fX;
  IncomingState fIncomingState;
  OutgoingState fOutgoingState;
  /**
   * @brief \f$s\f$, squared centre of mass energy of the incoming particles' system, in \f$\mathrm{GeV}^2\f$
   */
  double fS;
  /**
   * @brief \f$\sqrt s\f$, centre of mass energy of the incoming particles' system, in \f$\mathrm{GeV}\f$
   */
  double fSqS;
  /**
   * @brief Number of dimensions on which the integration has to be performed.
   */
  unsigned int fNumDimensions;
  /**
   * @brief Set of cuts to apply on the final phase space
   */
  Kinematics fCuts;
  /**
   * @brief Event object containing all the information on the in- and outgoing
   * particles
   */
  Event* fEvent;
  /**
   * @brief Is the phase space point set ?
   */
  bool fIsPointSet;
  /**
   * @brief Are the event's incoming particles set ?
   */
  bool fIsInStateSet;
  /**
   * @brief Are the event's outgoing particles set ?
   */
  bool fIsOutStateSet;
  /**
   * @brief Is the full event's kinematic set ?
   */
  bool fIsKinematicSet;
  /**
   * @brief Name of the process (useful for logging and debugging)
   */
  std::string fName;
  
 private:

};

#endif
