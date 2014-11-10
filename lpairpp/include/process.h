#ifndef _PROCESS_H
#define _PROCESS_H

#include "kinematics.h"
#include "event.h"
#include "physics.h"

/**
 * Class template to define any process to compute using this MC
 * integrator/events generator
 * @author Laurent Forthomme <laurent.forthomme@uclouvain.be>
 * @date January 2014
 */
class Process
{
 public:
  Process();
  virtual ~Process();
  /**
   * @brief Returns the weight for this point in the phase-space
   */
  inline virtual double ComputeWeight() { std::cout << "***WARNING*** Calling ComputeWeight on a non-process!" << std::endl; return -1; };
  /**
   * Specifies the incoming particles' kinematics as well as their properties
   * using two Particle objects.
   * @brief Sets the momentum and PDG id for the incoming particles
   * @param[in] ip1_ Information on the first incoming particle
   * @param[in] ip2_ Information on the second incoming particle
   * @return A boolean stating whether or not the incoming kinematics is
   * properly set for this event
   */
  inline virtual bool SetIncomingParticles(Particle ip1_,Particle ip2_) { 
    double k = 0., *p1 = ip1_.P4(), *p2 = ip2_.P4();
    ip1_.role=(ip1_.Pz()>0.)?1:2; _ev->AddParticle(ip1_);
    ip2_.role=(ip2_.Pz()>0.)?1:2; _ev->AddParticle(ip2_);
    for (int i=0; i<3; i++) k += p1[i]*p2[i];
    _s = ip1_.M2()+ip2_.M2()+2.*(ip1_.E()*ip2_.E()-k);
    _ecm = sqrt(_s);
    _setin=true; return _setin;
  };
  /**
   * @brief Sets the PDG id for the outgoing particles
   * @param[in] part_ Role of the particle in the process
   * @param[in] pdgId_ Particle ID according to the PDG convention
   * @param[in] mothRole_ Integer role of the outgoing particle's mother
   * @return A boolean stating whether or not the outgoing kinematics is
   * properly set for this event
   */
  inline virtual bool SetOutgoingParticles(int part_, ParticleId pdgId_, int mothRole_=-1) {
    _ev->AddParticle(Particle(part_, pdgId_));
    if (mothRole_!=-1) _ev->GetOneByRole(part_)->SetMother(_ev->GetOneByRole(mothRole_));
    return true;
  };
  /**
   * Fills the private Event object with all the Particle object contained
   * in this event.
   * @param[in] symmetrise_ Do we have to symmetrise the event (randomise the
   * production of the positively- and negatively-charged lepton) ?
   * @brief Fills the Event object with the particles' kinematics
   */
  inline virtual void FillKinematics(bool symmetrise_=false) { if (symmetrise_) std::cout << "symmetrised" << std::endl; };
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
  void DumpPoint();
  /**
   * @brief Sets the list of kinematic cuts to apply on the outgoing particles'
   * final state
   * @param[in] cuts_ The Cuts object containing the kinematic parameters
   */
  inline virtual void SetKinematics(Kinematics cuts_) { _cuts=cuts_; };
  /**
   * Is the system's kinematics well defined and compatible with the process ?
   * This check is mandatory to perform the (@a _ndim)-dimensional point's
   * cross-section computation.
   * @brief Is the system's kinematics well defined?
   * @return A boolean stating if the input kinematics and the final states are
   * well defined
   */
  inline bool IsKinematicsDefined() {
    if (_ev->GetByRole(1).size()!=0 and _ev->GetByRole(1).size()!=0) _setin = true;
    if (_ev->GetByRole(3).size()!=0 and _ev->GetByRole(5).size()!=0 and (_ev->GetByRole(6).size()!=0 or _ev->GetByRole(7).size()!=0)) _setout = true;
    _setkin = _setin and _setout;
    return _setkin;
  }
  /**
   * Returns the complete list of Particle with their role in the process for
   * the point considered in the phase space as an Event object.
   * @brief Returns the event content (list of particles with an assigned role)
   * @return The Event object containing all the generated Particle objects
   */
  inline Event* GetEvent() { return this->_ev; };
  /**
   * @brief Returns the number of dimensions on which the integration is performed
   */
  inline unsigned int ndim() const { return this->_ndim; };
  /**
   * @brief Returns the value of a component of the @a _ndim -dimensional point considered
   */
  inline double x(const unsigned int idx_) { return (idx_>=this->_ndim)?-1.:this->_x[idx_]; }
  /**
   * @brief Returns the human-readable name of the process considered
   */
  inline std::string GetName() { return this->_name; };
 protected:
  /**
   * @brief Array of @a _ndim components representing the point on which the
   * weight in the cross-section is computed
   */
  double* _x;
  /**
   * @brief \f$s\f$, squared centre of mass energy of the incoming particles' system, in \f$\mathrm{GeV}^2\f$
   */
  double _s;
  /**
   * @brief \f$\sqrt s\f$, centre of mass energy of the incoming particles' system, in \f$\mathrm{GeV}\f$
   */
  double _ecm;
  /**
   * @brief Number of dimensions on which the integration has to be performed.
   */
  unsigned int _ndim;
  /**
   * @brief Set of cuts to apply on the final phase space
   */
  Kinematics _cuts;
  /**
   * @brief Event object containing all the information on the in- and outgoing
   * particles
   */
  Event* _ev;
  /**
   * @brief Is the phase space point set ?
   */
  bool _point_set;
  /**
   * @brief Are the event's incoming particles set ?
   */
  bool _setin;
  /**
   * @brief Are the event's outgoing particles set ?
   */
  bool _setout;
  /**
   * @brief Is the full event's kinematic set ?
   */
  bool _setkin;
  /**
   * @brief Name of the process (useful for logging and debugging)
   */
  std::string _name;
};

#endif
