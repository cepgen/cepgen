#ifndef _PROCESS_H
#define _PROCESS_H

#include "kinematics.h"
#include "event.h"

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
  ~Process();
  /**
   * @brief Returns the weight for this point in the phase-space
   */
  virtual double ComputeWeight()=0;
  /**
   * Sets the phase space point to compute the weight associated to it.
   * @brief Sets the phase space point to compute
   * @param[in] ndim_ The number of dimensions of the point in the phase space
   * @param[in] x_[] The (@a ndim_)-dimensional point in the phase space on
   * which the kinematics and the cross-section are computed
   */
  virtual void SetPoint(const unsigned int ndim_,double x_[]);
  virtual void DumpPoint();
  /**
   * @brief Sets the list of kinematic cuts to apply on the outgoing particles'
   * final state
   * @param[in] cuts_ The Cuts object containing the kinematic parameters
   */
  virtual inline void SetKinematics(Kinematics cuts_) { _cuts=cuts_;std::cout<<"aaa"<<std::endl; };
  /**
   * Is the system's kinematics well defined and compatible with the process ?
   * This check is mandatory to perform the (@a _ndim)-dimensional point's
   * cross-section computation.
   * @brief Is the system's kinematics well defined?
   * @return A boolean stating if the input kinematics and the final states are
   * well defined
   */
  inline bool IsKinematicsDefined() { return _setkin; }
  /**
   * Specifies the incoming particles' kinematics as well as their properties
   * using two Particle objects.
   * @brief Sets the momentum and PDG id for the incoming particles
   * @param[in] ip1_ Information on the first incoming particle
   * @param[in] ip2_ Information on the second incoming particle
   * @return A boolean stating whether or not the incoming kinematics is
   * properly set for this event
   */
  virtual inline bool SetIncomingParticles(Particle ip1_,Particle ip2_) { 
    ip1_.role=(ip1_.pz>0.)?1:2; _ev->AddParticle(&ip1_);
    ip2_.role=(ip2_.pz>0.)?1:2; _ev->AddParticle(&ip2_);
    _setin=true; return _setin;
  };
  /**
   * @brief Sets the PDG id for the outgoing particles
   * @param[in] part_ Role of the particle in the process
   * @param[in] pdgId_ Particle ID according to the PDG convention
   * @return A boolean stating whether or not the outgoing kinematics is
   * properly set for this event
   */
  virtual inline bool SetOutgoingParticles(int part_,int pdgId_) {
    _ev->AddParticle(new Particle(part_, pdgId_));
    return true;
  };
  /**
   * Fills the private Event object with all the Particle object contained
   * in this event.
   * @param[in] symmetrise_ Do we have to symmetrise the event (randomise the
   * production of the positively- and negatively-charged lepton) ?
   * @brief Fills the Event object with the particles' kinematics
   */
  virtual inline void FillKinematics(bool symmetrise_=false) { if (symmetrise_) std::cout << "symmetrised" << std::endl; };
  /**
   * Returns the complete list of Particle with their role in the process for
   * the point considered in the phase space as an Event object.
   * @brief Returns the event content (list of particles with an assigned role)
   * @return The Event object containing all the generated Particle objects
   */
  inline Event* GetEvent() { return this->_ev; };
  inline unsigned int ndim() const { return this->_ndim; };
  inline double x(const unsigned int idx_) { return (idx_>=this->_ndim)?-1.:this->_x[idx_]; }
  inline std::string GetName() { return this->_name; };
 protected:
  double* _x;
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
  std::string _name;
};

#endif
