#ifndef GenericKTProcess_h
#define GenericKTProcess_h

#include "GenericProcess.h"

/**
 * A generic kT-factorisation process.
 * First 4 dimensions of the phase space are required for the incoming partons'
 * virtualities (radial and azimuthal coordinates)
 * \brief Class template to define any kT-factorisation process
 * \author Laurent Forthomme <laurent.forthomme@cern.ch>
 * \date Apr 2016
 */
class GenericKTProcess : public GenericProcess
{
 public:
  /**
   * \brief Class constructor
   * \param[in] name_ Human-readable kT-factorised process name
   * \param[in] num_user_dimensions_ Number of additional dimensions required for the user process
   * \param[in] ip1_ First incoming parton
   * \param[in] ip2_ Second incoming parton (if undefined, same as the first)
   * \param[in] op1_ First produced final state particle
   * \param[in] op2_ Second produced final state particle (if undefined, same as the first)
   */
  GenericKTProcess(std::string name_="<generic process>",
                   unsigned int num_user_dimensions_=0,
                   Particle::ParticleCode ip1_=Particle::Photon,
                   Particle::ParticleCode op1_=Particle::Muon,
                   Particle::ParticleCode ip2_=Particle::invalidParticle,
                   Particle::ParticleCode op2_=Particle::invalidParticle);
  ~GenericKTProcess();
  
  void AddEventContent();
  int GetNdim(Kinematics::ProcessMode) const;
  double ComputeWeight();
  void FillKinematics(bool);
  
 protected:
  inline void SetKinematics(const Kinematics& kin_) {
    fCuts = kin_;
    fLogQmin = -10.; // FIXME //lqmin = std::log(std::sqrt(fCuts.q2min));
    fLogQmax = std::log(fCuts.qtmax);
  }
  /// Set the kinematics of the central system before any point computation
  inline virtual void PrepareKTKinematics() { DebugInsideLoop("Dummy kinematics prepared!"); }
  /// Jacobian weight of the point in the phase space for integration
  inline virtual double ComputeJacobian() {
    DebugInsideLoop("Dummy Jacobian returned!"); return 0.;
  }
  /// kT-factorised matrix element (event weight)
  /// \return Weight of the point in the phase space to the integral
  inline virtual double ComputeKTFactorisedMatrixElement() {
    DebugInsideLoop("Dummy matrix element returned!"); return 0.;
  }
  /// Compute the invariant masses of the outgoing protons (or remnants)
  void ComputeOutgoingPrimaryParticlesMasses();
  /// Set the kinematics of the incoming and outgoing protons (or remnants)
  void FillPrimaryParticlesKinematics();
  /// Set the kinematics of the outgoing central system
  inline virtual void FillCentralParticlesKinematics() { DebugInsideLoop("Dummy central particles list filled!"); }
 
  /// Get the elastic flux to be expected at a given x_bjorken / kT
  double ElasticFlux(double x_, double kt2_) const;
  /// Get the inelastic flux to be expected at a given x_bjorken / kT
  double InelasticFlux(double x_, double kt2_, double mx_) const;
  
  /// Minimal log-virtuality of the intermediate parton
  double fLogQmin;
  /// Maximal log-virtuality of the intermediate parton
  double fLogQmax;
  /// Virtuality of the first intermediate parton (photon, pomeron, ...)
  double fQT1;
  /// Azimuthal rotation of the first intermediate parton's transverse virtuality
  double fPhiQT1;
  /// Virtuality of the second intermediate parton (photon, pomeron, ...)
  double fQT2;
  /// Azimuthal rotation of the second intermediate parton's transverse virtuality
  double fPhiQT2;
  
  /// First outgoing proton
  Particle::Momentum fPX;
  /// Second outgoing proton
  Particle::Momentum fPY;
  
 private:
  void AddPartonContent();
  const unsigned int kNumRequiredDimensions;
  /// Number of additional dimensions required for the user process
  /// (in addition to the 4 required for the two partons' transverse momenta)
  unsigned int kNumUserDimensions;
  /// First intermediate parton (photon, pomeron, ...)
  Particle::ParticleCode kIntermediatePart1;
  /// Second intermediate parton (photon, pomeron, ...)
  Particle::ParticleCode kIntermediatePart2;
  /// Type of particle produced in the final state
  Particle::ParticleCode kProducedPart1;
  /// Type of particle produced in the final state
  Particle::ParticleCode kProducedPart2;
  
};

#endif
