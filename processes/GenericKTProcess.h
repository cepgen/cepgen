#ifndef GenericKTProcess_h
#define GenericKTProcess_h

#include "processes/GenericProcess.h"
#include "physics/FormFactors.h"
#include "physics/PhotonFluxes.h"

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
  GenericKTProcess( const std::string& name_="<generic process>",
                    const unsigned int& num_user_dimensions_=0,
                    const Particle::ParticleCode& ip1_=Particle::Photon,
                    const Particle::ParticleCode& op1_=Particle::Muon,
                    const Particle::ParticleCode& ip2_=Particle::invalidParticle,
                    const Particle::ParticleCode& op2_=Particle::invalidParticle);
  ~GenericKTProcess();

  /// Populate the event content with the generated process' topology
  void AddEventContent();
  /// Retrieve the total number of dimensions on which the integration is being performet
  /// \param[in] proc_mode_ Kinematics case considered
  unsigned int GetNdim( const Kinematics::ProcessMode& proc_mode_ ) const;
  /// Retrieve the event weight in the phase space
  double ComputeWeight();
  /// Populate the event content with the generated process' kinematics  
  void FillKinematics( bool );
  
 protected:
  inline void SetKinematics( const Kinematics& kin_ ) {
    fCuts = kin_;
    fLogQmin = -10.; // FIXME //lqmin = std::log(std::sqrt(fCuts.q2min));
    fLogQmax = std::log(fCuts.qtmax);
  }
  /// Set the kinematics of the central system before any point computation
  inline virtual void PrepareKTKinematics() { DebuggingInsideLoop("Dummy kinematics prepared!"); }
  /// Minimal Jacobian weight of the point considering a kT factorisation
  double MinimalJacobian() const;
  /// Jacobian weight of the point in the phase space for integration
  inline virtual double ComputeJacobian() {
    DebuggingInsideLoop("Dummy Jacobian returned!"); return 0.;
  }
  /// kT-factorised matrix element (event weight)
  /// \return Weight of the point in the phase space to the integral
  inline virtual double ComputeKTFactorisedMatrixElement() {
    DebuggingInsideLoop("Dummy matrix element returned!"); return 0.;
  }
  /// Compute the invariant masses of the outgoing protons (or remnants)
  void ComputeOutgoingPrimaryParticlesMasses();
  void ComputeIncomingFluxes( double, double, double, double );
  /// Set the kinematics of the incoming and outgoing protons (or remnants)
  void FillPrimaryParticlesKinematics();
  /// Set the kinematics of the outgoing central system
  inline virtual void FillCentralParticlesKinematics() { DebuggingInsideLoop("Dummy central particles list filled!"); }
  
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
  /// First incoming parton's flux
  double fFlux1;
  /// Second incoming parton's flux
  double fFlux2;
  
 private:
  void AddPartonContent();
  static const unsigned int kNumRequiredDimensions = 4;
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
