#include "GenericProcess.h"

/**
 * Class template to define any kT-factorisation process
 * \author Laurent Forthomme <laurent.forthomme@cern.ch>
 * \date Apr 2016
 */
class GenericKTProcess : public GenericProcess
{
 public:
  GenericKTProcess(std::string, Particle::ParticleCode, Particle::ParticleCode);
  ~GenericKTProcess();
  
  void AddEventContent();
  int GetNdim(ProcessMode) const;
  double ComputeWeight();
  
 protected:
  inline void SetKinematics(const Kinematics& kin_) {
    fCuts = kin_;
    fLogQmin = -10.; // FIXME //lqmin = std::log(std::sqrt(fCuts.q2min));
    fLogQmax = std::log(fCuts.qtmax);
  std::cout << fLogQmin << " /// " << fLogQmax << " /// " << fCuts.qtmax  << std::endl;
  }
  /// Set the kinematics of the central system before any point computation
  inline virtual void PrepareKTKinematics() { DebugInsideLoop("Dummy kinematics prepared!"); }
  /// Jacobian weight of the point in the phase space for integration
  inline virtual double ComputeJacobian() {
    DebugInsideLoop("Dummy Jacobian returned!"); return 0.;
  }
  /// kT-factorised matrix element (event weight)
  inline virtual double ComputeKTFactorisedMatrixElement() {
    DebugInsideLoop("Dummy matrix element returned!"); return 0.;
  }
  /// Compute the invariant masses of the outgoing protons (or remnants)
  void ComputeOutgoingPrimaryParticlesMasses();
  /// Set the kinematics of the incoming and outgoing protons (or remnants)
  void FillPrimaryParticlesKinematics();
 
  /// Get the elastic flux to be expected at a given x_bjorken / kT
  double ElasticFlux(double x_, double kt2_) const;
  /// Get the inelastic flux to be expected at a given x_bjorken / kT
  double InelasticFlux(double x_, double kt2_, double mx_) const;
  
  /// Minimal log-virtuality of the intermediate parton
  double fLogQmin;
  /// Maximal log-virtuality of the intermediate parton
  double fLogQmax;
  
  /// First outgoing proton
  Particle::Momentum fPX;
  /// Second outgoing proton
  Particle::Momentum fPY;
  
 private:
  /// First intermediate parton (photon, pomeron, ...)
  Particle::ParticleCode kIntermediatePart1;
  /// Second intermediate parton (photon, pomeron, ...)
  Particle::ParticleCode kIntermediatePart2;
  
};
