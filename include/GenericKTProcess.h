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
  inline virtual void PrepareKTKinematics() { DebugInsideLoop("Dummy kinematics prepared!"); }
  inline virtual double ComputeJacobian() {
    DebugInsideLoop("Dummy Jacobian returned!");
    return 0.;
  }
  inline virtual double ComputeKTFactorisedMatrixElement() {
    DebugInsideLoop("Dummy matrix element returned!");
    return 0.;
  }
  void ComputeOutgoingPrimaryParticlesMasses();
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
  /// Invariant mass of the first outgoing proton (or remnant), in GeV
  double fMX;
  /// Second outgoing proton
  Particle::Momentum fPY;
  /// Invariant mass of the second outgoing proton (or remnant), in GeV
  double fMY;
  
 private:
  /// First intermediate parton (photon, pomeron, ...)
  Particle::ParticleCode kIntermediatePart1;
  /// Second intermediate parton (photon, pomeron, ...)
  Particle::ParticleCode kIntermediatePart2;
  
};
