#ifndef CepGen_Physics_IncomingBeams_h
#define CepGen_Physics_IncomingBeams_h

#include "CepGen/Core/ParametersList.h"

#include "CepGen/Physics/Modes.h"
#include "CepGen/Physics/Momentum.h"

#include <iosfwd>
#include <vector>
#include <memory>

namespace cepgen {
  enum class KTFlux;
  namespace strfun {
    class Parameterisation;
  }
  namespace formfac {
    class Parameterisation;
  }

  /// Incoming beams characteristics
  struct Beam {
    Beam();             ///< Default constructor
    Momentum momentum;  ///< Incoming particle momentum
    pdgid_t pdg;        ///< PDG identifier for the beam
    mode::Beam mode;    ///< Beam treatment mode
    KTFlux kt_flux;     ///< Type of \f$k_{\rm T}\f$-factorised flux to be considered (if any)
  };
  /// Human-readable description of a beam particle/system
  std::ostream& operator<<(std::ostream&, const Beam&);

  /// Beam/primary particle's kinematics
  class IncomingBeams : private std::pair<Beam, Beam> {
  public:
    IncomingBeams() = default;
    explicit IncomingBeams(const ParametersList&);

    /// List containing all parameters handled
    ParametersList parameters() const;

    const Beam& positive() const { return this->first; }
    Beam& positive() { return this->first; }
    const Beam& negative() const { return this->second; }
    Beam& negative() { return this->second; }

    /// Set the beams for the type of kinematics to consider
    void setMode(const mode::Kinematics&);
    /// Type of kinematics to consider for the phase space
    mode::Kinematics mode() const;

    /// Set the beams momenta according to the centre of mass energy
    void setSqrtS(double sqrts);
    /// Incoming beams squared centre of mass energy (in GeV^2)
    double s() const;
    /// Incoming beams centre of mass energy (in GeV)
    double sqrtS() const;

    /// Form factors evaluator
    formfac::Parameterisation* formFactors() const { return form_factors_.get(); }
    /// Set a form factors evaluator object
    void setFormFactors(std::unique_ptr<formfac::Parameterisation>);

    /// Structure functions evaluator
    strfun::Parameterisation* structureFunctions() const { return str_fun_.get(); }
    /// Set a structure functions evaluator object
    void setStructureFunctions(std::unique_ptr<strfun::Parameterisation>);
    /// Set the integer-type of structure functions evaluator to build
    void setStructureFunctions(int, int);

  private:
    /// Type of form factors to consider
    std::shared_ptr<formfac::Parameterisation> form_factors_;
    /// Type of structure functions to consider
    std::shared_ptr<strfun::Parameterisation> str_fun_;
  };
}  // namespace cepgen

#endif
