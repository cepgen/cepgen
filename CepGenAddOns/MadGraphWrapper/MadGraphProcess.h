#ifndef CepGenAddOns_MadGraphWrapper_MadGraphProcess_h
#define CepGenAddOns_MadGraphWrapper_MadGraphProcess_h

#include <memory>
#include <string>

#include "CepGen/Physics/Momentum.h"

// forward-declaration of base MadGraph standalone_cpp process
class CPPProcess;

namespace cepgen {
  /// Wrapper around a generic MadGraph CPPProcess definition
  class MadGraphProcess {
  public:
    MadGraphProcess();
    ~MadGraphProcess();

    const std::string& name() const { return name_; }
    const std::string& description() const { return descr_; }

    void initialise(const std::string&);
    double eval();

    MadGraphProcess& setMomentum(size_t, const Momentum&);
    const std::vector<Momentum>& momenta();
    const std::vector<double>& masses() const;

  private:
    std::unique_ptr<CPPProcess> proc_;
    std::vector<double*> mom_;

    const std::string name_;
    const std::string descr_;
    const std::array<int, 2> incoming_pdgids_;
    const std::vector<int> central_pdgids_;

    std::vector<Momentum> momenta_;
  };
}  // namespace cepgen

#endif
