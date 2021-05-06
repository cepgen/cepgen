#ifndef CepGenAddOns_MadGraphWrapper_MadGraphInterface_h
#define CepGenAddOns_MadGraphWrapper_MadGraphInterface_h

#include "CepGen/Core/ParametersList.h"
#include "CepGen/Event/Particle.h"

#include <string>
#include <memory>

// forward-declaration of base MadGraph standalone_cpp process
class CPPProcess;

namespace cepgen {
  class MadGraphInterface {
  public:
    MadGraphInterface(const ParametersList&);

    std::string run() const;

  private:
    static std::string runCommand(const std::string&);
    static std::string generateLibrary(const std::string&, const std::string&, const std::string&);
    static std::string generateProcess(const std::string&);

    void prepareCard() const;
    std::string prepareMadGraphProcess() const;
    void linkCards() const;

    const std::string proc_;
    const std::string model_;
    const std::string card_path_;
    const std::string standalone_cpp_path_;
    const std::string tmp_dir_;
    const std::string log_filename_;

    const std::unordered_map<std::string, pdgid_t> mg5_parts_ = {
        {"d", (pdgid_t)1},     {"d~", (pdgid_t)1},    {"u", (pdgid_t)2},    {"u~", (pdgid_t)2},   {"s", (pdgid_t)3},
        {"s~", (pdgid_t)3},    {"c", (pdgid_t)4},     {"c~", (pdgid_t)4},   {"b", (pdgid_t)5},    {"b~", (pdgid_t)5},
        {"t", (pdgid_t)6},     {"t~", (pdgid_t)6},    {"e+", (pdgid_t)11},  {"e-", (pdgid_t)11},  {"ve", (pdgid_t)12},
        {"ve~", (pdgid_t)12},  {"mu+", (pdgid_t)13},  {"mu-", (pdgid_t)13}, {"vm", (pdgid_t)14},  {"vm~", (pdgid_t)14},
        {"tau+", (pdgid_t)15}, {"tau-", (pdgid_t)15}, {"vt", (pdgid_t)16},  {"vt~", (pdgid_t)16}, {"g", (pdgid_t)21},
        {"a", (pdgid_t)22},    {"z", (pdgid_t)23},    {"w+", (pdgid_t)24},  {"w-", (pdgid_t)24},  {"h", (pdgid_t)25},
    };
  };

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
