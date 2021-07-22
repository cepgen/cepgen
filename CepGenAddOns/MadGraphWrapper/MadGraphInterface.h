#ifndef CepGenAddOns_MadGraphWrapper_MadGraphInterface_h
#define CepGenAddOns_MadGraphWrapper_MadGraphInterface_h

#include <memory>
#include <string>

#include "CepGen/Core/ParametersList.h"
#include "CepGen/Physics/ParticleProperties.h"

// forward-declaration of base MadGraph standalone_cpp process
class CPPProcess;

namespace cepgen {
  class MadGraphInterface {
  public:
    MadGraphInterface(const ParametersList&);

    std::string run() const;

  private:
    static constexpr size_t cmd_buffer_size_ = 256;
    static const std::unordered_map<std::string, pdgid_t> mg5_parts_;

    static std::string runCommand(const std::string&);
    static std::string generateLibrary(const std::string&, const std::string&, const std::string&);
    static std::string generateProcess(const std::string&);

    using ProcessParticles = std::pair<std::vector<pdgid_t>, std::vector<pdgid_t> >;
    static ProcessParticles unpackProcessParticles(const std::string&);

    void prepareCard() const;
    void linkCards() const;
    std::string prepareMadGraphProcess() const;

    const std::string proc_;
    const std::string model_;
    const std::string card_path_;
    const std::string standalone_cpp_path_;
    const std::string tmp_dir_;
    const std::string log_filename_;
  };
}  // namespace cepgen

#endif
