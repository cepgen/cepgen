#include "CepGen/Modules/CardsHandlerFactory.h"
#include "CepGen/Cards/Handler.h"

#include "CepGen/Core/ParametersList.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/Filesystem.h"

#include "CepGen/Parameters.h"

namespace cepgen {
  namespace card {
    Handler::Handler(const ParametersList& params)
        : NamedModule(params), filename_(params.get<std::string>("filename")), rt_params_(new Parameters) {
      if (!filename_.empty())
        parse(filename_, rt_params_);
    }

    Parameters* Handler::parse(const std::string& filename) {
      try {
        auto parser = CardsHandlerFactory::get().build(utils::fileExtension(filename));
        return parser->parse(filename, new Parameters);
      } catch (const std::invalid_argument& err) {
        throw CG_FATAL("Cards:handler") << "Failed to parse the steering card at \"" << filename << "\"!\n"
                                        << err.what();
      }
    }

    void Handler::write(const Parameters* params, const std::string& filename) {
      try {
        auto writer = CardsHandlerFactory::get().build(utils::fileExtension(filename));
        writer->pack(params);
        return writer->write(filename);
      } catch (const std::invalid_argument& err) {
        throw CG_FATAL("Cards:handler") << "Failed to write the configuration to \"" << filename << "\"!\n"
                                        << err.what();
      }
    }
  }  // namespace card
}  // namespace cepgen
