#ifndef CepGen_Cards_Handler_h
#define CepGen_Cards_Handler_h

#include "CepGen/Modules/NamedModule.h"

namespace cepgen {
  class Parameters;
  class ParametersList;
  /// Location for all steering card parsers/writers
  namespace card {
    /// Base steering card module
    class Handler : public NamedModule<std::string> {
    public:
      /// Build a configuration from an external steering card
      explicit Handler(const ParametersList&);
      virtual ~Handler() = default;

      /// Get the list of runtime parameters parsed
      Parameters* parameters() { return params_; }
      /// Specify runtime parameters to the handler
      virtual void pack(const Parameters*){};
      /// Retrieve a configuration from a parsed steering card
      virtual Parameters* parse(const std::string&, Parameters* params) { return params; }
      /// Build a configuration from a steering card
      static Parameters* parse(const std::string&);
      /// Write the current configuration into a steering card
      virtual void write(const std::string&) const {}
      /// Write a steering card from a configuration
      static void write(const Parameters*, const std::string&);

    protected:
      /// Input filename
      const std::string filename_;
      /// List of parameters parsed from a card handler
      Parameters* params_;
    };
  }  // namespace card
}  // namespace cepgen

#endif
