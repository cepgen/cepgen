#include "CepGen/Cards/Handler.h"
#include "CepGen/Generator.h"  // for library loading
#include "CepGen/Parameters.h"

#include "CepGen/Event/Event.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Core/EventModifier.h"
#include "CepGen/Core/ExportModule.h"
#include "CepGen/Processes/Process.h"
#include "CepGen/StructureFunctions/Parameterisation.h"

#include "CepGen/Modules/CardsHandlerFactory.h"
#include "CepGen/Modules/ProcessesFactory.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Modules/EventModifierFactory.h"
#include "CepGen/Modules/ExportModuleFactory.h"

#include "CepGen/Utils/TimeKeeper.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <boost/property_tree/xml_parser.hpp>

namespace pt = boost::property_tree;

namespace cepgen {
  namespace card {
    /// Boost tree configuration cards reader/writer
    class BoostTreeHandler : public Handler {
    public:
      /// Boost tree parser from a configuration card
      explicit BoostTreeHandler(const ParametersList&);
      static std::string description() { return "Boost tree parser/writer"; }

      Parameters* parse(const std::string&, Parameters*) override;
      void pack(const Parameters* params) override;

    protected:
      /// Read and cast a file into the property tree
      virtual void read(const std::string&) = 0;

      /// The BOOST property tree translated by this configuration
      pt::ptree tree_;

    private:
      static constexpr const char* DAUGH_KEY = "DAUGHTER";
      static constexpr const char* MIN_KEY = "min";
      static constexpr const char* MAX_KEY = "max";

      static constexpr const char* ADDONS_NAME = "addons";
      static constexpr const char* PROCESS_NAME = "process";
      static constexpr const char* KIN_NAME = "kinematics";
      static constexpr const char* INTEGR_NAME = "integrator";
      static constexpr const char* GENERAL_NAME = "general";
      static constexpr const char* GENERATOR_NAME = "generator";
      static constexpr const char* EVT_MOD_SEQ_NAME = "eventSequence";
      static constexpr const char* OUTPUT_NAME = "output";
      static constexpr const char* TIMER_NAME = "timer";
      static constexpr const char* LOGGER_NAME = "logger";

      template <typename T>
      static pt::ptree pack(const std::vector<T>&);
      static pt::ptree pack(const ParametersList&);
      static pt::ptree pack(const Limits&);
      static ParametersList unpack(const pt::ptree&);

      static void add(ParametersList&, const std::string&, const pt::ptree&);

      ParametersList proc_, gen_, log_;
      ParametersList evt_mod_, evt_out_;
    };

    BoostTreeHandler::BoostTreeHandler(const ParametersList& params) : Handler(params) {}

    Parameters* BoostTreeHandler::parse(const std::string& filename, Parameters* params) {
      params_ = params;
      read(filename);

      if (tree_.count(ADDONS_NAME))
        for (const auto& lib : unpack(tree_.get_child(ADDONS_NAME)).keys())
          loadLibrary(lib);

      try {
        proc_ = unpack(tree_.get_child(PROCESS_NAME));
        params_->setProcess(proc::ProcessesFactory::get().build(proc_));
      } catch (const boost::exception&) {
        throw CG_FATAL("BoostTreeHandler") << "Failed to retrieve a valid \"" << PROCESS_NAME << "\" block"
                                           << " in the steering card!";
      }
      try {
        //----- phase space definition
        if (tree_.count(KIN_NAME))
          params_->kinematics = Kinematics(unpack(tree_.get_child(KIN_NAME)));
        if (tree_.count(INTEGR_NAME))
          *params_->integrator += unpack(tree_.get_child(INTEGR_NAME));
        if (tree_.count(GENERAL_NAME))
          *params_->general += unpack(tree_.get_child(GENERAL_NAME));
        if (tree_.count(GENERATOR_NAME))
          params_->generation() = Parameters::Generation(unpack(tree_.get_child(GENERATOR_NAME)));
        if (tree_.count(EVT_MOD_SEQ_NAME)) {
          evt_mod_ = unpack(tree_.get_child(EVT_MOD_SEQ_NAME));
          for (const auto& name : evt_mod_.keys()) {
            const auto& mod = evt_mod_.get<ParametersList>(name);
            if (!mod.empty())
              params_->addModifier(EventModifierFactory::get().build(name, mod));
          }
        }
        if (tree_.count(OUTPUT_NAME)) {
          evt_out_ = unpack(tree_.get_child(OUTPUT_NAME));
          for (const auto& name : evt_out_.keys()) {
            const auto& mod = evt_out_.get<ParametersList>(name);
            if (!mod.empty())
              params_->addOutputModule(io::ExportModuleFactory::get().build(name, mod));
          }
        }
        if (tree_.count(TIMER_NAME))
          params_->setTimeKeeper(new utils::TimeKeeper);
        if (tree_.count(LOGGER_NAME)) {
          log_ = unpack(tree_.get_child(LOGGER_NAME));
          utils::Logger::get().level =
              (utils::Logger::Level)log_.get<int>("level", (int)utils::Logger::Level::information);
          for (const auto& mod : log_.get<std::vector<std::string> >("enabledModules"))
            utils::Logger::get().addExceptionRule(mod);
        }
      } catch (const boost::exception&) {
      } catch (const Exception&) {
      }

      return params_;
    }

    template <>
    pt::ptree BoostTreeHandler::pack<ParametersList>(const std::vector<ParametersList>& vec) {
      pt::ptree out;
      for (const auto& elem : vec)
        out.push_back(std::make_pair("", pack(elem)));
      return out;
    }

    template <typename T>
    pt::ptree BoostTreeHandler::pack(const std::vector<T>& vec) {
      pt::ptree out;
      for (const auto& elem : vec) {
        pt::ptree elem_tree;
        elem_tree.put("", elem);
        out.push_back(std::make_pair("", elem_tree));
      }
      return out;
    }

    template <>
    pt::ptree BoostTreeHandler::pack<double>(const std::vector<double>& vec) {
      pt::ptree out;
      for (const auto& elem : vec) {
        pt::ptree elem_tree;
        std::ostringstream os;  // ensure floating point storage
        os << std::scientific << elem;
        elem_tree.put("", elem);
        out.push_back(std::make_pair("", elem_tree));
      }
      return out;
    }

    pt::ptree BoostTreeHandler::pack(const ParametersList& params) {
      pt::ptree out;
      for (const auto& key : params.keys()) {
        if (params.has<ParametersList>(key))
          out.add_child(key, pack(params.get<ParametersList>(key)));
        else if (params.has<int>(key))
          out.put(key, params.get<int>(key));
        else if (params.has<double>(key)) {
          std::ostringstream os;  // ensure floating point storage
          os << std::scientific << params.get<double>(key);
          out.put(key, os.str());
        } else if (params.has<std::string>(key))
          out.put(key, params.get<std::string>(key));
        else if (params.has<Limits>(key))
          out.add_child(key, pack(params.get<Limits>(key)));
        else if (params.has<std::vector<ParametersList> >(key))
          out.add_child(key, pack(params.get<std::vector<ParametersList> >(key)));
        else if (params.has<std::vector<int> >(key))
          out.add_child(key, pack(params.get<std::vector<int> >(key)));
        else if (params.has<std::vector<double> >(key))
          out.add_child(key, pack(params.get<std::vector<double> >(key)));
        else if (params.has<std::vector<std::string> >(key))
          out.add_child(key, pack(params.get<std::vector<std::string> >(key)));
        else
          throw CG_FATAL("BoostTreeHandler") << "Failed to recast the key \"" << key << "\" "
                                             << "with value \"" << params.getString(key) << "\"!";
      }
      return out;
    }

    pt::ptree BoostTreeHandler::pack(const Limits& lim) {
      pt::ptree out;
      if (lim.hasMin()) {
        pt::ptree min;
        std::ostringstream os;  // ensure floating point storage
        os << std::scientific << lim.min();
        min.put("", os.str());
        out.push_back(std::make_pair(MIN_KEY, min));
      }
      if (lim.hasMax()) {
        pt::ptree max;
        std::ostringstream os;  // ensure floating point storage
        os << std::scientific << lim.max();
        max.put("", os.str());
        out.push_back(std::make_pair(MAX_KEY, max));
      }
      return out;
    }

    ParametersList BoostTreeHandler::unpack(const pt::ptree& tree) {
      ParametersList out;
      if (tree.empty())
        throw NullStream();
      for (const auto& it : tree) {
        try {
          if (it.first.empty())  // this might be a vector
            try {
              out.operator[]<std::vector<ParametersList> >(DAUGH_KEY).emplace_back(unpack(it.second));
            } catch (const boost::exception&) {
              try {
                out.operator[]<std::vector<double> >(DAUGH_KEY).emplace_back(it.second.get_value<double>());
              } catch (const boost::exception&) {
                try {
                  out.operator[]<std::vector<int> >(DAUGH_KEY).emplace_back(it.second.get_value<int>());
                } catch (const boost::exception&) {
                  out.operator[]<std::vector<std::string> >(DAUGH_KEY).emplace_back(it.second.get_value<std::string>());
                }
              }
            }
          else
            add(out, it.first, it.second);
        } catch (const Exception&) {
          if (it.second.get_value<std::string>().find('.') != std::string::npos)
            try {
              out.set<double>(it.first, it.second.get_value<double>());
            } catch (const boost::exception&) {
              out.set<std::string>(it.first, it.second.get_value<std::string>());
            }
          else
            try {
              out.set<int>(it.first, it.second.get_value<int>());
            } catch (const boost::exception&) {
              out.set<std::string>(it.first, it.second.get_value<std::string>());
            }
        }
      }
      return out;
    }

    void BoostTreeHandler::add(ParametersList& base, const std::string& name, const pt::ptree& tree) {
      auto plist = unpack(tree);
      //--- first check if we have a limits set
      if (plist.keys().size() <= 2 && (plist.has<double>(MIN_KEY) || plist.has<double>(MAX_KEY))) {
        Limits lim;
        plist.fill<double>(MIN_KEY, lim.min());
        plist.fill<double>(MAX_KEY, lim.max());
        base.set<Limits>(name, lim);
      }
      //--- then check if daughter is a vector; if true, skip one hierarchy level
      else if (plist.has<std::vector<int> >(DAUGH_KEY))
        base.set<std::vector<int> >(name, plist.get<std::vector<int> >(DAUGH_KEY));
      else if (plist.has<std::vector<double> >(DAUGH_KEY)) {
        auto vec = plist.get<std::vector<double> >(DAUGH_KEY);
        base.set<std::vector<double> >(name, vec);
      } else if (plist.has<std::vector<std::string> >(DAUGH_KEY))
        base.set<std::vector<std::string> >(name, plist.get<std::vector<std::string> >(DAUGH_KEY));
      else
        base.set<ParametersList>(name, plist);
    }

    void BoostTreeHandler::pack(const Parameters* params) {
      params_ = const_cast<Parameters*>(params);
      tree_.add_child(PROCESS_NAME, pack(params_->process().parameters()));
      if (params_->integrator && !params_->integrator->keys().empty())
        tree_.add_child(INTEGR_NAME, pack(*params_->integrator));
      if (params_->general && !params_->general->keys().empty())
        tree_.add_child(GENERAL_NAME, pack(*params_->general));

      //----- kinematics block
      tree_.add_child(KIN_NAME, pack(params_->kinematics.parameters()));

      //----- generation block
      tree_.add_child(GENERATOR_NAME, pack(params_->generation().parameters()));

      //----- event modification and output
      if (!params_->eventModifiersSequence().empty()) {
        auto evt_mod_tree = pack(evt_mod_);
        for (const auto& mod : params_->eventModifiersSequence())
          evt_mod_tree.put("", mod->name());
        tree_.add_child(EVT_MOD_SEQ_NAME, evt_mod_tree);
      }
      if (!params_->outputModulesSequence().empty()) {
        auto out_mod_tree = pack(evt_out_);
        for (const auto& mod : params_->outputModulesSequence())
          out_mod_tree.add_child(mod->name(), pack(mod->parameters()));
        tree_.add_child(OUTPUT_NAME, out_mod_tree);
      }

      //----- timing and logging
      if (params_->timeKeeper())
        tree_.add_child(TIMER_NAME, pack(ParametersList()));
      log_.set<int>("level", (int)utils::Logger::get().level);
      //FIXME not yet implemented
      //for ( const auto& mod :  utils::Logger::get().exceptionRules() )
      //  log_.operator[]<std::vector<std::string> >( "enabledModules" ).emplace_back( mod );
      tree_.add_child(LOGGER_NAME, pack(log_));
    }

    //------------------------------------------------------------------
    // class specialisations for each Boost format to be handled
    //------------------------------------------------------------------

    /// A JSON configuration file parser
    class JsonHandler : public BoostTreeHandler {
      using BoostTreeHandler::BoostTreeHandler;
      void read(const std::string& filename) override { pt::read_json(filename, tree_); }
      void write(const std::string& filename) const override { pt::write_json(filename, tree_); }
    };

    /// An INFO configuration file parser
    class InfoHandler : public BoostTreeHandler {
      using BoostTreeHandler::BoostTreeHandler;
      void read(const std::string& filename) override { pt::read_info(filename, tree_); }
      void write(const std::string& filename) const override { pt::write_info(filename, tree_); }
    };

    /// An XML configuration file parser
    class XmlHandler : public BoostTreeHandler {
      using BoostTreeHandler::BoostTreeHandler;
      void read(const std::string& filename) override { pt::read_xml(filename, tree_); }
      void write(const std::string& filename) const override {
        std::ofstream file(filename);
        pt::write_xml(file, tree_);
      }
    };
  }  // namespace card
}  // namespace cepgen

REGISTER_CARD_HANDLER(".json", JsonHandler)
REGISTER_CARD_HANDLER(".info", InfoHandler)
REGISTER_CARD_HANDLER(".xml", XmlHandler)
