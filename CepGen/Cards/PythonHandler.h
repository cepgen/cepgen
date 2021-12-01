/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2021  Laurent Forthomme
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CepGen_Cards_PythonHandler_h
#define CepGen_Cards_PythonHandler_h

#include <Python.h>

#include "CepGen/Cards/Handler.h"

namespace cepgen {
  namespace strfun {
    class Parameterisation;
  }
  class Limits;
  class ParametersList;
  namespace card {
    /// CepGen Python configuration cards reader/writer
    class PythonHandler final : public Handler {
    public:
      /// Read a standard configuration card
      explicit PythonHandler(const ParametersList&);

      static std::string description() { return "Python 2/3 cards parser"; }
      static ParametersDescription parametersDescription();

      Parameters* parse(const std::string&, Parameters*) override;

    private:
      static constexpr const char* ADDONS_NAME = "addons";
      static constexpr const char* TIMER_NAME = "timer";
      static constexpr const char* PROCESS_NAME = "process";
      static constexpr const char* HADR_NAME = "hadroniser";
      static constexpr const char* EVT_MOD_SEQ_NAME = "eventSequence";
      static constexpr const char* LOGGER_NAME = "logger";
      static constexpr const char* INTEGRATOR_NAME = "integrator";
      static constexpr const char* GENERATOR_NAME = "generator";
      static constexpr const char* OUTPUT_NAME = "output";

      static constexpr const char* PDGLIST_NAME = "PDG";
      static constexpr const char* MCD_NAME = "mcdFile";

      void throwPythonError(const std::string&) const;
      std::string pythonPath(const std::string&) const;
      PyObject* element(PyObject*, const std::string&) const;
      PyObject* encode(const char* str) const;
      std::string decode(PyObject* obj) const;

      template <typename T>
      bool is(PyObject* obj) const;
      template <typename T>
      T get(PyObject* obj) const;
      template <typename T>
      bool isVector(PyObject* obj) const;
      template <typename T>
      std::vector<T> getVector(PyObject* obj) const;

      void fillParameter(PyObject* parent, const char* key, bool& out);
      void fillParameter(PyObject* parent, const char* key, int& out);
      void fillParameter(PyObject* parent, const char* key, unsigned long& out);
      void fillParameter(PyObject* parent, const char* key, unsigned int& out);
      void fillParameter(PyObject* parent, const char* key, double& out);
      void fillParameter(PyObject* parent, const char* key, std::string& out);
      void fillParameter(PyObject* parent, const char* key, Limits& out);
      void fillParameter(PyObject* parent, const char* key, std::vector<int>& out);
      void fillParameter(PyObject* parent, const char* key, std::vector<double>& out);
      void fillParameter(PyObject* parent, const char* key, std::vector<std::string>& out);
      void fillParameter(PyObject* parent, const char* key, ParametersList& out);
      void fillParameter(PyObject* parent, const char* key, std::vector<ParametersList>& out);

      void parseLogging(PyObject*);
      void parseIntegrator(PyObject*);
      void parseGenerator(PyObject*);
      void parseHadroniser(PyObject*);
      void parseEventModifiers(PyObject*);
      void parseOutputModule(PyObject*);
      void parseOutputModules(PyObject*);
      void parseExtraParticles(PyObject*);
    };
    template <>
    bool PythonHandler::is<bool>(PyObject* obj) const;
    template <>
    bool PythonHandler::is<int>(PyObject* obj) const;
    template <>
    bool PythonHandler::is<long>(PyObject* obj) const;
    template <>
    int PythonHandler::get<int>(PyObject* obj) const;
    template <>
    unsigned long PythonHandler::get<unsigned long>(PyObject* obj) const;
    template <>
    bool PythonHandler::is<ParametersList>(PyObject* obj) const;
    template <>
    ParametersList PythonHandler::get<ParametersList>(PyObject* obj) const;
    template <>
    bool PythonHandler::is<double>(PyObject* obj) const;
    template <>
    double PythonHandler::get<double>(PyObject* obj) const;
    template <>
    bool PythonHandler::is<std::string>(PyObject* obj) const;
    template <>
    std::string PythonHandler::get<std::string>(PyObject* obj) const;
    template <>
    bool PythonHandler::is<Limits>(PyObject* obj) const;
    template <>
    Limits PythonHandler::get<Limits>(PyObject* obj) const;
  }  // namespace card
}  // namespace cepgen

#endif
