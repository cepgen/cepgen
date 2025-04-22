/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2025  Laurent Forthomme
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
 */

#ifndef CepGen_Generator_h
#define CepGen_Generator_h

#include <functional>
#include <memory>

#include "CepGen/Utils/Value.h"

/// Common namespace for this Monte Carlo generator
namespace cepgen {
  class Event;
  class Integrator;
  class GeneratorWorker;
  class RunParameters;
}  // namespace cepgen
namespace cepgen::proc {
  class Process;
}

namespace cepgen {
  static std::vector<std::string> loaded_libraries;   ///< Collection of libraries loaded in RTE
  static std::vector<std::string> invalid_libraries;  ///< Collection of libraries tested not to work with RTE
  static std::vector<std::string> search_paths;       ///< Collection of search paths to build RTE
  /// Execute an action on a path if found in search paths collection
  bool callPath(const std::string&, bool (*callback)(const std::string&));
  bool loadLibrary(const std::string&, bool match = false);  ///< Import a shared library in RTE
  /// Launch the initialisation procedure
  /// \param[in] safe_mode Drop libraries initialisation?
  void initialise(bool safe_mode = false);
  void printHeader();  ///< Dump this program's header into the standard output stream

  /// Core generator object allowing for process definition, cross-section computation, and event generation
  /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
  /// \date Feb 2013
  class Generator {
  public:
    /// Initialise the Monte Carlo integrator and event generator
    /// \param[in] safe_mode Load the generator without external libraries?
    explicit Generator(bool safe_mode = false);
    explicit Generator(RunParameters*);  ///< Build an MC generator object
    ~Generator();

    void parseRunParameters(const std::string&);  ///< Read a steering card to populate the run parameters block
    const RunParameters& runParameters() const;   ///< Pointer to the parameters block
    RunParameters& runParameters();               ///< Run parameters block
    void setRunParameters(std::unique_ptr<RunParameters>&);  ///< Feed the generator with a RunParameters object

    void setIntegrator(std::unique_ptr<Integrator>);  ///< Specify an integrator algorithm configuration
    Integrator& integrator() const;                   ///< Retrieve the integrator object
    void integrate();                                 ///< Integrate the functional over the phase space of interest

    Value computeXsection();  ///< Compute the cross-section and uncertainty, in pb, for the run parameters
    /// Compute the cross-section for the run parameters
    /// \param[out] cross_section The computed cross-section, in pb
    /// \param[out] err The absolute integration error on the computed cross-section, in pb
    [[deprecated("Please use the parameters-less version")]] void computeXsection(double& cross_section, double& err);
    double crossSection() const { return cross_section_; }  ///< Last cross-section computed by the generator
    double crossSectionError() const {
      return cross_section_.uncertainty();
    }  ///< Last error on the cross-section computed

    void generate(size_t num_events, const std::function<void(const Event&, size_t)>&);            ///< Generate events
    void generate(size_t num_events, const std::function<void(const proc::Process&)>& = nullptr);  ///< Generate events
    const Event& next();  ///< Generate one single event

    /// Compute one single point from the total phase space
    /// \param[in] coordinates the n-dimensional point to compute
    /// \return the function value for the given point
    double computePoint(const std::vector<double>& coordinates);

  private:
    void initialise();       ///< Initialise event generation
    void clearRun();         ///< Remove all references to a previous generation/run
    void resetIntegrator();  ///< Reset integrator algorithm from the user-specified configuration

    std::unique_ptr<RunParameters> parameters_;  ///< Run parameters for event generation and cross-section computation
    std::unique_ptr<GeneratorWorker> worker_;    ///< Generator worker instance
    std::unique_ptr<Integrator> integrator_;     ///< Integration algorithm
    bool initialised_{false};                    ///< Has the event generator already been initialised?
    Value cross_section_{-1., -1.};              ///< Cross-section value computed at the last integration
  };
}  // namespace cepgen

#endif
