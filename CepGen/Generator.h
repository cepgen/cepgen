/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2024  Laurent Forthomme
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

#include <memory>

#include "CepGen/Event/Event.h"
#include "CepGen/Utils/Value.h"

////////////////////////////////////////////////////////////////////////////////

/**
 * \mainpage Foreword
 * This Monte Carlo generator was developed as a modern version of the LPAIR code introduced
 * in the early 1990s by J. Vermaseren *et al*\cite Baranov:1991yq\cite Vermaseren:1982cz. This
 * latter allows to compute the cross-section and to generate events for the
 * \f$\gamma\gamma\to\ell^{+}\ell^{-}\f$ process for ee, ep, and pp collisions.
 *
 * Soon after the integration of its matrix element, it was extended as a tool to compute and
 * generate events for any generic 2\f$\rightarrow\f$ 3 central exclusive process.
 * To do so, the main operation performed here is the integration of the matrix element (given
 * as a subset of a Process object) over the full available phase space.
 *
 */

////////////////////////////////////////////////////////////////////////////////

/// Common namespace for this Monte Carlo generator
namespace cepgen {
  class Integrator;
  class GeneratorWorker;
  class RunParameters;
  namespace proc {
    class Process;
  }

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
  void dumpModules();  ///< List the modules registered in RTE database

  ////////////////////////////////////////////////////////////////////////////////

  /**
   * This object represents the core of this Monte Carlo generator, with its
   * capability to generate the events (using the embedded Vegas object) and to
   * study the phase space in term of the variation of resulting cross section
   * while scanning the various parameters (point \f${\bf x}\f$ in the
   * multi-dimensional phase space).
   *
   * The phase space is constrained using the RunParameters object given as an
   * argument to the constructor, and the differential cross-sections for each
   * value of the array \f${\bf x}\f$ are computed in the \a f-function defined
   * outside (but populated inside) this object.
   *
   * This f-function embeds a Process-inherited object which defines all the
   * methods to compute this differential cross-section as well as the in- and outgoing
   * kinematics associated to each particle.
   *
   * \author Laurent Forthomme <laurent.forthomme@cern.ch>
   * \date Feb 2013
   * \brief Core of the Monte-Carlo generator
   *
   */
  class Generator {
  public:
    /// Core of the Monte Carlo integrator and events generator
    /// \param[in] safe_mode Load the generator without external libraries?
    explicit Generator(bool safe_mode = false);
    explicit Generator(RunParameters*);  ///< Build a MC generator object
    ~Generator();

    const RunParameters& runParameters() const;  ///< Pointer to the parameters block
    RunParameters& runParameters();              ///< Run parameters block
    void setRunParameters(RunParameters* ip);    ///< Feed the generator with a RunParameters object

    void setIntegrator(std::unique_ptr<Integrator>);  ///< Specify an integrator algorithm configuration
    void integrate();                                 ///< Integrate the functional over the phase space of interest

    /// Compute the cross section for the run parameters
    /// \return The computed cross-section and uncertainty, in pb
    Value computeXsection();
    /// Compute the cross section for the run parameters
    /// \param[out] xsec The computed cross-section, in pb
    /// \param[out] err The absolute integration error on the computed cross-section, in pb
    [[deprecated("Please use the parameters-less version")]] void computeXsection(double& cross_section, double& err);
    double crossSection() const { return xsect_; }                     ///< Last cross section computed by the generator
    double crossSectionError() const { return xsect_.uncertainty(); }  ///< Last error on the cross section computed

    void generate(size_t num_events, const std::function<void(const Event&, size_t)>&);  ///< Launch event generation
    /// Launch event generation
    void generate(size_t num_events, const std::function<void(const proc::Process&)>& = nullptr);
    const Event& next();  //</ Generate one event

    /// Compute one single point from the total phase space
    /// \param[in] x the n-dimensional point to compute
    /// \return the function value for the given point
    double computePoint(const std::vector<double>& x);

  private:
    void initialise();       ///< Initialise event generation
    void clearRun();         ///< Remove all references to a previous generation/run
    void resetIntegrator();  ///< Reset integrator algorithm from the user-specified configuration

    std::unique_ptr<RunParameters> parameters_;  ///< Run parameters for event generation and cross-section computation
    std::unique_ptr<GeneratorWorker> worker_;    ///< Generator worker instance
    std::unique_ptr<Integrator> integrator_;     ///< Integration algorithm
    bool initialised_{false};                    ///< Has the event generator already been initialised?
    Value xsect_{-1., -1.};                      ///< Cross section value computed at the last integration
  };
}  // namespace cepgen

#endif
