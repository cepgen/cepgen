#ifndef CepGen_Generator_h
#define CepGen_Generator_h

#include "CepGen/Event/Event.h"

#include <iosfwd>
#include <memory>
#include <functional>

////////////////////////////////////////////////////////////////////////////////

/**
 * \mainpage Foreword
 * This Monte Carlo generator was developed as a modern version of the LPAIR code introduced
 * in the early 1990s by J. Vermaseren *et al*\cite Baranov:1991yq\cite Vermaseren:1982cz. This latter allows to
 * compute the cross-section and to generate events for the \f$\gamma\gamma\to\ell^{+}\ell^{-}\f$
 * process in the scope of high energy physics.
 *
 * Soon after the integration of its matrix element, it was extended as a tool to compute and
 * generate events for any generic 2\f$\rightarrow\f$ 3 central exclusive process.
 * To do so, the main operation performed here is the integration of the matrix element (given as a
 * subset of a Process object) over the full available phase space.
 *
 */

////////////////////////////////////////////////////////////////////////////////

/// Common namespace for this Monte Carlo generator
namespace cepgen
{
  class Integrator;
  class GeneratorWorker;
  class Parameters;

  /// Collection of libraries loaded in the runtime environment
  static std::vector<std::string> loaded_libraries;
  /// Collection of search paths to build the runtime environment
  static std::vector<std::string> search_paths;
  /// Import a shared library in the runtime environment
  void loadLibrary( const std::string&, bool match = false );
  /// Launch the initialisation procedure
  /// \param[in] safe_mode Drop libraries initialisation?
  void initialise( bool safe_mode = false );
  /// Dump this program's header into the standard output stream
  void printHeader();
  /// List the modules registered in the runtime database
  void dumpModules();

  ////////////////////////////////////////////////////////////////////////////////

  /**
   * This object represents the core of this Monte Carlo generator, with its
   * capability to generate the events (using the embedded Vegas object) and to
   * study the phase space in term of the variation of resulting cross section
   * while scanning the various parameters (point \f${\bf x}\f$ in the
   * multi-dimensional phase space).
   *
   * The phase space is constrained using the Parameters object given as an
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
      Generator( bool safe_mode = false );
      /// Core of the Monte Carlo integrator and events generator
      /// \param[in] ip List of input parameters defining the phase space on which to perform the integration
      Generator( Parameters *ip );
      ~Generator();

      /// Pointer to the parameters block
      const Parameters* parameters() const { return parameters_.get(); }
      /// Extracted pointer to the parameters block
      Parameters* parametersPtr() { return parameters_.release(); }
      /// Getter to the run parameters block
      Parameters& parametersRef();
      /// Feed the generator with a Parameters object
      void setParameters( Parameters* ip );
      /// Specify an integrator algorithm configuration
      void setIntegrator( std::unique_ptr<Integrator> );
      /// Remove all references to a previous generation/run
      void clearRun();
      /// Integrate the functional over the whole phase space
      void integrate();
      /**
       * Compute the cross section for the run parameters defined by this object.
       * This returns the cross section as well as the absolute error computed along.
       * \brief Compute the cross-section for the given process
       * \param[out] xsec The computed cross-section, in pb
       * \param[out] err The absolute integration error on the computed cross-section, in pb
       */
      void computeXsection( double& xsec, double& err );
      /// Last cross section computed by the generator
      double crossSection() const { return result_; }
      /// Last error on the cross section computed by the generator
      double crossSectionError() const { return result_error_; }

      //void terminate();
      /// \deprecated Replaced by the generic method \a generate.
      [[deprecated("Please use generate instead")]]
      const Event& generateOneEvent( Event::callback callback = nullptr );
      /// Launch the generation of events
      void generate( size_t num_events = 0, Event::callback callback = nullptr );
      /// Compute one single point from the total phase space
      /// \param[in] x the n-dimensional point to compute
      /// \return the function value for the given point
      double computePoint( const std::vector<double>& x );

   private:
      /// Physical Parameters used in the events generation and cross-section computation
      std::unique_ptr<Parameters> parameters_;
      /// Generator worker instance
      std::unique_ptr<GeneratorWorker> generator_;
      /// Integration algorithm
      std::unique_ptr<Integrator> integrator_;
      /// Cross section value computed at the last integration
      double result_;
      /// Error on the cross section as computed in the last integration
      double result_error_;
  };
}

#endif

