#ifndef CepGen_Integration_Integrand_h
#define CepGen_Integration_Integrand_h

#include <vector>
#include <memory>

namespace cepgen
{
  class Parameters;
  class Event;
  namespace utils { class Timer; }
  namespace proc { class Process; }
  /// Wrapper to the function to be integrated
  class Integrand
  {
    public:
      explicit Integrand( const Parameters* );
      ~Integrand();

      /// Compute the integrand for a given phase space point (or "event")
      /// \param[in] x Phase space point coordinates
      /// \note This weight includes the matrix element of the process
      ///  considered, along with all the kinematic factors, and the cut
      ///  restrictions imposed on this phase space.
      ///  \f${\bf x}=\{x_1,\ldots,x_N\}\f$ is therefore an array of random
      ///  numbers defined inside its boundaries (as normalised so that
      ///  \f$\forall i=1,\ldots,N\f$, \f$0<x_i<1\f$).
      double eval( const std::vector<double>& x );
      /// Phase space dimension
      size_t size() const;
      /// Thread-local physics process
      const proc::Process& process() const { return *process_; }

      /// Specify if the generated events are to be stored
      void setStorage( bool store ) { storage_ = store; }
      /// Are the events currently generated in this run to be stored?
      bool storage() const { return storage_; }

    private:
      const Parameters* params_; ///< Generator-owned runtime parameters
      const std::unique_ptr<utils::Timer> tmr_; ///< A precious timekeeper for event timing
      std::unique_ptr<proc::Process> process_; ///< Local instance of the physics process
      Event* event_; ///< Process-owned event
      bool storage_; ///< Is the next event to be generated to be stored?
  };
}

#endif
