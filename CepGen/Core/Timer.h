#ifndef CepGen_Core_Timer_h
#define CepGen_Core_Timer_h

#include <chrono>

namespace cepgen
{
  /**
   * A generic timer to extract the processing time between two steps in this software's flow
   * \author Laurent Forthomme <laurent.forthomme@cern.ch>
   */
  class Timer
  {
    public:
      inline Timer() { reset(); }
      /**
       * Get the time elapsed since the last @a reset call (or class construction)
       * \return Elapsed time (since the last reset), in seconds
       */
      inline double elapsed() const {
        auto end = std::chrono::high_resolution_clock::now();
        return std::chrono::duration<double>( end-beg_ ).count();
      }
      /// Reset the clock counter
      inline void reset() {
        beg_ = std::chrono::high_resolution_clock::now();
      }

    private:
      /// Timestamp marking the beginning of the counter
      std::chrono::time_point<std::chrono::high_resolution_clock> beg_;
  };
}

#endif
