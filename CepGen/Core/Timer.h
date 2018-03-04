#ifndef CepGen_Core_Timer_h
#define CepGen_Core_Timer_h

#include <time.h>

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
  inline double elapsed() {
    timespec end;
    clock_gettime( clock_, &end );
    return end.tv_sec -beg_.tv_sec+( end.tv_nsec-beg_.tv_nsec )/1.e9;
  }
  /// Reset the clock counter
  inline void reset() {
    clock_gettime( clock_, &beg_ );
  }
 private:
  static constexpr clockid_t clock_ = CLOCK_REALTIME;
  /// Timestamp marking the beginning of the counter
  timespec beg_;
};

#endif
