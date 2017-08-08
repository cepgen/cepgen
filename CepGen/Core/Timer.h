#ifndef Timer_h
#define Timer_h

#include <time.h>

/**
 * A generic timer to extract the processing time between two steps in this software's flow
 * \author Laurent Forthomme <laurent.forthomme@cern.ch>
 */
class Timer
{
 public:
  inline Timer() { clock_gettime( CLOCK_REALTIME, &beg_ ); }
  /**
   * Get the time elapsed since the last @a reset call (or class construction)
   * @return Elapsed time (since the last reset), in seconds
   */
  inline double elapsed() {
    clock_gettime( CLOCK_REALTIME, &end_ );
    return end_.tv_sec -beg_.tv_sec+( end_.tv_nsec - beg_.tv_nsec )/1.e9;
  }
  /// Reset the clock counter
  inline void reset() {
    clock_gettime( CLOCK_REALTIME, &beg_ );
  }
 private:
  /// Timestamp marking the beginning of the counter
  timespec beg_;
  /// Timestamp marking the end of the counter
  timespec end_;
};

#endif
