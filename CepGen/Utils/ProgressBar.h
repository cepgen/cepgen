#ifndef CepGen_Utils_ProgressBar_h
#define CepGen_Utils_ProgressBar_h

#include <string>

namespace cepgen
{
  namespace utils
  {
    /// A simple progress indicator
    class ProgressBar
    {
      public:
        /// Progress bar constructor
        /// \param[in] tot Total number of steps before completion
        /// \param[in] freq Frequency at which the tick is refreshed
        ProgressBar( size_t tot, size_t freq = 10 );
        /// Broadcast the current progress to the bar
        /// \param[in] iter Current iteration
        void update( size_t iter ) const;

      private:
        const size_t bar_length_;
        const std::string bar_pattern_;
        size_t total_, frequency_;
    };
  }
}

#endif
