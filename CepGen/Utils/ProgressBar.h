#ifndef CepGen_Utils_ProgressBar_h
#define CepGen_Utils_ProgressBar_h

#include <string>

namespace cepgen
{
  namespace utils
  {
    class ProgressBar
    {
      public:
        ProgressBar( size_t tot, size_t freq = 10 );
        void update( size_t iter ) const;

      private:
        static constexpr size_t BAR_LENGTH = 50;
        const std::string bar_pattern_;
        size_t total_, frequency_;
    };
  }
}

#endif
