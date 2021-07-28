#include <cmath>

#include "CepGen/Utils/ProgressBar.h"
#include "CepGen/Utils/String.h"

namespace cepgen {
  namespace utils {
    ProgressBar::ProgressBar(size_t tot, size_t freq)
        : bar_length_(std::stoi(env::get("COLUMNS", "60")) - 10),
          bar_pattern_(bar_length_, '='),
          total_(tot),
          frequency_(freq) {}

    void ProgressBar::update(size_t iter) const {
      const size_t percent = iter * 100. / total_;
      if (percent % frequency_ == 0 || iter == total_) {
        int lpad = int(percent / 100. * bar_length_);
        int rpad = bar_length_ - lpad;
        fprintf(stderr, "\r%3zu%% [%.*s%*s]", percent, lpad, bar_pattern_.c_str(), rpad, "");
        fflush(stderr);
        if (iter == total_)
          fprintf(stderr, "\n");
      }
    }
  }  // namespace utils
}  // namespace cepgen
