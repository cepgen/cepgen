/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2017-2024  Laurent Forthomme
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
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CepGen_Utils_AbortHandler_h
#define CepGen_Utils_AbortHandler_h

#include <atomic>
#include <csignal>

#include "CepGen/Core/Exception.h"

namespace cepgen::utils {
  extern std::atomic<int> gSignal;
  /// Exception raised when the user terminates the process
  struct RunAbortedException : std::runtime_error {
    RunAbortedException() : std::runtime_error("CepGen run aborted") {}
    ~RunAbortedException() noexcept override { CG_INFO("RunAbortedException") << "Run aborted by user interaction."; }

    const char* what() const noexcept override { return "User abort through C-c."; }
  };

  /// Object handling a user-driven process abortion
  class AbortHandler {
  public:
    /// Define a process abortion procedure
    explicit AbortHandler(int flags = SA_SIGINFO) {
      action_.sa_sigaction = handle_ctrl_c;
      sigemptyset(&action_.sa_mask);
      action_.sa_flags = flags;
      init();
    }

  private:
    static void handle_ctrl_c(int signal, siginfo_t* si, void*) {
      gSignal = signal;
      if (abs(si->si_code) != SIGABRT)
        throw RunAbortedException();
    }
    void init() const {
      if (sigaction(SIGINT, &action_, nullptr) != 0 || sigaction(SIGTERM, &action_, nullptr) != 0)
        throw CG_FATAL("AbortHandler") << "Failed to initialise the C-c handler!";
    }
    struct sigaction action_;
  };
}  // namespace cepgen::utils

#endif
