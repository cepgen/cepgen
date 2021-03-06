#ifndef CepGen_Utils_AbortHandler_h
#define CepGen_Utils_AbortHandler_h

#include "CepGen/Core/Exception.h"

#include <csignal>
#include <atomic>

namespace cepgen {
  namespace utils {
    extern std::atomic<int> gSignal;
    /// Exception raised when the user terminates the process
    struct RunAbortedException : NullStream {
      using NullStream::NullStream;
      ~RunAbortedException() override {
        CG_LOG("RunAbortedException");
        CG_INFO("RunAbortedException") << "User abort through C-c.";
      }
    };

    /// Object handling an user-driven process abortion
    class AbortHandler {
    public:
      /// Define a process abortion procedure
      AbortHandler(int flags = SA_SIGINFO) {
        action_.sa_sigaction = handle_ctrl_c;
        sigemptyset(&action_.sa_mask);
        action_.sa_flags = flags;
        init();
      }
      /// Switch on/off multithreading capabilities
      void setMT(bool mt_on = true) {
        if (mt_on)
          action_.sa_sigaction = handle_ctrl_c_mt;
        else
          action_.sa_sigaction = handle_ctrl_c;
        init();
      }

    private:
      static void handle_ctrl_c_mt(int signal, siginfo_t*, void*) { gSignal = signal; }
      static void handle_ctrl_c(int signal, siginfo_t*, void*) {
        gSignal = signal;
        throw RunAbortedException();
      }
      void init() {
        if (sigaction(SIGINT, &action_, nullptr) != 0 || sigaction(SIGTERM, &action_, nullptr) != 0)
          throw CG_FATAL("AbortHandler") << "Failed to initialise the C-c handler!";
      }
      struct sigaction action_;
    };
  }  // namespace utils
}  // namespace cepgen

#endif
