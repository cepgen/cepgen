/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2024  Laurent Forthomme
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

#ifndef CepGen_Utils_StreamCollector_h
#define CepGen_Utils_StreamCollector_h

#ifdef _MSC_VER
#include <io.h>
#define popen _popen
#define pclose _pclose
#define stat _stat
#define dup _dup
#define dup2 _dup2
#define fileno _fileno
#define close _close
#define pipe _pipe
#define read _read
#define eof _eof
#else
#include <unistd.h>
#endif
#include <fcntl.h>
#include <stdio.h>

#include <array>
#include <mutex>
#include <string>
#include <thread>

#ifndef STD_OUT_FD
#define STD_OUT_FD (fileno(stdout))
#endif

#ifndef STD_ERR_FD
#define STD_ERR_FD (fileno(stderr))
#endif

namespace cepgen::utils {
  class StreamCollector {
  public:
    explicit StreamCollector(std::string&);
    virtual ~StreamCollector();

  private:
    int secure_dup(int src);
    void secure_pipe(int* pipes);
    static void secure_dup2(int src, int dest);
    static void secure_close(int& fd);

    std::array<int, 2> pipes_{0, 0};
    int old_stdout_{0}, old_stderr_{0};
    std::mutex mutex_;
    std::string& captured_stream_;
  };
}  // namespace cepgen::utils

#endif
