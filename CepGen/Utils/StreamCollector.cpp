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

#include <fcntl.h>
#include <unistd.h>

#include <cstdio>

#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/StreamCollector.h"

namespace cepgen::utils {
  enum PIPES { READ = 0, WRITE = 1 };

  StreamCollector::StreamCollector(std::string& captured_stream) : captured_stream_(captured_stream) {
    std::lock_guard<std::mutex> lock(mutex_);
    setvbuf(stdout, nullptr, _IONBF, 0);
    setvbuf(stderr, nullptr, _IONBF, 0);
    secure_pipe(pipes_.data());
    old_stdout_ = secure_dup(STD_OUT_FD);
    old_stderr_ = secure_dup(STD_ERR_FD);
    secure_dup2(pipes_[WRITE], STD_OUT_FD);
    secure_dup2(pipes_[WRITE], STD_ERR_FD);
#ifndef _MSC_VER
    secure_close(pipes_[WRITE]);
#endif
  }

  StreamCollector::~StreamCollector() {
    std::lock_guard<std::mutex> lock(mutex_);

    captured_stream_.clear();
    secure_dup2(old_stdout_, STD_OUT_FD);
    secure_dup2(old_stderr_, STD_ERR_FD);

    const int bufSize = 1025;
    char buf[bufSize];
    int bytesRead = 0;
    bool fd_blocked(false);
    do {
      bytesRead = 0;
      fd_blocked = false;
#ifdef _MSC_VER
      if (!eof(pipes_[READ]))
        bytesRead = ::read(pipes_[READ], buf, bufSize - 1);
#else
      bytesRead = ::read(pipes_[READ], buf, bufSize - 1);
#endif
      if (bytesRead > 0) {
        buf[bytesRead] = 0;
        captured_stream_ += buf;
      } else if (bytesRead < 0) {
        fd_blocked = (errno == EAGAIN || errno == EWOULDBLOCK || errno == EINTR);
        if (fd_blocked)
          std::this_thread::sleep_for(std::chrono::milliseconds(10));
      }
    } while (fd_blocked || bytesRead == (bufSize - 1));

    secure_close(old_stdout_);
    secure_close(old_stderr_);
    secure_close(pipes_[READ]);
#ifdef _MSC_VER
    secure_close(pipes_[WRITE]);
#endif
  }

  int StreamCollector::secure_dup(int src) {
    int ret = -1;
    bool fd_blocked = false;
    do {
      ret = ::dup(src);
      fd_blocked = (errno == EINTR || errno == EBUSY);
      if (fd_blocked)
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    } while (ret < 0);
    return ret;
  }
  void StreamCollector::secure_pipe(int* pipes) {
    int ret = -1;
    bool fd_blocked = false;
    do {
#ifdef _MSC_VER
      ret = ::pipe(pipes, 65536, O_BINARY);
#else
      ret = ::pipe(pipes) == -1;
#endif
      if (fd_blocked = errno == EINTR || errno == EBUSY; fd_blocked)
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    } while (ret < 0);
  }
  void StreamCollector::secure_dup2(int src, int dest) {
    int ret = -1;
    bool fd_blocked = false;
    do {
      ret = ::dup2(src, dest);
      if (fd_blocked = errno == EINTR || errno == EBUSY; fd_blocked)
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    } while (ret < 0);
  }

  void StreamCollector::secure_close(int& fd) {
    int ret = -1;
    bool fd_blocked = false;
    do {
      ret = close(fd);
      if (fd_blocked = errno == EINTR; fd_blocked)
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    } while (ret < 0);

    fd = -1;
  }
}  // namespace cepgen::utils
