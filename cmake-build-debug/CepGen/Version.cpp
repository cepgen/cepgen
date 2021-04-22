#include "CepGen/Version.h"

namespace cepgen {
  const std::string version::tag = "1.0.alpha1";
  const std::string version::extended = "62ec946(devel)";
  const std::string version::banner = "CepGen version "+version::tag+" ("+version::extended+")\n"
    "Copyright (c) 2020 L. Forthomme.\n"
    "License GPLv3+: GNU GPL version 3 or later <https://gnu.org/licenses/gpl.html>.\n"
    "This is free software: you are free to change and redistribute it.\n"
    "There is NO WARRANTY, to the extent permitted by law.";
}