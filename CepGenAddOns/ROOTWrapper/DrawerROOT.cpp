#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGenAddOns/ROOTWrapper/ROOTCanvas.h"

namespace cepgen {
  namespace utils {
    class DrawerROOT : public Drawer {
    public:
      explicit DrawerROOT(const ParametersList&);

      const DrawerROOT& draw(const Drawable&) const override;

    private:
    };

    DrawerROOT::DrawerROOT(const ParametersList& params) : Drawer(params) {}

    const DrawerROOT& DrawerROOT::draw(const Drawable& obj) const { return *this; }
  }  // namespace utils
}  // namespace cepgen

REGISTER_DRAWER("root", DrawerROOT)
