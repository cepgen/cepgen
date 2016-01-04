#include "Kinematics.h"

Kinematics::Kinematics() :
  ptmin(3.), ptmax(-1.),
  emin(0.), emax(-1.),
  etamin(-999.), etamax(999.),
  mxmin(1.07), mxmax(320.),
  q2min(0.), q2max(1.e5),
  wmin(0.), wmax(-1.),
  ptdiffmin(0.), ptdiffmax(300.) //FIXME need to load this from somewhere else
{}

Kinematics::~Kinematics() {}

void
Kinematics::Dump()
{
  std::string s_mode;
  switch (this->mode) {
    case 0: s_mode = "none"; break;      
    case 1: s_mode = "single"; break;      
    case 2: s_mode = "both"; break;
    default: s_mode = "unknown"; break;
  }
  std::cout
    << std::setfill(' ')
    << __PRETTY_FUNCTION__ << " Dump" << std::endl
    << std::setw(25) << "Cuts mode :" << std::setw(2) << this->mode << "->" << std::setw(4) << s_mode << std::endl    
    << "===== Single leptons" << std::endl
    << std::setw(25) << "Minimal pT :" << std::setw(8) << this->ptmin << std::endl
    << std::setw(25) << "Maximal pT :" << std::setw(8) << this->ptmax << std::endl
    << std::setw(25) << "Minimal energy :" << std::setw(8) << this->emin << std::endl
    << std::setw(25) << "Maximal energy :" << std::setw(8) << this->emax << std::endl
    << std::setw(25) << "Minimal pseudorapidity :" << std::setw(8) << this->etamin << std::endl
    << std::setw(25) << "Maximal pseudorapidity :" << std::setw(8) << this->etamax << std::endl
    << "===== Central kinematics" << std::endl
    << std::setw(25) << "Minimal Q**2 :" << std::setw(8) << this->q2min << std::endl
    << std::setw(25) << "Maximal Q**2 :" << std::setw(8) << this->q2max << std::endl
    << std::setw(25) << "Minimal W :" << std::setw(8) << this->wmin << std::endl
    << std::setw(25) << "Maximal W :" << std::setw(8) << this->wmax << std::endl;
}

