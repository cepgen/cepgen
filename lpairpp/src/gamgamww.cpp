#include "gamgamww.h"

GamGamWW::GamGamWW() :
  _sqs(-1.),
  setp1(false), setp2(false), setp3(false), setp5(false), setll(false)
{
  _name = "gamma,gamma->W+,W-";
}

GamGamWW::~GamGamWW()
{}

bool
GamGamWW::SetOutgoingParticles(int part_, int pdgId_)
{
  double mass_, outm, dm;

  if (!_point_set) return false;

  mass_ = GetMassFromPDGId(pdgId_);

  return true;
}

bool
GamGamWW::SetIncomingParticles(Particle ip1_, Particle ip2_)
{
  Particle *p1, *p2;
  int role1, role2;

  role1 = (ip1_.pz>0.) ? 1:2;
  role2 = (ip2_.pz>0.) ? 1:2;
  if (role1==role2) return false;
  ip1_.role = role1;
  ip2_.role = role2;

  this->_ev->AddParticle(&ip1_);
  this->_ev->AddParticle(&ip2_);

  p1 = this->_ev->GetOneByRole(1);
  p2 = this->_ev->GetOneByRole(2);

  _etot = p1->E()+p2->E();
  _ptot = std::sqrt(std::pow(p1->px+p2->px, 2)
                   +std::pow(p1->py+p2->py, 2)
                   +std::pow(p1->pz+p2->pz, 2));

  _setin = p1->Valid() and p2->Valid();
  _setkin = _setin && _setout;
  return _setkin;
}

void
GamGamWW::ComputeCMenergy()
{
  double k;

  k = 0.;
  /*for (int i=0; i<3; i++) {
    k += _p3_p1[i]*_p3_p2[i];
  }
  _s = std::pow(_mp1,2)+std::pow(_mp2,2)+2.*(_ep1*_ep2-k);
  _sqs = sqrt(_s);*/

#ifdef DEBUG
  std::cout << "[GamGamWW::ComputeCMenergy] [DEBUG] Centre of mass energy : " << _sqs << " GeV" << std::endl;
#endif

}

double
GamGamWW::ComputeMX(double x_, double outmass_, double *dw_)
{
  double wx2min, wx2max;
  double mx2, dmx2;

  if (_sqs<0.) {
    this->ComputeCMenergy();
  }
  
  /*wx2min = std::pow(GetMassFromPDGId(2212)+GetMassFromPDGId(211), 2);
  wx2max = std::pow(_sqs-_mp2-2.*outmass_, 2);
  Map(x_, wx2min, wx2max, &mx2, &dmx2);*/

#ifdef DEBUG
  std::cout << "[GamGamWW::ComputeMX] [DEBUG]" << std::endl
	    << "\tMX**2 in range [" << wx2min << ", " << wx2max << "]" << std::endl
	    << "\tx = " << x_ << std::endl
	    << "\tMX**2 = " << mx2 << ", dMX**2 = " << dmx2 << std::endl
	    << "\tMX = " << sqrt(mx2) << ", dMX = " << sqrt(dmx2) << std::endl;
#endif

  *dw_ = sqrt(dmx2);
  return sqrt(mx2);
}

double
GamGamWW::ComputeWeight()
{
  return -1.;
}

void
GamGamWW::FillKinematics(bool symmetrise_)
{
  if (symmetrise_) std::cout << "symmetrise" << std::endl;
}

void
GamGamWW::SetKinematics(Kinematics cuts_)
{
  _cuts = cuts_;
}
