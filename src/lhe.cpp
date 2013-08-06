#include "lhe.h"

std::string GetLHEvent(GamGam *g_, bool symm_)
{
  std::stringstream s;
  
  /*Particle ga1, ga2;
  ga1 = g_->GetParticle(41);
  ga2 = g_->GetParticle(42);*/
  
  s << "<event>" << std::endl;
  s << g_->GetParticle(1)->GetLHEline() << std::endl;
  s << g_->GetParticle(2)->GetLHEline() << std::endl;
  s << g_->GetParticle(3)->GetLHEline() << std::endl;
  s << g_->GetParticle(5)->GetLHEline() << std::endl;
  s << g_->GetParticle(6)->GetLHEline(symm_) << std::endl;
  s << g_->GetParticle(7)->GetLHEline(symm_) << std::endl;
  s << "</event>" << std::endl;
  //std::cout << s.str() << std::endl;
  return s.str();
}
