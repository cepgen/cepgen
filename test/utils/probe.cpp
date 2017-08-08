#include "core/MCGen.h"

using namespace std;

int main()
{
  MCGen g;
  Parameters* p = g.parameters;
  //p->process = new GamGamLL;
  p->process = new PPtoLL;
  p->process_mode = Kinematics::ElasticElastic;
  //p->process_mode = Kinematics::InelasticElastic;
  //p->process_mode = Kinematics::ElasticInelastic;
  p->minpt = 5.;
  p->mineta = -2.5; p->maxeta = 2.5;
  p->minmx = 1.07;
  p->maxmx = 320.;
  
  p->Dump();
  Logger::GetInstance()->Level = Logger::DebugInsideLoop;

  const unsigned int ndim = g.GetNdim();
  double x[ndim];
  for (unsigned int i=0; i<ndim; i++) { x[i] = 0.3;}
  
  cout << g.ComputePoint(x) << endl;
  
  return 0;
}
