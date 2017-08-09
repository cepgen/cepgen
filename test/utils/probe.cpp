#include "CepGen/Generator.h"

using namespace std;

int main()
{
  CepGen::Generator g;
  CepGen::Parameters* p = g.parameters.get();
  //p->setProcess( new GamGamLL );
  p->setProcess( new CepGen::Process::PPtoLL );
  p->process_mode = CepGen::Kinematics::ElasticElastic;
  //p->process_mode = CepGen::Kinematics::InelasticElastic;
  //p->process_mode = CepGen::Kinematics::ElasticInelastic;
  p->minpt = 5.;
  p->mineta = -2.5; p->maxeta = 2.5;
  p->minmx = 1.07;
  p->maxmx = 320.;
  
  p->dump();
  CepGen::Logger::get().level = CepGen::Logger::DebugInsideLoop;

  const unsigned short ndim = g.numDimensions();
  double x[ndim];
  for (unsigned int i=0; i<ndim; i++) { x[i] = 0.3;}
  
  cout << g.computePoint(x) << endl;
  
  return 0;
}
