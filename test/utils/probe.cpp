#include "CepGen/Generator.h"
#include "CepGen/Processes/PPtoLL.h"

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
  p->kinematics.pt_min = 5.;
  p->kinematics.eta_min = -2.5; p->kinematics.eta_max = 2.5;
  p->kinematics.mx_min = 1.07;
  p->kinematics.mx_max = 320.;
  
  p->dump();
  CepGen::Logger::get().level = CepGen::Logger::DebugInsideLoop;

  const unsigned short ndim = g.numDimensions();
  double x[12];
  for ( unsigned int i=0; i<ndim; i++ ) { x[i] = 0.3; }
  
  cout << g.computePoint( x ) << endl;
  
  return 0;
}
