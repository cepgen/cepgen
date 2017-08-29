#include "CepGen/Generator.h"
#include "CepGen/Processes/PPtoLL.h"

using namespace std;

int main()
{
  CepGen::Generator g;
  CepGen::Parameters* p = g.parameters.get();
  //p->setProcess( new GamGamLL );
  p->setProcess( new CepGen::Process::PPtoLL );
  p->kinematics.mode = CepGen::Kinematics::ElasticElastic;
  //p->kinematics.mode = CepGen::Kinematics::InelasticElastic;
  //p->kinematics.mode = CepGen::Kinematics::ElasticInelastic;
  p->kinematics.central_cuts[CepGen::Cuts::pt_single] = 5.;
  p->kinematics.central_cuts[CepGen::Cuts::eta_single] = { -2.5, 2.5 };
  p->kinematics.remnant_cuts[CepGen::Cuts::mass] = { 1.07, 320. };
  
  p->dump();
  CepGen::Logger::get().level = CepGen::Logger::DebugInsideLoop;

  const unsigned short ndim = g.numDimensions();
  double x[12];
  for ( unsigned int i=0; i<ndim; i++ ) { x[i] = 0.3; }
  
  cout << g.computePoint( x ) << endl;
  
  return 0;
}
