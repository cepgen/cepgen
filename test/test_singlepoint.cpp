#include "CepGen/Generator.h"
#include "CepGen/Parameters.h"
#include "CepGen/Processes/GamGamLL.h"

using namespace std;

int main()
{
  CepGen::Generator g;
  CepGen::Parameters* p = g.parameters.get();
  //p->setProcess( new GamGamLL );
  p->setProcess( new CepGen::Process::GamGamLL );
  p->kinematics.mode = CepGen::Kinematics::Mode::ElasticElastic;
  //p->kinematics.mode = CepGen::Kinematics::Mode::InelasticElastic;
  //p->kinematics.mode = CepGen::Kinematics::Mode::ElasticInelastic;
  p->kinematics.cuts.central.pt_single = 5.;
  p->kinematics.cuts.central.eta_single = { -2.5, 2.5 };
  p->kinematics.cuts.remnants.mass_single = { 1.07, 320. };

  p->dump();
  CepGen::Logger::get().level = CepGen::Logger::Level::debugInsideLoop;

  const unsigned short ndim = g.numDimensions();
  double x[12];
  for ( unsigned int i = 0; i < ndim; ++i )
    x[i] = 0.3;

  cout << g.computePoint( x ) << endl;

  return 0;
}
