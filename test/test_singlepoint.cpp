#include "CepGen/Generator.h"
#include "CepGen/Parameters.h"
#include "CepGen/Processes/GamGamLL.h"
#include "CepGen/Core/Exception.h"

using namespace std;

int main()
{
  cepgen::Generator g;
  cepgen::Parameters* p = g.parameters.get();
  //p->setProcess( new GamGamLL );
  p->setProcess( new cepgen::Process::GamGamLL );
  p->kinematics.mode = cepgen::Kinematics::Mode::ElasticElastic;
  //p->kinematics.mode = cepgen::Kinematics::Mode::InelasticElastic;
  //p->kinematics.mode = cepgen::Kinematics::Mode::ElasticInelastic;
  p->kinematics.cuts.central.pt_single = 5.;
  p->kinematics.cuts.central.eta_single = { -2.5, 2.5 };
  p->kinematics.cuts.remnants.mass_single = { 1.07, 320. };

  CG_INFO( "main" ) << p;
  cepgen::Logger::get().level = cepgen::Logger::Level::debugInsideLoop;

  const unsigned short ndim = g.numDimensions();
  double x[12];
  for ( unsigned int i = 0; i < ndim; ++i )
    x[i] = 0.3;

  cout << g.computePoint( x ) << endl;

  return 0;
}
