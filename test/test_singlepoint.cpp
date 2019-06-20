#include "CepGen/Generator.h"
#include "CepGen/Parameters.h"
#include "CepGen/Processes/ProcessesHandler.h"
#include "CepGen/Core/Exception.h"

using namespace std;

int main()
{
  cepgen::Generator gen;
  auto& params = gen.parameters();
  params.setProcess( cepgen::proc::ProcessesHandler::get().build( "lpair", cepgen::ParametersList()
    .set<int>( "mode", (int)cepgen::KinematicsMode::ElasticElastic )
  ) );
  params.kinematics.cuts.central.pt_single = 5.;
  params.kinematics.cuts.central.eta_single = { -2.5, 2.5 };
  params.kinematics.cuts.remnants.mass_single = { 1.07, 320. };

  CG_INFO( "main" ) << &params;
  cepgen::utils::Logger::get().level = cepgen::utils::Logger::Level::debugInsideLoop;

  double x[12];
  for ( unsigned int i = 0; i < gen.numDimensions(); ++i )
    x[i] = 0.3;

  cout << gen.computePoint( x ) << endl;

  return 0;
}
