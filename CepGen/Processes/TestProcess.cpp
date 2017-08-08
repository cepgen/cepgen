#include "TestProcess.h"

TestProcess::TestProcess() : GenericProcess("<test process>")
{}

TestProcess::~TestProcess()
{}

unsigned int
TestProcess::GetNdim( const Kinematics::ProcessMode& ) const
{
  return 3;
}

double
TestProcess::ComputeWeight()
{
  
  double A = 1./(M_PI*M_PI*M_PI);
  return A/(1.-cos(x(0)*M_PI)*cos(x(1)*M_PI)*cos(x(2)*M_PI));
}

void
TestProcess::FillKinematics(bool)
{
  return;
}
