#include "test.h"

TestProcess::TestProcess() : Process("<test process>")
{}

TestProcess::~TestProcess()
{}

int
TestProcess::GetNdim(int process_mode_) const
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
TestProcess::FillKinematics()
{
  return;
}
