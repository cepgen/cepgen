#ifndef TEST_H
#define TEST_H

#include "../include/GenericProcess.h"

class TestProcess : public GenericProcess
{
 public:
  TestProcess();
  ~TestProcess();
  
  int GetNdim(int) const;
  double ComputeWeight();
  void FillKinematics();
  
 private:

};

#endif

