#ifndef TEST_H
#define TEST_H

#include "../include/process.h"

class TestProcess : public Process
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

