#ifndef Test_h
#define Test_h

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

