/**
 * Class template to define any process to compute using this MC integrator/generator
 * @author Laurent Forthomme <laurent.forthomme@uclouvain.be>
 * @date January 2014
 */
class Process
{
 public:
  Process();
  ~Process();
  /**
   * @brief Returns the weight for this point in the phase-space
   */
  inline double ComputeWeight() { return -1.; };
 private:
  double a;
};
