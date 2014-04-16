#ifndef _GAMGAMWW_H
#define _GAMGAMWW_H

#include "process.h"

/**
 * @brief Computes the matrix element for a CE \f$\gamma\gamma\rightarrow W^+W^-\f$ process
 */
class GamGamWW : public Process
{
 public:
  GamGamWW();
  ~GamGamWW();
  bool SetIncomingParticles(Particle, Particle);
  bool SetOutgoingParticles(int, int);
  void FillKinematics(bool);
  void SetKinematics(Kinematics);
  void ComputeCMenergy();
  double ComputeMX(double x_, double outmass_, double* dw_);
  double ComputeWeight();
  void StoreEvent(std::ofstream*,double);
  void PrepareHadronisation(Particle *part_);
 private:
  /** @brief \f$s\f$, squared centre of mass energy of the incoming particles' system */
  double _s;
  /** @brief \f$\sqrt{s}\f$, centre of mass energy of the incoming particles' system */
  double _sqs;
  /** @brief Total energy provided by the two incoming proton-like particles */
  double _etot;
  /** @brief Total momentum provided by the two incoming proton-like particles (along the \f$z\f$-axis) */
  double _ptot;

  /** @brief Is the first incoming proton-like particle's kinematic set ? */
  bool setp1;
  /** @brief Is the second incoming proton-like particle's kinematic set ? */
  bool setp2;
  /** @brief Is the first outgoing proton-like particle's kinematic set ? */
  bool setp3;
  /** @brief Is the second outgoing proton-like particle's kinematic set ? */
  bool setp5;
  /** @brief Is the outgoing leptons' state set ? */
  bool setll;

};

#endif

