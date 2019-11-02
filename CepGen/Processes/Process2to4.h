#ifndef CepGen_Processes_Process2to4_h
#define CepGen_Processes_Process2to4_h

#include "CepGen/Processes/KTProcess.h"

namespace cepgen
{
  namespace proc
  {
    class Process2to4 : public KTProcess
    {
      public:
        Process2to4( const ParametersList& params, const std::string& name, const std::string& desc, std::array<pdgid_t,2> partons, pdgid_t cs_id );

      protected:
        void setKinematics( const Kinematics& kin ) override;
        void setCuts( const Cuts& single );

        void preparePhaseSpace() override;
        void fillCentralParticlesKinematics() override;
        double computeKTFactorisedMatrixElement() override;

        virtual void prepareProcessKinematics() = 0;
        virtual double computeCentralMatrixElement() const = 0;

        ParticleProperties cs_prop_; ///< PDG id of the central particles

        Cuts single_limits_;

        Momentum p1_, p2_, q1_, q2_;
        Momentum p_c1_; ///< Momentum of the first central particle
        Momentum p_c2_; ///< Momentum of the second central particle
        double y_c1_; ///< Rapidity of the first central particle
        double y_c2_; ///< Rapidity of the second central particle
        double pt_diff_; ///< Transverse momentum difference for the two central particle
        double phi_pt_diff_; ///< Azimuthal angle difference for the two central particles

      private:
        double ww_;
    };
  }
}

#endif

