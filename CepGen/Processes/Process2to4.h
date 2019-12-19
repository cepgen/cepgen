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
        void setCuts( const Cuts& single );

        void preparePhaseSpace() override;
        void fillCentralParticlesKinematics() override;
        double computeKTFactorisedMatrixElement() override;

        virtual void prepareProcessKinematics() = 0;
        virtual double computeCentralMatrixElement() const = 0;

        //--- Mandelstam variables
        double shat() const;
        double that() const;
        double uhat() const;

        static const Limits x_limits_;
        ParticleProperties cs_prop_; ///< PDG id of the central particles

        Cuts single_limits_;

        Momentum pA_, pB_, q1_, q2_;
        Momentum p_c1_; ///< Momentum of the first central particle
        Momentum p_c2_; ///< Momentum of the second central particle
        double y_c1_; ///< Rapidity of the first central particle
        double y_c2_; ///< Rapidity of the second central particle
        double pt_diff_; ///< Transverse momentum difference for the two central particle
        double phi_pt_diff_; ///< Azimuthal angle difference for the two central particles
        double amt1_, amt2_;

      private:
        double ww_;
    };
  }
}

#endif

