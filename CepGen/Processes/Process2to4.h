#ifndef CepGen_Processes_Process2to4_h
#define CepGen_Processes_Process2to4_h

#include "CepGen/Processes/KTProcess.h"

namespace cepgen
{
  namespace proc
  {
    /// A 2-to-4 (or 2-to-2 central) process
    class Process2to4 : public KTProcess
    {
      public:
        /// Initialise a 2-to-4 process
        /// \param[in] params Collection of user-defined steering parameters
        /// \param[in] partons Incoming hard scattering particles
        /// \param[in] cs_id Central particles PDG id
        Process2to4( const ParametersList& params, std::array<pdgid_t,2> partons, pdgid_t cs_id );

      protected:
        /// Set all cuts for the single outgoing particle phase space definition
        void setCuts( const cuts::Central& single );

        void preparePhaseSpace() override;
        void fillCentralParticlesKinematics() override;
        double computeKTFactorisedMatrixElement() override;

        /// Conform all kinematics variables to the user-defined phase space
        virtual void prepareProcessKinematics() = 0;
        /// Computation rule for the central matrix element
        virtual double computeCentralMatrixElement() const = 0;

        //--- Mandelstam variables
        double shat() const; ///< \f$\hat s=(p_1+p_2)^2=(p_3+p_4)^2\f$
        double that() const; ///< \f$\hat t=\frac{1}{2}\left[(p_1-p_3)^2+(p_2-p_4)^2\right]\f$
        double uhat() const; ///< \f$\hat u=\frac{1}{2}\left[(p_1-p_4)^2+(p_2-p_3)^2\right]\f$

        static const Limits x_limits_; ///< Standard [0,1] limits for input variables
        ParticleProperties cs_prop_; ///< PDG id of the central particles

        cuts::Central single_limits_; ///< Limits to be applied on single central system's particles

        Momentum pA_; ///< Momentum of the positive-z incoming beam particle
        Momentum pB_; ///< Momentum of the negative-z incoming beam particle
        Momentum q1_; ///< Momentum of the first hard scattering particle
        Momentum q2_; ///< Momentum of the second hard scattering particle
        Momentum p_c1_; ///< Momentum of the first central particle
        Momentum p_c2_; ///< Momentum of the second central particle
        double y_c1_; ///< Rapidity of the first central particle
        double y_c2_; ///< Rapidity of the second central particle
        double pt_diff_; ///< Transverse momentum difference for the two central particle
        double phi_pt_diff_; ///< Azimuthal angle difference for the two central particles
        double amt1_; ///< Transverse mass of the first central particle
        double amt2_; ///< Transverse mass of the second central particle

      private:
        double ww_;
    };
  }
}

#endif

