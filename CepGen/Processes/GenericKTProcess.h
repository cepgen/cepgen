#ifndef GenericKTProcess_h
#define GenericKTProcess_h

#include "GenericProcess.h"
#include "CepGen/Physics/FormFactors.h"
#include "CepGen/Physics/PhotonFluxes.h"

namespace CepGen
{
  namespace Process
  {
    /**
     * A generic kT-factorisation process.
     * First 4 dimensions of the phase space are required for the incoming partons' virtualities (radial and azimuthal coordinates)
     * \brief Class template to define any kT-factorisation process
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date Apr 2016
     */
    class GenericKTProcess : public GenericProcess
    {
      public:
        /**
         * \brief Class constructor
         * \param[in] name Human-readable kT-factorised process name
         * \param[in] num_user_dimensions_ Number of additional dimensions required for the user process
         * \param[in] ip1 First incoming parton
         * \param[in] ip2 Second incoming parton (if undefined, same as the first)
         * \param[in] op1 First produced final state particle
         * \param[in] op2 Second produced final state particle (if undefined, same as the first)
         */
        GenericKTProcess( const std::string& name="<generic process>",
                          const unsigned int& num_user_dimensions=0,
                          const Particle::ParticleCode& ip1=Particle::Photon,
                          const Particle::ParticleCode& op1=Particle::Muon,
                          const Particle::ParticleCode& ip2=Particle::invalidParticle,
                          const Particle::ParticleCode& op2=Particle::invalidParticle);
        ~GenericKTProcess();

        /// Populate the event content with the generated process' topology
        void addEventContent();
        /// Retrieve the total number of dimensions on which the integration is being performet
        /// \param[in] proc_mode_ Kinematics case considered
        unsigned int numDimensions( const Kinematics::ProcessMode& proc_mode_ ) const;
        /// Retrieve the event weight in the phase space
        double computeWeight();
        /// Populate the event content with the generated process' kinematics  
        void fillKinematics( bool );

      protected:
        /// Set the kinematics associated to the phase space definition
        inline void setKinematics( const Kinematics& kin ) {
          cuts_ = kin;
          log_qmin_ = -10.; // FIXME //log_qmin_ = std::log( std::sqrt( cuts_.q2min ) );
          log_qmax_ = log( cuts_.qt_max );
        }
        /// Set the kinematics of the central system before any point computation
        inline virtual void prepareKTKinematics() { DebuggingInsideLoop("Dummy kinematics prepared!"); }
        /// Minimal Jacobian weight of the point considering a kT factorisation
        double minimalJacobian() const;
        /// Jacobian weight of the point in the phase space for integration
        virtual double computeJacobian() = 0;
        /// kT-factorised matrix element (event weight)
        /// \return Weight of the point in the phase space to the integral
        virtual double computeKTFactorisedMatrixElement() = 0;
        /// Compute the invariant masses of the outgoing protons (or remnants)
        void computeOutgoingPrimaryParticlesMasses();
        void computeIncomingFluxes( double, double, double, double );
        /// Set the kinematics of the incoming and outgoing protons (or remnants)
        void fillPrimaryParticlesKinematics();
        /// Set the kinematics of the outgoing central system
        virtual void fillCentralParticlesKinematics() = 0;
  
        /// Minimal log-virtuality of the intermediate parton
        double log_qmin_;
        /// Maximal log-virtuality of the intermediate parton
        double log_qmax_;
        /// Virtuality of the first intermediate parton (photon, pomeron, ...)
        double qt1_;
        /// Azimuthal rotation of the first intermediate parton's transverse virtuality
        double phi_qt1_;
        /// Virtuality of the second intermediate parton (photon, pomeron, ...)
        double qt2_;
        /// Azimuthal rotation of the second intermediate parton's transverse virtuality
        double phi_qt2_;

        /// First outgoing proton
        Particle::Momentum PX_;
        /// Second outgoing proton
        Particle::Momentum PY_;
        /// First incoming parton's flux
        double flux1_;
        /// Second incoming parton's flux
        double flux2_;

      private:
        void addPartonContent();
        static const unsigned int kNumRequiredDimensions = 4;
        /// Number of additional dimensions required for the user process
        /// (in addition to the 4 required for the two partons' transverse momenta)
        unsigned int kNumUserDimensions;
        /// First intermediate parton (photon, pomeron, ...)
        Particle::ParticleCode kIntermediatePart1;
        /// Second intermediate parton (photon, pomeron, ...)
        Particle::ParticleCode kIntermediatePart2;
        /// Type of particle produced in the final state
        Particle::ParticleCode kProducedPart1;
        /// Type of particle produced in the final state
        Particle::ParticleCode kProducedPart2;
    };
  }
}

#endif
