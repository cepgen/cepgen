#ifndef CepGen_Processes_GenericKTProcess_h
#define CepGen_Processes_GenericKTProcess_h

#include "GenericProcess.h"
#include "CepGen/Physics/FormFactors.h"

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
         * \param[in] name Generic process name
         * \param[in] description Human-readable kT-factorised process name
         * \param[in] num_user_dimensions Number of additional dimensions required for the user process
         * \param[in] partons First and second incoming parton
         * \param[in] output Produced final state particles
         */
        GenericKTProcess( const std::string& name,
                          const std::string& description,
                          const unsigned int& num_user_dimensions,
                          const std::array<Particle::ParticleCode,2>& partons,
                          const std::vector<Particle::ParticleCode>& output );
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
        void setKinematics( const Kinematics& kin );
        /// Set the kinematics of the central system before any point computation
        virtual void setExtraContent() {}
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
        /// Compute the unintegrated photon fluxes (for inelastic distributions, interpolation on double logarithmic grid)
        void computeIncomingFluxes( double, double, double, double );
        /// Set the kinematics of the incoming and outgoing protons (or remnants)
        void fillPrimaryParticlesKinematics();
        /// Set the kinematics of the outgoing central system
        virtual void fillCentralParticlesKinematics() = 0;


        /// Get the elastic flux to be expected at a given x_bjorken / kT
        /// \param[in] x Bjorken x
        /// \param[in] kt2 Transverse 2-momentum \f$\mathbf{q}_{\mathrm{T}}^2\f$ of the incoming photon
        double elasticFlux( double x, double kt2 );

        /// Get the inelastic flux to be expected at a given x_bjorken / kT
        /// \param[in] x Bjorken x
        /// \param[in] kt2 Transverse 2-momentum \f$\mathbf{q}_{\mathrm{T}}^2\f$ of the incoming photon
        /// \param[in] mx Outgoing diffractive proton mass
        double inelasticFlux( double x, double kt2, double mx );

        /// Retrieve a component of the phase space point for the kT-factorised process
        inline double xkt( const unsigned int i ) const { return x( kNumRequiredDimensions + i ); }
  
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
        static constexpr unsigned int kNumRequiredDimensions = 4;
        /// Number of additional dimensions required for the user process
        /// (in addition to the 4 required for the two partons' transverse momenta)
        unsigned int kNumUserDimensions;
        /// First and second intermediate parton (photon, pomeron, ...)
        std::array<Particle::ParticleCode,2> kIntermediateParts;
        /// Type of particles produced in the final state
        std::vector<Particle::ParticleCode> kProducedParts;
    };
  }
}

#endif
