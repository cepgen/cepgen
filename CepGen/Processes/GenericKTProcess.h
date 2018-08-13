#ifndef CepGen_Processes_GenericKTProcess_h
#define CepGen_Processes_GenericKTProcess_h

#include "GenericProcess.h"

namespace CepGen
{
  namespace Process
  {
    /**
     * A generic kT-factorisation process.
     * \note
     * - First 4 dimensions of the phase space are required for the
     *    incoming partons' virtualities (radial and azimuthal coordinates).
     * - Last 0-2 dimensions may be used for the scattered diffractive
     *    system(s)' invariant mass definition.
     * \brief Class template to define any kT-factorisation process
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date Apr 2016
     */
    class GenericKTProcess : public GenericProcess
    {
      public:
        /// Class constructor
        /// \param[in] name Generic process name
        /// \param[in] description Human-readable kT-factorised process name
        /// \param[in] partons First and second incoming parton
        /// \param[in] output Produced final state particles
        GenericKTProcess( const std::string& name,
                          const std::string& description,
                          const std::array<PDG,2>& partons,
                          const std::vector<PDG>& output );

        /// Populate the event content with the generated process' topology
        void addEventContent() override;
        /// Retrieve the total number of dimensions on which the integration is being performet
        unsigned int numDimensions( const Kinematics::Mode& ) const override;
        /// Retrieve the event weight in the phase space
        double computeWeight() override;
        /// Populate the event content with the generated process' kinematics
        void fillKinematics( bool ) override;

        void dumpVariables() const;
        unsigned short computationMethod() const { return method_; }
        void setComputationMethod( unsigned short i ) { method_ = i; }

        /// Type of incoming partons fluxes
        enum class Flux
        {
          Elastic = 0,
          ElasticBudnev = 10,
          Inelastic = 1,
          InelasticBudnev = 11
        };
        friend std::ostream& operator<<( std::ostream&, const Flux& );
        /// Get the flux at a given parton x/kT
        /// \param[in] kt2 Transverse 2-momentum \f$\mathbf{q}_{\mathrm{T}}^2\f$ of the incoming parton
        /// \param[in] x Parton momentum fraction
        /// \param[in] mx Outgoing diffractive proton mass
        static double flux( const Flux& type, double x, double kt2, StructureFunctions& sf, double mx = 0. );

      protected:
        /// Set the kinematics associated to the phase space definition
        void setKinematics( const Kinematics& kin ) override;
        /// Set the kinematics of the central system before any point computation
        virtual void setExtraContent() {}
        /// Prepare the central part of the Jacobian (only done once, as soon as the kinematics is set)
        virtual void preparePhaseSpace() = 0;
        /// kT-factorised matrix element (event weight)
        /// \return Weight of the point in the phase space to the integral
        virtual double computeKTFactorisedMatrixElement() = 0;
        /// Compute the unintegrated photon fluxes (for inelastic distributions, interpolation on double logarithmic grid)
        void computeIncomingFluxes( double, double, double, double );
        /// Set the kinematics of the incoming and outgoing protons (or remnants)
        void fillPrimaryParticlesKinematics();
        /// Set the kinematics of the outgoing central system
        virtual void fillCentralParticlesKinematics() = 0;

        /// Type of mapping to apply on the variable
        enum class Mapping
        {
          /// a linear \f$\textrm dx\f$ mapping
          linear = 0,
          /// a logarithmic \f$\frac{\textrm dx}{x} = \textrm d(\log x)\f$ mapping
          logarithmic,
          /// a square \f$\textrm dx^2=2x\cdot\textrm dx\f$ mapping
          square
        };
        friend std::ostream& operator<<( std::ostream&, const Mapping& );
        /// Register a variable to be handled and populated whenever
        ///  a new phase space point weight is to be calculated.
        /// \note To be run once per generation (before any point computation)
        /// \param[out] out Reference to the variable to be mapped
        /// \param[in] type Type of mapping to apply
        /// \param[in] in Integration limits
        /// \param[in] default_limits Limits to apply if none retrieved from the user configuration
        /// \param[in] description Human-readable description of the variable
        void registerVariable( double& out, const Mapping& type, const Limits& in, Limits default_limits, const char* description );
        /// Generate and initialise all variables handled by this process
        /// \return Phase space point-dependent component of the Jacobian weight of the point in the phase space for integration
        /// \note To be run at each point computation (therefore, to be optimised!)
        double generateVariables() const;

        unsigned short num_dimensions_;
        unsigned short method_;

        /// Phase space point-independant component of the Jacobian weight of the point in the phase space for integration
        double kt_jacobian_;

        /// Log-virtuality range of the intermediate parton
        Limits log_qt_limits_;
        /// Intermediate azimuthal angle range
        Limits phi_qt_limits_;
        /// Invariant mass range for the scattered excited system
        Limits mx_limits_;

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

        /// Handler to a variable mapped by this process
        struct MappingVariable
        {
          /// Human-readable description of the variable
          std::string description;
          /// Kinematic limits to apply on the variable
          Limits limits;
          /// Reference to the process variable to generate/map
          double& variable;
          /// Interpolation type
          Mapping type;
          /// Corresponding integration variable
          unsigned short index;
        };
        /// Collection of variables to be mapped at the weight generation stage
        std::vector<MappingVariable> mapped_variables_;

      private:
        static const double kMinFlux;
        /// First and second intermediate parton (photon, pomeron, ...)
        std::array<PDG,2> kIntermediateParts;
        /// Type of particles produced in the final state
        std::vector<PDG> kProducedParts;
    };
  }
}

#endif

