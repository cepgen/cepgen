#ifndef CepGen_Processes_GamGamLL_h
#define CepGen_Processes_GamGamLL_h

#include "CepGen/Processes/GenericProcess.h"

namespace CepGen
{

  namespace Process
  {

    /**
     * Full class of methods and objects to compute the full analytic matrix element
     * \cite Vermaseren1983347 for the \f$\gamma\gamma\to\ell^{+}\ell^{-}\f$ process
     * according to a set of kinematic constraints provided for the incoming and
     * outgoing particles (the Kinematics object).
     * The particle roles in this process are defined as following: \n
     * \image latex lpair_kinematics.pdf Detailed particle roles in the two-photon process as defined by the @a GamGamLL object. The incoming protons/electrons are denoted by a role 1, and 2, as the outgoing protons/protons remnants/ electrons carry the indices 3 and 5. The two outgoing leptons have the roles 6 and 7, while the lepton/antilepton distinction is done randomly (thus, the arrow convention is irrelevant here).
     *
     * The \a f function created by this Process child has its \a _ndim -dimensional
     * coordinates mapped as :
     * - 0 = \f$t_1\f$, first incoming photon's virtuality
     * - 1 = \f$t_2\f$, second incoming photon's virtuality
     * - 2 = \f$s_2\f$ mapping
     * - 3 = yy4 = \f$\cos\left(\pi x_3\right)\f$ definition
     * - 4 = \f$w_4\f$, the two-photon system's invariant mass
     * - 5 = xx6 = \f$\frac{1}{2}\left(1-\cos\theta^\text{CM}_6\right)\f$ definition (3D rotation of the first outgoing lepton with respect to the two-photon centre-of-mass system). If the @a nm_ optimisation flag is set this angle coefficient value becomes
     *   \f[\frac{1}{2}\left(\frac{a_\text{map}}{b_\text{map}}\frac{\beta-1}{\beta+1}+1\right)\f]
     *   with \f$a_\text{map}=\frac{1}{2}\left(w_4-t_1-t_2\right)\f$, \f$b_\text{map}=\frac{1}{2}\sqrt{\left(\left(w_4-t_1-t_2\right)^2-4t_1t_2\right)\left(1-4\frac{w_6}{w_4}\right)}\f$, and \f$\beta=\left(\frac{a_\text{map}+b_\text{map}}{a_\text{map}-b_\text{map}}\right)^{2x_5-1}\f$
     *   and the \a fJacobian element is scaled by a factor \f$\frac{1}{2}\frac{\left(a_\text{map}^2-b_\text{map}^2\cos^2\theta^\text{CM}_6\right)}{a_\text{map}b_\text{map}}\log\left(\frac{a_\text{map}+b_\text{map}}{a_\text{map}-b_\text{map}}\right)\f$
     * - 6 = _phicm6_, or \f$\phi_6^\text{CM}\f$ the rotation angle of the dilepton system in the centre-of-mass
     *   system
     * - 7 = \f$x_q\f$, \f$w_X\f$ mappings, as used in the single- and double-dissociative
     *   cases only
     * \brief Compute the matrix element for a CE \f$\gamma\gamma\to\ell^{+}\ell^{-}\f$
     *  process
     */
    class GamGamLL : public GenericProcess
    {
      public:
        /// Class constructor ; set the mandatory parameters before integration and events generation
        /// \param[in] nopt Optimisation (legacy from LPAIR)
        explicit GamGamLL( const ParametersList& params = ParametersList() );
        ProcessPtr clone() const override { return ProcessPtr( new GamGamLL( *this ) ); }

        void addEventContent() override;
        void beforeComputeWeight() override;
        /// Compute the process' weight for the given point
        /// \return \f$\mathrm d\sigma(\mathbf x)(\gamma\gamma\to\ell^{+}\ell^{-})\f$,
        ///   the differential cross-section for the given point in the phase space.
        double computeWeight() override;
        unsigned int numDimensions() const override;
        void setKinematics( const Kinematics& cuts ) override;
        void fillKinematics( bool ) override;
        /// Compute the ougoing proton remnant mass
        /// \param[in] x A random number (between 0 and 1)
        /// \param[in] outmass The maximal outgoing particles' invariant mass
        /// \param[in] lepmass The outgoing leptons' mass
        /// \param[out] dw The size of the integration bin
        /// \return Mass of the outgoing proton remnant
        double computeOutgoingPrimaryParticlesMasses( double x, double outmass, double lepmass, double& dw );
        /// Set all the kinematic variables for the outgoing proton remnants, and prepare the hadronisation
        /// \param[in] part_ Particle to "prepare" for the hadronisation to be performed
        void prepareHadronisation( Particle *part_ );

      private:
        /**
         * Calculate energies and momenta of the
         *  1st, 2nd (resp. the "proton-like" and the "electron-like" incoming particles),
         *  3rd (the "proton-like" outgoing particle),
         *  4th (the two-photons central system), and
         *  5th (the "electron-like" outgoing particle) particles in the overall centre-of-mass frame.
         * \brief Energies/momenta computation for the various particles, in the CM system
         * \return Success state of the operation
         */
        bool orient();
        /**
         * Compute the expression of the matrix element squared for the \f$\gamma\gamma\rightarrow\ell^{+}\ell^{-}\f$ process.
         * It returns the value of the convolution of the form factor or structure functions with the central two-photons matrix element squared.
         * \brief Computes the matrix element squared for the requested process
         * \return Full matrix element for the two-photon production of a pair of spin\f$-\frac{1}{2}-\f$point particles.
         *  It is noted as \f[
         *  M = \frac{1}{4bt_1 t_2}\sum_{i=1}^2\sum_{j=1}^2 u_i v_j t_{ij} = \frac{1}{4}\frac{u_1 v_1 t_{11}+u_2 v_1 t_{21}+u_1 v_2 t_{12}+u_2 v_2 t_{22}}{t_1 t_2 b}
         * \f] where \f$b\f$ = \a bb_ is defined in \a ComputeWeight as : \f[
         *  b = t_1 t_2+\left(w_{\gamma\gamma}\sin^2{\theta^\text{CM}_6}+4m_\ell\cos^2{\theta^\text{CM}_6}\right) p_g^2
         * \f]
         */
        double periPP( int, int );
        /**
         * Describe the kinematics of the process \f$p_1+p_2\to p_3+p_4+p_5\f$ in terms of Lorentz-invariant variables.
         * These variables (along with others) will then be fed into the \a PeriPP method (thus are essential for the evaluation of the full matrix element).
         * \return Success state of the operation
         */
        bool pickin();

        /// Internal switch for the optimised code version (LPAIR legacy ; unimplemented here)
        int n_opt_;
        int pair_;

        Limits w_limits_;
        Limits q2_limits_;
        Limits mx_limits_;
        struct Masses
        {
          Masses();
          /// squared mass of the first proton-like outgoing particle
          double MX2_;
          /// squared mass of the second proton-like outgoing particle
          double MY2_;
          /// squared mass of the outgoing leptons
          double Ml2_;
          /// \f$\delta_2=m_1^2-m_2^2\f$ as defined in Vermaseren's paper
          /// \cite Vermaseren1983347 for the full definition of this quantity
          double w12_;

          /// \f$\delta_1=m_3^2-m_1^2\f$ as defined in Vermaseren's paper
          /// \cite Vermaseren1983347 for the full definition of this quantity
          double w31_;
          double dw31_;
          /// \f$\delta_4=m_5^2-m_2^2\f$ as defined in Vermaseren's paper
          /// \cite Vermaseren1983347 for the full definition of this quantity
          double w52_;
          double dw52_;
        };
        Masses masses_;

        /// energy of the first proton-like incoming particle
        double ep1_;
        /// energy of the second proton-like incoming particle
        double ep2_;
        double p_cm_;

        /// energy of the two-photon central system
        double ec4_;
        /// 3-momentum norm of the two-photon central system
        double pc4_;
        /// mass of the two-photon central system
        double mc4_;
        /// squared mass of the two-photon central system
        double w4_;

        /// \f$p_{12} = \frac{1}{2}\left(s-m_{p_1}^2-m_{p_2}^2\right)\f$
        double p12_;
        double p1k2_, p2k1_;
        /// \f$p_{13} = -\frac{1}{2}\left(t_1-m_{p_1}^2-m_{p_3}^2\right)\f$
        double p13_;
        double p14_, p25_;

        double q1dq_, q1dq2_;

        double s1_, s2_;

        double epsi_;
        double g5_, g6_;
        double a5_, a6_;
        double bb_;

        double gram_;
        double dd1_, dd2_, dd3_;
        /// \f$\delta_5=m_4^2-t_1\f$ as defined in Vermaseren's paper
        /// \cite Vermaseren1983347 for the full definition of this quantity
        double dd4_;
        double dd5_;
        /**
         * Invariant used to tame divergences in the matrix element computation. It is defined as
         * \f[\Delta = \left(p_1\cdot p_2\right)\left(q_1\cdot q_2\right)-\left(p_1\cdot q_2\right)\left(p_2\cdot q_1\right)\f]
         * with \f$p_i, q_i\f$ the 4-momenta associated to the incoming proton-like particle and to the photon emitted from it.
         */
        double delta_;
        double g4_;
        double sa1_, sa2_;

        double sl1_;

        /// cosine of the polar angle for the two-photons centre-of-mass system
        double cos_theta4_;
        /// sine of the polar angle for the two-photons centre-of-mass system
        double sin_theta4_;

        double al4_;
        double be4_;
        double de3_, de5_;
        double pt4_;

        /// Kinematics of the first incoming proton
        Particle::Momentum p1_lab_;
        /// Kinematics of the second incoming proton
        Particle::Momentum p2_lab_;
        /// Kinematics of the first outgoing proton
        Particle::Momentum p3_lab_;
        /// Kinematics of the two-photon system (in the two-proton CM)
        Particle::Momentum p4_lab_;
        /// Kinematics of the second outgoing proton
        Particle::Momentum p5_lab_;
        /// Kinematics of the first outgoing lepton (in the two-proton CM)
        Particle::Momentum p6_cm_;
        /// Kinematics of the second outgoing lepton (in the two-proton CM)
        Particle::Momentum p7_cm_;
        double jacobian_;

      private:
        /**
         * Define modified variables of integration to avoid peaks integrations (see @cite Vermaseren1983347 for details)
         * Return a set of two modified variables of integration to maintain the stability of the integrant. These two new variables are :
         * - \f$y_{out} = x_{min}\left(\frac{x_{max}}{x_{min}}\right)^{exp}\f$ the new variable
         * - \f$\mathrm dy_{out} = x_{min}\left(\frac{x_{max}}{x_{min}}\right)^{exp}\log\frac{x_{min}}{x_{max}}\f$, the new variable's differential form
         * @brief Redefine the variables of integration in order to avoid the strong peaking of the integrant.
         * @param[in] expo Exponant
         * @param[in] xmin Minimal value of the variable
         * @param[in] xmax Maximal value of the variable
         * @param[out] out The new variable definition
         * @param[out] dout The differential variant of the new variable definition
         * @param[in] var_name The variable name
         * @note This method overrides the set of `mapxx` subroutines in ILPAIR, with a slight difference according to the sign of the
         *  \f$\mathrm dy_{out}\f$ parameter :
         *  - left unchanged :
         * > `mapw2`, `mapxq`, `mapwx`, `maps2`
         *  - opposite sign :
         * > `mapt1`, `mapt2`
         */
        void map( double expo, const Limits& lim, double& out, double& dout, const std::string& var_name="" );
        void mapla( double y, double z, int u, double xm, double xp, double& x, double& d );
        /// Compute the electric/magnetic form factors for the two considered \f$Q^{2}\f$ momenta transfers
        void formFactors( double q1, double q2, FormFactors& fp1, FormFactors& fp2 ) const;
    };
  }
}

#endif

