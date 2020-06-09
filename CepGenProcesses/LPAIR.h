#ifndef CepGenProcesses_LPAIR_h
#define CepGenProcesses_LPAIR_h

#include "CepGen/Processes/Process.h"

namespace cepgen
{
  namespace proc
  {
    /**
     * Full class of methods and objects to compute the full analytic matrix element
     * \cite Vermaseren:1982cz for the \f$\gamma\gamma\to\ell^{+}\ell^{-}\f$ process
     * according to a set of kinematic constraints provided for the incoming and
     * outgoing particles (the Kinematics object).
     * The \a f function created by this Process child has its \a _ndim -dimensional
     * coordinates mapped as :
     * - 0 = \f$t_1\f$, first incoming photon's virtuality
     * - 1 = \f$t_2\f$, second incoming photon's virtuality
     * - 2 = \f$s_2\f$ mapping
     * - 3 = yy4 = \f$\cos\left(\pi x_3\right)\f$ definition
     * - 4 = \f$w_4\f$, the two-photon system's invariant mass
     * - 5 = xx6 = \f$\frac{1}{2}\left(1-\cos\theta^{\rm CM}_6\right)\f$ definition (3D rotation of the first outgoing lepton with respect to the two-photon centre-of-mass system). If the \a nm_ optimisation flag is set this angle coefficient value becomes
     *   \f[\frac{1}{2}\left(\frac{a_{\rm map}}{b_{\rm map}}\frac{\beta-1}{\beta+1}+1\right)\f]
     *   with \f$a_{\rm map}=\frac{1}{2}\left(w_4-t_1-t_2\right)\f$, \f$b_{\rm map}=\frac{1}{2}\sqrt{\left(\left(w_4-t_1-t_2\right)^2-4t_1t_2\right)\left(1-4\frac{w_6}{w_4}\right)}\f$, and \f$\beta=\left(\frac{a_{\rm map}+b_{\rm map}}{a_{\rm map}-b_{\rm map}}\right)^{2x_5-1}\f$
     *   and the Jacobian element is scaled by a factor \f$\frac{1}{2}\frac{\left(a_{\rm map}^2-b_{\rm map}^2\cos^2\theta^{\rm CM}_6\right)}{a_{\rm map}b_{\rm map}}\log\left(\frac{a_{\rm map}+b_{\rm map}}{a_{\rm map}-b_{\rm map}}\right)\f$
     * - 6 = _phicm6_, or \f$\phi_6^{\rm CM}\f$ the rotation angle of the dilepton system in the centre-of-mass
     *   system
     * - 7 = \f$x_q\f$, \f$w_X\f$ mappings, as used in the single- and double-dissociative
     *   cases only
     * \brief Compute the matrix element for a CE \f$\gamma\gamma\to\ell^{+}\ell^{-}\f$
     *  process
     */
    class LPAIR : public Process
    {
      public:
        /// \brief Class constructor: set the mandatory parameters before integration and events generation
        /// \param[in] params General process parameters (nopt = Optimisation, legacy from LPAIR)
        explicit LPAIR( const ParametersList& params = ParametersList() );
        ProcessPtr clone( const ParametersList& params ) const override {
          return ProcessPtr( new LPAIR( *this ) );
        }

        void addEventContent() override;
        double computeWeight() override;
        void prepareKinematics() override;
        void fillKinematics( bool ) override;

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
         *  b = t_1 t_2+\left(w_{\gamma\gamma}\sin^2{\theta^{\rm CM}_6}+4m_\ell\cos^2{\theta^{\rm CM}_6}\right) p_g^2
         * \f]
         */
        double periPP() const;
        /**
         * Describe the kinematics of the process \f$p_1+p_2\to p_3+p_4+p_5\f$ in terms of Lorentz-invariant variables.
         * These variables (along with others) will then be fed into the \a PeriPP method (thus are essential for the evaluation of the full matrix element).
         * \return Success state of the operation
         */
        bool pickin();

        /// Internal switch for the optimised code version (LPAIR legacy ; unimplemented here)
        const int n_opt_;
        pdgid_t pair_;
        const bool symmetrise_;

        std::vector<double> x_tmp_;
        double theta4_;
        double phi6_cm_;
        double x6_;

        Limits w_limits_;
        struct Masses
        {
          /// squared mass of the outgoing leptons
          double Ml2 = 0.;
          /// \f$\delta_2=m_1^2-m_2^2\f$ as defined in \cite Vermaseren:1982cz
          double w12 = 0.;
          /// \f$\delta_1=m_3^2-m_1^2\f$ as defined in \cite Vermaseren:1982cz
          double w31 = 0.;
          /// \f$\delta_4=m_5^2-m_2^2\f$ as defined in \cite Vermaseren:1982cz
          double w52 = 0.;
        } masses_;

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
        /// \cite Vermaseren:1982cz for the full definition of this quantity
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
        Momentum p1_lab_;
        /// Kinematics of the second incoming proton
        Momentum p2_lab_;
        /// Kinematics of the first outgoing proton
        Momentum p3_lab_;
        /// Kinematics of the two-photon system (in the two-proton CM)
        Momentum p4_lab_;
        /// Kinematics of the second outgoing proton
        Momentum p5_lab_;
        /// Kinematics of the first outgoing lepton (in the two-proton CM)
        Momentum p6_cm_;
        /// Kinematics of the second outgoing lepton (in the two-proton CM)
        Momentum p7_cm_;
        double jacobian_;

      private:
        /**
         * Define modified variables of integration to avoid peaks integrations (see \cite Vermaseren:1982cz for details)
         * Return a set of two modified variables of integration to maintain the stability of the integrant. These two new variables are :
         * - \f$y_{out} = x_{min}\left(\frac{x_{max}}{x_{min}}\right)^{exp}\f$ the new variable
         * - \f$\mathrm dy_{out} = x_{min}\left(\frac{x_{max}}{x_{min}}\right)^{exp}\log\frac{x_{min}}{x_{max}}\f$, the new variable's differential form
         * \brief Redefine the variables of integration in order to avoid the strong peaking of the integrant
         * \param[in] expo Exponant
         * \param[in] lim Min/maximal value of the variable
         * \param[in] var_name The variable name
         * \return A pair containing the value and the bin width the new variable definition
         * \note This method overrides the set of `mapxx` subroutines in ILPAIR, with a slight difference according to the sign of the
         *  \f$\mathrm dy_{out}\f$ parameter :
         *  - left unchanged :
         * > `mapw2`, `mapxq`, `mapwx`, `maps2`
         *  - opposite sign :
         * > `mapt1`, `mapt2`
         */
        std::pair<double,double> map( double expo, const Limits& lim, const std::string& var_name = "" );
        std::pair<double,double> mapla( double y, double z, int u, const Limits& lim );
    };
  }
}

#endif
