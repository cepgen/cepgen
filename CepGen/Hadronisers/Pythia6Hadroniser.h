#ifndef CepGen_Hadronisers_Pythia6Hadroniser_h
#define CepGen_Hadronisers_Pythia6Hadroniser_h

#include <algorithm>

#include "GenericHadroniser.h"

/// Maximal number of characters to fetch for the particle's name
#define NAME_CHR 16

extern "C"
{
  /// Get the particle's mass in GeV from the Pythia6 module
  extern double pymass_( int& );
  /// Get the resonant particle's width in GeV from the Pythia6 module
  //extern double pywidt_( int& );
  /// Launch the Pythia6 fragmentation
  extern void pyexec_();
  /// Set a parameter value to the Pythia6 module
  extern void pygive_( const char*, int );
  extern void pyckbd_();
  /// List all the particles in the event in a human-readable format
  extern void pylist_( int& );
  /// Join two coloured particles in a colour singlet
  extern void pyjoin_( int&, int& );
  /// Get a particle's human-readable name from the Pythia6 module
  extern void pyname_( int&, char*, int );
  /// Get kinematic information on a particle from the Pythia6 module
  extern double pyp_( int&, int& );
  /// Store one parton/particle in the PYJETS common block
  extern void py1ent_( int&, int&, double&, double&, double& );

  /// Particles content of the event
  extern struct
  {
    /// Number of particles in the event
    int n;
    int npad;
    /// Particles' general information (status, PDG id, mother, daughter 1, daughter 2)
    int k[5][4000];
    /// Particles' kinematics, in GeV (px, py, pz, E, M)
    double p[5][4000];
    /// Primary vertex for the particles
    double v[5][4000];
  } pyjets_;
}

namespace CepGen
{
  namespace Hadroniser
  {
    /**
     * Full interface to the Pythia6 @cite Sjostrand:2006za algorithm. It can be used in a single particle decay mode as well as a full event hadronisation using the string model, as in Jetset.
     * @brief Pythia6 hadronisation algorithm
     */
    class Pythia6Hadroniser : public GenericHadroniser
    {
      public:
        Pythia6Hadroniser();
        ~Pythia6Hadroniser();

        bool hadronise( const Particle* );
        bool hadronise( Event* );

      private:
        inline static double pymass(int pdgid_) { return pymass_(pdgid_); }
        //inline static double pywidt(int pdgid_) { return pywidt_(pdgid_); }
        inline static void pyexec() { pyexec_(); }
        inline static void pyckbd() { pyckbd_(); }
        inline static void pygive( const std::string& line ) { pygive_( line.c_str(), line.length() ); }
        inline static void pylist( int mlist ) { pylist_( mlist ); }
        inline static double pyp( int role, int qty ) { return pyp_( role, qty ); }
        //inline static void py1ent( int* kf, double* pe, double theta, double phi ) { int one=1; py1ent_( &one, kf, pe, theta, phi ); }
        //inline static void py1ent( int* kf, double* pe, double theta, double phi ) { py1ent_( 1, kf, pe, theta, phi ); }
        inline static std::string pyname( int pdgid ) {
          char out[NAME_CHR];
          std::string s;
          pyname_( pdgid, out, NAME_CHR );
          s = std::string( out, NAME_CHR );
          //s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
          s.erase( remove( s.begin(), s.end(), ' ' ), s.end() );
          return s;
        }
        /**
         * \brief Connect entries with colour flow information
         * \param[in] njoin_ Number of particles to join in the colour flow
         * \param[in] ijoin_ List of particles unique identifier to join in the colour flow
         */
        inline static void pyjoin( int njoin, int ijoin[2] ) { return pyjoin_( njoin, *ijoin ); }
        bool prepareHadronisation( Event* );
    };
  }
}

#endif
