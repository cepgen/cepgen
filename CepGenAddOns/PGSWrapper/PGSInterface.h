/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2019-2022  Laurent Forthomme
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CepGenAddOns_EventInterfaces_PGSInterface_h
#define CepGenAddOns_EventInterfaces_PGSInterface_h

namespace pgs {
  enum object { photon = 0, electron = 1, muon = 2, tau_had = 3, jet = 4, heavy_charged = 5 };
  enum printmask : int { hepevt = 1 << 0, calo_sum = 1 << 1, calo_clus = 1 << 2, trg_obj = 1 << 3, off_obj = 1 << 4 };
}  // namespace pgs

extern "C" {
void pgs_user_args_();
//void pgs_user_init_();
void pgs_initialize_();
void pgs_load_param_();
bool pgs_next_event_();
void pgs_dump_event_(int&, double&, double&);
void pgs_trigger_();
void pgs_recon_();
void pgs_write_event_(const char* cmd, int cmd_size);

// callback routines
void pgs_user_event_(int& done) { done = 0; }
void pgs_user_herwig_() {}
void pgs_user_pythia_() {}

/// PGS event header and control information
extern struct {
  int numarg;                 ///< number of arguments supplied to program
  char pgs_args[10][40];      ///< list of arguments (char*40)
  int nevpgs;                 ///< number of events to generate/read
  int target_lum;             ///< target luminosity (in pb-1)
  int nprpgs;                 ///< number of events to print out
  int pgs_iseed, pgs_jseed;   ///< seeds for pgs_ranmar
  int pgs_log_unit;           ///< log file unit number
  char optpgs[6];             ///< type of run: 'PYTHIA', 'ISAJET', 'FILE', ...
  char evtlum[6];             ///< number of events ('events') or luminosity ('pb-1')
  char pgs_input_file[80];    ///< input file
  char pgs_output_file[80];   ///< output file
  char pgs_log_file[80];      ///< log file
  char pgs_param_file[80];    ///< detector parameter file
  char pgs_isajet_decay[80];  ///< ISAJET decay table file name
  char pgs_isajet_cards[80];  ///< ISAJET card file name
  char pgs_pythia_cards[80];  ///< PYTHIA card file name
  int pgs_herwig_proc;        ///< HERWIG process to generate
  char pgs_herwig_susy[80];   ///< HERWIG SUSY data file
  char pgs_alpgen_stem[80];   ///< ALPGEN unweighted events file stem
} pgsevt_;
const int ntrkmx = 500;

/// PGS track list
extern struct {
  int numtrk;  ///< number of tracks
  int dumtrk;
  int indtrk[ntrkmx];      ///< index to HEPEVT particle
  double ptrk[ntrkmx][3];  ///< track 3-vector
  double qtrk[ntrkmx];     ///< track charge
} pgstrk_;
const int nphimax = 600, netamax = 600;

/// PGS calorimeter tower arrays
extern struct {
  double ecal[nphimax][netamax];  ///< electromagnetic energy in each tower (phi,eta)
  double hcal[nphimax][netamax];  ///< hadronic energy in each tower
  double met_cal;                 ///< calorimeter missing ET
  double phi_met_cal;             ///< calorimeter missing ET phi
  double met_cor;                 ///< missing ET corrected for muons
  double phi_met_cor;             ///< corrected missing ET phi
} pgscal_;
const int nmxobj = 500;

/// PGS reconstructed object candidate list
extern struct {
  int numobj;  ///< number of reconstructed objects
  int dumobj;
  int indobj[nmxobj];          ///< index to HEPEVT particle (where relevant)
  pgs::object typobj[nmxobj];  ///< reconstructed type
  double pobj[nmxobj][4];      ///< four vector of reconstructed object
  double qobj[nmxobj];         ///< charge of reconstructed object
  double vecobj[nmxobj][10];   ///< interesting object quantities (see below)
  bool unique[nmxobj];         ///< true for object if it is uniquely identified and passes cuts in pgs_object_cuts
} pgsrec_;
const int nmxhep = 4000;
extern struct {
  int nevhep;              ///< event number
  int nhep;                ///< number of particles in event
  int isthep[nmxhep];      ///< particle status
  int idhep[nmxhep];       ///< particle PDG id
  int jmohep[nmxhep][2];   ///< particle parents
  int jdahep[nmxhep][2];   ///< particle secondary products
  double phep[nmxhep][5];  ///< particle 4-momentum and mass
  double vhep[nmxhep][4];  ///< particle production 4-vector
} hepevt_;

// disable all HERWIG6 routines
void hwbgen_() {}
void hwcdec_() {}
void hwcfor_() {}
void hwdhad_() {}
void hwdhob_() {}
void hwdhvy_() {}
void hwefin_() {}
void hweini_() {}
void hwepro_() {}
void hwigin_() {}
void hwissp_() {}
void hwmevt_() {}
void hwufne_() {}
void hwuinc_() {}
void hwuine_() {}
void hwusta_() {}
void hwwarn_() {}
// disable all TAUOLA routines
void tauola_init_() {}
// disable all PYTHIA6 routines
void pyevnt_() {}
void pygive_() {}
void pyinit_() {}
void pylist_() {}
void lunhep_() {}
// disable all other routines
void stdflisxsec_() {}
void stdflpyxsec_() {}
void stdchg_() {}
void stdxropen_() {}
void stdxrd_() {}
void stdxwinit_() {}
void stdxwrt_() {}
void stdxend_() {}
void upveto_() {}
void pdgrdtb_() {}
}

#endif
