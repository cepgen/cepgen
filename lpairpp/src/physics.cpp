#include "physics.h"

PhysicsBoundaries::PhysicsBoundaries() :
  wmin(20.), wmax(0.),
  q2min(4.), q2max(100.),
  zmin(0.), zmax(1.)
{}

PhysicsBoundaries::~PhysicsBoundaries()
{}

// values of a, b, c provided from the fits on ep data and retrieved from
// http://dx.doi.org/10.1016/0550-3213(76)90231-5 with 1.110 <= w2 <=1.990

double abrass[56] = {5.045,5.126,5.390,5.621,5.913,5.955,6.139,6.178,6.125,5.999,
                     5.769,5.622,5.431,5.288,5.175,5.131,5.003,5.065,5.045,5.078,
                     5.145,5.156,5.234,5.298,5.371,5.457,5.543,5.519,5.465,5.384,
                     5.341,5.320,5.275,5.290,5.330,5.375,5.428,5.478,5.443,5.390,
                     5.333,5.296,5.223,5.159,5.146,5.143,5.125,5.158,5.159,5.178,
                     5.182,5.195,5.160,5.195,5.163,5.172};
double bbrass[56] = {0.798,1.052,1.213,1.334,1.397,1.727,1.750,1.878,1.887,1.927,
                     2.041,2.089,2.148,2.205,2.344,2.324,2.535,2.464,2.564,2.610,
                     2.609,2.678,2.771,2.890,2.982,3.157,3.183,3.315,3.375,3.450,
                     3.477,3.471,3.554,3.633,3.695,3.804,3.900,4.047,4.290,4.519,
                     4.709,4.757,4.840,5.017,5.015,5.129,5.285,5.322,5.545,5.623,
                     5.775,5.894,6.138,6.151,6.301,6.542};
double cbrass[56] = { 0.043, 0.024, 0.000,-0.013,-0.023,-0.069,-0.060,-0.080,-0.065,-0.056,
                     -0.065,-0.056,-0.043,-0.034,-0.054,-0.018,-0.046,-0.015,-0.029,-0.048,
                     -0.032,-0.045,-0.084,-0.115,-0.105,-0.159,-0.164,-0.181,-0.203,-0.223,
                     -0.245,-0.254,-0.239,-0.302,-0.299,-0.318,-0.383,-0.393,-0.466,-0.588,
                     -0.622,-0.568,-0.574,-0.727,-0.665,-0.704,-0.856,-0.798,-1.048,-0.980,
                     -1.021,-1.092,-1.313,-1.341,-1.266,-1.473};

bool PSF(double q2_, double mX2_, double* sigT_, double* w1_, double* w2_)
{
  int nBin;
  double xBin, dx, nu2, logqq0, gd2;
  double sigLow, sigHigh;
  double mX = std::sqrt(mX2_);
  double mP = GetMassFromPDGId(2212);
  double mPI = 0.135; //FIXME pi0 mass ???

  if (mX>=mP+mPI && mX<1.99) {
    if (mX<1.11) {
      nBin = 0;
      xBin = mX-mP-mPI;
      dx = 1.11-mP-mPI; // Delta w bin sizes
    }
    else if (mX<1.77) { // w in [1.11, 1.77[
      dx = 0.015; // Delta w bin sizes
      nBin = (mX-1.11)/dx+1;
      xBin = fmod(mX-1.11, dx);
    }
    else { // w in [1.77, 1.99[
      dx = 0.02; // Delta w bin sizes
      nBin = (mX-1.77)/dx+45;
      xBin = fmod(mX-1.77, dx);
    }
  }
  else {
    *sigT_ = 0.;
    *w1_ = 0.;
    *w2_ = 0.;
    return false;
  }
  nu2 = std::pow((mX2_-q2_-std::pow(mP, 2))/(2.*mP), 2);
  logqq0 = log((nu2-q2_)/std::pow((mX2_-std::pow(mP, 2))/(2.*mP), 2))/2.;
  gd2 = std::pow(1./(1-q2_/.71), 4); // dipole form factor of the proton

  sigLow = (nBin!=0) ?
    exp(abrass[nBin-1]+bbrass[nBin-1]*logqq0+cbrass[nBin-1]*std::pow(fabs(logqq0), 3))*gd2 :
    0.;
  sigHigh =
    exp(abrass[nBin]  +bbrass[nBin]  *logqq0+cbrass[nBin]  *std::pow(fabs(logqq0), 3))*gd2;

  *sigT_ = sigLow+xBin*(sigHigh-sigLow)/dx;
  *w1_ = (mX2_-std::pow(mP, 2))/(8.*std::pow(pi, 2)*mP*alphaF)*muBarn*(*sigT_);
  *w2_ = (*w1_)*q2_/(q2_-nu2);

  return true;
}

Particles EPA(Particle el_, Particle pr_, int mode_, PhysicsBoundaries b_)
{
  Particles op;
  el_.id = 0;
  el_.role = 1;
  el_.Dump();
  op.push_back(el_);
  op[0].Dump();

  int isum, irnd;
  double psum;
  double dsum, qsum, dsumt, qsumt, dsuml, qsuml;
  double r, epamax, epa, epat, epal, eparho;
  double sthe, cthe, phi;
  double eqe, q2, w, lf;
  double pesc, eesc;
  double elpr, s, esmp2, eel, wmin2, w12, y, ysqr;
  double dymax, dymin;
  double q2min, q2max;
  double emsqr, emqe2;
  double emy, exy;
  double epsil, ftrans;
  int ierr1, ierr2;

  ierr1 = ierr2 = 0;

  bool _gephot_first = true;

  if (_gephot_first) {
    isum = irnd = 0;
    dsum = qsum = dsumt = qsumt = dsuml = qsuml = 0.;
    w12 = dymin = 0.;
    dymax = 1.;

    // Calculate CMS s=(P+K)**2

    psum = 0.;
    for (int i=0; i<3; i++) psum+= pr_.P(i)*el_.P(i);
    elpr = pr_.E()*el_.E()-psum; // 4-vector product

    esmp2 = std::pow(2.*elpr+el_.M2(), 2);
    s = (el_+pr_).E2();
    //Particle tot = el_+pr_;

    // Evaluate photon flux in proton rest frame: set EEL to approx. 50TeV

    if (mode_>3) eel = elpr/pr_.M();
    else eel = el_.E();

    wmin2 = std::pow(b_.wmin, 2);
    w12 = wmin2-pr_.M2();

    // Calculate Y - bounds from
    // ALI, A. et al. (1987): Heavy quark physics at HERA. -
    //   Proc. HERA workshop, Hamburg 1987 (ed. R.D. PECCEI), 395-494.

    ysqr = std::sqrt(std::pow(s-w12, 2)-4.*w12*el_.M2());
    dymax = (s+w12+ysqr)/(2.*(s+el_.M2()));

    // Use trick for quadratic equations. See
    // W.H. PRESS et al. (1988): Numerical Recipes in C. Cambridge (Cambridge
    // Univ. Press), p. 156.

    dymin = std::max(w12/(dymax*(s+el_.M2())), b_.zmin);

    // Calculate absolute maximum of y, irrespective of final state

    dymax = std::min(s/(s+el_.M2()), b_.zmax);
    std::cout << "dymax = " << s/(s+el_.M2()) << " <-> " << b_.zmax << std::endl;
    dymax = std::min((std::pow(b_.wmax, 2)-pr_.M2()+b_.q2max)/(2.*elpr), dymax);
    std::cout << "final dymax = " << dymax << std::endl;

    // Set max. photon weight for efficient rejection plane
    q2min = std::pow(el_.M()+dymin, 2)/(1.-dymin);
    if (q2min<b_.q2min) q2min = b_.q2min;
    q2max = dymax*s;
    std::cout << "q2max = " << q2max << ", boundary = " << b_.q2max << std::endl;
    if (q2max>b_.q2max) q2max = b_.q2max;

    if (mode_==1) { // WWA - approximation
      epamax = alphared*(4.*(1.-dymin)+std::pow(dymin, 2));
      std::cout << "alphared = " << alphared << ", dymin = " << dymin << " -> epamax = " << epamax << std::endl;
    }
    else { // Full transversal spectrum (2) or full longitudinal and transversal (3) spectrum

      eqe = q2min/std::pow(eel, 2);
      emqe2 = std::pow(dymin-eqe/4., 2);
      emsqr = (std::pow(dymin*elpr, 2)+q2min*pr_.M2())/(std::pow(elpr, 2)+el_.M2()*pr_.M2());

      if (emsqr<0.) {
	std::cerr << "[EPA] ERROR: problem with sqrt(emsqr)=" << emsqr << " at EPAMAX determination." << std::endl;
	exit(0);
      }

      if (mode_==2) {
	// Transversal spectrum
	epamax = alphared*dymin*std::sqrt(emsqr)*(2.*(1.-dymin)+emqe2+eqe)/(emqe2+eqe);
      }
      else { // Longitudinal & transversal spectrum
	epamax = alphared*dymin*std::sqrt(emsqr)/(emqe2+eqe)*(4.*(1.-dymin)+emqe2+eqe);
      }
    }
    std::cout << "dymax = " << dymax << ", dymin = " << dymin << ", q2max = " << q2max << ", q2min = " << q2min << std::endl;
    epamax *= log(dymax/dymin)*log(q2max/q2min);
    std::cout << "mode = " << mode_ << ", epamax = " << epamax << std::endl;

    _gephot_first = false;
  }

  // Reset rndloop counter every new event
  irnd = 0;

  // Begin main loop over Y,Q2 random production
  do {
    do {
      do {
	isum++; irnd++;
	// Produce Y spect. ( 1/y weighted shape )
	y = dymin*(std::pow(dymax/dymin, drand()));
	ysqr = std::pow(y, 2);
	
	// Calculate actual Q2_min, Q2_max from Y
	q2min = el_.M2()*ysqr/(1.-y);
	q2max = y*s;
	
	// Take Q2_cut from steering, if it is kinematicly reachable.
	if (q2min<b_.q2min) q2min = b_.q2min;
	if (q2max>b_.q2max) q2max = b_.q2max;
	//isum++;
      } while (q2min>q2max);
      
      // Produce Q2 spect. (1/x weighted shape )
      q2 = q2min*std::pow(q2max/q2min, drand());
      
      //---------------------------------------------------------------
      // EPA - WWA spectrum
      //---------------------------------------------------------------
      
      // Calc. photon weight
      if (mode_==1) { // WWA - approximation
	r = alphared/(y*q2);
	epat = r*(2.*(1.-y)*(1.-el_.M2()*ysqr/(1.-y)*q2))+ysqr;
	epal = r*2.*(1.-y);
      }
      else {
	// Full transversal spectrum (2) or full longitudinal and transversal (3)
	// spectrum from:
	// ABT, I. & J.R. SMITH (1992): MC Upgrades to Study Untagged Events. -
	//    H1-10/92-249.
	// See also:
	// SMITH, J.R. (1992): An Experimentalist's Guide to Photon Flux
	//    Calculations. - H1-12/92-259.
	// SMITH, J.R. (1993): Polarization Decomposition of Fluxes and
	//    Kinematics in ep Reactions. - H1-04/93-282.
	
	eqe = q2/std::pow(eel, 2);
	emqe2 = std::pow(y-eqe/4., 2);
	emsqr = (std::pow(y*elpr, 2)+q2*pr_.M2())/(std::pow(elpr, 2)+el_.M2()*pr_.M2());
	
	if (emsqr<0.) {
	  std::cerr << "[EPA] WARNING: problem with sqrt(emsqr)= " << emsqr << ": y, Q2 pair rejected" << std::endl;
	  //CALL ERRLOG (280, 'S: GEPHOT: EMSQR<0!')
	  ierr1++;
	  if (ierr1>10) {
	    std::cerr << "[EPA] ERROR: too many sqrt problems: try WWA" << std::endl;
	    //CALL ERRLOG (281, 'F: GEPHOT: EMSQR<0 too often!')
	    exit(0);
	  }
	}
	
	if (mode_==2) { // Transversal spectrum
	  epat = alphared/q2*std::sqrt(emsqr)*(2.*(1.-y)+emqe2+eqe)/(emqe2+eqe);
	  epal = 0.;
	}
	else { // Longitudinal & transversal spectrum
	  r = alphared/q2*std::sqrt(emsqr)/(emqe2+eqe);
	  epat = r*(2.*(1.-y)+emqe2+eqe);
	  epat = r*2.*(1.-y);
	}
      }
      
      epa = epat+epal;
      lf = epat/epa;
      
      // Unweight MC
      r = y*q2*log(dymax/dymin)*log(q2max/q2min);
      w = std::sqrt(y*2.*elpr-q2+pr_.M2());
      // Check if W > W_min, else reject photon
      if (w<b_.wmin) r = 0.;
      // Check if W < W_max, else reject photon
      if (w>b_.wmax) r = 0.;
      epa *= r;
      epat *= r;
      epal *= r;
      
      // Update upper EPA bound
      if (epa>epamax) {
	if (epa>1.1*epamax) std::cout << "[EPA] INFO: EPA > 1.1*EPAMAX !" << std::endl;
	else if (epa>1.01*epamax) std::cout << "[EPA] INFO: EPA > 1.01*EPAMAX !" << std::endl;
	else std::cout << "[EPA] INFO: EPA > EPAMAX !" << std::endl;
	epamax = epa;
	std::cout << "[EPA] INFO: update of maximal weight : " << epamax << std::endl;
      }
      
      // Global counter for overall integration
      dsum += epa;
      qsum += std::pow(epa, 2);
      dsumt += epat;
      qsumt += std::pow(epat, 2);
      dsuml += epal;
      qsuml += std::pow(epal, 2);
      
      if (irnd>10000) { // Kin. loop failed
	std::cerr << "[EPA] ERROR: Kinematic loop failed after " << irnd << " trials." << std::endl
		  << "  EPAMAX too high for efficient mc! EPAMAX=" << epamax << std::endl;
	//CALL ERRLOG (285, 'F: GEPHOT: More than 10000 iterations!')
	exit(0);
      }
      
      // End rnd loop: rejection method
      eparho = drand()*epamax;
      
    } while (eparho>epa);
    
    // Continue with unweighted kinematics
    
    // Scattering angle of electron in LAB.: E.Lohrmann DESY HERA 83/08
    // x = Q2 / (y s)
    // E_sc = E_e(1-y) + E_p x y
    // cos t = [E_e(1-y) - E_p x y] / E_sc
    
    emy = el_.E()*(1.-y);
    exy = pr_.E()*q2/s;
    
    eesc = emy+exy;
    cthe = (emy-exy)/eesc;
    sthe = 2.*std::sqrt(emy*exy)/eesc;
    
    // Control scattering angle 
    ierr2++;

    if (ierr2>100) {
      std::cerr << "[EPA] ERROR: too many problems for CTHE or STHE:" << std::endl
		<< "  CTHE=" << cthe << ", STHE=" << sthe << std::endl;
    }
  } while (fabs(cthe)>1. or fabs(sthe)>1.);

  phi = 2.*pi*drand();

  pesc = -sqrt(std::pow(eesc, 2)-el_.M2());
  Particle outele(2, op[0].pdgId);
  outele.P(pesc*sthe*cos(phi),
	   pesc*sthe*sin(phi),
	   pesc*cthe,
	   eesc);
  outele.id = op.size();
  outele.SetMother(&(op[0]));
  op.push_back(outele);

  Particle outgam = op[0]-op[1];
  outgam.role = 3;
  outgam.pdgId = 22;
  outgam.helicity = Heli(lf);
  outgam.id = op.size();
  outgam.SetMother(&(op[0]));
  op.push_back(outgam);

  ftrans = epat;
  epsil = epal/epat;

  return op;
}
