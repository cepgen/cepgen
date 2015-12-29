#include "PPtoLL.h"

PPtoLL::PPtoLL() : GenericProcess("gamma,gamma->l+,l- (kT-factorization approach)")
{}

PPtoLL::~PPtoLL()
{}

void
PPtoLL::AddEventContent()
{
  IncomingState is; OutgoingState os;
  is.insert(ParticleWithRole(Particle::IncomingBeam1,    Particle::Proton));
  is.insert(ParticleWithRole(Particle::IncomingBeam2,    Particle::Proton));
  is.insert(ParticleWithRole(Particle::Parton1,          Particle::Photon));
  is.insert(ParticleWithRole(Particle::Parton2,          Particle::Photon));
  os.insert(ParticleWithRole(Particle::OutgoingBeam1,    Particle::Proton));
  os.insert(ParticleWithRole(Particle::OutgoingBeam2,    Particle::Proton));
  os.insert(ParticleWithRole(Particle::CentralParticle1, Particle::Muon));
  os.insert(ParticleWithRole(Particle::CentralParticle2, Particle::Muon));
  GenericProcess::SetEventContent(is, os);
}

int
PPtoLL::GetNdim(ProcessMode process_mode_) const
{
  switch (process_mode_) {
    default:
    case ElasticElastic:     return 8;
    case ElasticInelastic:
    case InelasticElastic:   return 9;
    case InelasticInelastic: return 10;
  }
}

double
PPtoLL::ComputeWeight()
{
  double lqmin, lqmax;
  double jac, weight;
  
  //lqmin = std::log(std::sqrt(fCuts.q2min));
  lqmin = -10.; //FIXME
  lqmax = std::log(std::sqrt(fCuts.q2max));
  
  // Incoming photons
  _q1t = std::exp(lqmin+(lqmax-lqmin)*x(0));
  _q2t = std::exp(lqmin+(lqmax-lqmin)*x(1));
  _phiq1t = 2.*pi*x(2);
  _phiq2t = 2.*pi*x(3);
  
  // Outgoing leptons
  //const double ymin = EtaToY(fCuts.etamin, GetParticle(Particle::CentralParticle1)->M(), pt),
  //             ymax = EtaToY(fCuts.etamax);
  ////////////////////////////////////
  const double ymin = fCuts.etamin, //
               ymax = fCuts.etamax; //
  ///////////// FIXME ////////////////
  
  _y1 = ymin+(ymax-ymin)*x(4);
  _y2 = ymin+(ymax-ymin)*x(5);
  _ptdiff = fCuts.ptdiffmin+(fCuts.ptdiffmax-fCuts.ptdiffmin)*x(6);
  _phiptdiff = 2.*pi*x(7);

  // Outgoing protons (or remnants)
  switch (fCuts.kinematics) {
    case 0: default: { Error("PPtoLL is intended for p-on-p collisions! Aborting!"); exit(0); break; }
    case 1: 
      _mx = GetParticle(Particle::IncomingBeam1)->M();
      _my = GetParticle(Particle::IncomingBeam2)->M();
      break;
    case 2:
      _mx = fCuts.mxmin+(fCuts.mxmax-fCuts.mxmin)*x(8);
      _my = GetParticle(Particle::IncomingBeam2)->M();
      break;
    case 3:
      _mx = GetParticle(Particle::IncomingBeam1)->M();
      _my = fCuts.mxmin+(fCuts.mxmax-fCuts.mxmin)*x(8);
      break;
    case 4:
      _mx = fCuts.mxmin+(fCuts.mxmax-fCuts.mxmin)*x(8);
      _my = fCuts.mxmin+(fCuts.mxmax-fCuts.mxmin)*x(9);
      break;
  }
  
  // Jacobian computation
  jac = 1.;
  jac *= (lqmax-lqmin)*_q1t; // d(q1t) . q1t
  jac *= (lqmax-lqmin)*_q2t; // d(q2t) . q2t
  jac *= 2.*pi; // d(phi1)
  jac *= 2.*pi; // d(phi2)
  jac *= (ymax-ymin); // d(y1)
  jac *= (ymax-ymin); // d(y2)
  if (fCuts.kinematics>1) {
    jac *= (fCuts.mxmax-fCuts.mxmin); // d(mx/y)
    if (fCuts.kinematics>3) {
      jac *= (fCuts.mxmax-fCuts.mxmin); // d(my/x)
    }
  }
  jac *= (fCuts.ptdiffmax-fCuts.ptdiffmin); // d(Dpt)
  jac *= 2.*pi; // d(phiDpt)
  
  weight = jac*INCqqbar();
  
  return weight;
}

double
PPtoLL::INCqqbar()
{
  int iterm11, iterm12, itermtt, iterm22;
  int imethod, imat1, imat2;
  int idif, idely;
  double pdif, dely_min, dely_max;
  const double alpha_em = 1./137.035;
  const double mp = Particle::GetMassFromPDGId(Particle::Proton), mp2 = pow(mp, 2);
  const double ml = GetParticle(Particle::CentralParticle1)->M(), ml2 = pow(ml, 2);
  double units;

  iterm11 = 1; // Long-long
  iterm22 = 1; // Trans-trans
  iterm12 = 1; // Long-trans
  itermtt = 1; // Trans-trans(')

  //=================================================================
  //     How matrix element is calculated
  //         imethod = 0: on-shell formula
  //         imethod = 1: off-shell formula
  //=================================================================
  imethod = 1;
  //imethod = 0;
  
  //=================================================================
  //     two terms in Wolfgang's formula for 
  //     off-shell gamma gamma --> l^+ l^-
  //=================================================================
  imat1 = 2;
  imat2 = 0;

  //=================================================================
  //     extra cuts on the p1t(l) and p2t(l) plane
  //         idif = 0: no extra cut
  //         idif = 1: extra cut
  //=================================================================
  idif = 0; // 0 is a standard
  pdif = 2.5;
  
  //=================================================================
  //     the distance in rapidity between l^+ and l^-
  //=================================================================
  idely = 0; // 0 or 1 
  dely_min = 4.0;
  dely_max = 5.0;
  
  //=================================================================
  //     conversion factor
  //     1/GeV^2 --> pb
  //=================================================================
  units = 1.e4*pow(197.3271, 2);
  
  //=================================================================
  //     matrix element computation
  //=================================================================
  double q1tx, q1ty, q2tx, q2ty;
  double ptsumx, ptsumy, ptsum;
  double ptdiffx, ptdiffy;
  double pt1x, pt1y, pt1, pt2x, pt2y, pt2;
  double amt1, amt2;
  double dely;
  double alpha1, alpha2, beta1, beta2;
  double q1t, q1t2, q2t, q2t2;
  double z1p, z1m, z2p, z2m;
  double f1, f2;
  
  //const double stild = fS/2.*(1+sqrt(1.-(4*pow(mp2, 2))/pow(fS, 2)));
  
  // Inner photons
  q1tx = _q1t*cos(_phiq1t);
  q1ty = _q1t*sin(_phiq1t);
  q2tx = _q2t*cos(_phiq2t);
  q2ty = _q2t*sin(_phiq2t);
  
  // Two-photon system
  ptsumx = q1tx+q2tx;
  ptsumy = q1ty+q2ty;
  ptsum = sqrt(pow(ptsumx, 2)+pow(ptsumy, 2));
  
  ptdiffx = _ptdiff*cos(_phiptdiff);
  ptdiffy = _ptdiff*sin(_phiptdiff);
  
  // Outgoing leptons
  pt1x = (ptsumx+ptdiffx)/2.;
  pt1y = (ptsumy+ptdiffy)/2.;
  pt1 = sqrt(pow(pt1x, 2)+pow(pt1y, 2));
  
  pt2x = (ptsumx-ptdiffx)/2.;
  pt2y = (ptsumy-ptdiffy)/2.;
  pt2 = sqrt(pow(pt2x, 2)+pow(pt2y, 2));
  
  if (pt1<fCuts.ptmin or pt2<fCuts.ptmin) {
    return 0.;
  }
  /*if (pt1>fCuts.ptmax or pt2>fCuts.ptmax) {
    return 0.;
  }*/
  amt1 = sqrt(pow(pt1, 2)+ml2);
  amt2 = sqrt(pow(pt2, 2)+ml2);
  
  //=================================================================
  //     a window in transverse momentum difference
  //=================================================================
  if (idif==1 and fabs(pt1-pt2)>pdif) return 0.;

  //const double pcaptx = pt1x+pt2x, pcapty = pt1y+pt2y;
  dely = fabs(_y1-_y2);
  //=================================================================
  //     a window in rapidity distance
  //=================================================================
  if (idely==1 and (dely<dely_min or dely>dely_max)) return 0.;
  
  //=================================================================
  //     auxiliary quantities
  //=================================================================

  alpha1 = amt1/fSqS*exp( _y1);
  alpha2 = amt2/fSqS*exp( _y2);
  beta1  = amt1/fSqS*exp(-_y1);
  beta2  = amt2/fSqS*exp(-_y2);

  q1t2 = pow(q1tx, 2)+pow(q1ty, 2);
  q2t2 = pow(q2tx, 2)+pow(q2ty, 2);

  //x2 = 0.; //FIXME figure out where this comes from
  //const double delta_x1 = (pow(_mx, 2)+q2t2)/((1.-x2)*fS);

  //x1 = alpha1+alpha2+delta_x1;
  const double x1 = alpha1+alpha2;
  const double x2 = beta1+beta2;

  /*const double xi_x1 = log10(x1);
  const double xi_x2 = log10(x2);*/

  z1p = alpha1/x1;
  z1m = alpha2/x1;
  z2p = beta1/x2;
  z2m = beta2/x2;
  
  if (x1>1. or x2>1.) return 0.; // sanity check

  double px_plus, px_minus;
  double py_plus, py_minus;
  double p10, p1x, p1y, p1z, p20, p2x, p2y, p2z;
  
  // FIXME FIXME FIXME
  const double ak10 = GetParticle(Particle::IncomingBeam1)->E(),
               ak1z = fabs(GetParticle(Particle::IncomingBeam1)->GetMomentum().Pz()),
               ak20 = GetParticle(Particle::IncomingBeam2)->E(),
               ak2z = fabs(GetParticle(Particle::IncomingBeam2)->GetMomentum().Pz());
  
  //=================================================================
  //     additional conditions for energy-momentum conservation
  //=================================================================
  
  const double s1_eff = x1*fS-pow(_q1t,2), s2_eff = x2*fS-pow(_q2t,2);
  const double invm = sqrt(pow(amt1,2)+pow(amt2,2)+2.*amt1*amt2*cosh(_y1-_y2)-pow(ptsum,2));
  
  switch (fCuts.kinematics) {
    case 2:
      if (sqrt(s1_eff)<=(_my+invm)) return 0.;
    case 3:
      if (sqrt(s2_eff)<=(_mx+invm)) return 0.;
    case 4:
      if (sqrt(s1_eff)<=(_my+invm)) return 0.;
      if (sqrt(s2_eff)<=(_mx+invm)) return 0.;
  }
  
  //const double qcaptx = pcaptx;
  //const double qcapty = pcapty;
  
  //=================================================================
  //     four-momenta of the outgoing protons (or remnants)
  //=================================================================

  px_plus = (1.-x1)*ak1z*sqrt(2.);
  px_minus = (pow(_mx, 2)+pow(q1tx, 2)+pow(q1ty, 2))/2./px_plus;
      
  _px_0 = (px_plus+px_minus)/sqrt(2.);
  _px_z = (px_plus-px_minus)/sqrt(2.);
  _px_x = -q1tx;
  _px_y = -q1ty;
      
  py_minus = (1.-x2)*ak2z*sqrt(2.); // warning! sign of pz??
  py_plus = (pow(_my, 2)+pow(q2tx, 2)+pow(q2ty, 2))/2./py_minus;
      
  _py_0 = (py_plus+py_minus)/sqrt(2.);
  _py_z = (py_plus-py_minus)/sqrt(2.);
  _py_x = -q2tx;
  _py_y = -q2ty;
  
  q1t = sqrt(q1t2);
  q2t = sqrt(q2t2);
  
  //=================================================================
  //     four-momenta of the outgoing l^+ and l^-
  //=================================================================

  p10 = alpha1*ak10+beta1*ak20;
  p1x = pt1x;
  p1y = pt1y;
  p1z = alpha1*ak1z+beta1*ak2z;

  _pl1_0 = sqrt(pow(pt1, 2)+ml2)*cosh(_y1);
  _pl1_x = pt1x;
  _pl1_y = pt1y;
  _pl1_z = sqrt(pow(pt1, 2)+ml2)*sinh(_y1);

  p20 = alpha2*ak10+beta2*ak20;
  p2x = pt2x;
  p2y = pt2y;
  p2z = alpha2*ak1z+beta2*ak2z;
  
  _pl2_0 = sqrt(pow(pt2, 2)+ml2)*cosh(_y2);
  _pl2_x = pt2x;
  _pl2_y = pt2y;
  _pl2_z = sqrt(pow(pt2, 2)+ml2)*sinh(_y2);

  //=================================================================
  //     four-momenta squared of the virtual photons
  //=================================================================

  //FIXME FIXME FIXME ////////////
  const double q10 = ak10, q20 = ak20;
  const double q1z = ak1z, q2z = ak2z;
  ////////////////////////////////
  
  /*const double q12 = pow(q10, 2)-pow(q1tx, 2)-pow(q1ty, 2)-pow(q1z, 2);
  const double q22 = pow(q20, 2)-pow(q2tx, 2)-pow(q2ty, 2)-pow(q2z, 2);*/

  //=================================================================
  //     Mendelstam variables
  //=================================================================
  double that1, that2, that, uhat1, uhat2, uhat;
  
  //shat = fS*x1*x2; // ishat = 1 (approximation)
  //const double shat = pow(q10+q20, 2)-pow(q1tx+q2tx, 2)-pow(q1ty+q2ty, 2)-pow(q1z+q2z, 2); // ishat = 2 (exact formula)

  that1 = pow(q10-p10, 2)-pow(q1tx-p1x, 2)-pow(q1ty-p1y, 2)-pow(q1z-p1z, 2);
  uhat1 = pow(q10-p20, 2)-pow(q1tx-p2x, 2)-pow(q1ty-p2y, 2)-pow(q1z-p2z, 2);
  that2 = pow(q20-p20, 2)-pow(q2tx-p2x, 2)-pow(q2ty-p2y, 2)-pow(q2z-p2z, 2);
  uhat2 = pow(q20-p10, 2)-pow(q2tx-p1x, 2)-pow(q2ty-p1y, 2)-pow(q2z-p1z, 2);

  //const double mll = sqrt(shat);

  that = (that1+that2)/2.;
  uhat = (uhat1+uhat2)/2.;

  //=================================================================
  //     matrix elements
  //=================================================================
  double amat2 = 0.;
  if (imethod==0) {
  
    //=================================================================
    //     on-shell formula for M^2
    //=================================================================
    
    const double term1 = 6.*pow(ml, 8),
                 term2 = -3.*pow(ml, 4)*pow(that, 2),
                 term3 = -14.*pow(ml, 4)*that*uhat,
                 term4 = -3.*pow(ml, 4)*pow(uhat, 2),
                 term5 = ml2*pow(that, 3),
                 term6 = 7.*ml2*pow(that, 2)*uhat,
                 term7 = 7.*ml2*that*pow(uhat, 2),
                 term8 = ml2*pow(uhat, 3),
                 term9  = -pow(that, 3)*uhat,
                 term10 = -that*pow(uhat, 3);

    const double auxil_gamgam = -2.*(term1+term2+term3+term4+term5+term6+term7+term8+term9+term10)/(pow(ml2-that, 2)*pow(ml2-uhat, 2));
    const double g_em = sqrt(4.*pi*alpha_em);
    amat2 = pow(g_em, 4)*auxil_gamgam;
  }
  else if (imethod==1) {
  
    //=================================================================
    //     Wolfgang's formulae
    //=================================================================

    double Phi10, Phi11_x, Phi11_y, Phi102, Phi11_dot_e, Phi11_cross_e;
    double Phi20, Phi21_x, Phi21_y, Phi202, Phi21_dot_e, Phi21_cross_e;
    double aux2_1, aux2_2;

    const double ak1_x = z1m*pt1x-z1p*pt2x, ak1_y = z1m*pt1y-z1p*pt2y;
    const double ak2_x = z2m*pt1x-z2p*pt2x, ak2_y = z2m*pt1y-z2p*pt2y;

    const double t1abs = (q1t2+x1*(pow(_mx, 2)-mp2)+pow(x1, 2)*mp2)/(1.-x1);
    const double t2abs = (q2t2+x2*(pow(_my, 2)-mp2)+pow(x2, 2)*mp2)/(1.-x2);

    const double eps12 = ml2+z1p*z1m*t1abs;
    const double eps22 = ml2+z2p*z2m*t2abs;

    Phi10 = 1./(pow(ak1_x+z1p*q2tx, 2)+pow(ak1_y+z1p*q2ty, 2)+eps12)
           -1./(pow(ak1_x-z1m*q2tx, 2)+pow(ak1_y-z1m*q2ty, 2)+eps12);
    Phi11_x = (ak1_x+z1p*q2tx)/(pow(ak1_x+z1p*q2tx, 2)+pow(ak1_y+z1p*q2ty, 2)+eps12)
             -(ak1_x-z1m*q2tx)/(pow(ak1_x-z1m*q2tx, 2)+pow(ak1_y-z1m*q2ty, 2)+eps12);
    Phi11_y = (ak1_y+z1p*q2ty)/(pow(ak1_x+z1p*q2tx, 2)+pow(ak1_y+z1p*q2ty, 2)+eps12)
             -(ak1_y-z1m*q2ty)/(pow(ak1_x-z1m*q2tx, 2)+pow(ak1_y-z1m*q2ty, 2)+eps12);

    Phi102 = Phi10*Phi10;
    //const double Phi112 = pow(Phi11_x, 2)+pow(Phi11_y, 2);

    Phi20 = 1./(pow(ak2_x+z2p*q1tx, 2)+pow(ak2_y+z2p*q1ty, 2)+eps22)
           -1./(pow(ak2_x-z2m*q1tx, 2)+pow(ak2_y-z2m*q1ty, 2)+eps22);
    Phi21_x = (ak2_x+z2p*q1tx)/(pow(ak2_x+z2p*q1tx, 2)+pow(ak2_y+z2p*q1ty, 2)+eps22)
             -(ak2_x-z2m*q1tx)/(pow(ak2_x-z2m*q1tx ,2)+pow(ak2_y-z2m*q1ty, 2)+eps22);
    Phi21_y = (ak2_y+z2p*q1ty)/(pow(ak2_x+z2p*q1tx, 2)+pow(ak2_y+z2p*q1ty, 2)+eps22)
             -(ak2_y-z2m*q1ty)/(pow(ak2_x-z2m*q1tx, 2)+pow(ak2_y-z2m*q1ty, 2)+eps22);

    Phi202 = Phi20*Phi20;
    //const double Phi212 = pow(Phi21_x, 2)+pow(Phi21_y, 2);

    /*aux2_1 = iterm11*(ml2+4.*z1p*z1m*t1abs)*Phi102
            +iterm22*(pow(z1p, 2)+pow(z1m, 2))*Phi112
            -iterm12*4.*z1p*z1m*(z1p-z1m)*Phi10*(q1tx*Phi11_x+q1ty*Phi11_y);
    aux2_2 = iterm11*(ml2+4.*z2p*z2m*t2abs)*Phi202
            +iterm22*(pow(z2p, 2)+pow(z2m, 2))*Phi212
            -iterm12*4.*z2p*z2m*(z2p-z2m)*Phi20*(q2tx*Phi21_x+q2ty*Phi21_y);*/

    Phi11_dot_e = (Phi11_x*q1tx + Phi11_y*q1ty)/sqrt(q1t2);
    Phi11_cross_e = (Phi11_x*q1ty - Phi11_y*q1tx)/sqrt(q1t2);

    Phi21_dot_e = (Phi21_x*q2tx +Phi21_y*q2ty)/sqrt(q2t2);
    Phi21_cross_e = (Phi21_x*q2ty -Phi21_y*q2tx)/sqrt(q2t2);

    aux2_1 = iterm11*(ml2+4.*pow(z1p, 2)*pow(z1m, 2)*t1abs)*Phi102
            +iterm22*((pow(z1p, 2)+pow(z1m, 2))*(pow(Phi11_dot_e, 2)+pow(Phi11_cross_e, 2)))
            +itermtt*(pow(Phi11_cross_e, 2)-pow(Phi11_dot_e, 2))
            -iterm12*4.*z1p*z1m*(z1p-z1m)*Phi10*(q1tx*Phi11_x+q1ty*Phi11_y);

    aux2_2 = iterm11*(ml2+4.*pow(z2p, 2)*pow(z2m, 2)*t2abs)*Phi202
            +iterm22*((pow(z2p, 2)+pow(z2m, 2))*(pow(Phi21_dot_e, 2)+pow(Phi21_cross_e, 2)))
            +itermtt*(pow(Phi21_cross_e, 2)-pow(Phi21_dot_e, 2))
            -iterm12*4.*z2p*z2m*(z2p-z2m)*Phi20*(q2tx*Phi21_x+q2ty*Phi21_y);


    //=================================================================
    //     convention of matrix element as in our kt-factorization
    //     for heavy flavours
    //=================================================================
    double amat2_1, amat2_2;
    
    amat2_1 = pow(4.*pi*alpha_em, 2)*pow(x1*x2*fS, 2)*aux2_1*2.*z1p*z1m*t1abs/(q1t2*q2t2)*t2abs/q2t2;
    amat2_2 = pow(4.*pi*alpha_em, 2)*pow(x1*x2*fS, 2)*aux2_2*2.*z2p*z2m*t2abs/(q1t2*q2t2);

    //=================================================================
    //     symmetrization
    //=================================================================

    amat2 = (imat1*amat2_1+imat2*amat2_2)/2.;

    /*const double xx1 = alpha1+alpha2, xx2 = beta1+beta2;

    const double sudakov_2 = (pow(_mx, 2)-mp2+q2t2+xx2*mp2)/((1.-xx2)*fS);
    const double sudakov_1 = (q1t2 + xx1*mp2)/((1.-xx1)*fS);
    const double ratio1 = sudakov_1 / xx1,
                 ratio2 = sudakov_2 / xx2;*/

    //if (ratio1>0.01) return 0.;
  }

  //============================================
  //     unintegrated photon distributions
  //     interpolation on double logarithmic grid
  //     of inelastic distributions
  //============================================
  
  switch (fCuts.kinematics) {
    case 1: // elastic-elastic
      f1 = ElasticFlux(x1, q1t2);
      f2 = ElasticFlux(x2, q2t2);
      break;
    case 2: // elastic-inelastic
      f1 = ElasticFlux(x1, q1t2);
      f2 = InelasticFlux(x2, q2t2, _my);
      break;
    case 3: // inelastic-elastic
      f1 = InelasticFlux(x1, q1t2, _mx);
      f2 = ElasticFlux(x2, q2t2);
      break;
    case 4: // inelastic-inelastic
      f1 = InelasticFlux(x1, q1t2, _mx);
      f2 = InelasticFlux(x2, q2t2, _my);
      break;
  }
  if (f1<1.e-20) f1 = 0.;
  if (f2<1.e-20) f2 = 0.;
  
  //=================================================================
  //     factor 2.*pi below from integration over phi_sum
  //     factor 1/4 below from jacobian of transformations
  //     factors 1/pi and 1/pi due to integration
  //     over d^2 kappa_1 d^2 kappa_2 instead d kappa_1^2 d kappa_2^2
  //=================================================================

  const double aintegral = (2.*pi)*1./(16.*pow(pi, 2)*pow(x1*x2*fS, 2)) * amat2
                         * f1/pi*f2/pi*(1./4.)*units
                         * 0.5*4./(4.*pi);

  //=================================================================
  return aintegral*q1t*q2t*_ptdiff;
  //=================================================================

  /*
  //=================================================================
  //     four-momenta squared of virtual photons
  //=================================================================
  q12 = pow(q10, 2)-pow(q1tx, 2)-pow(q1ty, 2)-pow(q1z, 2);
  q22 = pow(q20, 2)-pow(q2tx, 2)-pow(q2ty, 2)-pow(q2z, 2);
  */
}

void
PPtoLL::FillKinematics(bool)
{
  //=================================================================
  //     outgoing protons
  //=================================================================
  Particle *op1 = GetParticle(Particle::OutgoingBeam1),
           *op2 = GetParticle(Particle::OutgoingBeam2);
  switch (fCuts.kinematics) {
    case 1:
      op1->status = Particle::FinalState;
      op2->status = Particle::FinalState;
      break;
    case 2:
      op1->status = Particle::Undecayed; op1->SetM();
      op2->status = Particle::FinalState;
      break;
    case 3:
      op1->status = Particle::FinalState;
      op2->status = Particle::Undecayed; op2->SetM();
      break;
    case 4:
      op1->status = Particle::Undecayed; op1->SetM();
      op2->status = Particle::Undecayed; op2->SetM();
      break;    
  }
  
  if (!op1->SetMomentum(_px_x, _px_y, _px_z, _px_0)) {
    Error(Form("Invalid outgoing proton 1: energy: %.2f", _px_0));
  }
  
  if (!op2->SetMomentum(_py_x, _py_y, _py_z, _py_0)) {
    Error(Form("Invalid outgoing proton 2: energy: %.2f", _py_0));
  }

  // randomise the charge of the outgoing leptons
  int sign = (drand()>.5) ? +1 : -1;

  //=================================================================
  //     first outgoing lepton
  //=================================================================
  Particle* ol1 = GetParticle(Particle::CentralParticle1);
  ol1->SetPDGId(ol1->GetPDGId(), sign);
  if (!ol1->SetMomentum(_pl1_x, _pl1_y, _pl1_z, _pl1_0)) {
    Error("Invalid outgoing lepton 1");
  }

  //=================================================================
  //     second outgoing lepton
  //=================================================================
  Particle* ol2 = GetParticle(Particle::CentralParticle2);
  ol2->SetPDGId(ol2->GetPDGId(), -sign);
  if (!ol2->SetMomentum(_pl2_x, _pl2_y, _pl2_z, _pl2_0)) {
    Error("Invalid outgoing lepton 2");
  }
}
