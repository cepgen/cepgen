#include "PPtoLL.h"
#include <assert.h>

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
  lqmax = std::log(fCuts.qtmax);

  // Incoming photons
  _q1t = std::exp(lqmin+(lqmax-lqmin)*x(0));
  _q2t = std::exp(lqmin+(lqmax-lqmin)*x(1));
  _phiq1t = 2.*pi*x(2);
  _phiq2t = 2.*pi*x(3);
  DebugInsideLoop(Form("photons transverse virtualities:\n\t  mag = %f / %f (%.2f < log(qt) < %.2f)\n\t  phi = %f / %f", _q1t, _q2t, lqmin, lqmax, _phiq1t, _phiq2t));
  
  // Outgoing leptons
  //const double ymin = EtaToY(fCuts.etamin, GetParticle(Particle::CentralParticle1)->M(), pt),
  //             ymax = EtaToY(fCuts.etamax);
  ////////////////////////////////////
  const double ymin = fCuts.etamin, //
               ymax = fCuts.etamax; //
  ///////////// FIXME ////////////////
  
  _y1 = ymin+(ymax-ymin)*x(4);
  _y2 = ymin+(ymax-ymin)*x(5);
  DebugInsideLoop(Form("leptons rapidities (%.2f < y < %.2f): %f / %f", ymin, ymax, _y1, _y2));
 
  if (fCuts.ptdiffmax<0.) fCuts.ptdiffmax = 400.; //FIXME
  _ptdiff = fCuts.ptdiffmin+(fCuts.ptdiffmax-fCuts.ptdiffmin)*x(6);
  _phiptdiff = 2.*pi*x(7);
  DebugInsideLoop(Form("leptons pt difference:\n\t  mag = %f (%.2f < Dpt < %.2f)\n\t  phi = %f", _ptdiff, fCuts.ptdiffmin, fCuts.ptdiffmax, _phiptdiff));

  // Outgoing protons (or remnants)
  switch (fCuts.kinematics) {
    case 0: default: { Error("PPtoLL is intended for p-on-p collisions! Aborting!"); exit(0); break; }
    case 1: 
      _mx = GetParticle(Particle::IncomingBeam1)->M();
      _my = GetParticle(Particle::IncomingBeam2)->M();
      break;
    case 2:
      _mx = GetParticle(Particle::IncomingBeam1)->M();
      _my = fCuts.mxmin+(fCuts.mxmax-fCuts.mxmin)*x(8);
      break;
    case 3:
      _mx = fCuts.mxmin+(fCuts.mxmax-fCuts.mxmin)*x(8);
      _my = GetParticle(Particle::IncomingBeam2)->M();
      break;
    case 4:
      _mx = fCuts.mxmin+(fCuts.mxmax-fCuts.mxmin)*x(8);
      _my = fCuts.mxmin+(fCuts.mxmax-fCuts.mxmin)*x(9);
      break;
  }
  DebugInsideLoop(Form("outgoing remnants invariant mass: %f / %f (%.2f < M(X/Y) < %.2f)", _mx, _my, fCuts.mxmin, fCuts.mxmax));
  
  // Jacobian computation
  jac = 1.;
  jac *= (lqmax-lqmin)*_q1t; // d(q1t) . q1t
  jac *= (lqmax-lqmin)*_q2t; // d(q2t) . q2t
  jac *= 2.*pi; // d(phi1)
  jac *= 2.*pi; // d(phi2)
  jac *= (ymax-ymin); // d(y1)
  jac *= (ymax-ymin); // d(y2)
  switch (fCuts.kinematics) {
    case 1: default: break;
    case 2: jac *= (fCuts.mxmax-fCuts.mxmin)*2.*_my; break;
    case 3: jac *= (fCuts.mxmax-fCuts.mxmin)*2.*_mx; break;
    case 4: jac *= (fCuts.mxmax-fCuts.mxmin)*2.*_mx;
            jac *= (fCuts.mxmax-fCuts.mxmin)*2.*_my; break;
  } // d(mx/y**2)
  jac *= (fCuts.ptdiffmax-fCuts.ptdiffmin); // d(Dpt)
  jac *= 2.*pi; // d(phiDpt)
 
  const double integrand = INCqqbar();

  weight = jac*integrand;
  DebugInsideLoop(Form("Jacobian = %f\n\tIntegrand = %f\n\tdW = %f", jac, integrand, weight));
  
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
  double f1, f2;
  
  //const double stild = fS/2.*(1+sqrt(1.-(4*pow(mp2, 2))/pow(fS, 2)));
  
  // Inner photons
  const double q1tx = _q1t*cos(_phiq1t), q1ty = _q1t*sin(_phiq1t),
               q2tx = _q2t*cos(_phiq2t), q2ty = _q2t*sin(_phiq2t);
  DebugInsideLoop(Form("q1t(x/y) = %e / %e\n\t"
                       "q2t(x/y) = %e / %e", q1tx, q1ty, q2tx, q2ty));
  
  // Two-photon system
  const double ptsumx = q1tx+q2tx,
               ptsumy = q1ty+q2ty,
               ptsum = sqrt(pow(ptsumx, 2)+pow(ptsumy, 2));
  
  const double ptdiffx = _ptdiff*cos(_phiptdiff),
               ptdiffy = _ptdiff*sin(_phiptdiff);
  
  // Outgoing leptons
  const double pt1x = (ptsumx+ptdiffx)/2., pt1y = (ptsumy+ptdiffy)/2., pt1 = sqrt(pow(pt1x, 2)+pow(pt1y, 2)),
               pt2x = (ptsumx-ptdiffx)/2., pt2y = (ptsumy-ptdiffy)/2., pt2 = sqrt(pow(pt2x, 2)+pow(pt2y, 2));
  
  if (pt1<fCuts.ptmin or pt2<fCuts.ptmin) {
    return 0.;
  }
  /*if (pt1>fCuts.ptmax or pt2>fCuts.ptmax) {
    return 0.;
  }*/
  // transverse mass for the two leptons
  const double amt1 = sqrt(pow(pt1, 2)+ml2),
               amt2 = sqrt(pow(pt2, 2)+ml2);
  
  //=================================================================
  //     a window in transverse momentum difference
  //=================================================================
  if (idif==1 and fabs(pt1-pt2)>pdif) return 0.;

  //const double pcaptx = pt1x+pt2x, pcapty = pt1y+pt2y;
  // rapidity difference
  const double dely = fabs(_y1-_y2);
  //=================================================================
  //     a window in rapidity distance
  //=================================================================
  if (idely==1 and (dely<dely_min or dely>dely_max)) return 0.;
  
  //=================================================================
  //     auxiliary quantities
  //=================================================================

  const double alpha1 = amt1/fSqS*exp( _y1),
               alpha2 = amt2/fSqS*exp( _y2),
               beta1  = amt1/fSqS*exp(-_y1),
               beta2  = amt2/fSqS*exp(-_y2);
  DebugInsideLoop(Form("Sudakov parameters:\n\t"
                       "  alpha1/2 = %f / %f\n\t"
                       "   beta1/2 = %f / %f", alpha1, alpha2, beta1, beta2));

  const double q1t2 = pow(q1tx, 2)+pow(q1ty, 2),
               q2t2 = pow(q2tx, 2)+pow(q2ty, 2);

  //const double old_x2 = 0.; //FIXME figure out where this comes from
  //const double delta_x1 = (pow(_mx, 2)+q2t2)/((1.-old_x2)*fS);

  //x1 = alpha1+alpha2+delta_x1;
  const double x1 = alpha1+alpha2,
               x2 = beta1 +beta2;

  /*const double xi_x1 = log10(x1);
  const double xi_x2 = log10(x2);*/

  const double z1p = alpha1/x1, z1m = alpha2/x1,
               z2p = beta1 /x2, z2m = beta2 /x2;
  DebugInsideLoop(Form("z(1/2)p = %f / %f\n\t"
                       "z(1/2)m = %f / %f", z1p, z2p, z1m, z2m));
  
  if (x1>1. or x2>1.) return 0.; // sanity check

  // FIXME FIXME FIXME
  const double ak10 = GetParticle(Particle::IncomingBeam1)->E(),
               ak1z = GetParticle(Particle::IncomingBeam1)->GetMomentum().Pz(),
               ak20 = GetParticle(Particle::IncomingBeam2)->E(),
               ak2z = GetParticle(Particle::IncomingBeam2)->GetMomentum().Pz();
  DebugInsideLoop(Form("incoming particles: p1: %f / %f\n\t"
                       "                    p2: %f / %f", ak1z, ak10, ak2z, ak20));
  
  //=================================================================
  //     additional conditions for energy-momentum conservation
  //=================================================================
  
  const double s1_eff = x1*fS-pow(_q1t,2), s2_eff = x2*fS-pow(_q2t,2);
  const double invm = sqrt(pow(amt1,2)+pow(amt2,2)+2.*amt1*amt2*cosh(_y1-_y2)-pow(ptsum,2));
  DebugInsideLoop(Form("s(1/2)_eff = %f / %f GeV^2\n\tdilepton invariant mass = %f GeV", s1_eff, s2_eff, invm));

  switch (fCuts.kinematics) {
    case 2: if (sqrt(s1_eff)<=(_my+invm)) return 0.;
    case 3: if (sqrt(s2_eff)<=(_mx+invm)) return 0.;
    case 4: if (sqrt(s1_eff)<=(_my+invm)) return 0.;
            if (sqrt(s2_eff)<=(_mx+invm)) return 0.;
  }
  
  //const double qcaptx = pcaptx, qcapty = pcapty;
  
  //=================================================================
  //     four-momenta of the outgoing protons (or remnants)
  //=================================================================

  const double px_plus  = (1.-x1)*fabs(ak1z)*sqrt(2.),
               px_minus = (pow(_mx, 2)+pow(q1tx, 2)+pow(q1ty, 2))/2./px_plus;
  
  const double py_minus = (1.-x2)*fabs(ak2z)*sqrt(2.), // warning! sign of pz??
               py_plus  = (pow(_my, 2)+pow(q2tx, 2)+pow(q2ty, 2))/2./py_minus;

  DebugInsideLoop(Form("px_(+/-) = %f / %f\n\t"
                       "py_(+/-) = %f / %f", px_plus, px_minus, py_plus, py_minus));
  
  fPX = Particle::Momentum(-q1tx, -q1ty, (px_plus-px_minus)/sqrt(2.), (px_plus+px_minus)/sqrt(2.));
  fPY = Particle::Momentum(-q2tx, -q2ty, (py_plus-py_minus)/sqrt(2.), (py_plus+py_minus)/sqrt(2.));

  DebugInsideLoop(Form("First remnant:  (E,p) = (%f, %f, %f, %f)\n\t"
                       "Second remnant: (E,p) = (%f, %f, %f, %f)",
                       fPX.Px(), fPX.Py(), fPX.Pz(), fPX.E(),
                       fPY.Px(), fPY.Py(), fPY.Pz(), fPY.E()));
  
  assert(fabs(fPX.M()-_mx)<1.e-6);
  assert(fabs(fPY.M()-_my)<1.e-6);

  //=================================================================
  //     four-momenta of the outgoing l^+ and l^-
  //=================================================================

  Particle::Momentum p1(pt1x, pt1y, alpha1*ak1z+beta1*ak2z, alpha1*ak10+beta1*ak20),
                     p2(pt2x, pt2y, alpha2*ak1z+beta2*ak2z, alpha2*ak10+beta2*ak20);
  DebugInsideLoop(Form("unboosted first lepton:  (E,p), m = (%f, %f, %f, %f), %f\n\t"
                       "          second lepton: (E,p), m = (%f, %f, %f, %f), %f",
                       p1.Px(), p1.Py(), p1.Pz(), p1.E(), p1.M(),
                       p2.Px(), p2.Py(), p2.Pz(), p2.E(), p2.M()));

  fPl1 = Particle::Momentum(pt1x, pt1y, sqrt(pow(pt1, 2)+ml2)*sinh(_y1), sqrt(pow(pt1, 2)+ml2)*cosh(_y1));
  fPl2 = Particle::Momentum(pt2x, pt2y, sqrt(pow(pt2, 2)+ml2)*sinh(_y2), sqrt(pow(pt2, 2)+ml2)*cosh(_y2));

  DebugInsideLoop(Form("First lepton:  (E,p), m = (%f, %f, %f, %f), %f\n\t"
                       "Second lepton: (E,p), m = (%f, %f, %f, %f), %f",
                       fPl1.Px(), fPl1.Py(), fPl1.Pz(), fPl1.E(), fPl1.M(),
                       fPl2.Px(), fPl2.Py(), fPl2.Pz(), fPl2.E(), fPl2.M()));

  assert(fabs(fPl1.M()-GetParticle(Particle::CentralParticle1)->M())<1.e-6);
  assert(fabs(fPl2.M()-GetParticle(Particle::CentralParticle2)->M())<1.e-6);

  //=================================================================
  //     four-momenta squared of the virtual photons
  //=================================================================

  // FIXME FIXME FIXME /////////////////////
  Particle::Momentum q1(q1tx, q1ty, 0., 0.),
                     q2(q2tx, q2ty, 0., 0.);
  //////////////////////////////////////////
  
  DebugInsideLoop(Form("First photon*:  (E,p), m2 = (%f, %f, %f, %f), %e\n\t"
                       "Second photon*: (E,p), m2 = (%f, %f, %f, %f), %e",
                       q1.Px(), q1.Py(), q1.Pz(), q1.E(), q1.M2(),
                       q2.Px(), q2.Py(), q2.Pz(), q2.E(), q2.M2()));
  //const double q12 = q1.M2(), q22 = q2.M2();

  //=================================================================
  //     Mendelstam variables
  //=================================================================
  
  //const double shat = fS*x1*x2; // ishat = 1 (approximation)
  //const double shat = (q1+q2).M2(); // ishat = 2 (exact formula)

  const double that1 = (q1-p1).M2(), that2 = (q2-p2).M2(),
               uhat1 = (q1-p2).M2(), uhat2 = (q2-p1).M2();
  DebugInsideLoop(Form("that(1/2) = %f / %f\n\t"
                       "uhat(1/2) = %f / %f",
                       that1, that2, uhat1, uhat2));

  //const double mll = sqrt(shat);

  const double that = (that1+that2)/2., uhat = (uhat1+uhat2)/2.;

  //=================================================================
  //     matrix elements
  //=================================================================
  double amat2 = 0.;
  if (imethod==0) {
  
    //=================================================================
    //     on-shell formula for M^2
    //=================================================================
    const double ml4 = ml2*ml2, ml8 = ml4*ml4;
 
    const double term1  =  6. *ml8,
                 term2  = -3. *ml4*pow(that, 2),
                 term3  = -14.*ml4*that*uhat,
                 term4  = -3. *ml4*pow(uhat, 2),
                 term5  =      ml2*pow(that, 3),
                 term6  =  7.* ml2*pow(that, 2)*uhat,
                 term7  =  7.* ml2*that*pow(uhat, 2),
                 term8  =      ml2*pow(uhat, 3),
                 term9  =         -pow(that, 3)*uhat,
                 term10 =         -that*pow(uhat, 3);

    const double auxil_gamgam = -2.*(term1+term2+term3+term4+term5+term6+term7+term8+term9+term10)/(pow(ml2-that, 2)*pow(ml2-uhat, 2));
    const double g_em = sqrt(4.*pi*alpha_em);
    amat2 = pow(g_em, 4)*auxil_gamgam;
  }
  else if (imethod==1) {
  
    //=================================================================
    //     Wolfgang's formulae
    //=================================================================

    double aux2_1, aux2_2;

    const double ak1_x = z1m*pt1x-z1p*pt2x, ak1_y = z1m*pt1y-z1p*pt2y,
                 ak2_x = z2m*pt1x-z2p*pt2x, ak2_y = z2m*pt1y-z2p*pt2y;

    const double t1abs = (q1t2+x1*(pow(_mx, 2)-mp2)+pow(x1, 2)*mp2)/(1.-x1),
                 t2abs = (q2t2+x2*(pow(_my, 2)-mp2)+pow(x2, 2)*mp2)/(1.-x2);

    const double eps12 = ml2+z1p*z1m*t1abs,
                 eps22 = ml2+z2p*z2m*t2abs;

    const double Phi10 = 1./(pow(ak1_x+z1p*q2tx, 2)+pow(ak1_y+z1p*q2ty, 2)+eps12)
                        -1./(pow(ak1_x-z1m*q2tx, 2)+pow(ak1_y-z1m*q2ty, 2)+eps12),
                 Phi11_x = (ak1_x+z1p*q2tx)/(pow(ak1_x+z1p*q2tx, 2)+pow(ak1_y+z1p*q2ty, 2)+eps12)
                          -(ak1_x-z1m*q2tx)/(pow(ak1_x-z1m*q2tx, 2)+pow(ak1_y-z1m*q2ty, 2)+eps12),
                 Phi11_y = (ak1_y+z1p*q2ty)/(pow(ak1_x+z1p*q2tx, 2)+pow(ak1_y+z1p*q2ty, 2)+eps12)
                          -(ak1_y-z1m*q2ty)/(pow(ak1_x-z1m*q2tx, 2)+pow(ak1_y-z1m*q2ty, 2)+eps12),
                 Phi102 = pow(Phi10, 2);
    //const double Phi112 = pow(Phi11_x, 2)+pow(Phi11_y, 2);

    const double Phi20 = 1./(pow(ak2_x+z2p*q1tx, 2)+pow(ak2_y+z2p*q1ty, 2)+eps22)
                        -1./(pow(ak2_x-z2m*q1tx, 2)+pow(ak2_y-z2m*q1ty, 2)+eps22),
                 Phi21_x = (ak2_x+z2p*q1tx)/(pow(ak2_x+z2p*q1tx, 2)+pow(ak2_y+z2p*q1ty, 2)+eps22)
                          -(ak2_x-z2m*q1tx)/(pow(ak2_x-z2m*q1tx ,2)+pow(ak2_y-z2m*q1ty, 2)+eps22),
                 Phi21_y = (ak2_y+z2p*q1ty)/(pow(ak2_x+z2p*q1tx, 2)+pow(ak2_y+z2p*q1ty, 2)+eps22)
                          -(ak2_y-z2m*q1ty)/(pow(ak2_x-z2m*q1tx, 2)+pow(ak2_y-z2m*q1ty, 2)+eps22),
                 Phi202 = pow(Phi20, 2);
    //const double Phi212 = pow(Phi21_x, 2)+pow(Phi21_y, 2);

    /*aux2_1 = iterm11*(ml2+4.*z1p*z1m*t1abs)*Phi102
            +iterm22*(pow(z1p, 2)+pow(z1m, 2))*Phi112
            -iterm12*4.*z1p*z1m*(z1p-z1m)*Phi10*(q1tx*Phi11_x+q1ty*Phi11_y);
    aux2_2 = iterm11*(ml2+4.*z2p*z2m*t2abs)*Phi202
            +iterm22*(pow(z2p, 2)+pow(z2m, 2))*Phi212
            -iterm12*4.*z2p*z2m*(z2p-z2m)*Phi20*(q2tx*Phi21_x+q2ty*Phi21_y);*/

    const double Phi11_dot_e = (Phi11_x*q1tx+Phi11_y*q1ty)/_q1t, Phi11_cross_e = (Phi11_x*q1ty-Phi11_y*q1tx)/_q1t;
    const double Phi21_dot_e = (Phi21_x*q2tx+Phi21_y*q2ty)/_q2t, Phi21_cross_e = (Phi21_x*q2ty-Phi21_y*q2tx)/_q2t;
    DebugInsideLoop(Form("Phi1: E, px, py = %e, %e, %e\n\t"
                         "Phi2: E, px, py = %e, %e, %e\n\t"
                         "(dot):   %e / %e\n\t"
                         "(cross): %e / %e",
                         Phi10, Phi11_x, Phi11_y, Phi20, Phi21_x, Phi21_y,
                         Phi11_dot_e, Phi21_dot_e, Phi11_cross_e, Phi21_cross_e));

    aux2_1 = iterm11*(ml2+4.*pow(z1p*z1m, 2)*t1abs)*Phi102
            +iterm22*((pow(z1p, 2)+pow(z1m, 2))*(pow(Phi11_dot_e, 2)+pow(Phi11_cross_e, 2)))
            +itermtt*(pow(Phi11_cross_e, 2)-pow(Phi11_dot_e, 2))
            -iterm12*4.*z1p*z1m*(z1p-z1m)*Phi10*(q1tx*Phi11_x+q1ty*Phi11_y);

    aux2_2 = iterm11*(ml2+4.*pow(z2p*z2m, 2)*t2abs)*Phi202
            +iterm22*((pow(z2p, 2)+pow(z2m, 2))*(pow(Phi21_dot_e, 2)+pow(Phi21_cross_e, 2)))
            +itermtt*(pow(Phi21_cross_e, 2)-pow(Phi21_dot_e, 2))
            -iterm12*4.*z2p*z2m*(z2p-z2m)*Phi20*(q2tx*Phi21_x+q2ty*Phi21_y);


    //=================================================================
    //     convention of matrix element as in our kt-factorization
    //     for heavy flavours
    //=================================================================
    double amat2_1, amat2_2;
    
    amat2_1 = pow(4.*pi*alpha_em, 2)*pow(x1*x2*fS, 2)*aux2_1*2.*z1p*z1m*t1abs/(q1t2*q2t2)*t2abs/q2t2;
if (amat2_1<0.) {
//  std::cout << "amat21=" << amat2_1 << "\t" << x1 << "\t" << x2 << "\t" << z1p << "\t" << z1m << "\t" << t1abs << "\t" << t2abs << "\t" << q1t2 << "\t" << q2t2 << "\t" << aux2_1 << std::endl;
//return 0.;
}
    amat2_2 = pow(4.*pi*alpha_em, 2)*pow(x1*x2*fS, 2)*aux2_2*2.*z2p*z2m*t2abs/(q1t2*q2t2);

    //=================================================================
    //     symmetrization
    //=================================================================

    amat2 = (imat1*amat2_1+imat2*amat2_2)/2.;

    DebugInsideLoop(Form("aux2(1/2) = %e / %e\n\t"
                         "amat2(1/2), amat2 = %e / %e / %e", aux2_1, aux2_2, amat2_1, amat2_2, amat2));
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
  DebugInsideLoop(Form("Form factors: %e / %e", f1, f2));
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
  if (aintegral*_q1t*_q2t*_ptdiff!=0.) {
//std::cout << amat2 << "\t" << f1 << "\t" << f2 << "\t" << aintegral << "\t" << _ptdiff << std::endl;
    //GenericProcess::DumpPoint(Information);
    //Information(Form("matrix element: %E", aintegral*_q1t*_q2t*_ptdiff));
  }

  //=================================================================
  return aintegral*_q1t*_q2t*_ptdiff;
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
  
  if (!op1->SetMomentum(fPX)) { Error(Form("Invalid outgoing proton 1: energy: %.2f", fPX.E())); }
  if (!op2->SetMomentum(fPY)) { Error(Form("Invalid outgoing proton 2: energy: %.2f", fPY.E())); }

  // randomise the charge of the outgoing leptons
  int sign = (drand()>.5) ? +1 : -1;

  //=================================================================
  //     first outgoing lepton
  //=================================================================
  Particle* ol1 = GetParticle(Particle::CentralParticle1);
  ol1->SetPDGId(ol1->GetPDGId(), sign);
  if (!ol1->SetMomentum(fPl1)) { Error("Invalid outgoing lepton 1"); }

  //=================================================================
  //     second outgoing lepton
  //=================================================================
  Particle* ol2 = GetParticle(Particle::CentralParticle2);
  ol2->SetPDGId(ol2->GetPDGId(), -sign);
  if (!ol2->SetMomentum(fPl2)) { Error("Invalid outgoing lepton 2"); }
}
