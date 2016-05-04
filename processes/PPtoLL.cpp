#include "PPtoLL.h"
#include <assert.h>

PPtoLL::PPtoLL() : GenericKTProcess("gamma,gamma->l+,l-", 4, Particle::Photon, Particle::Muon)
{}

void
PPtoLL::PrepareKTKinematics()
{
  ////////////////////////////////////
  fYmin = fCuts.etamin;             //
  //fYmin = EtaToY(fCuts.etamin, GetParticle(Particle::CentralParticle1)->M(), pt);
  fYmax = fCuts.etamax;             //
  //fYmax = EtaToY(fCuts.etamax);
  ///////////// FIXME ////////////////
  
  // Outgoing leptons  
  fY1 = fYmin+(fYmax-fYmin)*x(4);
  fY2 = fYmin+(fYmax-fYmin)*x(5);
  DebugInsideLoop(Form("leptons rapidities (%.2f < y < %.2f): %f / %f", fYmin, fYmax, fY1, fY2));
 
  if (fCuts.ptdiffmax<0.) fCuts.ptdiffmax = 400.; //FIXME
  fPtDiff = fCuts.ptdiffmin+(fCuts.ptdiffmax-fCuts.ptdiffmin)*x(6);
  fPhiPtDiff = 2.*Constants::Pi*x(7);
  DebugInsideLoop(Form("leptons pt difference:\n\t"
                       "  mag = %f (%.2f < Dpt < %.2f)\n\t"
                       "  phi = %f",
                       fPtDiff, fCuts.ptdiffmin, fCuts.ptdiffmax, fPhiPtDiff));
}

double
PPtoLL::ComputeJacobian()
{
  double jac = GenericKTProcess::MinimalJacobian();
  jac *= (fYmax-fYmin); // d(y1)
  jac *= (fYmax-fYmin); // d(y2)
  jac *= (fCuts.ptdiffmax-fCuts.ptdiffmin); // d(Dpt)
  jac *= 2.*Constants::Pi; // d(phiDpt)
  
  return jac;
}

double
PPtoLL::ComputeKTFactorisedMatrixElement()
{
  const double mp = Particle::GetMassFromPDGId(Particle::Proton), mp2 = pow(mp, 2);
  const double ml = GetParticle(Particle::CentralParticle1)->M(), ml2 = pow(ml, 2);

  const unsigned int iterm11 = 1, // Long-long
                     iterm22 = 1, // Trans-trans
                     iterm12 = 1, // Long-trans
                     itermtt = 1; // Trans-trans(')

  //=================================================================
  //     How matrix element is calculated
  //=================================================================
  const bool off_shell = true;
  
  //=================================================================
  //     two terms in Wolfgang's formula for 
  //     off-shell gamma gamma --> l^+ l^-
  //=================================================================
  const unsigned int imat1 = 2,
                     imat2 = 0;

  //=================================================================
  //     extra cuts on the p1t(l) and p2t(l) plane
  //=================================================================
  const bool pt_window = false;
  const double pdif = 2.5;
  
  //=================================================================
  //     the distance in rapidity between l^+ and l^-
  //=================================================================
  const bool delta_y_window = false;
  const double dely_min = 4.0,
               dely_max = 5.0;
    
  //=================================================================
  //     matrix element computation
  //=================================================================
  //const double stild = fS/2.*(1+sqrt(1.-(4*pow(mp2, 2))/pow(fS, 2)));
  
  // Inner photons
  const double q1tx = fQT1*cos(fPhiQT1), q1ty = fQT1*sin(fPhiQT1),
               q2tx = fQT2*cos(fPhiQT2), q2ty = fQT2*sin(fPhiQT2);
  DebugInsideLoop(Form("q1t(x/y) = %e / %e\n\t"
                       "q2t(x/y) = %e / %e", q1tx, q1ty, q2tx, q2ty));
  
  // Two-photon system
  const double ptsumx = q1tx+q2tx,
               ptsumy = q1ty+q2ty,
               ptsum = sqrt(pow(ptsumx, 2)+pow(ptsumy, 2));
  
  const double ptdiffx = fPtDiff*cos(fPhiPtDiff),
               ptdiffy = fPtDiff*sin(fPhiPtDiff);
  
  // Outgoing leptons
  const double pt1x = (ptsumx+ptdiffx)/2., pt1y = (ptsumy+ptdiffy)/2., pt1 = sqrt(pow(pt1x, 2)+pow(pt1y, 2)),
               pt2x = (ptsumx-ptdiffx)/2., pt2y = (ptsumy-ptdiffy)/2., pt2 = sqrt(pow(pt2x, 2)+pow(pt2y, 2));
  
  if (pt1<fCuts.ptmin or pt2<fCuts.ptmin) return 0.;
  //if (pt1>fCuts.ptmax or pt2>fCuts.ptmax) return 0.;
  // transverse mass for the two leptons
  const double amt1 = sqrt(pow(pt1, 2)+ml2),
               amt2 = sqrt(pow(pt2, 2)+ml2);
  
  //=================================================================
  //     a window in transverse momentum difference
  //=================================================================
  if (pt_window and fabs(pt1-pt2)>pdif) return 0.;

  //const double pcaptx = pt1x+pt2x, pcapty = pt1y+pt2y;
  // rapidity difference
  const double dely = fabs(fY1-fY2);
  //=================================================================
  //     a window in rapidity distance
  //=================================================================
  if (delta_y_window and (dely<dely_min or dely>dely_max)) return 0.;
  
  //=================================================================
  //     auxiliary quantities
  //=================================================================

  const double alpha1 = amt1/fSqS*exp( fY1),
               alpha2 = amt2/fSqS*exp( fY2),
               beta1  = amt1/fSqS*exp(-fY1),
               beta2  = amt2/fSqS*exp(-fY2);
  DebugInsideLoop(Form("Sudakov parameters:\n\t"
                       "  alpha1/2 = %f / %f\n\t"
                       "   beta1/2 = %f / %f", alpha1, alpha2, beta1, beta2));

  const double q1t2 = pow(q1tx, 2)+pow(q1ty, 2),
               q2t2 = pow(q2tx, 2)+pow(q2ty, 2);

  //const double old_x2 = 0.; //FIXME figure out where this comes from
  //const double delta_x1 = (pow(fMX, 2)+q2t2)/((1.-old_x2)*fS);

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
  
  const double s1_eff = x1*fS-pow(fQT1,2), s2_eff = x2*fS-pow(fQT2,2);
  const double invm = sqrt(pow(amt1,2)+pow(amt2,2)+2.*amt1*amt2*cosh(fY1-fY2)-pow(ptsum,2));
  DebugInsideLoop(Form("s(1/2)_eff = %f / %f GeV^2\n\t"
                       "dilepton invariant mass = %f GeV", s1_eff, s2_eff, invm));

  switch (fCuts.kinematics) {
    case Kinematics::ElasticInelastic:   if (sqrt(s1_eff)<=(fMY+invm)) return 0.;
    case Kinematics::InelasticElastic:   if (sqrt(s2_eff)<=(fMX+invm)) return 0.;
    case Kinematics::InelasticInelastic: if (sqrt(s1_eff)<=(fMY+invm)) return 0.;
                                         if (sqrt(s2_eff)<=(fMX+invm)) return 0.;
    default: break;
  }
  
  //const double qcaptx = pcaptx, qcapty = pcapty;
  
  //=================================================================
  //     four-momenta of the outgoing protons (or remnants)
  //=================================================================

  const double px_plus  = (1.-x1)*fabs(ak1z)*sqrt(2.),
               px_minus = (pow(fMX, 2)+pow(q1tx, 2)+pow(q1ty, 2))/2./px_plus;
  
  const double py_minus = (1.-x2)*fabs(ak2z)*sqrt(2.), // warning! sign of pz??
               py_plus  = (pow(fMY, 2)+pow(q2tx, 2)+pow(q2ty, 2))/2./py_minus;

  DebugInsideLoop(Form("px_(+/-) = %f / %f\n\t"
                       "py_(+/-) = %f / %f", px_plus, px_minus, py_plus, py_minus));
  
  fPX = Particle::Momentum(-q1tx, -q1ty, (px_plus-px_minus)/sqrt(2.), (px_plus+px_minus)/sqrt(2.));
  fPY = Particle::Momentum(-q2tx, -q2ty, (py_plus-py_minus)/sqrt(2.), (py_plus+py_minus)/sqrt(2.));

  DebugInsideLoop(Form("First remnant:  (E,p) = (%f, %f, %f, %f)\n\t"
                       "Second remnant: (E,p) = (%f, %f, %f, %f)",
                       fPX.Px(), fPX.Py(), fPX.Pz(), fPX.E(),
                       fPY.Px(), fPY.Py(), fPY.Pz(), fPY.E()));
  
  assert(fabs(fPX.M()-fMX)<1.e-6);
  assert(fabs(fPY.M()-fMY)<1.e-6);

  //=================================================================
  //     four-momenta of the outgoing l^+ and l^-
  //=================================================================

  Particle::Momentum p1(pt1x, pt1y, alpha1*ak1z+beta1*ak2z, alpha1*ak10+beta1*ak20),
                     p2(pt2x, pt2y, alpha2*ak1z+beta2*ak2z, alpha2*ak10+beta2*ak20);
  DebugInsideLoop(Form("unboosted first lepton:  (E,p), m = (%f, %f, %f, %f), %f\n\t"
                       "          second lepton: (E,p), m = (%f, %f, %f, %f), %f",
                       p1.Px(), p1.Py(), p1.Pz(), p1.E(), p1.M(),
                       p2.Px(), p2.Py(), p2.Pz(), p2.E(), p2.M()));

  fPl1 = Particle::Momentum(pt1x, pt1y, sqrt(pow(pt1, 2)+ml2)*sinh(fY1), sqrt(pow(pt1, 2)+ml2)*cosh(fY1));
  fPl2 = Particle::Momentum(pt2x, pt2y, sqrt(pow(pt2, 2)+ml2)*sinh(fY2), sqrt(pow(pt2, 2)+ml2)*cosh(fY2));

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
  if (!off_shell) {
  
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
    const double g_em = sqrt(4.*Constants::Pi*Constants::AlphaEM);
    amat2 = pow(g_em, 4)*auxil_gamgam;
  }
  else if (off_shell) {
  
    //=================================================================
    //     Wolfgang's formulae
    //=================================================================

    double aux2_1, aux2_2;

    const double ak1_x = z1m*pt1x-z1p*pt2x, ak1_y = z1m*pt1y-z1p*pt2y,
                 ak2_x = z2m*pt1x-z2p*pt2x, ak2_y = z2m*pt1y-z2p*pt2y;

    const double t1abs = (q1t2+x1*(pow(fMX, 2)-mp2)+pow(x1, 2)*mp2)/(1.-x1),
                 t2abs = (q2t2+x2*(pow(fMY, 2)-mp2)+pow(x2, 2)*mp2)/(1.-x2);

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

    const double Phi11_dot_e = (Phi11_x*q1tx+Phi11_y*q1ty)/fQT1, Phi11_cross_e = (Phi11_x*q1ty-Phi11_y*q1tx)/fQT1;
    const double Phi21_dot_e = (Phi21_x*q2tx+Phi21_y*q2ty)/fQT2, Phi21_cross_e = (Phi21_x*q2ty-Phi21_y*q2tx)/fQT2;
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
    
    amat2_1 = pow(4.*Constants::Pi*Constants::AlphaEM, 2)*pow(x1*x2*fS, 2)*aux2_1*2.*z1p*z1m*t1abs/(q1t2*q2t2)*t2abs/q2t2;
    amat2_2 = pow(4.*Constants::Pi*Constants::AlphaEM, 2)*pow(x1*x2*fS, 2)*aux2_2*2.*z2p*z2m*t2abs/(q1t2*q2t2);

    //=================================================================
    //     symmetrization
    //=================================================================

    amat2 = (imat1*amat2_1+imat2*amat2_2)/2.;

    DebugInsideLoop(Form("aux2(1/2) = %e / %e\n\t"
                         "amat2(1/2), amat2 = %e / %e / %e", aux2_1, aux2_2, amat2_1, amat2_2, amat2));
    /*const double xx1 = alpha1+alpha2, xx2 = beta1+beta2;

    const double sudakov_2 = (pow(fMX, 2)-mp2+q2t2+xx2*mp2)/((1.-xx2)*fS);
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
  
  double f1, f2;
  switch (fCuts.kinematics) {
    case Kinematics::ElasticElastic:
      f1 = GenericKTProcess::ElasticFlux(x1, q1t2);
      f2 = GenericKTProcess::ElasticFlux(x2, q2t2);
      break;
    case Kinematics::ElasticInelastic:
      f1 = GenericKTProcess::ElasticFlux(x1, q1t2);
      f2 = GenericKTProcess::InelasticFlux(x2, q2t2, fMY);
      break;
    case Kinematics::InelasticElastic:
      f1 = GenericKTProcess::InelasticFlux(x1, q1t2, fMX);
      f2 = GenericKTProcess::ElasticFlux(x2, q2t2);
      break;
    case Kinematics::InelasticInelastic:
      f1 = GenericKTProcess::InelasticFlux(x1, q1t2, fMX);
      f2 = GenericKTProcess::InelasticFlux(x2, q2t2, fMY);
      break;
    default: return 0.;
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

  const double aintegral = (2.*Constants::Pi)
                         *1./(16.*pow(Constants::Pi, 2)
                         *pow(x1*x2*fS, 2)) * amat2
                         * f1/Constants::Pi
                         * f2/Constants::Pi
                         *(1./4.)*Constants::GeV2toBarn
                         * 0.5*4./(4.*Constants::Pi);
  if (aintegral*fQT1*fQT2*fPtDiff!=0.) {
    //GenericProcess::DumpPoint(Information);
    //Information(Form("matrix element: %E", aintegral*fQT1*fQT2*fPtDiff));
  }

  //=================================================================
  return aintegral*fQT1*fQT2*fPtDiff;
  //=================================================================
}

void
PPtoLL::FillCentralParticlesKinematics()
{
  // randomise the charge of the outgoing leptons
  int sign = (drand()>.5) ? +1 : -1;

  //=================================================================
  //     first outgoing lepton
  //=================================================================
  Particle* ol1 = GetParticle(Particle::CentralParticle1);
  ol1->SetPDGId(ol1->GetPDGId(), sign);
  ol1->status = Particle::FinalState;
  if (!ol1->SetMomentum(fPl1)) { Error("Invalid outgoing lepton 1"); }

  //=================================================================
  //     second outgoing lepton
  //=================================================================
  Particle* ol2 = GetParticle(Particle::CentralParticle2);
  ol2->SetPDGId(ol2->GetPDGId(), -sign);
  ol2->status = Particle::FinalState;
  if (!ol2->SetMomentum(fPl2)) { Error("Invalid outgoing lepton 2"); }
}
