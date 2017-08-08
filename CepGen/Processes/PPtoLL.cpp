#include "PPtoLL.h"
#include <assert.h>

using namespace CepGen::Process;

PPtoLL::PPtoLL() : GenericKTProcess("gamma,gamma->l+,l-", 4, Particle::Photon, Particle::Muon)
{}

void
PPtoLL::prepareKTKinematics()
{
  ////////////////////////////////////
  y_min_ = cuts_.etamin;             //
  //y_min_ = EtaToY(cuts_.etamin, particlePtr(Particle::CentralParticle1)->mass(), pt);
  y_max_ = cuts_.etamax;             //
  //y_max_ = EtaToY(cuts_.etamax);
  ///////////// FIXME ////////////////
  
  // Outgoing leptons  
  y1_ = y_min_+(y_max_-y_min_)*x(4);
  y2_ = y_min_+(y_max_-y_min_)*x(5);
  DebuggingInsideLoop(Form("leptons rapidities (%.2f < y < %.2f): %f / %f", y_min_, y_max_, y1_, y2_));
 
  if (cuts_.ptdiffmax<0.) cuts_.ptdiffmax = 400.; //FIXME
  pt_diff_ = cuts_.ptdiffmin+(cuts_.ptdiffmax-cuts_.ptdiffmin)*x(6);
  phi_pt_diff_ = 2.*M_PI*x(7);
  DebuggingInsideLoop(Form("leptons pt difference:\n\t"
                           "  mag = %f (%.2f < Dpt < %.2f)\n\t"
                           "  phi = %f",
                           pt_diff_, cuts_.ptdiffmin, cuts_.ptdiffmax, phi_pt_diff_));
}

double
PPtoLL::computeJacobian()
{
  double jac = GenericKTProcess::minimalJacobian();
  jac *= ( y_max_-y_min_ ); // d(y1)
  jac *= ( y_max_-y_min_ ); // d(y2)
  jac *= ( cuts_.ptdiffmax-cuts_.ptdiffmin ); // d(Dpt)
  jac *= 2.*M_PI; // d(phiDpt)
  
  return jac;
}

double
PPtoLL::computeKTFactorisedMatrixElement()
{
  const double mp = Particle::massFromPDGId(Particle::Proton), mp2 = mp*mp;
  const double ml = particlePtr( Particle::CentralParticle1 )->mass(), ml2 = ml*ml;

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
  //const double stild = s_/2.*(1+sqrt(1.-(4*pow(mp2, 2))/s_*s_));
  
  // Inner photons
  const double q1tx = qt1_*cos(phi_qt1_), q1ty = qt1_*sin(phi_qt1_),
               q2tx = qt2_*cos(phi_qt2_), q2ty = qt2_*sin(phi_qt2_);
  DebuggingInsideLoop(Form("q1t(x/y) = %e / %e\n\t"
                           "q2t(x/y) = %e / %e", q1tx, q1ty, q2tx, q2ty));
  
  // Two-photon system
  const double ptsumx = q1tx+q2tx,
               ptsumy = q1ty+q2ty,
               ptsum = sqrt(ptsumx*ptsumx+ptsumy*ptsumy);
  
  const double ptdiffx = pt_diff_*cos(phi_pt_diff_),
               ptdiffy = pt_diff_*sin(phi_pt_diff_);
  
  // Outgoing leptons
  const double pt1x = (ptsumx+ptdiffx)/2., pt1y = (ptsumy+ptdiffy)/2., pt1 = sqrt(pt1x*pt1x+pt1y*pt1y),
               pt2x = (ptsumx-ptdiffx)/2., pt2y = (ptsumy-ptdiffy)/2., pt2 = sqrt(pt2x*pt2x+pt2y*pt2y);
  
  if (pt1<cuts_.ptmin or pt2<cuts_.ptmin) return 0.;
  //if (pt1>cuts_.ptmax or pt2>cuts_.ptmax) return 0.;
  // transverse mass for the two leptons
  const double amt1 = sqrt(pt1*pt1+ml2),
               amt2 = sqrt(pt2*pt2+ml2);
  
  //=================================================================
  //     a window in transverse momentum difference
  //=================================================================
  if (pt_window and fabs(pt1-pt2)>pdif) return 0.;

  //const double pcaptx = pt1x+pt2x, pcapty = pt1y+pt2y;
  // rapidity difference
  const double dely = fabs(y1_-y2_);
  //=================================================================
  //     a window in rapidity distance
  //=================================================================
  if (delta_y_window and (dely<dely_min or dely>dely_max)) return 0.;
  
  //=================================================================
  //     auxiliary quantities
  //=================================================================

  const double alpha1 = amt1/sqs_*exp( y1_),
               alpha2 = amt2/sqs_*exp( y2_),
               beta1  = amt1/sqs_*exp(-y1_),
               beta2  = amt2/sqs_*exp(-y2_);
  DebuggingInsideLoop(Form("Sudakov parameters:\n\t"
                           "  alpha1/2 = %f / %f\n\t"
                           "   beta1/2 = %f / %f", alpha1, alpha2, beta1, beta2));

  const double q1t2 = q1tx*q1tx+q1ty*q1ty,
               q2t2 = q2tx*q2tx+q2ty*q2ty;

  //const double old_x2 = 0.; //FIXME figure out where this comes from
  //const double delta_x1 = (MX_*MX_+q2t2)/((1.-old_x2)*s_);

  //x1 = alpha1+alpha2+delta_x1;
  const double x1 = alpha1+alpha2,
               x2 = beta1 +beta2;

  /*const double xi_x1 = log10(x1);
  const double xi_x2 = log10(x2);*/

  const double z1p = alpha1/x1, z1m = alpha2/x1,
               z2p = beta1 /x2, z2m = beta2 /x2;
  DebuggingInsideLoop(Form("z(1/2)p = %f / %f\n\t"
                           "z(1/2)m = %f / %f", z1p, z2p, z1m, z2m));
  
  if (x1>1. or x2>1.) return 0.; // sanity check

  // FIXME FIXME FIXME
  const double ak10 = particlePtr(Particle::IncomingBeam1)->energy(),
               ak1z = particlePtr(Particle::IncomingBeam1)->momentum().pz(),
               ak20 = particlePtr(Particle::IncomingBeam2)->energy(),
               ak2z = particlePtr(Particle::IncomingBeam2)->momentum().pz();
  DebuggingInsideLoop(Form("incoming particles: p1: %f / %f\n\t"
                           "                    p2: %f / %f", ak1z, ak10, ak2z, ak20));
  
  //=================================================================
  //     additional conditions for energy-momentum conservation
  //=================================================================
  
  const double s1_eff = x1*s_-qt1_*qt1_, s2_eff = x2*s_-qt2_*qt2_;
  const double invm = sqrt(amt1*amt1+amt2*amt2+2.*amt1*amt2*cosh(y1_-y2_)-ptsum*ptsum);
  DebuggingInsideLoop(Form("s(1/2)_eff = %f / %f GeV^2\n\t"
                           "dilepton invariant mass = %f GeV", s1_eff, s2_eff, invm));

  switch (cuts_.kinematics) {
    case Kinematics::ElasticInelastic:   if (sqrt(s1_eff)<=(MY_+invm)) return 0.;
    case Kinematics::InelasticElastic:   if (sqrt(s2_eff)<=(MX_+invm)) return 0.;
    case Kinematics::InelasticInelastic: if (sqrt(s1_eff)<=(MY_+invm)) return 0.;
                                         if (sqrt(s2_eff)<=(MX_+invm)) return 0.;
    default: break;
  }
  
  //const double qcaptx = pcaptx, qcapty = pcapty;
  
  //=================================================================
  //     four-momenta of the outgoing protons (or remnants)
  //=================================================================

  const double px_plus  = (1.-x1)*fabs(ak1z)*sqrt(2.),
               px_minus = (MX_*MX_+q1tx*q1tx+q1ty*q1ty)/2./px_plus;
  
  const double py_minus = (1.-x2)*fabs(ak2z)*sqrt(2.), // warning! sign of pz??
               py_plus  = (MY_*MY_+q2tx*q2tx+q2ty*q2ty)/2./py_minus;

  DebuggingInsideLoop(Form("px_(+/-) = %f / %f\n\t"
                           "py_(+/-) = %f / %f", px_plus, px_minus, py_plus, py_minus));
  
  PX_ = Particle::Momentum(-q1tx, -q1ty, (px_plus-px_minus)/sqrt(2.), (px_plus+px_minus)/sqrt(2.));
  PY_ = Particle::Momentum(-q2tx, -q2ty, (py_plus-py_minus)/sqrt(2.), (py_plus+py_minus)/sqrt(2.));

  DebuggingInsideLoop(Form("First remnant:  (E,p) = (%f, %f, %f, %f)\n\t"
                           "Second remnant: (E,p) = (%f, %f, %f, %f)",
                           PX_.px(), PX_.py(), PX_.pz(), PX_.energy(),
                           PY_.px(), PY_.py(), PY_.pz(), PY_.energy()));
  
  assert(fabs(PX_.mass()-MX_)<1.e-6);
  assert(fabs(PY_.mass()-MY_)<1.e-6);

  //=================================================================
  //     four-momenta of the outgoing l^+ and l^-
  //=================================================================

  Particle::Momentum p1(pt1x, pt1y, alpha1*ak1z+beta1*ak2z, alpha1*ak10+beta1*ak20),
                     p2(pt2x, pt2y, alpha2*ak1z+beta2*ak2z, alpha2*ak10+beta2*ak20);
  DebuggingInsideLoop(Form("unboosted first lepton:  (E,p), m = (%f, %f, %f, %f), %f\n\t"
                           "          second lepton: (E,p), m = (%f, %f, %f, %f), %f",
                           p1.px(), p1.py(), p1.pz(), p1.energy(), p1.mass(),
                           p2.px(), p2.py(), p2.pz(), p2.energy(), p2.mass()));

  Pl1_ = Particle::Momentum(pt1x, pt1y, sqrt(pt1*pt1+ml2)*sinh(y1_), sqrt(pt1*pt1+ml2)*cosh(y1_));
  Pl2_ = Particle::Momentum(pt2x, pt2y, sqrt(pt2*pt2+ml2)*sinh(y2_), sqrt(pt2*pt2+ml2)*cosh(y2_));

  DebuggingInsideLoop(Form("First lepton:  (E,p), m = (%f, %f, %f, %f), %f\n\t"
                           "Second lepton: (E,p), m = (%f, %f, %f, %f), %f",
                           Pl1_.px(), Pl1_.py(), Pl1_.pz(), Pl1_.energy(), Pl1_.mass(),
                           Pl2_.px(), Pl2_.py(), Pl2_.pz(), Pl2_.energy(), Pl2_.mass()));

  assert(fabs(Pl1_.mass()-particlePtr(Particle::CentralParticle1)->mass())<1.e-6);
  assert(fabs(Pl2_.mass()-particlePtr(Particle::CentralParticle2)->mass())<1.e-6);

  //=================================================================
  //     four-momenta squared of the virtual photons
  //=================================================================

  // FIXME FIXME FIXME /////////////////////
  Particle::Momentum q1(q1tx, q1ty, 0., 0.),
                     q2(q2tx, q2ty, 0., 0.);
  //////////////////////////////////////////
  
  DebuggingInsideLoop(Form("First photon*:  (E,p), m2 = (%f, %f, %f, %f), %e\n\t"
                           "Second photon*: (E,p), m2 = (%f, %f, %f, %f), %e",
                           q1.px(), q1.py(), q1.pz(), q1.energy(), q1.mass2(),
                           q2.px(), q2.py(), q2.pz(), q2.energy(), q2.mass2()));
  //const double q12 = q1.mass2(), q22 = q2.mass2();

  //=================================================================
  //     Mendelstam variables
  //=================================================================
  
  //const double shat = s_*x1*x2; // ishat = 1 (approximation)
  //const double shat = (q1+q2).mass2(); // ishat = 2 (exact formula)

  const double that1 = (q1-p1).mass2(), that2 = (q2-p2).mass2(),
               uhat1 = (q1-p2).mass2(), uhat2 = (q2-p1).mass2();
  DebuggingInsideLoop(Form("that(1/2) = %f / %f\n\t"
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
    const double g_em = sqrt(4.*M_PI*Constants::alphaEM);
    amat2 = pow(g_em, 4)*auxil_gamgam;
  }
  else if (off_shell) {
  
    //=================================================================
    //     Wolfgang's formulae
    //=================================================================

    double aux2_1, aux2_2;

    const double ak1_x = z1m*pt1x-z1p*pt2x, ak1_y = z1m*pt1y-z1p*pt2y,
                 ak2_x = z2m*pt1x-z2p*pt2x, ak2_y = z2m*pt1y-z2p*pt2y;

    const double t1abs = (q1t2+x1*(MX_*MX_-mp2)+x1*x1*mp2)/(1.-x1),
                 t2abs = (q2t2+x2*(MY_*MY_-mp2)+x2*x2*mp2)/(1.-x2);

    const double eps12 = ml2+z1p*z1m*t1abs,
                 eps22 = ml2+z2p*z2m*t2abs;

    const double Phi10 = 1./(pow(ak1_x+z1p*q2tx, 2)+pow(ak1_y+z1p*q2ty, 2)+eps12)
                        -1./(pow(ak1_x-z1m*q2tx, 2)+pow(ak1_y-z1m*q2ty, 2)+eps12),
                 Phi11_x = (ak1_x+z1p*q2tx)/(pow(ak1_x+z1p*q2tx, 2)+pow(ak1_y+z1p*q2ty, 2)+eps12)
                          -(ak1_x-z1m*q2tx)/(pow(ak1_x-z1m*q2tx, 2)+pow(ak1_y-z1m*q2ty, 2)+eps12),
                 Phi11_y = (ak1_y+z1p*q2ty)/(pow(ak1_x+z1p*q2tx, 2)+pow(ak1_y+z1p*q2ty, 2)+eps12)
                          -(ak1_y-z1m*q2ty)/(pow(ak1_x-z1m*q2tx, 2)+pow(ak1_y-z1m*q2ty, 2)+eps12),
                 Phi102 = Phi10*Phi10;
    //const double Phi112 = pow(Phi11_x, 2)+pow(Phi11_y, 2);

    const double Phi20 = 1./(pow(ak2_x+z2p*q1tx, 2)+pow(ak2_y+z2p*q1ty, 2)+eps22)
                        -1./(pow(ak2_x-z2m*q1tx, 2)+pow(ak2_y-z2m*q1ty, 2)+eps22),
                 Phi21_x = (ak2_x+z2p*q1tx)/(pow(ak2_x+z2p*q1tx, 2)+pow(ak2_y+z2p*q1ty, 2)+eps22)
                          -(ak2_x-z2m*q1tx)/(pow(ak2_x-z2m*q1tx ,2)+pow(ak2_y-z2m*q1ty, 2)+eps22),
                 Phi21_y = (ak2_y+z2p*q1ty)/(pow(ak2_x+z2p*q1tx, 2)+pow(ak2_y+z2p*q1ty, 2)+eps22)
                          -(ak2_y-z2m*q1ty)/(pow(ak2_x-z2m*q1tx, 2)+pow(ak2_y-z2m*q1ty, 2)+eps22),
                 Phi202 = Phi20*Phi20;
    //const double Phi212 = pow(Phi21_x, 2)+pow(Phi21_y, 2);

    /*aux2_1 = iterm11*(ml2+4.*z1p*z1m*t1abs)*Phi102
            +iterm22*(pow(z1p, 2)+pow(z1m, 2))*Phi112
            -iterm12*4.*z1p*z1m*(z1p-z1m)*Phi10*(q1tx*Phi11_x+q1ty*Phi11_y);
    aux2_2 = iterm11*(ml2+4.*z2p*z2m*t2abs)*Phi202
            +iterm22*(pow(z2p, 2)+pow(z2m, 2))*Phi212
            -iterm12*4.*z2p*z2m*(z2p-z2m)*Phi20*(q2tx*Phi21_x+q2ty*Phi21_y);*/

    const double Phi11_dot_e = (Phi11_x*q1tx+Phi11_y*q1ty)/qt1_, Phi11_cross_e = (Phi11_x*q1ty-Phi11_y*q1tx)/qt1_;
    const double Phi21_dot_e = (Phi21_x*q2tx+Phi21_y*q2ty)/qt2_, Phi21_cross_e = (Phi21_x*q2ty-Phi21_y*q2tx)/qt2_;
    DebuggingInsideLoop(Form("Phi1: E, px, py = %e, %e, %e\n\t"
                             "Phi2: E, px, py = %e, %e, %e\n\t"
                             "(dot):   %e / %e\n\t"
                             "(cross): %e / %e",
                             Phi10, Phi11_x, Phi11_y, Phi20, Phi21_x, Phi21_y,
                             Phi11_dot_e, Phi21_dot_e, Phi11_cross_e, Phi21_cross_e));

    aux2_1 = iterm11*(ml2+4.*pow(z1p*z1m, 2)*t1abs)*Phi102
            +iterm22*((z1p*z1p+z1m*z1m)*(pow(Phi11_dot_e, 2)+pow(Phi11_cross_e, 2)))
            +itermtt*(pow(Phi11_cross_e, 2)-pow(Phi11_dot_e, 2))
            -iterm12*4.*z1p*z1m*(z1p-z1m)*Phi10*(q1tx*Phi11_x+q1ty*Phi11_y);

    aux2_2 = iterm11*(ml2+4.*pow(z2p*z2m, 2)*t2abs)*Phi202
            +iterm22*((z2p*z2p+z2m*z2m)*(pow(Phi21_dot_e, 2)+pow(Phi21_cross_e, 2)))
            +itermtt*(pow(Phi21_cross_e, 2)-pow(Phi21_dot_e, 2))
            -iterm12*4.*z2p*z2m*(z2p-z2m)*Phi20*(q2tx*Phi21_x+q2ty*Phi21_y);


    //=================================================================
    //     convention of matrix element as in our kt-factorization
    //     for heavy flavours
    //=================================================================
    double amat2_1, amat2_2;
    
    amat2_1 = pow(4.*M_PI*Constants::alphaEM, 2)*pow(x1*x2*s_, 2)*aux2_1*2.*z1p*z1m*t1abs/(q1t2*q2t2)*t2abs/q2t2;
    amat2_2 = pow(4.*M_PI*Constants::alphaEM, 2)*pow(x1*x2*s_, 2)*aux2_2*2.*z2p*z2m*t2abs/(q1t2*q2t2);

    //=================================================================
    //     symmetrization
    //=================================================================

    amat2 = (imat1*amat2_1+imat2*amat2_2)/2.;

    DebuggingInsideLoop(Form("aux2(1/2) = %e / %e\n\t"
                             "amat2(1/2), amat2 = %e / %e / %e", aux2_1, aux2_2, amat2_1, amat2_2, amat2));
    /*const double xx1 = alpha1+alpha2, xx2 = beta1+beta2;

    const double sudakov_2 = (MX_*MX_-mp2+q2t2+xx2*mp2)/((1.-xx2)*s_);
    const double sudakov_1 = (q1t2 + xx1*mp2)/((1.-xx1)*s_);
    const double ratio1 = sudakov_1 / xx1,
                 ratio2 = sudakov_2 / xx2;*/

    //if (ratio1>0.01) return 0.;
  }

  //============================================
  //     unintegrated photon distributions
  //     interpolation on double logarithmic grid
  //     of inelastic distributions
  //============================================
  
  GenericKTProcess::computeIncomingFluxes( x1, q1t2, x2, q2t2 );

  //=================================================================
  //     factor 2.*pi below from integration over phi_sum
  //     factor 1/4 below from jacobian of transformations
  //     factors 1/pi and 1/pi due to integration
  //     over d^2 kappa_1 d^2 kappa_2 instead d kappa_1^2 d kappa_2^2
  //=================================================================

  const double aintegral = (2.*M_PI)
                         *1./(16.*M_PI*M_PI
                         *pow(x1*x2*s_, 2)) * amat2
                         * flux1_/M_PI
                         * flux2_/M_PI
                         *(1./4.)*Constants::GeV2toBarn
                         * 0.5*4./(4.*M_PI);
  if (aintegral*qt1_*qt2_*pt_diff_!=0.) {
    //GenericProcess::DumpPoint(Information);
    //Information(Form("matrix element: %E", aintegral*qt1_*qt2_*pt_diff_));
  }

  //=================================================================
  return aintegral*qt1_*qt2_*pt_diff_;
  //=================================================================
}

void
PPtoLL::fillCentralParticlesKinematics()
{
  // randomise the charge of the outgoing leptons
  int sign = ( drand()>.5 ) ? +1 : -1;

  //=================================================================
  //     first outgoing lepton
  //=================================================================
  Particle* ol1 = particlePtr( Particle::CentralParticle1 );
  ol1->setPdgId( ol1->pdgId(), sign );
  ol1->status = Particle::FinalState;
  if ( !ol1->setMomentum( Pl1_ ) ) { InError( "Invalid outgoing lepton 1" ); }

  //=================================================================
  //     second outgoing lepton
  //=================================================================
  Particle* ol2 = particlePtr( Particle::CentralParticle2 );
  ol2->setPdgId( ol2->pdgId(), -sign );
  ol2->status = Particle::FinalState;
  if ( !ol2->setMomentum( Pl2_ ) ) { InError( "Invalid outgoing lepton 2" ); }
}
