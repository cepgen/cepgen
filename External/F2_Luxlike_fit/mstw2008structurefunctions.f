C--   G.W. 02/12/2008 Standalone code for structure functions used in
C--   MSTW 2008 fits.  Code copied directly from fitting program.
C--   Only neutral-current structure functions in this version,
C--   not charged-current structure function (to be done).
C--   Comments to Graeme Watt <watt(at)hep.ucl.ac.uk>
C--   or Robert Thorne <thorne(at)hep.ucl.ac.uk>.

C--   Input variables for MSTWNC:
C--     x = Bjorken-x value.
C--     q = sqrt(Q^2) in GeV.
C--     IPN = 1, 2 for p or n.
C--     IORD = 0, 1, 2 for LO, NLO, NNLO (pass in COMMON block).
C--   Output variables: f2,f2c,f2b,fl,flc,flb.
      SUBROUTINE MSTWNC(x,q,ipn,f2,f2c,f2b,fl,flc,flb)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION FTEMP(5)
      INTEGER alphaSorder,alphaSnfmax
      DOUBLE PRECISION distance,tolerance,
     &     mCharm,mBottom,alphaSQ0,alphaSMZ
      COMMON/mstwCommon/distance,tolerance,
     &     mCharm,mBottom,alphaSQ0,alphaSMZ,alphaSorder,alphaSnfmax
      COMMON/GRPTHY/FLAVOR
      COMMON/DYLAMB/xlam,S0
      COMMON/iordCommon/iord
      COMMON/GAUS96/XI(96),WI(96),XX(97),NTERMS ! G.W. 15/02/2007
      DATA PI,PI2/3.14159,9.8696/

      IF (IORD.EQ.2) THEN
         CALL MSTWNCnnlo(x,q,ipn,f2,f2c,f2b,fl,flc,flb)
         RETURN
      ELSE IF (IORD.NE.0.AND.IORD.NE.1) THEN
         WRITE(6,*) "Error in MSTWNC, IORD = ",IORD
         STOP
      END IF

      Q2=q**2                   ! photon virtuality
      WW2=(1.-X)*Q2/X+0.88      ! 0.88 is square of proton mass

      xlam = 0.3D0
      Q02 = 1.D0
      S0=LOG(Q02/xlam**2)
      S=LOG(LOG(Q2/xlam**2)/S0)
      qsdt = 4.D0*mCharm**2
      qsct = 4.D0*mBottom**2
      epsc4=qsdt/q2
      fpsc4=qsdt/ww2
      epsc=epsc4/4.D0
      epsb4=qsct/q2
      fpsb4=qsct/ww2
      epsb=epsb4/4.D0
      FLAVOR=3.D0
      FAC=FLAVOR
      CF=4.D0/3.D0
      enf=flavor
      AL1=LOG(1.-X)
      T=S0*EXP(S)
      AL=ALPHA(T)/(4.D0*PI)
      Schm=LOG(LOG(0.25*Qsdt/xlam**2)/S0)
      Sbot=LOG(LOG(0.25*Qsct/xlam**2)/S0)
      tchm=s0*exp(schm)
      tbot=s0*exp(sbot)
      alchm=alpha(tchm)/(4.D0*PI)
      albot=alpha(tbot)/(4.D0*PI)
      ca=3.D0
      al39=al

      CALL FETCH(X,S,IPN,FTEMP)
      theory=ftemp(1)
      fx=ftemp(1)
      fg=ftemp(2)
      fc=ftemp(3)
      fb=ftemp(4)
      fsing=ftemp(5)
      sfac=2./9.
      ffx=0.
      ffc=0.
      ffb=0.
      fflx=0.
      fflc=0. 
      fflb=0.
      if(epsc.gt.1.) fc=0.
      if(epsb.gt.1.) fb=0.

      IF3=0

      FAC=6./9.
      facc=4./9.
      facb=1./9.

      IF (IORD.EQ.0) THEN
      ffx=FX
      ELSE
      ffx=FX+FX*AL39*CF*(-9.-2.*PI2/3.+AL1*(-3.+2.*AL1))
c$$$      fflx=fflx+fx*al**2*(-0.012)
      fflx=fflx+fx*al**2*CLNN2C(x) ! G.W. 02/11/2007
      END IF

      DO 23 I=1,NTERMS
      Y=0.5*(1.-X)*XI(I)+0.5*(1.+X)
      XY=X/Y
      AL1=LOG(1.-Y)
      CALL FETCH(XY,S,IPN,FTEMP)
      if(epsc.gt.1.) ftemp(3)=0.
      if(epsb.gt.1.) ftemp(4)=0.
      IF (IORD.NE.0) THEN
      C22=CF*(6.+4.*Y-2.*(1.+Y*Y)/(1.-Y)*LOG(Y)-2.*(1.+Y)*AL1
     2-IF3*2.*(1.+Y))
      C23=CF*(-3.+4.*AL1)/(1.-Y)
      CG2=2.*FAC*(-1.+8.*Y*(1.-Y)+(1.-2.*Y+2.*Y*Y)*LOG(1./Y-1.))
      END IF
      f1lq=4.*cf*y
      f1lg=8.*fac*y*(1.-y)
      IF (IORD.NE.0) THEN
      ffx=ffx+.5*(1.-X)*WI(I)*AL39*
     2(C22*FTEMP(1)+C23*(FTEMP(1)-FX))
      ffx=ffx+.5*(1.-X)*WI(I)*AL39*CG2*FTEMP(2)
      END IF
      fflx=fflx+.5*(1.-x)*wi(i)*al*f1lq*ftemp(1)
      fflx=fflx+.5*(1.-x)*wi(i)*al*f1lg*ftemp(2)

      IF (IORD.NE.0) THEN
      Y1=1.-Y
      DL=LOG(Y)
      DL2=DL*DL
      DLM1=LOG(Y1)
      DLM2=DLM1*DLM1
      DLM3=DLM2*DLM1
      DLM4=DLM3*DLM1
c$$$      FNS2LQ=128./9.d0*y*DLM1**2-46.50*y*DLM1-84.094*DL*DLM1
c$$$     x-37.338 +89.53*y
c$$$     X+33.82*y**2+y*DL*(32.90+18.41*DL)-128./9.d0*DL
c$$$     X+16./27.d0*flavor*(6.*y*DLM1-12*y*DL-25.*y+6.)
      FNS2LQ = CLNN2A(Y,INT(flavor)) ! G.W. 02/11/2007
      FS2LQ=fac*((15.94-5.212*y)*Y1*Y1*DLM1+(0.421+1.520*y)*DL*DL
     X+28.09*Y1*DL-(2.370/Y-19.27)*Y1**3)
      F2LG=fac*((94.74-49.20*y)*y1*DLM2+864.8*Y1*DLM1
     x+1161.*y*DLM1*DL
     X+60.06*y*DL*DL+39.66*Y1*DL-5.333*(1./Y-1.))
      fflx=fflx+0.5*(1.-x)*WI(I)*AL*AL*
     X(FNS2LQ*ftemp(1)+FS2LQ*ftemp(5)+F2LG*ftemp(2))
      END IF

   23 CONTINUE
   21 CONTINUE

      xcmax=1./(1.+epsc4)
      if(xcmax.le.x) go to 321
      xcmup=x/xcmax
      CALL FETCH(Xcmup,S,IPN,FTEMP)
      if(epsc.gt.1.) ftemp(3)=0.
      fc=ftemp(3)
      IF (IORD.EQ.0) THEN
      ffc=fc
      ELSE
      AL1c=LOG(1.-xcmup)
      ffc=Fc+Fc*AL*CF*(-9.-2.*PI2/3.+AL1c*(-3.+2.*AL1c))
      END IF

      DO 323 I=1,NTERMS
      Y=0.5*(xcmax-X)*XI(I)+0.5*(xcmax+X)
      XY=X/Y
      CALL FETCH(XY,S,IPN,FTEMP)
      if(epsc.gt.1.) ftemp(3)=0.
      gluxy=ftemp(2)
      fcxy=ftemp(3)
      ycmup=y/xcmax
c      if(ycmup.gt.0.9999) ycmup=0.9999
      IF (IORD.NE.0) THEN
      c0c=1.
      p0qg=ycmup**2+(1.-ycmup)**2
      C22c=CF*(6.+4.*Ycmup-2.*(1.+Ycmup*Ycmup)/(1.-Ycmup)*
     .LOG(Ycmup)-2.*(1.+Ycmup)*log(1.-ycmup)
     2-IF3*2.*(1.+Ycmup))
      C23c=CF*(-3.+4.*log(1.-ycmup))/(1.-Ycmup)
      if(epsc.gt.1.d0) c0c=0.d0
      cg21c=2.*facc*cheavy(1,y,epsc)
      cg22c=2.*facc*c0c*p0qg*log(1./epsc)
      END IF
      clg2c=2.*facc*cheavy(7,y,epsc)
      f1lq=cheavy(8,ycmup,epsc)
      IF (IORD.NE.0) THEN
      ffc=ffc+0.5/xcmax*(xcmax-X)*WI(I)*AL*(C22c*fcxy+C23c*(fcxy-fc))      
      ffc=ffc+0.5*(xcmax-x)*wi(i)*al*(cg21c-cg22c/xcmax)*gluxy
      END IF
      fflc=fflc+0.5/xcmax*(xcmax-x)*wi(i)*al*f1lq*fcxy
      fflc=fflc+0.5*(xcmax-x)*wi(i)*al*clg2c*gluxy

  323 CONTINUE

      if(epsc.gt.1.) then 
      xcmax=1./(1.+epsc4)
      eps=epsc
      if(xcmax.le.x) go to 321
      DO 324 I=1,NTERMS
      Y=0.5*(xcmax-X)*XI(I)+0.5*(xcmax+X)
      XY=X/Y
      CALL FETCH(XY,S,IPN,FTEMP)
      gluxy=ftemp(2)
      IF (IORD.EQ.0) THEN
      cg21c=2.*facc*cheavy(1,y,eps)
      ffc=ffc+0.5*(xcmax-x)*wi(i)*al*(cg21c)*gluxy
      ELSE
      singxy=ftemp(5)
      cgff2=facc*(c2gffnsl(y,eps)*(1-0.5*exp(1-eps**2))+
     .c2gffnsh(y,eps)*0.5*exp(1-eps**2))
      cqff2=facc*(c2qffnsl(y,eps)*(1-0.5*exp(1-eps**2))+
     .c2qffnsh(y,eps)*0.5*exp(1-eps**2))
c$$$      ffc=ffc+0.5*(xcmax-x)*wi(i)*al**2.*(cgff2*gluxy+cqff2*singxy)
C--   G.W. 25/07/2007 Add charm and bottom contributions to singlet.
      ffc=ffc+0.5*(xcmax-x)*wi(i)*al**2.*(cgff2*gluxy+cqff2*
     &     ( singxy + 9./4.*ftemp(3) + 9.*ftemp(4) ))
      END IF

  324 CONTINUE

      else 
c$$$      xcmax=1./(1.+ 4.)
      xcmax=1.D0/(1.D0+ 4.D0)
      eps=1.
      if(xcmax.le.x) go to 321
      DO 325 I=1,NTERMS
      Y=0.5*(xcmax-X)*XI(I)+0.5*(xcmax+X)
      XY=X/Y
      CALL FETCH(XY,Schm,IPN,FTEMP)
      gluxy=ftemp(2)
      IF (IORD.EQ.0) THEN
      cg21c=2.*facc*cheavy(1,y,eps)
      ffc=ffc+0.5*(xcmax-x)*wi(i)*alchm*(cg21c)*gluxy
      ELSE
      singxy=ftemp(5)
      cgff2=facc*(c2gffnsl(y,eps)*(1-0.5*exp(1-eps**2))+
     .c2gffnsh(y,eps)*0.5*exp(1-eps**2))
      cqff2=facc*(c2qffnsl(y,eps)*(1-0.5*exp(1-eps**2))+
     .c2qffnsh(y,eps)*0.5*exp(1-eps**2))
c$$$      ffc=ffc+0.5*(xcmax-x)*wi(i)*alchm**2.*(cgff2*gluxy+cqff2*singxy)
C--   G.W. 25/07/2007 Add charm and bottom contributions to singlet.
      ffc=ffc+0.5*(xcmax-x)*wi(i)*alchm**2.*(cgff2*gluxy+cqff2*
     &     ( singxy + 9./4.*ftemp(3) + 9.*ftemp(4) ))
      END IF

  325 CONTINUE
      endif


      IF (IORD.NE.0) THEN
      xcmax=1./(1.+epsc4)
      if(xcmax.le.x) go to 321
      xcmup=x/xcmax
      CALL FETCH(XCMUP,S,IPN,FTEMP)
      if(epsc.gt.1.) ftemp(3)=0. ! G.W. 05/11/2007
c$$$      fflc=fflc+ftemp(3)*(AL**2*(-0.0012))
c$$$     .*1.25*(1/(1+4.*epsc)-0.2)
      fflc=fflc+ftemp(3)*(AL**2*CLNN2C(xcmup))
     &     *1.25*(1/(1+4.*epsc)-0.2) ! G.W. 05/11/2007

 529  continue      
      if(epsc.gt.1.) then 
      xcmax=1./(1.+epsc4)
      eps=epsc
      if(xcmax.le.x) go to 321
      DO 524 I=1,NTERMS
      Y=0.5*(xcmax-X)*XI(I)+0.5*(xcmax+X)
      XY=X/Y
      CALL FETCH(XY,S,IPN,FTEMP)
      gluxy=ftemp(2)
      singxy=ftemp(5)
      cgffl=facc*(clgffnsl(y,eps)*(1-0.5*exp(1-eps**2))+
     .clgffnsh(y,eps)*0.5*exp(1-eps**2)) 
      cqffl=facc*(clqffnsl(y,eps)*(1-0.5*exp(1-eps**2))+
     .clqffnsh(y,eps)*0.5*exp(1-eps**2))
c$$$      fflc=fflc+0.5*(xcmax-x)*wi(i)*al**2.*(cgffl*gluxy+cqffl*singxy)
C--   G.W. 25/07/2007 Add charm and bottom contributions to singlet.
      fflc=fflc+0.5*(xcmax-x)*wi(i)*al**2.*(cgffl*gluxy+cqffl*
     &     ( singxy + 9./4.*ftemp(3) + 9.*ftemp(4) ))

  524 CONTINUE
      else 
      xcmax=(1./(1.+epsc4)) 
      eps=epsc
      if(xcmax.le.x) go to 321
      DO 525 I=1,NTERMS
      Y=0.5*(xcmax-X)*XI(I)+0.5*(xcmax+X)
      XY=X/Y
      CALL FETCH(XY,S,IPN,FTEMP)
      if(epsc.gt.1.) ftemp(3)=0.
      gluxy=ftemp(2)
      fcxy=ftemp(3)
      singxy=ftemp(5)
      ymul=y*(1+epsc4)
      Y1mul=1.-Ymul
      DL=LOG(Ymul)
      DL2=DL*DL
      DLM1=LOG(Y1mul)
      DLM2=DLM1*DLM1
      DLM3=DLM2*DLM1
      DLM4=DLM3*DLM1
c$$$      FNS2LQmul=128./9.d0*ymul*DLM1**2-46.50*ymul*DLM1-84.094*DL*DLM1
c$$$     x-37.338 +89.53*ymul
c$$$     X+33.82*ymul**2+ymul*DL*(32.90+18.41*DL)-128./9.d0*DL
c$$$     X+16./27.d0*enf*(6.*ymul*DLM1-12*ymul*DL-25.*ymul+6.)
      FNS2LQmul = CLNN2A(YMUL,INT(enf))   ! G.W. 02/11/2007
      FS2LQmul=((15.94-5.212*ymul)*Y1mul*Y1mul*DLM1+(0.421+1.520*ymul)
     x*DL*DL+28.09*Y1mul*DL-(2.370/Ymul-19.27)*Y1mul**3)
      cgvfl=facc*((clgffnsh(y,eps)*(1-0.5*exp(1-1/eps**2))+
     .clgffnsl(y,eps)*0.5*exp(1-1/eps**2))-clgvfsub(ymul,eps)/xcmax)
      cqvfl=facc*((clqffnsh(y,eps)*(1-0.5*exp(1-1/eps**2))+
     .clqffnsl(y,eps)*0.5*exp(1-1/eps**2)))
      fflc=fflc+0.5/xcmax*(xcmax-x)*wi(i)*al**2.*(fcxy*FNS2LQmul+
     xfcxy*FS2LQmul)*1.25*(1/(1+4.*eps)-0.2)
c$$$      fflc=fflc+0.5*(xcmax-x)*wi(i)*al**2.*(cgvfl*gluxy+cqvfl*singxy)
C--   G.W. 25/07/2007 Add charm and bottom contributions to singlet.
      fflc=fflc+0.5*(xcmax-x)*wi(i)*al**2.*(cgvfl*gluxy+cqvfl*
     &     ( singxy + 9./4.*ftemp(3) + 9.*ftemp(4) ))

  525 CONTINUE
      endif
      END IF


  321 CONTINUE

      IF (IORD.NE.0) THEN
      if(ffc.lt.0.) ffc=0.
      END IF

      xbmax=1./(1.+epsb4)
      if(xbmax.le.x+0.00001) go to 421
      xbmup=x/xbmax
      CALL FETCH(Xbmup,S,IPN,FTEMP)
      if(epsb.gt.1.) ftemp(4)=0.
      fb=ftemp(4)
      IF (IORD.EQ.0) THEN       ! G.W. 05/11/2008
      ffb=fb                    ! G.W. 05/11/2008
      ELSE                      ! G.W. 05/11/2008
      AL1b=LOG(1.-xbmup)
      ffb=Fb+Fb*AL*CF*(-9.-2.*PI2/3.+AL1b*(-3.+2.*AL1b))
      END IF

      DO 423 I=1,NTERMS
      Y=0.5*(xbmax-X)*XI(I)+0.5*(xbmax+X)
      XY=X/Y
      CALL FETCH(XY,S,IPN,FTEMP)
      if(epsb.gt.1.) ftemp(4)=0.
      gluxy=ftemp(2)
      fbxy=ftemp(4)
      ybmup=y/xbmax
      if(ybmup.gt.0.99999d0) ybmup=0.99999d0
      IF (IORD.NE.0) THEN
      c0b=1.
      p0qg=ybmup**2+(1.-ybmup)**2
      C22b=CF*(6.+4.*Ybmup-2.*(1.+Ybmup*Ybmup)/(1.-Ybmup)*
     .LOG(Ybmup)-2.*(1.+Ybmup)*log(1.-ybmup)
     2-IF3*2.*(1.+Ybmup))
      C23b=CF*(-3.+4.*log(1.-ybmup))/(1.-Ybmup)
      if(epsb.gt.1.d0) c0b=0.d0
      cg21b=2.*facb*cheavy(1,y,epsb)
      cg22b=2.*facb*c0b*p0qg*log(1./epsb)
      ffb=ffb+0.5/xbmax*(xbmax-X)*WI(I)*AL*(C22b*fbxy+C23b*(fbxy-fb))      
      ffb=ffb+0.5*(xbmax-x)*wi(i)*al*(cg21b-cg22b/xbmax)*gluxy
      END IF
      clg2b=2.*facb*cheavy(7,y,epsb)
      f1lq=cheavy(8,ybmup,epsb)
      fflb=fflb+0.5/xbmax*(xbmax-x)*wi(i)*al*f1lq*fbxy
      fflb=fflb+0.5*(xbmax-x)*wi(i)*al*clg2b*gluxy


  423 CONTINUE

      if(epsb.gt.1.) then 
      xbmax=1./(1.+epsb4)
      eps=epsb
      if(xbmax.le.x) go to 421
      DO 424 I=1,NTERMS
      Y=0.5*(xbmax-X)*XI(I)+0.5*(xbmax+X)
      XY=X/Y
      CALL FETCH(XY,S,IPN,FTEMP)
      gluxy=ftemp(2)
      IF (IORD.EQ.0) THEN
      cg21b=2.*facb*cheavy(1,y,eps)
      ffb=ffb+0.5*(xbmax-x)*wi(i)*al*(cg21b)*gluxy
      ELSE
c$$$      singxy=ftemp(5)+9./8.*ftemp(3) ! Why 9/8 and not 9/4?
      singxy=ftemp(5)           ! G.W. 27/07/2007
      cgff2=facb*(c2gffnsl(y,eps)*(1-0.5*exp(1-eps**2))+
     .c2gffnsh(y,eps)*0.5*exp(1-eps**2))
      cqff2=facb*(c2qffnsl(y,eps)*(1-0.5*exp(1-eps**2))+
     .c2qffnsh(y,eps)*0.5*exp(1-eps**2))
c$$$      ffb=ffb+0.5*(xbmax-x)*wi(i)*al**2.*(cgff2*gluxy+cqff2*singxy)
C--   G.W. 25/07/2007 Add charm and bottom contributions to singlet.
      ffb=ffb+0.5*(xbmax-x)*wi(i)*al**2.*(cgff2*gluxy+cqff2*
     &     ( singxy + 9./4.*ftemp(3) + 9.*ftemp(4) ))
      END IF
      
  424 CONTINUE

      else 
c$$$      xbmax=1./(1.+ 4.)
      xbmax=1.D0/(1.D0+ 4.D0)
      eps=1.
      if(xbmax.le.x) go to 421
      DO 425 I=1,NTERMS
      Y=0.5*(xbmax-X)*XI(I)+0.5*(xbmax+X)
      XY=X/Y
      CALL FETCH(XY,Sbot,IPN,FTEMP)
      gluxy=ftemp(2)
      IF (IORD.EQ.0) THEN
      cg21b=2.*facb*cheavy(1,y,eps)
      ffb=ffb+0.5*(xbmax-x)*wi(i)*albot*(cg21b)*gluxy
      ELSE
c$$$      singxy=ftemp(5)+9./8.*ftemp(3) ! Why 9/8 and not 9/4?
      singxy=ftemp(5)           ! G.W. 27/07/2007
      cgff2=facb*(c2gffnsl(y,eps)*(1-0.5*exp(1-eps**2))+
     .c2gffnsh(y,eps)*0.5*exp(1-eps**2))
      cqff2=facb*(c2qffnsl(y,eps)*(1-0.5*exp(1-eps**2))+
     .c2qffnsh(y,eps)*0.5*exp(1-eps**2))
c$$$      ffb=ffb+0.5*(xbmax-x)*wi(i)*albot**2.*(cgff2*gluxy+cqff2*singxy)
C--   G.W. 25/07/2007 Add charm and bottom contributions to singlet.
      ffb=ffb+0.5*(xbmax-x)*wi(i)*albot**2.*(cgff2*gluxy+cqff2*
     &     ( singxy + 9./4.*ftemp(3) + 9.*ftemp(4) ))
      END IF

  425 CONTINUE
      endif

      IF (IORD.NE.0) THEN

C--   G.W. 05/11/2007 This contribution was missing for b, but included for c.
      xbmax=1./(1.+epsb4)
      if(xbmax.le.x) go to 421
      xbmup=x/xbmax
      CALL FETCH(XBMUP,S,IPN,FTEMP)
      if(epsb.gt.1.) ftemp(4)=0. ! G.W. 05/11/2007
      fflb=fflb+ftemp(4)*(AL**2*CLNN2C(xbmup))
     &     *1.25*(1/(1+4.*epsb)-0.2) ! G.W. 05/11/2007

 629  continue
      if(epsb.gt.1.) then 
      xbmax=1./(1.+epsb4)
      eps=epsb
      if(xbmax.le.x) go to 421
      DO 624 I=1,NTERMS
      Y=0.5*(xbmax-X)*XI(I)+0.5*(xbmax+X)
      XY=X/Y
      CALL FETCH(XY,S,IPN,FTEMP)
      gluxy=ftemp(2)
      singxy=ftemp(5)
      cgffl=facb*(clgffnsl(y,eps)*(1-0.5*exp(1-eps**2))+
     .clgffnsh(y,eps)*0.5*exp(1-eps**2)) 
      cqffl=facb*(clqffnsl(y,eps)*(1-0.5*exp(1-eps**2))+
     .clqffnsh(y,eps)*0.5*exp(1-eps**2))
c$$$      fflb=fflb+0.5*(xbmax-x)*wi(i)*al**2.*(cgffl*gluxy+cqffl*singxy)
C--   G.W. 25/07/2007 Add charm and bottom contributions to singlet.
      fflb=fflb+0.5*(xbmax-x)*wi(i)*al**2.*(cgffl*gluxy+cqffl*
     &     ( singxy + 9./4.*ftemp(3) + 9.*ftemp(4) ))

  624 CONTINUE
      else 
      xbmax=(1./(1.+epsb4)) 
      eps=epsb
      if(xbmax.le.x) go to 421
      DO 625 I=1,NTERMS
      Y=0.5*(xbmax-X)*XI(I)+0.5*(xbmax+X)
      XY=X/Y
      CALL FETCH(XY,S,IPN,FTEMP)
      if(epsb.gt.1.) ftemp(3)=0.
      gluxy=ftemp(2)
      fbxy=ftemp(4)
      singxy=ftemp(5)
      ymul=y*(1+epsb4)
      Y1mul=1.-Ymul
      DL=LOG(Ymul)
      DL2=DL*DL
      DLM1=LOG(Y1mul)
      DLM2=DLM1*DLM1
      DLM3=DLM2*DLM1
      DLM4=DLM3*DLM1
c$$$      FNS2LQmul=128./9.d0*ymul*DLM1**2-46.50*ymul*DLM1-84.094*DL*DLM1
c$$$     x-37.338 +89.53*ymul
c$$$     X+33.82*ymul**2+ymul*DL*(32.90+18.41*DL)-128./9.d0*DL
c$$$     X+16./27.d0*enf*(6.*ymul*DLM1-12*ymul*DL-25.*ymul+6.)
      FNS2LQmul = CLNN2A(YMUL,INT(ENF)) ! G.W. 02/11/2007
      FS2LQmul=((15.94-5.212*ymul)*Y1mul*Y1mul*DLM1+(0.421+1.520*ymul)
     x*DL*DL+28.09*Y1mul*DL-(2.370/Ymul-19.27)*Y1mul**3)
      cgvfl=facb*((clgffnsh(y,eps)*(1-0.5*exp(1-1/eps**2))+
     .clgffnsl(y,eps)*0.5*exp(1-1/eps**2))-clgvfsub(ymul,eps)/xcmax)
      cqvfl=facb*((clqffnsh(y,eps)*(1-0.5*exp(1-1/eps**2))+
     .clqffnsl(y,eps)*0.5*exp(1-1/eps**2)))
      fflb=fflb+0.5/xbmax*(xbmax-x)*wi(i)*al**2.*(fbxy*FNS2LQmul+
     xfbxy*FS2LQmul)*1.25*(1/(1+4.*eps)-0.2)
c$$$      fflb=fflb+0.5*(xbmax-x)*wi(i)*al**2.*(cgvfl*gluxy+cqvfl*singxy)
C--   G.W. 25/07/2007 Add charm and bottom contributions to singlet.
      fflb=fflb+0.5*(xbmax-x)*wi(i)*al**2.*(cgvfl*gluxy+cqvfl*
     &     ( singxy + 9./4.*ftemp(3) + 9.*ftemp(4) ))

  625 CONTINUE
      endif
      END IF

  421 CONTINUE

      if(ffb.lt.0.) ffb=0.

      F2 = ffx+ffc+ffb
      F2C = ffc
      F2B = ffb
      FL = fflx+fflc+fflb
      FLC = fflc
      FLB = fflb

      RETURN
      END
C--   End of MSTWNC subroutine.



      SUBROUTINE MSTWNCnnlo(x,q,ipn,f2,f2c,f2b,fl,flc,flb)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION FTEMP(5)
      INTEGER alphaSorder,alphaSnfmax
      DOUBLE PRECISION distance,tolerance,
     &     mCharm,mBottom,alphaSQ0,alphaSMZ
      COMMON/mstwCommon/distance,tolerance,
     &     mCharm,mBottom,alphaSQ0,alphaSMZ,alphaSorder,alphaSnfmax
      COMMON/GRPTHY/FLAVOR
      COMMON/DYLAMB/xlam,S0
      COMMON/iordCommon/iord
      COMMON/GAUS96/XI(96),WI(96),XX(97),NTERMS ! G.W. 15/02/2007
      DATA PI,PI2/3.14159,9.8696/

      IF (IORD.NE.2) THEN
         WRITE(6,*) "Error in MSTWNCnnlo, IORD = ",IORD
         STOP
      END IF

      Q2=q**2                   ! photon virtuality
      WW2=(1.-X)*Q2/X+0.88      ! 0.88 is square of proton mass

      xlam = 0.3D0
      Q02 = 1.D0
      S0=LOG(Q02/xlam**2)
      S=LOG(LOG(Q2/xlam**2)/S0)
      qsdt = 4.D0*mCharm**2
      qsct = 4.D0*mBottom**2
      epsc4=qsdt/q2
      fpsc4=qsdt/ww2
      epsc=epsc4/4.D0
      epsb4=qsct/q2
      fpsb4=qsct/ww2
      epsb=epsb4/4.D0
      FLAVOR=3.D0
      FAC=FLAVOR
      CF=4.D0/3.D0
      enf=flavor
      AL1=LOG(1.D0-X)
      T=S0*EXP(S)
      AL=ALPHA(T)/(4.D0*PI)
      Schm=LOG(LOG(0.25*Qsdt/xlam**2)/S0)
      Sbot=LOG(LOG(0.25*Qsct/xlam**2)/S0)
      tchm=s0*exp(schm)
      tbot=s0*exp(sbot)
      alchm=alpha(tchm)/(4.D0*PI)
      albot=alpha(tbot)/(4.D0*PI)
      ca=3.D0
      zz2=1.6449D0
      zz3=1.20205D0
      al39=al

      CALL FETCH(X,S,IPN,FTEMP)
      theory=ftemp(1)
      fx=ftemp(1)
      fg=ftemp(2)
      fc=ftemp(3)
      fb=ftemp(4)
      fsing=ftemp(5)
      sfac=2./9.
      ffx=0.
      ffc=0.
      ffb=0.
      fflx=0.
      fflc=0. 
      fflb=0.
      if(epsc.gt.1.) fc=0.
      if(epsb.gt.1.) fb=0.

      IF3=0

      FAC=6./9.
      facc=4./9.
      facb=1./9.

      ffx=FX+FX*AL39*CF*(-9.-2.*PI2/3.+AL1*(-3.+2.*AL1))
c$$$      fflx=fflx+fx*al**2*(-0.012)+fx*al**3*(0.113+enf*0.006)
      fflx=fflx+fx*al**2*CLNN2C(x)+fx*al**3*(0.113+enf*0.006) ! G.W. 02/11/2007
      if(iord.gt.1) then
      ffx=ffx+fx*al39*al39*c2nn2c(x,3)
      ffx=ffx+fg*al39*al39*sfac*c2g2c(x,3)
      endif

      if(epsc.gt.1.) go to 281
      if(epsc.lt.1.) go to 282
 

 281  ffx=ffx+fx*2./3.*al*al*1./epsc*(352./225.+16./45.*log(epsc))
      go to 207

 282  if(epsc.lt.0.001) go to 283
      AL1=log(1.-X)
      ffx=ffx+fx*al*al*2./3.*(+8./9.*(log(epsc))**3.
     .+76./9.*(log(epsc))**2.-(32./3.*zz2-1060/27.)*log(epsc)
     .+6710./81-304./9.*zz2-32./3.*zz3+
     .(78.7+26.32*log(epsc))*(1.-epsc*(1.-log(epsc)))
     .+2.*epsc*(-6.44-3.147*log(epsc)+7.266*(log(epsc))**2.
     .+0.1064*(log(epsc))**4))
      ffx=ffx-fx*2./3.*al*al*((8./3.*(log(epsc))**2.
     .+80./9.*log(epsc)+224./27.)*al1+(2*(log(epsc))**2.
     .+(16./3.*zz2+2./3.)*log(epsc)
     .+(73./18.+40./9.*zz2-8./3.*zz3)))
      ffx=ffx+fx*2./3.*log(epsc)*al*al*CF*(-9.-2.*PI2/3.
     .+AL1*(-3.+2.*AL1))
      DO 201 I=1,NTERMS
      z=0.5*(1.-X)*XI(I)+0.5*(1.+X)
      aL1=log(1.-Z)
      aLZ=log(Z)
      xz=x/z 
      CALL FETCH(XZ,S,IPN,FTEMP)
      C22=4./3.*(6.+4.*z-2.*(1.+z*z)/(1.-z)*log(z)-2.*(1.+z)*AL1)
      C23=4./3.*(-3.+4.*AL1)/(1.-z)
      ffx=ffx-0.5*(1.-x)*al*al*wi(i)*2./3.*
     .((8./3.*(log(epsc))**2.+80./9.*log(epsc)+224./27.)*
     .(ftemp(1)-fx)/(1.-z)+((8./3.*log(z)*log(epsc)+(2./3.*alz*alz
     .+20./9.*alz))*(1.+z*z)/(1.-z)-4./3.*(1.+z)*(log(epsc))**2.
     .+(8./9.-88./9.*z)*log(epsc)
     .+(8./3.*(1.-z)*alz+44./27.-268./27.*z))*ftemp(1))
      ffx=ffx+2./3.*log(epsc)*0.5*(1.-X)*WI(I)*AL*al*
     .(C22*ftemp(1)+C23*(ftemp(1)-fx))
 201  continue
      go to 207

 283  continue
      ffx=ffx+fx*al*al*c2nn2ch(x)  
      DO 202 I=1,NTERMS
      z=0.5*(1.-X)*XI(I)+0.5*(1.+X)
      xz=x/z 
      CALL FETCH(XZ,S,IPN,FTEMP)
      ffx=ffx+0.5*(1.-x)*wi(i)*al*al*
     .(ftemp(1)*c2nn2ah(z) + (ftemp(1)-fx)*c2ns2bh(z))
 202  continue
      go to 207
       
  
 207  continue

      if(epsc.gt.1.) go to 291 
      if(epsc.lt.1) go to 292

 291  xcmax=1./(1.+epsc4)
      eps=epsc
      if(xcmax.le.x) go to 221
      DO 224 I=1,NTERMS
      Y=0.5*(xcmax-X)*XI(I)+0.5*(xcmax+X)
      XY=X/Y
      CALL FETCH(XY,S,IPN,FTEMP)
      cl2ff2l=cl2ffnsl(y,eps)
      ffx=ffx+0.5*(xcmax-x)*wi(i)*al**2.*(cl2ff2l*ftemp(1))
  224 CONTINUE
      go to 221
      
 292  if(epsc.lt.0.001) go to 293
      xcmax=(1./(1.+epsc4)) 
      eps=epsc
      if(xcmax.le.x) go to 221
      DO 225 I=1,NTERMS
      Y=0.5*(xcmax-X)*XI(I)+0.5*(xcmax+X)
      XY=X/Y
      CALL FETCH(XY,S,IPN,FTEMP)
      cl2ff2h=cl2ffnsh(y,eps)
      ffx=ffx+0.5*(xcmax-x)*wi(i)*al**2.*(cl2ff2h*ftemp(1))
  225 CONTINUE
      go to 221 

  293 ffx=ffx+0.
      go to 221
      

  221 continue

      if(epsb.gt.1.) go to 2811
      if(epsb.lt.1.) go to 2821
 

 2811 ffx=ffx+fx*2./3.*al*al*1./epsb*(352./225.+16./45.*log(epsb))
      go to 2071

 2821 if(epsb.lt.0.001) go to 2831
      AL1=log(1.-X)
      ffx=ffx+fx*al*al*2./3.*(+8./9.*(log(epsb))**3.
     .+76./9.*(log(epsb))**2.-(32./3.*zz2-1060/27.)*log(epsb)
     .+6710./81-304./9.*zz2-32./3.*zz3+
     .(78.7+26.32*log(epsb))*(1.-epsb*(1.-log(epsb)))
     .+2.*epsb*(-6.44-3.147*log(epsb)+7.266*(log(epsb))**2.
     .+0.1064*(log(epsb))**4))
      ffx=ffx-fx*2./3.*al*al*((8./3.*(log(epsb))**2.
     .+80./9.*log(epsb)+224./27.)*al1+(2*(log(epsb))**2.
     .+(16./3.*zz2+2./3.)*log(epsb)
     .+(73./18.+40./9.*zz2-8./3.*zz3)))
      ffx=ffx+fx*2./3.*log(epsb)*al*al*CF*(-9.-2.*PI2/3.
     .+AL1*(-3.+2.*AL1))
      DO 2011 I=1,NTERMS
      z=0.5*(1.-X)*XI(I)+0.5*(1.+X)
      aL1=log(1.-Z)
      aLZ=log(Z)
      xz=x/z 
      CALL FETCH(XZ,S,IPN,FTEMP)
      C22=4./3.*(6.+4.*z-2.*(1.+z*z)/(1.-z)*log(z)-2.*(1.+z)*AL1)
      C23=4./3.*(-3.+4.*AL1)/(1.-z)
      fpxz=ftemp(1)
      ffx=ffx-0.5*(1.-x)*al*al*wi(i)*2./3.*
     .((8./3.*(log(epsb))**2.+80./9.*log(epsb)+224./27.)*
     .(fpxz-fx)/(1.-z)+((8./3.*log(z)*log(epsb)+(2./3.*alz*alz
     .+20./9.*alz))*(1.+z*z)/(1.-z)-4./3.*(1.+z)*(log(epsb))**2.
     .+(8./9.-88./9.*z)*log(epsb)
     .+(8./3.*(1.-z)*alz+44./27.-268./27.*z))*fpxz)
      ffx=ffx+2./3.*log(epsb)*0.5*(1.-X)*WI(I)*AL*al*
     .(C22*fpxz+C23*(fpxz-fx))
 2011 continue
      go to 2071

 2831 continue
      ffx=ffx+fx*al*al*c2nn2ch(x)  
      DO 2021 I=1,NTERMS
      z=0.5*(1.-X)*XI(I)+0.5*(1.+X)
      xz=x/z 
      CALL FETCH(XZ,S,IPN,FTEMP)
      fpxz=ftemp(1)
      ffx=ffx+0.5*(1.-x)*wi(i)*al*al*
     .(fpxz*c2nn2ah(z) + (fpxz-fx)*c2ns2bh(z))
 2021 continue
      go to 2071

  
 2071 continue

      if(epsb.gt.1.) go to 2911 
      if(epsb.lt.1.) go to 2921

 2911 xbmax=1./(1.+epsb4)
      eps=epsb
      if(xbmax.le.x) go to 2211
      DO 2241 I=1,NTERMS
      Y=0.5*(xbmax-X)*XI(I)+0.5*(xbmax+X)
      XY=X/Y
      CALL FETCH(XY,S,IPN,FTEMP)
      fpxy=ftemp(1)
      cl2ff2l=cl2ffnsl(y,eps)
      ffx=ffx+0.5*(xbmax-x)*wi(i)*al**2.*(cl2ff2l*fpxy)
 2241 CONTINUE
      go to 2211
      
 2921 if(epsb.lt.0.001) go to 2931
      xbmax=(1./(1.+epsb4)) 
      eps=epsb
      if(xbmax.le.x) go to 2211
      DO 2251 I=1,NTERMS
      Y=0.5*(xbmax-X)*XI(I)+0.5*(xbmax+X)
      XY=X/Y
      CALL FETCH(XY,S,IPN,FTEMP)
      fpxy=ftemp(1)
      cl2ff2h=cl2ffnsh(y,eps)
      ffx=ffx+0.5*(xbmax-x)*wi(i)*al**2.*(cl2ff2h*fpxy)
 2251 CONTINUE
      go to 2211 

 2931 ffx=ffx+0.
      go to 2211
      

 2211 continue



      DO 23 I=1,NTERMS
      Y=0.5*(1.-X)*XI(I)+0.5*(1.+X)
      XY=X/Y
      AL1=LOG(1.-Y)
      CALL FETCH(XY,S,IPN,FTEMP)
      if(epsc.gt.1.) ftemp(3)=0.
      if(epsb.gt.1.) ftemp(4)=0.
      C22=CF*(6.+4.*Y-2.*(1.+Y*Y)/(1.-Y)*LOG(Y)-2.*(1.+Y)*AL1
     2-IF3*2.*(1.+Y))
      C23=CF*(-3.+4.*AL1)/(1.-Y)
      CG2=2.*FAC*(-1.+8.*Y*(1.-Y)+(1.-2.*Y+2.*Y*Y)*LOG(1./Y-1.))
      f1lq=4.*cf*y
      f1lg=8.*fac*y*(1.-y)
      ffx=ffx+.5*(1.-X)*WI(I)*AL39*
     2(C22*FTEMP(1)+C23*(FTEMP(1)-FX))
      ffx=ffx+.5*(1.-X)*WI(I)*AL39*CG2*FTEMP(2)
      fflx=fflx+.5*(1.-x)*wi(i)*al*f1lq*ftemp(1)
      fflx=fflx+.5*(1.-x)*wi(i)*al*f1lg*ftemp(2)
      if(iord.gt.1) then
c$$$      ffx=ffx+0.5*(1.-x)*wi(i)*al39*al39*sfac*ftemp(5)*c2s2a(y,3)
C--   G.W. 25/07/2007 Add charm and bottom contributions to singlet.
      ffx=ffx+0.5*(1.-x)*wi(i)*al39*al39*sfac*c2s2a(y,3)*
     &     ( ftemp(5) + 9./4.*ftemp(3) + 9.*ftemp(4) )
      termm=+0.5*(1.-x)*wi(i)*al39*al39*
     .(ftemp(1)*c2nn2a(y,3) + (ftemp(1)-fx)*c2ns2b(y,3))
      ffx=ffx+termm
      ffx=ffx+0.5*(1.-x)*wi(i)*al39*al39*sfac*ftemp(2)*c2g2a(y,3)
      endif

      Y1=1.-Y
      DL=log(Y)
      DL2=DL*DL
      DLM1=log(Y1)
      DLM2=DLM1*DLM1
      DLM3=DLM2*DLM1
      DLM4=DLM3*DLM1
c$$$      FNS2LQ=128./9.d0*y*DLM1**2-46.50*y*DLM1-84.094*DL*DLM1
c$$$     x-37.338 +89.53*y
c$$$     X+33.82*y**2+y*DL*(32.90+18.41*DL)-128./9.d0*DL
c$$$     X+16./27.d0*flavor*(6.*y*DLM1-12*y*DL-25.*y+6.)
      FNS2LQ = CLNN2A(Y,INT(flavor)) ! G.W. 02/11/2007
      FS2LQ=fac*((15.94-5.212*y)*Y1*Y1*DLM1+(0.421+1.520*y)*DL*DL
     X+28.09*Y1*DL-(2.370/Y-19.27)*Y1**3)
      F2LG=fac*((94.74-49.20*y)*y1*DLM2+864.8*Y1*DLM1
     x+1161.*y*DLM1*DL
     X+60.06*y*DL*DL+39.66*Y1*DL-5.333*(1./Y-1.))
      fflx=fflx+0.5*(1.-x)*WI(I)*AL*AL*
     X(FNS2LQ*ftemp(1)+FS2LQ*ftemp(5)+F2LG*ftemp(2))

      c3lg1=((144.*DLM4-47024./27.*DLM3+6319.*DLM2
     x+53160.*DLM1)*Y1
     X+72549*DL*DLM1+88238*DL*DL*DLM1+(3709.-33514.*Y
     x-9533.*Y*Y)*Y1
     X+66773.*Y*DL*DL -1117.*DL+45.37*DL*DL-5360./27.*DL**3
     X-(2044.7*Y1+409.506*DL)/Y)
      c3lg2=enf*((32./3.*DLM3-1216./9.*DLM2-592.3*DLM1
     X+1511.*y*DLM1)*Y1+311.3*DL*DLM1+14.24*DL*DL*DLM1
     X+(577.3-729.0*y)*Y1+30.78*y*DL**3+366.0*DL+1000./9*DL*DL
     X+160./9.*DL**3+88.5037*y1/y)

      c3lq1=512./27.*DLM4-177.4*DLM3+650.6*DLM2-2729.*DLM1
     x-2220.5-7884.*y+4168*y*Y-(844.7*DL+517.3*DLM1)*DL*DLM1
     x+(195.6*DLM1-125.3)*y1*DLM3+208.3*y*DL**3-1355.7*DL
     x-7456./27.*DL*DL-1280./81*DL**3
      c3lq2=enf*(1024./81.*DLM3-112.35*DLM2+344.1*DLM1+408.4
     x-9.345*y-919.3*Y*Y+(239.7+20.63*DLM1)*y1*DLM2+(887.3
     x+294.5*DL-59.14*DLM1)*DL*DLM1-1792./81.*y*DL**3+200.73*DL
     x+64./3.*DL*DL )
      c3lq3=enf*enf*(3.*y*DLM2+(6.-25*y)*DLM1-19.+(317./6.
     x-12*1.645)*y-6.*y*DL*DLM1+6.*y*(y+y**2/4.+y**3/9.+y**4/16.
     x+y**5/25.+y**6/36.)+9.*y*DL*DL-(6.-50.*y)*DL)*64./81.
      c3lq4=((1568./27.*DLM3-3968./9.*DLM2+5124*DLM1)*y1*y1
     X+(2184.*DL+6059.*y1)*DL*DLM1-(795.6+1036.*y)*y1*y1
     x-143.6*y1*DL+2848./9.*DL*DL-1600./27.*DL**3
     x-(885.53*y1+182.0*DL)*y1/y)
      c3lq5=enf*((-32./9.*DLM2+29.52*DLM1)*y1*y1+(35.18*DL
     X+73.06*y1)*DL*DLM1-35.24*y*DL*DL-(14.16-69.84*y)*y1*y1
     X-69.41*y1*DL-128./9.*DL*DL+40.239*y1*y1/y)
      c3lg=c3lg1+c3lg2
      c3lqns=c3lq1+c3lq2+c3lq3
      c3lqps=c3lq4+c3lq5
c$$$      fflx=fflx+0.5*(1.-x)*wi(i)*al*al*al*
c$$$     .(c3lqns*ftemp(1)+fac*c3lqps*ftemp(5)+c3lg*fac*ftemp(2))
C--   G.W. 25/07/2007 Add charm and bottom contributions to singlet.
      fflx=fflx+0.5*(1.-x)*wi(i)*al*al*al*(c3lqns*ftemp(1)+
     &     fac*c3lqps*(ftemp(5)+9./4.*ftemp(3)+9.*ftemp(4))+
     &     c3lg*fac*ftemp(2))


   23 CONTINUE
   21 CONTINUE


      xcmax=1./(1.+epsc4)
      if(xcmax.le.x) go to 321
      xcmup=x/xcmax
      CALL FETCH(Xcmup,S,IPN,FTEMP)
      if(epsc.gt.1.) ftemp(3)=0.
      fc=ftemp(3)
      AL1c=log(1.-xcmup)
      ffc=Fc+Fc*AL*CF*(-9.-2.*PI2/3.+AL1c*(-3.+2.*AL1c))
      ffc=ffc+fc*al*al*c2nn2c(xcmup,3)

      DO 323 I=1,NTERMS
      Y=0.5*(xcmax-X)*XI(I)+0.5*(xcmax+X)
      XY=X/Y
      CALL FETCH(XY,S,IPN,FTEMP)
      if(epsc.gt.1.) ftemp(3)=0.
      gluxy=ftemp(2)
      fcxy=ftemp(3)
      ycmup=y/xcmax
c      if(ycmup.gt.0.9999d0) ycmup=0.9999d0
      c0c=1.
      p0qg=ycmup**2+(1.-ycmup)**2
      C22c=CF*(6.+4.*Ycmup-2.*(1.+Ycmup*Ycmup)/(1.-Ycmup)*
     .log(Ycmup)-2.*(1.+Ycmup)*log(1.-ycmup)
     2-IF3*2.*(1.+Ycmup))
      C23c=CF*(-3.+4.*log(1.-ycmup))/(1.-Ycmup)
      if(epsc.gt.1.d0) c0c=0.d0
      cg21c=2.*facc*cheavy(1,y,epsc)
      cg22c=2.*facc*c0c*p0qg*log(1./epsc)
      clg2c=2.*facc*cheavy(7,y,epsc)
      f1lq=cheavy(8,ycmup,epsc)
      ffc=ffc+0.5/xcmax*(xcmax-X)*WI(I)*AL*(C22c*fcxy+C23c*(fcxy-fc))      
      ffc=ffc+0.5*(xcmax-x)*wi(i)*al*(cg21c-cg22c/xcmax)*gluxy
      fflc=fflc+0.5/xcmax*(xcmax-x)*wi(i)*al*f1lq*fcxy
      fflc=fflc+0.5*(xcmax-x)*wi(i)*al*clg2c*gluxy

  323 CONTINUE

      if(epsc.gt.1.) then 
      xcmax=1./(1.+epsc4)
      eps=epsc
      if(xcmax.le.x) go to 321
      DO 324 I=1,NTERMS
      Y=0.5*(xcmax-X)*XI(I)+0.5*(xcmax+X)
      XY=X/Y
      CALL FETCH(XY,S,IPN,FTEMP)
      gluxy=ftemp(2)
      singxy=ftemp(5)
      cgff2=facc*(c2gffnsl(y,eps)*(1-0.5*exp(1-eps**2))+
     .c2gffnsh(y,eps)*0.5*exp(1-eps**2))
      cqff2=facc*(c2qffnsl(y,eps)*(1-0.5*exp(1-eps**2))+
     .c2qffnsh(y,eps)*0.5*exp(1-eps**2))
c$$$      ffc=ffc+0.5*(xcmax-x)*wi(i)*al**2.*(cgff2*gluxy+cqff2*singxy)
C--   G.W. 25/07/2007 Add charm and bottom contributions to singlet.
      ffc=ffc+0.5*(xcmax-x)*wi(i)*al**2.*(cgff2*gluxy+cqff2*
     &     ( singxy + 9./4.*ftemp(3) + 9.*ftemp(4) ))

  324 CONTINUE

      else 
c      xcmax=1./(1.+ 4.)
      xcmax=(1./(1.+epsc4)) 
c      eps=1.
      eps=epsc
      if(xcmax.le.x) go to 321
      DO 325 I=1,NTERMS
      Y=0.5*(xcmax-X)*XI(I)+0.5*(xcmax+X)
      XY=X/Y
      CALL FETCH(XY,S,IPN,FTEMP)
      if(epsc.gt.1.) ftemp(3)=0.
      gluxy=ftemp(2)
      fcxy=ftemp(3)
      singxy=ftemp(5)
      ymul=y*(1+epsc4)
      cgvf2=facc*((c2gffnsh(y,eps)*(1-0.5*exp(1-1/eps**2))+
     .c2gffnsl(y,eps)*0.5*exp(1-1/eps**2))-c2gvfsub(ymul,eps)/xcmax)
      cqvf2=facc*((c2qffnsh(y,eps)*(1-0.5*exp(1-1/eps**2))+
     .c2qffnsl(y,eps)*0.5*exp(1-1/eps**2))-c2qvfsub(ymul,eps)/xcmax)
      ffc=ffc+0.5/xcmax*(xcmax-x)*wi(i)*al*al*
     .(fcxy*c2nn2a(ymul,3) + (fcxy-fc)*c2ns2b(ymul,3))
c$$$      ffc=ffc+0.5*(xcmax-x)*wi(i)*al**2.*(cgvf2*gluxy+cqvf2*singxy)
C--   G.W. 25/07/2007 Add charm and bottom contributions to singlet.
      ffc=ffc+0.5*(xcmax-x)*wi(i)*al**2.*(cgvf2*gluxy+cqvf2*
     &     ( singxy + 9./4.*ftemp(3) + 9.*ftemp(4) ))

  325 CONTINUE
      endif

      if(epsc.gt.1.) then 
      xcmax=1./(1.+epsc4)
      eps=epsc
      if(xcmax.le.x) go to 321
      DO 326 I=1,NTERMS
      Y=0.5*(xcmax-X)*XI(I)+0.5*(xcmax+X)
      XY=X/Y
      CALL FETCH(XY,S,IPN,FTEMP)
      gluxy=ftemp(2)
      singxy=ftemp(5)
      cgff23=facc*c2gffns3(y,eps)
      cqff23=facc*c2qffns3(y,eps)
c$$$      ffc=ffc+0.5*(xcmax-x)*wi(i)*al**3.*(cgff23*gluxy+cqff23*singxy)
C--   G.W. 25/07/2007 Add charm and bottom contributions to singlet.
      ffc=ffc+0.5*(xcmax-x)*wi(i)*al**3.*(cgff23*gluxy+cqff23*
     &     ( singxy + 9./4.*ftemp(3) + 9.*ftemp(4) ))

  326 CONTINUE

      else 
      xcmax=1./(1.+ 4.)
      eps=1.
      if(xcmax.le.x+0.00001) go to 321
      DO 327 I=1,NTERMS
      Y=0.5*(xcmax-X)*XI(I)+0.5*(xcmax+X)
      XY=X/Y
      CALL FETCH(XY,Schm,IPN,FTEMP)
      gluxy=ftemp(2)
      singxy=ftemp(5)
      cgff23=facc*c2gffns3(y,eps)
      cqff23=facc*c2qffns3(y,eps)
c$$$      ffc=ffc+0.5*(xcmax-x)*wi(i)*alchm**3*(cgff23*gluxy+cqff23*singxy)
C--   G.W. 25/07/2007 Add charm and bottom contributions to singlet.
      ffc=ffc+0.5*(xcmax-x)*wi(i)*alchm**3*(cgff23*gluxy+cqff23*
     &     ( singxy + 9./4.*ftemp(3) + 9.*ftemp(4) ))

  327 CONTINUE
      endif

      xcmax=1./(1.+epsc4)
      if(xcmax.le.x) go to 321
C--   G.W. 02/11/2007 Only add these terms if above threshold.
      xcmup=x/xcmax
      CALL FETCH(XCMUP,S,IPN,FTEMP)
      if(epsc.gt.1.) ftemp(3)=0. ! G.W. 05/11/2007
c$$$      fflc=fflc+ftemp(3)*(AL**2*(-0.0012)+al**3*(0.113+enf*0.006))
c$$$     .*1.25*(1/(1+4.*epsc)-0.2) ! N.B. -0.0012 not 0.012?
      fflc=fflc+ftemp(3)*(AL**2*CLNN2C(xcmup)+al**3*(0.113+enf*0.006))
     &     *1.25*(1/(1+4.*epsc)-0.2) ! G.W. 02/11/2007

 529  continue      
      if(epsc.gt.1.) then 
      xcmax=1./(1.+epsc4)
      eps=epsc
      if(xcmax.le.x) go to 321
      DO 524 I=1,NTERMS
      Y=0.5*(xcmax-X)*XI(I)+0.5*(xcmax+X)
      XY=X/Y
      CALL FETCH(XY,S,IPN,FTEMP)
      gluxy=ftemp(2)
      singxy=ftemp(5)
      cgffl=facc*(clgffnsl(y,eps)*(1-0.5*exp(1-eps**2))+
     .clgffnsh(y,eps)*0.5*exp(1-eps**2)) 
      cqffl=facc*(clqffnsl(y,eps)*(1-0.5*exp(1-eps**2))+
     .clqffnsh(y,eps)*0.5*exp(1-eps**2))
c$$$      fflc=fflc+0.5*(xcmax-x)*wi(i)*al**2.*(cgffl*gluxy+cqffl*singxy)
C--   G.W. 25/07/2007 Add charm and bottom contributions to singlet.
      fflc=fflc+0.5*(xcmax-x)*wi(i)*al**2.*(cgffl*gluxy+cqffl*
     &     ( singxy + 9./4.*ftemp(3) + 9.*ftemp(4) ))

  524 CONTINUE
      else 
      xcmax=(1./(1.+epsc4)) 
      eps=epsc
      if(xcmax.le.x) go to 321
      DO 525 I=1,NTERMS
      Y=0.5*(xcmax-X)*XI(I)+0.5*(xcmax+X)
      XY=X/Y
      CALL FETCH(XY,S,IPN,FTEMP)
      if(epsc.gt.1.) ftemp(3)=0.
      gluxy=ftemp(2)
      fcxy=ftemp(3)
      singxy=ftemp(5)
      ymul=y*(1+epsc4)
      Y1mul=1.-Ymul
      DL=log(Ymul)
      DL2=DL*DL
      DLM1=log(Y1mul)
      DLM2=DLM1*DLM1
      DLM3=DLM2*DLM1
      DLM4=DLM3*DLM1
c$$$      FNS2LQmul=128./9.d0*ymul*DLM1**2-46.50*ymul*DLM1-84.094*DL*DLM1
c$$$     x-37.338 +89.53*ymul
c$$$     X+33.82*ymul**2+ymul*DL*(32.90+18.41*DL)-128./9.d0*DL
c$$$     X+16./27.d0*enf*(6.*ymul*DLM1-12*ymul*DL-25.*ymul+6.)
      FNS2LQmul = CLNN2A(YMUL,INT(ENF)) ! G.W. 02/11/2007
      FS2LQmul=((15.94-5.212*ymul)*Y1mul*Y1mul*DLM1+(0.421+1.520*ymul)
     x*DL*DL+28.09*Y1mul*DL-(2.370/Ymul-19.27)*Y1mul**3)
      c3lg1mul=((144.*DLM4-47024./27.*DLM3+6319.*DLM2
     x+53160.*DLM1)*Y1mul
     X+72549*DL*DLM1+88238*DL*DL*DLM1+(3709.-33514.*Ymul
     x-9533.*Ymul*Ymul)*Y1mul
     X+66773.*Ymul*DL*DL -1117.*DL+45.37*DL*DL-5360./27.*DL**3
     X-(2044.7*Y1mul+409.506*DL)/Ymul)
      c3lg2mul=enf*((32./3.*DLM3-1216./9.*DLM2-592.3*DLM1
     X+1511.*ymul*DLM1)*Y1mul+311.3*DL*DLM1+14.24*DL*DL*DLM1
     X+(577.3-729.0*ymul)*Y1mul+30.78*ymul*DL**3+366.0*DL+1000./9*DL*DL
     X+160./9.*DL**3+88.5037*y1mul/ymul)
      c3lq1mul=512./27.*DLM4-177.4*DLM3+650.6*DLM2-2729.*DLM1
     x-2220.5-7884.*ymul+4168*ymul*Ymul-(844.7*DL+517.3*DLM1)*DL*DLM1
     x+(195.6*DLM1-125.3)*y1mul*DLM3+208.3*ymul*DL**3-1355.7*DL
     x-7456./27.*DL*DL-1280./81*DL**3
      c3lq2mul=enf*(1024./81.*DLM3-112.35*DLM2+344.1*DLM1+408.4
     x-9.345*ymul-919.3*Ymul*Ymul+(239.7+20.63*DLM1)*y1mul*DLM2+(887.3
     x+294.5*DL-59.14*DLM1)*DL*DLM1-1792./81.*ymul*DL**3+200.73*DL
     x+64./3.*DL*DL )
      c3lq3mul=enf*enf*(3.*ymul*DLM2+(6.-25*ymul)*DLM1-19.+(317./6.
     x-12*1.645)*ymul-6.*ymul*DL*DLM1+6.*ymul*(ymul+ymul**2/4.
     x+ymul**3/9.+ymul**4/16.
     x+ymul**5/25.+ymul**6/36.)+9.*ymul*DL*DL-(6.-50.*ymul)*DL)*64./81.
      c3lq4mul=((1568./27.*DLM3-3968./9.*DLM2+5124*DLM1)*y1mul*y1mul
     X+(2184.*DL+6059.*y1mul)*DL*DLM1-(795.6+1036.*ymul)*y1mul*y1mul
     x-143.6*y1mul*DL+2848./9.*DL*DL-1600./27.*DL**3
     x-(885.53*y1mul+182.0*DL)*y1mul/ymul)
      c3lq5mul=enf*((-32./9.*DLM2+29.52*DLM1)*y1mul*y1mul+(35.18*DL
     X+73.06*y1mul)*DL*DLM1-35.24*ymul*DL*DL
     x-(14.16-69.84*ymul)*y1mul*y1mul
     X-69.41*y1mul*DL-128./9.*DL*DL+40.239*y1mul*y1mul/ymul)
      c3lgmul=c3lg1mul+c3lg2mul
      c3lqnsmul=c3lq1mul+c3lq2mul+c3lq3mul
      c3lqpsmul=c3lq4mul+c3lq5mul
      cgvfl=facc*((clgffnsh(y,eps)*(1-0.5*exp(1-1/eps**2))+
     .clgffnsl(y,eps)*0.5*exp(1-1/eps**2))-clgvfsub(ymul,eps)/xcmax)
      cqvfl=facc*((clqffnsh(y,eps)*(1-0.5*exp(1-1/eps**2))+
     .clqffnsl(y,eps)*0.5*exp(1-1/eps**2)))
      fflc=fflc+0.5/xcmax*(xcmax-x)*wi(i)*al**2.*(fcxy*FNS2LQmul+
     xfcxy*FS2LQmul)*1.25*(1/(1+4.*eps)-0.2)
c$$$      fflc=fflc+0.5*(xcmax-x)*wi(i)*al**2.*(cgvfl*gluxy+cqvfl*singxy)
c$$$      fflc=fflc+0.5*(xcmax-x)*wi(i)*al*al*al*facc*
c$$$     .(c3lqnsmul*fcxy+c3lqpsmul*singxy
c$$$     .+c3lgmul*gluxy)*1.25*(1/(1+4.*eps)-0.2)
C--   G.W. 25/07/2007 Add charm and bottom contributions to singlet.
      fflc=fflc+0.5*(xcmax-x)*wi(i)*al**2.*(cgvfl*gluxy+cqvfl*
     &     ( singxy + 9./4.*ftemp(3) + 9.*ftemp(4) ))
      fflc=fflc+0.5*(xcmax-x)*wi(i)*al*al*al*facc*
     &     (c3lqnsmul*fcxy+c3lqpsmul*
     &     (singxy+9./4.*ftemp(3)+9.*ftemp(4))
     &     +c3lgmul*gluxy)*1.25*(1/(1+4.*eps)-0.2)

  525 CONTINUE
      endif


      if(epsc.gt.1.) then 
      xcmax=1./(1.+epsc4)
      eps=epsc
      if(xcmax.le.x) go to 321
      DO 526 I=1,NTERMS
      Y=0.5*(xcmax-X)*XI(I)+0.5*(xcmax+X)
      XY=X/Y
      CALL FETCH(XY,S,IPN,FTEMP)
      gluxy=ftemp(2)
      singxy=ftemp(5)
      cgffl3=facc*clgffns3(y,eps)
      cqffl3=facc*clqffns3(y,eps)
c$$$      fflc=fflc+0.5*(xcmax-x)*wi(i)*al**3.*(cgffl3*gluxy+cqffl3*singxy)
C--   G.W. 25/07/2007 Add charm and bottom contributions to singlet.
      fflc=fflc+0.5*(xcmax-x)*wi(i)*al**3.*(cgffl3*gluxy+cqffl3*
     &     ( singxy + 9./4.*ftemp(3) + 9.*ftemp(4) ))

  526 CONTINUE

      else 
      xcmax=1./(1.+ 4.)
      eps=1.
      if(xcmax.le.x+0.00001) go to 321
      DO 527 I=1,NTERMS
      Y=0.5*(xcmax-X)*XI(I)+0.5*(xcmax+X)
      XY=X/Y
      CALL FETCH(XY,Schm,IPN,FTEMP)
      gluxy=ftemp(2)
      singxy=ftemp(5)
      cgffl3=facc*clgffns3(y,eps)
      cqffl3=facc*clqffns3(y,eps)
      clgcvagg=facc*clgconagg(y,eps)
c$$$      fflc=fflc+0.5*(xcmax-x)*wi(i)*alchm**3*((-clgcvagg+cgffl3)*gluxy
c$$$     .+cqffl3*singxy)*1.25*(1-1./(1.+4.*epsc))
C--   G.W. 25/07/2007 Add charm and bottom contributions to singlet.
      fflc=fflc+0.5*(xcmax-x)*wi(i)*alchm**3*((-clgcvagg+cgffl3)*gluxy
     &     +cqffl3*(singxy+9./4.*ftemp(3)+9.*ftemp(4)))*
     &     1.25*(1-1./(1.+4.*epsc))

  527 CONTINUE
      endif

  321 CONTINUE

      if(ffc.lt.0.) ffc=0.

      xbmax=1./(1.+epsb4)
      if(xbmax.le.x+0.00001) go to 421
      xbmup=x/xbmax
      CALL FETCH(Xbmup,S,IPN,FTEMP)
      if(epsb.gt.1.) ftemp(4)=0.
      fb=ftemp(4)
      AL1b=log(1.-xbmup)
      ffb=Fb+Fb*AL*CF*(-9.-2.*PI2/3.+AL1b*(-3.+2.*AL1b))
      ffb=ffb+fb*al*al*c2nn2c(xbmup,3)

      DO 423 I=1,NTERMS
      Y=0.5*(xbmax-X)*XI(I)+0.5*(xbmax+X)
      XY=X/Y
      CALL FETCH(XY,S,IPN,FTEMP)
      if(epsb.gt.1.) ftemp(4)=0.
      gluxy=ftemp(2)
      fbxy=ftemp(4)
      ybmup=y/xbmax
      if(ybmup.gt.0.99999d0) ybmup=0.99999d0
      c0b=1.
      p0qg=ybmup**2+(1.-ybmup)**2
      C22b=CF*(6.+4.*Ybmup-2.*(1.+Ybmup*Ybmup)/(1.-Ybmup)*
     .log(Ybmup)-2.*(1.+Ybmup)*log(1.-ybmup)
     2-IF3*2.*(1.+Ybmup))
      C23b=CF*(-3.+4.*log(1.-ybmup))/(1.-Ybmup)
      if(epsb.gt.1.d0) c0b=0.d0
      cg21b=2.*facb*cheavy(1,y,epsb)
      cg22b=2.*facb*c0b*p0qg*log(1./epsb)
      ffb=ffb+0.5/xbmax*(xbmax-X)*WI(I)*AL*(C22b*fbxy+C23b*(fbxy-fb))      
      ffb=ffb+0.5*(xbmax-x)*wi(i)*al*(cg21b-cg22b/xbmax)*gluxy
      clg2b=2.*facb*cheavy(7,y,epsb)
      f1lq=cheavy(8,ybmup,epsb)
      fflb=fflb+0.5/xbmax*(xbmax-x)*wi(i)*al*f1lq*fbxy
      fflb=fflb+0.5*(xbmax-x)*wi(i)*al*clg2b*gluxy


  423 CONTINUE

      if(epsb.gt.1.) then 
      xbmax=1./(1.+epsb4)
      eps=epsb
      if(xbmax.le.x) go to 421
      DO 424 I=1,NTERMS
      Y=0.5*(xbmax-X)*XI(I)+0.5*(xbmax+X)
      XY=X/Y
      CALL FETCH(XY,S,IPN,FTEMP)
      gluxy=ftemp(2)
c$$$      singxy=ftemp(5)+9./8.*ftemp(3) ! Why 9/8 and not 9/4?
      singxy=ftemp(5)           ! G.W. 27/07/2007
      cgff2=facb*(c2gffnsl(y,eps)*(1-0.5*exp(1-eps**2))+
     .c2gffnsh(y,eps)*0.5*exp(1-eps**2))
      cqff2=facb*(c2qffnsl(y,eps)*(1-0.5*exp(1-eps**2))+
     .c2qffnsh(y,eps)*0.5*exp(1-eps**2))
c$$$      ffb=ffb+0.5*(xbmax-x)*wi(i)*al**2.*(cgff2*gluxy+cqff2*singxy)
C--   G.W. 25/07/2007 Add charm and bottom contributions to singlet.
      ffb=ffb+0.5*(xbmax-x)*wi(i)*al**2.*(cgff2*gluxy+cqff2*
     &     ( singxy + 9./4.*ftemp(3) + 9.*ftemp(4) ))

  424 CONTINUE

      else 
c      xbmax=1./(1.+ 4.)
      xbmax=(1./(1.+epsb4)) 
c      eps=1.
      eps=epsb
      if(xbmax.le.x) go to 421
      DO 425 I=1,NTERMS
      Y=0.5*(xbmax-X)*XI(I)+0.5*(xbmax+X)
      XY=X/Y
      CALL FETCH(XY,S,IPN,FTEMP)
      if(epsb.gt.1.)  ftemp(4)=0.
      gluxy=ftemp(2)
      fbxy=ftemp(4)
c$$$      singxy=ftemp(5)+9./8.*ftemp(3) ! Why 9/8 and not 9/4?
      singxy=ftemp(5)           ! G.W. 27/07/2007
      ymul=y*(1+epsb4)
      cgvf2=facb*((c2gffnsh(y,eps)*(1-0.5*exp(1-1/eps**2))+
     .c2gffnsl(y,eps)*0.5*exp(1-1/eps**2))-c2gvfsub(ymul,eps)/xbmax)
      cqvf2=facb*((c2qffnsh(y,eps)*(1-0.5*exp(1-1/eps**2))+
     .c2qffnsl(y,eps)*0.5*exp(1-1/eps**2))-c2qvfsub(ymul,eps)/xbmax)
      ffb=ffb+0.5/xbmax*(xbmax-x)*wi(i)*al*al*
     .(fbxy*c2nn2a(ymul,3) + (fbxy-fb)*c2ns2b(ymul,3))
c$$$      ffb=ffb+0.5*(xbmax-x)*wi(i)*al**2.*(cgvf2*gluxy+cqvf2*singxy)
C--   G.W. 25/07/2007 Add charm and bottom contributions to singlet.
      ffb=ffb+0.5*(xbmax-x)*wi(i)*al**2.*(cgvf2*gluxy+cqvf2*
     &     ( singxy + 9./4.*ftemp(3) + 9.*ftemp(4) ))

  425 CONTINUE
      endif

      if(epsb.gt.1.) then 
      xbmax=1./(1.+epsb4)
      eps=epsb
      if(xbmax.le.x) go to 421
      DO 426 I=1,NTERMS
      Y=0.5*(xbmax-X)*XI(I)+0.5*(xbmax+X)
      XY=X/Y
      CALL FETCH(XY,S,IPN,FTEMP)
      gluxy=ftemp(2)
c$$$      singxy=ftemp(5)+9./8.*ftemp(3)! Why 9/8 and not 9/4?
      singxy=ftemp(5)           ! G.W. 27/07/2007
      cgff23=facb*c2gffns3(y,eps)
      cqff23=facb*c2qffns3(y,eps)
c$$$      ffb=ffb+0.5*(xbmax-x)*wi(i)*al**3.*(cgff23*gluxy+cqff23*singxy)
C--   G.W. 25/07/2007 Add charm and bottom contributions to singlet.
      ffb=ffb+0.5*(xbmax-x)*wi(i)*al**3.*(cgff23*gluxy+cqff23*
     &     ( singxy + 9./4.*ftemp(3) + 9.*ftemp(4) ))

  426 CONTINUE

      else 
      xbmax=1./(1.+ 4.)
      eps=1.
      if(xbmax.le.x+0.00001) go to 421
      DO 427 I=1,NTERMS
      Y=0.5*(xbmax-X)*XI(I)+0.5*(xbmax+X)
      XY=X/Y
      CALL FETCH(XY,Sbot,IPN,FTEMP)
      gluxy=ftemp(2)
c$$$      singxy=ftemp(5)+9./8.*ftemp(3) ! Why 9/8 and not 9/4?
      singxy=ftemp(5)           ! G.W. 27/07/2007
      cgff23=facb*c2gffns3(y,eps)
      cqff23=facb*c2qffns3(y,eps)
c$$$      ffb=ffb+0.5*(xbmax-x)*wi(i)*albot**3*(cgff23*gluxy+cqff23*singxy)
C--   G.W. 25/07/2007 Add charm and bottom contributions to singlet.
      ffb=ffb+0.5*(xbmax-x)*wi(i)*albot**3*(cgff23*gluxy+cqff23*
     &     ( singxy + 9./4.*ftemp(3) + 9.*ftemp(4) ))

  427 CONTINUE
      endif

C--   G.W. 02/11/2007 This contribution was missing for b, but included for c.
      xbmax=1./(1.+epsb4)
      if(xbmax.le.x) go to 421
      xbmup=x/xbmax
      CALL FETCH(XBMUP,S,IPN,FTEMP)
      if(epsb.gt.1.) ftemp(4)=0. ! G.W. 05/11/2007
      fflb=fflb+ftemp(4)*(AL**2*CLNN2C(xbmup)+al**3*(0.113+enf*0.006))
     &     *1.25*(1/(1+4.*epsb)-0.2) ! G.W. 02/11/2007

 629  continue
      if(epsb.gt.1.) then 
      xbmax=1./(1.+epsb4)
      eps=epsb
      if(xbmax.le.x) go to 421
      DO 624 I=1,NTERMS
      Y=0.5*(xbmax-X)*XI(I)+0.5*(xbmax+X)
      XY=X/Y
      CALL FETCH(XY,S,IPN,FTEMP)
      gluxy=ftemp(2)
      singxy=ftemp(5)
      cgffl=facb*(clgffnsl(y,eps)*(1-0.5*exp(1-eps**2))+
     .clgffnsh(y,eps)*0.5*exp(1-eps**2)) 
      cqffl=facb*(clqffnsl(y,eps)*(1-0.5*exp(1-eps**2))+
     .clqffnsh(y,eps)*0.5*exp(1-eps**2))
c$$$      fflb=fflb+0.5*(xbmax-x)*wi(i)*al**2.*(cgffl*gluxy+cqffl*singxy)
C--   G.W. 25/07/2007 Add charm and bottom contributions to singlet.
      fflb=fflb+0.5*(xbmax-x)*wi(i)*al**2.*(cgffl*gluxy+cqffl*
     &     ( singxy + 9./4.*ftemp(3) + 9.*ftemp(4) ))

  624 CONTINUE
      else 
      xbmax=(1./(1.+epsb4)) 
      eps=epsb
      if(xbmax.le.x) go to 421
      DO 625 I=1,NTERMS
      Y=0.5*(xbmax-X)*XI(I)+0.5*(xbmax+X)
      XY=X/Y
      CALL FETCH(XY,S,IPN,FTEMP)
      if(epsb.gt.1.) ftemp(3)=0.
      gluxy=ftemp(2)
      fbxy=ftemp(4)
      singxy=ftemp(5)
      ymul=y*(1+epsb4)
      Y1mul=1.-Ymul
      DL=log(Ymul)
      DL2=DL*DL
      DLM1=log(Y1mul)
      DLM2=DLM1*DLM1
      DLM3=DLM2*DLM1
      DLM4=DLM3*DLM1
c$$$      FNS2LQmul=128./9.d0*ymul*DLM1**2-46.50*ymul*DLM1-84.094*DL*DLM1
c$$$     x-37.338 +89.53*ymul
c$$$     X+33.82*ymul**2+ymul*DL*(32.90+18.41*DL)-128./9.d0*DL
c$$$     X+16./27.d0*enf*(6.*ymul*DLM1-12*ymul*DL-25.*ymul+6.)
      FNS2LQmul = CLNN2A(YMUL,INT(ENF)) ! G.W. 02/11/2007
      FS2LQmul=((15.94-5.212*ymul)*Y1mul*Y1mul*DLM1+(0.421+1.520*ymul)
     x*DL*DL+28.09*Y1mul*DL-(2.370/Ymul-19.27)*Y1mul**3)
      c3lg1mul=((144.*DLM4-47024./27.*DLM3+6319.*DLM2
     x+53160.*DLM1)*Y1mul
     X+72549*DL*DLM1+88238*DL*DL*DLM1+(3709.-33514.*Ymul
     x-9533.*Ymul*Ymul)*Y1mul
     X+66773.*Ymul*DL*DL -1117.*DL+45.37*DL*DL-5360./27.*DL**3
     X-(2044.7*Y1mul+409.506*DL)/Ymul)
      c3lg2mul=enf*((32./3.*DLM3-1216./9.*DLM2-592.3*DLM1
     X+1511.*ymul*DLM1)*Y1mul+311.3*DL*DLM1+14.24*DL*DL*DLM1
     X+(577.3-729.0*ymul)*Y1mul+30.78*ymul*DL**3+366.0*DL+1000./9*DL*DL
     X+160./9.*DL**3+88.5037*y1mul/ymul)
      c3lq1mul=512./27.*DLM4-177.4*DLM3+650.6*DLM2-2729.*DLM1
     x-2220.5-7884.*ymul+4168*ymul*Ymul-(844.7*DL+517.3*DLM1)*DL*DLM1
     x+(195.6*DLM1-125.3)*y1mul*DLM3+208.3*ymul*DL**3-1355.7*DL
     x-7456./27.*DL*DL-1280./81*DL**3
      c3lq2mul=enf*(1024./81.*DLM3-112.35*DLM2+344.1*DLM1+408.4
     x-9.345*ymul-919.3*Ymul*Ymul+(239.7+20.63*DLM1)*y1mul*DLM2+(887.3
     x+294.5*DL-59.14*DLM1)*DL*DLM1-1792./81.*ymul*DL**3+200.73*DL
     x+64./3.*DL*DL )
      c3lq3mul=enf*enf*(3.*ymul*DLM2+(6.-25*ymul)*DLM1-19.+(317./6.
     x-12*1.645)*ymul-6.*ymul*DL*DLM1+6.*ymul*(ymul+ymul**2/4.
     x+ymul**3/9.+ymul**4/16.
     x+ymul**5/25.+ymul**6/36.)+9.*ymul*DL*DL-(6.-50.*ymul)*DL)*64./81.
      c3lq4mul=((1568./27.*DLM3-3968./9.*DLM2+5124*DLM1)*y1mul*y1mul
     X+(2184.*DL+6059.*y1mul)*DL*DLM1-(795.6+1036.*ymul)*y1mul*y1mul
     x-143.6*y1mul*DL+2848./9.*DL*DL-1600./27.*DL**3
     x-(885.53*y1mul+182.0*DL)*y1mul/ymul)
      c3lq5mul=enf*((-32./9.*DLM2+29.52*DLM1)*y1mul*y1mul+(35.18*DL
     X+73.06*y1mul)*DL*DLM1-35.24*ymul*DL*DL
     x-(14.16-69.84*ymul)*y1mul*y1mul
     X-69.41*y1mul*DL-128./9.*DL*DL+40.239*y1mul*y1mul/ymul)
      c3lgmul=c3lg1mul+c3lg2mul
      c3lqnsmul=c3lq1mul+c3lq2mul+c3lq3mul
      c3lqpsmul=c3lq4mul+c3lq5mul
      cgvfl=facb*((clgffnsh(y,eps)*(1-0.5*exp(1-1/eps**2))+
     .clgffnsl(y,eps)*0.5*exp(1-1/eps**2))-clgvfsub(ymul,eps)/xcmax)
      cqvfl=facb*((clqffnsh(y,eps)*(1-0.5*exp(1-1/eps**2))+
     .clqffnsl(y,eps)*0.5*exp(1-1/eps**2)))
      fflb=fflb+0.5/xbmax*(xbmax-x)*wi(i)*al**2.*(fbxy*FNS2LQmul+
     xfbxy*FS2LQmul)*1.25*(1/(1+4.*eps)-0.2)
c$$$      fflb=fflb+0.5*(xbmax-x)*wi(i)*al**2.*(cgvfl*gluxy+cqvfl*singxy)
c$$$      fflb=fflb+0.5*(xbmax-x)*wi(i)*al*al*al*facb*
c$$$     .(c3lqnsmul*fbxy+c3lqpsmul*singxy
c$$$     .+c3lgmul*gluxy)*1.25*(1/(1+4.*eps)-0.2)
C--   G.W. 25/07/2007 Add charm and bottom contributions to singlet.
      fflb=fflb+0.5*(xbmax-x)*wi(i)*al**2.*(cgvfl*gluxy+cqvfl*
     &     ( singxy + 9./4.*ftemp(3) + 9.*ftemp(4) ))
      fflb=fflb+0.5*(xbmax-x)*wi(i)*al*al*al*facb*
     &     (c3lqnsmul*fbxy+c3lqpsmul*
     &     (singxy+9./4.*ftemp(3)+9.*ftemp(4))
     &     +c3lgmul*gluxy)*1.25*(1/(1+4.*eps)-0.2)

  625 CONTINUE
      endif


      if(epsb.gt.1.) then 
      xbmax=1./(1.+epsb4)
      eps=epsb
      if(xbmax.le.x) go to 421
      DO 626 I=1,NTERMS
      Y=0.5*(xbmax-X)*XI(I)+0.5*(xbmax+X)
      XY=X/Y
      CALL FETCH(XY,S,IPN,FTEMP)
      gluxy=ftemp(2)
      singxy=ftemp(5)
      cgffl3=facb*clgffns3(y,eps)
      cqffl3=facb*clqffns3(y,eps)
c$$$      fflb=fflb+0.5*(xbmax-x)*wi(i)*al**3.*(cgffl3*gluxy+cqffl3*singxy)
C--   G.W. 25/07/2007 Add charm and bottom contributions to singlet.
      fflb=fflb+0.5*(xbmax-x)*wi(i)*al**3.*(cgffl3*gluxy+cqffl3*
     &     ( singxy + 9./4.*ftemp(3) + 9.*ftemp(4) ))

  626 CONTINUE

      else 
      xbmax=1./(1.+ 4.)
      eps=1.
      if(xbmax.le.x+0.00001) go to 421
      DO 627 I=1,NTERMS
      Y=0.5*(xbmax-X)*XI(I)+0.5*(xbmax+X)
      XY=X/Y
      CALL FETCH(XY,Sbot,IPN,FTEMP)
      gluxy=ftemp(2)
      singxy=ftemp(5)
      cgffl3=facb*clgffns3(y,eps)
      cqffl3=facb*clqffns3(y,eps)
      clgcvagg=facb*clgconagg(y,eps)
c$$$      fflb=fflb+0.5*(xbmax-x)*wi(i)*albot**3*((-clgcvagg+cgffl3)*gluxy
c$$$     .+cqffl3*singxy)*1.25*(1-1./(1.+4.*epsb))
C--   G.W. 25/07/2007 Add charm and bottom contributions to singlet.
      fflb=fflb+0.5*(xbmax-x)*wi(i)*albot**3*((-clgcvagg+cgffl3)*gluxy
     &     +cqffl3*(singxy+9./4.*ftemp(3)+9.*ftemp(4)))*
     &     1.25*(1-1./(1.+4.*epsb))

  627 CONTINUE
      endif

  421 CONTINUE

      if(ffb.lt.0.) ffb=0.

      F2 = ffx+ffc+ffb
      F2C = ffc
      F2B = ffb
      FL = fflx+fflc+fflb
      FLC = fflc
      FLB = fflb

      RETURN
      END
C--   End of MSTWNCnnlo subroutine.


C--   TO DO: charged current version.
C--   Only W+,p and (p + n)/2,(W+ + W-)/2 currently in fitting code.
C--   Don't include nuclear corrections: allow user to add separately.
C--   Input variables for MSTWCC:
C--     x = Bjorken-x value.
C--     q = sqrt(Q^2) in GeV.
C--     ipn = 1, 2 for p or n.
C--     iwpm = 1 for W+, 2 for W-.
C--     iord = 0, 1, 2 for LO, NLO, NNLO (pass in COMMON block).
C--   Output variables: f2,f2c,fl,flc,xf3,xf3c.
C      SUBROUTINE MSTWCC(x,q,ipn,f2,f2c,fl,flc,xf3,xf3c)


C--   Rewrite FETCH to return PDFs from grids.
      SUBROUTINE FETCH(X,S,IPN,FTEMP)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) ! G.W. 15/02/2007
      DIMENSION FTEMP(5)
      COMMON/GRPTHY/FLAVOR
      COMMON/DYLAMB/xlam,S0 ! G.W. 17/07/2007

C--   Structure functions.
      CALL GINTRP(7,X,S,IPN,SING)
      CALL GINTRP(3,X,S,IPN,CPLUS)
      CALL GINTRP(6,X,S,IPN,BPLUS)
      ftemp(5)=sing
   20 CONTINUE
      CALL GINTRP(8,X,S,IPN,FTEMP(2))
      CALL GINTRP(4,X,S,IPN,SPLUS)
      CALL GINTRP(5,X,S,IPN,UPLUS)
      SPLUS=-SPLUS+SING/FLAVOR
      UPLUS=UPLUS+SING/FLAVOR
      DPLUS=SING-UPLUS-SPLUS
C--   Proton structure functions.
      FTEMP(1)=(4.*(UPLUS)+DPLUS+SPLUS)/9.
      ftemp(3)=4.*cplus/9.
      ftemp(4)=bplus/9.

      RETURN
      END


      SUBROUTINE GINTRP(I,X,S,IPN,ANS)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON/DYLAMB/XLAM,S0
      COMMON/GRPTHY/FLAVOR
      CHARACTER prefix*50,prefix1*50,cl*4
      INTEGER ISET
      COMMON/mstwfiles/iset,prefix,cl

      IF (iset.EQ.0) THEN
         prefix1 = prefix
      ELSE
         prefix1 = prefix(1:len_trim(prefix))//'.'//cl
      END IF

      XMUF = sqrt(XLAM**2*exp(S0*exp(S)))
      CALL GetAllPDFs(prefix1,iset,x,xmuf,upv,dnv,usea,dsea,str,sbar,
     &     chm,cbar,bot,bbar,glu,phot)

      IF (IPN.EQ.2) THEN        ! neutron
C--   Swap upv and dnv.
         upvOrig = upv
         upv = dnv
         dnv = upvOrig
C--   Swap usea and dsea.
         useaOrig = usea
         usea = dsea
         dsea = useaOrig
      ELSE IF (IPN.NE.1) THEN
         WRITE(6,*) "Error in GINTRP, IPN = ",IPN
         STOP
      END IF

      sing = upv + 2.D0*usea + dnv + 2.D0*dsea + str + sbar
      IF (I.EQ.1) THEN
         ANS = upv+dnv
      ELSE IF (I.EQ.2) THEN
         ANS = dnv - (upv+dnv)/2.D0
      ELSE IF (I.EQ.3) THEN
         ANS = chm+cbar
      ELSE IF (I.EQ.4) THEN
         ANS = -(str+sbar) + sing/flavor
      ELSE IF (I.EQ.5) THEN
         ANS = upv + 2.D0*usea - sing/flavor
      ELSE IF (I.EQ.6) THEN
         ANS = bot + bbar
      ELSE IF (I.EQ.7) THEN
         ANS = sing
      ELSE IF (I.EQ.8) THEN
         ANS = glu
      ELSE IF (I.EQ.9) THEN
         ANS = str-sbar - (upv+dnv)/2.D0
      ELSE IF (I.EQ.10) THEN
         ANS = chm-cbar - (upv+dnv)/2.D0
      ELSE IF (I.EQ.11) THEN
         ANS = bot-bbar - (upv+dnv)/2.D0
      ELSE
         WRITE(6,*) "Error in GINTRP: I = ",I
         STOP
      END IF
         
      RETURN
      END


      FUNCTION ALPHA(T)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON/DYLAMB/XLAM,S0
C--   G.W. 15/06/2007 Use new routine from PEGASUS.
      QS = XLAM**2*EXP(T)
      ALPHA = ALPHAS(sqrt(QS))
      RETURN
      END


      SUBROUTINE WATE96
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) ! G.W. 15/02/2007
C*******************************************************************
C*****							       *****
C***** THE X(I)	AND W(I) ARE THE DIRECT	OUTPUT FROM A PROGRAM  *****
C***** USING NAG ROUTINE D01BCF	TO CALCULATE THE	       *****
C***** GAUSS-LEGENDRE WEIGHTS FOR 96 POINT INTEGRATION.	       *****
C***** THEY AGREE TO TYPICALLY 14 DECIMAL PLACES WITH THE      *****
C***** TABLE IN	ABRAMOWITZ & STEGUN, PAGE 919.		       *****
C*****							       *****
C***** ---->   PETER HARRIMAN, APRIL 3RD 1990.		       *****
C*****							       *****
C*******************************************************************
      DIMENSION	X(48),W(48)
c$$$      COMMON/GAUS96/XI(96),WI(96),NTERMS,XX(97)
      COMMON/GAUS96/XI(96),WI(96),XX(97),NTERMS ! G.W. 15/02/2007
      NTERMS=96

      X( 1)=   0.01627674484960183561
      X( 2)=   0.04881298513604856015
      X( 3)=   0.08129749546442434360
      X( 4)=   0.11369585011066471632
      X( 5)=   0.14597371465489567682
      X( 6)=   0.17809688236761733390
      X( 7)=   0.21003131046056591064
      X( 8)=   0.24174315616383866556
      X( 9)=   0.27319881259104774468
      X(10)=   0.30436494435449495954
      X(11)=   0.33520852289262397655
      X(12)=   0.36569686147231213885
      X(13)=   0.39579764982890709712
      X(14)=   0.42547898840729897474
      X(15)=   0.45470942216774136446
      X(16)=   0.48345797392059470382
      X(17)=   0.51169417715466604391
      X(18)=   0.53938810832435567233
      X(19)=   0.56651041856139533470
      X(20)=   0.59303236477757022282
      X(21)=   0.61892584012546672523
      X(22)=   0.64416340378496526886
      X(23)=   0.66871831004391424358
      X(24)=   0.69256453664216964528
      X(25)=   0.71567681234896561582
      X(26)=   0.73803064374439816819
      X(27)=   0.75960234117664555964
      X(28)=   0.78036904386743123629
      X(29)=   0.80030874413913884180
      X(30)=   0.81940031073792957139
      X(31)=   0.83762351122818502758
      X(32)=   0.85495903343459936363
      X(33)=   0.87138850590929436968
      X(34)=   0.88689451740241818933
      X(35)=   0.90146063531585023110
      X(36)=   0.91507142312089592706
      X(37)=   0.92771245672230655266
      X(38)=   0.93937033975275308073
      X(39)=   0.95003271778443564022
      X(40)=   0.95968829144874048809
      X(41)=   0.96832682846326217918
      X(42)=   0.97593917458513455843
      X(43)=   0.98251726356301274934
      X(44)=   0.98805412632962202890
      X(45)=   0.99254390032376081654
      X(46)=   0.99598184298720747465
      X(47)=   0.99836437586317963722
      X(48)=   0.99968950388322870559
      W( 1)=   0.03255061449236316962
      W( 2)=   0.03251611871386883307
      W( 3)=   0.03244716371406427668
      W( 4)=   0.03234382256857594104
      W( 5)=   0.03220620479403026124
      W( 6)=   0.03203445623199267876
      W( 7)=   0.03182875889441101874
      W( 8)=   0.03158933077072719007
      W( 9)=   0.03131642559686137819
      W(10)=   0.03101033258631386231
      W(11)=   0.03067137612366917839
      W(12)=   0.03029991542082762553
      W(13)=   0.02989634413632842385
      W(14)=   0.02946108995816795100
      W(15)=   0.02899461415055528410
      W(16)=   0.02849741106508543861
      W(17)=   0.02797000761684838950
      W(18)=   0.02741296272602931385
      W(19)=   0.02682686672559184485
      W(20)=   0.02621234073567250055
      W(21)=   0.02557003600534944960
      W(22)=   0.02490063322248370695
      W(23)=   0.02420484179236479915
      W(24)=   0.02348339908592633665
      W(25)=   0.02273706965832950717
      W(26)=   0.02196664443874448477
      W(27)=   0.02117293989219144572
      W(28)=   0.02035679715433347898
      W(29)=   0.01951908114014518992
      W(30)=   0.01866067962741165898
      W(31)=   0.01778250231604547316
      W(32)=   0.01688547986424539715
      W(33)=   0.01597056290256253144
      W(34)=   0.01503872102699521608
      W(35)=   0.01409094177231515264
      W(36)=   0.01312822956696188190
      W(37)=   0.01215160467108866759
      W(38)=   0.01116210209983888144
      W(39)=   0.01016077053500880978
      W(40)=   0.00914867123078384552
      W(41)=   0.00812687692569928101
      W(42)=   0.00709647079115442616
      W(43)=   0.00605854550423662775
      W(44)=   0.00501420274292825661
      W(45)=   0.00396455433844564804
      W(46)=   0.00291073181793626202
      W(47)=   0.00185396078894924657
      W(48)=   0.00079679206555731759
      DO 1 I=1,48
      XI(I)=-X(49-I)
      WI(I)=W(49-I)
      XI(I+48)=X(I)
      WI(I+48)=W(I)
    1 CONTINUE
      DO 2 I=1,96
    2 XX(I)=0.5*(XI(I)+1.)
      XX(97)=1.0
c$$$      EXPON=1.5D0
      EXPON=2.D0                ! G.W. 04/07/2007
c$$$      EXPON=3.D0
      DO 3 I=1,96
      YI=2.*(0.5*(1.+XI(I)))**EXPON-1.
      WI(I)=WI(I)/(1.+XI(I))*(1.+YI)*EXPON
      XI(I)=YI
      XX(I)=0.5*(1.+YI)
    3 CONTINUE
      RETURN
      END


C----------------------------------------------------------------------


      function cheavy(i,z,eps)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z) ! G.W. 15/02/2007

c     this function returns the values of C_g(z,Q^2) 
c     and the deriv. wrt log Q^2. Here eps=m^2/Q^2.
c     If i=1,3  C_g for F2.  If i=2,4 deriv of C_g for F2 
c     i=1,2 refer to photon current
c     i=3,4 refer to W current for F2
c     i=5,6 refer to W current for xF3
c     i=7,8 refer to FL 

C--   G.W. 22/08/2006 Added error message.
      if(i.gt.8) then
         write(6,*) "Error in cheavy: i = ",i
         stop
      end if
      z1=1.-z
      z2=z*z
      z3=z2*z
      zr=z1/z
      eps1=1.+eps
      eps2=eps*eps
      beta2=1.-4.*eps*z/z1
      if(i.gt.2) beta2=1.-eps*z/z1
      if(i.gt.6) beta2=1.-4.*eps*z/z1
c      if(beta2.lt.0.) go to 10
      beta=sqrt(beta2)
      a=z2+z1*z1
      b=4.*z*(1.-3.*z)
      c=-8.*z2
      aa=8.*z*z1-1.
      bb=-4.*z*z1
      arg=(1.+beta)/(1.-beta)
      fac=log(arg)
      alam=1./(1.+eps)
      zp=z/alam
      alz=alam*z
      alz1=1.-alz
      al1=1.-alam
      al2=alam*alam
      al3=al2*alam
      al4=al2*al2
      z1=1.-z
      z2=z*z
      zp1=1.-zp
      cf=4./3.
      go to (1,2,3,4,5,6,7,8) i
    1 cheavy=(a+b*eps+c*eps2)*fac+(aa+bb*eps)*beta
      return
    2 cheavy=(-b*eps-2.*c*eps2)*fac+(a+b*eps+c*eps2)/beta
     .      +(-bb*eps)*beta +(aa+bb*eps)*2.*z*eps/z1/beta
      return
    3 pqg=0.5*(zp*zp+zp1*zp1)
      arg0=alam*zp1*zp1/z1/z
      c0=pqg*log(arg0)
      arg=alam*z1/al1/z
      bigl=log(arg)
      f1=c0+(8.-18.*al1+12.*al1*al1)*zp*zp1+(al1/z1-1.)
      f3=al1*z*bigl*6.*(1.-2.*z)
      f2=pqg*( bigl-log(alam))
      cheavy=2.*(f1+f3+f2)/alam
      return
    4 arg=alam/al1/z2
      argbig=alam*z1/al1/z
      bigl=log(argbig)
      f1=(1.-2.*zp+2.*zp*zp)/al2
      f2=(4.*z*alam-6.*z2-al2)/al4*log(alam/z1/z)
      f3=4.*z/al4*(3.*z-2.*(1.+3.*z)*alam+3.*(1.+2.*z)*al2)-2.*z/al2/z1
      f4=12.*z*(1.-2.*z)/al2*(1.-bigl)
      f5=(1.-2.*zp+2.*zp*zp)/alam/al1
     .+(4.*z*alam-6.*z2-al2)/al4*log(z1/al1/z)
      f6=2./al4*(4.*z*alam-6.*z2-al2)*log(zp1)
      cheavy=(f1+f2+f3+f4+f5+f6)*eps*al2
      return
    5 alam=1./(1.+eps)
      zp=z/alam
      alz=alam*z
      alz1=1.-alz
      al1=1.-alam
      z1=1.-z
      zp1=1.-zp
      pqg=0.5*(zp*zp+zp1*zp1)
      arg0=alam*zp1*zp1/z1/z
      c0=pqg*log(arg0)
      arg=alam*z1/al1/z
      bigl=log(arg)
      f1=c0+2.*al1*zp*zp1+al1*zp*bigl*(-2.*zp1+2.*z)
      f2=pqg*(-bigl-log(alam))
      cheavy=2.*(f1+f2)
      return
    6 alam=1./(1.+eps)
      zp=z/alam
      alz=alam*z
      alz1=1.-alz
      al1=1.-alam
      al2=alam*alam
      al3=al2*alam
      z1=1.-z
      zp1=1.-zp
      arg=alam/al1/z2
      f1=(1.-2.*zp+2.*zp*zp)*(-1./alam/al1)
      f2=2.*z/al3*(alam-2.*z)*log(arg)
      f3=4.*z/al3*(3.*z-2.*alam)
      f4=4.*z/al3*(alam-2.*z)*log(zp1)
      cheavy=(f1+f2+f3+f4)*eps*al2
      return
    7 cheavy=-bb*beta+c*eps*fac
      return
    8 cheavy=4.*cf*z*1.25*(1/(1+4.*eps)-0.2)
      return
   10 print 99
   99 format(1x,'x > x0')
      print 98,i,z,eps,beta2
   98 format(1x,i3,' z=',f10.6,' eps=',f10.6,' beta2=',f10.6)
      stop
      end


c Simplified version of NLO clg in FFNS Q^2<M^2
      function clgffnsl(z,eps)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) ! G.W. 15/02/2007
      beta2=1.-4.*eps*z/(1.-z)
      beta=sqrt(beta2)
      xi=1./eps
      pi = 3.14159265359d0
 
      term1 =beta**2*xi*(-22.06*z**2*(1+
     .2.24*log(xi)+1.298*(log(xi))**2+6.15*(log(xi))**3)
     .+2.091*(z*log(z))
     .*(1+3.175*log(xi)+0.7*(log(xi))**2)
     .+1566.8*(z**3)
     .*(1-1.923*log(xi)+1.931*(log(xi))**2)
     .-10.41*(z)
     .*(1-2.791*log(xi)-0.181*(log(xi))**2) 
     .+0.679*(1-0.00156*log(1./z))*(1+
     .0.2368*log(xi)-0.1357*(log(xi))**2-0.0484*(log(xi))**3))/z
 
      clgffnsl=term1 
  
      return
      end

c Simplified version of NLO clq in FFNS Q^2<M^2
      function clqffnsl(z,eps)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) ! G.W. 15/02/2007
      beta2=1.-4.*eps*z/(1.-z)
      beta=sqrt(beta2)
      xi=1./eps

      term1=beta**2*xi*(24.83*z**2*(1+
     .0.324*log(xi)+0.429*(log(xi))**2+0.0923*(log(xi))**3)
     .+1.091*(z*log(z))
     .*(1-1.251*log(xi)-0.3127*(log(xi))**2)
     .-5.358*(z**3)
     .*(1-1.684*log(xi)-5.652*(log(xi))**2)
     .-3.111*(z)
     .*(1+1.553*log(xi)+0.233*(log(xi))**2) 
     .+0.2443*(1+0.02402*log(1./z))*(1+
     .0.2075*log(xi)-0.1234*(log(xi))**2-0.0314*(log(xi))**3))/z
 
      clqffnsl=term1 
  
      return
      end

c Subtraction term for NLO clg in VFNS
      function clgvfsub(z,eps)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) ! G.W. 15/02/2007
      COMMON/iordCommon/iord
      beta2=1.-4.*eps*z/(1.-z)
      beta=sqrt(beta2)
      xi=1./eps

      term1 = 2./3.*(32.*z*log(z)+16.+16.*z-32.*z**2)*(log(xi))

      IF (IORD.EQ.2) THEN
         clgvfsub=(term1)*1.25*(1/(1+4.*eps)-0.2)
      ELSE
         clgvfsub=(term1)*(1.-eps**0.5)
      END IF

      return
      end

c Simplified version of NLO clg in FFNS Q^2>M^2
      function clgffnsh(z,eps)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) ! G.W. 15/02/2007
      beta2=1.-4.*eps*z/(1.-z)
      beta=sqrt(beta2)
      xi=1./eps

      term1 = 2./3.*(32.*z*log(z)+16.+16.*z-32.*z**2)*beta*(log(xi))

       DL  = log(z)
       DL1 = log(1.-z)

      term2 = ((94.74-49.2*z)*(1.-z)*DL1**2. + 864.8*(1.-z)*DL1 
     1 +1161.*z*DL1*DL+60.06*z*DL**2.+39.66*(1.-z)*DL
     2 -5.333*(1./z-1))*beta**3.

      term3 =beta**2.3/(z*xi)*(-155.41*z**2*(1+
     .2.057*log(xi)-0.017*(log(xi))**2-0.0639*(log(xi))**3)
     .-14.67*(z*log(z))
     .*(1-1.378*log(xi)-2.755*(log(xi))**2)
     .+948.72*(z**3)
     .*(1+1.144*log(xi)-0.400*(log(xi))**2)
     .+35.24*(z)
     .*(1+3.634*log(xi)+2.395*(log(xi))**2) 
     .+7.134*(1-0.0199*log(1./z))*(1+
     .1.144*log(xi)+0.8324*(log(xi))**2+0.4316*(log(xi))**3))

      clgffnsh=(term1+term2+term3)

      return
      end

c Simplified version of NLO clq in FFNS Q^2>M^2
      function clqffnsh(z,eps)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) ! G.W. 15/02/2007

      beta2=1.-4.*eps*z/(1.-z)
      beta=sqrt(beta2)
      xi=1/eps

       DL  = LOG (z)
       DL1 = LOG (1.-z)

      term1 = ((15.94-5.212*z)*(1.-z)**2.*DL1 + 
     1(0.421+1.520*z)*DL**2 +28.09*(1.-z)*DL
     2-(2.370/z-19.27)*(1.-z)**3. )*beta**3.

      term2 =beta*(-1.198*z**2*(1+
     .4.146*log(xi)-8.22*(log(xi))**2+1.236*(log(xi))**3)
     .+0.00335*(z*(log(z)))
     .*(1-7.135*log(xi)+12.282*(log(xi))**2)
     .-0.188*(z*(log(z))**2)
     .*(1+6.095*log(xi)-17.989*(log(xi))**2)
     .-4.8214*(z)
     .*(1-1.686*log(xi)+3.069*(log(xi))**2) 
     .+4.214*(1-0.03673*log(1./z))*(1+
     .1.1373*log(xi)+0.53361*(log(xi))**2
     .+0.37647*(log(xi))**3))/(z*xi)

      clqffnsh=(term1+term2)

      return
      end


c Simplified version of NLO c2g in FFNS Q^2<M^2
      function c2gffnsl(z,eps)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) ! G.W. 15/02/2007

      beta2=1.-4.*eps*z/(1.-z)
      beta=sqrt(beta2)
      xi=1./eps
      pi = 3.14159265359d0
      rho=1.-beta2 
 
      term1=rho*(1.5*0.25d0/pi*1.d0/(1.d0 
     .+ 0.25d0*(xi))*
     .(beta*(log(8.d0*beta*beta))**2- 5.d0*beta*log(8.d0*beta*beta) 
     .- 0.25d0*pi*pi)+2./3.*0.25d0/pi*1.d0/(1.d0 
     .+ 0.25d0*(xi))*pi*pi/2.d0)*xi*16.*pi/z

      term2 = beta*(-0.747*(xi)*(z*(log(z))**4)*(1+
     .1.626*log(xi)+1.008*(log(xi))**2-0.205*(log(xi))**3)
     .-39.34*(xi)*((log(1.-z)))
     .*(1-5.521*log(xi)+1.934*(log(xi))**2)
     .+867.6*(xi)*((log(1.-z))**2)
     .*(1-0.6569*log(xi)+0.1751*(log(xi))**2)
     .+1.978*(xi)*(z*(log(z))**3)
     .*(1-5.107*log(xi)-1.874*(log(xi))**2) 
     .+8.008*(xi)*(1+0.0333*log(1./z))*beta**14.238*(1-
     .0.5042*log(xi)-0.1053*(log(xi))**2+0.01987*(log(xi))**3)
     .+6.541*(xi)**0.4252*beta)/z
 
      c2gffnsl=term1+term2 
  
      return
      end

c Simplified version of NLO c2q in FFNS Q^2<M^2
      function c2qffnsl(z,eps)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) ! G.W. 15/02/2007

      beta2=1.-4.*eps*z/(1.-z)
      beta=sqrt(beta2)
      xi=1./eps

      term1=beta*(320.54*(xi)*z**2.*(1+
     .2.257*log(xi)+3.104*(log(xi))**2-0.5681*(log(xi))**3)
     .+11.37*(xi)*(z*(log(z)))
     .*(1+4.427*log(xi)-18.199*(log(xi))**2)
     .-9.518*(xi)*(z*(log(z))**2)
     .*(1-2.815*log(xi)+2.935*(log(xi))**2)
     .-12.684*(xi)*(z)
     .*(1+5.0125*log(xi)+34.086*(log(xi))**2) 
     .+4.654*(xi)*(1+0.05894*log(1./z))*(1-
     .0.6226*log(xi)-0.0299*(log(xi))**2
     .-0.00108*(log(xi))**3))/z
 
      c2qffnsl=term1 
  
      return
      end

c Simplified version of NLO c2g in FFNS Q^2>M^2
      function c2gffnsh(z,eps)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) ! G.W. 15/02/2007

      beta2=1.-4.*eps*z/(1.-z)
      beta=sqrt(beta2)
      xi=1./eps

      term1 = (2.d0/3.d0*((8.d0-16.*z+16.*z**2)*log(1.-z)
     x-(4.d0-8.*z+16.*z**2)*log(z)-2.d0+8.*z)
     x-3.d0/2.d0*((8.d0-16.*z+16.*z**2)*log(1.-z)
     x+(8.d0+32.*z)*log(z)+16./3./z+4.d0+32.*z-124./3.*z**2))
     x*beta*(log(xi))**2.

      term2 = (80./3./z +304.4*z-378.3*z**2+5.96*log(1.-z)
     x+33.9*log(z)
     x-0.277*(log(1.-z))**2-4.47*(log(z))**2+80.15*z**3
     x-153.2*log(z)*log(1.-z))*beta*(log(xi))

       DL  = log(z)
       DL1 = log(1.-z)

      term3 =   ( 1./z * (11.90 - 42.77* DL1) + 6.10 * DL**3  
     1         - 49.92 * DL**2 - 251.02 * DL - 1123.1 - 665.7* DL1
     2         + (214.4 - 215.1*z) * DL1**3 - 145.75 * DL1**2
     3         - 770.5 * DL**2 * DL1 - 935.04 * DL * DL1**2 )*beta

      term4 = (-224./9./z -10./9.*(log(1.-z))**3-316.15*z
     x+200.0*z**2-27.24*log(1.-z)-14.52*log(z)
     x-2.28*(log(1.-z))**2+13.21*(log(z))**2+96.77*z**3
     x+217.06*log(z)*log(1.-z))*beta

      term5 =(9.540*((log(1.-z))**4)*(1+
     .0.8336*log(xi)+1.072*(log(xi))**2
     .-0.06285*(log(xi))**3)
     .+3.853*(-(log(1.-z)))
     .*(1+9.5591*log(xi)+14.252*(log(xi))**2)
     .+365.18*((log(1.-z))**2)
     .*(1-0.202*log(xi)-0.1031*(log(xi))**2)
     .+3.125*((log(1.-z))**6)
     .*(1-1.578*log(xi)+0.2167*(log(xi))**2) 
     .+31.652*(1-0.002965*log(1./z))*beta**1.2823*(1+
     .0.31946*log(xi)+0.13524*(log(xi))**2
     .+0.06488*(log(xi))**3)-0.937*beta*xi**(-0.9556))/(z*xi)

      c2gffnsh=(term1+term2+term3+term4+term5)

      return
      end

c Simplified version of NLO c2q in FFNS Q^2>M^2
      function c2qffnsh(z,eps)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) ! G.W. 15/02/2007

      beta2=1.-4.*eps*z/(1.-z)
      beta=sqrt(beta2)
      xi=1/eps

      term1 = -(2.d0/3.d0*(8.*(1.+z)*log(z)+16./3./z+4.-4.*z
     .-16./3.*z**2))*beta*(log(xi))**2.

      term2 = -2./3.*(8.*(1.+z)*(log(z))**2
     .-(8.+40.*z+64./3.*z**2)*log(z)-160./9./z+16.-48.*z
     .+448./9.*z**2)*beta*(log(xi))

      term3 = (-896./81./z -88.76*z
     x+59.07*z**2+3.41*log(z)+8.078*(log(z))**2
     x+16.42*log(z)*(log(1.-z))**2+40.98*z**3 
     x+97.04*log(z)*log(1.-z))*beta

       DL  = LOG (z)
       DL1 = LOG (1.-z)

       term4 =   ( 5.290 * (1./z-1.) + 4.310 * DL**3   
     1         - 2.086 * DL**2 + 39.78 * DL - 0.101 * (1.-z) * DL1**3 
     2- (24.75 - 13.80 * z) * DL**2 * DL1 + 30.23 * DL * DL1)*beta

      term5 =(0.219*(1+
     .15.57*log(xi)+13.89*(log(xi))**2
     .+18.17*(log(xi))**3)
     .+0.0522*(-(log(z)))
     .*(1-12.47*log(xi)-30.14*(log(xi))**2)
     .-0.0075*((log(z))**2)
     .*(1-4.404*log(xi)-14.44*(log(xi))**2)
     .+1.86*(z)
     .*(1-5.697*log(xi)-8.0884*(log(xi))**2) 
     .+13.738*(1-0.00555*log(1./z))*beta*(1+
     .0.3379*log(xi)+0.3241*(log(xi))**2
     .-0.2286*(log(xi))**3))/(z*xi)

      c2qffnsh=(term1+term2+term3+term4+term5)

      return
      end
*
*
* ..This is the regular non-singlet piece for the electromagnetic F2, 
*    corresponding to CLNSP+C2NSM in W. van Neerven's program. The 
*    8 numerical coefficients are fitted to his results, using x values 
*    between 10^-6 and 1-10^-6. 
*
       FUNCTION CLNN2A (Y, NF)
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
*
       CLNN2A = 
     1          - 40.41 + 97.48 * Y
     2          + (26.56 * Y - 0.031) * DL**2 - 14.85 * DL 
     3          + 13.62 * DL1**2 - 55.79 * DL1 - 150.5 * DL * DL1 
     4        + NF * 16./27.D0 * ( 6.* Y*DL1 - 12.* Y*DL - 25.* Y + 6.)
*
       RETURN
       END

* ---------------------------------------------------------------------
*
*
* ..This is the 'local' NS piece for the e.m. FL, with no counterpart 
*    in WvN's program, as it does not exist in the exact expressions.
*    The value is fixed from the lowest integer moments.
*
       FUNCTION CLNN2C (Y)
       IMPLICIT REAL*8 (A-Z)
*
       CLNN2C = -0.164
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..This is the regular non-singlet piece for the electromagnetic F2, 
*    corresponding to C2NSP+C2NSN in W. van Neerven's program. The 
*    (10+8) numerical coefficients are fitted to his results, using x 
*    values between 10^-6 and 1-10^-6. 
*
      FUNCTION C2NN2A (Y, NF)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
*
       C2NN2A = 
     1          - 69.59 - 1008.* Y
     2          - 2.835 * DL**3 - 17.08 * DL**2 + 5.986 * DL 
     3          - 17.19 * DL1**3 + 71.08 * DL1**2 - 660.7 * DL1
     4          - 174.8 * DL * DL1**2 + 95.09 * DL**2 * DL1
     5        + NF * ( - 5.691 - 37.91 * Y 
     6          + 2.244 * DL**2 + 5.770 * DL 
     7          - 1.707 * DL1**2  + 22.95 * DL1
     8          + 3.036 * DL**2 * DL1 + 17.97 * DL * DL1 )     
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..This is the singular NS piece, denoted by SOFT2 in WvN's program. 
*    It is the same for all F2 and F3 cases. The numerical coefficients 
*    are exact, but truncated.
*
       FUNCTION C2NS2B (Y, NF)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

*
       DL1 = LOG (1.-Y)
       DM  = 1./(1.-Y)
*
       C2NS2B = 
     1          + 14.2222 * DL1**3 - 61.3333 * DL1**2 - 31.105 * DL1 
     2          + 188.64 
     3        + NF * ( 1.77778 * DL1**2 - 8.5926 * DL1 + 6.3489 ) 
       C2NS2B = DM * C2NS2B
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..This is the 'local' NS piece for the e.m. F2, denoted by COR2 in 
*    WvN's program. The numerical coefficients of the logs are exact,
*    but truncated, the constant one (from the delta-function) is 
*    slightly adjusted (+ 0.485 - 0.0035 NF) using the lowest moments.
*
       FUNCTION C2NN2C (Y, NF)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

       DL1 = LOG (1.-Y)
*
       C2NN2C = 
     1          + 3.55555 * DL1**4 - 20.4444 * DL1**3 - 15.5525 * DL1**2
     2          + 188.64 * DL1 - 338.531 + 0.485 
     3        + NF * (0.592593 * DL1**3 - 4.2963 * DL1**2 
     4          + 6.3489 * DL1 + 46.844 - 0.0035)
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..This is the pure singlet piece, denoted by C2S in WvN's program. 
*    Seven numerical coefficients (all but the one of 1/y, which is 
*    exact up to truncation) are fitted to his results, using x values
*    between 10^-6 and 1-10^-6.
*
      FUNCTION C2S2A (Y, NF)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
*
       C2S2A =   NF * ( 5.290 * (1./Y-1.) + 4.310 * DL**3   
     1         - 2.086 * DL**2 + 39.78 * DL - 0.101 * (1.-Y) * DL1**3 
     2         - (24.75 - 13.80 * Y) * DL**2 * DL1 + 30.23 * DL * DL1 )
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..This is the regular gluon piece, denoted by C2G2 in WvN's program. 
*    Nine numerical coefficients are fitted as above, the ones of 1/y, 
*    ln^3(1-y), and ln^2(1-y) are exact up to truncation.
*
      FUNCTION C2G2A (Y, NF)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
*
       C2G2A =   NF * ( 1./Y * (11.90 + 1494.* DL1) + 5.319 * DL**3  
     1         - 59.48 * DL**2 - 284.8 * DL + 392.4 - 1483.* DL1
     2         + (6.445 + 209.4 * (1.-Y)) * DL1**3 - 24.00 * DL1**2
     3         - 724.1 * DL**2 * DL1 - 871.8 * DL * DL1**2 )
*
       RETURN
       END
* 
* ---------------------------------------------------------------------
*
*
* ..This is the 'local' gluon piece, which has no counterpart in WvN's
*    program, as it does not exist in the exact expressions. Here it 
*    is, however, relevant for achieving a high accuracy of the convo-
*    lution, as are the adjustments of the constant in the non-singlet
*    quark coefficient functions. The value is fixed from the lowest 
*    even-integer moments. 
*
      FUNCTION C2G2C (Y, NF)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

       C2G2C = - NF * 0.28  
*
       RETURN
       END

c Approx version of NNLO c2g in FFNS Q^2<M^2
      function c2gffns3(z,eps)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      beta2=1.-4.*eps*z/(1.-z)
      beta=sqrt(beta2)
      xi=1./eps
      eta=xi*(1.-z)/(4.*z)-1
      xmax=(1./(1.+eps*4))
      pi=3.14159265359d0

      term1 = 96*beta*(13.073*xi-23.827*xi**2+24.107*xi**3
     .-9.173*xi**4)/z*(log(xmax/z)-4.)*(1.-z/xmax)**20

      term2=256./z*pi**3*xi*(1./(1.+xi/4))*2.1/(1.+eta)*
     .((0.0144-0.0273*log(eta)+0.00235*(log(eta))**2
     .-0.0001033*(log(eta))**3-0.0000478*(log(eta))**4)
     .+log(xi)*(0.0205-0.00373*log(eta)+0.00339*(log(eta))**2
     .+0.000128*(log(eta))**3-0.000044*(log(eta))**4)  
     .+(log(xi))**2*(-0.00065-0.0003*log(eta)
     .+0.000178*(log(eta))**2+0.0000206*(log(eta))**3))

      c2gffns3=(term1+term2)

      return
      end

c Approx version of NNLO c2q in FFNS Q^2<M^2
      function c2qffns3(z,eps)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      beta2=1.-4.*eps*z/(1.-z)
      beta=sqrt(beta2)
      xi=1./eps
      eta=1.-beta2
      xmax=(1./(1.+eps*4.))
      pi=3.14159265359d0

      term1 = 4./9.*96*beta*(13.073*xi-23.827*xi**2+24.107*xi**3
     .-9.173*xi**4)/z*(log(xmax/z)-4.)*(1.-z/xmax)**20

      c2qffns3=term1

      return
      end

c Approx version of NNLO clg in FFNS Q^2<M^2
      function clgffns3(z,eps)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      beta2=1.-4.*eps*z/(1.-z)
      beta=sqrt(beta2)
      xi=1./eps
      eta=xi*(1.-z)/(4.*z)-1
      xmax=(1./(1.+eps*4))
      pi=3.14159265359d0

      term1 = 96*beta**3*(0.484*xi**2-0.567*xi**3
     .+0.239*xi**4)/z*(log(xmax/z)-4.)*(1.-z/xmax)**20

      clgffns3=term1

      return
      end

c Approx version of NNLO clq in FFNS Q^2<M^2
      function clqffns3(z,eps)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      beta2=1.-4.*eps*z/(1.-z)
      beta=sqrt(beta2)
      xi=1./eps
      eta=1.-beta2
      xmax=(1./(1.+eps*4.))
      pi=3.14159265359d0

      term1 = 4./9.*96*beta**3*(0.484*xi**2-0.567*xi**3
     .+0.239*xi**4)/z*(log(xmax/z)-4.)*(1.-z/xmax)**20

      clqffns3=term1

      return
      end



* =======================================================================
*
* ..This is the regular non-singlet piece for the electromagnetic F2, 
*    corresponding to C2NSP+C2NSN in W. van Neerven's program. The 
*    (10+8) numerical coefficients are fitted to his results, using x 
*    values between 10^-6 and 1-10^-6. Heavy flvour part. 
*
      FUNCTION C2NN2AH (Y)
      
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

*
       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
*
       C2NN2AH = ( - 5.691 - 37.91 * Y 
     .          + 2.244 * DL**2 + 5.770 * DL 
     .          - 1.707 * DL1**2  + 22.95 * DL1
     .          + 3.036 * DL**2 * DL1 + 17.97 * DL * DL1 )     
*
       RETURN
       END

* ---------------------------------------------------------------------
*
*
* ..This is the 'local' NS piece for the e.m. F2, denoted by COR2 in 
*    WvN's program. The numerical coefficients of the logs are exact,
*    but truncated, the constant one (from the delta-function) is 
*    slightly adjusted (+ 0.485 - 0.0035 NF) using the lowest moments.
*  Heavy flavour piece.
      FUNCTION C2NN2CH (Y)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

       DL1 = LOG (1.-Y)
*
       C2NN2CH = (0.592593 * DL1**3 - 4.2963 * DL1**2 
     .          + 6.3489 * DL1 + 46.844 - 0.0035)
*
       RETURN
       END
*

* ---------------------------------------------------------------------
*
*
* ..This is the singular NS piece, denoted by SOFT2 in WvN's program. 
*    It is the same for all F2 and F3 cases. The numerical coefficients 
*    are exact, but truncated. Heavy flavour piece
*
      FUNCTION C2NS2BH (Y)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

*
       DL1 = LOG (1.-Y)
       DM  = 1./(1.-Y)
*
       C2NS2BH = ( 1.77778 * DL1**2 - 8.5926 * DL1 + 6.3489 ) 
       C2NS2BH = DM * C2NS2BH
*
       RETURN
       END
*



c Simplified version of NNLO l2q in FFNS Q^2<M^2
      function cl2ffnsl(z,eps)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      xi=1./eps
      beta2=1.-4.*eps*z/(1.-z)
      beta=sqrt(beta2)

      term1=beta**2*xi*(-78.24*z**2*(1+
     .0.388*log(xi)-5.232*(log(xi))**2-0.197*(log(xi))**3)
     .-1.458*(z*log(z))
     .*(1-6.136*log(xi)+8.047*(log(xi))**2)
     .+153.1*(z**3)
     .*(1-1.203*log(xi)-7.212*(log(xi))**2)
     .+2.8375*(z)
     .*(1+9.749*log(xi)-19.823*(log(xi))**2) 
     .+0.262*(1-0.1015*log(1./z))*(1-
     .0.0239*log(xi)+0.07244*(log(xi))**2+0.0435*(log(xi))**3))/z
 
      cl2ffnsl=term1 
  
      return
      end

c Simplified version of NLO l2 in FFNS Q^2>M^2
      function cl2ffnsh(z,eps)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      pi = 3.14159265359d0
      xi=1/eps

      term1 = 2.d0/3.d0*(4./3.*(1.+z**2)/(1.-z))*(log(xi))**2.

      term2 =2.d0/3.d0*((1.+z**2)/(1.-z)*(8./3.*log(1.-z)
     .-16./3.*log(z)-58./9.)+2./3.+26./3.*z)*(log(xi))

      term3 =2.d0/3.d0*((1.+z**2)/(1.-z)*(-8./3.*(1.0379*(1.-z)
     .+0.509*(1.-z)**3.-0.2713*(1.-z)**7.+0.3676*(1.-z)**10.)
     .-8./3.*pi**2./6.-16./3.*log(z)*log(1.-z)
     .+4./3.*(log(1.-z))**2.+4.*(log(z))**2.-58./9.*log(1.-z)
     .+134./9.*log(z)+359./27.)+(2./3.+26./3.*z)*log(1.-z)
     .-(2.+46./3.*z)*log(z)+29./9.-295./9.*z) 

      term4 =1/(z*xi)*(7.455*(log(1.-z)/(1.-z))*(1-
     .18.32*log(xi)+1.741*(log(xi))**2-0.235*(log(xi))**3)
     .+0.4824*(z/(1.-z))
     .*(1+49.627*log(xi)-11.47*(log(xi))**2)
     .-15.76*((log(1.-z)))
     .*(1-10.16*log(xi)+0.843*(log(xi))**2)
     .-9.023*(((log(1.-z))**2)/(1.-z))
     .*(1+5.997*log(xi)-0.84*(log(xi))**2+0.0688*(log(xi))**3) 
     .-17.11*((log(1.-z))**2.)*(1+
     .2.768*log(xi)+0.142*(log(xi))**2+0.104*(log(xi))**3))

      cl2ffnsh=(term1+term2+term3+term4)

      return
      end


c Subtraction term for NNLO c2g in VFNS
      function c2gvfsub(z,eps)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      beta2=1.-4.*eps*z/(1.-z)
      beta=sqrt(beta2)
      xi=1./eps

      term1 = (2.d0/3.d0*((8.d0-16.*z+16.*z**2)*log(1.-z)
     x-(4.d0-8.*z+16.*z**2)*log(z)-2.d0+8.*z)
     x-3.d0/2.d0*((8.d0-16.*z+16.*z**2)*log(1.-z)
     x+(8.d0+32.*z)*log(z)+16./3./z+4.d0+32.*z-124./3.*z**2))
     x*(log(xi))**2.

      term2 = (80./3./z +304.4*z-378.3*z**2+5.96*log(1.-z)
     x+33.9*log(z)
     x-0.277*(log(1.-z))**2-4.47*(log(z))**2+80.15*z**3
     x-153.2*log(z)*log(1.-z))*(log(xi))

c$$$ C--  Fit to (B.3) of hep-ph/9612398, but typo fixed in journal!
c$$$      term3 = (-224./9./z -10./9.*(log(1.-z))**3-316.15*z
c$$$     x+200.0*z**2-27.24*log(1.-z)-14.52*log(z)
c$$$     x-2.28*(log(1.-z))**2+13.21*(log(z))**2+96.77*z**3
c$$$     x+217.06*log(z)*log(1.-z))
      term3 = A2HGA(Z) ! G.W. 12/06/2008 Use A.Vogt's parameterisation.

      c2gvfsub=(term1+term2+term3)

      return
      end

c Subtraction term for NNLO c2q in VFNS
      function c2qvfsub(z,eps)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      beta2=1.-4.*eps*z/(1.-z)
      beta=sqrt(beta2)
      xi=1/eps

      term1 = -(2.d0/3.d0*(8.*(1.+z)*log(z)+16./3./z+4.-4.*z
     .-16./3.*z**2))*(log(xi))**2.

      term2 = -2./3.*(8.*(1.+z)*(log(z))**2
     .-(8.+40.*z+64./3.*z**2)*log(z)-160./9./z+16.-48.*z
     .+448./9.*z**2)*(log(xi))

      term3 = (-896./81./z -88.76*z
     x+59.07*z**2+3.41*log(z)+8.078*(log(z))**2
     x+16.42*log(z)*(log(1.-z))**2+40.98*z**3 
     x+97.04*log(z)*log(1.-z))


      c2qvfsub=(term1+term2+term3)

      return
      end


c Simplified version of LO clg convoluted with a2gg
      function clgconagg(z,eps)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      beta2=1.-4.*eps*z/(1.-z)
      beta=sqrt(beta2)
      xi=1./eps

      term1 = 3.85*(0.2-z)+16525*z*(z-0.2)*log(5*z)
     .+457.82*log(5*z)*log(1-5*z)-54222*(1-5.*z)*z**3
     .+17745*z**2*(log(5*z))**2*log(1-5*z)+5619*z*log(5.*z)
     .+282565*z**2*(0.2-z)+138.7*z*(1-5*z)-18965*z**2*(1-5*z)
     .-11346*z**2*log(5*z)*log(1-5*z)

      clgconagg=(term1)/z

      return
      end


! =====================================================================
!
! ..File: xa2hgp.f (provided by A.Vogt for G.Salam's HOPPET)
!
!
! ..Calculation of the alpha_s^2 heavy-quark singlet operator matrix
!    element (OME) a^2_Hg(x) in the MS(bar) scheme for mu_f = m_H via
!    a compact parametrization involving only logarithms.
!    The coupling constant is normalized as  a_s = alpha_s/(4*pi).
!
! ..This quantity, presented in Appendix B of M. Buza, Y. Matiounine,
!    J. Smith and W.L. van Neerven, Eur. Phys. J. C1 (1998) 301 (BSMN),
!    is required for the N_f matching of the NNLO parton densities.
!
!  ..The distributions (in the mathematical sense) are given as in eq.
!    (B.26) of Floratos, Kounnas, Lacaze: Nucl. Phys. B192 (1981) 417.
!    The name-endings A, B, and C of the functions below correspond to
!    the kernel superscripts [2], [3], and [1] in that equation.
!
!  ..The relative accuracy of the OME, as well as of its convolutions
!    with gluon distributions, amounts to a few thousandth.
!
!  ..The user should also cite the original calculation by BSMN.
!
! =====================================================================
!
!
! ..This is the regular piece.
!
      FUNCTION A2HGA (Y)
      IMPLICIT REAL*8 (A-Z)

      DL  = LOG (Y) 
      DL1 = LOG (1.-Y)
      
      A2HGA = - 24.89d0 / Y - 187.8d0 + 249.6d0 * Y
     &     - 146.8d0 * DL**2 * DL1
     &     - 1.556d0 * DL**3  - 3.292d0 * DL**2  - 93.68d0 * DL
     &     - 1.111d0 * DL1**3 - 0.400d0 * DL1**2 - 2.770d0 * DL1

      RETURN 
      END                                           
!
! ---------------------------------------------------------------------
!
! ..This is the 'local' piece, which has no counterpart in W. van
!    Neerven's program, as it does not exist in the exact expressions.
      
      FUNCTION A2HGC (Y) 
      IMPLICIT REAL*8 (A-Z) 
      
      A2HGC = - 0.006d0 
      
      RETURN 
      END
! =====================================================================
