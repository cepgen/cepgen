        subroutine F2_fit_luxlike(xbj,q2,F2,FL)
c       input: x,q2
c       output: F2,F1
        implicit real*8 (a-h,o-z)
        common/luxlike_params/amp,am_pi,alpha_em,
     &     q2_cut,w2_lo,w2_hi,
     &     ires_model,icont_model
c       -----------------------------
        w_thr = amp+am_pi
        w2 = amp**2 + q2*(1.d0-xbj)/xbj
        w = dsqrt(w2)
c       -----------------------------
        omega = (w2-w2_lo)/(w2_hi-w2_lo)
        rho = 2.d0*omega**2 - omega**4 

        if(q2.ge.q2_cut) then
           if(w2.gt.w2_hi) then ! MSTW grid
              call CepGen_Structure_Functions(205,Q2,xbj,F2p,FLp)
              F2 = F2p
              FL = FLp
           elseif(w2.lt.w2_hi) then
              call F2_cont(xbj,q2,F2c,FLc)
              F2 = F2c
              FL = FLc
           endif
        else ! Q2 < Q2cut
           if(w2.le.w2_lo) then
              if(ires_model.eq.1) then     ! Christy-Bosted
                 call CepGen_Structure_Functions(102,Q2,xbj,F2r,FLr)
              elseif(ires_model.eq.2) then ! Fiore-Brasse
                 call CepGen_Structure_Functions(101,Q2,xbj,F2r,FLr)
              endif
              F2 = F2r
              FL = FLr
           elseif(w2.gt.w2_lo.and.w2.lt.w2_hi) then
              if(ires_model.eq.1) then     ! Christy-Bosted
                 call CepGen_Structure_Functions(102,Q2,xbj,F2r,FLr)
              elseif(ires_model.eq.2) then ! Fiore-Brasse
                 call CepGen_Structure_Functions(101,Q2,xbj,F2r,FLr)
              endif
              call F2_cont(xbj,q2,F2c,FLc)
              F2 = (1.d0-rho)*F2r + rho*F2c
              FL = (1.d0-rho)*FLr + rho*FLc
           elseif(w2.ge.w2_hi) then
              call F2_cont(xbj,q2,F2c,FLc)
              F2 = F2c
              FL = FLc
           endif
        endif

        return
        end
c       ---------------------------------------------------------
        subroutine F2_CONT(xbj,q2,F2c,FLc)
        implicit real*8 (a-h,o-z)
        common/luxlike_params/amp,am_pi,alpha_em,
     &     q2_cut,w2_lo,w2_hi,
     &     ires_model,icont_model
        if(icont_model.eq.1) then     ! GD11p
           call CepGen_Structure_Functions(204,Q2,xbj,F2,FL)
        elseif(icont_model.eq.2) then ! ALLM91
           call CepGen_Structure_Functions(201,Q2,xbj,F2,FL)
        elseif(icont_model.eq.3) then ! ALLM97
           call CepGen_Structure_Functions(202,Q2,xbj,F2,FL)
        endif
        F2c = F2
        R = R_1998(xbj,Q2)
        tau = 4.d0*xbj**2*amp**2/q2
        FLc = F2c*(1+tau)*R/(1.d0+R)

        return
        end
c       -------------------------------------------------
        double precision function R_1998(x,q2)
        implicit real*8 (a-h,o-z)

c       from hep-ex/9808028v1, K.Abe et al.(SLAC)

        data a1,a2,a3/0.0485d0,0.5470d0,2.0621d0/
        data a4,a5,a6/-0.3804d0,0.5090d0,-0.0285d0/
        data b1,b2,b3/0.0481d0,0.6114d0,-0.3509d0/
        data b4,b5,b6/-0.4611d0,0.7172d0,-0.0317d0/ 
        data c1,c2,c3/0.0577d0,0.4644d0,1.8288d0/
        data c4,c5,c6/12.3708d0,-43.1043d0,41.7415d0/

        q2_b = 0.34
        u = q2/q2_b
        xl = dlog(q2/0.04d0)

        pa = (1.d0 + a4*x +a5*x**2)*x**a6
        pb = (1.d0 + b4*x +b5*x**2)*x**b6
        tt = theta(x,q2)
        q2_thr = c4*x +c5*x**2+c6*x**3

        ra = a1/xl*tt+a2/(q2**4.d0+a3**4.d0)**.25d0*pa
        rb = b1/xl*tt+(b2/q2+b3/(q2**2+.3d0**2))*pb
        rc = c1/xl*tt+c2*((q2-q2_thr)**2 +c3**2)**(-.5d0)

       
        R =(ra+rb+rc)/3.d0

        if(q2.gt.q2_b) then
                R_1998 = R
        else
                R_1998 = R*(3.d0*u - u**3)/2.d0
        endif

        return
        end

c       -------------------------------------------------
        double precision function theta(x,q2)
c       
c       function needed for the 1998 parametrisation of R=sigma_L/sigma_T
c
        implicit real*8(a-h,o-z)
        theta = 1.d0 + 12.d0*q2/(q2+1.d0)
     >                  *(0.125d0**2/(.125d0**2 +x**2))
        return
        end
c       --------------------------------------------------
