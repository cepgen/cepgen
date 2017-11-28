c       --------------------------------------------------
        subroutine F2_fit_luxlike(xbj,q2,F2,FL)
c       --------------------------------------------------
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
c in the Lux-paper a "higher-twist" correction is applied to F2
              F2 = F2p*(1.d0+5.5d0/Q2)
              FL = FLp
           else
              call F2_cont(xbj,q2,icont_model,F2c,FLc)
              F2 = F2c
              FL = FLc
           endif
        else ! Q2 < Q2cut
           if(w2.le.w2_lo) then
              call F2_res(xbj,q2,ires_model,F2r,FLr)
              F2 = F2r
              FL = FLr
           elseif(w2.gt.w2_lo.and.w2.le.w2_hi) then
              call F2_res(xbj,q2,ires_model,F2r,FLr)
              call F2_cont(xbj,q2,icont_model,F2c,FLc)
              F2 = (1.d0-rho)*F2r + rho*F2c
              FL = (1.d0-rho)*FLr + rho*FLc
           elseif(w2.gt.w2_hi) then
              call F2_cont(xbj,q2,icont_model,F2c,FLc)
              F2 = F2c
              FL = FLc
           endif
        endif

        return
        end
c       ---------------------------------------------------------
        subroutine F2_RES(xbj,q2,irm,F2r,FLr)
        implicit real*8 (a-h,o-z)
        if(irm.eq.1) then     ! Christy-Bosted
           call CepGen_Structure_Functions(102,Q2,xbj,F2,FL)
        elseif(irm.eq.2) then ! Fiore-Brasse
           call CepGen_Structure_Functions(101,Q2,xbj,F2,FL)
        endif
        F2r = F2
        FLr = FL

        return
        end
c       ---------------------------------------------------------
        subroutine F2_CONT(xbj,q2,icm,F2c,FLc)
        implicit real*8 (a-h,o-z)
        if(icm.eq.1) then     ! GD11p
           call CepGen_Structure_Functions(204,Q2,xbj,F2,FL)
        elseif(icm.eq.2) then ! ALLM91
           call CepGen_Structure_Functions(201,Q2,xbj,F2,FL)
        elseif(icm.eq.3) then ! ALLM97
           call CepGen_Structure_Functions(202,Q2,xbj,F2,FL)
        endif
        F2c = F2
        FLc = FL

        return
        end
c       --------------------------------------------------

