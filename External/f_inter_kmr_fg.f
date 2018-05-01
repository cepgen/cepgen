      subroutine f_inter_kmr_fg(rx,rkt2,rmu2,iread,parton)
      implicit double precision (a-h,o-z)
c     =================================================================
c     interpolation routine of the KMR UGDF from grids
c     =================================================================
      integer iread
      double precision rx,rkt2,rmu2,parton
      common/transferg/fmap_kmr_fg(140,140,140)

c     -----------------------------------------------------------------
c     rx = log(x)
c     -----------------------------------------------------------------
      rx_min = -7.d0
      rx_max =  0.d0
      nrx = 140
      drx = (rx_max-rx_min)/nrx
c     -----------------------------------------------------------------
c     rkt2 = log(kt2)
c     -----------------------------------------------------------------
      rkt2_min = -1.d0
      rkt2_max =  6.d0
      nrkt2 = 140
      drkt2 = (rkt2_max-rkt2_min)/nrkt2

c     -----------------------------------------------------------------
c     rmu2 = log(mu2)
c     -----------------------------------------------------------------
      rmu2_min = -1.d0
      rmu2_max =  6.d0
      nrmu2 = 140
      drmu2 = (rmu2_max-rmu2_min)/nrmu2

c     =================================================================
c     first time, only reading
c     =================================================================

      if(iread.eq.0) then

      do i = 1,140
      do j = 1,140
      do k = 1,140

      fmap_kmr_fg(i,j,k) = 0.d0

      enddo
      enddo
      enddo

c     =================================================================
c     reading file with the grid
c     =================================================================
c      open(unit=11,file='KMRgrid_gluon_mstw08.dat',
c     &     status='unknown')

cc      open(unit=11,
cc     &   file='gluon_mmht2014lo_Watt.dat',
cc     &     status='unknown')

      open(unit=11,
     &   file='gluon_mmht2014nlo_Watt.dat',
     &     status='unknown')


c      open(unit=11,
c     &   file='KMRWatt_gluon_mmht2014_LO.dat',
c     &     status='unknown')

c     =================================================================

      do 100 irx=1,nrx

      rx = rx_min + irx*drx - drx/2.

      do 100 irkt2=1,nrkt2

      rkt2 = rkt2_min + irkt2*drkt2 - drkt2/2.

      do 100 irmu2=1,nrmu2

      rmu2 = rmu2_min + irmu2*drmu2 - drmu2/2.

      read(11,*) aln_x, aln_kt2, aln_mu2, rglue

      fmap_kmr_fg(irx,irkt2,irmu2) = rglue

  100 continue

      close(unit=11)

      elseif(iread.eq.1) then

c    ==================================================================
c     outside of the grid
c    ==================================================================

      if(rmu2.lt.(rmu2_min+drmu2/2.).or.rmu2.gt.(rmu2_max-drmu2/2.))
     & then
      parton = 0.d0
      goto 200
      endif

      if(rkt2.lt.(rkt2_min+drkt2/2.).or.rkt2.gt.(rkt2_max-drkt2/2.))
     & then
      parton = 0.d0
      goto 200
      endif

      if(rx.lt.(rx_min+drx/2.).or.rx.gt.(rx_max-drx/2.))
     & then
      parton = 0.d0
      goto 200
      endif

c     =================================================================
c     interpolation variables
c     =================================================================

      delrx = (rx-rx_min)
      delrkt2 = (rkt2-rkt2_min)
      delrmu2 = (rmu2-rmu2_min)

      srx = (delrx-drx/2.) / drx
      irx = int(srx)
      irx_lo = irx+1
      irx_up = irx_lo+1

      srkt2 = (delrkt2-drkt2/2.) / drkt2
      irkt2 = int(srkt2)
      irkt2_lo = irkt2+1
      irkt2_up = irkt2_lo+1

      srmu2 = (delrmu2-drmu2/2.) / drmu2
      irmu2 = int(srmu2)
      irmu2_lo = irmu2+1
      irmu2_up = irmu2_lo+1

c     =================================================================
c     for testing the code only
c     =================================================================

      rx_lo = rx_min+irx_lo*drx-drx/2.
      rx_up = rx_min+irx_up*drx-drx/2.
      rkt2_lo = rkt2_min+irkt2_lo*drkt2-drkt2/2.
      rkt2_up = rkt2_min+irkt2_up*drkt2-drkt2/2.
      rmu2_lo = rmu2_min+irmu2_lo*drmu2-drmu2/2.
      rmu2_up = rmu2_min+irmu2_up*drmu2-drmu2/2.

c     =================================================================
c     neighbouring grid points
c     =================================================================

      f111 = fmap_kmr_fg(irx_lo,irkt2_lo,irmu2_lo)
      f112 = fmap_kmr_fg(irx_lo,irkt2_lo,irmu2_up)
      f121 = fmap_kmr_fg(irx_lo,irkt2_up,irmu2_lo)
      f122 = fmap_kmr_fg(irx_lo,irkt2_up,irmu2_up)
      f211 = fmap_kmr_fg(irx_up,irkt2_lo,irmu2_lo)
      f212 = fmap_kmr_fg(irx_up,irkt2_lo,irmu2_up)
      f221 = fmap_kmr_fg(irx_up,irkt2_up,irmu2_lo)
      f222 = fmap_kmr_fg(irx_up,irkt2_up,irmu2_up)

c     =================================================================
c     weights for neighbouring points
c     =================================================================

      ddrx = (srx - float(irx))*drx
      ddrkt2  = (srkt2 - float(irkt2))*drkt2
      ddrmu2  = (srmu2 - float(irmu2))*drmu2

      w111 =   (drx-ddrx)      /drx
     *       * (drkt2-ddrkt2)  /drkt2
     *       * (drmu2-ddrmu2)/drmu2
      w112 =   (drx-ddrx)      /drx
     *       * (drkt2-ddrkt2)    /drkt2
     *       *   ddrmu2        /drmu2
      w121 =   (drx-ddrx)      /drx
     *       *   ddrkt2             /drkt2
     *       * (drmu2-ddrmu2)/drmu2
      w122 =   (drx-ddrx)      /drx
     *       *   ddrkt2             /drkt2
     *       *   ddrmu2         /drmu2
      w211 =     ddrx            /drx
     *       * (drkt2-ddrkt2)        /drkt2
     *       * (drmu2-ddrmu2)/drmu2
      w212 =     ddrx            /drx
     *       * (drkt2-ddrkt2)        /drkt2
     *       *   ddrmu2         /drmu2
      w221 =     ddrx            /drx
     *       *   ddrkt2             /drkt2
     *       * (drmu2-ddrmu2)/drmu2
      w222 =     ddrx            /drx
     *       *   ddrkt2             /drkt2
     *       *   ddrmu2         /drmu2
c     =================================================================
      f_int =  w111*f111 + w112*f112
     *       + w121*f121 + w122*f122
     *       + w211*f211 + w212*f212
     *       + w221*f221 + w222*f222
c     =================================================================

      parton = f_int

  200 continue

      endif

      return

      end

