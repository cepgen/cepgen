      program main

      implicit none
      integer i,niter,nsf
      double precision f2,fl
      double precision xbj,min_xbj,max_xbj
      double precision q2

c      nsf=204 ! GD11p
      nsf=301 ! LUXlike
      q2=10.225
      min_xbj=1.0d-3
      max_xbj=1.0
      niter=101

      do i=1,niter
         xbj=min_xbj+(max_xbj-min_xbj)*(i-1)/(niter-1)
         call CepGen_Structure_Functions(nsf,q2,xbj,f2,fl)
         print *,xbj,f2
      enddo
      end
