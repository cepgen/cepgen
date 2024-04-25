      program main
      implicit none
      integer i, niter
      double precision f2_ll, fl_ll, f2_sy, fl_sy
      double precision xbj, min_xbj, max_xbj
      double precision q2

      call CepGen_init

      q2 = 10.225
      min_xbj = 1.0d - 3
      max_xbj = 1.0
      niter = 101

      do i = 1, niter
        xbj = min_xbj + (max_xbj - min_xbj) * (i - 1) / (niter - 1)
        call CepGen_Structure_Functions(301, xbj, q2, f2_ll, fl_ll)
        call CepGen_Structure_Functions(11, xbj, q2, f2_sy, fl_sy)
        print *, q2, xbj, f2_ll, fl_ll, f2_sy, fl_sy
      enddo
      call exit(0)
      end
