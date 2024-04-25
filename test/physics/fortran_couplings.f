      program main
      implicit none
      integer i, niter
      double precision min_q, max_q, q
      double precision CepGen_AlphaEM, CepGen_AlphaS

      call CepGen_init

      min_q = 1.0d0
      max_q = 1.0d2
      niter = 101

      do i = 1, niter
        q = min_q + (max_q - min_q) * (i - 1) / (niter - 1)
        print *, q, CepGen_AlphaEM('burkhardt', q),
     &              CepGen_AlphaS('pegasus', q)
      enddo
      call exit(0)
      end
