      subroutine cepgen_print
      implicit none
      include 'cepgen_blocks.inc'
      logical params_shown
      data params_shown/.false./
      save params_shown

      if(params_shown) return

      print *,'========================================================'
      print *,'List of cuts retrieved from CepGen:'
      print *,'========================================================'
      print *,'Cut                               minimum     maximum'
      print 100,'pt',ipt,pt_min,pt_max
      print 100,'eta',ieta,eta_min,eta_max
      print 100,'delta(y)',idely,dely_min,dely_max
      print *,'========================================================'

      params_shown=.true.

100   format(A25,' (',L1,'):',f12.4,f12.4)

      end

