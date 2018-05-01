      subroutine cepgen_print
      implicit none
      include 'cepgen_blocks.inc'
      logical params_shown
      data params_shown/.false./
      save params_shown

      if(params_shown) return

      print *,'List of cuts retrieved from CepGen:'
      print *,'pt cuts: ',pt_min,pt_max
      print *,'delta(y) cuts: ',dely_min,dely_max

      params_shown=.true.

      end

