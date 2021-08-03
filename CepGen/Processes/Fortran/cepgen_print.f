      ! CepGen: a central exclusive processes event generator
      ! Copyright (C) 2013-2021  Laurent Forthomme
      !
      ! This program is free software: you can redistribute it and/or modify
      ! it under the terms of the GNU General Public License as published by
      ! the Free Software Foundation, either version 3 of the License, or
      ! any later version.
      !
      ! This program is distributed in the hope that it will be useful,
      ! but WITHOUT ANY WARRANTY; without even the implied warranty of
      ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
      ! GNU General Public License for more details.
      !
      ! You should have received a copy of the GNU General Public License
      ! along with this program.  If not, see <http://www.gnu.org/licenses/>.

      !> \file cepgen_print.f
      subroutine cepgen_print
      !> Print useful run information in standard stream
      implicit none
      include 'CepGen/Processes/Fortran/KTBlocks.inc'
      logical params_shown
      data params_shown/.false./
      save params_shown

      if(params_shown) return

      print *,'========================================================'
      print *,'Parameter                                  value(s)'
      print *,'--------------------------------------------------------'
      print 101,'Process mode:',icontri
      print 101,'Structure functions:',sfmod
      print 103,'Beams momenta:',inp1,inp2
      print 102,'Fluxes modes:',iflux1,iflux2
      print 104,'Beams (A,Z):',a_nuc1,z_nuc1,a_nuc2,z_nuc2
      print *,'========================================================'
      print *,'Cut                        enabled   minimum     maximum'
      print *,'--------------------------------------------------------'
      print 100,'pt(single)',ipt,pt_min,pt_max
      print 100,'energy(single)',iene,ene_min,ene_max
      print 100,'eta(single)',ieta,eta_min,eta_max
      print 100,'m(sum)',iinvm,invm_min,invm_max
      print 100,'pt(sum)',iptsum,ptsum_min,ptsum_max
      print 100,'delta(y)',idely,dely_min,dely_max
      print *,'========================================================'
      print *,'Process-specific parameters'
      call cepgen_list_params
      print *,'========================================================'

      params_shown=.true.

100   format(A26,'     ',L2,f12.4,f12.4)
101   format(A33,I12)
102   format(A33,I12,I12)
103   format(A33,f12.2,f12.2)
104   format(A33,'   (',I3,',',I3,'),  (',I3,',',I3,')')

      end

