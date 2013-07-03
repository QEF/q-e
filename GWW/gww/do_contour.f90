!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

!this subroutine add to the integral part of the self-energy the poles part
  SUBROUTINE do_contour(sr,wp,options)
!NOT_TO_BE_INCLUDED_START
    USE  contour, ONLY : w_poles, w_poles_value 
    USE kinds, ONLY : DP
    USE self_energy_storage, ONLY : self_on_real
    USE basic_structures, ONLY : wannier_u,free_memory, initialize_memory
    USE input_gw,          ONLY : input_options

    implicit none

    TYPE(self_on_real), INTENT(inout) :: sr
    TYPE(w_poles), INTENT(in) :: wp
    TYPE(input_options), INTENT(in) :: options


    TYPE(wannier_u) :: uu
    INTEGER :: ie,jj,ii,is
    COMPLEX(kind=DP) :: energy

!reads KS eigen-energies

    call read_data_pw_u(uu,options%prefix)

!loop on spin
    do is=1,sr%nspin
!loop on real energy grid
       do ie=1,sr%n
          energy=sr%grid(ie)
!divide by in valence and in conduction case
          if(dble(energy) <= uu%ene(uu%nums_occ(is),is)) then
!consider valece states
             do jj=sr%i_min,uu%nums_occ(is)!ATTENZIONE
!loop on poles
!for selected poles add terms
                if(uu%ene(jj,is)>=dble(energy) )then
                   do ii=sr%i_min,sr%i_max
                      sr%diag(ie,ii,is)=sr%diag(ie,ii,is)-w_poles_value(uu%ene(jj,is)-energy,wp,jj,ii,is)!GIUSTO CUSSI'
                   enddo
                endif
             enddo
          else
             do jj=uu%nums_occ(is)+1,sr%i_max
!loop on poles
!for selected poles add terms 
                if(uu%ene(jj,is)<=dble(energy) )then
                   do ii=sr%i_min,sr%i_max
                      sr%diag(ie,ii,is)=sr%diag(ie,ii,is)+w_poles_value(uu%ene(jj,is)-energy,wp,jj,ii,is)
                   enddo
                endif
             enddo

          endif
       enddo
    enddo
    
    call free_memory(uu)
    return

!NOT_TO_BE_INCLUDED_END
  END SUBROUTINE do_contour
