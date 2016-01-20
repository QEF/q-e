!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------------
subroutine set_irr_sym_new ( t, tmq, npertx )
!---------------------------------------------------------------------
!
!     This subroutine computes: 
!     1) the matrices which represent the small group of q on the
!        pattern basis.
!
  USE kinds, ONLY : DP
  USE constants, ONLY: tpi
  USE ions_base, ONLY : nat
  USE cell_base, ONLY : at, bg
  USE symm_base, ONLY : s, irt
  USE modes,     ONLY : u, nirr, npert
  USE control_flags, ONLY : modenum
  USE mp,        ONLY : mp_bcast
  USE mp_images, ONLY : intra_image_comm
  USE io_global, ONLY : ionode_id

  USE qpoint,       ONLY : xq
  USE lr_symm_base, ONLY : nsymq, irotmq, rtau, minus_q

  implicit none
!
!   first the dummy variables
!
  integer, intent(in) ::  npertx
! input: maximum dimension of the irreducible representations 
!
  complex(DP), intent(out) :: t(npertx, npertx, 48, 3*nat), &
                              tmq (npertx, npertx, 3*nat)
! output: the symmetry matrices
! output: the matrice sending q -> -q+G
!
!   here the local variables
!
  integer :: na, imode, jmode, ipert, jpert, kpert, nsymtot, imode0, &
       irr, ipol, jpol, isymq, irot, sna
  ! counters and auxiliary variables

  real(DP) :: arg
! the argument of the phase

  complex(DP) :: wrk_u (3, nat), wrk_ru (3, nat), fase, wrk
! pattern
! rotated pattern
! the phase factor

!
!   We compute the matrices which represent the symmetry transformation
!   in the basis of the displacements
!
  t(:,:,:,:) = (0.d0, 0.d0)
  tmq(:,:,:) = (0.d0, 0.d0)
  if (minus_q) then
     nsymtot = nsymq + 1
  else
     nsymtot = nsymq
  endif
  do isymq = 1, nsymtot
     if (isymq.le.nsymq) then
        irot = isymq
     else
        irot = irotmq
     endif
     imode0 = 0
     do irr = 1, nirr
        do ipert = 1, npert (irr)
           if (modenum /= 0 .AND. modenum /= irr) CYCLE
           imode = imode0 + ipert
           do na = 1, nat
              do ipol = 1, 3
                 jmode = 3 * (na - 1) + ipol
                 wrk_u (ipol, na) = u (jmode, imode)
              enddo
           enddo
!
!     transform this pattern to crystal basis
!
           do na = 1, nat
              call trnvecc (wrk_u (1, na), at, bg, - 1)
           enddo
!
!     the patterns are rotated with this symmetry
!
           wrk_ru(:,:) = (0.d0, 0.d0)
           do na = 1, nat
              sna = irt (irot, na)
              arg = 0.d0
              do ipol = 1, 3
                 arg = arg + xq (ipol) * rtau (ipol, irot, na)
              enddo
              arg = arg * tpi
              if (isymq.eq.nsymtot.and.minus_q) then
                 fase = CMPLX (cos (arg), sin (arg) )
              else
                 fase = CMPLX (cos (arg), - sin (arg) )
              endif
              do ipol = 1, 3
                 do jpol = 1, 3
                    wrk_ru (ipol, sna) = wrk_ru (ipol, sna) + s (jpol, ipol, irot) &
                         * wrk_u (jpol, na) * fase
                 enddo
              enddo
           enddo
!
!    Transform back the rotated pattern
!
           do na = 1, nat
              call trnvecc (wrk_ru (1, na), at, bg, 1)
           enddo
            
!
!     Computes the symmetry matrices on the basis of the pattern
!
           do jpert = 1, npert (irr)
              imode = imode0 + jpert
              do na = 1, nat
                 do ipol = 1, 3
                    jmode = ipol + (na - 1) * 3
                    if (isymq.eq.nsymtot.and.minus_q) then
                       tmq (jpert, ipert, irr) = tmq (jpert, ipert, irr) + CONJG(u ( &
                            jmode, imode) * wrk_ru (ipol, na) )
                    else
                       t (jpert, ipert, irot, irr) = t (jpert, ipert, irot, irr) &
                            + CONJG(u (jmode, imode) ) * wrk_ru (ipol, na)
                    endif
                 enddo
              enddo
           enddo
        enddo
        imode0 = imode0 + npert (irr)

!
! If the representations are irreducible, the rotations should be unitary matrices
! if this is not the case, the way the representations have been chosen has failed
! for some reasons (check set_irr.f90)
!
        
        if(isymq<=nsymq) then 
        do ipert = 1, npert (irr)
           IF (modenum /= 0 .AND. modenum /= irr) CYCLE
           do jpert = 1, npert (irr)
              wrk = cmplx(0.d0,0.d0)
              do kpert = 1, npert (irr)
                 wrk = wrk + t (ipert,kpert,irot,irr) * conjg( t(jpert,kpert,irot,irr))
              enddo
              if (jpert.ne.ipert .and. abs(wrk)> 1.d-6 ) &
                     call errore('set_irr_sym_new','wrong representation',100*irr+10*jpert+ipert)
              if (jpert.eq.ipert .and. abs(wrk-1.d0)> 1.d-6 ) &
                     call errore('set_irr_sym_new','wrong representation',100*irr+10*jpert+ipert)
           enddo
        enddo
        endif

     enddo
  enddo

!
! parallel stuff: first node broadcasts everything to all nodes
!
  call mp_bcast (t, ionode_id, intra_image_comm)
  call mp_bcast (tmq, ionode_id, intra_image_comm)
  return
end subroutine set_irr_sym_new
