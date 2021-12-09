
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine stres_loc (sigmaloc)
  !----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE atom,                 ONLY : msh, rgrid
  USE m_gth,                ONLY : dvloc_gth
  USE ions_base,            ONLY : ntyp => nsp
  USE cell_base,            ONLY : omega, tpiba2
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : fwfft
  USE gvect,                ONLY : ngm, gstart, g, ngl, gl, igtongl
  USE scf,                  ONLY : rho
  USE vlocal,               ONLY : strf, vloc
  USE control_flags,        ONLY : gamma_only
  USE wavefunctions, ONLY : psic
  USE uspp_param,           ONLY : upf
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE Coul_cut_2D,          ONLY : do_cutoff_2D, cutoff_stres_evloc, cutoff_stres_sigmaloc 
  !
  implicit none
  !
  real(DP) :: sigmaloc (3, 3)
  real(DP) , allocatable :: dvloc(:)
  real(DP) :: evloc, fact
  integer :: ng, nt, l, m
  ! counter on g vectors
  ! counter on atomic type
  ! counter on angular momentum
  ! counter on spin components

  allocate(dvloc(ngl))
  sigmaloc(:,:) = 0.d0
  psic(:) = CMPLX(rho%of_r(:,1), KIND=dp)

  CALL fwfft ('Rho', psic, dfftp)
  ! psic contains now the charge density in G space
  if (gamma_only) then
     fact = 2.d0
  else
     fact = 1.d0
  end if
  evloc = 0.0d0
  do nt = 1, ntyp
     if (gstart==2) evloc = evloc + &
          psic (dfftp%nl (1) ) * strf (1, nt) * vloc (igtongl (1), nt)
     do ng = gstart, ngm
        evloc = evloc +  DBLE (CONJG(psic (dfftp%nl (ng) ) ) * strf (ng, nt) ) &
             * vloc (igtongl (ng), nt) * fact
     enddo
  enddo
  ! 2D:  add contribution from cutoff long-range part of Vloc
  IF (do_cutoff_2D)  call cutoff_stres_evloc( psic, strf, evloc )
  !
  !      WRITE( 6,*) ' evloc ', evloc, evloc*omega   ! DEBUG
  !
  do nt = 1, ntyp
     IF ( upf(nt)%is_gth ) THEN
        !
        ! special case: GTH pseudopotential
        !
        call dvloc_gth (nt, upf(nt)%zp, tpiba2, ngl, gl, omega, dvloc)
        !
     ELSE IF ( upf(nt)%tcoulombp ) THEN
        !
        ! special case: pseudopotential is coulomb 1/r potential
        !
        call dvloc_coul (upf(nt)%zp, tpiba2, ngl, gl, omega, dvloc)
        !
     ELSE
        !
        ! normal case: dvloc contains dV_loc(G)/dG
        !
        call dvloc_of_g (rgrid(nt)%mesh, msh (nt), rgrid(nt)%rab, rgrid(nt)%r,&
          upf(nt)%vloc(:), upf(nt)%zp, tpiba2, ngl, gl, omega, dvloc)
        !
     END IF
     ! no G=0 contribution
     do ng = 1, ngm
        do l = 1, 3
           do m = 1, l
              sigmaloc(l, m) = sigmaloc(l, m) +  DBLE( CONJG( psic(dfftp%nl(ng) ) ) &
                    * strf (ng, nt) ) * 2.0d0 * dvloc (igtongl (ng) ) &
                    * tpiba2 * g (l, ng) * g (m, ng) * fact
           enddo
        enddo
     enddo
  enddo
  IF (do_cutoff_2D)  call cutoff_stres_sigmaloc( psic, strf, sigmaloc) ! 2D: re-add LR Vloc to sigma here
  !
  do l = 1, 3
     sigmaloc (l, l) = sigmaloc (l, l) + evloc
     do m = 1, l - 1
        sigmaloc (m, l) = sigmaloc (l, m)
     enddo
  enddo
  !
  call mp_sum(  sigmaloc, intra_bgrp_comm )
  !
  deallocate(dvloc)
  return
end subroutine stres_loc

