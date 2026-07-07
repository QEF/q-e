!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
MODULE lr_magnons_routines
  !---------------------------------------------------------------------------
  !
 CONTAINS
  !
  !
  !---------------------------------------------------------------------------
  SUBROUTINE pauli( evc_g, ipol )
    !---------------------------------------------------------------------------
    !
    ! Applies a Pauli matrix or (i sigma_y) to a wavefunction in G-space
    ! 
    ! if ipol=1 output is
    !
    !   |  0     1  |  | evc_1 |   | evc_2 |
    !   |           |  |       | = |       |
    !   |  1     0  |  | evc_2 |   | evc_1 |
    !
    ! 
    ! if ipol=2 output is
    !
    !   |  0    -i  |  | evc_1 |   | -i*evc_2 |
    !   |           |  |       | = |          |
    !   |  i     0  |  | evc_2 |   |  i*evc_1 |
    !
    ! 
    ! if ipol=3 output is
    !
    !   |  1     0  |  | evc_1 |   |  evc_1 |
    !   |           |  |       | = |        |
    !   |  0    -1  |  | evc_2 |   | -evc_2 |
    !
    ! 
    ! if ipol=4 output is
    !
    !   |  0     1  |  | evc_1 |   |  evc_2 |
    !   |           |  |       | = |        |
    !   | -1     0  |  | evc_2 |   | -evc_1 |
    !
    USE kinds,               ONLY : DP
    USE wvfct,               ONLY : npwx
    USE noncollin_module,    ONLY : noncolin, npol
    USE io_global,            ONLY : stdout  

    IMPLICIT NONE
    
    complex(DP), intent(inout)  ::  evc_g (npwx*npol) 
    integer, intent(in)         ::  ipol
   
    complex(DP)                 ::  psic (npwx) 
    integer                     ::  ig
  
    IF (.not. noncolin) &
        & call errore('lr_magnons_rountines: pauli','Pauli matrices can be used only&
        &              in non collinear calculations.',1)
  
    ! sigma_x

    IF ( ipol == 1 ) THEN
       !
       psic(:) = evc_g(1:npwx)
       evc_g(1:npwx) =  evc_g(npwx+1:npol*npwx)
       evc_g(npwx+1:npol*npwx) = psic(:)
       !
    ! sigma_y
    ELSEIF ( ipol == 2 ) THEN
       !
       psic(:) = evc_g(1:npwx)                  
       evc_g(1:npwx) =  - evc_g(npwx+1:npol*npwx)
       evc_g(npwx+1:npol*npwx) = psic(:)
       evc_g(:) = (0.0d0,1.0d0)*evc_g(:)
       !
    ! sigma_z
    ELSEIF ( ipol == 3 ) THEN
       !
       evc_g(npwx+1:npol*npwx) = - evc_g(npwx+1:npol*npwx)
       !
    ! (i sigma_y)
    ELSEIF ( ipol == 4 ) THEN
       !
       psic(:) = evc_g(1:npwx)                  
       evc_g(1:npwx) =  evc_g(npwx+1:npol*npwx)
       evc_g(npwx+1:npol*npwx) = -psic(:)
       !
    ELSE
       call errore('lr_magnons_routines: pauli','ipol must be either 1, 2, 3 or 4.',1) 
    ENDIF
    !
    RETURN
    !
  END SUBROUTINE pauli
  !
  !---------------------------------------------------------------------------
  SUBROUTINE T_rev( evc_g, npw_k, ig_k, npw_mk, ig_mk, c_evc_g, spinorial )
    !---------------------------------------------------------------------------
    !
    ! Input:    wfc in G-space evc_g
    !
    ! Performs complex conjugation in real space --> evc_r^*, then Fourier
    ! transforms back giving the 
    !
    ! Output: c_evc_g, wfc in G-space
    !
    ! If spinorial (and noncolin) = .true., it also applies (i sigma_y) 
    !
    !
    ! Notice that the Fourier transform back to G-space is done with the -k
    ! cutoff, since u_k^* has Fourier components such that |-k+G|^2 < e_cut
    !
    USE kinds,                   ONLY : DP
    USE wvfct,                   ONLY : npwx
    USE fft_base,                ONLY : dffts
    USE fft_interfaces,          ONLY : fwfft, invfft
    USE mp_global,               ONLY : inter_pool_comm, intra_bgrp_comm
    USE noncollin_module,        ONLY : noncolin, npol
    USE io_global,               ONLY : stdout
    USE mp,                      ONLY : mp_sum
    !
    IMPLICIT NONE

    !-------------------------INPUT-----------------------------
    ! 
    ! input wavefunction u_k(G)
    !
    COMPLEX(DP), INTENT(IN)  ::   evc_g(npol*npwx)
    !
    ! number of plane-waves of u_k and u_{-k}
    !
    INTEGER, INTENT(IN)      ::   npw_k, npw_mk
    !
    ! mapping G--> k+G and G--> -k+G, output of gk_sort
    !
    INTEGER, INTENT(IN)      ::   ig_k(npwx), ig_mk(npwx)
    ! 
    ! if .true. applies also (i sigma_y) to evc_g
    !
    LOGICAL,INTENT(IN)       ::   spinorial
   
    !-------------------------OUTPUT----------------------------
    !
    !             F[u_k^*](G)    if spinorial = .false.
    ! (i sigma_y) F[u_k^*](G)    if spinorial = .true.
    !
    COMPLEX(DP), INTENT(OUT) ::   c_evc_g(npol*npwx)
    !
    !-------------------INTERNAL-VARIABLES----------------------
    !
    ! u_k(r)
    !
    COMPLEX(DP)              ::   evc_r (dffts%nnr,npol)
    !
    ! counter over Gs
    !
    INTEGER                  ::   ig
    !
    !! COMPLEX(DP)              ::   norm
   

    evc_r   = (0.d0, 0.d0)
    c_evc_g = (0.d0, 0.d0)
    ! To real space
    do ig = 1, npw_k
       evc_r (dffts%nl (ig_k (ig) ),1 ) = evc_g (ig)
    enddo
    CALL invfft ('Wave', evc_r(:,1), dffts)
    IF (noncolin) THEN
       DO ig = 1, npw_k
         evc_r (dffts%nl(ig_k(ig)),2) = evc_g (ig+npwx)
       ENDDO
       CALL invfft ('Wave', evc_r(:,2), dffts)
    ENDIF
 
    ! Complex conjg
    evc_r = conjg(evc_r)

   
    ! Back to G-space
    CALL fwfft ('Wave', evc_r(:,1), dffts)
    do ig = 1, npw_mk
       c_evc_g (ig) =  evc_r (dffts%nl (ig_mk (ig) ), 1 )
    enddo
    IF (noncolin) THEN
       CALL fwfft ('Wave', evc_r(:,2), dffts)
       DO ig = 1, npw_mk
          c_evc_g (ig+npwx) =  evc_r (dffts%nl(ig_mk(ig)),2)
       ENDDO
    ENDIF
    !
    ! Apply (i sigma_y) if spinorial and noncolin = .true.
    !
    if (spinorial) then
       if (noncolin) then
          call pauli( c_evc_g, 4)
       else
          call errore('T-rev','in the non-spinorial case T_rev is &
                       & only complex conjugation',1)
       endif
    endif
    !

    RETURN
    !
  END SUBROUTINE T_rev
  !
 END MODULE lr_magnons_routines
