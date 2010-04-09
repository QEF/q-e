!
! Copyright (C) 2008 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!-----------------------------------------------------------------------
! Modified by Osman Baris Malcioglu (2009)
SUBROUTINE lr_ortho(dvpsi, evq, ikk, ikq, sevc, inverse)
!subroutine lr_ortho(sv)
  !
  !
  ! This routine ortogonalizes dvpsi to the valence states: ps = <evq|dvpsi>
  ! It should be quite general. It works for metals and insulators, with 
  ! NC as well as with US PP, both SR or FR.
  ! Note that on output it changes sign. So it applies -P^+_c.
  !
  !OBM!! evc0 ->evq sevc0 -> sevc dvpsi -> input/output
  !
USE kinds, ONLY : DP
use gvect,                only : gstart
USE klist, ONLY : lgauss, degauss, ngauss
USE noncollin_module, ONLY : noncolin, npol
USE wvfct, ONLY : npwx, nbnd, et
USE ener, ONLY : ef
!USE qpoint, ONLY : npwq
USE control_ph,  ONLY : alpha_pv, nbnd_occ
!USE becmod,      ONLY : becp, becp_nc, calbec
USE uspp,        ONLY : vkb, okvan
USE mp_global,   ONLY : intra_pool_comm
USE mp,          ONLY : mp_sum
!use lr_variables, ONLY : lr_alpha_pv, nbnd_occ, 
use lr_variables, ONLY : lr_verbosity
use realus,       ONLY : npw_k
use control_flags,        only : gamma_only
USE io_global,      ONLY : stdout
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: ikk, ikq   ! the index of the k and k+q points
COMPLEX(DP), INTENT(IN) :: evq(npwx*npol,nbnd)
COMPLEX(DP), INTENT(INOUT) :: dvpsi(npwx*npol,nbnd)
COMPLEX(DP), INTENT(IN) :: sevc(npwx*npol,nbnd) ! work space allocated by
                                                   ! the calling routine (was called dpsi)
!real(kind=dp), intent(IN) :: lr_alpha_pv !This is calculated manually in tddfpt
logical, intent(in):: inverse !if .true. |dvspi> =  |dvpsi> - |evq><sevc|dvpsi>  instead of |dvspi> =  |dvpsi> - |sevc><evq|dvpsi> 


logical::  inverse_mode

! functions computing the delta and theta function

CALL start_clock ('lr_ortho')

If (lr_verbosity > 5)   WRITE(stdout,'("<lr_ortho>")')
  !
  !if (.not. present(inverse)) then 
  !  inverse_mode=.false.
  !else
    inverse_mode=inverse
  !endif
  if (gamma_only) then
     !
     call lr_ortho_gamma()
     !
  else if (noncolin) then
    !
    call lr_ortho_noncolin()
    !
  else
     !
     call lr_ortho_k()
     !
  end if

CALL stop_clock ('lr_ortho')
RETURN
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!multiple K point specific
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine lr_ortho_k()
  IMPLICIT NONE
  COMPLEX(DP), ALLOCATABLE :: ps(:,:)
  INTEGER :: ibnd, jbnd, nbnd_eff
  REAL(DP) :: wg1, w0g, wgp, wwg, deltae, theta
  REAL(DP), EXTERNAL :: w0gauss, wgauss 
  
   ALLOCATE(ps(nbnd,nbnd))
   !
   if (lgauss) then
       !
       !  metallic case
       !
       ps = (0.d0, 0.d0)
       if (inverse_mode) then
       CALL ZGEMM( 'C', 'N', nbnd, nbnd_occ (ikk), npw_k(ikk), (1.d0,0.d0), &
                     sevc, npwx, dvpsi, npwx, (0.d0,0.d0), ps, nbnd )
       else
       CALL ZGEMM( 'C', 'N', nbnd, nbnd_occ (ikk), npw_k(ikk), (1.d0,0.d0), &
                     evq, npwx, dvpsi, npwx, (0.d0,0.d0), ps, nbnd )
       endif
       !
       DO ibnd = 1, nbnd_occ (ikk)
          wg1 = wgauss ((ef-et(ibnd,ikk)) / degauss, ngauss)
          w0g = w0gauss((ef-et(ibnd,ikk)) / degauss, ngauss) / degauss
          DO jbnd = 1, nbnd
             wgp = wgauss ( (ef - et (jbnd, ikq) ) / degauss, ngauss)
             deltae = et (jbnd, ikq) - et (ibnd, ikk)
             theta = wgauss (deltae / degauss, 0)
             wwg = wg1 * (1.d0 - theta) + wgp * theta
             IF (jbnd <= nbnd_occ (ikq) ) THEN
                IF (abs (deltae) > 1.0d-5) THEN
                    wwg = wwg + alpha_pv * theta * (wgp - wg1) / deltae
                ELSE
                   !
                   !  if the two energies are too close takes the limit
                   !  of the 0/0 ratio
                   !
                   wwg = wwg - alpha_pv * theta * w0g
                ENDIF
             ENDIF
             !
             ps(jbnd,ibnd) = wwg * ps(jbnd,ibnd)
             !
          ENDDO
          call DSCAL (2*npw_k(ikk), wg1, dvpsi(1,ibnd), 1)
       END DO
       nbnd_eff=nbnd
   ELSE
      !
      !  insulators
      !
      ps = (0.d0, 0.d0)
      !OBM!!!
      !ps = <evq|dvpsi>
      ! in the old version it was <sevc|dvpsi> 
      if (inverse_mode) then
      CALL ZGEMM( 'C', 'N', nbnd_occ(ikq), nbnd_occ (ikk), npw_k(ikk), &
                (1.d0,0.d0), sevc, npwx, dvpsi, npwx, &
                (0.d0,0.d0), ps, nbnd )
      else
      CALL ZGEMM( 'C', 'N', nbnd_occ(ikq), nbnd_occ (ikk), npw_k(ikk), &
                (1.d0,0.d0), evq, npwx, dvpsi, npwx, &
                (0.d0,0.d0), ps, nbnd )
      endif
      nbnd_eff=nbnd_occ(ikk)
   END IF
#ifdef __PARA
   call mp_sum(ps(:,1:nbnd_eff),intra_pool_comm)
#endif
   !!
   !! |dvspi> =  -(|dvpsi> - |sevc><evq|dvpsi>)
   !!
   !OBM!!! changed to |dvspi> =  |dvpsi> - |sevc><evq|dvpsi> 
   if (lgauss) then
      !
      !  metallic case
      !
      if (inverse_mode) then
      CALL ZGEMM( 'N', 'N', npw_k(ikk), nbnd_occ(ikk), nbnd, &
                (-1.d0,0.d0), evq, npwx, ps, nbnd, (1.0d0,0.d0), &
                 dvpsi, npwx )
      else
      CALL ZGEMM( 'N', 'N', npw_k(ikk), nbnd_occ(ikk), nbnd, &
                (-1.d0,0.d0), sevc, npwx, ps, nbnd, (1.0d0,0.d0), &
                 dvpsi, npwx )
      endif
   ELSE
      !
      !  Insulators: note that nbnd_occ(ikk)=nbnd_occ(ikq) in an insulator
      !
      if (inverse_mode) then
      CALL ZGEMM( 'N', 'N', npw_k(ikk), nbnd_occ(ikk), nbnd_occ(ikk), &
                (-1.d0,0.d0), evq, npwx, ps, nbnd, (1.0d0,0.d0), &
                 dvpsi, npwx )
      else
      CALL ZGEMM( 'N', 'N', npw_k(ikk), nbnd_occ(ikk), nbnd_occ(ikk), &
                (-1.d0,0.d0), sevc, npwx, ps, nbnd, (1.0d0,0.d0), &
                 dvpsi, npwx )
      endif
   ENDIF
   DEALLOCATE(ps)
  end subroutine lr_ortho_k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!Gamma point specific
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine lr_ortho_gamma()
  IMPLICIT NONE
  COMPLEX(DP), ALLOCATABLE :: ps_c(:,:)
  REAL(DP), ALLOCATABLE :: ps(:,:)
  INTEGER :: ibnd, jbnd, nbnd_eff
  REAL(DP) :: wg1, w0g, wgp, wwg, deltae, theta
  REAL(DP), EXTERNAL :: w0gauss, wgauss 
  
   ALLOCATE(ps(nbnd,nbnd))
   ALLOCATE(ps_c(nbnd,nbnd))
   !
   if (lgauss) then
       call errore ('lr_ortho', "degauss with gamma point algorithms",1)
   ELSE
      !
      !  insulators
      !  ps = <evq|dvpsi>
      ! in old version it was ps = <S evc0|sv> 
     ps = 0.d0
     if (inverse_mode) then
     CALL DGEMM( 'C', 'N', nbnd, nbnd ,2*npw_k(1), &
               2.d0, sevc, 2*npwx, dvpsi, 2*npwx, &
               0.d0, ps, nbnd )                          
               !ps = 2*<sevc|dvpsi>
     else
     CALL DGEMM( 'C', 'N', nbnd, nbnd ,2*npw_k(1), &
               2.d0, evq, 2*npwx, dvpsi, 2*npwx, &
               0.d0, ps, nbnd )
               !ps = 2*<evq|dvpsi>
     endif
     nbnd_eff=nbnd
     if (gstart == 2) then
     if (inverse_mode) then
         CALL DGER( nbnd, nbnd, -1.D0, sevc, 2*npwx, dvpsi, 2*npwx, ps, nbnd )
         !PS = PS - sevc*dvpsi
     else
         CALL DGER( nbnd, nbnd, -1.D0, evq, 2*npwx, dvpsi, 2*npwx, ps, nbnd )
         !PS = PS - evc*dvpsi
     endif
     endif
  END IF
#ifdef __PARA
   call mp_sum(ps(:,:),intra_pool_comm)
#endif
  ! in the original dpsi was used as a storage for sevc, since in
  ! tddfpt we have it stored in memory as sevc0 this part is obsolote
  !!
  !! dpsi is used as work space to store S|evc>
  !!
  !IF (noncolin) THEN
  !   IF (okvan) CALL calbec ( npw_k(ikk), vkb, evq, becp_nc, nbnd_eff )
  !ELSE
  !   IF (okvan) CALL calbec ( npwq, vkb, evq, becp, nbnd_eff)
  !ENDIF
  !CALL s_psi (npwx, npwq, nbnd_eff, evq, dpsi)


  ps_c = cmplx(ps, 0.d0, dp) 
  !!
  !! |dvspi> =  -(|dvpsi> - S|evq><evq|dvpsi>)
  !!
  !OBM!!! changed to |dvspi> =  |dvpsi>  - |sevc><evq|dvpsi> 
  if (lgauss) then
          !errore ? 
  ELSE
     !
     !  Insulators: note that nbnd_occ(ikk)=nbnd_occ(ikq) in an insulator
     !
     if (inverse_mode) then
     CALL ZGEMM( 'N', 'N', npw_k(1), nbnd, nbnd, &
               (-1.d0,0.d0), evq, npwx, ps_c, nbnd, (1.0d0,0.d0), &
                dvpsi, npwx )
                !dvpsi=dvpsi-|evq><sevc|dvpsi>
     else
     CALL ZGEMM( 'N', 'N', npw_k(1), nbnd, nbnd, &
               (-1.d0,0.d0), sevc, npwx, ps_c, nbnd, (1.0d0,0.d0), &
                dvpsi, npwx )
                !dvpsi=dvpsi-|sevc><evq|dvpsi>
     endif
  ENDIF
   DEALLOCATE(ps)
   DEALLOCATE(ps_c)
  end subroutine lr_ortho_gamma
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!noncolin specific
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine lr_ortho_noncolin()
  IMPLICIT NONE
  COMPLEX(DP), ALLOCATABLE :: ps(:,:)
  INTEGER :: ibnd, jbnd, nbnd_eff
  REAL(DP) :: wg1, w0g, wgp, wwg, deltae, theta
  REAL(DP), EXTERNAL :: w0gauss, wgauss 
  
   ALLOCATE(ps(nbnd,nbnd))
   !
   if (lgauss) then
       !
       !  metallic case
       !
       ps = (0.d0, 0.d0)
          CALL ZGEMM( 'C', 'N', nbnd, nbnd_occ (ikk), npwx*npol, (1.d0,0.d0), &
                     evq, npwx*npol, dvpsi, npwx*npol, (0.d0,0.d0), ps, nbnd )
       !
       DO ibnd = 1, nbnd_occ (ikk)
          wg1 = wgauss ((ef-et(ibnd,ikk)) / degauss, ngauss)
          w0g = w0gauss((ef-et(ibnd,ikk)) / degauss, ngauss) / degauss
          DO jbnd = 1, nbnd
             wgp = wgauss ( (ef - et (jbnd, ikq) ) / degauss, ngauss)
             deltae = et (jbnd, ikq) - et (ibnd, ikk)
             theta = wgauss (deltae / degauss, 0)
             wwg = wg1 * (1.d0 - theta) + wgp * theta
             IF (jbnd <= nbnd_occ (ikq) ) THEN
                IF (abs (deltae) > 1.0d-5) THEN
                    wwg = wwg + alpha_pv * theta * (wgp - wg1) / deltae
                ELSE
                   !
                   !  if the two energies are too close takes the limit
                   !  of the 0/0 ratio
                   !
                   wwg = wwg - alpha_pv * theta * w0g
                ENDIF
             ENDIF
             !
             ps(jbnd,ibnd) = wwg * ps(jbnd,ibnd)
             !
          ENDDO
             CALL DSCAL (2*npwx*npol, wg1, dvpsi(1,ibnd), 1)
       END DO
       nbnd_eff=nbnd
   ELSE
      !
      !  insulators
      !
      ps = (0.d0, 0.d0)
         CALL ZGEMM( 'C', 'N',nbnd_occ(ikq), nbnd_occ(ikk), npwx*npol, &
                (1.d0,0.d0), evq, npwx*npol, dvpsi, npwx*npol, &
                (0.d0,0.d0), ps, nbnd )
      nbnd_eff=nbnd_occ(ikk)
   END IF
#ifdef __PARA
   call mp_sum(ps(:,1:nbnd_eff),intra_pool_comm)
#endif
   ! in the original dpsi was used as a storage for sevc, since in
   ! tddfpt we have it stored in memory as sevc0 this part is obsolote
   !!
   !! dpsi is used as work space to store S|evc>
   !!
   !IF (noncolin) THEN
   !   IF (okvan) CALL calbec ( npw_k(ikk), vkb, evq, becp_nc, nbnd_eff )
   !ELSE
   !   IF (okvan) CALL calbec ( npwq, vkb, evq, becp, nbnd_eff)
   !ENDIF
   !CALL s_psi (npwx, npwq, nbnd_eff, evq, dpsi)
   !!
   !! |dvspi> =  -(|dvpsi> - S|evq><evq|dvpsi>)
   !OBM!!! changed to |dvspi> =  |dvpsi> - S|evq><evq|dvpsi> using this
   !!
   if (lgauss) then
      !
      !  metallic case
      !
         CALL ZGEMM( 'N', 'N', npwx*npol, nbnd_occ(ikk), nbnd, &
                   (-1.d0,0.d0), sevc, npwx*npol, ps, nbnd, (1.0d0,0.d0), &
                   dvpsi, npwx*npol )
   ELSE
      !
      !  Insulators: note that nbnd_occ(ikk)=nbnd_occ(ikq) in an insulator
      !
         CALL ZGEMM( 'N', 'N', npwx*npol, nbnd_occ(ikk), nbnd_occ(ikk), &
                   (-1.d0,0.d0),sevc,npwx*npol,ps,nbnd,(1.0d0,0.d0), &
                   dvpsi, npwx*npol )
   ENDIF
   DEALLOCATE(ps)
  end subroutine lr_ortho_noncolin


end subroutine lr_ortho
!-----------------------------------------------------------------------
