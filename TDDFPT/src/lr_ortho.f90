!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE lr_ortho(dvpsi, evq, ikk, ikq, sevc, inverse)
  !---------------------------------------------------------------------
  !
  ! Inspired by LR_Modules/orthogonalize.f90
  ! This routine orthogonalizes dvpsi to the valence states.
  ! It should be quite general. It works for metals and insulators, with 
  ! NC as well as with US PP, both SR or FR.
  !
  ! If inverse=.true.  then apply P_c   = 1 - |evq><sevc| = 1 -  |psi><psi|S
  ! If inverse=.false. then apply P_c^+ = 1 - |sevc><evq| = 1 - S|psi><psi| 
  !
  ! See definitions of P_c and P_c^+ eq.(29) in 
  ! B. Walker and R. Gebauer, J. Chem. Phys. 127, 164106 (2007)
  !
  ! NB: IN/OUT is dvpsi ; sevc is used as a workspace
  !
  ! evc0 -> evq 
  ! sevc0 -> sevc
  !
  ! Note: S^{-1} P_c^+(k) = P_c(k) S^{-1}
  !
  ! This subroutine is a modified version of the the subroutine
  ! orthogonalize.f90. This subroutine should be removed in the future,
  ! and orthogonalize.f90 should be used instead. Currently this 
  ! subroutine is used only in several places in lr_dav_routines.f90.
  ! TODO: Check if lr_ortho can be replaced by orthogonalize in
  ! lr_dav_routines.
  !
  ! Modified by Osman Baris Malcioglu (2009)
  ! Modified by Iurii Timrov (2013)
  !
  USE kinds,             ONLY : DP
  use gvect,             only : gstart
  USE klist,             ONLY : ngk, lgauss, degauss, ngauss
  USE noncollin_module,  ONLY : noncolin, npol
  USE wvfct,             ONLY : npwx, nbnd, et
  USE ener,              ONLY : ef
  USE control_lr,        ONLY : alpha_pv, nbnd_occ
  USE uspp,              ONLY : vkb, okvan
  USE mp_global,         ONLY : intra_bgrp_comm
  USE mp,                ONLY : mp_sum
  use lr_variables,      ONLY : lr_verbosity
  use control_flags,     only : gamma_only
  USE io_global,         ONLY : stdout
#if defined(__CUDA)
  USE cublas
#endif  
  !
  IMPLICIT NONE
  INTEGER, INTENT(in) :: ikk, ikq   
  ! the index of the k and k+q points
  ! (in the optical case ikq = ikk)
  COMPLEX(DP), INTENT(in) ::    evq(npwx*npol,nbnd)
  COMPLEX(DP), INTENT(inout) :: dvpsi(npwx*npol,nbnd)
  COMPLEX(DP), INTENT(in) ::    sevc(npwx*npol,nbnd)
  ! sevc is a work space allocated by the calling routine (was called dpsi)
  LOGICAL, INTENT(in):: inverse 
  LOGICAL::  inverse_mode
  !
  CALL start_clock ('lr_ortho')
  !
  IF (lr_verbosity > 5)   WRITE(stdout,'("<lr_ortho>")')
  !
  !if (.not. present(inverse)) then
  !  inverse_mode=.false.
  !else
    inverse_mode = inverse
  !endif
  !$acc data present_or_copyin(evq(1:npwx*npol,1:nbnd), sevc(1:npwx*npol,1:nbnd)) present_or_copy(dvpsi(1:npwx*npol,1:nbnd))
  !
  IF (gamma_only) THEN
     !
     CALL lr_ortho_gamma()
     !
  ELSEIF (noncolin) THEN
    !
    CALL lr_ortho_noncolin()
    !
  ELSE
     !
     CALL lr_ortho_k()
     !
  ENDIF
  !
  !$acc end data
  !
  CALL stop_clock ('lr_ortho')
  !
  RETURN
  !
CONTAINS

SUBROUTINE lr_ortho_k()
  !
  ! General K points algorithm
  !
  IMPLICIT NONE

  COMPLEX(DP), ALLOCATABLE :: ps(:,:)
  INTEGER :: ibnd, jbnd, nbnd_eff
  REAL(DP) :: wg1, w0g, wgp, wwg(nbnd), deltae, theta
  REAL(DP), EXTERNAL :: w0gauss, wgauss

  ALLOCATE(ps(nbnd,nbnd))
  !
  !$acc enter data create(ps(:,:))
  !
  IF (lgauss) THEN
     !
     !  metallic case
     !
     !$acc kernels
     ps(:,:) = (0.d0, 0.d0)
     !$acc end kernels
     !
     !$acc host_data use_device(sevc, evq, dvpsi, ps)
     !
     IF (inverse_mode) THEN
        CALL ZGEMM( 'C', 'N', nbnd, nbnd_occ(ikk), ngk(ikk), (1.d0,0.d0), &
                    sevc, npwx, dvpsi, npwx, (0.d0,0.d0), ps, nbnd )
     ELSE
        CALL ZGEMM( 'C', 'N', nbnd, nbnd_occ(ikk), ngk(ikk), (1.d0,0.d0), &
                    evq, npwx, dvpsi, npwx, (0.d0,0.d0), ps, nbnd )
     ENDIF
     !
     !$acc end host_data
     !
     DO ibnd = 1, nbnd_occ (ikk)
        wg1 = wgauss ((ef-et(ibnd,ikk)) / degauss, ngauss)
        w0g = w0gauss((ef-et(ibnd,ikk)) / degauss, ngauss) / degauss
        DO jbnd = 1, nbnd
           wgp = wgauss ( (ef - et (jbnd, ikq) ) / degauss, ngauss)
           deltae = et (jbnd, ikq) - et (ibnd, ikk)
           theta = wgauss (deltae / degauss, 0)
           wwg(jbnd) = wg1 * (1.d0 - theta) + wgp * theta
           IF (jbnd <= nbnd_occ (ikq) ) THEN
              IF (abs (deltae) > 1.0d-5) THEN
                 wwg(jbnd) = wwg(jbnd) + alpha_pv * theta * (wgp - wg1) / deltae
              ELSE
                 !
                 !  if the two energies are too close takes the limit
                 !  of the 0/0 ratio
                 !
                 wwg(jbnd) = wwg(jbnd) - alpha_pv * theta * w0g
              ENDIF
           ENDIF
           !
        ENDDO
        !
        !$acc parallel loop copyin(wwg)
        DO jbnd = 1, nbnd
           ps(jbnd,ibnd) = wwg(jbnd) * ps(jbnd,ibnd)
        ENDDO
        !$acc kernels
        dvpsi(1:ngk(ikk), ibnd) = wg1 * dvpsi(1:ngk(ikk), ibnd)
        !$acc end kernels
     ENDDO
     !
     nbnd_eff = nbnd
     !
  ELSE
     !
     !  insulators
     !
     !$acc kernels
     ps(:,:) = (0.d0, 0.d0)
     !$acc end kernels
     !
     !$acc host_data use_device(sevc, evq, dvpsi, ps)
     !
     IF (inverse_mode) THEN
        !
        ! ps = <sevc|dvpsi>
        !
        CALL ZGEMM( 'C', 'N', nbnd_occ(ikq), nbnd_occ (ikk), ngk(ikk), &
                  (1.d0,0.d0), sevc, npwx, dvpsi, npwx, &
                  (0.d0,0.d0), ps, nbnd )
        !
     ELSE
        !
        ! ps = <evq|dvpsi>
        !
        CALL ZGEMM( 'C', 'N', nbnd_occ(ikq), nbnd_occ (ikk), ngk(ikk), &
                  (1.d0,0.d0), evq, npwx, dvpsi, npwx, &
                  (0.d0,0.d0), ps, nbnd )
     ENDIF
     !
     !$acc end host_data
     !
     nbnd_eff = nbnd_occ(ikk)
     !
  ENDIF
  !
#if defined(__MPI)
  !$acc host_data use_device(ps) 
  CALL mp_sum(ps(:,1:nbnd_eff),intra_bgrp_comm)
  !$acc end host_data
#endif
  !
  !$acc host_data use_device(sevc, evq, dvpsi, ps)
  !
  IF (lgauss) THEN
     !
     !  Metallic case
     !
     if (inverse_mode) then
        CALL ZGEMM( 'N', 'N', ngk(ikk), nbnd_occ(ikk), nbnd, &
                  (-1.d0,0.d0), evq, npwx, ps, nbnd, (1.0d0,0.d0), &
                  dvpsi, npwx )
     else
        CALL ZGEMM( 'N', 'N', ngk(ikk), nbnd_occ(ikk), nbnd, &
                  (-1.d0,0.d0), sevc, npwx, ps, nbnd, (1.0d0,0.d0), &
                  dvpsi, npwx )
     endif
     !
  ELSE
     !
     !  Insulators: note that nbnd_occ(ikk)=nbnd_occ(ikq) in an insulator
     !
     if (inverse_mode) then
        !
        ! |dvspi> =  |dvpsi> - |evq><sevc|dvpsi>
        !
        CALL ZGEMM( 'N', 'N', ngk(ikk), nbnd_occ(ikk), nbnd_occ(ikk), &
                  (-1.d0,0.d0), evq, npwx, ps, nbnd, (1.0d0,0.d0), &
                  dvpsi, npwx )
        !
     else
        !
        ! |dvspi> =  |dvpsi> - |sevc><evq|dvpsi>
        !
        CALL ZGEMM( 'N', 'N', ngk(ikk), nbnd_occ(ikk), nbnd_occ(ikk), &
                  (-1.d0,0.d0), sevc, npwx, ps, nbnd, (1.0d0,0.d0), &
                  dvpsi, npwx )
        !
     endif
     !
  ENDIF
  !
  !$acc end host_data
  !
  !$acc exit data delete(ps)
  !
  DEALLOCATE(ps)
  !
END SUBROUTINE lr_ortho_k

SUBROUTINE lr_ortho_gamma()
  !
  ! gamma_only algorithm
  !
  IMPLICIT NONE

  COMPLEX(DP), ALLOCATABLE :: ps_c(:,:)
  REAL(DP), ALLOCATABLE :: ps(:,:)
  INTEGER :: nbnd_eff

  ALLOCATE(ps(nbnd,nbnd))
  ALLOCATE(ps_c(nbnd,nbnd))
  !
  !$acc enter data create(ps(:,:), ps_c(:,:))
  !
  IF (lgauss) THEN
      !
      ! Metals
      !
      CALL errore ('lr_ortho', "degauss with gamma point algorithm is not allowed",1)
      !
  ELSE
      !
      ! Insulators
      !
      !$acc kernels
      ps(:,:) = 0.d0
      !$acc end kernels
      !
      !$acc host_data use_device(sevc,evq,dvpsi,ps)
      !
      IF (inverse_mode) THEN
         ! 
         ! ps = 2 * <sevc|dvpsi>
         !
         CALL DGEMM( 'C', 'N', nbnd, nbnd ,2*ngk(1), &
               2.d0, sevc, 2*npwx, dvpsi, 2*npwx, 0.d0, ps, nbnd )
         !
      ELSE
         !
         ! ps = 2 * <evq|dvpsi>
         !
         CALL DGEMM( 'C', 'N', nbnd, nbnd ,2*ngk(1), &
               2.d0, evq, 2*npwx, dvpsi, 2*npwx, 0.d0, ps, nbnd )
         !
      ENDIF
      !
      !$acc end host_data
      !
      nbnd_eff = nbnd
      !
      ! G=0 has been accounted twice, therefore we subtruct one contribution.
      !
      !$acc host_data use_device(sevc,evq,dvpsi,ps)
      !
      IF (gstart == 2) THEN
         !
         IF (inverse_mode) THEN
            !
            ! ps = ps - <sevc|dvpsi>
            !
            CALL MYDGER( nbnd, nbnd, -1.D0, sevc, 2*npwx, dvpsi, 2*npwx, ps, nbnd)
            !
         ELSE
            !
            ! ps = ps - <evq|dvpsi>
            !
            CALL MYDGER( nbnd, nbnd, -1.D0, evq, 2*npwx, dvpsi, 2*npwx, ps, nbnd )
            !
         ENDIF
         !
      ENDIF
      !
      !$acc end host_data
      !
  ENDIF
  !
#if defined(__MPI)
  !$acc host_data use_device(ps)
  CALL mp_sum(ps(:,:),intra_bgrp_comm)
  !$acc end host_data
#endif
  !
  ! In the original code dpsi was used as a storage for sevc, since in
  ! tddfpt we have it stored in memory as sevc0 this part is obsolote
  !
  !$acc kernels
  ps_c(:,:) = cmplx(ps(:,:), 0.d0, dp)
  !$acc end kernels
  !
  IF (lgauss) THEN
      !
      ! Metals
      !
      CALL errore ('lr_ortho', "degauss with gamma point algorithm is not allowed",1)
      !
  ELSE
     !
     !  Insulators: note that nbnd_occ(ikk) = nbnd_occ(ikq)
     !
     !$acc host_data use_device(sevc,evq,dvpsi,ps_c)
     !
     IF (inverse_mode) THEN
        !
        ! |dvpsi> = |dvpsi> - |evq><sevc|dvpsi>
        !
        CALL ZGEMM( 'N', 'N', ngk(1), nbnd, nbnd, &
                  (-1.d0,0.d0), evq, npwx, ps_c, nbnd, (1.0d0,0.d0), dvpsi, npwx)
        !
     ELSE
        !
        ! |dvpsi> = |dvpsi> - |sevc><evq|dvpsi>
        !
        CALL ZGEMM( 'N', 'N', ngk(1), nbnd, nbnd, &
                  (-1.d0,0.d0), sevc, npwx, ps_c, nbnd, (1.0d0,0.d0), dvpsi, npwx )
        !
     ENDIF
     !
     !$acc end host_data
     !
  ENDIF
  !
  !$acc exit data delete(ps, ps_c)
  !
  DEALLOCATE(ps)
  DEALLOCATE(ps_c)
  !
  RETURN
  !
END SUBROUTINE lr_ortho_gamma

SUBROUTINE lr_ortho_noncolin()
  !
  ! Noncollinear case 
  ! Warning: the inverse mode is not implemented!
  !
  IMPLICIT NONE

  COMPLEX(DP), ALLOCATABLE :: ps(:,:)
  INTEGER :: ibnd, jbnd, nbnd_eff
  REAL(DP) :: wg1, w0g, wgp, wwg(nbnd), deltae, theta
  REAL(DP), EXTERNAL :: w0gauss, wgauss
  
  IF (inverse_mode) CALL errore ('lr_ortho', "The inverse mode is not implemented!",1)
  !
  ALLOCATE(ps(nbnd,nbnd))
  !
  !$acc enter data create(ps(:,:))
             
  !
  IF (lgauss) THEN
       !
       !  metallic case
       !
       !$acc kernels
       ps(:,:) = (0.d0, 0.d0)
       !$acc end kernels
       !
       !$acc host_data use_device(evq,dvpsi,ps)
       CALL ZGEMM( 'C', 'N', nbnd, nbnd_occ (ikk), npwx*npol, (1.d0,0.d0), &
                    evq, npwx*npol, dvpsi, npwx*npol, (0.d0,0.d0), ps, nbnd )
       !$acc end host_data
       !
       DO ibnd = 1, nbnd_occ (ikk)
          !
          wg1 = wgauss ((ef-et(ibnd,ikk)) / degauss, ngauss)
          w0g = w0gauss((ef-et(ibnd,ikk)) / degauss, ngauss) / degauss
          !
          DO jbnd = 1, nbnd
             !
             wgp = wgauss ( (ef - et (jbnd, ikq) ) / degauss, ngauss)
             deltae = et (jbnd, ikq) - et (ibnd, ikk)
             theta = wgauss (deltae / degauss, 0)
             wwg = wg1 * (1.d0 - theta) + wgp * theta
             !
             IF (jbnd <= nbnd_occ (ikq) ) THEN
                IF (abs (deltae) > 1.0d-5) THEN
                    wwg(jbnd) = wwg(jbnd) + alpha_pv * theta * (wgp - wg1) / deltae
                ELSE
                   !
                   !  if the two energies are too close takes the limit
                   !  of the 0/0 ratio
                   !
                   wwg(jbnd) = wwg(jbnd) - alpha_pv * theta * w0g
                ENDIF
             ENDIF
             !
          ENDDO
          !
          !$acc parallel loop copyin(wwg)
          DO jbnd = 1, nbnd
             ps(jbnd,ibnd) = wwg(jbnd) * ps(jbnd,ibnd)
          ENDDO
          !
          !$acc kernels
          dvpsi(1:npwx*npol, ibnd) = wg1 * dvpsi(1:npwx*npol, ibnd)
          !$acc end kernels
       ENDDO
       !
       nbnd_eff = nbnd
       !
   ELSE
      !
      !  insulators
      !
      !$acc kernels
      ps(:,:) = (0.d0, 0.d0)
      !$acc end kernels
      !
      ! ps = <evq|dvpsi>
      !
      !
      !$acc host_data use_device(evq,dvpsi,ps)
      CALL ZGEMM( 'C', 'N',nbnd_occ(ikq), nbnd_occ(ikk), npwx*npol, &
                (1.d0,0.d0), evq, npwx*npol, dvpsi, npwx*npol, (0.d0,0.d0), ps, nbnd )
      !$acc end host_data  
      !
      nbnd_eff = nbnd_occ(ikk)
      !
   ENDIF
   !
#if defined(__MPI)
   !$acc host_data use_device(ps)
   CALL mp_sum(ps(:,1:nbnd_eff),intra_bgrp_comm)
   !$acc end host_data
#endif
   !
   ! In the original dpsi was used as a storage for sevc, since in
   ! tddfpt we have it stored in memory as sevc0 this part is obsolote.
   !
   !$acc host_data use_device(sevc,dvpsi,ps)
   IF (lgauss) THEN
      !
      !  metallic case
      !
      ! |dvpsi> = |dvpsi> - |sevc><evq|dvpsi>
      !
      CALL ZGEMM( 'N', 'N', npwx*npol, nbnd_occ(ikk), nbnd, &
                (-1.d0,0.d0), sevc, npwx*npol, ps, nbnd, (1.0d0,0.d0), dvpsi,npwx*npol )
      !
   ELSE
      !
      ! Insulators: note that nbnd_occ(ikk)=nbnd_occ(ikq)
      !  
      ! |dvpsi> = |dvpsi> - |sevc><evq|dvpsi>
      !
      CALL ZGEMM( 'N', 'N', npwx*npol, nbnd_occ(ikk), nbnd_occ(ikk), &
                (-1.d0,0.d0),sevc,npwx*npol,ps,nbnd,(1.0d0,0.d0), dvpsi,npwx*npol )
      !
   ENDIF
   !$acc end host_data
   !
   !$acc exit data delete(ps)
   !
   DEALLOCATE(ps)
   !
   RETURN
   !   
END SUBROUTINE lr_ortho_noncolin

END SUBROUTINE lr_ortho
!-----------------------------------------------------------------------
