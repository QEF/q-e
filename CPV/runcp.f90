!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

! Written and Revised by Carlo Cavazzoni

!=----------------------------------------------------------------------------------=!


   SUBROUTINE runcp_uspp_x &
      ( nfi, fccc, ccc, ema0bg, dt2bye, rhos, bec, c0, cm, fromscra, restart )
      !
      !  This subroutine performs a Car-Parrinello or Steepest-Descent step
      !  on the electronic variables, computing forces on electrons
      ! 
      !  on input:
      !  c0  wave functions at time t
      !  cm  wave functions at time t - dt 
      !
      !  on output:
      !  cm  wave functions at time t + dt, not yet othogonalized 
      !
      USE parallel_include
      USE kinds,               ONLY : DP
      USE mp_global,           ONLY : nogrp, ogrp_comm, me_image
      USE fft_base,            ONLY : dffts
      use wave_base,           only : wave_steepest, wave_verlet
      use control_flags,       only : lwf, tsde, use_task_groups, program_name
      use uspp,                only : deeq, vkb
      use reciprocal_vectors,  only : gstart
      use electrons_base,      only : n=>nbsp, ispin, f, nspin, nupdwn, iupdwn
      use wannier_subroutines, only : ef_potential
      use efield_module,       only : dforce_efield, tefield, dforce_efield2, tefield2
      USE task_groups,         ONLY : nolist
      use gvecw,               only : ngw
      USE cp_interfaces,       ONLY : dforce
      USE ldaU
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: nfi
      REAL(DP) :: fccc, ccc
      REAL(DP) :: ema0bg(:), dt2bye
      REAL(DP) :: rhos(:,:)
      REAL(DP) :: bec(:,:)
      COMPLEX(DP) :: c0(:,:), cm(:,:)
      LOGICAL, OPTIONAL, INTENT(IN) :: fromscra
      LOGICAL, OPTIONAL, INTENT(IN) :: restart
      !
      !
     real(DP) ::  verl1, verl2, verl3
     real(DP),    allocatable :: emadt2(:)
     real(DP),    allocatable :: emaver(:)
     complex(DP), allocatable :: c2(:), c3(:)
     INTEGER,     ALLOCATABLE :: recv_cnt(:), recv_displ(:)
     REAL(DP),    ALLOCATABLE :: tg_rhos(:,:)
     integer :: i, nsiz, incr, idx, idx_in, ierr
     integer :: iwfc, nwfc, is, ii, tg_rhos_siz, c2_siz
     integer :: iflag
     logical :: ttsde



     iflag = 0
     !
     IF( PRESENT( fromscra ) ) THEN
       IF( fromscra ) iflag = 1
     END IF
     IF( PRESENT( restart ) ) THEN
       IF( restart ) iflag = 2
     END IF

     IF( use_task_groups ) THEN
        tg_rhos_siz = (nogrp+1)*dffts%nr1x * dffts%nr2x*MAXVAL( dffts%npp )
        c2_siz      = (nogrp+1)*ngw 
     ELSE
        tg_rhos_siz = 1
        c2_siz      = ngw 
     END IF

     !
     ! ...  set verlet variables 
     !
     verl1 = 2.0d0 * fccc
     verl2 = 1.0d0 - verl1
     verl3 = 1.0d0 * fccc

     ALLOCATE( emadt2( ngw ) )
     ALLOCATE( emaver( ngw ) )

     ccc    = fccc * dt2bye
     emadt2 = dt2bye * ema0bg
     emaver = emadt2 * verl3

     IF( iflag == 0 ) THEN
       ttsde  = tsde
     ELSE IF( iflag == 1 ) THEN
       ttsde = .TRUE.
     ELSE IF( iflag == 2 ) THEN
       ttsde = .FALSE.
     END IF

     IF( lwf ) THEN

        call ef_potential( nfi, rhos, bec, deeq, vkb, c0, cm, emadt2, emaver, verl1, verl2 )

     ELSE

        allocate( c2( c2_siz ), c3( c2_siz ) )
        allocate( tg_rhos( tg_rhos_siz, nspin ) )

        tg_rhos = 0D0
        c2      = 0D0
        c3      = 0D0

        IF( use_task_groups ) THEN
           !
           !  The potential in rhos is distributed accros all processors
           !  We need to redistribute it so that it is completely contained in the
           !  processors of an orbital TASK-GROUP
           !
           ALLOCATE( recv_cnt( nogrp ), recv_displ( nogrp ) )
           
           recv_cnt(1)   = dffts%npp( nolist(1) + 1 ) * dffts%nr1x * dffts%nr2x
           recv_displ(1) = 0
           DO i = 2, NOGRP
              recv_cnt(i) = dffts%npp( nolist(i) + 1 ) * dffts%nr1x * dffts%nr2x
              recv_displ(i) = recv_displ(i-1) + recv_cnt(i-1)  ! "i-i" _NON_ "i", CAZZOOOO!!!
           ENDDO

#if defined (__PARA) && defined (__MPI)
           DO i = 1, nspin
              nsiz = dffts%npp(me_image+1)* dffts%nr1x * dffts%nr2x
              CALL MPI_Allgatherv( rhos(1,i), nsiz, MPI_DOUBLE_PRECISION, &
                   tg_rhos(1,i), recv_cnt, recv_displ, MPI_DOUBLE_PRECISION, ogrp_comm, IERR)
           ENDDO
#endif
           DEALLOCATE( recv_cnt, recv_displ )

           incr = 2 * nogrp

        ELSE

           incr = 2

        END IF

        

        DO i = 1, n, incr

           IF( use_task_groups ) THEN
              !
              !The input coefficients to dforce cover eigenstates i:i+2*NOGRP-1
              !Thus, in dforce the dummy arguments for c0(1,i) and
              !c0(1,i+1) hold coefficients for eigenstates i,i+2*NOGRP-2,2
              !and i+1,i+2*NOGRP...for example if NOGRP is 4 then we would have
              !1-3-5-7 and 2-4-6-8
              !

              if( tefield .OR. tefield2 ) then
                 CALL errore( ' runcp_uspp ', ' electric field with task group not implemented yet ', 1 )
              end if
              if( lda_plus_u ) then
                 CALL errore( ' runcp_uspp ', ' lda_plus_u with task group not implemented yet ', 1 )
              end if

              CALL dforce( i, bec, vkb, c0, c2, c3, tg_rhos, tg_rhos_siz, ispin, f, n, nspin )

           ELSE

              CALL dforce( i, bec, vkb, c0, c2, c3, rhos, SIZE(rhos,1), ispin, f, n, nspin )

           END IF


           IF ( lda_plus_u ) THEN
              c2(:) = c2(:) - vupsi(:,i)
              c3(:) = c3(:) - vupsi(:,i+1)
           END IF

           IF( tefield ) THEN
             CALL dforce_efield ( bec, i, c0, c2, c3, rhos)
           END IF

           IF( tefield2 ) THEN
             CALL dforce_efield2 ( bec, i, c0, c2, c3, rhos)
           END IF

           IF( iflag == 2 ) THEN
              DO idx = 1, incr, 2
                 IF( i + idx - 1 <= n ) THEN
                    cm( :, i+idx-1) = c0(:,i+idx-1)
                    cm( :, i+idx  ) = c0(:,i+idx  )
                 END IF
              ENDDO
           END IF

           idx_in = 1
           DO idx = 1, incr, 2
              IF( i + idx - 1 <= n ) THEN
                 IF (tsde) THEN
                    CALL wave_steepest( cm(:, i+idx-1 ), c0(:, i+idx-1 ), emaver, c2, ngw, idx_in )
                    CALL wave_steepest( cm(:, i+idx   ), c0(:, i+idx   ), emaver, c3, ngw, idx_in )
                 ELSE
                    CALL wave_verlet( cm(:, i+idx-1 ), c0(:, i+idx-1 ), verl1, verl2, emaver, c2, ngw, idx_in )
                    CALL wave_verlet( cm(:, i+idx   ), c0(:, i+idx   ), verl1, verl2, emaver, c3, ngw, idx_in )
                 ENDIF
                 IF ( gstart == 2 ) THEN
                    cm(1,i+idx-1) = cmplx(real(cm(1,i+idx-1)),0.0d0)
                    cm(1,i+idx  ) = cmplx(real(cm(1,i+idx  )),0.0d0)
                 END IF
              END IF
              !
              idx_in = idx_in + 1
              !
           END DO

        end do

        DEALLOCATE( c2, c3, tg_rhos )

     END IF

     DEALLOCATE( emadt2 )
     DEALLOCATE( emaver )
!
   END SUBROUTINE runcp_uspp_x
!
!
!=----------------------------------------------------------------------------=!
!
!

!=----------------------------------------------------------------------------=!

    SUBROUTINE runcp_uspp_force_pairing_x  &
       ( nfi, fccc, ccc, ema0bg, dt2bye, rhos, bec, c0, cm, intermed, fromscra, &
         restart )
  !
!  same as runcp, except that electrons are paired forcedly
!  i.e. this handles a state dependant Hamiltonian for the paired and unpaired electrons
!  unpaired is assumed to exist, to be unique, and located in highest index band

      USE kinds,               ONLY : DP
      USE wave_base,           ONLY : wave_steepest, wave_verlet
      USE control_flags,       ONLY : lwf, tsde, program_name
      USE uspp,                ONLY : deeq, vkb
      USE reciprocal_vectors,  ONLY : gstart
      USE wannier_subroutines, ONLY : ef_potential
      USE efield_module,       ONLY : dforce_efield, tefield
      USE electrons_base,      ONLY : ispin, nspin, f, n=>nbsp
      USE cp_interfaces,       ONLY : dforce
  !
      USE gvecw, ONLY: ngw
  !
  !
      USE electrons_base,   ONLY: nx=>nbnd, nupdwn, iupdwn, nbspx, nbsp
      USE mp, ONLY: mp_sum 
      USE mp_global, ONLY: intra_image_comm 
!@@@@
      USE ldaU
!@@@@
  !
      IMPLICIT NONE
      INTEGER, INTENT(in) :: nfi
      REAL(DP) :: fccc, ccc
      REAL(DP) :: ema0bg(:), dt2bye
      REAL(DP) :: rhos(:,:)
      REAL(DP) :: bec(:,:)
      COMPLEX(DP) :: c0(:,:), cm(:,:)
      REAL(DP)    :: intermed
      LOGICAL, OPTIONAL, INTENT(in) :: fromscra
      LOGICAL, OPTIONAL, INTENT(in) :: restart
!
      REAL(DP) ::  verl1, verl2, verl3
      REAL(DP), ALLOCATABLE:: emadt2(:)
      REAL(DP), ALLOCATABLE:: emaver(:)
      COMPLEX(DP), ALLOCATABLE:: c2(:), c3(:)
      INTEGER :: i
      INTEGER :: iflag
      LOGICAL :: ttsde
!
       INTEGER     :: ierr,  nb, np_dw, is_dw, npair, n_unp, n_dwn, n_pair 
       REAL(DP)    :: ei_unp_mem, ei_unp_wfc
       COMPLEX(DP) :: intermed3
       REAL(DP),    ALLOCATABLE :: occ(:)
       COMPLEX(DP), ALLOCATABLE :: c4(:), c5(:)
!
! ... Controlling on sic applicability
!
       IF( lwf ) CALL errore('runcp_uspp_force_pairing', &
                           'Wannier function and sic are not compatibile',1)
       IF( tefield ) CALL errore('runcp_uspp_force_pairing', &
                           'Electric field and sic are not implemented',2)
       IF( nspin == 1 ) CALL errore(' runcp_force_pairing ',' inconsistent nspin ', 1)

!       
       ALLOCATE( emadt2( ngw ) )
       ALLOCATE( emaver( ngw ) )      
!
       iflag = 0
       IF( PRESENT( fromscra ) ) THEN
          IF( fromscra ) iflag = 1
       END IF
       IF( PRESENT( restart ) ) THEN
          IF( restart ) iflag = 2
       END IF
!       
       IF( iflag == 0 ) THEN
          ttsde  = tsde
       ELSE IF( iflag == 1 ) THEN
          ttsde = .TRUE.
       ELSE IF( iflag == 2 ) THEN
          ttsde = .FALSE.
       END IF
!
       ALLOCATE( c2(ngw), c3(ngw), c4(ngw), c5(ngw) )
       !
       ! ...  set verlet variables
       !
       verl1 = 2.0d0 * fccc
       verl2 = 1.0d0 - verl1
       verl3 = 1.0d0 * fccc 
!
       ccc    = fccc * dt2bye
       emadt2 = dt2bye * ema0bg
       emaver = emadt2 * verl3
!
       IF( nupdwn(1) /= (nupdwn(2) + 1) ) &
          CALL errore(' runcp_force_pairing ',' inconsistent number of states ', 1)

       n_unp = nupdwn(1)
       n_dwn = nupdwn(2)
       is_dw = iupdwn(2) 
       np_dw = nbsp 
!
       ALLOCATE( occ( nbspx ) )
!
       occ( 1:np_dw )  = 1.0d0
       occ( nbspx   )  = 0.0d0
!
! c0(dwn_paired) == c0(up_paired)
! cm(dwn_paired) == cm(up_paired)
! the nbspx dwn state has to be empty
!
!
      c0(:, is_dw:np_dw ) = c0(:, 1:n_dwn )
      cm(:, is_dw:np_dw ) = cm(:, 1:n_dwn )
!
      c0(:, nbspx ) = (0.d0, 0.d0)
      cm(:, nbspx ) = (0.d0, 0.d0)
!
     IF( MOD(n_unp, 2) == 0 ) npair = n_unp - 2
     IF( MOD(n_unp, 2) /= 0 ) npair = n_unp - 1

      DO i = 1, npair, 2 
         !
         CALL dforce(i,bec,vkb,c0,c2,c3,rhos(:,1:1),SIZE(rhos,1),ispin,f,n,nspin)
         CALL dforce(i,bec,vkb,c0,c4,c5,rhos(:,2:2),SIZE(rhos,1),ispin,f,n,nspin)
         !
         c2 = occ( i )*(c2 + c4)  
         c3 = occ(i+1)*(c3 + c5) 
         !
         IF( iflag == 2 ) THEN
            cm(:,i)        = c0(:,i)
            cm(:,i+1)      = c0(:,i+1)
         END IF
         !
         IF( ttsde ) THEN
             CALL wave_steepest( cm(:, i  ), c0(:, i  ), emaver, c2 )
             CALL wave_steepest( cm(:, i+1), c0(:, i+1), emaver, c3 )
         ELSE
             CALL wave_verlet( cm(:, i  ), c0(:, i  ), verl1, verl2, emaver, c2 )
             CALL wave_verlet( cm(:, i+1), c0(:, i+1), verl1, verl2, emaver, c3 )
         END IF
         !
         IF ( gstart == 2 ) THEN
                cm(1,  i)    = CMPLX(DBLE(cm(1,  i)),0.d0)
                cm(1, i+1)   = CMPLX(DBLE(cm(1,  i+1)),0.d0)
         END IF
      !
      END DO
      !
      IF( MOD(n_unp, 2) == 0 ) THEN

         npair = n_unp - 1 
!
         CALL dforce(npair,bec,vkb,c0,c2,c3,rhos(:,1:1),SIZE(rhos,1),ispin,f,n,nspin)
         CALL dforce(npair,bec,vkb,c0,c4,c5,rhos(:,2:2),SIZE(rhos,1),ispin,f,n,nspin)
!
         c2 = c2 + c4
         !
         IF( iflag == 2 ) cm( :, npair ) = c0( :, npair )
!
         IF( ttsde ) THEN
           CALL wave_steepest( cm(:, npair  ), c0(:, npair  ), emaver, c2 )
         ELSE
           CALL wave_verlet( cm(:, npair), c0(:, npair), verl1, verl2, emaver, c2 )
         ENDIF
!
         IF ( gstart == 2 ) cm(1, npair) = CMPLX(DBLE(cm(1, npair)),0.d0)

      ENDIF
!
      c0(:, is_dw:np_dw ) = c0(:, 1:n_dwn )
      cm(:, is_dw:np_dw ) = cm(:, 1:n_dwn )
!
      c0(:, nbspx ) = (0.d0, 0.d0)
      cm(:, nbspx ) = (0.d0, 0.d0)
!

!
! The electron unpaired is signed by n_unp and spin up 
! for the unpaired electron the ei_unp is the value of lambda
! "TRUE" ONLY WHEN THE POT is NORM_CONSERVING
!

      CALL dforce( n_unp, bec, vkb, c0, c2, c3, rhos, SIZE(rhos,1), ispin,f,n,nspin )
      !
      intermed  = - 2.d0 * sum(c2 * conjg(c0(:,n_unp)))
      IF ( gstart == 2 ) THEN
        intermed  = intermed + 1.d0 * c2(1) * conjg(c0(1,n_unp))
      END IF
      CALL mp_sum ( intermed, intra_image_comm )
      !           
      IF( iflag == 2 ) cm(:, n_unp) = c0(:, n_unp) 
      !
      IF( ttsde ) THEN
        CALL wave_steepest( cm(:, n_unp), c0(:, n_unp), emaver, c2 )
      ELSE
        CALL wave_verlet( cm(:, n_unp), c0(:, n_unp), verl1, verl2, emaver, c2 )
      ENDIF 
      !
      IF ( gstart == 2 ) cm(1, n_unp) = CMPLX(DBLE(cm(1, n_unp)),0.d0)
      !
      DEALLOCATE( occ )
      DEALLOCATE( emadt2 )
      DEALLOCATE( emaver )
      DEALLOCATE(c2, c4)
      DEALLOCATE(c3, c5)

   END SUBROUTINE runcp_uspp_force_pairing_x

