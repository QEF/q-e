!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!




   SUBROUTINE runcp_x &
      ( ttprint, tortho, tsde, cm, c0, cp, vpot, vkb, fi, ekinc, ht, ei, bec, fccc )

      !     This subroutine performs a Car-Parrinello or Steepest-Descent step
      !     on the electronic variables, computing forces on electrons and,
      !     when required, the eigenvalues of the Hamiltonian 
      !
      !     On output "cp" contains the new plave waves coefficients, while
      !     "cm" and "c0" are not changed
      !  ----------------------------------------------

      ! ...   declare modules
      USE kinds
      USE electrons_module,   ONLY : nb_l
      USE electrons_base,     ONLY : nupdwn, iupdwn, nspin, nbsp, nudx
      USE cp_electronic_mass, ONLY : emass
      USE cp_main_variables,  ONLY : ema0bg
      USE wave_base,          ONLY : hpsi
      USE cell_module,        ONLY : boxdimensions
      USE time_step,          ONLY : delt
      USE uspp,               ONLY : nkb
      USE gvecw,              ONLY : ngw
      USE cp_interfaces,      ONLY : runcp_ncpp, eigs, ortho, elec_fakekine, update_lambda

      IMPLICIT NONE

      ! ...   declare subroutine arguments

      LOGICAL :: ttprint, tortho, tsde
      COMPLEX(DP) :: cm(:,:), c0(:,:), cp(:,:)
      COMPLEX(DP)  ::  vkb(:,:)
      REAL(DP), INTENT(IN)  ::  fi(:)
      REAL(DP), INTENT(IN)  ::  bec(:,:)
      TYPE (boxdimensions), INTENT(IN)  ::  ht
      REAL (DP) ::  vpot(:,:)
      REAL(DP) :: ei(:,:)
      REAL(DP) :: ekinc
      REAL(DP), INTENT(IN) :: fccc

      ! ...   declare other variables
      INTEGER :: ierr, is

      REAL(DP),    ALLOCATABLE :: gam(:,:,:)

      ! ...   end of declarations

      ngw = SIZE( cp, 1 )
      !
      ALLOCATE( gam( nudx, nudx, nspin ), STAT=ierr)
      IF( ierr /= 0 ) CALL errore(' runcp ', ' allocating gam, prod ', ierr)

      ekinc    = 0.0d0

      !  Compute electronic forces and move electrons

      CALL runcp_ncpp( cm, c0, cp, vpot, vkb, fi, bec, fccc, gam, lambda = ttprint )

      !  Compute eigenstate
      !
      IF( ttprint ) THEN
        DO is = 1, nspin
          CALL eigs( nupdwn(is), gam(:,:,is), tortho, fi(iupdwn(is):iupdwn(is)+nupdwn(is)-1), ei(:,is) )
        END DO
      END IF

      !  Orthogonalize the new wave functions "cp"

      IF( tortho ) THEN
         CALL ortho( c0, cp, nupdwn, iupdwn, nspin )
      ELSE
         DO is = 1, nspin
            CALL gram( vkb, bec, nkb, cp(1,iupdwn(is)), SIZE(cp,1), nupdwn(is) )
         END DO
      END IF

      !  Compute fictitious kinetic energy of the electrons at time t

      CALL elec_fakekine( ekinc, ema0bg, emass, cp, cm, ngw, nbsp, 1, 2.0d0 * delt )

      DEALLOCATE( gam, STAT=ierr)
      IF( ierr /= 0 ) CALL errore(' runcp ', ' deallocating 1 ', ierr)


      RETURN
    END SUBROUTINE runcp_x


!=----------------------------------------------------------------------------------=!


    SUBROUTINE runcp_ncpp_x &
       ( cm, c0, cp, vpot, vkb, fi, bec, fccc, gam, lambda, fromscra, diis, restart )

       !     This subroutine performs a Car-Parrinello or Steepest-Descent step
       !     on the electronic variables, computing forces on electrons and,
       !     when required, the eigenvalues of the Hamiltonian 
       !
       !     On output "cp" contains the new plave waves coefficients, while
       !     "cm" and "c0" are not changed
       !  ----------------------------------------------

       ! ...   declare modules
      USE kinds
      USE electrons_base,     ONLY: nupdwn, iupdwn, nspin, nudx
      USE cp_electronic_mass, ONLY: emass
      USE cp_main_variables,  ONLY: ema0bg
      USE wave_base,          ONLY: wave_steepest, wave_verlet
      USE time_step,          ONLY: delt
      USE cp_interfaces,      ONLY: dforce, update_lambda
      USE control_flags,      ONLY: tsde
      USE gvecw,              ONLY: ngw
      USE reciprocal_vectors, ONLY: gzero

      IMPLICIT NONE

! ...   declare subroutine arguments

      COMPLEX(DP) :: cm(:,:), c0(:,:), cp(:,:)
      REAL(DP)    :: gam(:,:,:)
      COMPLEX(DP) :: vkb(:,:)
      REAL(DP), INTENT(IN)  ::  fi(:)
      REAL (DP) ::  vpot(:,:)
      REAL (DP), INTENT(IN) ::  bec(:,:)
      REAL(DP), INTENT(IN) :: fccc
      LOGICAL, OPTIONAL, INTENT(IN) :: lambda, fromscra, diis, restart

! ...   declare other variables
      REAL(DP) ::  svar1, svar2
      INTEGER :: i, ig, nx, nb, ierr, is
      INTEGER :: iflag, iwfc, nwfc

      COMPLEX(DP), ALLOCATABLE :: c2(:), c3(:)
      REAL(DP),    ALLOCATABLE :: svar3(:)
      LOGICAL :: tlam, ttsde


! ...   end of declarations
!  ----------------------------------------------

      IF( PRESENT( lambda ) ) THEN
        tlam = lambda
      ELSE
        tlam = .FALSE.
      END IF

      iflag = 0
      IF( PRESENT( fromscra ) ) THEN
        IF( fromscra ) iflag = 1
      END IF
      IF( PRESENT( restart ) ) THEN
        IF( restart ) iflag = 2
      END IF


      ALLOCATE( c2(ngw), c3(ngw), svar3(ngw), STAT = ierr )
      IF( ierr /= 0 ) CALL errore(' runcp_ncpp ', ' allocating c2, c3, svar3 ', ierr)

      !
      ! ...   determines friction dynamically according to the Nose' dynamics
      !

      IF( iflag == 0 ) THEN
        ttsde   = tsde
      ELSE IF ( iflag == 1 ) THEN
        ttsde   = .TRUE.
      ELSE IF ( iflag == 2 ) THEN
        ttsde   = .FALSE.
      END IF

      svar1   = 2.d0 * fccc
      svar2   = 1.d0 - svar1
      svar3( 1:ngw ) = delt * delt * ema0bg / emass * fccc


      DO is = 1, nspin

          nx   = nupdwn( is )

          nb = nx - MOD(nx, 2)

          DO i = 1, nb, 2

            iwfc = iupdwn(is)
            nwfc = nupdwn(is)

            CALL dforce( i, c0, fi, c2, c3, vpot(:,is), vkb, bec, nupdwn(is), iupdwn(is) )

            IF( tlam ) THEN
               CALL update_lambda( i  , gam( :, :,is), c0, c2, nwfc, iwfc )
               CALL update_lambda( i+1, gam( :, :,is), c0, c3, nwfc, iwfc )
            END IF

            IF( iflag == 2 ) THEN
              c0(:,i+iwfc-1) = cp(:,i+iwfc-1)
              c0(:,i+1+iwfc-1) = cp(:,i+1+iwfc-1)
            END IF

            IF ( ttsde ) THEN
              CALL wave_steepest( cp(:,i+iwfc-1), c0(:,i+iwfc-1), svar3, c2 )
              CALL wave_steepest( cp(:,i+1+iwfc-1), c0(:,i+1+iwfc-1), svar3, c3 )
            ELSE
              cp(:,i+iwfc-1) = cm(:,i+iwfc-1)
              cp(:,i+1+iwfc-1) = cm(:,i+1+iwfc-1)
              CALL wave_verlet( cp(:,i+iwfc-1), c0(:,i+iwfc-1), svar1, svar2, svar3, c2 )
              CALL wave_verlet( cp(:,i+1+iwfc-1), c0(:,i+1+iwfc-1), svar1, svar2, svar3, c3 )
            END IF

            IF( gzero ) cp(1,i+iwfc-1)  = DBLE( cp(1,i+iwfc-1) )
            IF( gzero ) cp(1,i+1+iwfc-1)= DBLE( cp(1,i+1+iwfc-1) )

          END DO

          IF( MOD(nx,2) /= 0) THEN

            nb = nx

            iwfc = iupdwn(is)
            nwfc = nupdwn(is)

            CALL dforce( nx, c0, fi, c2, c3, vpot(:,is), vkb, bec, nupdwn(is), iupdwn(is) )

            IF( tlam ) THEN
               CALL update_lambda( nb, gam( :, :,is), c0, c2, nwfc, iwfc )
            END IF

            IF( iflag == 2 ) THEN
              c0(:,nb+iwfc-1) = cp(:,nb+iwfc-1)
            END IF

            IF ( ttsde ) THEN
              CALL wave_steepest( cp(:,nb+iwfc-1), c0(:,nb+iwfc-1), svar3, c2 )
            ELSE
              cp(:,nb+iwfc-1) = cm(:,nb+iwfc-1)
              CALL wave_verlet( cp(:,nb+iwfc-1), c0(:,nb+iwfc-1), svar1, svar2, svar3, c2 )
            END IF
            IF( gzero ) cp(1,nb+iwfc-1) = DBLE( cp(1,nb+iwfc-1) )

          END IF

      END DO

      DEALLOCATE(svar3, c2, c3, STAT=ierr)
      IF( ierr /= 0 ) CALL errore(' runcp_ncpp ', ' deallocating 1 ', ierr)

      RETURN
    END SUBROUTINE runcp_ncpp_x


!=----------------------------------------------------------------------------------=!


!eigr==e^ig*r f is the occupation number
!fnl if the factor non local

    SUBROUTINE runcp_force_pairing_x &
       (ttprint, tortho, tsde, cm, c0, cp, vpot, vkb, fi, ekinc, ht, ei, bec, fccc)

!  same as runcp, except that electrons are paired forcedly
!  i.e. this handles a state dependant Hamiltonian for the paired and unpaired electrons
!  unpaired is assumed to exist, to be unique, and located in highest index band
!  ----------------------------------------------
!  END manual

! ...   declare modules
      USE kinds
      USE mp_global,          ONLY: intra_image_comm
      USE mp,                 ONLY: mp_sum
      USE electrons_module,   ONLY: nb_l
      USE electrons_base,     ONLY: iupdwn, nupdwn, nspin
      USE cp_electronic_mass, ONLY: emass
      USE cp_main_variables,  ONLY: ema0bg
      USE wave_base,          ONLY: wave_steepest, wave_verlet
      USE wave_base,          ONLY: hpsi
      USE cell_module,        ONLY: boxdimensions
      USE time_step,          ONLY: delt
      USE cp_interfaces,      ONLY: dforce, eigs, ortho, elec_fakekine, update_lambda
      USE constants,          ONLY: autoev
      USE io_global,          ONLY: ionode
      USE uspp,               ONLY: nkb
      use reciprocal_vectors, only: gzero, gstart
      USE gvecw,              ONLY: ngw

        IMPLICIT NONE

! ...   declare subroutine arguments

      LOGICAL :: ttprint, tortho, tsde
      COMPLEX(DP) :: cm(:,:), c0(:,:), cp(:,:)
      COMPLEX(DP)  ::  vkb(:,:)
      REAL(DP), INTENT(INOUT) ::  fi(:)
      TYPE (boxdimensions), INTENT(IN)  ::  ht
      REAL (DP) ::  vpot(:,:)
      REAL(DP) :: ei(:,:)
      REAL(DP), INTENT(IN) :: bec(:,:)
      REAL(DP) :: ekinc
      REAL(DP), INTENT(IN) :: fccc

! ...   declare other variables
      REAL(DP) :: s3, s4
      REAL(DP) :: svar1, svar2, etmp
      INTEGER :: i, ig, nx, nb, j, nb_g, nb_lx, ierr, ibl
      INTEGER :: ispin_wfc, n_unp, iwfc, nwfc 
      REAL(DP), ALLOCATABLE :: occup(:), occdown(:), occsum(:)
      REAL(DP) :: intermed, intermed2, ei_unp_mem, ei_unp_wfc
      COMPLEX(DP) ::  intermed3, intermed4


      COMPLEX(DP), ALLOCATABLE :: c2(:)
      COMPLEX(DP), ALLOCATABLE :: c3(:)
      COMPLEX(DP), ALLOCATABLE :: c4(:)
      COMPLEX(DP), ALLOCATABLE :: c5(:)
      REAL(DP),    ALLOCATABLE :: svar3(:)
      REAL(DP),    ALLOCATABLE :: gam(:,:)
      REAL(DP),    ALLOCATABLE :: ei_t(:,:)

! ...   end of declarations
!  ----------------------------------------------

      IF( nspin == 1 ) &
        CALL errore(' runcp_forced_pairing ',' inconsistent nspin ', 1)

      ALLOCATE(c2(ngw), c3(ngw), c4(ngw), c5(ngw), svar3(ngw), STAT=ierr)
      IF( ierr /= 0 ) CALL errore(' runcp_forced_pairing ', ' allocating c2, c3, svar3 ', ierr)


      svar1   = 2.d0 * fccc
      svar2   = 1.d0 - svar1
      svar3(1:ngw) = delt * delt * ema0bg / emass * fccc

      ekinc    = 0.0d0

      nx    = nupdwn(1)
      n_unp = nupdwn(1)

      IF( nupdwn(1) /= (nupdwn(2) + 1) ) &
        CALL errore(' runcp_forced_pairing ',' inconsistent spin numbers ', 1)


      nb_g = nupdwn(1)
      nb_lx = MAX( nb_l(1), nb_l(2) )
      nb_lx = MAX( nb_lx, 1 )

      ALLOCATE( gam(nb_lx,nb_g), STAT=ierr)
      IF( ierr /= 0 ) CALL errore(' runcp_forced_pairing ', ' allocating gam, prod ', ierr)

      ALLOCATE( occup(nx), occdown(nx), STAT=ierr )
      if ( ierr/=0 ) CALL errore(' runcp_forced_pairing ', 'allocating occup, occdown', ierr)

      ALLOCATE (ei_t(nx,2), STAT=ierr)
      IF( ierr /= 0 ) CALL errore(' runcp_forced_pairing ', 'allocating iei_t', ierr)

      occup   = 0.D0
      occdown = 0.D0
      occup(  1:nupdwn(1) )  = fi( 1:nupdwn(1) )
      occdown( 1:nupdwn(2) ) = fi( iupdwn(2):iupdwn(2)+nupdwn(2)-1 ) 


        IF( MOD( n_unp, 2 ) == 0 ) nb =  n_unp - 1
        IF( MOD( n_unp, 2 ) /= 0 ) nb =  n_unp - 2

        DO i = 1, nb, 2
          !
          CALL dforce( i, c0, fi, c2, c3, vpot(:,1), vkb, bec, nupdwn(2), iupdwn(1) )
          CALL dforce( i, c0, fi, c4, c5, vpot(:,2), vkb, bec, nupdwn(2), iupdwn(1) )
          !
          c2 = occup(i  )* (c2 + c4)
          c3 = occup(i+1)* (c3 + c5)

          IF( ttprint ) then
            !
            CALL update_lambda( i  , gam( :, :), c0, c2, nupdwn(1), 1 )
            CALL update_lambda( i+1, gam( :, :), c0, c3, nupdwn(1), 1 )

          END IF

          IF ( tsde ) THEN
             CALL wave_steepest( cp(:,i)  , c0(:,i)  , svar3, c2 )
             CALL wave_steepest( cp(:,i+1), c0(:,i+1), svar3, c3 )
          ELSE
            cp(:,i  ) = cm(:,i  )
            cp(:,i+1) = cm(:,i+1)
            CALL wave_verlet( cp(:,i), c0(:,i), svar1, svar2, svar3, c2 )
            CALL wave_verlet( cp(:,i+1), c0(:,i+1), svar1, svar2, svar3, c3 )
          END IF

          IF( gzero ) cp(1,i  )  = DBLE( cp(1,i) )
          IF( gzero ) cp(1,i+1)  = DBLE( cp(1,i+1) )

        END DO ! bande


        IF( MOD( n_unp, 2 ) /= 0 .and. n_unp > 1 ) THEN
          !
          nb = n_unp - 1
          !
          CALL dforce( nb, c0, fi, c2, c3, vpot(:,1), vkb, bec, nupdwn(2), iupdwn(1) )
          CALL dforce( nb, c0, fi, c4, c5, vpot(:,2), vkb, bec, nupdwn(2), iupdwn(1) )

          c2 = occup(nb)* (c2 + c4)

          IF( ttprint ) THEN
            CALL update_lambda( nb, gam( :, :), c0, c2, nupdwn(2), iupdwn(1) )
          END IF

          IF ( tsde ) THEN
             CALL wave_steepest( cp(:,nb), c0(:,nb), svar3, c2 )
          ELSE
             cp(:,nb) = cm(:,nb)
             CALL wave_verlet( cp(:,nb), c0(:,nb), svar1, svar2, svar3, c2 )
          END IF
          IF( gzero ) cp(1,nb) = DBLE( cp(1,nb) )
        END IF

        !
        CALL dforce( n_unp, c0, fi, c2, c3, vpot(:,1), vkb, bec, nupdwn(1), iupdwn(1) )

        intermed  = -2.d0 * sum( c2 * conjg( c0(:, n_unp ) ) )
        IF ( gstart == 2 ) THEN
           intermed  = intermed + 1.d0 * c2(1) * conjg( c0( 1, n_unp ) )
         END IF
        intermed3 = sum(c0(:,n_unp) * conjg( c0(:, n_unp)))

        CALL mp_sum ( intermed, intra_image_comm )
        CALL mp_sum ( intermed3, intra_image_comm )
        !  Eigenvalue of unpaired
        ei_unp_mem = intermed

        !  <Phiunpaired|Phiunpaired>
        ei_unp_wfc = intermed3
        !write(6,*) '  <psi|psi> = ', intermed3, '  ei_unp(au) = ', intermed

        IF ( tsde ) THEN
           CALL wave_steepest( cp( :, n_unp ), c0( :, n_unp ), svar3, c2 )
        ELSE
          cp( :, n_unp ) = cm( :, n_unp )
          CALL wave_verlet( cp( :, n_unp ), c0( :, n_unp ), svar1, svar2, svar3, c2 )
        END IF
        IF( gzero ) cp( 1, n_unp ) = DBLE( cp( 1, n_unp ) )


        IF( ttprint ) THEN

            ALLOCATE( occsum( nupdwn(1) ) )
            occsum(:) = occup(:) + occdown(:)

            CALL eigs( nupdwn(2), gam, tortho, occsum, ei(:,1) )

            DEALLOCATE( occsum )
            DO i = 1, nupdwn(2)
               ei( i, 2 ) = ei( i , 1)
            END DO

            ei( nupdwn(1), 1) = ei_unp_mem
            ei( nupdwn(1), 2) = 0.d0

            WRITE(6,*) 'SIC EIGENVALUES(eV), dwn and up electrons Kpoint',1
            WRITE(6,1004) ( ei( i, 2 ) * autoev, i = 1, nupdwn(2) )
            WRITE(6,1005) ( ei( i, 1 ) * autoev, i = 1, nupdwn(1) )

1004        FORMAT(/,3X,'SIC EIGENVALUES DW=',3X,10F8.2)
1005        FORMAT(/,3X,'SIC EIGENVALUES UP=',3X,10F8.2)

        ENDIF

      IF( tortho ) THEN
         CALL ortho( c0, cp, nupdwn, iupdwn, nspin )
      ELSE
         CALL gram( vkb, bec, nkb, cp(1,1), SIZE(cp,1), nupdwn(1) )
         cp(:,iupdwn(2):iupdwn(2)+nupdwn(2)-1) = cp(:,1:nupdwn(2))
      END IF

      !  Compute fictitious kinetic energy of the electrons at time t

      CALL elec_fakekine( etmp , ema0bg, emass, cp, cm, ngw, nupdwn(1), 1, 2.0d0 * delt )
      CALL elec_fakekine( ekinc, ema0bg, emass, cp, cm, ngw, nupdwn(2), 1, 2.0d0 * delt )

      ekinc = ekinc + etmp

      DEALLOCATE( ei_t, svar3, c2, c3, c4, c5, gam, occup, occdown, STAT=ierr)
      IF( ierr /= 0 ) CALL errore(' runcp_force_pairing ', ' deallocating ', ierr)

      RETURN
    END SUBROUTINE runcp_force_pairing_x


!=----------------------------------------------------------------------------=!


   SUBROUTINE runcp_uspp_x &
      ( nfi, fccc, ccc, ema0bg, dt2bye, rhos, bec, c0, cm, fromscra, restart )
      !
      USE kinds,             ONLY: DP
      use wave_base, only: wave_steepest, wave_verlet
      use control_flags, only: lwf, tsde
      !use uspp, only : nhsa=> nkb, betae => vkb, rhovan => becsum, deeq
      use uspp, only : deeq, betae => vkb
      use reciprocal_vectors, only : gstart
      use electrons_base, only : n=>nbsp, ispin, f, nspin
      use wannier_subroutines, only: ef_potential
      use efield_module, only: dforce_efield, tefield, dforce_efield2, tefield2

      use gvecw, only: ngw
     !
     IMPLICIT NONE
     integer, intent(in) :: nfi
     real(DP) :: fccc, ccc
     real(DP) :: ema0bg(:), dt2bye
     real(DP) :: rhos(:,:)
     real(DP) :: bec(:,:)
     complex(DP) :: c0(:,:), cm(:,:)
     logical, optional, intent(in) :: fromscra
     logical, optional, intent(in) :: restart
     !
     !
     real(DP) ::  verl1, verl2, verl3
     real(DP), allocatable:: emadt2(:)
     real(DP), allocatable:: emaver(:)
     complex(DP), allocatable:: c2(:), c3(:)
     integer :: i
     integer :: iflag
     logical :: ttsde

     iflag = 0
     IF( PRESENT( fromscra ) ) THEN
       IF( fromscra ) iflag = 1
     END IF
     IF( PRESENT( restart ) ) THEN
       IF( restart ) iflag = 2
     END IF

     !
     ! ...  set verlet variables 
     !
     verl1 = 2.0d0 * fccc
     verl2 = 1.0d0 - verl1
     verl3 = 1.0d0 * fccc

     allocate(c2(ngw))
     allocate(c3(ngw))
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

      if( lwf ) then
        call ef_potential( nfi, rhos, bec, deeq, betae, c0, cm, emadt2, emaver, verl1, verl2, c2, c3 )
      else
        do i=1,n,2
           call dforce(bec,betae,i,c0(1,i),c0(1,i+1),c2,c3,rhos,ispin,f,n,nspin)
           if( tefield ) then
             CALL dforce_efield ( bec, i, c0, c2, c3, rhos)
           end if
           if( tefield2 ) then
             CALL dforce_efield2 ( bec, i, c0, c2, c3, rhos)
           end if
           IF( iflag == 2 ) THEN
             cm(:,i)   = c0(:,i)
             cm(:,i+1) = c0(:,i+1)
           END IF
           if( ttsde ) then
              CALL wave_steepest( cm(:, i  ), c0(:, i  ), emaver, c2 )
              CALL wave_steepest( cm(:, i+1), c0(:, i+1), emaver, c3 )
           else
              CALL wave_verlet( cm(:, i  ), c0(:, i  ), verl1, verl2, emaver, c2 )
              CALL wave_verlet( cm(:, i+1), c0(:, i+1), verl1, verl2, emaver, c3 )
           endif
           if ( gstart == 2) THEN
              cm(1,  i)=CMPLX(DBLE(cm(1,  i)),0.d0)
              cm(1,i+1)=CMPLX(DBLE(cm(1,i+1)),0.d0)
           end if
        end do
      end if

     DEALLOCATE( emadt2 )
     DEALLOCATE( emaver )
     deallocate(c2)
     deallocate(c3)
!
   END SUBROUTINE runcp_uspp_x
!
!
!=----------------------------------------------------------------------------=!
!
!

   SUBROUTINE runcp_uspp_bgl_x &
      ( nfi, fccc, ccc, ema0bg, dt2bye, rhos, bec, c0, cm, fromscra, restart )
     !
      USE kinds,             ONLY: DP
     use wave_base,              only: my_wave_steepest, my_wave_verlet
     use control_flags,          only: lwf, tsde
     use uspp,                   only: deeq, betae => vkb
     use reciprocal_vectors,     only: gstart
     use electrons_base,         only: n=>nbsp, nspin
     use wannier_subroutines,    only: ef_potential
     use efield_module,          only: dforce_efield, tefield, dforce_efield2, tefield2
     use gvecw,                  only: ngw
     use smooth_grid_dimensions, only: nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx, nnrsx
     USE fft_base,               ONLY: dffts
     USE mp_global,              ONLY: me_image, nogrp, me_ogrp
     USE parallel_include
     use task_groups

     !
     IMPLICIT NONE
     integer, intent(in) :: nfi
     real(DP) :: fccc, ccc
     real(DP) :: ema0bg(:), dt2bye
     real(DP) :: rhos(:,:)
     real(DP) :: bec(:,:)
     complex(DP) :: c0(:,:), cm(:,:)
     logical, optional, intent(in) :: fromscra
     logical, optional, intent(in) :: restart
     !
     !
     !
     real(8) ::  verl1, verl2, verl3
     real(8), allocatable:: emadt2(:)
     real(8), allocatable:: emaver(:)
     complex(8), allocatable:: c2(:), c3(:)
     integer :: i, index_in, index
     integer :: iflag, ierr
     logical :: ttsde

     iflag = 0
     IF( PRESENT( fromscra ) ) THEN
       IF( fromscra ) iflag = 1
     END IF
     IF( PRESENT( restart ) ) THEN
       IF( restart ) iflag = 2
     END IF

     !
     ! ...  set verlet variables 
     !
     verl1 = 2.0d0 * fccc
     verl2 = 1.0d0 - verl1
     verl3 = 1.0d0 * fccc

     allocate(c2(ngw))
     allocate(c3(ngw))
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

     if( lwf ) then
        !
        call ef_potential( nfi, rhos, bec, deeq, betae, c0, cm, emadt2, emaver, verl1, verl2, c2, c3 )
        !
     else

        IF (.NOT.(ALLOCATED(tg_c2))) ALLOCATE(tg_c2((NOGRP+1)*ngw))
        IF (.NOT.(ALLOCATED(tg_c3))) ALLOCATE(tg_c3((NOGRP+1)*ngw))

        !---------------------------------------------------------------
        !This loop is parallelized accross the eigenstates as well as
        !in the FFT, similar to rhoofr
        !---------------------------------------------------------------

        !------------------------------------
        !The potential in rhos
        !is distributed accros all processors
        !We need to redistribute it so that
        !it is completely contained in the
        !processors of an orbital TASK-GROUP
        !------------------------------------
        !
        recv_cnt(1)   = dffts%npp(NOLIST(1)+1)*nr1sx*nr2sx
        recv_displ(1) = 0
        DO i = 2, NOGRP
           recv_cnt(i) = dffts%npp(NOLIST(i)+1)*nr1sx*nr2sx
           recv_displ(i) = recv_displ(i-1) + recv_cnt(i)
        ENDDO
        IF (.NOT.ALLOCATED(tg_rhos)) ALLOCATE(tg_rhos( (NOGRP+1)*nr1sx*nr2sx*maxval(dffts%npp),nspin))

        tg_c3(:) = 0D0
        tg_c3(:) = 0D0
        tg_rhos(:,:) = 0D0

#if defined (__PARA) && defined (__MPI)
        DO i = 1, nspin
           CALL MPI_Allgatherv(rhos(1,i), dffts%npp(me_image+1)*nr1sx*nr2sx, MPI_DOUBLE_PRECISION, &
                tg_rhos(1,i), recv_cnt, recv_displ, MPI_DOUBLE_PRECISION, ME_OGRP, IERR)
        ENDDO
#endif
        do i = 1, n, 2*NOGRP ! 2*NOGRP eigenvalues are treated at each iteration

           !----------------------------------------------------------------
           !The input coefficients to dforce cover eigenstates i:i+2*NOGRP-1
           !Thus, in dforce the dummy arguments for c0(1,i,1,1) and
           !c0(1,i+1,1,1) hold coefficients for eigenstates i,i+2*NOGRP-2,2
           !and i+1,i+2*NOGRP...for example if NOGRP is 4 then we would have
           !1-3-5-7 and 2-4-6-8
           !----------------------------------------------------------------

           call dforce_bgl( bec, betae, i, c0(1,i), c0(1,i+1), tg_c2, tg_c3, tg_rhos)

           !-------------------------------------------------------
           !C. Bekas: This is not implemented yet! I need to see it
           !-------------------------------------------------------

           if( tefield ) then
             CALL errore( ' runcp_uspp ', ' electric field on BGL not implemented yet ', 1 )
           end if

           IF( iflag == 2 ) THEN
             DO index = 1, 2 * NOGRP, 2
                cm(:,i+index-1) = c0(:,i+index-1)
                cm(:,i+index) = c0(:,i+index)
             ENDDO
           END IF

           index_in = 1
           DO index = 1, 2*NOGRP, 2
              IF (tsde) THEN
                 CALL my_wave_steepest( cm(:, i+index-1 ), c0(:, i+index-1 ), emaver, tg_c2, ngw, index_in )
                 CALL my_wave_steepest( cm(:, i+index), c0(:, i+index), emaver, tg_c3, ngw, index_in )
              ELSE
                 CALL my_wave_verlet( cm(:, i+index-1 ), c0(:, i+index-1 ), &
                      verl1, verl2, emaver, tg_c2, ngw, index_in )
                 CALL my_wave_verlet( cm(:, i+index), c0(:, i+index ), &
                      verl1, verl2, emaver, tg_c3, ngw, index_in )

              ENDIF
              if ( gstart == 2 ) then
                 cm(1,  i+index-1)=cmplx(real(cm(1,  i+index-1)),0.0)
                 cm(1,i+index)=cmplx(real(cm(1,i+index)),0.0)
              end if
             index_in = index_in+1

           ENDDO ! End loop accross 2*NOGRP current eigenstates

        end do ! End loop accross eigenstates

     end if

     DEALLOCATE( emadt2 )
     DEALLOCATE( emaver )
     deallocate(c2)
     deallocate(c3)


   END SUBROUTINE runcp_uspp_bgl_x
!
!=----------------------------------------------------------------------------=!
!=----------------------------------------------------------------------------=!

    SUBROUTINE runcp_uspp_force_pairing_x  &
       ( nfi, fccc, ccc, ema0bg, dt2bye, rhos, bec, c0, cm, intermed, fromscra, restart )
  !
      USE kinds,               ONLY : DP
      USE wave_base,           ONLY : wave_steepest, wave_verlet
      USE control_flags,       ONLY : lwf, tsde
  !   use uspp,                only : nhsa=> nkb, betae => vkb, rhovan => becsum, deeq
      USE uspp,                ONLY : deeq, betae => vkb
      USE reciprocal_vectors,  ONLY : gstart
      USE wannier_subroutines, ONLY : ef_potential
      USE efield_module,       ONLY : dforce_efield, tefield
      USE electrons_base,      ONLY : ispin, nspin, f, n=>nbsp
  !
      USE gvecw, ONLY: ngw
  !
  !
      USE electrons_base,   ONLY: nx=>nbnd, nupdwn, iupdwn, nbspx, nbsp
      USE mp, ONLY: mp_sum 
      USE mp_global, ONLY: intra_image_comm 
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
         CALL dforce(bec,betae,i,c0(1,i),c0(1,i+1),c2,c3,rhos(1,1),ispin,f,n,nspin)
         CALL dforce(bec,betae,i,c0(1,i),c0(1,i+1),c4,c5,rhos(1,2),ispin,f,n,nspin)
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
         CALL dforce(bec,betae,npair,c0(1,npair),c0(1,nbspx),c2,c3,rhos(1,1),ispin,f,n,nspin)
         CALL dforce(bec,betae,npair,c0(1,npair),c0(1,nbspx),c4,c5,rhos(1,2),ispin,f,n,nspin)
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

      CALL dforce( bec, betae, n_unp, c0(1,n_unp), c0(1,n_unp), c2, c3, rhos(1,1),ispin,f,n,nspin )
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

