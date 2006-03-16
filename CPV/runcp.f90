!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!  AB INITIO COSTANT PRESSURE MOLECULAR DYNAMICS

!=----------------------------------------------------------------------------=!
  MODULE runcp_module
!=----------------------------------------------------------------------------=!

        IMPLICIT NONE
        PRIVATE
        SAVE

        PUBLIC :: runcp, runcp_force_pairing
        PUBLIC :: runcp_uspp, runcp_uspp_force_pairing, runcp_ncpp

!=----------------------------------------------------------------------------=!
        CONTAINS
!=----------------------------------------------------------------------------=!


!  ----------------------------------------------
!  BEGIN manual

    SUBROUTINE runcp( ttprint, tortho, tsde, cm, c0, cp, cdesc, &
      vpot, eigr, fi, ekinc, timerd, timeorto, ht, ei, bec, fccc )

!     This subroutine performs a Car-Parrinello or Steepest-Descent step
!     on the electronic variables, computing forces on electrons and,
!     when required, the eigenvalues of the Hamiltonian 
!
!     On output "cp" contains the new plave waves coefficients, while
!     "cm" and "c0" are not changed
!  ----------------------------------------------
!  END manual

! ...   declare modules
      USE kinds
      USE mp_global, ONLY: mpime, nproc
      USE mp, ONLY: mp_sum
      USE electrons_module, ONLY:  pmss, eigs, nb_l
      USE cp_electronic_mass, ONLY: emass
      USE wave_functions, ONLY : cp_kinetic_energy
      USE wave_base, ONLY: hpsi
      USE cell_module, ONLY: boxdimensions
      USE time_step, ONLY: delt
      USE forces, ONLY: dforce
      USE orthogonalize, ONLY: ortho
      USE wave_types, ONLY: wave_descriptor
      USE pseudo_projector, ONLY: projector
      USE wave_constrains, ONLY: update_lambda
      USE reciprocal_space_mesh, ONLY: gkmask_l
      USE uspp,             ONLY : vkb, nkb

      IMPLICIT NONE

! ...   declare subroutine arguments

      LOGICAL :: ttprint, tortho, tsde
      COMPLEX(DP) :: cm(:,:,:,:), c0(:,:,:,:), cp(:,:,:,:)
      TYPE (wave_descriptor), INTENT(IN) :: cdesc
      COMPLEX(DP)  ::  eigr(:,:)
      REAL(DP), INTENT(IN)  ::  fi(:,:,:)
      REAL(DP), INTENT(IN)  ::  bec(:,:)
      TYPE (boxdimensions), INTENT(IN)  ::  ht
      REAL (DP) ::  vpot(:,:)
      REAL(DP) :: ei(:,:,:)
      REAL(DP) :: timerd, timeorto
      REAL(DP) :: ekinc(:)
      REAL(DP), INTENT(IN) :: fccc

! ...   declare other variables
      REAL(DP) :: s1, s2, s3, s4
      INTEGER :: ik, nx, nb_lx, ierr, nkl, is

      COMPLEX(DP), ALLOCATABLE :: cgam(:,:,:)
      REAL(DP),    ALLOCATABLE :: gam(:,:,:)

      REAL(DP), EXTERNAL :: cclock

! ...   end of declarations
!  ----------------------------------------------

      s1 = cclock()

      nb_lx = MAX( nb_l(1), nb_l(2) )
      nb_lx = MAX( nb_lx, 1 )
      IF( cdesc%gamma ) THEN
        ALLOCATE( cgam(1,1,1), gam( nb_lx, SIZE( c0, 2 ), cdesc%nspin ), STAT=ierr)
      ELSE
        ALLOCATE( cgam( nb_lx, SIZE( c0, 2 ), cdesc%nspin ), gam(1,1,1), STAT=ierr)
      END IF
      IF( ierr /= 0 ) CALL errore(' runcp ', ' allocating gam, prod ', ierr)

      ekinc    = 0.0d0
      timerd   = 0.0d0
      timeorto = 0.0d0

      !  Compute electronic forces and move electrons

      CALL runcp_ncpp( cm, c0, cp, cdesc, vpot, eigr, fi, bec, fccc, &
           gam, cgam, lambda = ttprint )

      !  Compute eigenstate
      !
      IF( ttprint ) THEN
        DO is = 1, cdesc%nspin
          nx = cdesc%nbt( is )
          nkl  = cdesc%nkl
          DO ik = 1, nkl
              CALL eigs( nx, gam(:,:,is), cgam(:,:,is), tortho, fi(:,ik,is), ei(:,ik,is), cdesc%gamma )
          END DO
        END DO
      END IF

      s2 = cclock()
      timerd = s2 - s1

      !  Orthogonalize the new wave functions "cp"

      IF( tortho ) THEN
         CALL ortho(c0, cp, cdesc, pmss, emass)
      ELSE
         DO is = 1, cdesc%nspin
            CALL gram( vkb, bec, nkb, cp(1,1,1,is), SIZE(cp,1), cdesc%nbt( is ) )
         END DO
      END IF

      s3 = cclock()
      timeorto = s3 - s2

      !  Compute fictitious kinetic energy of the electrons at time t

      DO is = 1, cdesc%nspin
        ekinc(is) = cp_kinetic_energy( is, cp(:,:,:,is), cm(:,:,:,is), cdesc, pmss, delt)
      END DO

      DEALLOCATE( cgam, gam, STAT=ierr)
      IF( ierr /= 0 ) CALL errore(' runcp ', ' deallocating 1 ', ierr)


      RETURN
    END SUBROUTINE runcp


!=----------------------------------------------------------------------------------=!


!  ----------------------------------------------
!  BEGIN manual

    SUBROUTINE runcp_ncpp( cm, c0, cp, cdesc, &
      vpot, eigr, fi, bec, fccc, gam, cgam, lambda, fromscra, diis, restart )

!     This subroutine performs a Car-Parrinello or Steepest-Descent step
!     on the electronic variables, computing forces on electrons and,
!     when required, the eigenvalues of the Hamiltonian 
!
!     On output "cp" contains the new plave waves coefficients, while
!     "cm" and "c0" are not changed
!  ----------------------------------------------
!  END manual

! ...   declare modules
      USE kinds
      USE mp_global, ONLY: mpime, nproc
      USE mp, ONLY: mp_sum
      USE electrons_module, ONLY:  pmss
      USE cp_electronic_mass, ONLY: emass
      USE wave_base, ONLY: wave_steepest, wave_verlet
      USE time_step, ONLY: delt
      USE forces, ONLY: dforce
      USE wave_types, ONLY: wave_descriptor
      USE wave_constrains, ONLY: update_lambda
      USE control_flags, ONLY: tsde
      USE pseudo_projector, ONLY: projector
      USE reciprocal_space_mesh, ONLY: gkmask_l

      IMPLICIT NONE

! ...   declare subroutine arguments

      COMPLEX(DP) :: cm(:,:,:,:), c0(:,:,:,:), cp(:,:,:,:)
      COMPLEX(DP) :: cgam(:,:,:)
      REAL(DP)    :: gam(:,:,:)
      TYPE (wave_descriptor), INTENT(IN) :: cdesc
      COMPLEX(DP) :: eigr(:,:)
      REAL(DP), INTENT(IN)  ::  fi(:,:,:)
      REAL (DP) ::  vpot(:,:)
      REAL (DP), INTENT(IN) ::  bec(:,:)
      REAL(DP), INTENT(IN) :: fccc
      LOGICAL, OPTIONAL, INTENT(IN) :: lambda, fromscra, diis, restart

! ...   declare other variables
      REAL(DP) ::  svar1, svar2
      INTEGER :: i, ig, nx, ngw, nb, ierr, is
      INTEGER :: iflag

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


      ngw  = cdesc%ngwl

      ALLOCATE( c2(ngw), c3(ngw), svar3(ngw), STAT = ierr )
      IF( ierr /= 0 ) CALL errore(' runcp_ncpp ', ' allocating c2, c3, svar3 ', ierr)

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
      svar3( 1:ngw ) = delt * delt / pmss( 1:ngw ) * fccc


      DO is = 1, cdesc%nspin

        nx   = cdesc%nbt( is )
        IF( nx > SIZE( fi, 1 ) ) &
          CALL errore(' runcp ',' inconsistent occupation numbers ', 1)

          nb = nx - MOD(nx, 2)

          DO i = 1, nb, 2

            !WRITE( 6, * ) 'DEBUG = ', fi(i,1,is), fi(i+1,1,is)
            CALL dforce( i, is, c0(:,:,1,is), cdesc, fi(:,1,is), c2, c3, vpot(:,is), eigr, bec )

            IF( tlam ) THEN
               CALL update_lambda( i, gam( :, :,is), c0(:,:,1,is), cdesc, c2 )
               CALL update_lambda( i+1, gam( :, :,is), c0(:,:,1,is), cdesc, c3 )
            END IF

            IF( iflag == 2 ) THEN
              c0(:,i,1,is) = cp(:,i,1,is)
              c0(:,i+1,1,is) = cp(:,i+1,1,is)
            END IF

            IF ( ttsde ) THEN
              CALL wave_steepest( cp(:,i,1,is), c0(:,i,1,is), svar3, c2 )
              CALL wave_steepest( cp(:,i+1,1,is), c0(:,i+1,1,is), svar3, c3 )
            ELSE
              cp(:,i,1,is) = cm(:,i,1,is)
              cp(:,i+1,1,is) = cm(:,i+1,1,is)
              CALL wave_verlet( cp(:,i,1,is), c0(:,i,1,is), svar1, svar2, svar3, c2 )
              CALL wave_verlet( cp(:,i+1,1,is), c0(:,i+1,1,is), svar1, svar2, svar3, c3 )
            END IF

            IF( cdesc%gzero ) cp(1,i,1,is)  = DBLE( cp(1,i,1,is) )
            IF( cdesc%gzero ) cp(1,i+1,1,is)= DBLE( cp(1,i+1,1,is) )

          END DO

          IF( MOD(nx,2) /= 0) THEN

            nb = nx

            CALL dforce( nx, is, c0(:,:,1,is), cdesc, fi(:,1,is), c2, c3, vpot(:,is), eigr, bec )

            IF( tlam ) THEN
               CALL update_lambda( nb, gam( :, :,is), c0(:,:,1,is), cdesc, c2 )
            END IF

            IF( iflag == 2 ) THEN
              c0(:,nb,1,is) = cp(:,nb,1,is)
            END IF

            IF ( ttsde ) THEN
              CALL wave_steepest( cp(:,nb,1,is), c0(:,nb,1,is), svar3, c2 )
            ELSE
              cp(:,nb,1,is) = cm(:,nb,1,is)
              CALL wave_verlet( cp(:,nb,1,is), c0(:,nb,1,is), svar1, svar2, svar3, c2 )
            END IF
            IF( cdesc%gzero ) cp(1,nb,1,is) = DBLE( cp(1,nb,1,is) )

          END IF

      END DO

      DEALLOCATE(svar3, c2, c3, STAT=ierr)
      IF( ierr /= 0 ) CALL errore(' runcp_ncpp ', ' deallocating 1 ', ierr)

      RETURN
    END SUBROUTINE runcp_ncpp


!=----------------------------------------------------------------------------------=!


!cdesc is the desciptor for the wf
!eigr==e^ig*r f is the occupation number
!fnl if the factor non local

    SUBROUTINE runcp_force_pairing(ttprint, tortho, tsde, cm, c0, cp, cdesc, &
        vpot, eigr, fi, ekinc, timerd, timeorto, ht, ei, bec, fccc)

!  same as runcp, except that electrons are paired forcedly
!  i.e. this handles a state dependant Hamiltonian for the paired and unpaired electrons
!  unpaired is assumed to exist, to be unique, and located in highest index band
!  ----------------------------------------------
!  END manual

! ...   declare modules
      USE kinds
      USE mp_global, ONLY: mpime, nproc, group
      USE mp, ONLY: mp_sum
      USE electrons_module, ONLY: pmss, eigs, nb_l, nupdwn, nspin
      USE cp_electronic_mass, ONLY: emass
      USE wave_functions, ONLY : cp_kinetic_energy
      USE wave_base, ONLY: wave_steepest, wave_verlet
      USE wave_base, ONLY: hpsi
      USE cell_module, ONLY: boxdimensions
      USE time_step, ONLY: delt
      USE forces, ONLY: dforce
      USE orthogonalize, ONLY: ortho
      USE wave_types, ONLY: wave_descriptor
      USE constants, ONLY: au
      USE io_global, ONLY: ionode
      USE wave_constrains, ONLY: update_lambda
      USE reciprocal_space_mesh, ONLY: gkmask_l
      USE uspp,             ONLY : vkb, nkb

        IMPLICIT NONE

! ...   declare subroutine arguments

      LOGICAL :: ttprint, tortho, tsde
      COMPLEX(DP) :: cm(:,:,:,:), c0(:,:,:,:), cp(:,:,:,:)
      TYPE (wave_descriptor), INTENT(IN) :: cdesc
      COMPLEX(DP)  ::  eigr(:,:)
      REAL(DP), INTENT(INOUT) ::  fi(:,:,:)
      TYPE (boxdimensions), INTENT(IN)  ::  ht
      REAL (DP) ::  vpot(:,:)
      REAL(DP) :: ei(:,:,:)
      REAL(DP), INTENT(IN) :: bec(:,:)
      REAL(DP) :: timerd, timeorto
      REAL(DP) :: ekinc(:)
      REAL(DP), INTENT(IN) :: fccc

! ...   declare other variables
      REAL(DP) :: s3, s4
      REAL(DP) ::  svar1, svar2
      INTEGER :: i, ik,ig, nx, ngw, nb, j, nb_g, nb_lx, ierr, nkl, ibl
      INTEGER :: ispin_wfc, n_unp 
      REAL(DP), ALLOCATABLE :: occup(:,:), occdown(:,:), occsum(:)
      REAL(DP) :: intermed, intermed2, ei_unp_mem, ei_unp_wfc
      COMPLEX(DP) ::  intermed3, intermed4


      COMPLEX(DP), ALLOCATABLE :: c2(:)
      COMPLEX(DP), ALLOCATABLE :: c3(:)
      COMPLEX(DP), ALLOCATABLE :: c4(:)
      COMPLEX(DP), ALLOCATABLE :: c5(:)
      COMPLEX(DP), ALLOCATABLE :: cgam(:,:)
      REAL(DP),    ALLOCATABLE :: svar3(:)
      REAL(DP),    ALLOCATABLE :: gam(:,:)
      REAL(DP),    ALLOCATABLE :: ei_t(:,:,:)

      REAL(DP), EXTERNAL :: cclock

! ...   end of declarations
!  ----------------------------------------------

      IF( nspin == 1 ) &
        CALL errore(' runcp_forced_pairing ',' inconsistent nspin ', 1)

      nkl  = cdesc%nkl
      IF( nkl /= SIZE( fi, 2 ) ) &
        CALL errore(' runcp_forced_pairing ',' inconsistent number of kpoints ', 1)

      ngw  = cdesc%ngwl

      ALLOCATE(c2(ngw), c3(ngw), c4(ngw), c5(ngw), svar3(ngw), STAT=ierr)
      IF( ierr /= 0 ) CALL errore(' runcp_forced_pairing ', ' allocating c2, c3, svar3 ', ierr)


      svar1   = 2.d0 * fccc
      svar2   = 1.d0 - svar1
      svar3(1:ngw) = delt * delt / pmss(1:ngw) * fccc

      ekinc    = 0.0d0
      timerd   = 0.0d0
      timeorto = 0.0d0

      nx    = cdesc%nbt( 1 )
      n_unp = nupdwn(1)

      IF( nx /= SIZE( fi, 1 ) ) &
        CALL errore(' runcp_forced_pairing ',' inconsistent occupation numbers ', 1)

      IF( nupdwn(1) /= (nupdwn(2) + 1) ) &
        CALL errore(' runcp_forced_pairing ',' inconsistent spin numbers ', 1)


      nb_g = cdesc%nbt( 1 )
      nb_lx = MAX( nb_l(1), nb_l(2) )
      nb_lx = MAX( nb_lx, 1 )

      IF( cdesc%gamma ) THEN
        ALLOCATE(cgam(1,1), gam(nb_lx,nb_g), STAT=ierr)
      ELSE
        ALLOCATE(cgam(nb_lx,nb_g), gam(1,1), STAT=ierr)
      END IF
      IF( ierr /= 0 ) CALL errore(' runcp_forced_pairing ', ' allocating gam, prod ', ierr)

      ALLOCATE( occup(nx,nkl), occdown(nx,nkl), STAT=ierr )
      if ( ierr/=0 ) CALL errore(' runcp_forced_pairing ', 'allocating occup, occdown', ierr)

      ALLOCATE (ei_t(nx,nkl,2), STAT=ierr)
      IF( ierr /= 0 ) CALL errore(' runcp_forced_pairing ', 'allocating iei_t', ierr)

      occup   = 0.D0
      occdown = 0.D0
      occup(  1:nupdwn(1), 1:nkl )  = fi( 1:nupdwn(1), 1:nkl, 1 )
      occdown( 1:nupdwn(2), 1:nkl ) = fi( 1:nupdwn(2), 1:nkl, 2 ) 

      !  ciclo sui punti K

      KAPPA: DO ik = 1, nkl

        s4 = cclock()

        IF( MOD( n_unp, 2 ) == 0 ) nb =  n_unp - 1
        IF( MOD( n_unp, 2 ) /= 0 ) nb =  n_unp - 2

        DO i = 1, nb, 2
          !
          CALL dforce( i, 2, c0(:,:,1,1), cdesc, fi(:,1,1), c2, c3, vpot(:,1), eigr, bec )
          CALL dforce( i, 2, c0(:,:,1,1), cdesc, fi(:,1,1), c4, c5, vpot(:,2), eigr, bec )
          !
          c2 = occup(i  , ik)* (c2 + c4)
          c3 = occup(i+1, ik)* (c3 + c5)

          IF( ttprint ) then
            !
            CALL update_lambda( i  , gam( :, :), c0(:,:,ik,1), cdesc, c2 )
            CALL update_lambda( i+1, gam( :, :), c0(:,:,ik,1), cdesc, c3 )

          END IF ! ttprint

          IF ( tsde ) THEN
             CALL wave_steepest( cp(:,i,ik,1)  , c0(:,i,ik,1)  , svar3, c2 )
             CALL wave_steepest( cp(:,i+1,ik,1), c0(:,i+1,ik,1), svar3, c3 )
          ELSE
            cp(:,i  ,ik,1) = cm(:,i  ,ik,1)
            cp(:,i+1,ik,1) = cm(:,i+1,ik,1)
            CALL wave_verlet( cp(:,i  ,ik,1), c0(:,i  ,ik,1), svar1, svar2, svar3, c2 )
            CALL wave_verlet( cp(:,i+1,ik,1), c0(:,i+1,ik,1), svar1, svar2, svar3, c3 )
          END IF

          IF( cdesc%gzero ) cp(1,i  ,ik,1)  = DBLE( cp(1,i  ,ik,1) )
          IF( cdesc%gzero ) cp(1,i+1,ik,1)  = DBLE( cp(1,i+1,ik,1) )

        END DO ! bande


        IF( MOD( n_unp, 2 ) /= 0 .and. n_unp > 1 ) THEN
          !
          nb = n_unp - 1
          !
          CALL dforce( nb, 2, c0(:,:,1,1), cdesc, fi(:,1,1), c2, c3, vpot(:,1), eigr, bec )
          CALL dforce( nb, 2, c0(:,:,1,1), cdesc, fi(:,1,2), c4, c5, vpot(:,2), eigr, bec )

          c2 = occup(nb , ik)* (c2 + c4)

          IF( ttprint ) THEN
            CALL update_lambda( nb, gam( :, :), c0(:,:,ik,1), cdesc, c2 )
          END IF

          IF ( tsde ) THEN
             CALL wave_steepest( cp(:,nb,ik,1), c0(:,nb,ik,1), svar3, c2 )
          ELSE
             cp(:,nb,ik,1) = cm(:,nb,ik,1)
             CALL wave_verlet( cp(:,nb,ik,1), c0(:,nb,ik,1), svar1, svar2, svar3, c2 )
          END IF
          IF( cdesc%gzero ) cp(1,nb,ik,1) = DBLE( cp(1,nb,ik,1) )
        END IF

        !
        CALL dforce( n_unp, 1, c0(:,:,1,1), cdesc, fi(:,1,1), c2, c3, vpot(:,1), eigr, bec )

        intermed  = -2.d0 * sum( c2 * conjg( c0(:, n_unp, ik, 1 ) ) )
        intermed3 = sum(c0(:,n_unp, ik, 1) * conjg( c0(:, n_unp, ik, 1)))

        CALL mp_sum ( intermed, group )
        CALL mp_sum ( intermed3, group )
        !  Eigenvalue of unpaired
        ei_unp_mem = intermed
        !  <Phiunpaired|Phiunpaired>
        ei_unp_wfc = intermed3
        !write(6,*) '  <psi|psi> = ', intermed3, '  ei_unp(au) = ', intermed

        IF ( tsde ) THEN
           CALL wave_steepest( cp( :, n_unp, ik, 1 ), c0( :, n_unp, ik, 1 ), svar3, c2 )
        ELSE
          cp( :, n_unp, ik, 1 ) = cm( :, n_unp, ik, 1 )
          CALL wave_verlet( cp( :, n_unp, ik, 1 ), c0( :, n_unp, ik, 1 ), svar1, svar2, svar3, c2 )
        END IF
        IF( cdesc%gzero ) cp( 1, n_unp, ik, 1 ) = DBLE( cp( 1, n_unp, ik, 1 ) )


        IF( ttprint ) THEN

            ALLOCATE( occsum( SIZE( fi, 1 ) ) )
            occsum(:) = occup(:,ik) + occdown(:,ik)

            if( cdesc%gamma ) then
               CALL eigs( nupdwn(2), gam, cgam, tortho, occsum, ei(:,ik,1), cdesc%gamma )
            else
               CALL eigs( nupdwn(2), gam, cgam, tortho, occsum, ei(:,ik,1), cdesc%gamma )
            endif
            DEALLOCATE( occsum )
            DO i = 1, nupdwn(2)
               ei( i, ik, 2 ) = ei( i , ik, 1)
            END DO

            ei( nupdwn(1), ik, 1) = ei_unp_mem
            ei( nupdwn(1), ik, 2) = 0.d0

            WRITE(6,*) 'SIC EIGENVALUES(eV), dwn and up electrons Kpoint',ik
            WRITE(6,1004) ( ei( i, ik, 2 ) * au, i = 1, nupdwn(2) )
            WRITE(6,1005) ( ei( i, ik, 1 ) * au, i = 1, nupdwn(1) )

1004        FORMAT(/,3X,'SIC EIGENVALUES DW=',3X,10F8.2)
1005        FORMAT(/,3X,'SIC EIGENVALUES UP=',3X,10F8.2)

        ENDIF

      END DO KAPPA

      s3 = cclock()
      timerd = timerd + s3 - s4

      IF( tortho ) THEN
         CALL ortho( 1, c0(:,:,:,1), cp(:,:,:,1), cdesc, pmss, emass )
      ELSE
         CALL gram( vkb, bec, nkb, cp(1,1,1,1), SIZE(cp,1), cdesc%nbt( 1 ) )
      END IF


      s4 = cclock()
      timeorto = timeorto + s4 - s3

      !  Compute fictitious kinetic energy of the electrons at time t

      ekinc(1) = cp_kinetic_energy( 1, cp(:,:,:,1), cm(:,:,:,1), cdesc, pmss, delt)
      ekinc(2) = cp_kinetic_energy( 2, cp(:,:,:,1), cm(:,:,:,1), cdesc, pmss, delt)


      DEALLOCATE( ei_t, svar3, c2, c3, c4, c5, cgam, gam, occup, occdown, STAT=ierr)
      IF( ierr /= 0 ) CALL errore(' runcp_force_pairing ', ' deallocating ', ierr)

      RETURN
    END SUBROUTINE runcp_force_pairing


!=----------------------------------------------------------------------------=!


   SUBROUTINE runcp_uspp( nfi, fccc, ccc, ema0bg, dt2bye, rhos, bec, c0, cm, &
              fromscra, restart )
     !
     use wave_base, only: wave_steepest, wave_verlet
     use control_flags, only: tbuff, lwf, tsde
     !use uspp, only : nhsa=> nkb, betae => vkb, rhovan => becsum, deeq
     use uspp, only : deeq, betae => vkb
     use reciprocal_vectors, only : gstart
     use electrons_base, only : n=>nbsp
     use wannier_subroutines, only: ef_potential
     use efield_module, only: dforce_efield, tefield, dforce_efield2, tefield2

     use gvecw, only: ngw
     !
     IMPLICIT NONE
     integer, intent(in) :: nfi
     real(8) :: fccc, ccc
     real(8) :: ema0bg(:), dt2bye
     real(8) :: rhos(:,:)
     real(8) :: bec(:,:)
     complex(8) :: c0(:,:), cm(:,:)
     logical, optional, intent(in) :: fromscra
     logical, optional, intent(in) :: restart
     !
     real(8) ::  verl1, verl2, verl3
     real(8), allocatable:: emadt2(:)
     real(8), allocatable:: emaver(:)
     complex(8), allocatable:: c2(:), c3(:)
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
           call dforce(bec,betae,i,c0(1,i),c0(1,i+1),c2,c3,rhos)
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
!==== end of loop which updates electronic degrees of freedom
!
!     buffer for wavefunctions is unit 21
!
     if(tbuff) rewind 21

   END SUBROUTINE runcp_uspp
!
!=----------------------------------------------------------------------------=!
    SUBROUTINE runcp_uspp_force_pairing( nfi, fccc, ccc, ema0bg, dt2bye, rhos, bec, c0, cm, &
                                         intermed, fromscra, restart )
  !
      USE wave_base, ONLY: wave_steepest, wave_verlet
      USE control_flags, ONLY: tbuff, lwf, tsde
  !   use uspp, only : nhsa=> nkb, betae => vkb, rhovan => becsum, deeq
      USE uspp, ONLY : deeq, betae => vkb
      USE reciprocal_vectors, ONLY : gstart
      USE wannier_subroutines, ONLY: ef_potential
      USE efield_module, ONLY: dforce_efield, tefield
  !
      USE gvecw, ONLY: ngw
  !
  !
      USE electrons_base,   ONLY: nx=>nbnd, nupdwn, iupdwn, nbspx, nbsp
      USE mp, ONLY: mp_sum 
  !
      IMPLICIT NONE
      INTEGER, INTENT(in) :: nfi
      REAL(8) :: fccc, ccc
      REAL(8) :: ema0bg(:), dt2bye
      REAL(8) :: rhos(:,:)
      REAL(8) :: bec(:,:)
      COMPLEX(8) :: c0(:,:), cm(:,:)
      LOGICAL, OPTIONAL, INTENT(in) :: fromscra
      LOGICAL, OPTIONAL, INTENT(in) :: restart
!
      REAL(8) ::  verl1, verl2, verl3
      REAL(8), ALLOCATABLE:: emadt2(:)
      REAL(8), ALLOCATABLE:: emaver(:)
      COMPLEX(8), ALLOCATABLE:: c2(:), c3(:)
      INTEGER :: i
      INTEGER :: iflag
      LOGICAL :: ttsde
!
       INTEGER    :: ierr,  nb, np_dw, is_dw, npair, n_unp, n_dwn, n_pair 
       REAL(8)    :: intermed, ei_unp_mem, ei_unp_wfc
       COMPLEX(8) ::  intermed3
       INTEGER(8), ALLOCATABLE:: occ(:)
       COMPLEX(8), ALLOCATABLE:: c4(:), c5(:)
!
! ... Controlling on sic applicanility
!
       IF( lwf ) &
            &  CALL errore('Wannier function and sic are not compatibile',1)
       IF( tefield ) &
            &  CALL errore('Electric field and sic are not implemeted',2)
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
       ALLOCATE(c2(ngw),c3(ngw),c4(ngw),c5(ngw) )
!
  !
  ! ...  set verlet variables
  !
!
       verl1 = 2.0d0 * fccc
       verl2 = 1.0d0 - verl1
       verl3 = 1.0d0 * fccc !delt * delt / pmss(1:ngw) * fccc
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
       ALLOCATE( occ(nbspx))
!
       occ(1:np_dw)  = 1
       occ(nbspx)  = 0
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
            CALL dforce(bec,betae,i,c0(1,i),c0(1,i+1),c2,c3,rhos(1,1))
            CALL dforce(bec,betae,i,c0(1,i),c0(1,i+1),c4,c5,rhos(1,2))
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
         CALL dforce(bec,betae,i,c0(1,npair),c0(1,npair+1),c2,c3,rhos(1,1))
         CALL dforce(bec,betae,i,c0(1,nbsp), c0(1,nbspx)  ,c4,c5,rhos(1,2))
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
      CALL dforce(bec,betae,n_unp,c0(1,n_unp),c0(1,nbspx),c2,c3,rhos(1,1))
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
      intermed  = -2.d0 * sum(c2 * conjg(c0(:,n_unp)))
      CALL mp_sum ( intermed )
!
!           write(6,*) 'Debug:: ei_unp(au) = ', intermed
!
       DEALLOCATE( occ )
       DEALLOCATE( emadt2 )
       DEALLOCATE( emaver )
       DEALLOCATE(c2, c4)
       DEALLOCATE(c3, c5)
!
!==== end of loop which updates electronic degrees of freedom
!
!     buffer for wavefunctions is unit 21
!
       IF(tbuff) REWIND 21

     END SUBROUTINE runcp_uspp_force_pairing

!=----------------------------------------------------------------------------=!

!=----------------------------------------------------------------------------=!
   END MODULE runcp_module
!=----------------------------------------------------------------------------=!
