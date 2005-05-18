!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"

!
!     EMPTY STATES
!

! ---------------------------------------------------------------------- !
      MODULE empty_states
! ---------------------------------------------------------------------- !

        USE kinds
        USE io_files, ONLY: empty_file, emptyunit

        IMPLICIT NONE
        SAVE

        PRIVATE
 
        INTEGER   :: max_emp = 0    !  maximum number of iterations
        REAL(dbl) :: ethr_emp       !  threshold for convergence
        REAL(dbl) :: delt_emp       !  delt for the empty states updating
        REAL(dbl) :: emass_emp      !  fictitious mass for the empty states

        LOGICAL :: prn_emp       = .FALSE.

        CHARACTER(LEN=256) :: fileempty
        LOGICAL :: first = .TRUE.

        INTERFACE EMPTY
          MODULE PROCEDURE EMPTY_SD
        END INTERFACE

        PUBLIC :: empty , empty_init , empty_print_info

! ---------------------------------------------------------------------- !
      CONTAINS
! ---------------------------------------------------------------------- !

        SUBROUTINE empty_print_info(iunit)
          !
          USE electrons_module, ONLY: n_emp
          INTEGER, INTENT(IN) :: iunit
          !
          IF ( n_emp > 0 ) WRITE (iunit,620) n_emp, max_emp, delt_emp
620       FORMAT(3X,'Empty states minimization : states = ',I4, &
             ' maxiter = ',I8,' delt = ',F8.4)
          !
          RETURN
        END SUBROUTINE empty_print_info

!----------------------------------------------------------------------

        SUBROUTINE empty_init( max_emp_ , delt_emp_ , ethr_emp_ )

          INTEGER, INTENT(IN) :: max_emp_
          REAL(dbl), INTENT(IN) :: delt_emp_ , ethr_emp_

          delt_emp  = delt_emp_
          emass_emp = 200.0d0
          max_emp   = max_emp_
          ethr_emp  = ethr_emp_

          RETURN
        END SUBROUTINE empty_init

!
!=----------------------------------------------------------------------------=!
!

      LOGICAL FUNCTION readempty( c_emp, wempt )

! ...   This subroutine reads empty states from unit emptyunit

        USE wave_types, ONLY: wave_descriptor
        USE mp_global, ONLY: mpime, nproc, group, root
        USE io_global, ONLY: stdout
        USE mp, ONLY: mp_bcast
        USE mp_wave, ONLY: splitwf
        USE io_files, ONLY: scradir
        USE reciprocal_vectors, ONLY: ig_l2g

        IMPLICIT none

        COMPLEX(dbl), INTENT(INOUT) :: c_emp(:,:,:,:)
        TYPE (wave_descriptor), INTENT(IN) :: wempt

        LOGICAL :: exst
        INTEGER :: ierr, ig, i, ik, nl, ispin
        INTEGER :: ngw_rd, ne_rd(2), nk_rd, nspin_rd
        INTEGER :: nk, ne(2), ngwm_g, nspin

        COMPLEX(dbl), ALLOCATABLE :: ctmp(:)
!
! ... Subroutine Body
!

        IF( wempt%nspin < 1 .OR. wempt%nspin > 2 ) &
          CALL errore( ' readempty ', ' nspin out of range ', 1 )

        nk        = wempt%nkl
        nspin     = wempt%nspin
        ngwm_g    = wempt%ngwt 
        ne        = 0
        ne( 1 : nspin ) = wempt%nbl( 1 : nspin )

        ALLOCATE( ctmp(ngwm_g) )

        IF( LEN( TRIM( scradir ) ) > 1 ) THEN
          nl = INDEX(scradir,' ')
          fileempty = scradir(1:nl-1) // '/' // TRIM( empty_file )
        ELSE
          fileempty = TRIM( empty_file )
        END IF

        IF ( mpime == 0 ) THEN

          INQUIRE( FILE = TRIM(fileempty), EXIST = EXST )

          IF ( EXST ) THEN
            OPEN( UNIT=emptyunit, FILE=TRIM(fileempty), STATUS='OLD', FORM='UNFORMATTED' )
            READ(emptyunit) ngw_rd, ne_rd(1), ne_rd(2), nk_rd, nspin_rd
            IF( ( ne(1) > ne_rd(1) ) .OR. ( ne(2) > ne_rd(2) ) .OR. &
                (nk_rd /= nk ) .OR. (nspin_rd /= nspin) ) THEN
              EXST = .false.
              WRITE( stdout,10)  ngw_rd, ne_rd(1), ne_rd(2), nk_rd, nspin_rd
            END IF
          END IF

        END IF

 10     FORMAT('*** EMPTY STATES : wavefunctions dimensions changed  ', &
          /,'*** NGW = ', I8, ' NE1 = ', I4, ' NE2 = ', I4, ' NK = ', I4, ' NSPIN = ', I2)

        CALL mp_bcast(exst, 0, group)
        CALL mp_bcast(ne_rd, 0, group)
        CALL mp_bcast(ngw_rd, 0, group)

        IF (exst) THEN

          DO ispin = 1, nspin
            DO ik = 1, nk
              DO i = 1, ne_rd(ispin)
                IF (mpime == 0) THEN
                  READ(emptyunit) ( ctmp(ig), ig = 1, MIN( SIZE(ctmp), ngw_rd ) )
                END IF
                IF( i <= ne(ispin) ) THEN
                  CALL splitwf(c_emp(:,i,ik,ispin), ctmp, wempt%ngwl, ig_l2g, mpime, nproc, root)
                END IF
              END DO
            END DO
          END DO

        END IF

        IF (mpime == 0 .AND. EXST ) THEN
          CLOSE(emptyunit)
        END IF

        readempty = exst
        DEALLOCATE(ctmp)

        RETURN
      END FUNCTION readempty

!
!=----------------------------------------------------------------------------=!
!

      SUBROUTINE writeempty( c_emp, wempt )

! ...   This subroutine writes empty states to unit emptyunit

        USE wave_types, ONLY: wave_descriptor
        USE mp_global, ONLY: mpime, nproc, group, root
        USE mp_wave, ONLY: mergewf
        USE io_files, ONLY: scradir
        USE reciprocal_vectors, ONLY: ig_l2g

        COMPLEX(dbl), INTENT(IN) :: c_emp(:,:,:,:)
        TYPE (wave_descriptor), INTENT(IN) :: wempt

        INTEGER :: ig, i, ik, nl, ne(2), ngwm_g, nk, ispin, nspin, ngw
        LOGICAL :: exst
        COMPLEX(dbl), ALLOCATABLE :: ctmp(:)
!
! ... Subroutine Body
!
        IF( wempt%nspin < 1 .OR. wempt%nspin > 2 ) &
          CALL errore( ' writeempty ', ' nspin out of range ', 1 )

        nk        = wempt%nkl
        nspin     = wempt%nspin
        ngwm_g    = wempt%ngwt 
        ngw       = wempt%ngwl
        ne        = 0
        ne( 1 : nspin ) = wempt%nbl( 1 : nspin )

        ALLOCATE( ctmp( ngwm_g ) )

        IF( LEN( TRIM( scradir ) ) > 1  ) THEN
          nl = INDEX(scradir,' ')
          fileempty = scradir(1:nl-1) // '/' // TRIM( empty_file )
        ELSE
          fileempty = TRIM( empty_file )
        END IF

        IF( mpime == 0 ) THEN
          OPEN(UNIT=emptyunit, FILE=TRIM(fileempty), status='unknown', FORM='UNFORMATTED')
          WRITE (emptyunit)  ngwm_g, ne(1), ne(2), nk, nspin
        END IF

        DO ispin = 1, nspin
          DO ik = 1, nk
            DO i = 1, ne(ispin)
              ctmp = 0.0d0
              CALL MERGEWF( c_emp(:,i,ik,ispin), ctmp(:), ngw, ig_l2g, mpime, nproc, root)
              IF( mpime == 0 ) THEN
                WRITE (emptyunit) ( ctmp(ig), ig=1, ngwm_g )
              END IF
            END DO
          END DO
        END DO

        IF( mpime == 0 ) THEN
          CLOSE (emptyunit)
        END IF

        DEALLOCATE(ctmp)

        RETURN
      END SUBROUTINE writeempty

!
!======================================================================
!

      SUBROUTINE gram_empty( ispin, tortho, cf, wfill, ce,  wempt )

! This subroutine orthogonalize the empty states CE to the
! filled states CF using gram-shmitd . 
! If TORTHO is FALSE the subroutine orthonormalizes the 
! empty states CE and orthogonalize them to the CF states.

      USE wave_types, ONLY: wave_descriptor
      USE mp, ONLY: mp_sum
      USE mp_global, ONLY: nproc, mpime, group

      REAL(dbl) SQRT, DNRM2
      
! ... ARGUMENTS
      LOGICAL, INTENT(IN) :: TORTHO
      COMPLEX(dbl), INTENT(INOUT) :: CF(:,:), CE(:,:)
      type (wave_descriptor), INTENT(IN) :: wfill, wempt
      INTEGER, INTENT(IN) :: ispin

! ... LOCALS
      INTEGER      :: i, j, ig, NF, NE, NGW, ldw
      REAL(dbl)    :: ANORM
      REAL(dbl)   , ALLOCATABLE :: SF(:),  SE(:),  TEMP(:)
      COMPLEX(dbl), ALLOCATABLE :: CSF(:), CSE(:), CTEMP(:)
      COMPLEX(dbl) :: czero, cone, cmone
!
! ... SUBROUTINE BODY
!

      NF    = wfill%nbl( ispin )
      NE    = wempt%nbl( ispin )
      NGW   = wfill%ngwl
      czero  = 0.0d0
      cone   = 1.0d0
      cmone  = -1.0d0

      IF( wfill%gamma ) THEN

        ldw   = 2 * SIZE( cf, 1 )
        ALLOCATE( SF(nf), SE(ne), TEMP(2*ngw) )
        DO I = 1, NE
          SF = 0.D0
          IF( wfill%gzero ) THEN
            CALL DAXPY( nf, - REAL( ce(1,i) ), cf(1,1), 2*ngw, sf, 1 )
          END IF
          CALL DGEMV( 'T', 2*ngw, nf, 2.0d0, cf(1,1), ldw, ce(1,i), 1, 1.d0, sf, 1 )
          CALL mp_sum( SF, group )
          temp = 0.0d0  
          CALL DGEMV( 'N', 2*ngw, nf, 1.d0, cf(1,1), ldw, sf, 1, 1.d0, TEMP, 1 )
          IF(.NOT.TORTHO) THEN
            IF( I > 1 ) THEN
              SE = 0.D0
              IF( wfill%gzero ) THEN
                CALL DAXPY( i-1, -REAL( ce(1,i) ), ce(1,1), 2*ngw, se, 1 )
              END IF
              CALL DGEMV( 'T', 2*ngw, i-1, 2.0d0, ce(1,1), ldw, ce(1,i), 1, 1.d0, se, 1 )
              CALL mp_sum( SE, group )
              CALL DGEMV( 'N', 2*ngw, i-1, 1.d0, ce(1,1), ldw, se, 1, 1.d0, temp, 1 )
            END IF
          END IF
          CALL DAXPY(2*NGW,-1.D0,TEMP,1,CE(1,I),1)
          IF(.NOT.TORTHO) THEN
            IF(wempt%gzero) THEN
              ANORM = DNRM2(2*(NGW-1),CE(2,I),1)
              ANORM=2.D0*ANORM*ANORM+CE(1,I)*CE(1,I)
            ELSE
              ANORM = DNRM2(2*NGW,CE(1,I),1)
              ANORM=2.D0*ANORM*ANORM
            END IF
            CALL mp_sum(ANORM,group)
            ANORM=SQRT(ANORM)
            CALL DSCAL(2*NGW,1.0D0/ANORM,CE(1,I),1)
          END IF
        END DO
        DEALLOCATE( SF, SE, TEMP )
    
      ELSE

        ldw   = SIZE( cf, 1 )
        ALLOCATE( CSF(nf), CSE(ne), CTEMP(ngw) )
        DO i = 1, NE
          csf   = 0.0d0
          CALL ZGEMV('C', ngw, nf, cone, cf(1,1), ldw, &
            ce(1,i), 1, czero, csf(1), 1)
          CALL mp_sum(csf, group)
          CALL ZGEMV('N', ngw, nf, cone, cf(1,1), ldw, &
            csf(1), 1, czero, ctemp, 1)
          IF( .NOT. TORTHO ) THEN
            cse   = 0.0d0
            IF( i .GT. 1 ) THEN
              CALL ZGEMV('C', ngw, i-1, cone, ce(1,1), ldw, &
                ce(1,i), 1, czero, cse(1), 1)
              CALL mp_sum(cse, group)
              CALL ZGEMV('N', ngw, i-1, cone, ce(1,1), ldw, &
                cse(1), 1, cone, ctemp, 1)
            END IF
          END IF
          CALL ZAXPY(ngw, cmone, ctemp, 1, ce(1,i), 1)
          IF(.NOT.TORTHO) THEN
            anorm = 0.0d0
            do ig = 1, ngw
              anorm = anorm + REAL( ce(ig,i) * CONJG(ce(ig,i)) )
            enddo
            CALL mp_sum(anorm,group)
            anorm = 1.0d0 / MAX( sqrt(anorm), 1.d-14 )
            CALL ZDSCAL(NGW, anorm, CE(1,i), 1)
          END IF
        END DO
        DEALLOCATE( CSF, CSE, CTEMP )
    
      END IF

      RETURN
      END SUBROUTINE gram_empty

!     
!================================================================
!
      SUBROUTINE sodomizza(c_occ, wfill, ampre, c_emp, wempt, kp)
        USE wave_types, ONLY: wave_descriptor
        USE reciprocal_vectors, ONLY: ig_l2g
        USE mp_global, ONLY: mpime, nproc, root
        USE mp_wave, ONLY: splitwf
        USE brillouin, ONLY: kpoints
        USE control_flags, ONLY: force_pairing
        USE reciprocal_space_mesh, ONLY: gkmask_l
        USE wave_functions, ONLY: wave_rand_init

! ...   Arguments
        COMPLEX(dbl), INTENT(INOUT) :: c_occ(:,:,:,:), c_emp(:,:,:,:)
        TYPE (wave_descriptor), INTENT(IN) :: wfill, wempt
        TYPE (kpoints), INTENT(IN) :: kp
        REAL(dbl)                  :: ampre  

        REAL(dbl) :: rranf
        EXTERNAL rranf

! ...   Locals
        INTEGER   :: ig_local
        INTEGER   :: ngw, ngwt
        INTEGER   :: ib, ik, ispin, ispin_wfc
        LOGICAL   :: tortho = .FALSE.
        COMPLEX(dbl), ALLOCATABLE :: pwt( : )
!
! ...   Subroutine body

! ...   initialize the wave functions in such a way that the values
! ...   of the components are independent on the number of processors

        ngwt  = wfill%ngwt
        ngw   = wfill%ngwl

        ALLOCATE( pwt( ngwt ) )

        DO ispin = 1, wempt%nspin
          ispin_wfc = ispin
          IF( force_pairing ) ispin_wfc = 1
          DO ik = 1, wempt%nkl
            CALL wave_rand_init( c_emp( :, :, ik, ispin ) )
            IF ( .NOT. kp%gamma_only ) THEN
              ! ..  set to zero all elements outside the cutoff sphere
              DO ib = 1, wempt%nbl( ispin )
                c_emp(:,ib,ik,ispin) = c_emp(:,ib,ik,ispin) * gkmask_l(:,ik)
              END DO
            END IF
            IF ( wempt%gzero ) THEN
              c_emp(1,:,ik,ispin) = (0.0d0, 0.0d0)
            END IF
            CALL gram_empty( ispin, tortho, c_occ(:,:,ik,ispin_wfc), wfill, c_emp(:,:,ik,ispin), wempt )
          END DO
        END DO

        DEALLOCATE( pwt )

        RETURN
      END SUBROUTINE sodomizza

!
!=======================================================================
!

    SUBROUTINE empty_sd( tortho, atoms, c_occ, wfill, c_emp, wempt, kp, vpot, eigr, ps)

      USE wave_types, ONLY: wave_descriptor
      USE wave_functions, ONLY: cp_kinetic_energy, crot, dft_kinetic_energy, fixwave
      USE wave_base, ONLY: hpsi, converg_base, dotp
      USE pseudopotential, ONLY: nsanl
      USE constants, ONLY: au
      USE cell_base, ONLY: tpiba2
      USE electrons_module, ONLY: pmss, n_emp
      USE cp_electronic_mass, ONLY: emass
      USE time_step, ONLY: delt
      USE forces, ONLY: dforce_all
      USE brillouin, ONLY: kpoints
      USE pseudo_projector, ONLY: allocate_projector, projector, &
        deallocate_projector
      USE orthogonalize, ONLY: ortho
      USE nl, ONLY: nlsm1
      USE mp, ONLY: mp_sum
      USE mp_global, ONLY: mpime, nproc, group
      USE check_stop, ONLY: check_stop_now
      USE atoms_type_module, ONLY: atoms_type
      USE io_global, ONLY: ionode
      USE io_global, ONLY: stdout
      USE cp_types, ONLY: pseudo
      USE control_flags, ONLY: force_pairing, gamma_only
      USE reciprocal_space_mesh, ONLY: gkmask_l, gkx_l, gk_l
      USE reciprocal_vectors, ONLY: g, gx
      USE uspp_param, ONLY: nhm

      IMPLICIT NONE

! ... ARGUMENTS
!
      COMPLEX(dbl), INTENT(INOUT) ::  c_occ(:,:,:,:), c_emp(:,:,:,:)
      TYPE (wave_descriptor), INTENT(IN) ::  wfill, wempt
      TYPE (atoms_type), INTENT(INOUT) :: atoms ! ions structure
      REAL (dbl), INTENT(IN) ::   vpot(:,:,:,:)
      TYPE (kpoints), INTENT(IN) :: kp
      TYPE (pseudo), INTENT(IN) :: ps
      LOGICAL, INTENT(IN) :: tortho
      COMPLEX(dbl) :: eigr(:,:)
!
! ... LOCALS
!

      INTEGER   ::  i, k, j, iter, ik, nk
      INTEGER   ::  nspin, ispin, ispin_wfc
      INTEGER   ::  n_occ( wfill%nspin )
      INTEGER   ::  ig, iprinte, iks, nrl, jl, ngw
      REAL(dbl) ::  dek, ekinc, ekinc_old
      REAL(dbl) :: ampre
      REAL(dbl), ALLOCATABLE :: dt2bye( : )
      REAL(dbl), PARAMETER :: small    = 1.0d-14
      TYPE (projector) :: fnle( SIZE(c_emp, 3), SIZE(c_emp, 4) )
      COMPLEX(dbl), ALLOCATABLE :: eforce(:,:,:,:), cp_emp(:,:,:,:)
      REAL(dbl), ALLOCATABLE :: fi(:,:,:)
      LOGICAL       :: gamma
      LOGICAL, SAVE :: exst
!
! ... SUBROUTINE BODY
!    

      nk        = wfill%nkl
      nspin     = wfill%nspin
      ngw       = wfill%ngwl
      n_occ     = wfill%nbt
      gamma     = wfill%gamma
      ampre     = 0.001d0

      ekinc_old = 1.d+10
      ekinc     = 0.0d0

      CALL allocate_projector(fnle, nsanl, n_emp, nhm, gamma) 
      ALLOCATE( eforce( ngw, wempt%ldb, nk, nspin ) )
      ALLOCATE( fi( wempt%ldb, nk, nspin ) )
      ALLOCATE( cp_emp( SIZE(c_emp,1), SIZE(c_emp,2), SIZE(c_emp,3), SIZE(c_emp,4) ) )
      ALLOCATE( dt2bye( ngw ) )

      IF( ionode ) WRITE( stdout,56)

      exst = readempty( c_emp, wempt )
      ! .. WRITE( stdout, * )' DEBUG empty 1 ', exst
      IF( .NOT. exst ) THEN
        CALL sodomizza(c_occ, wfill, ampre, c_emp, wempt, kp)
      END IF

      dt2bye = delt * delt / pmss
      cp_emp = c_emp
      fi     = 2.0d0 / nspin

      ITERATIONS: DO iter = 1, max_emp

        ekinc = 0.0d0

        SPIN_LOOP: DO ispin = 1, nspin

          ispin_wfc = ispin
          IF( force_pairing ) ispin_wfc = 1

          IF( n_emp < 1 ) CYCLE SPIN_LOOP 

          DO ik = 1, kp%nkpt

            IF( gamma_only ) THEN
              CALL nlsm1( ispin, ps%wnl(:,:,:,ik), atoms, eigr, c_emp(:,:,ik,ispin), wempt, &
                g(:), gx(:,:), fnle(ik,ispin))
            ELSE
              CALL nlsm1( ispin, ps%wnl(:,:,:,ik), atoms, eigr, c_emp(:,:,ik,ispin), wempt, &
                gk_l(:,ik), gkx_l(:,:,ik), fnle(ik,ispin))
            END IF

            CALL dforce_all( ispin, c_emp(:,:,:,ispin), wempt, fi(:,:,ispin), eforce(:,:,:,ispin), &
              vpot(:,:,:,ispin), fnle(:,ispin), eigr, ps, ik)

! ...       Steepest descent
            DO i = 1, n_emp
              cp_emp(:,i,ik,ispin) = c_emp(:,i,ik,ispin) +  dt2bye(:) * eforce(:, i, ik, ispin)
            END DO

            CALL fixwave( ispin, cp_emp(:,:,ik,ispin), wempt, gkmask_l(:,ik) )

            IF (tortho) THEN

              CALL ortho( ispin, c_emp(:,:,ik,ispin), cp_emp(:,:,ik,ispin), wempt, pmss, emass)

              CALL gram_empty( ispin, tortho, c_occ(:,:,ik,ispin_wfc), wfill, cp_emp(:,:,ik,ispin), wempt)

            ELSE

              CALL gram_empty( ispin, .FALSE., c_occ(:,:,ik,ispin_wfc), wfill, cp_emp(:,:,ik,ispin), wempt)

            END IF

          END DO
          ekinc = ekinc + cp_kinetic_energy( ispin, cp_emp(:,:,:,ispin), c_emp(:,:,:,ispin), wempt, kp, pmss, delt)

        END DO SPIN_LOOP

        dek = ekinc - ekinc_old
        IF( ionode ) WRITE( stdout,113) ITER, dek, ekinc

        c_emp = cp_emp

! ...   check for exit
!
        IF ( check_stop_now() ) THEN
          EXIT ITERATIONS
        END IF

! ...   check for convergence
!
        IF( ( ekinc / n_emp ) <  ethr_emp ) THEN
          IF( ionode ) WRITE( stdout,112) 
          EXIT ITERATIONS
        END IF

        ekinc_old = ekinc

      END DO ITERATIONS

      CALL empty_eigs(tortho, atoms, c_emp, wempt, fi, vpot, eforce, fnle, ps, eigr, kp)

      CALL writeempty( c_emp, wempt )

      CALL deallocate_projector(fnle)
      DEALLOCATE( eforce )
      DEALLOCATE( cp_emp )
      DEALLOCATE( fi )
      DEALLOCATE( dt2bye )
              
 55   FORMAT(1X,I8,4F12.6)
 56   FORMAT(/,3X,'Empty states minimization starting ')
111   FORMAT(I5,2X,F12.8)
113   FORMAT(I5,2X,F22.18,F22.18)
112   FORMAT(/,3X,'Empty states: convergence achieved') 

    RETURN
    END SUBROUTINE empty_sd

!=----------------------------------------------------------------------------=!
!
!   Compute the eigenvalues of the empty states 
!
!=----------------------------------------------------------------------------=!

    SUBROUTINE empty_eigs(tortho, atoms, c_emp, wempt, fi, vpot, eforce, fnle, ps, eigr, kp)

      USE wave_types, ONLY: wave_descriptor
      USE wave_functions, ONLY: crot
      USE wave_constrains, ONLY: update_lambda
      USE wave_base, ONLY: hpsi, dotp
      USE constants, ONLY: au
      USE control_flags, ONLY: gamma_only
      USE cell_base, ONLY: tpiba2
      USE electrons_module, ONLY: eigs, ei_emp, pmss, emass, n_emp_l, n_emp
      USE forces, ONLY: dforce_all
      USE brillouin, ONLY: kpoints
      USE pseudo_projector, ONLY: projector
      USE orthogonalize, ONLY: ortho
      USE nl, ONLY: nlsm1
      USE mp, ONLY: mp_sum
      USE mp_global, ONLY: mpime, nproc, group
      USE atoms_type_module, ONLY: atoms_type
      USE cp_types, ONLY: pseudo
      USE reciprocal_space_mesh, ONLY: gkx_l, gk_l
      USE reciprocal_vectors,    ONLY: gx, g

      IMPLICIT NONE

! ... ARGUMENTS

      COMPLEX(dbl), INTENT(inout) ::  c_emp(:,:,:,:)
      TYPE (wave_descriptor), INTENT(IN) ::  wempt
      TYPE(atoms_type), INTENT(INOUT) :: atoms ! ions structure
      REAL (dbl), INTENT(in) ::  vpot(:,:,:,:), fi(:,:,:)
      COMPLEX (dbl) ::  eforce(:,:,:,:)
      TYPE (kpoints), INTENT(in) :: kp
      TYPE (pseudo), INTENT(in) :: ps
      LOGICAL, INTENT(IN) :: TORTHO
      COMPLEX(dbl) :: eigr(:,:)
      TYPE (projector) :: fnle(:,:)
!
! ... LOCALS
!

      INTEGER     kk, i, k, j, iopt, iter, nwh, ik, nk, ibl
      INTEGER     ngw, n_occ, nspin, ispin
      INTEGER     ig, iprinte, iks, nrl, jl
      LOGICAL     gamma

      REAL(dbl),    ALLOCATABLE :: gam(:,:)
      COMPLEX(dbl), ALLOCATABLE :: cgam(:,:)

!
! ... SUBROUTINE BODY
!
      nspin     = wempt%nspin
      nk        = wempt%nkl
      ngw       = wempt%ngwl
      gamma     = wempt%gamma

!
! ... empty state diagonalization ==

      SPIN_LOOP: DO ispin = 1, nspin

        IF( n_emp < 1 ) CYCLE SPIN_LOOP

        ALLOCATE(gam(n_emp_l(ispin), n_emp), cgam(n_emp_l(ispin), n_emp))

        DO ik = 1, nk
          IF( gamma_only ) THEN
            CALL nlsm1( ispin, ps%wnl(:,:,:,ik), atoms, eigr, c_emp(:,:,ik,ispin), wempt, &
              g(:), gx(:,:), fnle(ik,ispin))
          ELSE
            CALL nlsm1( ispin, ps%wnl(:,:,:,ik), atoms, eigr, c_emp(:,:,ik,ispin), wempt, &
              gk_l(:,ik), gkx_l(:,:,ik), fnle(ik,ispin))
          END IF
        END DO

! ...   Calculate | dH / dpsi(j) >
        CALL dforce_all( ispin, c_emp(:,:,:,ispin), wempt, fi(:,:,ispin), eforce(:,:,:,ispin), &
          vpot(:,:,:,ispin), fnle(:,ispin), eigr, ps)

        DO ik = 1, kp%nkpt

! ...     Calculate Eij = < psi(i) | H | psi(j) > = < psi(i) | dH / dpsi(j) >
          DO i = 1, n_emp
            IF( gamma ) THEN
              CALL update_lambda( i, gam, c_emp(:,:,ik, ispin), wempt, eforce(:,i,ik,ispin) )
            ELSE
              CALL update_lambda( i, cgam, c_emp(:,:,ik, ispin), wempt, eforce(:,i,ik,ispin) )
            END IF
          END DO

          CALL eigs( n_emp, gam, cgam, tortho, fi(:,ik,ispin), ei_emp(:,ik,ispin), gamma)
        END DO

        DEALLOCATE(gam, cgam)

      END DO SPIN_LOOP


      RETURN
    END SUBROUTINE empty_eigs

! ---------------------------------------------------------------------- !
      END MODULE empty_states
! ---------------------------------------------------------------------- !
