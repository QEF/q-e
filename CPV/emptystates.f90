!
! Copyright (C) 2002 FPMD group
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
        LOGICAL :: tconjgrad_emp = .FALSE.

        CHARACTER(LEN=256) :: fileempty
        LOGICAL :: first = .TRUE.

#define __EMPTY_SD

        INTERFACE EMPTY
#if defined __EMPTY_CG
          MODULE PROCEDURE EMPTY_CG
#elif defined __EMPTY_CG2
          MODULE PROCEDURE EMPTY_CG2
#elif defined __EMPTY_IT
          MODULE PROCEDURE EMPTY_IT
#else
          MODULE PROCEDURE EMPTY_SD
#endif
        END INTERFACE

        PUBLIC :: empty, empty_setup, empty_print_info

! ---------------------------------------------------------------------- !
      CONTAINS
! ---------------------------------------------------------------------- !

        SUBROUTINE empty_print_info(iunit)
          USE electrons_module, ONLY: n_emp
          INTEGER, INTENT(IN) :: iunit
          IF ( n_emp > 0 ) WRITE (iunit,620) n_emp, max_emp, delt_emp
620       FORMAT(3X,'Empty states minimization : states = ',I4, &
             ' maxiter = ',I8,' delt = ',F8.4)
          RETURN
        END SUBROUTINE empty_print_info

!----------------------------------------------------------------------

        SUBROUTINE empty_setup( max_emp_inp, delt_emp_inp, ethr_emp_inp )

          INTEGER, INTENT(IN) :: max_emp_inp
          REAL(dbl), INTENT(IN) :: delt_emp_inp, ethr_emp_inp

          delt_emp  = delt_emp_inp
          emass_emp = 200.0d0
          max_emp   = max_emp_inp
          ethr_emp  = ethr_emp_inp

          tconjgrad_emp = .TRUE.
          !tconjgrad_emp = .FALSE.

          RETURN
        END SUBROUTINE empty_setup

!
!=----------------------------------------------------------------------------=!
!

      LOGICAL FUNCTION readempty( c_emp, wempt, gv )

! ...   This subroutine reads empty states from unit emptyunit

        USE wave_types, ONLY: wave_descriptor
        USE mp_global, ONLY: mpime, nproc, group, root
        USE io_global, ONLY: stdout
        USE mp, ONLY: mp_bcast
        USE mp_wave, ONLY: splitwf
        USE environment, ONLY: tscra, scradir
        USE cp_types, ONLY: recvecs

        IMPLICIT none

        COMPLEX(dbl), INTENT(INOUT) :: c_emp(:,:,:,:)
        TYPE (wave_descriptor), INTENT(IN) :: wempt
        TYPE (recvecs) :: gv

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

        IF( tscra ) THEN
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
                  CALL splitwf(c_emp(:,i,ik,ispin), ctmp, gv%ngw_l, gv%ig, mpime, nproc, root)
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

      SUBROUTINE writeempty( c_emp, wempt, gv )

! ...   This subroutine writes empty states to unit emptyunit

        USE wave_types, ONLY: wave_descriptor
        USE mp_global, ONLY: mpime, nproc, group, root
        USE mp_wave, ONLY: mergewf
        USE environment, ONLY: tscra, scradir
        USE cp_types, ONLY: recvecs

        COMPLEX(dbl), INTENT(IN) :: c_emp(:,:,:,:)
        TYPE (wave_descriptor), INTENT(IN) :: wempt
        TYPE (recvecs) :: gv

        INTEGER :: ig, i, ik, nl, ne(2), ngwm_g, nk, ispin, nspin
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
        ne        = 0
        ne( 1 : nspin ) = wempt%nbl( 1 : nspin )

        ALLOCATE( ctmp( ngwm_g ) )

        IF( tscra ) THEN
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
              CALL MERGEWF( c_emp(:,i,ik,ispin), ctmp(:), gv%ngw_l, gv%ig, mpime, nproc, root)
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
      END SUBROUTINE GRAM_EMPTY

!     
!================================================================
!
      SUBROUTINE randomizza(c_occ, wfill, ampre, c_emp, wempt, gv, kp)
        USE wave_types, ONLY: wave_descriptor
        USE reciprocal_vectors, ONLY: ig_l2g, ngw, ngwt
        USE mp_global, ONLY: mpime, nproc, root
        USE mp_wave, ONLY: splitwf
        USE brillouin, ONLY: kpoints
        USE cp_types, ONLY: recvecs
        USE control_flags, ONLY: force_pairing

! ...   Arguments
        COMPLEX(dbl), INTENT(INOUT) :: c_occ(:,:,:,:), c_emp(:,:,:,:)
        TYPE (wave_descriptor), INTENT(IN) :: wfill, wempt
        TYPE (recvecs), INTENT(IN) :: gv 
        TYPE (kpoints), INTENT(IN) :: kp
        REAL(dbl)                  :: ampre  

        REAL(dbl) :: rranf
        EXTERNAL rranf

! ...   Locals
        INTEGER   :: ig_local
        INTEGER   :: ib, ik, ispin, i, ig, ntest, j, ispin_wfc
        LOGICAL   :: tortho = .FALSE.
        REAL(dbl) :: rranf1, rranf2
        COMPLEX(dbl), ALLOCATABLE :: pwt( : )
!
        ntest = gv%ngw_g / 4

! ...   Subroutine body

! ...   initialize the wave functions in such a way that the values
! ...   of the components are independent on the number of processors

        ALLOCATE( pwt( ngwt ) )

        DO ispin = 1, wempt%nspin
          ispin_wfc = ispin
          IF( force_pairing ) ispin_wfc = 1
          DO ik = 1, wempt%nkl
            c_emp(:,:,ik,ispin) = 0.0d0
            DO ib = 1, wempt%nbl( ispin )
              pwt( 1 ) = 0.0d0
              DO ig = 2, ntest
                rranf1 = 0.5d0 - rranf()
                rranf2 = rranf()
                pwt( ig ) = ampre * CMPLX(rranf1, rranf2)
              END DO
              CALL splitwf ( c_emp( :, ib, ik, ispin ), pwt, ngw, ig_l2g, mpime, nproc, 0 )
            END DO
            IF ( .NOT. kp%gamma_only ) THEN
! ..  .       set to zero all elements outside the cutoff sphere
              DO ib = 1, wempt%nbl( ispin )
                c_emp(:,ib,ik,ispin) = c_emp(:,ib,ik,ispin) * gv%kg_mask_l(:,ik)
              END DO
            END IF
            IF ( gv%gzero ) THEN
              c_emp(1,:,ik,ispin) = (0.0d0, 0.0d0)
            END IF
            CALL gram_empty( ispin, tortho, c_occ(:,:,ik,ispin_wfc), wfill, c_emp(:,:,ik,ispin), wempt )
          END DO
        END DO

        DEALLOCATE( pwt )

        RETURN
      END SUBROUTINE RANDOMIZZA

!
!=======================================================================
!

#if defined __EMPTY_SD

    SUBROUTINE EMPTY_SD( tortho, atoms, gv, c_occ, wfill, c_emp, wempt, kp, vpot, eigr, ps)

      USE wave_types, ONLY: wave_descriptor
      USE wave_functions, ONLY: cp_kinetic_energy, crot, dft_kinetic_energy, fixwave
      USE wave_base, ONLY: hpsi, converg_base, dotp
      USE pseudopotential, ONLY: nsanl, ngh
      USE constants, ONLY: au
      USE cell_base, ONLY: tpiba2
      USE electrons_module, ONLY: pmss
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
      USE gvecw, ONLY: tecfix
      USE atoms_type_module, ONLY: atoms_type
      USE io_global, ONLY: ionode
      USE io_global, ONLY: stdout
      USE cp_types, ONLY: recvecs, pseudo, phase_factors
      USE control_flags, ONLY: force_pairing

      IMPLICIT NONE

! ... ARGUMENTS
!
      COMPLEX(dbl), INTENT(INOUT) ::  c_occ(:,:,:,:), c_emp(:,:,:,:)
      TYPE (wave_descriptor), INTENT(IN) ::  wfill, wempt
      TYPE (atoms_type), INTENT(INOUT) :: atoms ! ions structure
      TYPE (recvecs), INTENT(IN) ::  gv
      REAL (dbl), INTENT(IN) ::   vpot(:,:,:,:)
      TYPE (kpoints), INTENT(IN) :: kp
      TYPE (pseudo), INTENT(IN) :: ps
      LOGICAL, INTENT(IN) :: tortho
      TYPE (phase_factors), INTENT(IN) :: eigr
!
! ... LOCALS
!

      INTEGER   ::  i, k, j, iter, ik, nk
      INTEGER   ::  ngw, ngw_g, nspin, ispin, ispin_wfc
      INTEGER   ::  n_occ( wfill%nspin )
      INTEGER   ::  n_emp( wempt%nspin )
      INTEGER   ::  ig, iprinte, iks, nrl, jl
      REAL(dbl) ::  dek, ekinc, ekinc_old
      REAL(dbl) :: ampre
      REAL(dbl), ALLOCATABLE :: dt2bye( : )
      REAL(dbl), PARAMETER :: small    = 1.0d-14
      TYPE (projector) :: fnle( SIZE(c_emp, 3), SIZE(c_emp, 4) )
      COMPLEX(dbl), ALLOCATABLE :: eforce(:,:,:,:), cp_emp(:,:,:,:)
      REAL(dbl), ALLOCATABLE :: fi(:,:,:)
      LOGICAL       :: gamma, gzero
      LOGICAL, SAVE :: exst
!
! ... SUBROUTINE BODY
!    

      nk        = wfill%nkl
      nspin     = wfill%nspin
      ngw       = wfill%ngwl
      ngw_g     = wfill%ngwt
      n_occ     = wfill%nbt
      n_emp     = wempt%nbt
      gzero     = wfill%gzero
      gamma     = wfill%gamma
      ampre     = 0.001d0

      ekinc_old = 1.d+10
      ekinc     = 0.0d0

      CALL allocate_projector(fnle, nsanl, MAXVAL( n_emp ), ngh, gamma) 
      ALLOCATE( eforce( ngw, wempt%ldb, nk, nspin ) )
      ALLOCATE( fi( wempt%ldb, nk, nspin ) )
      ALLOCATE( cp_emp( SIZE(c_emp,1), SIZE(c_emp,2), SIZE(c_emp,3), SIZE(c_emp,4) ) )
      ALLOCATE( dt2bye( ngw ) )

      IF( ionode ) WRITE( stdout,56)

      exst = readempty( c_emp, wempt, gv )
      ! .. WRITE( stdout, * )' DEBUG empty 1 ', exst
      IF( .NOT. exst ) THEN
        CALL randomizza(c_occ, wfill, ampre, c_emp, wempt, gv, kp)
      END IF

      dt2bye = delt * delt / pmss
      cp_emp = c_emp
      fi     = 2.0d0 / nspin

      ITERATIONS: DO iter = 1, max_emp

        ekinc = 0.0d0

        SPIN_LOOP: DO ispin = 1, nspin

          ispin_wfc = ispin
          IF( force_pairing ) ispin_wfc = 1

          IF( n_emp( ispin ) < 1 ) CYCLE SPIN_LOOP 

          DO ik = 1, kp%nkpt

            CALL nlsm1( ispin, ps%wnl(:,:,:,ik), atoms, eigr%xyz, c_emp(:,:,ik,ispin), wempt, &
                gv%khg_l(:,ik), gv%kgx_l(:,:,ik), fnle(ik,ispin))

            CALL dforce_all( ispin, c_emp(:,:,:,ispin), wempt, fi(:,:,ispin), eforce(:,:,:,ispin), &
              gv, vpot(:,:,:,ispin), fnle(:,ispin), eigr, ps, ik)

! ...       Steepest descent
            DO i = 1, n_emp( ispin )
              cp_emp(:,i,ik,ispin) = c_emp(:,i,ik,ispin) +  dt2bye(:) * eforce(:, i, ik, ispin)
            END DO

            CALL fixwave( ispin, cp_emp(:,:,ik,ispin), wempt, gv%kg_mask_l(:,ik) )

            IF (tortho) THEN

              CALL ortho( ispin, c_emp(:,:,ik,ispin), cp_emp(:,:,ik,ispin), wempt, pmss, emass)
              ! WRITE( stdout,*) ' EMPTY DEBUG ', c_emp(ik,ispin)%w(4,n_emp)
              ! WRITE( stdout,*) ' EMPTY DEBUG ', cp_emp(ik,ispin)%w(4,n_emp)

              CALL gram_empty( ispin, tortho, c_occ(:,:,ik,ispin_wfc), wfill, cp_emp(:,:,ik,ispin), wempt)

            ELSE

              CALL gram_empty( ispin, .FALSE., c_occ(:,:,ik,ispin_wfc), wfill, cp_emp(:,:,ik,ispin), wempt)

            END IF

          END DO
          ekinc = ekinc + cp_kinetic_energy( ispin, cp_emp(:,:,:,ispin), c_emp(:,:,:,ispin), wempt, &
            kp, gv%kg_mask_l, pmss, delt)

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
        IF( ( ekinc / MAXVAL( n_emp ) ) <  ethr_emp ) THEN
          IF( ionode ) WRITE( stdout,112) 
          EXIT ITERATIONS
        END IF

        ekinc_old = ekinc

      END DO ITERATIONS

      CALL empty_eigs(tortho, atoms, c_emp, wempt, fi, vpot, eforce, fnle, ps, eigr, gv, kp)

      CALL writeempty( c_emp, wempt, gv )

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
    END SUBROUTINE EMPTY_SD

#endif

!=----------------------------------------------------------------------------=!
!
!   Compute the eigenvalues of the empty states 
!
!=----------------------------------------------------------------------------=!

    SUBROUTINE empty_eigs(tortho, atoms, c_emp, wempt, fi, vpot, eforce, fnle, ps, eigr, gv, kp)

      USE wave_types, ONLY: wave_descriptor
      USE wave_functions, ONLY: crot
      USE wave_base, ONLY: hpsi, dotp
      USE constants, ONLY: au
      USE cell_base, ONLY: tpiba2
      USE electrons_module, ONLY: eigs, ei_emp, pmss, emass, emp_desc
      USE forces, ONLY: dforce_all
      USE brillouin, ONLY: kpoints
      USE pseudo_projector, ONLY: projector
      USE orthogonalize, ONLY: ortho
      USE nl, ONLY: nlsm1
      USE mp, ONLY: mp_sum
      USE mp_global, ONLY: mpime, nproc, group
      USE descriptors_module, ONLY: get_local_dims, owner_of, local_index
      USE atoms_type_module, ONLY: atoms_type
      USE cp_types, ONLY: recvecs, pseudo, phase_factors

      IMPLICIT NONE

! ... ARGUMENTS

      COMPLEX(dbl), INTENT(inout) ::  c_emp(:,:,:,:)
      TYPE (wave_descriptor), INTENT(IN) ::  wempt
      TYPE(atoms_type), INTENT(INOUT) :: atoms ! ions structure
      TYPE (recvecs), INTENT(in) ::  gv
      REAL (dbl), INTENT(in) ::  vpot(:,:,:,:), fi(:,:,:)
      COMPLEX (dbl) ::  eforce(:,:,:,:)
      TYPE (kpoints), INTENT(in) :: kp
      TYPE (pseudo), INTENT(in) :: ps
      LOGICAL, INTENT(IN) :: TORTHO
      TYPE (phase_factors), INTENT(IN) :: eigr
      TYPE (projector) :: fnle(:,:)
!
! ... LOCALS
!

      INTEGER     kk, i, k, j, iopt, iter, nwh, ik, nk, ibl
      INTEGER     ngw, ngw_g, n_occ, n_emp, n_emp_l, nspin, ispin
      INTEGER     ig, iprinte, iks, nrl, jl
      LOGICAL     gamma, gzero

      REAL(dbl),    ALLOCATABLE :: gam(:,:)
      COMPLEX(dbl), ALLOCATABLE :: cgam(:,:)
      REAL(dbl),    ALLOCATABLE :: prod(:)
      COMPLEX(dbl), ALLOCATABLE :: cprod(:)

!
! ... SUBROUTINE BODY
!
      nspin     = wempt%nspin
      nk        = wempt%nkl
      ngw       = wempt%ngwl
      ngw_g     = wempt%ngwt 
      gzero     = wempt%gzero
      gamma     = wempt%gamma

      CALL get_local_dims( emp_desc, n_emp_l )

!
! ... empty state diagonalization ==

      SPIN_LOOP: DO ispin = 1, nspin

        n_emp     = wempt%nbl( ispin )

        IF( n_emp < 1 ) CYCLE SPIN_LOOP

        ALLOCATE(gam(n_emp_l, n_emp), cgam(n_emp_l, n_emp))
        ALLOCATE(prod(n_emp), cprod(n_emp))

        DO ik = 1, nk
          CALL nlsm1( ispin, ps%wnl(:,:,:,ik), atoms, eigr%xyz, c_emp(:,:,ik,ispin), wempt, &
              gv%khg_l(:,ik), gv%kgx_l(:,:,ik), fnle(ik,ispin))
        END DO

! ...   Calculate | dH / dpsi(j) >
        CALL dforce_all( ispin, c_emp(:,:,:,ispin), wempt, fi(:,:,ispin), eforce(:,:,:,ispin), gv, &
          vpot(:,:,:,ispin), fnle(:,ispin), eigr, ps)

        DO ik = 1, kp%nkpt

! ...     Calculate Eij = < psi(i) | H | psi(j) > = < psi(i) | dH / dpsi(j) >
          DO i = 1, n_emp
            IF( gamma ) THEN
              prod = hpsi( gzero, c_emp(:,:,ik, ispin), eforce(:,i,ik,ispin) )
              CALL mp_sum( prod )
              IF( mpime == owner_of( i, emp_desc, 'R' ) ) THEN
                ibl = local_index( i, emp_desc, 'R' )
                gam(ibl,:) = prod(:)
              END IF
            ELSE
              cprod = hpsi(c_emp(:,:,ik, ispin), eforce(:,i,ik,ispin) )
              CALL mp_sum( cprod )
              IF( mpime == owner_of( i, emp_desc, 'R' ) ) THEN
                ibl = local_index( i, emp_desc, 'R' )
                cgam(ibl,:) = cprod(:)
              END IF
            END IF
          END DO

          CALL eigs( n_emp, gam, cgam, tortho, fi(:,ik,ispin), ei_emp(:,ik,ispin), gamma)
        END DO

        DEALLOCATE(gam, cgam, prod, cprod)

      END DO SPIN_LOOP


      RETURN
    END SUBROUTINE


#if ! defined __EMPTY_SD

! ---------------------------------------------------------------------- !
         REAL(dbl) FUNCTION eenergy(ik, tortho, atoms, cp_emp, c_emp, wempt, hstep, hacca, eforce, &
           c_occ, wfill, vpot, fnle, ps, eigr, gv, kp)

           USE wave_types, ONLY: wave_descriptor
           USE wave_functions, ONLY: crot, fixwave
           USE wave_base, ONLY: hpsi, dotp
           USE electrons_module, ONLY: eigs, ei, ei_emp, pmss, emass, &
             emp_desc
           USE forces, ONLY: dforce_all, dforce_all2
           USE brillouin, ONLY: kpoints
           USE pseudo_projector, ONLY: projector
           USE gvecw, ONLY: tecfix
           USE orthogonalize, ONLY: ortho
           USE nl, ONLY: nlsm1
           USE mp, ONLY: mp_sum
           USE mp_global, ONLY: mpime, nproc, group
           USE descriptors_module, ONLY: get_local_dims, owner_of, local_index
           USE atoms_type_module, ONLY: atoms_type
           USE cp_types, ONLY: recvecs, pseudo, phase_factors

           IMPLICIT NONE


! ...      ARGUMENTS
           COMPLEX(dbl), INTENT(inout) ::  c_emp(:), cp_emp(:), c_occ(:)
           TYPE (wave_descriptor), INTENT(in) ::  wempt, wfill
           TYPE(atoms_type), INTENT(INOUT) :: atoms ! ions structure
           TYPE (recvecs), INTENT(in) ::  gv
           COMPLEX (dbl) ::  eforce(:,:,:)
           COMPLEX (dbl) ::  hacca(:,:,:)
           REAL(dbl) :: hstep
           TYPE (kpoints), INTENT(in) :: kp
           LOGICAL, INTENT(IN) :: TORTHO
           REAL (dbl), INTENT(in) ::  vpot(:,:,:)
           TYPE (pseudo), INTENT(in) :: ps
           TYPE (phase_factors), INTENT(IN) :: eigr
           TYPE (projector) :: fnle(:)
           INTEGER, INTENT(IN) :: ik
!
! ... LOCALS
!

           INTEGER     kk, i, k, j, iopt, iter, nwh, nk
           INTEGER     ngw, ngw_g, n_occ, n_emp, n_emp_l, nspin, ispin
           INTEGER     ib, ig, iprinte, iks, nrl, jl, ibl
           LOGICAL     gamma, gzero

           REAL(dbl),    ALLOCATABLE :: gam(:,:)
           COMPLEX(dbl), ALLOCATABLE :: cgam(:,:)
           REAL(dbl),    ALLOCATABLE :: prod(:)
           COMPLEX(dbl), ALLOCATABLE :: cprod(:)
           REAL(dbl),    POINTER :: gmod2(:)
           REAL(dbl) :: sk1

!
! ...      SUBROUTINE BODY
!
           IF( ik < 1 .OR. ik > SIZE( c_emp ) ) THEN
             CALL errore(' empty_states: eenergy ', ' ik out of bonds ', ik)
           END IF
           ngw       = SIZE(c_emp(ik)%w, 1)
           n_emp     = SIZE(c_emp(ik)%w, 2)
           ngw_g     = c_emp(ik)%ngw_g
           gzero     = c_emp(ik)%gzero
           gamma = c_emp(ik)%gamma

           CALL get_local_dims( emp_desc, n_emp_l )

           ALLOCATE(gam(n_emp_l, n_emp), cgam(n_emp_l, n_emp))
           ALLOCATE(prod(n_emp), cprod(n_emp))

           DO i = 1, n_emp
             cp_emp(ik)%w(:,i) = c_emp(ik)%w(:,i) + hstep * hacca(:, i, ik)
           END DO
           CALL fixwave( cp_emp(ik), wempt, gv%kg_mask_l(:,ik) )
           CALL gram_empty( .FALSE., c_occ(ik), wfill, cp_emp(ik), wempt )
           CALL nlsm1( ps%wnl(:,:,:,ik), atoms, eigr%xyz, cp_emp(ik), wempt, &
               gv%khg_l(:,ik), gv%kgx_l(:,:,ik), fnle(ik))

           CALL dforce_all( cp_emp, wempt, eforce(:,:,:), gv, vpot(:,:,:), fnle, eigr, ps, ik)

           ei_emp(:,ik,1) = 0.0d0
           DO i = 1, n_emp
             IF( gamma ) THEN
               prod = hpsi( gzero, cp_emp(ik)%w(:,:), eforce(:,i,ik) )
               CALL mp_sum( prod )
               IF( mpime == owner_of( i, emp_desc, 'R' ) ) THEN
                 ibl = local_index( i, emp_desc, 'R' )
                 gam(ibl,:) = prod(:)
               END IF
             ELSE
               cprod = hpsi(cp_emp(ik)%w(:,:), eforce(:,i,ik) )
               CALL mp_sum( cprod )
               IF( mpime == owner_of( i, emp_desc, 'R' ) ) THEN
                 ibl = local_index( i, emp_desc, 'R' )
                 cgam(ibl,:) = cprod(:)
               END IF
             END IF
           END DO
           CALL eigs( n_emp, gam, cgam, tortho, ei_emp(:,ik,1), gamma)
           eenergy = SUM( ( ei_emp(:,ik,1)-ei(1,ik,1) )**2 )

           DEALLOCATE(gam, cgam, prod, cprod)

         END FUNCTION

! ---------------------------------------------------------------------- !
! ---------------------------------------------------------------------- !

      SUBROUTINE EMPTY_LINMIN(ik, ispin, atoms, tbad, emin, ediff, tortho, cp_emp, c_emp, &
        wempt, c_occ, wfill, vpot, eforce, hacca, fnle, ps, eigr, gv, kp)

        USE wave_types, ONLY: wave_descriptor
        USE brillouin, ONLY: kpoints
        USE pseudo_projector, ONLY: projector
        USE mp_global, ONLY: mpime, nproc, group
        USE atoms_type_module, ONLY: atoms_type
        USE cp_types, ONLY: recvecs, pseudo, phase_factors

        IMPLICIT NONE

! ...   ARGUMENTS
        REAL(dbl) :: ediff, emin
        COMPLEX(dbl), INTENT(inout) ::  c_emp(:)
        COMPLEX(dbl), INTENT(inout) ::  cp_emp(:)
        COMPLEX(dbl), INTENT(inout) ::  c_occ(:)
        TYPE (wave_descriptor), INTENT(in) ::  wempt, wfill
        TYPE(atoms_type), INTENT(INOUT) :: atoms ! ions structure
        TYPE (recvecs), INTENT(in) ::  gv
        REAL (dbl), INTENT(in) ::  vpot(:,:,:)
        COMPLEX (dbl) ::  hacca(:,:,:)
        COMPLEX (dbl) ::  eforce(:,:,:)
        TYPE (kpoints), INTENT(in) :: kp
        TYPE (pseudo), INTENT(in) :: ps
        LOGICAL, INTENT(IN)  :: TORTHO
        LOGICAL, INTENT(OUT) :: tbad
        TYPE (phase_factors), INTENT(IN) :: eigr
        TYPE (projector) :: fnle(:)
        INTEGER, INTENT(IN) :: ik, ispin
!
! ... LOCALS
!

        REAL(dbl) :: GOLD, GLIMIT, TINY, CGOLD, ZEPS
        INTEGER   :: itmax
        PARAMETER (GOLD=1.618034D0, GLIMIT=100.D0, TINY=1.D-20)
        PARAMETER (ITMAX=20,CGOLD=.3819660D0,ZEPS=1.0D-10)

        REAL(dbl) :: ax, bx, cx, fa, fb, fc, dum, u, fu ,r, q, ulim
        REAL(dbl) :: x, p, v, w, e, fw, fv, xm, tol1, tol2, a, b, etemp, d
        REAL(dbl) :: fx, xmin, brent, eold
        LOGICAL   :: tbrent
        INTEGER   :: iter

!
! ... SUBROUTINE BODY
!
        tbrent         = .TRUE.
        tbad           = .FALSE.

        ax = 0.0d0
        bx = 1.0d0

        ! FA=FUNC(AX)
        fa = eenergy(ik, tortho, atoms, cp_emp, c_emp, wempt, ax, hacca, eforce, c_occ, wfill, &
          vpot, fnle, ps, eigr, gv, kp)

        eold = fa

        ! FB=FUNC(BX)
        fb = eenergy(ik, tortho, atoms, cp_emp, c_emp, wempt, bx, hacca, eforce, c_occ, wfill, &
          vpot, fnle, ps, eigr, gv, kp)

        !WRITE( stdout,*) ' ### fa = ', fa
        !WRITE( stdout,*) ' ### fb = ', fb

        IF(FB .GT. FA)THEN
          tbad = .TRUE.
          DUM=AX; AX=BX; BX=DUM
          DUM=FB; FB=FA; FA=DUM
        ENDIF
        CX=BX+GOLD*(BX-AX)

        ! FC=FUNC(CX)
        fc = eenergy(ik, tortho, atoms, cp_emp, c_emp, wempt, cx, hacca, eforce, c_occ, wfill, &
          vpot, fnle, ps, eigr, gv, kp)

100     IF(FB.GE.FC)THEN
          R=(BX-AX)*(FB-FC)
          Q=(BX-CX)*(FB-FA)
          U=BX-((BX-CX)*Q-(BX-AX)*R)/(2.*SIGN(MAX(ABS(Q-R),TINY),Q-R))
          ULIM=BX+GLIMIT*(CX-BX)
          IF((BX-U)*(U-CX).GT.0.)THEN
            ! FU=FUNC(U)
            fu = eenergy(ik, tortho, atoms, cp_emp, c_emp, wempt, u, hacca, eforce, c_occ, wfill, &
              vpot, fnle, ps, eigr, gv, kp)
            IF(FU.LT.FC)THEN
              AX=BX; FA=FB; BX=U; FB=FU;
              GO TO 100
            ELSE IF(FU.GT.FB)THEN
              CX=U; FC=FU;
              GO TO 100
            ENDIF
            U=CX+GOLD*(CX-BX)
            ! FU=FUNC(U)
            fu = eenergy(ik, tortho, atoms, cp_emp, c_emp, wempt, u, hacca, eforce, c_occ, wfill, &
              vpot, fnle, ps, eigr, gv, kp)
          ELSE IF((CX-U)*(U-ULIM).GT.0.)THEN
            ! FU=FUNC(U)
            fu = eenergy(ik, tortho, atoms, cp_emp, c_emp, wempt, u, hacca, eforce, c_occ, wfill, &
              vpot, fnle, ps, eigr, gv, kp)
            IF(FU.LT.FC)THEN
              BX=CX; CX=U
              U=CX+GOLD*(CX-BX)
              FB=FC; FC=FU
              ! FU=FUNC(U)
              fu = eenergy(ik, tortho, atoms, cp_emp, c_emp, wempt, u, hacca, eforce, c_occ, wfill,  &
                vpot, fnle, ps, eigr, gv, kp)
            ENDIF
          ELSE IF((U-ULIM)*(ULIM-CX).GE.0.)THEN
            U=ULIM
            ! FU=FUNC(U)
            fu = eenergy(ik, tortho, atoms, cp_emp, c_emp, wempt, u, hacca, eforce, c_occ, wfill,  &
              vpot, fnle, ps, eigr, gv, kp)
          ELSE
            U=CX+GOLD*(CX-BX)
            ! FU=FUNC(U)
            fu = eenergy(ik, tortho, atoms, cp_emp, c_emp, wempt, u, hacca, eforce, c_occ, wfill,  &
              vpot, fnle, ps, eigr, gv, kp)
          ENDIF
          AX=BX; BX=CX; CX=U; FA=FB; FB=FC; FC=FU
          GO TO 100
        ENDIF

        IF( tbrent .AND. tbad ) THEN

          IF( mpime .EQ. 0 .AND. prn_emp ) WRITE( stdout,114) ax, bx, cx, fa, fb, fc

          A=MIN(AX,CX); B=MAX(AX,CX)
          V=BX; W=V; X=V; E=0.d0
          ! FX=F(X)
          fx = eenergy(ik, tortho, atoms, cp_emp, c_emp, wempt, x, hacca, eforce, c_occ, wfill, &
                vpot, fnle, ps, eigr, gv, kp)
          FV=FX; FW=FX

          DO ITER = 1, ITMAX
            XM = 0.5d0 * (A+B)
            ! TOL1=TOL*ABS(X)+ZEPS
            TOL1 = ethr_emp * ABS(X) + ZEPS
            TOL2 = 2.d0 * TOL1
            IF(ABS(X-XM).LE.(TOL2-.5d0*(B-A))) GOTO 103
            IF(ABS(E).GT.TOL1) THEN
              R=(X-W)*(FX-FV)
              Q=(X-V)*(FX-FW)
              P=(X-V)*Q-(X-W)*R
              Q=2.d0*(Q-R)
              IF(Q.GT.0.d0) P=-P
              Q=ABS(Q)
              ETEMP=E
              E=D
              IF(ABS(P).GE.ABS(.5d0*Q*ETEMP).OR.P.LE.Q*(A-X).OR. P.GE.Q*(B-X)) GOTO 101
              D=P/Q
              U=X+D
              IF(U-A.LT.TOL2 .OR. B-U.LT.TOL2) D=SIGN(TOL1,XM-X)
              GOTO 102
            ENDIF
101         IF(X.GE.XM) THEN
              E=A-X
            ELSE
              E=B-X
            ENDIF
            D =CGOLD*E
102         IF(ABS(D).GE.TOL1) THEN
              U=X+D
            ELSE
              U=X+SIGN(TOL1,D)
            ENDIF
            ! FU=F(U)
            fu = eenergy(ik, tortho, atoms, cp_emp, c_emp, wempt, u, hacca, eforce, c_occ, wfill, &
                vpot, fnle, ps, eigr, gv, kp)
            IF(FU.LE.FX) THEN
              IF(U.GE.X) THEN
                A=X
              ELSE
                B=X
              ENDIF
              V=W; FV=FW; W=X; FW=FX; X=U; FX=FU
            ELSE
              IF(U.LT.X) THEN
                A=U
              ELSE
                B=U
              ENDIF
              IF(FU.LE.FW .OR. W.EQ.X) THEN
                V=W; FV=FW; W=U; FW=FU
              ELSE IF(FU.LE.FV .OR. V.EQ.X .OR. V.EQ.W) THEN
                V=U; FV=FU
              ENDIF
            ENDIF
          END DO
          WRITE( stdout, fmt='(" EMPTY_LINMIN, WARNING: Brent exceed maximum iterations ")' )
          ! CALL ERROR('EMPTY_LINMIN', 'Brent exceed maximum iterations.',itmax)
103       XMIN=X
          BRENT=FX
  
        ELSE

          x = bx

        END IF

        emin = eenergy(ik, tortho, atoms, cp_emp, c_emp, wempt, x, hacca, eforce, c_occ, wfill, &
          vpot, fnle, ps, eigr, gv, kp)

        IF( mpime .EQ. 0 .AND. prn_emp ) WRITE( stdout,114) ax, x, cx, fa, emin, fc

        IF( tbad ) THEN
          ediff = ABS(emin - fa)
        ELSE
          ediff = ABS(emin - eold)
        END IF

113     FORMAT(3X,'lm',I5,2X,3F22.18,2X,2F10.6)
114     FORMAT(3X,'lm',3F10.5,3F12.6)


      END SUBROUTINE

#endif

! ---------------------------------------------------------------------- !
      END MODULE empty_states
! ---------------------------------------------------------------------- !
