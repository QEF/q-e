!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"

! ---------------------------------------------------------------------- !
      MODULE runcg_module
! ---------------------------------------------------------------------- !

        USE kinds

        IMPLICIT NONE
        SAVE

        PRIVATE

        LOGICAL :: tforce  = .FALSE.
        LOGICAL :: tstress = .FALSE.


        INTEGER :: nsteep = 5
        REAL(DP) :: cg_dt = 4.0d0
        REAL(DP) :: cg_dt2 = 25.0d0
        REAL(DP) :: cg_emass = 200.d0
        LOGICAL :: cg_prn = .FALSE.
        REAL(DP) :: old_clock_value = -1.0d0

        INTERFACE runcg
          MODULE PROCEDURE runcg_new
        END INTERFACE

        REAL(DP), EXTERNAL :: cclock

        PUBLIC :: runcg, runcg_info


! ---------------------------------------------------------------------- !
      CONTAINS
! ---------------------------------------------------------------------- !

        SUBROUTINE runcg_info( unit )
          INTEGER, INTENT(IN) :: unit
 100      FORMAT(/,3X,'Using Conjugate Gradient for electronic minimization')
          RETURN
        END SUBROUTINE runcg_info



!  -----------------------------------------------------------------------
!  BEGIN manual

   SUBROUTINE runcg_new(tortho, tprint, rhoe, atoms_0, &
                bec, becdr, eigr, ei1, ei2, ei3, sfac, c0, cm, cp, cdesc, tcel, ht0, occ, ei, &
                vpot, doions, edft, maxnstep, cgthr, tconv )

!  this routine computes the electronic ground state via ...
!  END manual

! ... declare modules
      USE energies, ONLY: dft_energy_type, print_energies
      USE electrons_module, ONLY: pmss, eigs, nb_l
      USE cp_electronic_mass, ONLY: emass
      USE wave_functions, ONLY: cp_kinetic_energy, proj, fixwave
      USE wave_base, ONLY: dotp, hpsi
      USE wave_constrains, ONLY: update_lambda
      USE check_stop, ONLY: check_stop_now
      USE forces
      USE io_global, ONLY: ionode
      USE io_global, ONLY: stdout
      USE orthogonalize
      USE cell_module, ONLY: boxdimensions
      USE wave_types
      USE potentials, ONLY: kspotential
      USE time_step, ONLY: delt
      USE atoms_type_module, ONLY: atoms_type
      USE control_flags, ONLY: force_pairing
      USE environment, ONLY: start_cclock_val
      USE reciprocal_space_mesh, ONLY: gkmask_l
      USE uspp,             ONLY : vkb, nkb

      IMPLICIT NONE

! ... declare subroutine arguments
      LOGICAL   :: tortho, tprint, tcel, doions, tconv
      TYPE (atoms_type) :: atoms_0
      COMPLEX(DP), INTENT(INOUT) :: c0(:,:,:,:), cm(:,:,:,:), cp(:,:,:,:)
      TYPE (wave_descriptor) :: cdesc
      REAL(DP) :: rhoe(:,:)
      COMPLEX(DP) :: eigr(:,:)
      COMPLEX(DP) :: ei1(:,:)
      COMPLEX(DP) :: ei2(:,:)
      COMPLEX(DP) :: ei3(:,:)
      COMPLEX(DP) :: sfac(:,:)
      TYPE (boxdimensions), INTENT(INOUT) :: ht0
      REAL(DP) :: occ(:,:,:)
      REAL(DP) :: bec(:,:)
      REAL(DP) :: becdr(:,:,:)
      TYPE (dft_energy_type) :: edft
      INTEGER :: maxnstep
      REAL(DP) :: cgthr

      REAL(DP)    :: ei(:,:,:)
      REAL(DP)    :: vpot(:,:)

! ... declare other variables
      LOGICAL :: ttsde, ttprint, ttforce, ttstress, gzero
      REAL(DP) :: timepre, s0, s1, s2, s3, s4, s5, s6, seconds_per_iter
      REAL(DP) :: dene, eold, timerd, timeorto, ekinc
      COMPLEX(DP), ALLOCATABLE :: cgam(:,:)
      REAL(DP),    ALLOCATABLE :: gam(:,:)
      REAL(DP), ALLOCATABLE :: dt2bye( : )


      REAL(DP)    :: gg, ggo, ekinc_old, emin, demin, dek, dt2fact
      COMPLEX(DP) :: lambda

      COMPLEX(DP), ALLOCATABLE :: hacca(:,:,:,:)

      INTEGER :: ib, ibl, ik, ispin, ngw, nfi_l, nspin, isteep, i
      INTEGER :: nk, iter, ierr
      LOGICAL :: gamma_symmetry
      LOGICAL :: tbad
      INTEGER :: nb  ( cdesc%nspin )
      INTEGER :: nb_g( cdesc%nspin )

! ... end of declarations
!  ----------------------------------------------

      nk          = cdesc%nkl
      nspin       = cdesc%nspin
      doions      = .FALSE.
      eold        = 1.0d10  ! a large number
      timerd      = 0
      timeorto    = 0
      isteep      = nsteep
      ttsde       = .TRUE.
      ttprint     = .FALSE.
      ttforce     = .FALSE.
      ttstress    = .FALSE.
      gzero       = cdesc%gzero
      gamma_symmetry = cdesc%gamma
      tbad        = .FALSE.
      dt2fact     = 1.0d0

      ngw         = cdesc%ngwl
      nb          = cdesc%nbl

      IF( force_pairing ) &
        CALL errore( ' runcg ', ' force pairing not implemented ', 1 )

      ALLOCATE(hacca( ngw, MAXVAL( nb ), nk, nspin ), STAT=ierr )
      IF( ierr/=0 ) CALL errore(' runcg ', ' allocating hacca ',ierr)

      ALLOCATE( dt2bye( ngw ) )
      dt2bye = delt * delt / pmss

      WRITE(stdout,100) cgthr, maxnstep
 100  FORMAT(/,3X,'Using Conjugate Gradient for electronic minimization', &
             /,3X,'energy threshold ........... = ',1D10.4, &
             /,3X,'maximum number of iterations = ', 1I6 )


      IF(ionode) THEN
        WRITE( stdout,'(/,3X,"Conjugate Gradient Optimizations, starting ...")' )
        WRITE( stdout,'(/,3X,"iter     erho          derho       ekinc      seconds")' )
      END IF

      CONJUGATE_GRADIENTS: DO iter = 1, maxnstep

        s1 = cclock()

        CALL kspotential( 1, ttprint, ttforce, ttstress, rhoe, &
          atoms_0, bec, becdr, eigr, ei1, ei2, ei3, sfac, c0, cdesc, tcel, ht0, occ, vpot, edft, timepre )

        s2 = cclock()

        DO ispin = 1, nspin

! ...     Calculate wave functions gradient (temporarely stored in cp)
! ...     |d H / dPsi_j > = H |Psi_j> - Sum{i} <Psi_i|H|Psi_j> |Psi_i>

          CALL dforce_all( ispin, c0(:,:,1,ispin), cdesc, occ(:,1,ispin), cp(:,:,1,ispin), &
            vpot(:,ispin), eigr, bec )
 
! ...     Project the gradient
          IF( gamma_symmetry ) THEN
            CALL proj( ispin, cp(:,:,1,ispin), cdesc, c0(:,:,1,ispin), cdesc )
          ELSE
            DO ik = 1, nk
              CALL proj( ispin, ik, cp(:,:,:,ispin), cdesc, c0(:,:,:,ispin), cdesc)
            END DO
          END IF
        END DO

        s3 = cclock()

! ...   Calculate new direction hacca for the line minimization
        DO ispin = 1, nspin
          DO ik = 1, nk
            DO i = 1, nb( ispin )
              cp(:, i, ik, ispin) = cp(:, i, ik, ispin) * dt2bye(:) * dt2fact
              IF( iter > 1 ) THEN
                IF( gamma_symmetry ) THEN
                  ggo = dotp( gzero,  cm(:, i, ik, ispin), cm(:, i, ik, ispin) )
                ELSE
                  ggo = dotp( cm(:, i, ik, ispin), cm(:, i, ik, ispin) )
                END IF
                cm(:, i, ik, ispin) = cp(:, i, ik, ispin) - cm(:, i, ik, ispin)
                IF( gamma_symmetry ) THEN
                  gg  = dotp( gzero,  cm(:, i, ik, ispin), cp(:, i, ik, ispin))
                ELSE
                  gg  = dotp( cm(:, i, ik, ispin), cp(:, i, ik, ispin))
                END IF
                lambda = gg / ggo
                hacca(:, i, ik, ispin) = cp(:, i, ik, ispin) + lambda * hacca(:, i, ik, ispin)
              ELSE
                hacca(:, i, ik, ispin) = cp(:, i, ik, ispin)
              END IF
            END DO
          END DO
        END DO

        !  save the gradient in "cm" for the next iteration

        cm = cp

        s4 = cclock()

        !  perform line minimization in the direction of "hacca"

        CALL CGLINMIN(emin, demin, tbad, edft, cp, c0, cdesc, occ, vpot, rhoe, hacca, &
          atoms_0, ht0, bec, becdr, eigr, ei1, ei2, ei3, sfac)

        ! CALL print_energies( edft )
        s5 = cclock()

        IF( tbad ) THEN

          !  if we find a bad direction slow down the move and ...

          IF( ionode ) WRITE( stdout, fmt='(3X,"bad step, advancing with steepest descent")')
          dt2fact = dt2fact * 0.5d0  

          !  ... with the up to date gradient "cm" perform a steepest descent step

          cp = c0 + cm

          CALL fixwave( cp, cdesc, gkmask_l )

          IF( tortho ) THEN
             CALL ortho( c0, cp, cdesc, pmss, emass )
          ELSE
             DO ispin = 1, nspin
                CALL gram( vkb, bec, nkb, cp(1,1,1,ispin), SIZE(cp,1), cdesc%nbt( ispin ) )
             END DO
          END IF

        END IF

        ekinc = 0.0d0
        DO ispin = 1, nspin
          ekinc = ekinc + cp_kinetic_energy( ispin, cp(:,:,:,ispin), c0(:,:,:,ispin), cdesc, pmss, delt)
        END DO
        IF( iter > 1 ) THEN
          dek   = ekinc - ekinc_old
        ELSE
          dek   = 1.0d0
        END IF

        IF( old_clock_value < 0.0d0 ) old_clock_value = start_cclock_val
        s0 = cclock()
        seconds_per_iter = (s0 - old_clock_value)
        old_clock_value = s0

        IF( ionode ) THEN
          WRITE( stdout,113) iter, emin, demin, ekinc, seconds_per_iter
113       FORMAT(1X,I5,2X,F14.6,2X,3D12.4)
        END IF

        c0 = cp

        s6 = cclock()

        ! ...   check for exit

        IF (check_stop_now()) THEN
          doions = .FALSE.
          tconv  = .FALSE.
          EXIT CONJUGATE_GRADIENTS
        END IF
        !
        IF( ABS( demin ) / MAXVAL( nb ) < cgthr ) THEN
          !
          IF(ionode) WRITE( stdout,*) "  convergence achieved successfully"
          !
          doions = .TRUE.
          tconv  = .TRUE.
          EXIT CONJUGATE_GRADIENTS
        END IF
        !
        ekinc_old = ekinc
        !
      END DO CONJUGATE_GRADIENTS

      !  set wave functions velocity to 0
      cm = c0

      IF( (iter .GT. maxnstep) .AND. ionode) THEN
        WRITE( stdout,*) "  convergence not achieved"
        WRITE( stdout,*) "  maximum number of iteration exceeded"
      END IF

      IF( tprint ) THEN
        DO ispin = 1, nspin

          CALL dforce_all( ispin, c0(:,:,1,ispin), cdesc, occ(:,1,ispin), hacca(:,:,1,ispin), &
            vpot(:,ispin), eigr, bec )

          nb_g( ispin ) = cdesc%nbt( ispin )

          IF( gamma_symmetry ) THEN
            ALLOCATE(cgam(1,1), gam( nb_l( ispin ), nb_g( ispin ) ), STAT=ierr)
          ELSE
            ALLOCATE(cgam(nb_l( ispin ),nb_g( ispin )), gam(1,1), STAT=ierr)
          END IF
          IF( ierr/=0 ) CALL errore(' runcg ', ' allocating gam ',ierr)
          DO ik = 1, nk
            DO i = 1, nb( ispin )
              IF( gamma_symmetry ) THEN
                CALL update_lambda( i,  gam, c0(:,:,ik,ispin), cdesc, hacca(:,i,ik,ispin) )
              ELSE
                CALL update_lambda( i, cgam, c0(:,:,ik,ispin), cdesc, hacca(:,i,ik,ispin) )
              END IF
            END DO
            CALL eigs( nb( ispin ), gam, cgam, tortho, occ(:,ik,ispin), ei(:,ik,ispin), gamma_symmetry)
          END DO
          DEALLOCATE( cgam, gam, STAT=ierr )
          IF( ierr/=0 ) CALL errore(' runcg ', ' deallocating gam ',ierr)
        END DO
      END IF

      DEALLOCATE( hacca, STAT=ierr )
      IF( ierr/=0 ) CALL errore(' runcg ', ' deallocating hacca ',ierr)
      DEALLOCATE( dt2bye, STAT=ierr )
      IF( ierr/=0 ) CALL errore(' runcg ', ' deallocating dt2bye ',ierr)

      RETURN
      END SUBROUTINE runcg_new


! ---------------------------------------------------------------------- !
!
!  The following subroutine performs the line minimizations along "hacca"
!
! ---------------------------------------------------------------------- !

    SUBROUTINE cglinmin(emin, ediff, tbad, edft, cp, c, cdesc, occ, vpot, rhoe, hacca, &
        atoms, ht, bec, becdr, eigr, ei1, ei2, ei3, sfac)

! ... declare modules

        USE wave_types
        USE energies, ONLY: dft_energy_type
        USE wave_functions, ONLY: fixwave
        USE io_global, ONLY: ionode
        USE io_global, ONLY: stdout
        USE cell_module, ONLY: boxdimensions
        USE potentials, ONLY: kspotential
        USE atoms_type_module, ONLY: atoms_type
        USE reciprocal_space_mesh, ONLY: gkmask_l
        USE uspp,             ONLY : vkb, nkb

        IMPLICIT NONE

! ...   ARGUMENTS
        REAL(DP) :: ediff, emin
        LOGICAL :: tbad
        TYPE (atoms_type), INTENT(INOUT) :: atoms
        COMPLEX(DP), INTENT(IN) :: c(:,:,:,:)
        COMPLEX(DP), INTENT(INOUT) :: cp(:,:,:,:)
        TYPE (wave_descriptor), INTENT(IN) :: cdesc
        REAL(DP) :: rhoe(:,:)
        COMPLEX(DP) :: sfac(:,:)
        COMPLEX(DP) :: eigr(:,:)
        COMPLEX(DP) :: ei1(:,:)
        COMPLEX(DP) :: ei2(:,:)
        COMPLEX(DP) :: ei3(:,:)
        TYPE (boxdimensions), INTENT(INOUT) ::  ht
        REAL(DP), INTENT(IN) :: occ(:,:,:)
        REAL(DP) :: bec(:,:)
        REAL(DP) :: becdr(:,:,:)
        TYPE (dft_energy_type) :: edft
        COMPLEX (DP) ::  hacca(:,:,:,:)
        REAL (DP), INTENT(in) ::  vpot(:,:)

!
! ... LOCALS
!

        REAL(DP) :: GOLD, GLIMIT, TINY, CGOLD, ZEPS
        INTEGER   :: itmax
        PARAMETER (GOLD=1.618034D0, GLIMIT=100.D0, TINY=1.D-20)
        PARAMETER (ITMAX=20, CGOLD=.3819660D0,ZEPS=1.0D-10)

        REAL(DP) :: ax, bx, cx, fa, fb, fc, dum, u, fu ,r, q, ulim
        REAL(DP) :: x, p, v, w, e, fw, fv, xm, tol1, tol2, a, b, etemp, d
        REAL(DP) :: fx, xmin, brent, eold, tol
        LOGICAL   :: tbrent
        INTEGER   :: iter, ispin

!
! ... SUBROUTINE BODY
!
        tbrent         = .FALSE.
        tol = 1.0d-8
        ax = 0.0d0
        bx = 1.0d0
        tbad = .FALSE.


        ! FA=FUNC(AX)
        fa =  cgenergy( ax )

        eold = fa

        ! FB=FUNC(BX)
        fb =  cgenergy( bx )

        IF(FB.GT.FA)THEN
          tbad = .TRUE.
          DUM=AX; AX=BX; BX=DUM
          DUM=FB; FB=FA; FA=DUM
        ENDIF
        CX=BX+GOLD*(BX-AX)

        ! FC=FUNC(CX)
        fc =  cgenergy( cx )

100     IF(FB.GE.FC)THEN
          R=(BX-AX)*(FB-FC)
          Q=(BX-CX)*(FB-FA)
          U=BX-((BX-CX)*Q-(BX-AX)*R)/(2.*SIGN(MAX(ABS(Q-R),TINY),Q-R))
          ULIM=BX+GLIMIT*(CX-BX)
          IF((BX-U)*(U-CX).GT.0.)THEN
            ! FU=FUNC(U)
            fu =  cgenergy( u )
            IF(FU.LT.FC)THEN
              AX=BX; FA=FB; BX=U; FB=FU;
              GO TO 100
            ELSE IF(FU.GT.FB)THEN
              CX=U; FC=FU;
              GO TO 100
            ENDIF
            U=CX+GOLD*(CX-BX)
            ! FU=FUNC(U)
            fu =  cgenergy( u )
          ELSE IF((CX-U)*(U-ULIM).GT.0.)THEN
            ! FU=FUNC(U)
            fu =  cgenergy( u )
            IF(FU.LT.FC)THEN
              BX=CX; CX=U
              U=CX+GOLD*(CX-BX)
              FB=FC; FC=FU
              ! FU=FUNC(U)
              fu =  cgenergy( u )
            ENDIF
          ELSE IF((U-ULIM)*(ULIM-CX).GE.0.)THEN
            U=ULIM
            ! FU=FUNC(U)
            fu =  cgenergy( u )
          ELSE
            U=CX+GOLD*(CX-BX)
            ! FU=FUNC(U)
            fu =  cgenergy( u )
          ENDIF
          AX=BX; BX=CX; CX=U; FA=FB; FB=FC; FC=FU
          GO TO 100
        ENDIF

        IF( tbrent .AND. tbad ) THEN

          IF( ionode .AND. cg_prn ) WRITE( stdout,114) ax, bx, cx, fa, fb, fc

          A=MIN(AX,CX); B=MAX(AX,CX)
          V=BX; W=V; X=V; E=0.d0
          ! FX=F(X)
          fx =  cgenergy( x )
          FV=FX; FW=FX

          DO ITER = 1, ITMAX
            XM = 0.5d0 * (A+B)
            ! TOL1=TOL*ABS(X)+ZEPS
            TOL1 = TOL * ABS(X) + ZEPS
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
            fu =  cgenergy( u )
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
          ! CALL errore('CGLINMIN', 'Brent exceed maximum iterations.',itmax)
          WRITE( stdout, fmt='(" CGLINMIN, WARNING: Brent exceed maximum iterations ")' )
103       XMIN=X
          BRENT=FX
  
        ELSE

          x = bx

        END IF

        emin =  cgenergy( x )

        IF( ionode .AND. cg_prn ) WRITE( stdout,114) ax, x, cx, fa, emin, fc

        IF( tbad ) THEN
          ediff = ABS(emin - fa)
        ELSE
          ediff = ABS(emin - eold)
        END IF

113     FORMAT(6X,'lm',I5,2X,3F22.18,2X,2F10.6)
114     FORMAT(6X,'lm',3F10.5,3F12.6)

    CONTAINS

      REAL(DP) FUNCTION cgenergy( hstep )
 
        ! ...   ARGUMENTS

        REAL(DP) :: hstep

        ! ... LOCALS

        LOGICAL      ttprint, ttforce, ttstress, tcel
        REAL(DP) :: timepre

        ! ...      SUBROUTINE BODY

        ttprint = .FALSE.
        ttforce = .FALSE.
        tcel     = .FALSE.
        ttstress = .FALSE.

        cp = c + hstep * hacca

        CALL fixwave( cp, cdesc, gkmask_l )

        DO ispin = 1, cdesc%nspin
           CALL gram( vkb, bec, nkb, cp(1,1,1,ispin), SIZE(cp,1), cdesc%nbt( ispin ) )
        END DO


        CALL kspotential( 1, ttprint, ttforce, ttstress, rhoe, &
            atoms, bec, becdr, eigr, ei1, ei2, ei3, sfac, cp, cdesc, tcel, ht, occ, vpot, edft, timepre )

        cgenergy = edft%etot

      END FUNCTION cgenergy

    END SUBROUTINE cglinmin

! ---------------------------------------------------------------------- !
      END MODULE runcg_module
! ---------------------------------------------------------------------- !
