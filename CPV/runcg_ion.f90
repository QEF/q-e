!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! ---------------------------------------------------------------------- !
      MODULE runcg_ion_module
! ---------------------------------------------------------------------- !

        USE kinds, ONLY: DP

        IMPLICIT NONE
        SAVE
 
        PRIVATE

        REAL(DP) :: old_clock_value = 0.0d0

        LOGICAL :: cg_prn = .FALSE. 

        PUBLIC :: runcg_ion

! ---------------------------------------------------------------------- !
      CONTAINS
! ---------------------------------------------------------------------- !

!  -----------------------------------------------------------------------
!  BEGIN manual

   SUBROUTINE runcg_ion(nfi, tortho, tprint, rhoe, atomsp, atoms0, atomsm, &
      bec, becdr, eigr, vkb, ei1, ei2, ei3, sfac, c0, cm, cp, cdesc, tcel, ht, fi, ei, &
      vpot, doions, edft, etol, ftol, maxiter, sdthr, maxnstep )

!  this routine computes the equilibrium ionic positions via conjugate gradient
!  END manual

      ! ... declare modules

      USE energies, ONLY: dft_energy_type, print_energies
      USE cp_interfaces, ONLY: update_wave_functions, printout
      USE wave_base, ONLY: dotp
      USE check_stop, ONLY: check_stop_now
      USE io_global, ONLY: ionode
      USE io_global, ONLY: stdout
      USE cell_base, ONLY: boxdimensions, s_to_r, r_to_s
      USE wave_types, ONLY: wave_descriptor
      USE time_step, ONLY: delt
      USE atoms_type_module, ONLY: atoms_type
      USE parameters, ONLY: nacx
      USE runsd_module, ONLY: runsd

      IMPLICIT NONE

      ! ... declare subroutine arguments

      INTEGER   :: nfi
      LOGICAL   :: tortho, tprint, tcel, doions
      TYPE (atoms_type) :: atomsp
      TYPE (atoms_type) :: atoms0
      TYPE (atoms_type) :: atomsm
      COMPLEX(DP), INTENT(INOUT) :: c0(:,:), cm(:,:), cp(:,:)
      TYPE (wave_descriptor) :: cdesc
      REAL(DP) :: rhoe(:,:)
      REAL(DP) :: bec(:,:)
      REAL(DP) :: becdr(:,:,:)
      COMPLEX(DP) :: eigr(:,:)
      COMPLEX(DP) :: vkb(:,:)
      COMPLEX(DP) :: ei1(:,:)
      COMPLEX(DP) :: ei2(:,:)
      COMPLEX(DP) :: ei3(:,:)
      COMPLEX(DP) :: sfac(:,:)
      TYPE (boxdimensions), INTENT(INOUT) ::  ht
      REAL(DP)  :: fi(:)
      TYPE (dft_energy_type) :: edft

      REAL(DP)    :: ei(:,:)
      REAL(DP)    :: vpot(:,:)

      INTEGER, INTENT(IN) :: maxnstep, maxiter
      REAL(DP), INTENT(IN) :: sdthr, etol, ftol

      ! ... declare other variables

      LOGICAL :: ttsde, ttprint, ttforce, ttstress, ttortho
      LOGICAL :: tbad
      REAL(DP) :: s0, s1, s2, s3, s4, s5, s6, seconds_per_iter
      REAL(DP) :: dene, eold, timerd, timeorto, ekinc

      REAL(DP)    :: gg, ggo, dgg, emin, demin, gam, fp, fret
      REAL(DP) :: lambda, fions(3), dumm

      REAL(DP), POINTER :: hacca(:,:)
      REAL(DP), POINTER :: gnew(:,:)
      REAL(DP), POINTER :: xi(:,:)

      INTEGER :: i, iter, ierr, is, ia, isa, k, j

      REAL(DP) :: displ, amtot

      REAL(DP) :: eps = 1.0d-20

      REAL(DP) ::  avgs(nacx)
      REAL(DP) ::  avgs_this_run(nacx)
      INTEGER   ::  nat

      REAL(DP), EXTERNAL :: cclock


! ... end of declarations
!  ----------------------------------------------

      doions      = .FALSE.
      eold        = 1.0d10  ! a large number
      timerd      = 0
      timeorto    = 0
      ttsde       = .TRUE.
      ttprint     = .FALSE.
      ttstress    = .FALSE.
      tcel    = .FALSE.
      ttortho = .TRUE.
      ttforce = .TRUE.
      tbad   = .FALSE.
      avgs   = 0.0d0
      avgs_this_run   =  0.0d0
      nat = atoms0%nat

!      maxnstep = 300
!      maxiter = 100

      amtot = 0.0d0
      DO is = 1, atoms0%nsp
        amtot = amtot + atoms0%m(is) * atoms0%na(is)
      END DO
      amtot = amtot / nat
      displ = delt * delt / amtot


      ALLOCATE(hacca( 3, nat ), gnew( 3, nat), xi( 3, nat), STAT=ierr )
      IF( ierr/=0 ) CALL errore(' runcg_ion ', ' allocating hacca ',ierr)

      IF(ionode) THEN
        WRITE( stdout,'(/,8X,"Conjugate Gradient Optimizations for Inos, starting ...")' )
        WRITE( stdout, fmt='(8X,"Displ = ",F10.6)' ) displ
        WRITE( stdout, fmt='(8X,"IonThr = ",D14.6," NstepIx = ",I5)' ) etol, maxiter
        WRITE( stdout, fmt='(8X,"ForThr = ",D14.6," NstepIx = ",I5)' ) ftol, maxiter
        WRITE( stdout, fmt='(8X,"EleThr = ",D14.6," NstepEx = ",I5)' ) sdthr, maxnstep
        WRITE( stdout, * )
      END IF

      s1 = cclock()
      old_clock_value = s1

      CALL runsd(ttortho, ttprint, ttforce, rhoe, atoms0, bec, becdr, eigr, vkb, &
                 ei1, ei2, ei3, sfac, c0, cm, cp, cdesc, tcel, ht, fi, ei, vpot, &
                 doions, edft, maxnstep, sdthr )

      IF( ionode .AND. cg_prn ) THEN
        DO j = 1, atoms0%nat
          WRITE( stdout,fmt="(6X,'F ',3D14.6)") (atoms0%for(i,j),i=1,3)
        END DO
      END IF

      xi(1:3,1:nat) = - atoms0%for(1:3,1:nat) 
      gnew  = -xi
      hacca = gnew
      xi    = hacca
      fp    = edft%etot 

      CONJUGATE_GRADIENTS: DO iter = nfi, nfi+maxiter

        s2 = cclock()

! ...   check for exit
        IF (check_stop_now()) THEN
          EXIT CONJUGATE_GRADIENTS
        END IF

        IF(ionode) &
          WRITE( stdout,fmt="(/,8X,'cgion: iter',I5,' line minimization along gradient starting')") iter

        CALL cglinmin(fret, edft, cp, c0, cm, cdesc, fi, ei, vpot, rhoe, xi, atomsp, atoms0, &
          ht, bec, becdr, eigr, vkb, ei1, ei2, ei3, sfac, maxnstep, sdthr, displ)

        IF( tbad ) THEN
!          displ = displ * 2.0d0
!          tbad = .FALSE.
        END IF

        s3 = cclock()

        IF( fp <= fret ) THEN

          IF( ionode ) WRITE( stdout, fmt='(8X,"cgion: bad step")')  ! perform steepest descent
          displ = displ / 2.0d0

          CALL runsd(ttortho, ttprint, ttforce, rhoe, atoms0, bec, becdr, eigr, vkb, ei1, ei2, ei3, &
                     sfac, c0, cm, cp, cdesc, tcel, ht, fi, ei, vpot, doions, edft, maxnstep, sdthr )
        
!          tbad = .TRUE.

          CYCLE CONJUGATE_GRADIENTS

        ELSE

          atoms0%taus(:,:) = atomsp%taus(:,:)
          atoms0%for(:,:) = atomsp%for(:,:)

          IF( ( 2.0d0 * ABS( fret - fp )                         < etol ) .AND. &
              ( MAXVAL( ABS( atoms0%for( 1:3, 1:atoms0%nat ) ) ) < ftol ) ) THEN
            IF(ionode) WRITE( stdout,fmt="(8X,'cgion:  convergence achieved successfully',/)")
            doions = .TRUE.
            EXIT CONJUGATE_GRADIENTS
          END IF

        END IF

        s0 = cclock()
        seconds_per_iter = (s0 - old_clock_value)
        old_clock_value = s0
        IF( ionode ) THEN
          WRITE( stdout,'(/,8X,"cgion:   iter     erho            derho      seconds")' )
          WRITE( stdout,113) iter, fret, ABS( fret-fp ), seconds_per_iter
113       FORMAT(8X,'cgion:',I5,2X,F14.6,2X,2D12.4)
        END IF

        CALL printout(iter, atoms0, 0.0d0, 0.0d0, .TRUE., ht, edft)

        fp  = fret
        xi(1:3,1:nat) = - atoms0%for(1:3,1:nat) 
        gg  = 0.0d0
        dgg = 0.0d0

        IF( ionode .AND. cg_prn ) THEN
          DO j = 1, SIZE(hacca, 2)
            WRITE( stdout,fmt="(6X,'F ',3D14.6)") (atoms0%for(i,j),i=1,3)
          END DO
        END IF

        DO isa = 1, nat
          DO k = 1, 3
            gg  = gg + gnew(k,isa)**2 
            dgg = dgg + xi(k,isa)**2
            dgg = dgg + ( xi(k,isa) + gnew(k,isa) ) * xi(k,isa)
          END DO
        END DO

        IF( gg == 0.0d0 ) THEN
          IF(ionode) WRITE( stdout,fmt="(8X,'cgion:  convergence achieved successfully',/)")
          doions = .TRUE.
          EXIT CONJUGATE_GRADIENTS
        END IF

        gam = dgg / gg

        gnew = -xi
        hacca = gnew + gam * hacca
        xi = hacca

      END DO CONJUGATE_GRADIENTS

! ... set wave functions velocity to 0
      atomsm%taus(:,:) = atoms0%taus(:,:)
      atomsm%for(:,:) = atoms0%for(:,:)

      IF( (iter .GT. maxiter) .AND. ionode) THEN
        WRITE( stdout,fmt="(8X,'cgion:  convergence not achieved')")
        WRITE( stdout,fmt="(8X,'cgion:  maximum number of iteration exceeded',/)")
      END IF

      CALL printout(iter, atoms0, 0.0d0, 0.0d0, .TRUE., ht, edft)

      DEALLOCATE( hacca, gnew, xi, STAT=ierr )
      IF( ierr/=0 ) CALL errore(' runcg_ion ', ' deallocating hacca ',ierr)

      nfi = iter

      RETURN
      END SUBROUTINE runcg_ion


! ---------------------------------------------------------------------- !
! ---------------------------------------------------------------------- !

      SUBROUTINE cglinmin(emin, edft, cp, c0, cm, cdesc, fi, ei, vpot, &
        rhoe, hacca, atomsp, atoms0, ht, bec, becdr, eigr, vkb, ei1, ei2, ei3, sfac, &
        maxnstep, sdthr, displ)

! ... declare modules

        USE wave_types, ONLY: wave_descriptor
        USE energies, ONLY: dft_energy_type
        USE cp_interfaces, ONLY: update_wave_functions
        USE io_global, ONLY: ionode
        USE io_global, ONLY: stdout
        USE cell_base, ONLY: boxdimensions, r_to_s
        USE atoms_type_module, ONLY: atoms_type
        USE check_stop, ONLY: check_stop_now
        USE runsd_module, ONLY: runsd

        IMPLICIT NONE

! ...   ARGUMENTS
        REAL(DP) :: emin
        TYPE (atoms_type) :: atomsp
        TYPE (atoms_type) :: atoms0
        COMPLEX(DP), INTENT(INOUT) :: c0(:,:)
        COMPLEX(DP), INTENT(INOUT) :: cp(:,:)
        COMPLEX(DP), INTENT(INOUT) :: cm(:,:)
        TYPE (wave_descriptor) :: cdesc
        REAL(DP) :: rhoe(:,:)
        COMPLEX(DP) :: eigr(:,:)
        COMPLEX(DP) :: vkb(:,:)
        COMPLEX(DP) :: ei1(:,:)
        COMPLEX(DP) :: ei2(:,:)
        COMPLEX(DP) :: ei3(:,:)
        COMPLEX(DP) :: sfac(:,:)
        TYPE (boxdimensions), INTENT(INOUT) ::  ht
        REAL(DP)  :: fi(:)
        TYPE (dft_energy_type) :: edft
        REAL (DP) ::  hacca(:,:)
        REAL (DP), INTENT(in) ::  vpot(:,:)
        REAL(DP) :: bec(:,:)
        REAL(DP) :: becdr(:,:,:)

        REAL(DP)    :: ei(:,:)

        INTEGER :: maxnstep
        REAL(DP) :: sdthr
        REAL(DP) :: displ

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
        INTEGER   :: iter, i, j

!
! ... SUBROUTINE BODY
!
        tbrent         = .FALSE.
        tol = 1.0d-8
        ax = 0.0d0
        bx = displ

        IF( ionode .AND. cg_prn ) THEN
          DO j = 1, SIZE(hacca, 2)
            WRITE( stdout,120) (hacca(i,j),i=1,3)
          END DO
 120      FORMAT(6X,'H ',3D14.6)
        END IF

        fa =  cgenergy( ax )

        eold = fa

        fb =  cgenergy( bx )

        IF(FB.GT.FA)THEN
          DUM=AX; AX=BX; BX=DUM
          DUM=FB; FB=FA; FA=DUM
        ENDIF
        CX=BX+GOLD*(BX-AX)

        fc =  cgenergy( cx )

100     IF(FB.GE.FC)THEN

          IF (check_stop_now()) THEN
            GO TO 300 
          END IF

          R=(BX-AX)*(FB-FC)
          Q=(BX-CX)*(FB-FA)
          U=BX-((BX-CX)*Q-(BX-AX)*R)/(2.d0*SIGN(MAX(ABS(Q-R),TINY),Q-R))
          ULIM=BX+GLIMIT*(CX-BX)
          IF((BX-U)*(U-CX).GT.0.)THEN
            fu =  cgenergy( u )
            IF(FU.LT.FC)THEN
              AX=BX; FA=FB; BX=U; FB=FU;
              GO TO 100
            ELSE IF(FU.GT.FB)THEN
              CX=U; FC=FU;
              GO TO 100
            ENDIF
            U=CX+GOLD*(CX-BX)
            fu =  cgenergy( u )
          ELSE IF((CX-U)*(U-ULIM).GT.0.)THEN
            fu =  cgenergy( u )
            IF(FU.LT.FC)THEN
              BX=CX; CX=U
              U=CX+GOLD*(CX-BX)
              FB=FC; FC=FU
              fu =  cgenergy( u )
            ENDIF
          ELSE IF((U-ULIM)*(ULIM-CX).GE.0.)THEN
            U=ULIM
            fu =  cgenergy( u )
          ELSE
            U=CX+GOLD*(CX-BX)
            fu =  cgenergy( u )
          ENDIF
          AX=BX; BX=CX; CX=U; FA=FB; FB=FC; FC=FU
          GO TO 100
        ENDIF

        IF( tbrent ) THEN

          IF( ionode .AND. cg_prn ) WRITE( stdout,114) ax, bx, cx, fa, fb, fc

          A=MIN(AX,CX); B=MAX(AX,CX)
          V=BX; W=V; X=V; E=0.d0
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
          WRITE( stdout, fmt='(" CGLINMIN, WARNING: Brent exceed maximum iterations ")' )
103       XMIN=X
          BRENT=FX
  
        ELSE

          x = bx

        END IF

300     continue
        IF (check_stop_now()) THEN
          x = 0.0d0
        END IF

        emin =  cgenergy( x )

        IF( ionode .AND. cg_prn ) WRITE( stdout,114) ax, x, cx, fa, emin, fc

113     FORMAT(6X,'lm',I5,2X,3F22.18,2X,2F10.6)
114     FORMAT(6X,'lm',3F10.5,3F12.6)

    CONTAINS

       REAL(DP) FUNCTION cgenergy( hstep )

         REAL(DP) :: hstep

         ! ... LOCALS

           INTEGER      ia, is, isa, k
           LOGICAL      ttprint, ttforce, ttstress, tcel, ttortho, doions
           REAL(DP) :: fions(3), dumm

         ! ... SUBROUTINE BODY

           tcel    = .FALSE.
           ttortho = .TRUE.
           ttprint = .FALSE.
           ttforce = .TRUE.
           doions  = .FALSE.

         ! ...  HERE UPDATE THE IONIC POSITION ALONG THE LINE

           isa = 0
           DO is = 1, atoms0%nsp
             DO ia = 1, atoms0%na(is)
               isa = isa + 1
               CALL r_to_s(hacca(:,isa), fions, ht)
               DO k = 1, 3
                 IF(  atoms0%mobile(k,isa) > 0 )THEN
                   dumm = hstep * fions(k)
                   atomsp%taus(k,isa) = atoms0%taus(k,isa) + dumm
                 ELSE
                   atomsp%taus(k,isa) = atoms0%taus(k,isa)
                 END IF
               END DO
             END DO
           END DO

         ! ...  Calculate Forces (fion) and DFT Total Energy (edft) for the new ionic
         ! ...  positions (atomsp)

           CALL runsd(ttortho, ttprint, ttforce, rhoe, atomsp, bec, becdr, eigr, &
                      vkb, ei1, ei2, ei3, sfac, c0, cm, cp, cdesc, tcel, ht, fi, ei, &
                      vpot, doions, edft, maxnstep, sdthr )

           cgenergy = edft%etot

         END FUNCTION

      END SUBROUTINE cglinmin

! ---------------------------------------------------------------------- !
    END MODULE runcg_ion_module
! ---------------------------------------------------------------------- !
