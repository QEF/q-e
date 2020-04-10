!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
SUBROUTINE init_us_b0
  !----------------------------------------------------------------------
  !! In this routine the beta_l(r) are smoothed.
  !
  USE kinds,        ONLY : DP
  USE gvecw,        ONLY : ecutwfc
  USE io_global,    ONLY : stdout
  USE constants,    ONLY : fpi
  USE atom,         ONLY : rgrid
  USE ions_base,    ONLY : ntyp => nsp
  USE us,           ONLY : dq
  USE uspp_param,   ONLY : upf, nbetam
  USE mp_bands,     ONLY : intra_bgrp_comm
  USE mp,           ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  LOGICAL, PARAMETER :: tprint=.FALSE.     ! Whether the beta_l(r) and its relatives are printed or not.
  INTEGER, PARAMETER :: nn=20              ! Smoothing parameter, order of the polynomial in the
                                           ! inverse gaussian approximant.
  REAL(DP), PARAMETER:: a=1.125d0*(2*nn+1) ! Smoothing parameter, exponent of the gaussian
                                           ! decaying factor, chosen so that the flex is at x=2/3.
  REAL(DP) :: rcut, drcut ! beta function cutoff radius and its estimated increase due to the filtering

  REAL(DP), PARAMETER :: eps = 1.d-9 ! error tollerance for intergrals, norms etc.
  !
  INTEGER :: nqx
  REAL(DP), ALLOCATABLE :: tab0(:,:), tab(:,:), beta(:,:), betas(:,:)
  REAL(DP), ALLOCATABLE :: power_r(:), power_q(:)
  !
  INTEGER :: nt, nb, l, ir, iq, startq, lastq, ndm
  ! various counters
  REAL(DP), ALLOCATABLE :: aux(:), besr(:), ffrr(:)
  ! various work space
  REAL(DP) :: qi, qmax, vqint
  ! the modulus of g for each shell
  ! q-point grid for interpolation
  ! interpolated value
  INTEGER :: kktest
  REAL(DP) :: test, safe_test, error_estimate
  LOGICAL :: first
  !
  CHARACTER(LEN=4) :: filename
  !
  CALL start_clock( 'init_us_b0' )
  !
  ! ... Initialization of variables
  !
  drcut = 3.D0 * ABS(LOG(eps)) / SQRT(ecutwfc)/ 2.D0  ! estimate (to be restricted later)
  qmax = 3.D0 * SQRT(ecutwfc) ! maximum q that avoids aliasing, higher Fourier components must be negligible
  nqx = INT( qmax / dq + 4)   ! Think about what happens in a variable cell calculations
  !
  ALLOCATE( tab0(nqx,nbetam), tab(nqx,nbetam) )
  ALLOCATE( power_r(nbetam), power_q(nbetam) )
  !
  WRITE( stdout, '(/5X,"PSEUDOPOTENTIAL Beta Smoothing: ==============================================")' )
  IF (tprint) THEN
     WRITE( stdout,'(/5X,"table grid dimensions:  NQX =",I7," NBETAM =", I4)'), nqx, nbetam
     WRITE( stdout,'(11X,"initial rcut indices:",10i7)') upf(1:ntyp)%kkbeta
     WRITE( stdout,'(19X,"rcut (ityp) :",10f7.3)'), (rgrid(nt)%r(upf(nt)%kkbeta),nt=1,ntyp)
  END IF
  !
  DO nt=1,ntyp
     rcut = rgrid(nt)%r(upf(nt)%kkbeta)
     DO ir = upf(nt)%kkbeta, upf(nt)%mesh
        IF ( rgrid(nt)%r(ir) < rcut + drcut ) upf(nt)%kkbeta=ir
     ENDDO
  ENDDO
  !
  ndm = MAXVAL( upf(:)%kkbeta )
  ALLOCATE( beta(ndm,nbetam), betas(ndm,nbetam) )
  ALLOCATE( aux(ndm), besr(ndm), ffrr(ndm) )
  !
  WRITE(6,'(/5X,a)') 'Pseudopotential projectors are smoothed with filter( q / qmax, a, nn )'
  WRITE(6,'(5X,a,f6.2,a,i4,2(a,f11.8))') 'FILTER : a=',a,', nn=',nn,', filter(1.0/3)=', filter(1.d0/3,a,nn), &
                                                                    ', filter(1.0)=', filter(1.d0,a,nn)
  WRITE(6,'(14X,a,1p,e12.4)') 'target error: eps =', eps
  !
  ! ... fill the interpolation table tab
  !
  CALL divide( intra_bgrp_comm, nqx, startq, lastq )
  !
  !- loop over pseudopotentials
  DO nt = 1, ntyp
     WRITE( stdout, '(/5X,a,i4)' ) 'Smoothing PSEUDO #', nt
     IF ( upf(nt)%is_gth ) CYCLE
     !
     tab0(:,:) = 0.d0  ; tab(:,:) = 0.d0
     beta(:,:) = 0.d0  ; betas(:,:) = 0.d0
     !
     IF (tprint) THEN
       filename = 'br_'                 !  the radial beta_l(r) as defined in the pseudopotential
       WRITE(filename( 4:4 ),'(i1)') nt
       OPEN(4, FILE=filename, FORM='formatted', STATUS='unknown')
       WRITE(4,*) '# the radial beta_l(r) as defined in the pseudopotential'
       WRITE(4,*) '# nbeta :', upf(nt)%nbeta,' kkbeta :',upf(nt)%kkbeta
       DO ir = 1, upf(nt)%kkbeta
          WRITE(4,'(12f16.10)') rgrid(nt)%r(ir), (upf(nt)%beta(ir,nb), nb=1,upf(nt)%nbeta)
       ENDDO
       CLOSE (4)
     ENDIF
     !- compute original beta normalization in real space where grid is sufficient by definition
     DO nb = 1, upf(nt)%nbeta
        l = upf(nt)%lll (nb)
        kktest = upf(nt)%kkbeta
        aux(1:kktest) = upf(nt)%beta(1:kktest,nb) * upf(nt)%beta(1:kktest,nb)
        CALL simpson( kktest, aux, rgrid(nt)%rab, power_r(nb) )
     ENDDO
     IF (tprint) THEN
       WRITE( stdout, '(5X,a)' ) "BETA FUNCTIONS NORMS:"
       WRITE( stdout, '(5X,a,1p,10e12.4)'), "norm2 R:", (power_r(nb), nb=1,upf(nt)%nbeta)
     END IF
     !-
     DO nb = 1, upf(nt)%nbeta
        l = upf(nt)%lll (nb)
        DO iq = startq, lastq
           qi = (iq - 1) * dq
           CALL sph_bes (upf(nt)%kkbeta, rgrid(nt)%r, qi, l, besr)
           DO ir = 1, upf(nt)%kkbeta
              aux (ir) = upf(nt)%beta (ir, nb) * besr (ir) * rgrid(nt)%r(ir)
           ENDDO
           CALL simpson( upf(nt)%kkbeta, aux, rgrid(nt)%rab, vqint )
           !
           tab0(iq, nb) = vqint
           tab (iq, nb) = vqint * filter( qi / qmax, a, nn )
           !
           DO ir = 1, upf(nt)%kkbeta
              beta (ir,nb) = beta (ir,nb) + qi*qi * dq * tab0(iq,nb) * besr (ir) * rgrid(nt)%r(ir)
              betas(ir,nb) = betas(ir,nb) + qi*qi * dq * tab (iq,nb) * besr (ir) * rgrid(nt)%r(ir)
           ENDDO
        ENDDO
     ENDDO
     CALL mp_sum( tab0, intra_bgrp_comm )
     CALL mp_sum( tab , intra_bgrp_comm )
     beta = beta * 8.d0/fpi ; CALL mp_sum( beta , intra_bgrp_comm )
     betas= betas* 8.d0/fpi ; CALL mp_sum( betas, intra_bgrp_comm )
     upf(nt)%beta(1:upf(nt)%kkbeta,1:upf(nt)%nbeta) = betas(1:upf(nt)%kkbeta,1:upf(nt)%nbeta) 
     !-
     !- compute the beta normalization in q space to see how much of the original is lost
     DO nb = 1, upf(nt)%nbeta
        l = upf(nt)%lll (nb)
        power_q(nb) =0.d0
        DO iq = startq, lastq
           qi = (iq - 1) * dq
           power_q(nb) = power_q(nb) + qi*qi * dq * tab0(iq,nb) * tab0(iq,nb)
        ENDDO
     ENDDO
     power_q  = power_q * 8.d0/fpi ; CALL mp_sum( power_q, intra_bgrp_comm )
     !
     IF (tprint) WRITE( stdout, '(5X,a,1p,10e12.4)'), "norm2 Q:", (power_q(nb), nb=1,upf(nt)%nbeta)
     IF (tprint) WRITE( stdout, '(5X,a,1p,10e12.4)'), "1-ratio:", (1.D0-power_q(nb)/power_r(nb), nb=1,upf(nt)%nbeta)
     !
     error_estimate= MAXVAL(1.d0-power_q(1:upf(nt)%nbeta)/power_r(1:upf(nt)%nbeta)) * filter(1.d0,a,nn)**2
     IF (error_estimate > eps ) CALL errore( 'init_us_b0','R and Q norms of beta too different',nt)
     WRITE ( stdout, '(5x,"Smoothing truncation error estimate in Q",1pe12.4," < eps")' ) error_estimate
     !
     !-Fix the core cutoff radius kkbeta so that no more than eps of the integral is lost.
     ! Parseval's identity is used to monitor the real space completeness of the integral of betas**2
     !
     ! Compute the integral in Fourier space of the smoothed projectors
     DO nb = 1, upf(nt)%nbeta
        l = upf(nt)%lll (nb)
        power_q(nb) =0.d0
        DO iq = startq, lastq
           qi = (iq - 1) * dq
           power_q(nb) = power_q(nb) + qi*qi * dq * tab(iq,nb) * tab(iq,nb)
        ENDDO
     ENDDO
     power_q  = power_q * 8.d0/fpi ; CALL mp_sum( power_q, intra_bgrp_comm )
     !
     ! Compute norms in real space up to the current value of kkbeta and reduce it if close enough to power_q
     kktest = upf(nt)%kkbeta + 1 ;  test = 0.d0 ; first =.TRUE. ! initialize mesh limit and test value
     DO WHILE ( test < eps ) 
        safe_test = test ; kktest = kktest - 1
        DO nb = 1, upf(nt)%nbeta
           l = upf(nt)%lll (nb)
           aux(1:kktest) = betas(1:kktest,nb) * betas(1:kktest,nb)
           CALL simpson( kktest, aux, rgrid(nt)%rab, power_r(nb) )
        ENDDO
        IF ( first ) THEN
           first = .FALSE.
           test = MAXVAL( 1.d0 - power_r(1:upf(nt)%nbeta)/power_q(1:upf(nt)%nbeta) )
           if (test>eps) WRITE (stdout,'(5X,"WARNING: R and Q norms disagree by",6X,1pe12.4," > eps !")') test
           power_q = power_r
        END IF
        test = MAXVAL( 1.d0 - power_r(1:upf(nt)%nbeta)/power_q(1:upf(nt)%nbeta) )
!        write( stdout, * ) kktest, rgrid(nt)%r(kktest), test
     ENDDO
     upf(nt)%kkbeta = kktest + 1 ! set the mesh limit to the last value fulfilling the test condition
     WRITE( stdout, '(5X,"kkbeta, r(kkbeta), error in R",i5,f7.3,1pe11.4," < eps")' ) &
                            upf(nt)%kkbeta, rgrid(nt)%r(upf(nt)%kkbeta), safe_test
     !
     IF (tprint) THEN
        !-
        filename(1:3) = 'bq_'                 ! the radial fourier transform of beta_l in reciprcal space up to qmax
        OPEN( 4, FILE=filename, FORM='formatted', STATUS='unknown' )
        WRITE(4,*) '# the radial fourier transform of beta_l in reciprcal space up to qmax'
        WRITE(4,*) '# nbeta :', upf(nt)%nbeta,' nqx :',nqx, dq*nqx
        DO iq = 1, nqx
           qi = (iq - 1) * dq
           WRITE(4,'(12f16.10)')  qi,(tab0(iq,nb), nb=1,upf(nt)%nbeta)
        ENDDO
        CLOSE(4)
        !-
        filename(1:3) = 'bqs'                 ! the smoothed radial fourier transform of beta_l in reciprcal space up to qmax
        OPEN (4,FILE=filename,FORM='formatted', STATUS='unknown')
        WRITE(4,*) '# the smoothed radial fourier transform of beta_l in reciprcal space up to qmax'
        WRITE(4,*) '# nbeta :', upf(nt)%nbeta,' nqx :',nqx, dq*nqx
        do iq=1,nqx
           qi = (iq - 1) * dq
           WRITE(4,'(12f16.10)')  qi,(tab(iq,nb), nb=1,upf(nt)%nbeta)
        ENDDO
        CLOSE(4)
        !-
        filename(1:3) = 'brq'                 ! the back radial fourier transform of beta_l in real space
        OPEN (4,FILE=filename, FORM='formatted', STATUS='unknown')
        WRITE(4,*) '# the back radial fourier transform of beta_l in real space'
        WRITE(4,*) '# nbeta :', upf(nt)%nbeta,' kkbeta :',upf(nt)%kkbeta
        do ir =1,upf(nt)%kkbeta
           WRITE(4,'(12f16.10)') rgrid(nt)%r(ir), (beta(ir,nb), nb=1,upf(nt)%nbeta)
        ENDDO
        CLOSE(4)
        !-
        filename(1:3) = 'brs'                 ! the back radial fourier transform of the smoothed beta_l in real space
        OPEN(4,FILE=filename, FORM='formatted', STATUS='unknown')
        WRITE(4,*) '# the back radial fourier transform of the smoothed beta_l in real space'
        WRITE(4,*) '# nbeta :', upf(nt)%nbeta,' kkbeta :',upf(nt)%kkbeta
        do ir = 1, upf(nt)%kkbeta
           WRITE(4,'(12f16.10)') rgrid(nt)%r(ir), (betas(ir,nb), nb=1,upf(nt)%nbeta)
        ENDDO
        CLOSE(4)
        !-
     ENDIF
     !
  ENDDO
  !- end of loop over pseudopotentials
  IF (tprint) THEN
     WRITE( stdout,'(/13X,"final rcut indices:",10i7)') upf(1:ntyp)%kkbeta
     WRITE( stdout,'(19X,"rcut (ityp) :",10f7.3)'), (rgrid(nt)%r(upf(nt)%kkbeta),nt=1,ntyp)
  END IF
  WRITE( stdout, '(/5X,"PSEUDOPOTENTIAL end Beta Smoothing: ==========================================")' )
  !
  IF (tprint) THEN
     filename(1:4) = 'ffqq'                ! the filter in reciprocal space
     OPEN (4,FILE=filename, FORM='formatted', STATUS='unknown')
     WRITE(4, *) '# the filter in reciprocal space'
     WRITE(4, *) '# fillter : a=',a,', nn=', nn
     DO iq = 1, nqx
        qi = (iq - 1) * dq
        WRITE(4,'(12f16.10)')  qi, filter( qi/qmax, a, nn)
     ENDDO
     CLOSE (4)
  ENDIF
  IF (tprint) THEN
     DO nt = 1, ntyp
        filename = 'ffr'                   !  the filter in real space
        WRITE(filename( 4:4 ),'(i1)') nt
        !- (spherically) fourier transform the filter in real space
        ffrr(:) = 0.d0
        DO iq = startq, lastq
           qi = (iq - 1) * dq
           CALL sph_bes (upf(nt)%kkbeta, rgrid(nt)%r, qi, 0, besr)
           DO ir = 1, upf(nt)%kkbeta
              ffrr (ir) = ffrr (ir) + qi*qi * dq * filter( qi / qmax, a, nn) * besr (ir)
           ENDDO
        END DO
        CALL mp_sum( ffrr, intra_bgrp_comm )
        ffrr(:) = 8.D0 / fpi * ffrr(:)
        aux(1:upf(nt)%kkbeta) = ffrr(1:upf(nt)%kkbeta) * rgrid(nt)%r(1:upf(nt)%kkbeta)**2
        CALL simpson( upf(nt)%kkbeta, aux, rgrid(nt)%rab, vqint )
        OPEN (4,FILE=filename, FORM='formatted', STATUS='unknown')
        WRITE(4, *) '# the filter in real space'
        WRITE(4, *) '# fillter : a=',a,', nn=', nn
        DO ir = 1, upf(nt)%kkbeta
           WRITE(4,'(12f16.10)')  rgrid(nt)%r(ir), ffrr(ir)
        ENDDO
        CLOSE (4)
        !-
     END DO
  END IF
  !
  DEALLOCATE( power_r, power_q )
  DEALLOCATE( tab0, tab )
  DEALLOCATE( beta, betas )
  DEALLOCATE( besr, aux, ffrr )
  !
  CALL stop_clock( 'init_us_b0' )
  !
  RETURN
  !
  CONTAINS
  !
  REAL(DP) FUNCTION filter( x, a, n )
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN) :: x, a
    INTEGER,  INTENT(IN) :: n
    REAL(DP) :: axx, ff
    INTEGER :: k
    !
    axx = a * x * x
    !
    ff = 1.d0
    DO k = n, 1, -1
       ff = (1.d0+axx/REAL(k)*ff)
    ENDDO
    filter = ff * EXP(-axx)
    !
    RETURN
    !
  END FUNCTION filter
  !
END SUBROUTINE init_us_b0

