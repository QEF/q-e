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
  USE constants,    ONLY : fpi, sqrt2, eps8
  USE atom,         ONLY : rgrid
  USE ions_base,    ONLY : ntyp => nsp
  USE us,           ONLY : dq
  USE uspp_param,   ONLY : upf, nbetam
  USE mp_bands,     ONLY : intra_bgrp_comm
  USE mp,           ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  LOGICAL, PARAMETER :: tprint=.FALSE.   ! Whether the beta_l(r) and 
                                         ! its relatives are printed or not.
  INTEGER, PARAMETER :: nn=16    ! Smoothing parameter, order of the polynomial
                                 ! inverse gaussian approximant.
  REAL(DP), PARAMETER :: a=22.0  ! Smoothing parameter, exponent of the gaussian
                                 ! decaying factor.
  REAL(DP) :: rcut, drcut ! beta function cutoff radius and estimated increase of
                          ! it due to the filtering
  REAL(DP), PARAMETER :: eps = 1.d-8
  !
  INTEGER :: nqx
  REAL(DP), ALLOCATABLE :: tab0(:,:), tab(:,:), beta(:,:), betas(:,:)
  REAL(DP), ALLOCATABLE :: power_r(:), power_q(:)
  !
  INTEGER :: nt, nb, l, ir, iq, startq, lastq, ndm
  ! various counters
  REAL(DP), ALLOCATABLE :: aux(:), besr(:)
  ! various work space
  REAL(DP) :: q, qi, qmax
  ! the modulus of g for each shell
  ! q-point grid for interpolation
  REAL(DP) ::  vqint 
  ! interpolated value
  !
  CHARACTER(LEN=4) :: filename
  !
  CALL start_clock( 'init_us_b0' )
  !
  ! ... Initialization of variables
  !
  drcut = ABS(LOG(eps8))/SQRT(ecutwfc)/2.d0 
  qmax = 3.d0 * SQRT(ecutwfc)
  nqx = INT( qmax / dq + 4)  ! Think about what happens in a variable cell calculations
  !
  ALLOCATE( tab0(nqx,nbetam), tab(nqx,nbetam) )
  ALLOCATE( power_r(nbetam), power_q(nbetam) )
  !
  IF (tprint) WRITE( stdout, * ) 'upf(nt)%kkbeta ', upf(1:ntyp)%kkbeta
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
  ALLOCATE( aux(ndm),besr(ndm) )
  !
  WRITE(6,'(a,f6.2,a,i4,4(a,f11.8))') 'FILTER : a=',a,', nn=',nn,', filter(1.1d0)=',   filter(1.1d0,a,nn),  &
                                                                 ', filter(1.1d0/3)=', filter(1.1d0/3,a,nn),&
                                                                 ', filter(1.2d0)=',   filter(1.2d0,a,nn),  &
                                                                 ', filter(1.2d0/3)=', filter(1.2d0/3,a,nn)
  IF (tprint) THEN
     WRITE( stdout, * ) " PSEUDOPOTENTIAL REPORT; initial values "
     WRITE( stdout, * ) ' NQX :', nqx, ' NBETAM :', nbetam
     WRITE( stdout, * ) ' NDM :', ndm, ' KKBETA :', upf(1:ntyp)%kkbeta
     WRITE( stdout, * ) ' r(ndm):', (rgrid(nt)%r(upf(nt)%kkbeta),nt=1,ntyp), ' drcut :', drcut
  ENDIF
  !
  ! ... fill the interpolation table tab
  !
  CALL divide( intra_bgrp_comm, nqx, startq, lastq )
  DO nt = 1, ntyp
     !     WRITE(*,*) ' ityp = ', nt
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
           ! tab (iq, nb) = vqint * filter( 1.1d0 * qi / qmax, a, nn )
           tab (iq, nb) = vqint * filter( 1.2d0 * qi / qmax, a, nn )
           !
           DO ir = 1, upf(nt)%kkbeta
              beta (ir,nb) = beta (ir,nb) + qi*qi * dq * tab0(iq,nb) * besr (ir) * rgrid(nt)%r(ir)
              betas(ir,nb) = betas(ir,nb) + qi*qi * dq * tab (iq,nb) * besr (ir) * rgrid(nt)%r(ir)
           ENDDO
        ENDDO
        !-
     ENDDO
     !
     CALL mp_sum( tab0, intra_bgrp_comm )
     CALL mp_sum( tab , intra_bgrp_comm )
     beta = beta * 8.0/fpi ; CALL mp_sum( beta , intra_bgrp_comm )
     betas= betas* 8.0/fpi ; CALL mp_sum( betas, intra_bgrp_comm )
     upf(nt)%beta(1:upf(nt)%kkbeta,1:upf(nt)%nbeta) = betas(1:upf(nt)%kkbeta,1:upf(nt)%nbeta) 
     !
     !-Fix the core cutoff radius kkbeta so that no more than eps of the integral is lost.
     ! Parseval's identity is used to monitor the real space completeness of the integral of the(squarred)  betas
     !
     ! Compute the integral in Fourier space
     DO nb = 1, upf(nt)%nbeta
        l = upf(nt)%lll (nb)
        power_q(nb) =0.d0
        DO iq = startq, lastq
           qi = (iq - 1) * dq
           power_q(nb) = power_q(nb) + qi*qi * dq * tab(iq,nb) * tab(iq,nb)
        ENDDO
     ENDDO
     power_q  = power_q * 8.0/fpi ; CALL mp_sum( power_q, intra_bgrp_comm )
     !
     ! Compute the integral in real space up to the current value of kkbeta and try to reduce it if close enough to power_q
     power_r = power_q                    ! initialize power_r so that the first test is passed. 
     upf(nt)%kkbeta = upf(nt)%kkbeta + 1  ! increase kkbeta because it must be decreased in a moment
     DO WHILE ( MAXVAL(1.d0-power_r(1:upf(nt)%nbeta)/power_q(1:upf(nt)%nbeta)) < eps ) 
        WRITE(*,*) upf(nt)%kkbeta,rgrid(nt)%r(upf(nt)%kkbeta),1.d0-power_r(1:upf(nt)%nbeta)/power_q(1:upf(nt)%nbeta)
        upf(nt)%kkbeta = upf(nt)%kkbeta - 1
        DO nb = 1, upf(nt)%nbeta
           l = upf(nt)%lll (nb)
           DO ir = 1, upf(nt)%kkbeta
              aux(ir) = betas(ir,nb) * betas(ir,nb)
           ENDDO
           CALL simpson( upf(nt)%kkbeta, aux, rgrid(nt)%rab, vqint )
           power_r(nb) = vqint
        ENDDO
     ENDDO
     !
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
  IF (tprint) THEN
     filename(1:4) = 'ffqq'                ! the filter in reciprocal space
     OPEN (4,FILE=filename, FORM='formatted', STATUS='unknown')
     WRITE(4, *) '# the filter in reciprocal space'
     WRITE(4, *) '# fillter : a=',a,', nn=', nn
     DO iq = 1, nqx
        qi = (iq - 1) * dq
        WRITE(4,'(12f16.10)')  qi, filter(1.1d0*qi/qmax, a, nn)
     ENDDO
     CLOSE (4)
  ENDIF
  IF (tprint) THEN
     WRITE( stdout, * ) " PSEUDOPOTENTIAL REPORT; updated values "
     WRITE( stdout, * ) ' NDM :', ndm, ' KKBETA :', upf(1:ntyp)%kkbeta
     WRITE( stdout, * ) ' r(ndm):', (rgrid(nt)%r(upf(nt)%kkbeta),nt=1,ntyp)
  ENDIF
  !
  DEALLOCATE( power_r, power_q )
  DEALLOCATE( tab0, tab )
  DEALLOCATE( beta, betas )
  DEALLOCATE( besr, aux )
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

