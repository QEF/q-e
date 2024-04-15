
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------------------
SUBROUTINE init_us_0(ecutrho,intra_bgrp_comm)
  !---------------------------------------------------------------------------------
  !! This routine performs the following task: for each uspp or paw pseudopotential
  !! the l-dependent aumentation charge \(\text{ q_nb_mb_l}\)(r), stored in
  !! \(\text{qfuncl}\)(ir,nmb,l), is:
  !
  !! * transformed in reciprocal space by bessel transform up to qmax = sqrt(ecutrho);
  !! * smoothed by multiplying with a filter function \(\textrm{filter}\)(q/qmax,a,nn);
  !! * brought back in real space,
  !
  !! where it overwrites the original array. The filter function is:
  !! \[ \text{filter}(x,a,\text{nn}) = e^{-\text{axx}} \sum_{k=0,\text{nn}} 
  !!                                   \frac{\text{axx}^k}{k!}\ . \]
  !
  USE upf_kinds,    ONLY: DP
  USE upf_io,       ONLY: stdout
  USE upf_const,    ONLY: fpi, sqrt2, eps8, eps6
  USE atom,         ONLY: rgrid
  USE uspp_param,   ONLY: nsp, upf, lmaxq, nbetam
  USE mp,           ONLY: mp_sum
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: ecutrho
  INTEGER,  INTENT(IN) :: intra_bgrp_comm
  !
  ! ... local variables
  !
  ! sdg 
  ! LOGICAL, PARAMETER :: tprint=.true. ! whether the q_l(r) and its relatives are printed or not
  LOGICAL, PARAMETER :: tprint=.FALSE.  ! whether the q_l(r) and its relatives are printed or not
  INTEGER, PARAMETER :: nn=16     ! smoothing parameter, order of the polynomial inverse gaussian
                                  ! approximant
  REAL(DP), PARAMETER:: a=22.0_DP ! smoothing parameter, exponent of the gaussian decaying factor
  REAL(DP), PARAMETER :: dq = 0.01_dp !nterpolation table 
                                  ! a=0.0 ; nn=0.0 would be no smoothing.
  !
  INTEGER :: nt, ih, jh, nb, mb, ijv, l, m, ir, ir0, iq, is, startq, lastq, ilast, ndm, ia
  ! various counters
  INTEGER :: nqxq
  REAL(DP), ALLOCATABLE :: qrad_q(:,:,:), qrad_r(:,:,:), qrad_rs(:,:,:), ffrr(:)
  REAL(DP), ALLOCATABLE :: power_0(:,:), power_r(:,:), power_q(:,:), power_rs(:,:), power_qs(:,:)
  REAL(DP), ALLOCATABLE :: aux(:), aux1(:)
  ! various work space
  REAL(DP) :: q, qi
  ! the modulus of g for each shell
  ! q-point grid for interpolation
  REAL(DP), ALLOCATABLE :: ylmk0(:)
  ! the spherical harmonics
  INTEGER :: lnb, lmb
  REAL(DP) :: qmax, rcut, drcut
  REAL(DP) :: target_ratio, ratio, ratio_s, fac
  !
  CHARACTER(LEN=6) :: filename
  !
  CALL start_clock( 'init_us_0' )
  !
  ! ... Initialization of the variables
  !
  drcut = ABS(LOG(eps8))/SQRT(ecutrho)
  !
  DO nt = 1, nsp
     !
     rcut = rgrid(nt)%r(upf(nt)%kkbeta)
     WRITE (stdout,*)  'RCUT:', rcut, drcut, rcut+drcut
     DO ir = upf(nt)%kkbeta, upf(nt)%mesh
        IF ( rgrid(nt)%r(ir) < rcut+drcut ) upf(nt)%kkbeta=ir
     ENDDO
     !
  ENDDO
  !
  ndm = MAXVAL( upf(:)%kkbeta )
  nqxq = INT( SQRT(ecutrho) / dq + 4 )
  !
  IF (tprint) THEN
     WRITE (stdout,*) " PSEUDOPOTENTIAL REPORT "
     WRITE (stdout,*) ' NDM :', ndm, '   ', upf(1:nsp)%kkbeta
     WRITE (stdout,*) ' LMAXQ :', lmaxq, ' NBETAM :', nbetam
  ENDIF
  !
  ALLOCATE( aux(ndm), aux1(ndm) )
  ALLOCATE( qrad_q(nqxq, nbetam*(nbetam+1)/2,lmaxq) )
  ALLOCATE( qrad_r(ndm , nbetam*(nbetam+1)/2,lmaxq) )
  ALLOCATE( qrad_rs(ndm, nbetam*(nbetam+1)/2,lmaxq) )
  ALLOCATE( ffrr(ndm) )    
  ALLOCATE( power_0(nbetam*(nbetam+1)/2 ,0:lmaxq) )
  ALLOCATE( power_r(nbetam*(nbetam+1)/2 ,0:lmaxq), power_q(nbetam*(nbetam+1)/2 ,0:lmaxq) )
  ALLOCATE( power_rs(nbetam*(nbetam+1)/2,0:lmaxq), power_qs(nbetam*(nbetam+1)/2,0:lmaxq) )
  ALLOCATE( ylmk0(lmaxq*lmaxq) )    
  !
  ! ... the following prevents an out-of-bound error: upf(nt)%nqlc=2*lmax+1
  ! but in some versions of the PP files lmax is not set to the maximum
  ! l of the beta functions but includes the l of the local potential
  !
  DO nt = 1, nsp
     !
     upf(nt)%nqlc = MIN( upf(nt)%nqlc, lmaxq )
     IF (upf(nt)%nqlc < 0)  upf(nt)%nqlc = 0
     !
  ENDDO
  !
  ! ... Here for the US types we compute the Fourier transform of the
  ! Q functions.
  !   
  CALL divide( intra_bgrp_comm, nqxq, startq, lastq )
  !
  qmax = SQRT(ecutrho)
  WRITE (stdout, *) ' qmax : sqrt(ecutrho) =', SQRT(ecutrho), dq*nqxq
  WRITE (stdout,'(a,f6.2,a,i4,4(a,f11.8))') 'FILTER : a=',a,', nn=',nn, &
                                ', filter(1.1d0)=', filter(1.1d0,a,nn), &
                                ', filter(1.0d0)=', filter(1.0d0,a,nn), &
                                ', filter(0.9d0)=', filter(0.9d0,a,nn), &
                                ', filter(0.8d0)=', filter(0.8d0,a,nn)
  !
  DO nt = 1, nsp
     WRITE (stdout,*) ' NT = ', nt
     !
     IF ( upf(nt)%tvanp ) THEN
        !-
        IF (tprint) THEN
           filename = 'qr_'   !  the radial q_l(r) as defined in the pseudopotential
           !
           DO nb = 1, upf(nt)%nbeta 
              DO mb = nb, upf(nt)%nbeta
                 !
                 ijv = mb * (mb - 1) / 2 + nb
                 lnb = upf(nt)%lll(nb) ; lmb = upf(nt)%lll(mb)
                 WRITE (filename(4:4),'(i1)') nb
                 WRITE (filename(5:5),'(i1)') mb
                 WRITE (filename(6:6),'(i1)') nt
                 OPEN (4,file=filename,form='formatted', status='unknown')
                 WRITE (4, *) '# the radial q_l(r) as defined in the pseudopotential'
                 WRITE (4, *) '# nb :', nb,lnb,' mb :',mb,lmb,' lmax :',lnb+lmb, &
                              ' kkbeta :',upf(nt)%kkbeta
                 !
                 DO ir =1,upf(nt)%kkbeta
                    write (4,'(12f16.10)') rgrid(nt)%r(ir), (upf(nt)%qfuncl(ir,ijv,l), l=0, lnb+lmb)
                 ENDDO
                 !
                 CLOSE (4)
                 !
              ENDDO
           ENDDO
           !
        ENDIF
        !-
        qrad_q(:,:,:) = 0.0_DP ; qrad_r(:,:,:) = 0.0_DP
        power_0(:,:)  = 0.0_DP ; power_r(:,:)  = 0.0_DP ; power_q(:,:) = 0.0_DP
        !
        DO l = 0, upf(nt)%nqlc-1
           !
           ! 1) first of all compute the integrated power spectrum of the Qs in real space.
           !
           DO nb = 1, upf(nt)%nbeta 
              DO mb = nb, upf(nt)%nbeta
                 !
                 ijv = mb * (mb - 1) / 2 + nb
                 aux1(1) = 0.0_DP
                 !
                 ir0 = 1
                 IF (rgrid(nt)%r(ir0) < eps8) ir0 = 2
                 !
                 DO ir = ir0, upf(nt)%kkbeta
                    aux1(ir) = (upf(nt)%qfuncl(ir,ijv,l)/rgrid(nt)%r(ir))**2
                 ENDDO
                 !
                 CALL simpson( upf(nt)%kkbeta, aux1, rgrid(nt)%rab, power_0(ijv,l+1) )
                 !
              ENDDO
           ENDDO
           !
           ! 2) compute the fourier transform of the Qs and their integrated power spectum in
           !    reciprocal space - FIXME: use routine init_tab_qrad
           !
           DO iq = startq, lastq
              !
              q = (iq - 1) * dq
              !
              ! ... here we compute the spherical bessel function for the given |q|
              !
              CALL sph_bes( upf(nt)%kkbeta, rgrid(nt)%r, q, l, aux )
              !
              ! .. and here we integrate with all the Q functions
              !
              DO nb = 1, upf(nt)%nbeta
                 DO mb = nb, upf(nt)%nbeta
                    !
                    ijv = mb * (mb - 1) / 2 + nb
                    lnb = upf(nt)%lll(nb)
                    lmb = upf(nt)%lll(mb)
                    !
                    IF ( (l >= ABS(lnb-lmb)) .AND. (l <= lnb+lmb) .AND. (MOD(l+lnb+lmb,2)==0) ) THEN
                       !
                       DO ir = 1, upf(nt)%kkbeta
                          aux1 (ir) = aux (ir) * upf(nt)%qfuncl(ir,ijv,l)
                       ENDDO
                       !
                       ! ... computes the Fourier transform of q_nb_mb_l(r) up to qmax
                       !
                       CALL simpson( upf(nt)%kkbeta, aux1, rgrid(nt)%rab, qrad_q(iq,ijv,l+1) )
                       !
                       ! ... and update the integrated power spectrum in reciprocal space
                       !
                       power_q(ijv,l+1) = power_q(ijv,l+1) + q*q * dq * qrad_q(iq,ijv,l+1)**2
                       !
                    ENDIF
                    !
                 ENDDO
              ENDDO
              !
              ! 3) back-fourier transform of the Qs to real space.
              !
              DO ir =1, upf(nt)%kkbeta
                 !
                 ! ... q_nb_mb_l(r) from the back fourier transform up to qmax of q_nb_mb_l(q)
                 !
                 qrad_r(ir,1:nbetam*(nbetam+1)/2,l+1) = qrad_r(ir,1:nbetam*(nbetam+1)/2,l+1) + &
                                               q*q * dq * aux(ir) * rgrid(nt)%r(ir)**2 * &
                                               qrad_q(iq,1:nbetam*(nbetam+1)/2,l+1)
              ENDDO
              !
           ENDDO
           !
        ENDDO
        !
        CALL mp_sum( qrad_q , intra_bgrp_comm )
        CALL mp_sum( power_q, intra_bgrp_comm ) ; power_q (:,:) = power_q (:,:) * 8.0_DP/fpi 
        CALL mp_sum( qrad_r , intra_bgrp_comm ) ; qrad_r(:,:,:) = qrad_r(:,:,:) * 8.0_DP/fpi 
        !
        ! 4) compute integrated power spectrum of the Qs in real space (completeness check).
        !
        DO l = 0, upf(nt)%nqlc-1
           !
           DO nb = 1, upf(nt)%nbeta 
              DO mb = nb, upf(nt)%nbeta
                 !
                 ijv = mb * (mb - 1)/2 + nb
                 aux1(1) = 0.0_DP
                 !
                 ir0 = 1
                 IF (rgrid(nt)%r(ir0) < eps8) ir0 = 2 
                 !
                 DO ir = ir0, upf(nt)%kkbeta
                    aux1(ir) = (qrad_r(ir,ijv,l+1)/rgrid(nt)%r(ir))**2
                 ENDDO
                 !
                 CALL simpson( upf(nt)%kkbeta, aux1, rgrid(nt)%rab, power_r(ijv,l+1) )
                 !
              ENDDO
           ENDDO
           !
        ENDDO
        !
        ! ... computes the ratio of power_0 and power_q to see how complete is the reciprocal
        ! space expansion.
        !
        power_0(:,0) = 0.0_DP
        power_q(:,0) = 0.0_DP
        power_r(:,0) = 0.0_DP
        target_ratio = 1.0_DP
        !
        DO nb = 1, upf(nt)%nbeta
           DO mb = nb, upf(nt)%nbeta
              !
              ijv = mb * (mb-1) / 2 + nb
              lnb = upf(nt)%lll(nb)
              lmb = upf(nt)%lll(mb)
              !
              DO l = 0, lnb+lmb
                 power_0(ijv,0) = power_0(ijv,0) + power_0(ijv,l+1)
                 power_q(ijv,0) = power_q(ijv,0) + power_q(ijv,l+1)
                 power_r(ijv,0) = power_r(ijv,0) + power_r(ijv,l+1)
              ENDDO
              !
              !write (stdout, *) ' nb :', nb,lnb,' mb :',mb,lmb,' lmax :',lnb+lmb
              !write (stdout,'(a,12f16.10)') 'power_0 ',(power_0 (ijv,l+1), l=0,lnb+lmb)
              !write (stdout,'(a,12f16.10)') 'power_r ',(power_r (ijv,l+1), l=0,lnb+lmb)
              !write (stdout,'(a,12f16.10)') 'power_q ',(power_q (ijv,l+1), l=0,lnb+lmb)
              !write (stdout,*) 'ratio   ',1.d0-(power_r (ijv,0)/power_q (ijv,0)), &
              !                                1.d0-(power_r (ijv,0)/power_0 (ijv,0))
              IF (power_0(ijv,0)>eps8) target_ratio = MIN(target_ratio, power_r(ijv,0)/power_0(ijv,0))
              !
           ENDDO
        ENDDO
        !
        WRITE (stdout,*) ' TARGET Qs SPILLOVER : 1.d0-target_ratio, target_ratio ', &
                                                    1.d0-target_ratio, target_ratio
        !
        fac = 1.2_DP
        !
 99     CONTINUE
        !
        ! ... smooth the Fourier transform of the Qs and compute its integrated power spectrum
        !
        power_qs(:,:) = 0.0_DP
        !
        DO nb = 1, upf(nt)%nbeta
           DO mb = nb, upf(nt)%nbeta
              !
              ijv = mb * (mb - 1) / 2 + nb
              lnb = upf(nt)%lll(nb)
              lmb = upf(nt)%lll(mb)
              !
              DO l = 0, upf(nt)%nqlc-1
                 !
                 IF ( (l >= ABS(lnb-lmb)) .AND. (l <= lnb+lmb) .AND. (MOD(l+lnb+lmb,2)==0) ) THEN
                    DO iq = startq, lastq
                       q = (iq - 1) * dq
                       power_qs(ijv,l+1) = power_qs(ijv,l+1) + q*q * dq * &
                                           (qrad_q(iq,ijv,l+1)*filter(fac*q/qmax,a,nn))**2
                    ENDDO
                 ENDIF
                 !
                 power_qs(ijv,0) = power_qs(ijv,0) + power_qs(ijv,l+1)
                 !
              ENDDO
              !
           ENDDO
        ENDDO
        !
        power_qs(:,:) = power_qs(:,:) * 8.0_DP/fpi
        !
        CALL mp_sum( power_qs, intra_bgrp_comm )
        !
        ratio = 1.0_DP ;  ratio_s = 1.0_DP
        !
        DO nb = 1, upf(nt)%nbeta
           DO mb = nb, upf(nt)%nbeta
              ijv = mb * (mb - 1) / 2 + nb
              IF (power_q(ijv,0)>eps8) ratio = MIN(ratio, power_qs(ijv,0)/power_q(ijv,0))
           ENDDO
        ENDDO
        !
        ! WRITE (stdout,*) ' filter factor and power_qs/power_q ratio:', fac, ratio
        ! 
        ! ... with a given fac the smoothed Qs qre built
        !        
        qrad_rs(:,:,:) = 0.0_DP ; ffrr(:) = 0.0_DP
        !
        DO l = 0, upf(nt)%nqlc-1
           !
           DO iq = startq, lastq
              !
              q = (iq - 1) * dq
              !
              ! ... here we compute the spherical bessel function for the given |q| ...
              !
              CALL sph_bes( upf(nt)%kkbeta, rgrid(nt)%r, q, l, aux )
              !
              ! ... and here we integrate with all the Q functions
              !
              ! 5) back-fourier transform of the smoothed Qs to real space.
              !
              DO ir = 1, upf(nt)%kkbeta
                 !
                 ! ... q_nb_mb_l(r) from the back fourier transform up to qmax of q_nb_mb_l(q)
                 !
                 qrad_rs(ir,1:nbetam*(nbetam+1)/2,l+1) = qrad_rs(ir,1:nbetam*(nbetam+1)/2,l+1) &
                                 + aux(ir) * q*q * dq * rgrid(nt)%r(ir)**2               &
                                 * qrad_q(iq,1:nbetam*(nbetam+1)/2,l+1) * filter(fac*q/qmax,a,nn)
                 !
                 ! ... build the filter function in real space from the back fourier transform up
                 ! to qmax 
                 !
                 IF (l==0) ffrr(ir) = ffrr(ir) + q*q * dq * aux(ir) * rgrid(nt)%r(ir)**2 &
                                                                      * filter(fac*q/qmax,a,nn)
                 !
              ENDDO
              !
           ENDDO
           !
        ENDDO
        !
        CALL mp_sum( qrad_rs, intra_bgrp_comm ) ; qrad_rs(:,:,:) = qrad_rs(:,:,:) * 8.0_DP/fpi 
        CALL mp_sum( ffrr, intra_bgrp_comm ) ; ffrr(:) = ffrr(:) * 8.0_DP/fpi
        !
        ! 6) compute intergrated power spectrum of the Qs in real space (completeness check).
        !
        DO l = 0, upf(nt)%nqlc-1
           !
           DO nb = 1, upf(nt)%nbeta 
              DO mb = nb, upf(nt)%nbeta
                 !
                 ijv = mb * (mb-1)/2 + nb
                 aux1(1) = 0.0_DP
                 !
                 ir0 = 1
                 IF (rgrid(nt)%r(ir0) < eps8) ir0 = 2
                 !
                 DO ir = ir0, upf(nt)%kkbeta
                    aux1(ir) = (qrad_rs(ir,ijv,l+1)/rgrid(nt)%r(ir))**2
                 ENDDO
                 !
                 CALL simpson( upf(nt)%kkbeta, aux1, rgrid(nt)%rab, power_rs(ijv,l+1) )
                 !
              ENDDO
           ENDDO
           !
        ENDDO
        !
        power_rs(:,0) = 0.0_DP
        !
        DO nb = 1, upf(nt)%nbeta
           DO mb = nb, upf(nt)%nbeta
              !
              ijv = mb * (mb-1)/2 + nb
              lnb = upf(nt)%lll(nb)
              lmb = upf(nt)%lll(mb)
              !
              DO l = 0, lnb+lmb
                 power_rs(ijv,0) = power_rs(ijv,0) + power_rs(ijv,l+1)
              ENDDO
              !
              !write (stdout, *) ' nb :', nb,lnb,' mb :',mb,lmb,' lmax :',lnb+lmb
              !write (stdout,'(a,12f16.10)') 'power_rs',(power_rs(ijv,l+1), l=0,lnb+lmb)
              !write (stdout,'(a,12f16.10)') 'power_qs',(power_qs(ijv,l+1), l=0,lnb+lmb)
              !write (stdout,*) 'ratio   ',1.d0-(power_rs(ijv,0)/power_qs(ijv,0)), &
              !                            1.d0-(power_qs(ijv,0)/power_q (ijv,0))
              IF (power_qs(ijv,0) > eps8) ratio_s = MIN(ratio_s, power_rs(ijv,0)/power_qs(ijv,0))
              !
           ENDDO
        ENDDO
        !
        WRITE (stdout,*) ' filter factor, 1.d0-power_rs/power_qs and power_qs/power_q ratio :', &
                                                                      fac, 1.d0-ratio_s, ratio
        fac = fac - 0.05_DP
        !sdg
        IF (ratio < target_ratio .AND. (1.d0-ratio_s) < eps8) GOTO 99
        !
        ! IF (ratio < target_ratio .AND. (1.d0-ratio_s) < eps6) GOTO 99
        !
        fac = fac + 0.05_DP ! reset the last successful value of fac 
        !
        !- save the smoothed real space qs in qfuncl
        !
        DO nb = 1, upf(nt)%nbeta
           DO mb = nb, upf(nt)%nbeta
              !
              ijv = mb * (mb-1)/2 + nb
              lnb = upf(nt)%lll(nb)
              lmb = upf(nt)%lll(mb)
              !
              DO l = 0, upf(nt)%nqlc-1
                 !
                 IF ( (l>=ABS(lnb-lmb)) .AND. (l<=lnb+lmb) .AND. (MOD(l+lnb+lmb,2)==0) ) THEN
                    !
                    upf(nt)%qfuncl(1:upf(nt)%kkbeta,ijv,l) = qrad_rs(1:upf(nt)%kkbeta,ijv,l+1)
                    !
                 ENDIF
                 !
              ENDDO
              !
           ENDDO
        ENDDO
        !
     ENDIF
     !
     IF (tprint) THEN
        !-
        filename = 'ffff'
        OPEN (4, FILE=filename, FORM='formatted', STATUS='unknown')
        WRITE (4, *) '# filter function : a=',a,', nn=',nn,', fac=', fac
        WRITE (4, *) '# nqxq :', nqxq,' dq :',dq, ' qmax :',qmax
        DO iq = 1, nqxq
           q = (iq-1)*dq
           WRITE (4,'(2f16.10)')  q, filter( fac*q/qmax, a, nn )
        ENDDO
        CLOSE (4)
        !-
        filename = 'ffrr'
        OPEN (4, FILE=filename, FORM='formatted', STATUS='unknown')
        WRITE (4, *) '# filter function : a=',a,', nn=',nn,', fac=', fac
        WRITE (4, *) '# kkbeta :', upf(nt)%kkbeta
        DO ir = 1, upf(nt)%kkbeta
           write (4,'(2f16.10)')  rgrid(nt)%r(ir), ffrr(ir)
        ENDDO
        CLOSE (4)
        !-
        DO nb = 1, upf(nt)%nbeta 
           DO mb = nb, upf(nt)%nbeta
              !
              ijv = mb * (mb - 1) / 2 + nb
              lnb = upf(nt)%lll(nb) ; lmb = upf(nt)%lll(mb)
              WRITE (filename(4:4),'(i1)') nb
              WRITE (filename(5:5),'(i1)') mb
              WRITE (filename(6:6),'(i1)') nt
              !-
              filename(1:3) = 'qq_'    ! the radial fourier transform of q_l in reciprocal space
              OPEN (4, FILE=filename, FORM='formatted', STATUS='unknown')
              WRITE (4,*) '# the radial fourier transform of q_l in reciprocal space'
              WRITE (4,*) '# nb :', nb, lnb,' mb :', mb, lmb,' lmax :', lnb+lmb, ' nqxq :', nqxq
              DO iq=1,nqxq
                 q = (iq-1)*dq
                 WRITE (4,'(12f16.10)')  q, (qrad_q(iq,ijv,l+1), l=0,lnb+lmb )
              ENDDO
              CLOSE (4)
              !-
              filename(1:3) = 'qqs'    ! the smoothed radial fourier transform of q_l in reciprocal space
              OPEN (4, FILE=filename, FORM='formatted', STATUS='unknown')
              WRITE (4,*) '# the smoothed radial fourier transform of q_l in reciprocal space'
              WRITE (4,*) '# nb :', nb,lnb,' mb :',mb,lmb,' lmax :',lnb+lmb, ' nqxq :',nqxq
              DO iq = 1, nqxq
                 q = (iq-1)*dq
                 WRITE (4,'(12f16.10)')  q,(qrad_q(iq,ijv,l+1)*filter(fac*q/qmax,a,nn), l=0,lnb+lmb )
              ENDDO
              CLOSE (4)
              !-
              filename(1:3) = 'qrq'    ! the radial q_l(r) as obtained back-transforming the q_l(q)
              OPEN (4, FILE=filename, FORM='formatted', STATUS='unknown')
              WRITE (4,*) '# the radial q_l(r) as obtained back-transforming the q_l(q)'
              WRITE (4,*) '# nb :', nb,lnb,' mb :',mb,lmb,' lmax :',lnb+lmb, ' kkbeta :',upf(nt)%kkbeta
              DO ir = 1, upf(nt)%kkbeta
                 WRITE (4,'(12f16.10)')  rgrid(nt)%r(ir),(qrad_r(ir,ijv,l+1), l=0,lnb+lmb )
              ENDDO
              CLOSE (4)
              !-
              filename(1:3) = 'qrs'    ! the radial q_l(r) as obtained back-transforming the smoothed q_l(q)
              OPEN (4, FILE=filename, FORM='formatted', STATUS='unknown')
              WRITE (4,*) '# the radial q_l(r) as obtained back-transforming the smoothed q_l(q)'
              WRITE (4,*) '# nb :', nb,lnb,' mb :',mb,lmb,' lmax :',lnb+lmb, ' kkbeta :',upf(nt)%kkbeta
              DO ir = 1, upf(nt)%kkbeta
                 WRITE (4,'(12f16.10)')  rgrid(nt)%r(ir),(qrad_rs(ir,ijv,l+1), l=0,lnb+lmb )
              ENDDO
              CLOSE (4)
              !-
           ENDDO
        ENDDO
        !
     ENDIF
     ! nsp
  ENDDO
  !
  DEALLOCATE( ylmk0 )
  DEALLOCATE( power_0, power_r, power_q, power_rs, power_qs )
  DEALLOCATE( qrad_q, qrad_r, qrad_rs, ffrr )
  DEALLOCATE( aux1, aux )
  !
  CALL stop_clock( 'init_us_0' )
  !
  RETURN
  !
  !
  !
  CONTAINS
  !
  REAL(DP) FUNCTION filter( x, a, n )
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN) :: x, a
    INTEGER, INTENT(IN) :: n
    REAL(DP) :: axx, ff
    INTEGER :: k
    !
    axx = a * x * x
    ff = 1.0_DP
    !
    DO k = n, 1, -1
       ff = (1._DP + axx/DBLE(k)*ff)
    ENDDO
    !
    filter = ff * EXP(-axx)
    !
    RETURN
    !
  END FUNCTION filter
  !
  !
END SUBROUTINE init_us_0
