
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine init_us_0
  !----------------------------------------------------------------------
  !
  !   This routine performs the following tasks:
  !   for each uspp or paw pseudopotential the l-dependent aumentation charge
  !   q_nb_mb_l(r), stored in qfuncl(ir,nmb,l), is 
  !   - transformed in reciprocal space by bessel transform up to qmax = sqrt(ecutrho), 
  !   - smoothed by multiplying with a filter function filter(q/qmax,a,nn) and
  !   - brought back in real space 
  !   where it overwrites the original array.
  !
  !   the filter function is     filter(x,a,nn) = exp(-axx) * \sum_{k=0,nn} (axx)**k/k!
  !
  USE kinds,        ONLY : DP
  USE gvect,        ONLY : ecutrho
  USE io_global,    ONLY : stdout
  USE constants,    ONLY : fpi, sqrt2, eps8
  USE atom,         ONLY : rgrid
  USE ions_base,    ONLY : ntyp => nsp
  USE cell_base,    ONLY : omega, tpiba
  USE us,           ONLY : nqxq, dq
  USE uspp_param,   ONLY : upf, lmaxq, nbetam
  USE mp_bands,     ONLY : intra_bgrp_comm
  USE mp,           ONLY : mp_sum
  !
  implicit none
  !
  !     here a few local variables
  !
  logical, parameter :: tprint=.true.    ! whether the q_l(r) and its relatives are printed or not
  integer, parameter :: nn=16   ! smoothing parameter, order of the polynomial inverse gaussian approximant
  real(DP), parameter:: a=22.0  ! smoothing parameter, exponent of the gaussian decaying factor
                                ! a=0.d0 ; nn=0 would be no smoothing.
  !
  integer :: nt, ih, jh, nb, mb, ijv, l, m, ir, ir0, iq, is, startq, lastq, ilast, ndm, ia
  ! various counters
  real(DP), allocatable :: qrad_q(:,:,:), qrad_r(:,:,:), qrad_rs(:,:,:), ffrr(:)
  real(DP), allocatable :: power_0(:,:), power_r(:,:), power_q(:,:), power_rs(:,:), power_qs(:,:)
  real(DP), allocatable :: aux (:), aux1 (:)
  ! various work space
  real(DP) :: q, qi
  ! the modulus of g for each shell
  ! q-point grid for interpolation
  real(DP), allocatable :: ylmk0 (:)
  ! the spherical harmonics
  integer, external :: sph_ind
  integer :: lnb, lmb
  real(DP) :: qmax, rcut, drcut
  real(DP) :: target_ratio, ratio, ratio_s, fac
  !
  character(LEN=6) :: filename
  !
  call start_clock ('init_us_0')
  !
  !    Initialization of the variables
  !
  drcut = abs(log(eps8))/sqrt(ecutrho)
  do nt=1,ntyp
     rcut = rgrid(nt)%r(upf(nt)%kkbeta)
     write (stdout,*)  'RCUT :', rcut, drcut, rcut+drcut
     do ir = upf(nt)%kkbeta, upf(nt)%mesh
        if ( rgrid(nt)%r(ir) < rcut + drcut ) upf(nt)%kkbeta=ir
     end do
  end do

  ndm = MAXVAL ( upf(:)%kkbeta )
  if (tprint) then
     write ( stdout, * ) " PSEUDOPOTENTIAL REPORT "
     write ( stdout, * ) ' NDM :', ndm, '   ', upf(1:ntyp)%kkbeta
     write ( stdout, * ) ' LMAXQ :', lmaxq, ' NBETAM :', nbetam
  end if

  allocate (aux ( ndm), aux1( ndm))    
  allocate (qrad_q (nqxq, nbetam*(nbetam+1)/2,lmaxq ) )
  allocate (qrad_r ( ndm, nbetam*(nbetam+1)/2,lmaxq ) )
  allocate (qrad_rs( ndm, nbetam*(nbetam+1)/2,lmaxq ) )
  allocate (ffrr( ndm))    
  allocate (power_0 ( nbetam*(nbetam+1)/2,0:lmaxq ) )
  allocate (power_r ( nbetam*(nbetam+1)/2,0:lmaxq ), power_q ( nbetam*(nbetam+1)/2,0:lmaxq ) )
  allocate (power_rs( nbetam*(nbetam+1)/2,0:lmaxq ), power_qs( nbetam*(nbetam+1)/2,0:lmaxq ) )
  allocate (ylmk0( lmaxq * lmaxq))    
  !
  ! the following prevents an out-of-bound error: upf(nt)%nqlc=2*lmax+1
  ! but in some versions of the PP files lmax is not set to the maximum
  ! l of the beta functions but includes the l of the local potential
  !
  do nt=1,ntyp
     upf(nt)%nqlc = MIN ( upf(nt)%nqlc, lmaxq )
     IF ( upf(nt)%nqlc < 0 )  upf(nt)%nqlc = 0
  end do

  !
  !   here for the US types we compute the Fourier transform of the
  !   Q functions.
  !   
  call divide (intra_bgrp_comm, nqxq, startq, lastq)
  qmax = sqrt(ecutrho)
  write (6, *) ' qmax : sqrt(ecutrho) =',sqrt(ecutrho), dq*nqxq*tpiba, tpiba 
  write (6,'(a,f6.2,a,i4,4(a,f11.8))') 'FILTER : a=',a,', nn=',nn,', filter(1.1d0)=', filter(1.1d0,a,nn), &
                                                                  ', filter(1.0d0)=', filter(1.0d0,a,nn), &
                                                                  ', filter(0.9d0)=', filter(0.9d0,a,nn), &
                                                                  ', filter(0.8d0)=', filter(0.8d0,a,nn)
  !
  do nt = 1, ntyp
     write (*,*) ' NT = ', nt
     !
     if ( upf(nt)%tvanp ) then
!-
        if (tprint) then
           filename = 'qr_'                 !  the radial q_l(r) as defined in the pseudopotential
           do nb=1, upf(nt)%nbeta 
              do mb=nb, upf(nt)%nbeta
                 ijv = mb * (mb - 1) / 2 + nb
                 lnb = upf(nt)%lll(nb) ; lmb = upf(nt)%lll(mb)
                 WRITE (filename( 4:4 ),'(i1)') nb ; WRITE (filename( 5:5 ),'(i1)') mb ; WRITE (filename( 6:6 ),'(i1)') nt
                 OPEN (4,file=filename,form='formatted', status='unknown')
                 write (4, *) '# the radial q_l(r) as defined in the pseudopotential'
                 write (4, *) '# nb :', nb,lnb,' mb :',mb,lmb,' lmax :',lnb+lmb, ' kkbeta :',upf(nt)%kkbeta
                 do ir =1,upf(nt)%kkbeta
                    write (4,'(12f16.10)') rgrid(nt)%r(ir), (upf(nt)%qfuncl(ir,ijv,l), l=0,lnb+lmb)
                 end do
                 CLOSE (4)
              end do
           end do
        end if
!-
        qrad_q(:,:,:)= 0.d0 ; qrad_r(:,:,:)= 0.d0 
        power_0(:,:) = 0.d0; power_r(:,:) = 0.d0; power_q(:,:) = 0.d0

        do l = 0, upf(nt)%nqlc -1
           !
           ! 1) first of all compute the integrated power spectrum of the Qs in real space
           !
           do nb=1, upf(nt)%nbeta 
              do mb=nb, upf(nt)%nbeta
                 ijv = mb * (mb - 1) / 2 + nb
                 aux1(1)=0.d0; ir0 =1; if (rgrid(nt)%r(ir0) < eps8) ir0=2 
                 do ir = ir0, upf(nt)%kkbeta
                    aux1 (ir) = (upf(nt)%qfuncl(ir,ijv,l)/rgrid(nt)%r(ir))**2
                 end do
                 call simpson ( upf(nt)%kkbeta, aux1, rgrid(nt)%rab, power_0(ijv,l+1) )
              end do
           end do
           !
           ! 2) compute the fourier transform of the Qs and their itegrated power spectum in reciprocal space
           !
           do iq = startq, lastq
              q = (iq - 1) * dq * tpiba
              ! here we compute the spherical bessel function for the given |q| ...
              call sph_bes ( upf(nt)%kkbeta, rgrid(nt)%r, q, l, aux)
              ! .. and here we integrate with all the Q functions
              do nb = 1, upf(nt)%nbeta
                 do mb = nb, upf(nt)%nbeta
                    ijv = mb * (mb - 1) / 2 + nb; lnb = upf(nt)%lll(nb) ; lmb = upf(nt)%lll(mb)
                    if ( ( l >= abs(lnb - lmb) ) .and. ( l <= lnb + lmb ) .and. (mod(l+lnb+lmb,2)==0) ) then
                       do ir = 1, upf(nt)%kkbeta
                          aux1 (ir) = aux (ir) * upf(nt)%qfuncl(ir,ijv,l)
                       enddo
                       ! compute the fourier transform of q_nb_mb_l(r) up to qmax ...
                       call simpson ( upf(nt)%kkbeta, aux1, rgrid(nt)%rab, qrad_q(iq,ijv,l+1) )
                       ! ... and update the integrated power spectrum in reciprocal space
                       power_q (ijv,l+1) = power_q (ijv,l+1) + q*q * dq * tpiba *  qrad_q(iq,ijv,l+1)**2
                    endif
                 enddo
              enddo
              !
              ! 3) back-fourier trasform of the Qs to real space
              !
              do ir =1, upf(nt)%kkbeta
                 ! q_nb_mb_l(r) from the back fourier transform up to qmax of q_nb_mb_l(q)
                 qrad_r(ir,1:nbetam*(nbetam+1)/2,l+1) = qrad_r(ir,1:nbetam*(nbetam+1)/2,l+1) + &
                                    q*q * dq*tpiba * aux(ir) * rgrid(nt)%r(ir)**2 * qrad_q( iq, 1:nbetam*(nbetam+1)/2,l+1) 
              enddo
           enddo
        enddo
        call mp_sum ( qrad_q , intra_bgrp_comm )
        call mp_sum ( power_q, intra_bgrp_comm ) ; power_q (:,:) = power_q (:,:) * 8.d0/fpi 
        call mp_sum ( qrad_r , intra_bgrp_comm ) ; qrad_r(:,:,:) = qrad_r(:,:,:) * 8.d0/fpi 
        !
        ! 4) compute intergrated power spectrum of the Qs in real space (completeness check)
        !
        do l = 0, upf(nt)%nqlc -1
           do nb=1, upf(nt)%nbeta 
              do mb=nb, upf(nt)%nbeta
                 ijv = mb * (mb - 1) / 2 + nb
                 aux1(1)=0.d0; ir0 =1; if (rgrid(nt)%r(ir0) < eps8) ir0=2 
                 do ir = ir0, upf(nt)%kkbeta
                    aux1 (ir) = (qrad_r(ir,ijv,l+1)/rgrid(nt)%r(ir))**2
                 end do
                 call simpson ( upf(nt)%kkbeta, aux1, rgrid(nt)%rab, power_r(ijv,l+1) )
              end do
           end do
           ! l
        enddo
!
! compute the ratio of power_0 and power_q to see how complete is the reciprocal space expansion
!
        power_0 (:,0) = 0.d0 ; power_q (:,0) = 0.d0 ; power_r(:,0) = 0.d0 ; target_ratio = 1.0d0
        do nb = 1, upf(nt)%nbeta
           do mb = nb, upf(nt)%nbeta
              ijv = mb * (mb-1) / 2 + nb
              lnb = upf(nt)%lll(nb) ; lmb = upf(nt)%lll(mb)
              do l=0,lnb+lmb
                 power_0 (ijv,0) = power_0 (ijv,0) + power_0 (ijv,l+1)
                 power_q (ijv,0) = power_q (ijv,0) + power_q (ijv,l+1)
                 power_r (ijv,0) = power_r (ijv,0) + power_r (ijv,l+1)
              end do
              write (*, *) ' nb :', nb,lnb,' mb :',mb,lmb,' lmax :',lnb+lmb
              write (*,'(a,12f16.10)') 'power_0 ',(power_0 (ijv,l+1), l=0,lnb+lmb)
              write (*,'(a,12f16.10)') 'power_r ',(power_r (ijv,l+1), l=0,lnb+lmb)
              write (*,'(a,12f16.10)') 'power_q ',(power_q (ijv,l+1), l=0,lnb+lmb)
              write (*,*) 'ratio   ',1.d0-(power_r (ijv,0)/power_q (ijv,0)), 1.d0-(power_r (ijv,0)/power_0 (ijv,0))
              if (power_0(ijv,0).gt.eps8) target_ratio = min(target_ratio,power_r(ijv,0)/power_0(ijv,0))
           enddo ! mb
        enddo ! nb
        write (*,*) ' TARGET Qs SPILLOVER :',1.d0-target_ratio, target_ratio

!
        fac = 1.2d0
 99     continue
        ! smooth the fourier trasform of the Qs and compute its integrated power spectrum
        power_qs(:,:) = 0.d0
        do nb = 1, upf(nt)%nbeta
           do mb = nb, upf(nt)%nbeta
              ijv = mb * (mb - 1) / 2 + nb; lnb = upf(nt)%lll(nb) ; lmb = upf(nt)%lll(mb)
              do l = 0, upf(nt)%nqlc -1
                 if ( ( l >= abs(lnb - lmb) ) .and. ( l <= lnb + lmb ) .and. (mod(l+lnb+lmb,2)==0) ) then
                    do iq = startq, lastq
                       q = (iq - 1) * dq * tpiba
                       power_qs (ijv,l+1) = power_qs (ijv,l+1) + &
                                            q*q * dq * tpiba *  (qrad_q(iq,ijv,l+1)*filter(fac*q/qmax,a,nn))**2
                    end do
                 end if
                 power_qs (ijv,0) = power_qs (ijv,0) + power_qs (ijv,l+1)
              end do
           end do
        end do
        power_qs(:,:) = power_qs(:,:) * 8.d0/fpi ; call mp_sum ( power_qs , intra_bgrp_comm )

        ratio = 1.0d0;  ratio_s = 1.0d0
        do nb = 1, upf(nt)%nbeta
           do mb = nb, upf(nt)%nbeta
              ijv = mb * (mb - 1) / 2 + nb
              if (power_q(ijv,0).gt.eps8) ratio = min (ratio, power_qs(ijv,0)/power_q(ijv,0))
           end do
        end do
        WRITE (*,*) ' filter factor and ratio :', fac, ratio
! 
! once fac is chosen the smoothed Qs qre built
!        
        qrad_rs(:,:,:) = 0.d0; ffrr(:) = 0.d0
        do l = 0, upf(nt)%nqlc -1
           do iq = startq, lastq
              q = (iq - 1) * dq * tpiba
              ! here we compute the spherical bessel function for the given |q| ...
              call sph_bes ( upf(nt)%kkbeta, rgrid(nt)%r, q, l, aux)
              ! .. and here we integrate with all the Q functions
              !
              ! 5) back-fourier trasform of the smoothed Qs to real space
              !
              do ir =1, upf(nt)%kkbeta
                 ! q_nb_mb_l(r) from the back fourier transform up to qmax of q_nb_mb_l(q)
                 qrad_rs(ir,1:nbetam*(nbetam+1)/2,l+1) = qrad_rs(ir,1:nbetam*(nbetam+1)/2,l+1) + aux(ir) *&
                   q*q * dq*tpiba * rgrid(nt)%r(ir)**2 * qrad_q( iq, 1:nbetam*(nbetam+1)/2,l+1) * filter(fac*q/qmax,a,nn)
                 
                 ! build the filter function in real space from the back fourier transform up to qmax 
                 if (l==0) ffrr(ir) = ffrr(ir) + q*q * dq*tpiba * aux(ir) * rgrid(nt)%r(ir)**2 * filter(fac*q/qmax,a,nn)
              enddo
           enddo
        enddo
        call mp_sum ( qrad_rs, intra_bgrp_comm ) ; qrad_rs(:,:,:) = qrad_rs(:,:,:) * 8.d0/fpi 
        call mp_sum ( ffrr, intra_bgrp_comm ) ; ffrr(:) = ffrr(:) * 8.d0/fpi
        !
        ! 6) compute intergrated power spectrum of the Qs in real space (completeness check)
        !
        do l = 0, upf(nt)%nqlc -1
           do nb=1, upf(nt)%nbeta 
              do mb=nb, upf(nt)%nbeta
                 ijv = mb * (mb - 1) / 2 + nb
                 aux1(1)=0.d0; ir0 =1; if (rgrid(nt)%r(ir0) < eps8) ir0=2 
                 do ir = ir0, upf(nt)%kkbeta
                    aux1 (ir) = (qrad_rs(ir,ijv,l+1)/rgrid(nt)%r(ir))**2
                 end do
                 call simpson ( upf(nt)%kkbeta, aux1, rgrid(nt)%rab, power_rs(ijv,l+1) )
              end do
           end do
           ! l
        enddo
        power_rs (:,0) = 0.d0
        do nb = 1, upf(nt)%nbeta
           do mb = nb, upf(nt)%nbeta
              ijv = mb * (mb-1) / 2 + nb
              lnb = upf(nt)%lll(nb) ; lmb = upf(nt)%lll(mb)
              do l=0,lnb+lmb
                 power_rs (ijv,0) = power_rs (ijv,0) + power_rs (ijv,l+1)
              end do
              write (*, *) ' nb :', nb,lnb,' mb :',mb,lmb,' lmax :',lnb+lmb
              write (*,'(a,12f16.10)') 'power_rs',(power_rs(ijv,l+1), l=0,lnb+lmb)
              write (*,'(a,12f16.10)') 'power_qs',(power_qs(ijv,l+1), l=0,lnb+lmb)
              write (*,*) 'ratio   ',1.d0-(power_rs(ijv,0)/power_qs(ijv,0)), 1.d0-(power_qs(ijv,0)/power_q (ijv,0))
              if (power_qs(ijv,0).gt.eps8) ratio_s = min (ratio_s, power_rs(ijv,0)/power_qs(ijv,0))
           enddo ! mb
        enddo ! nb

        WRITE (*,*) ' filter factor and ratio :', fac, ratio_s, ratio
        fac = fac -0.05d0
        if (ratio .lt. target_ratio .and. (1.d0-ratio_s).lt. eps8) go to 99
        fac = fac + 0.05d0

!- save the smoothed real space qs in qfuncl
        do nb = 1, upf(nt)%nbeta
           do mb = nb, upf(nt)%nbeta
              ijv = mb * (mb-1) / 2 + nb;  lnb = upf(nt)%lll(nb) ; lmb = upf(nt)%lll(mb)
              do l = 0, upf(nt)%nqlc -1
                  if ( (l>=abs(lnb-lmb)) .and. (l<=lnb+lmb) .and. (mod(l+lnb+lmb,2)==0) ) then
                     upf(nt)%qfuncl(1:upf(nt)%kkbeta,ijv,l) = qrad_r(1:upf(nt)%kkbeta,ijv,l+1) 
                  end if
              enddo
           enddo ! mb
        enddo ! nb
     end if

     if (tprint) then
!-
        filename = 'ffff'
        OPEN (4,file=filename,form='formatted', status='unknown')
        write (4, *) '# filter function : a=',a,', nn=',nn,', fac=', fac
        write (4, *) '# nqxq :', nqxq,' dq :',dq, ' qmax :',qmax
        do iq=1,nqxq
           q = (iq-1)*dq*tpiba
           write (4,'(2f16.10)')  q,filter(fac*q/qmax,a,nn)
        end do
        CLOSE (4)
!-
        filename = 'ffrr'
        OPEN (4,file=filename,form='formatted', status='unknown')
        write (4, *) '# filter function : a=',a,', nn=',nn,', fac=', fac
        write (4, *) '# kkbeta :', upf(nt)%kkbeta
        do ir=1,upf(nt)%kkbeta
           write (4,'(2f16.10)')  rgrid(nt)%r(ir),ffrr(ir)
        end do
        CLOSE (4)
!-
        do nb=1, upf(nt)%nbeta 
           do mb=nb, upf(nt)%nbeta
              ijv = mb * (mb - 1) / 2 + nb
              lnb = upf(nt)%lll(nb) ; lmb = upf(nt)%lll(mb)
              WRITE (filename( 4:4 ),'(i1)') nb ; WRITE (filename( 5:5 ),'(i1)') mb ; WRITE (filename( 6:6 ),'(i1)') nt
!-
              filename(1:3) = 'qq_'                 ! the radial fourier transform of q_l in reciprcal space
              OPEN (4,file=filename,form='formatted', status='unknown')
              write (4, *) '# the radial fourier transform of q_l in reciprcal space'
              write (4, *) '# nb :', nb,lnb,' mb :',mb,lmb,' lmax :',lnb+lmb, ' nqxq :',nqxq
              do iq=1,nqxq
                 q = (iq-1)*dq*tpiba
                 write (4,'(12f16.10)')  q,(qrad_q(iq,ijv,l+1), l=0,lnb+lmb )
              end do
              CLOSE (4)
!-
              filename(1:3) = 'qqs'                 ! the smoothed radial fourier transform of q_l in reciprcal space
              OPEN (4,file=filename,form='formatted', status='unknown')
              write (4, *) '# the smoothed radial fourier transform of q_l in reciprcal space'
              write (4, *) '# nb :', nb,lnb,' mb :',mb,lmb,' lmax :',lnb+lmb, ' nqxq :',nqxq
              do iq=1,nqxq
                 q = (iq-1)*dq*tpiba
                 write (4,'(12f16.10)')  q,(qrad_q(iq,ijv,l+1)*filter(fac*q/qmax,a,nn), l=0,lnb+lmb )
              end do
              CLOSE (4)
!-
              filename(1:3) = 'qrq'                 ! the radial q_l(r) as obtained back-transforming the q_l(q)
              OPEN (4,file=filename,form='formatted', status='unknown')
              write (4, *) '# the radial q_l(r) as obtained back-transforming the q_l(q)'
              write (4, *) '# nb :', nb,lnb,' mb :',mb,lmb,' lmax :',lnb+lmb, ' kkbeta :',upf(nt)%kkbeta
              do ir=1,upf(nt)%kkbeta
                 write (4,'(12f16.10)')  rgrid(nt)%r(ir),(qrad_r(ir,ijv,l+1), l=0,lnb+lmb )
              end do
              CLOSE (4)
!-
              filename(1:3) = 'qrs'                 ! the radial q_l(r) as obtained back-transforming the smoothed q_l(q)
              OPEN (4,file=filename,form='formatted', status='unknown')
              write (4, *) '# the radial q_l(r) as obtained back-transforming the smoothed q_l(q)'
              write (4, *) '# nb :', nb,lnb,' mb :',mb,lmb,' lmax :',lnb+lmb, ' kkbeta :',upf(nt)%kkbeta
              do ir=1,upf(nt)%kkbeta
                 write (4,'(12f16.10)')  rgrid(nt)%r(ir),(qrad_rs(ir,ijv,l+1), l=0,lnb+lmb )
              end do
              CLOSE (4)
!-
           end do
        end do
     end if
     ! ntyp
  enddo
  !
  deallocate (ylmk0)
  deallocate (power_0, power_r, power_q, power_rs, power_qs)
  deallocate (qrad_q, qrad_r, qrad_rs, ffrr)
  deallocate (aux1, aux)

  call stop_clock ('init_us_0')
  return

  contains

  REAL(DP) function filter ( x, a, n )
  implicit none
  REAL (DP), INTENT (IN) :: x, a
  INTEGER,   INTENT (IN) :: n
  REAL (DP) :: axx, ff
  integer :: k

  axx = a * x * x

  ff = 1.d0
  do k = n, 1, -1
     ff = (1.d0+axx/REAL(k)*ff)
  end do
  filter = ff * exp(-axx)

  return
  end function filter

end subroutine init_us_0
