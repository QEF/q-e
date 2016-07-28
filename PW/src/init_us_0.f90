
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
  USE constants,    ONLY : fpi, sqrt2
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
  logical, parameter :: tprint=.false.    ! whether the q_l(r) and its relatives are printed or not
  integer, parameter :: nn=16   ! smoothing parameter, order of the polynomial inverse gaussian approximant
  real(DP), parameter:: a=22.0  ! smoothing parameter, exponent of the gaussian decaying factor
                                ! a=0.d0 ; nn=0 would be no smoothing.
  !
  integer :: nt, ih, jh, nb, mb, ijv, l, m, ir, iq, is, startq, lastq, ilast, ndm, ia
  ! various counters
  real(DP), allocatable :: aux (:), aux1 (:), qrad0(:,:,:), qrad1(:,:,:), qrad2(:,:,:), ffrr(:)
  ! various work space
  real(DP) :: prefr, q, qi
  ! the prefactor of the q functions
  ! the modulus of g for each shell
  ! q-point grid for interpolation
  real(DP), allocatable :: ylmk0 (:)
  ! the spherical harmonics
  integer, external :: sph_ind
  integer :: lnb, lmb
  real(DP) ::  qmax
  !
  character(LEN=5) :: filename
  !
  call start_clock ('init_us_0')
  !
  !    Initialization of the variables
  !
  ndm = MAXVAL ( upf(:)%kkbeta )
  if (tprint) then
     write ( stdout, * ) " PSEUDOPOTENTIAL REPORT "
     write ( stdout, * ) ' NDM :', ndm, '   ', upf(1:ntyp)%kkbeta
     write ( stdout, * ) ' LMAXQ :', lmaxq, ' NBETAM :', nbetam
  end if

  allocate (aux ( ndm), aux1( ndm))    
  allocate (qrad0(nqxq, nbetam*(nbetam+1)/2,lmaxq ) )
  allocate (qrad1( ndm, nbetam*(nbetam+1)/2,lmaxq ) )
  allocate (qrad2( ndm, nbetam*(nbetam+1)/2,lmaxq ) )
  allocate (ffrr( ndm))    
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

  prefr = fpi / omega
  !
  !   here for the US types we compute the Fourier transform of the
  !   Q functions.
  !   
  call divide (intra_bgrp_comm, nqxq, startq, lastq)
  qmax = sqrt(ecutrho)
  write (6, *) ' qmax : sqrt(ecutrho) =',sqrt(ecutrho), dq*nqxq*tpiba 
  write (6,'(a,f6.2,a,i4,a,f16.10)') 'FILTER : a=',a,', nn=',nn,', filter(1.d0)=', filter(1.d0,a,nn)
  !
  do nt = 1, ntyp
     !
!-
     if (tprint) then
        filename = 'qr'                 !  the radial q_l(r) as defined in the pseudopotential
        do nb=1, upf(nt)%nbeta 
           do mb=nb, upf(nt)%nbeta
              ijv = mb * (mb - 1) / 2 + nb
              lnb = upf(nt)%lll(nb) ; lmb = upf(nt)%lll(mb)
              WRITE (filename( 3:3 ),'(i1)') nb ; WRITE (filename( 4:4 ),'(i1)') mb ; WRITE (filename( 5:5 ),'(i1)') nt
              OPEN (4,file=filename,form='formatted', status='unknown')
              write (4, *) ' nb :', nb,lnb,' mb :',mb,lmb,' lmax :',lnb+lmb, ' kkbeta :',upf(nt)%kkbeta
              do ir =1,upf(nt)%kkbeta
                 write (4,'(12f16.10)') rgrid(nt)%r(ir), (upf(nt)%qfuncl(ir,ijv,l), l=0,lnb+lmb)
              end do
              CLOSE (4)
           end do
        end do
     end if
!-

     if (lmaxq > 0) qrad0(:,:,:)= 0.d0 ; qrad1(:,:,:)= 0.d0 ; qrad2(:,:,:)= 0.d0

     if ( upf(nt)%tvanp ) then
        do l = 0, upf(nt)%nqlc -1
           !
           !     here we compute the spherical bessel function for each |g|
           !
           do iq = startq, lastq
              q = (iq - 1) * dq * tpiba
              call sph_bes ( upf(nt)%kkbeta, rgrid(nt)%r, q, l, aux)
              !
              !   and then we integrate with all the Q functions
              !
              do nb = 1, upf(nt)%nbeta
                 !
                 !    the Q are symmetric with respect to indices
                 !
                 do mb = nb, upf(nt)%nbeta
                    ijv = mb * (mb - 1) / 2 + nb
                    if ( ( l >= abs(upf(nt)%lll(nb) - upf(nt)%lll(mb)) ) .and. &
                         ( l <=     upf(nt)%lll(nb) + upf(nt)%lll(mb)  ) .and. &
                         (mod(l+upf(nt)%lll(nb)+upf(nt)%lll(mb),2)==0) ) then
                       do ir = 1, upf(nt)%kkbeta
                          aux1 (ir) = aux (ir) * upf(nt)%qfuncl(ir,ijv,l)
                       enddo
                       ! the fourier transform of q_nb_mb_l(r) up to qmax
                       call simpson ( upf(nt)%kkbeta, aux1, rgrid(nt)%rab, qrad0(iq,ijv,l+1) )
                    endif
                 enddo
              enddo
              ! igl
              do ir =1, upf(nt)%kkbeta
                 ! q_nb_mb_l(r) from the back fourier transform up to qmax of q_nb_mb_l(q)
                 qrad1(ir,1:nbetam*(nbetam+1)/2,l+1) = &
                 qrad1(ir,1:nbetam*(nbetam+1)/2,l+1) + q*q * dq*tpiba * aux(ir) * &
                    rgrid(nt)%r(ir)**2 *  qrad0( iq, 1:nbetam*(nbetam+1)/2,l+1) 

                 ! q_nb_mb_l(r) from the back fourier transform up to qmax of the smoothed q_nb_mb_l(q)
                 qrad2(ir,1:nbetam*(nbetam+1)/2,l+1) = &
                 qrad2(ir,1:nbetam*(nbetam+1)/2,l+1) + q*q * dq*tpiba * aux(ir) * &
                    rgrid(nt)%r(ir)**2 *  qrad0( iq, 1:nbetam*(nbetam+1)/2,l+1) * filter(q/qmax,a,nn)
              enddo
              if (l==0) then
                 do ir =1, upf(nt)%kkbeta
                    ! the filter function in real space from the back fourier transform up to qmax 
                    ffrr(ir) = ffrr(ir) + q*q * dq*tpiba * aux(ir) * rgrid(nt)%r(ir)**2 * filter(q/qmax,a,nn)
                 enddo
              end if
           enddo
           ! l
        enddo
        qrad0(:,:,:) = qrad0(:,:,:)*prefr
        call mp_sum ( qrad0 , intra_bgrp_comm )
        qrad1(:,:,:) = qrad1(:,:,:) * 8.d0/fpi
        call mp_sum ( qrad1 , intra_bgrp_comm )
        qrad2(:,:,:) = qrad2(:,:,:) * 8.d0/fpi
        call mp_sum ( qrad2 , intra_bgrp_comm )
        do l = 0, upf(nt)%nqlc -1
           do nb = 1, upf(nt)%nbeta
              do mb = nb, upf(nt)%nbeta
               if ( ( l >= abs(upf(nt)%lll(nb) - upf(nt)%lll(mb)) ) .and. &
                    ( l <=     upf(nt)%lll(nb) + upf(nt)%lll(mb)  ) .and. &
                    (mod (l+upf(nt)%lll(nb)+upf(nt)%lll(mb), 2) == 0) ) then
                 ijv = mb * (mb-1) / 2 + nb
                 upf(nt)%qfuncl(1:upf(nt)%kkbeta,ijv,l) = qrad2(1:upf(nt)%kkbeta,ijv,l+1) 
               endif 
              enddo ! mb
           enddo ! nb
        enddo
     endif

     if (tprint) then
!-
        filename = 'ffff'
        OPEN (4,file=filename,form='formatted', status='unknown')
        write (4, *) ' nqxq :', nqxq,' dq :',dq, ' qmax :',qmax
        do iq=1,nqxq
           q = (iq-1)*dq*tpiba
           write (4,'(2f16.10)')  q,filter(q/qmax,a,nn)
        end do
        CLOSE (4)
!-
        filename = 'ffrr'
        OPEN (4,file=filename,form='formatted', status='unknown')
        write (4, *) ' kkbeta :', upf(nt)%kkbeta
        do ir=1,upf(nt)%kkbeta
           write (4,'(2f16.10)')  rgrid(nt)%r(ir),ffrr(ir)
        end do
        CLOSE (4)
!-
        do nb=1, upf(nt)%nbeta 
           do mb=nb, upf(nt)%nbeta
              ijv = mb * (mb - 1) / 2 + nb
              lnb = upf(nt)%lll(nb) ; lmb = upf(nt)%lll(mb)
              WRITE (filename( 3:3 ),'(i1)') nb ; WRITE (filename( 4:4 ),'(i1)') mb ; WRITE (filename( 5:5 ),'(i1)') nt
!-
              filename(1:2) = 'qq'                 ! the radial fourier transform of q_l in reciprcal space
              OPEN (4,file=filename,form='formatted', status='unknown')
              write (4, *) ' nb :', nb,lnb,' mb :',mb,lmb,' lmax :',lnb+lmb, ' nqxq :',nqxq
              do iq=1,nqxq
                 q = (iq-1)*dq*tpiba
                 write (4,'(12f16.10)')  q,(qrad0(iq,ijv,l+1), l=0,lnb+lmb )
              end do
              CLOSE (4)
!-
              filename(1:2) = 'qs'                 ! the smoothed radial fourier transform of q_l in reciprcal space
              OPEN (4,file=filename,form='formatted', status='unknown')
              write (4, *) ' nb :', nb,lnb,' mb :',mb,lmb,' lmax :',lnb+lmb, ' nqxq :',nqxq
              do iq=1,nqxq
                 q = (iq-1)*dq*tpiba
                 write (4,'(12f16.10)')  q,(qrad0(iq,ijv,l+1)*filter(q/qmax,a,nn), l=0,lnb+lmb )
              end do
              CLOSE (4)
!-
              filename(1:2) = 'rq'                 ! the radial q_l(r) as obtained back-transforming the q_l(q)
              OPEN (4,file=filename,form='formatted', status='unknown')
              write (4, *) ' nb :', nb,lnb,' mb :',mb,lmb,' lmax :',lnb+lmb, ' kkbeta :',upf(nt)%kkbeta
              do ir=1,upf(nt)%kkbeta
                 write (4,'(12f16.10)')  rgrid(nt)%r(ir),(qrad1(ir,ijv,l+1), l=0,lnb+lmb )
              end do
              CLOSE (4)
!-
              filename(1:2) = 'rs'                 ! the radial q_l(r) as obtained back-transforming the smoothed q_l(q)
              OPEN (4,file=filename,form='formatted', status='unknown')
              write (4, *) ' nb :', nb,lnb,' mb :',mb,lmb,' lmax :',lnb+lmb, ' kkbeta :',upf(nt)%kkbeta
              do ir=1,upf(nt)%kkbeta
                 write (4,'(12f16.10)')  rgrid(nt)%r(ir),(qrad2(ir,ijv,l+1), l=0,lnb+lmb )
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
  deallocate (qrad0, qrad1, qrad2, ffrr)
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
