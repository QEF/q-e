!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
subroutine init_us_b0
  !----------------------------------------------------------------------
  !
  ! in this routine the beta_l(r) are smoothed 
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
  implicit none
  !
  !     here a few local variables
  !
  logical, parameter :: tprint=.false.     ! whether the beta_l(r) and its relatives are printed or not
  integer, parameter :: nn=16   ! smoothing parameter, order of the polynomial inverse gaussian approximant
  real(DP), parameter:: a=22.0  ! smoothing parameter, exponent of the gaussian decaying factor
  real (DP) :: rcut, drcut ! beta function cutoff radius and estimated increase of it due to the filtering
  real(DP), parameter :: eps = 1.d-8

  integer :: nqx
  real(DP),allocatable :: tab0(:,:), tab(:,:), beta(:,:), betas(:,:)
  real(DP),allocatable :: power_r(:), power_q(:)
  !
  integer :: nt, nb, l, ir, iq, startq, lastq, ndm
  ! various counters
  real(DP), allocatable :: aux (:), besr (:)
  ! various work space
  real(DP) :: q, qi, qmax
  ! the modulus of g for each shell
  ! q-point grid for interpolation
  real(DP) ::  vqint 
  ! interpolated value
  !
  character(LEN=4) :: filename

  !
  call start_clock ('init_us_b0')
  !
  !    Initialization of variables
  !
  drcut = abs(log(eps8))/sqrt(ecutwfc)/2.d0 
  qmax = 3.d0 * sqrt(ecutwfc)
  nqx = int( qmax / dq + 4)  ! Think about what happens in a variable cell calculations
  
  ALLOCATE (tab0( nqx , nbetam ), tab( nqx, nbetam ) )
  ALLOCATE (power_r( nbetam ), power_q( nbetam ) )

  if (tprint) write ( stdout, * ) 'upf(nt)%kkbeta ', upf(1:ntyp)%kkbeta
  do nt=1,ntyp
     rcut = rgrid(nt)%r(upf(nt)%kkbeta)
     do ir = upf(nt)%kkbeta, upf(nt)%mesh
        if ( rgrid(nt)%r(ir) < rcut + drcut ) upf(nt)%kkbeta=ir
     end do
  end do

  ndm = MAXVAL ( upf(:)%kkbeta )
  ALLOCATE (beta( ndm, nbetam ), betas( ndm, nbetam ) )
  allocate (aux ( ndm ), besr( ndm ))

  write (6,'(a,f6.2,a,i4,4(a,f11.8))') 'FILTER : a=',a,', nn=',nn,', filter(1.1d0)=', filter(1.1d0,a,nn), &
                                                                  ', filter(1.1d0/3)=', filter(1.1d0/3,a,nn), &
                                                                  ', filter(1.2d0)=', filter(1.2d0,a,nn), &
                                                                  ', filter(1.2d0/3)=', filter(1.2d0/3,a,nn)
  if (tprint) then
     write ( stdout, * ) " PSEUDOPOTENTIAL REPORT; initial values "
     write ( stdout, * ) ' NQX :', nqx, ' NBETAM :', nbetam
     write ( stdout, * ) ' NDM :', ndm, ' KKBETA :', upf(1:ntyp)%kkbeta
     write ( stdout, * ) ' r(ndm):', (rgrid(nt)%r(upf(nt)%kkbeta),nt=1,ntyp), ' drcut :', drcut
  end if

  !
  !     fill the interpolation table tab
  !
  call divide (intra_bgrp_comm, nqx, startq, lastq)
  do nt = 1, ntyp
!     write (*,*) ' ityp = ', nt
     if ( upf(nt)%is_gth ) cycle

     tab0(:,:) = 0.d0  ; tab (:,:) = 0.d0
     beta (:,:) = 0.d0 ; betas(:,:) = 0.d0

     if (tprint) then
       filename = 'br_'                 !  the radial beta_l(r) as defined in the pseudopotential
       WRITE (filename( 4:4 ),'(i1)') nt
       OPEN (4,file=filename,form='formatted', status='unknown')
       write (4, *) '# the radial beta_l(r) as defined in the pseudopotential'
       write (4, *) '# nbeta :', upf(nt)%nbeta,' kkbeta :',upf(nt)%kkbeta
       do ir =1,upf(nt)%kkbeta
          write (4,'(12f16.10)') rgrid(nt)%r(ir), (upf(nt)%beta(ir,nb), nb=1,upf(nt)%nbeta)
       end do
       CLOSE (4)
     end if
!-
     do nb = 1, upf(nt)%nbeta
        l = upf(nt)%lll (nb)
        do iq = startq, lastq
           qi = (iq - 1) * dq
           call sph_bes (upf(nt)%kkbeta, rgrid(nt)%r, qi, l, besr)
           do ir = 1, upf(nt)%kkbeta
              aux (ir) = upf(nt)%beta (ir, nb) * besr (ir) * rgrid(nt)%r(ir)
           enddo
           call simpson (upf(nt)%kkbeta, aux, rgrid(nt)%rab, vqint)

           tab0(iq, nb) = vqint
!           tab (iq, nb) = vqint * filter( 1.1d0 * qi / qmax, a, nn )
           tab (iq, nb) = vqint * filter( 1.2d0 * qi / qmax, a, nn )

           do ir = 1, upf(nt)%kkbeta
              beta (ir,nb) = beta (ir,nb) + qi*qi * dq * tab0(iq,nb) * besr (ir) * rgrid(nt)%r(ir)
              betas(ir,nb) = betas(ir,nb) + qi*qi * dq * tab (iq,nb) * besr (ir) * rgrid(nt)%r(ir)
           enddo
        enddo
!-
     enddo

     call mp_sum( tab0, intra_bgrp_comm )
     call mp_sum( tab , intra_bgrp_comm )
     beta = beta * 8.0/fpi ; call mp_sum( beta , intra_bgrp_comm )
     betas= betas* 8.0/fpi ; call mp_sum( betas, intra_bgrp_comm )
     upf(nt)%beta(1:upf(nt)%kkbeta,1:upf(nt)%nbeta) = betas(1:upf(nt)%kkbeta,1:upf(nt)%nbeta) 

!-Fix the core cutoff radius kkbeta so that no more than eps of the integral is lost.
! Parseval's identity is used to monitor the real space completeness of the integral of the(squarred)  betas

! Compute the integral in Fourier space
     do nb = 1, upf(nt)%nbeta
        l = upf(nt)%lll (nb)
        power_q(nb) =0.d0
        do iq = startq, lastq
           qi = (iq - 1) * dq
           power_q(nb) = power_q(nb) + qi*qi * dq * tab(iq,nb) * tab(iq,nb)
        end do
     end do
     power_q  = power_q * 8.0/fpi ; call mp_sum( power_q, intra_bgrp_comm )

! Compute the integral in real space up to the current value of kkbeta and try to reduce it if close enough to power_q
     power_r = power_q                    ! initialize power_r so that the first test is passed. 
     upf(nt)%kkbeta = upf(nt)%kkbeta + 1  ! increase kkbeta because it must be decreased in a moment
     do while ( MAXVAL (1.d0-power_r(1:upf(nt)%nbeta)/power_q(1:upf(nt)%nbeta)) .lt. eps ) 
        write (*,*) upf(nt)%kkbeta,rgrid(nt)%r(upf(nt)%kkbeta),1.d0-power_r(1:upf(nt)%nbeta)/power_q(1:upf(nt)%nbeta)
        upf(nt)%kkbeta = upf(nt)%kkbeta - 1
        do nb = 1, upf(nt)%nbeta
           l = upf(nt)%lll (nb)
           do ir = 1, upf(nt)%kkbeta
              aux (ir) = betas (ir, nb) * betas (ir, nb)
           enddo
           call simpson (upf(nt)%kkbeta, aux, rgrid(nt)%rab, vqint)
           power_r(nb) = vqint
        end do
     end do
!!!

     if (tprint) then
!-
        filename(1:3) = 'bq_'                 ! the radial fourier transform of beta_l in reciprcal space up to qmax
        OPEN (4,file=filename,form='formatted', status='unknown')
        write (4, *) '# the radial fourier transform of beta_l in reciprcal space up to qmax'
        write (4, *) '# nbeta :', upf(nt)%nbeta,' nqx :',nqx, dq*nqx
        do iq=1,nqx
           qi = (iq - 1) * dq
           write (4,'(12f16.10)')  qi,(tab0(iq,nb), nb=1,upf(nt)%nbeta)
        end do
        CLOSE (4)
!-
        filename(1:3) = 'bqs'                 ! the smoothed radial fourier transform of beta_l in reciprcal space up to qmax
        OPEN (4,file=filename,form='formatted', status='unknown')
        write (4, *) '# the smoothed radial fourier transform of beta_l in reciprcal space up to qmax'
        write (4, *) '# nbeta :', upf(nt)%nbeta,' nqx :',nqx, dq*nqx
        do iq=1,nqx
           qi = (iq - 1) * dq
           write (4,'(12f16.10)')  qi,(tab(iq,nb), nb=1,upf(nt)%nbeta)
        end do
        CLOSE (4)
!-
        filename(1:3) = 'brq'                 ! the back radial fourier transform of beta_l in real space
        OPEN (4,file=filename,form='formatted', status='unknown')
        write (4, *) '# the back radial fourier transform of beta_l in real space'
        write (4, *) '# nbeta :', upf(nt)%nbeta,' kkbeta :',upf(nt)%kkbeta
        do ir =1,upf(nt)%kkbeta
           write (4,'(12f16.10)') rgrid(nt)%r(ir), (beta(ir,nb), nb=1,upf(nt)%nbeta)
        end do
        CLOSE (4)
!-
        filename(1:3) = 'brs'                 ! the back radial fourier transform of the smoothed beta_l in real space
        OPEN (4,file=filename,form='formatted', status='unknown')
        write (4, *) '# the back radial fourier transform of the smoothed beta_l in real space'
        write (4, *) '# nbeta :', upf(nt)%nbeta,' kkbeta :',upf(nt)%kkbeta
        do ir =1,upf(nt)%kkbeta
           write (4,'(12f16.10)') rgrid(nt)%r(ir), (betas(ir,nb), nb=1,upf(nt)%nbeta)
        end do
        CLOSE (4)
!-
     end if

  enddo
  if (tprint) then
     filename(1:4) = 'ffqq'                ! the filter in reciprocal space
     OPEN (4,file=filename,form='formatted', status='unknown')
     write (4, *) '# the filter in reciprocal space'
     write (4, *) '# fillter : a=',a,', nn=', nn
     do iq = 1, nqx
        qi = (iq - 1) * dq
        write (4,'(12f16.10)')  qi,filter(1.1d0* qi / qmax, a, nn)
     end do
     CLOSE (4)
  end if
  if (tprint) then
     write ( stdout, * ) " PSEUDOPOTENTIAL REPORT; updated values "
     write ( stdout, * ) ' NDM :', ndm, ' KKBETA :', upf(1:ntyp)%kkbeta
     write ( stdout, * ) ' r(ndm):', (rgrid(nt)%r(upf(nt)%kkbeta),nt=1,ntyp)
  end if
  deallocate (power_r, power_q)
  deallocate (tab0, tab)
  deallocate (beta, betas)
  deallocate (besr, aux)

  call stop_clock ('init_us_b0')
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
  
end subroutine init_us_b0

