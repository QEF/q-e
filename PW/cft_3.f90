!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!     This is a collection of the fft routines which are used on differe
!     machines. The performance of the code depends heavily upon the
!     performance of these routines. Therefore it is convenient always
!     to use machine-specific routines
!
!     If __FFTW is defined, the fftw library is used.
!     Otherwise machine-specific routines are used (if available).
!     If a machine-specific routine is not available either a fake
!     routine issuing an error message is compiled.
!
#include "machine.h"
#ifdef __FFTW
!
!-----------------------------------------------------------------------

subroutine cft_3 (f, nr1, nr2, nr3, nr1x, nr2x, nr3x, igrid, isign)
  !-----------------------------------------------------------------------
  ! driver routine for 3d fft using fftw libraries (PG)
  !
  use parameters, only : DP
  use fftw_mod

  implicit none
  integer :: nr1, nr2, nr3, nr1x, nr2x, nr3x, isign, igrid
  complex(kind=DP) :: f (nr1x * nr2x * nr3x)
  !
  real(kind=DP) :: fac
  integer :: ibid
  !
  ! initialization variables
  !
  C_POINTER, save :: plan(2,2)
  data plan/0,0,0,0/
  !plan = reshape((/ 0, 0, 0, 0 /),(/2,2/))
  !
  !
  if (nr1.ne.nr1x.or.nr2.ne.nr2x) call errore ('cft_3', 'not implemented', 1)
  if (igrid.le.0.or.igrid.gt.2) call errore ('cft_3', 'which grid ?',1)
  if (isign.eq.1) then
     ibid = 1
  elseif (isign.eq. - 1) then
     ibid = 2
  else
     call errore ('cft_3', 'isign unexpected', isign)

  endif

  if (plan (igrid, ibid) .eq.0) call FFTW3D_F77_CREATE_PLAN (plan (igrid, ibid),&
       nr1, nr2, nr3, isign, FFTW_ESTIMATE+FFTW_IN_PLACE)
  call FFTWND_F77_ONE (plan (igrid, ibid), f, 0)
  if (isign.eq. - 1) then
     fac = 1.0 / float (nr1 * nr2 * nr3)
     call DSCAL (2 * nr1 * nr2 * nr3, fac, f, 1)
  endif
  return

end subroutine cft_3
#else
#undef PRESENT
#if defined(__SGI) || defined(__ORIGIN)
#define PRESENT
!----------------------------------------------------------------------

subroutine cft_3 (f, n1, n2, n3, nm1, nm2, nm3, igrid, sign)
  !     ===============
  !     silicon graphics  driver routine for 3d complex fft (complib)
  !
  !----------------------------------------------------------------------

  use parameters, only : DP
  implicit none
  integer :: n1, n2, n3, nm1, nm2, nm3, igrid, sign

  complex(kind=DP) :: f (nm1, nm2, nm3)
  integer :: ngrid, nmax
  parameter (ngrid = 2, nmax = 1000)
  logical, save :: first (ngrid)
  complex(kind=DP), save :: aux (nmax, ngrid)
  real(kind=DP) :: fac
  data first / ngrid * .true. /

  !save first, aux
  !first = (/(.true.,i=1,ngrid)/)
  if (n1 + n2 + n3 + 45.gt.nmax) call errore ('cft_3', 'nmax too small',&
       n1 + n2 + n3 + 45)
  if (first (igrid) ) then
     call zfft3di (n1, n2, n3, aux (1, igrid) )
     first (igrid) = .false.
  endif

  call zfft3d (sign, n1, n2, n3, f, nm1, nm2, aux (1, igrid) )
  if (sign.lt.0) then
     fac = 1.0d0 / dfloat (n1 * n2 * n3)
     call DSCAL (2 * nm1 * nm2 * n3, fac, f, 1)

  endif
  return

end subroutine cft_3
#endif
#if defined(EXEMPLAR)
#define PRESENT
!----------------------------------------------------------------------

subroutine cft_3 (f, n1, n2, n3, nm1, nm2, nm3, igrid, sign)
  !     ===============
  !     exemplar graphics  driver routine for 3d complex fft (veclib)
  !
  !----------------------------------------------------------------------

  use parameters, only : DP
  implicit none
  integer :: n1, n2, n3, nm1, nm2, nm3, igrid, sign

  complex(kind=DP) :: f (nm1, nm2, nm3)

  integer :: ier
  call z3dfft (f, n1, n2, n3, nm1, nm2, sign, ier)

  if (ier.ne.0) call errore ('cft_3', ' fft error ', abs (ier) )
  return

end subroutine cft_3
#endif
#ifdef CRAYY
#define PRESENT
!----------------------------------------------------------------------
subroutine cft_3 (ac, n1, n2, n3, nm1, nm2, nm3, igrid, isign)
  !----------------------------------------------------------------------
  !
  !      3d fft - cray scilib version
  !
  use parameters, only : DP
  implicit none
  integer :: n1, n2, n3, nm1, nm2, nm3, igrid, isign

  real :: ac (2, nm1, nm2, nm3)
  integer :: ngrid, nmax
  parameter (ngrid = 2, nmax = 256)
  integer, save :: ifax (19, 3, ngrid), np1
  integer       :: inc, lot, jump, j, k
  real, save    :: trig (2 * nmax, 3, ngrid)
  real ::  work (4 * nm1 * n2 * n3), fac
  logical, save :: first (ngrid)
  data first / ngrid * .true. /
  !save trig, ifax, first, np1
  external sscal
  !
  !
  !first = (/(.true.,i=1,ngrid)/)

  if (igrid.le.0.or.igrid.gt.ngrid) call errore ('cft_3', 'which grid?', 1)
  if (sign.ne. - 1.and.sign.ne.1) call errore ('cft_3', 'which fft ?', 2)

  if (n1.gt.nmax.or.n2.gt.nmax.or.n3.gt.nmax) call errore ('cft_3', &
       'increase nmax', 3)
  !
  if (first (igrid) ) then
     if (mod (n1, 2) .eq.0) then
        np1 = n1 + 1
     else
        np1 = n1
     endif
     if (np1.gt.nm1) call errore ('cft3', 'too large input dimension', np1)
     if (n2.gt.nm2) call errore ('cft3', 'too large input dimension', &
          n2)
     if (n3.gt.nm3) call errore ('cft3', 'too large input dimension', &
          n3)
     !
     call cftfax (n1, ifax (1, 1, igrid), trig (1, 1, igrid) )
     call cftfax (n2, ifax (1, 2, igrid), trig (1, 2, igrid) )
     call cftfax (n3, ifax (1, 3, igrid), trig (1, 3, igrid) )
     first (igrid) = .false.
  endif
  !

  if (np1.ne.nm1.or.n2.ne.nm2) call errore ('cft_3', 'no longer implemented', 1)
  !     & call mcpack(ac,nm1,nm2,nm3,ac,np1,n2,n3,1)
  if (n1.ne.np1) then
     do k = 1, n3
        do j = 1, n2
           ac (1, np1, j, k) = 0.0
           ac (2, np1, j, k) = 0.0
        enddo
     enddo
  endif
  !
  !     ... i-direction
  !
  inc = 2
  jump = 2 * np1
  lot = n2 * n3
  call cfftmlt (ac (1, 1, 1, 1), ac (2, 1, 1, 1), work, trig (1, 1, &
       igrid), ifax (1, 1, igrid), inc, jump, n1, lot, isign)
  !
  !     ... j-direction
  !
  inc = 2 * np1
  jump = 2
  lot = n1
  do k = 1, n3
     call cfftmlt (ac (1, 1, 1, k), ac (2, 1, 1, k), work, trig (1, 2, &
          igrid), ifax (1, 2, igrid), inc, jump, n2, lot, isign)
  enddo
  !
  !     ... k-direction
  !
  inc = 2 * np1 * n2
  jump = 2
  lot = np1 * n2

  call cfftmlt (ac (1, 1, 1, 1), ac (2, 1, 1, 1), work, trig (1, 3, &
       igrid), ifax (1, 3, igrid), inc, jump, n3, lot, isign)
  !      if( np1.ne.nm1 .or. n2.ne.nm2 )
  !     & call mcpack(ac,nm1,nm2,nm3,ac,np1,n2,n3,-1)
  !
  if (isign.lt.0) then
     fac = 1.0 / (n1 * n2 * n3)
     call sscal (2 * nm1 * nm2 * nm3, fac, ac, 1)
  endif
  !
  return

end subroutine cft_3
#endif

#ifdef __SX4
#define PRESENT
!----------------------------------------------------------------------
subroutine cft_3 (f, nr1, nr2, nr3, nrx1, nrx2, nrx3, igrid, sign)
  !----------------------------------------------------------------------
  !
  !      3d fft for NEC:
  !      uses GPFA routines
  !

  use parameters, only : DP
  implicit none

  integer :: nr1, nr2, nr3, nrx1, nrx2, nrx3, igrid, sign
  !
  ! input: the logical dimension of the FFT
  !
  !
  ! input: the physical dimension of the FFT
  !
  ! input: grid used (1=thick, 2=smooth)
  ! input: the sign of the transformation

  real(kind=DP) :: f (2, nrx1, nrx2, nrx3)
  ! inp/out: the function to transform
  integer :: ngrid, nmax
  ! max number of different grid allowed
  ! max value of n1, n2, n3 allowed

  parameter (ngrid = 2, nmax = 256)

  integer :: k, inc, jump, lot
  ! counter on z direction
  ! the increment between different values
  ! the jump between fft arrays
  ! how many fft

  real(kind=DP), save :: trig1 (2 * nmax, ngrid), trig2 (2 * nmax, ngrid), &
       trig3 (2 * nmax, ngrid), fact
  !
  !    trigonometric factors
  !
  !    the multiplication factor
  logical, save :: first (ngrid)
  ! is true at the first iteration

  data first / ngrid * .true. /
  !save first, trig1, trig2, trig3
  !
  !    test the sign and put the correct normalization on f
  !

  !first = (/(.true.,i=1,ngrid)/)

  if (sign.eq. - 1) then
     fact = 1.d0 / float (nr1 * nr2 * nr3)
     call sscal (2 * nrx1 * nrx2 * nrx3, fact, f, 1)
  elseif (sign.ne.1) then
     call errore ('cft_3', 'wrong isign', 1)
  endif
  if (igrid.le.0.or.igrid.gt.ngrid) call errore ('cft_3', 'which grid?', 1)

  if (nr1.gt.nmax.or.nr2.gt.nmax.or.nr3.gt.nmax) call errore ( &
       'cft_3', 'increase nmax', 3)
  !
  !   At the first iteration initialize
  !
  if (first (igrid) ) then
     call setgpfa (trig1 (1, igrid), nr1)
     call setgpfa (trig2 (1, igrid), nr2)
     call setgpfa (trig3 (1, igrid), nr3)
     first (igrid) = .false.
  endif
  !
  !     i-direction
  !
  inc = 2
  jump = 2 * nrx1
  lot = nr2 * nr3
  call gpfa (f (1, 1, 1, 1), f (2, 1, 1, 1), trig1 (1, igrid), &
       inc, jump, nr1, lot, sign)
  !
  !     ... j-direction ...
  !
  inc = 2 * nrx1
  jump = 2
  lot = nr1
  do k = 1, nr3
     call gpfa (f (1, 1, 1, k), f (2, 1, 1, k), trig2 (1, igrid), &
          inc, jump, nr2, lot, sign)
  enddo
  !
  !     ... k-direction
  !
  inc = 2 * nrx1 * nrx2
  jump = 2
  lot = nrx1 * nr2

  call gpfa (f (1, 1, 1, 1), f (2, 1, 1, 1), trig3 (1, igrid), &
       inc, jump, nr3, lot, sign)
  return

end subroutine cft_3
#endif

#ifdef __SX6
#define PRESENT
#define ASL

module afftnec

  use parameters
  implicit none

  integer, parameter :: ngrid=2
  integer, parameter :: dim_iw=60
  integer:: nrz1(ngrid),nrz2(ngrid),nrz3(ngrid)
  logical :: first(ngrid)        ! is true at the first iteration
  data first/ngrid*.true./
  real(kind=DP), dimension(ngrid) :: fact
  real(kind=DP), allocatable, target, dimension(:,:) :: auxp
  integer, target, dimension(dim_iw,ngrid) :: iw0

end module afftnec
!
!----------------------------------------------------------------------
subroutine cft_3(f,nr1,nr2,nr3,nrx1,nrx2,nrx3,igrid,sign)
  !----------------------------------------------------------------------
  !
  !      3d fft for NEC SX6 - uses ASL library routines
  !      contributed by Guido Roma
  !
  use parameters, only : DP
  use afftnec
  implicit none

  integer :: &
       &       nr1,&       !
       &       nr2,&       ! input: the logical dimension of the FFT
       &       nr3,&       !
       &       nrx1,&      !
       &       nrx2,&      ! input: the physical dimension of the FFT
       &       nrx3,&      !
       &       igrid,&     ! input: grid used (1=thick, 2=smooth)
       &       sign,&      ! input: the sign of the transformation
       &       ierr,&      ! 
       isw 
  complex(kind=DP) :: &       
       &       f(nrx1,nrx2,nrx3)    ! inp/out: the function to transform
  complex(kind=DP), allocatable, dimension(:,:,:) :: f1 ! for ASL Library FFT routines 
#ifdef ASL
  integer, pointer, dimension(:) :: iw
#if defined MICRO
  common/NEC_ASL_PARA/nbtasks
  integer :: nbtasks
#endif
#endif
  real(kind=DP), pointer, dimension(:) :: cw1
  complex(kind=DP), dimension(:), allocatable :: cw2   

  !     allocate auxp at the first call (independently of the grid)
  if (.not.allocated(auxp)) allocate(auxp(2*(nr1+nr2+nr3),ngrid))

  !
  !    test the sign and put the correct normalization on f
  !
  if (first(igrid)) then
     nrz1(igrid)=nrx1
     nrz2(igrid)=nrx2
     nrz3(igrid)=nrx3
     if (mod(nrx1,2)==0) nrz1(igrid)=nrx1+1
     if (mod(nrx2,2)==0) nrz2(igrid)=nrx2+1
     if (mod(nrx3,2)==0) nrz3(igrid)=nrx3+1
  end if
#ifdef ASL
  allocate(cw2(nrz1(igrid)*nrz2(igrid)*nrz3(igrid)))
#else
  allocate(cw2(3*nrz1(igrid)*nrz2(igrid)*nrz3(igrid)))
#endif
  allocate(f1(nrz1(igrid),nrz2(igrid),nrz3(igrid)))

  if ( sign.eq.-1 ) then
     fact(igrid)=1.0_8/dble(nr1*nr2*nr3)
     call DSCAL(2*nrx1*nrx2*nrx3,fact(igrid),f,1)
  else if ( sign.ne.1) then
     call errore('cft_3', 'wrong isign',1)
  endif
  if (igrid.le.0.or.igrid.gt.ngrid)&
       &  call errore('cft_3','which grid ?',1)

  !     copy f in the auxiliary f1 with odd dimensions
  !      call ZCOPY(nrx1*nrx2*nrx3,f,1,f1(1:nrx1,1:nrx2,1:nrx3),1)
  f1(1:nrx1,1:nrx2,1:nrx3)=f

#ifdef ASL
  call zfc3cl(f1,nr1,nr2,nr3,nrz1(igrid),nrz2(igrid),nrz3(igrid),ierr)
  call errore('cft_3', 'initialisation problem',ierr)
  iw=>iw0(:,igrid)
#endif
  cw1=>auxp(:,igrid)

  if (first(igrid)) then
     isw=0
     first(igrid)=.false.
     !         write(6,*)'________________________________________________________'
     !         write(6,*) 'igrid = ',igrid
     !         write(6,*) '  nrxs => ',nrx1,nrx2,nrx3
     !         write(6,*) '  nrzs => ',nrz1(igrid),nrz2(igrid),nrz3(igrid)
     !         write(6,*) '  nrs => ',nr1,nr2,nr3
     !      write(6,*)'size(auxp)',size(auxp,1),size(auxp,2)
     !      write(6,*)'size(cw1)',size(cw1)
     !      write(6,*)'size(iw)',size(iw)
     !         write(6,*)'________________________________________________________'
#ifdef ASL
#if defined MICRO
     call hfc3fb(nr1,nr2,nr3,f1,nrz1(igrid),nrz2(igrid),nrz3(igrid),&
          &            isw,iw,cw1,cw2,nbtasks,ierr)
#else
     call zfc3fb(nr1,nr2,nr3,f1,nrz1(igrid),nrz2(igrid),nrz3(igrid),&
          &            isw,iw,cw1,cw2,ierr)
#endif
     if (ierr.ne.0) call errore('cft_3','ierr=',ierr)
#else
     call ZZFFT3D(0,nr1,nr2,nr3,1.d0,f1,nrz1(igrid),nrz2(igrid),&
          &             f1,nrz1(igrid),nrz2(igrid),cw1,cw2,ierr)
#endif
  endif

  isw=sign
#ifdef ASL
#if defined MICRO
  call hfc3bf(nr1,nr2,nr3,f1,nrz1(igrid),nrz2(igrid),nrz3(igrid),&
       &            isw,iw,cw1,cw2,nbtasks,ierr)
#else
  call zfc3bf(nr1,nr2,nr3,f1,nrz1(igrid),nrz2(igrid),nrz3(igrid),&
       &            isw,iw,cw1,cw2,ierr)     
#endif
  if (ierr.ne.0) call errore('cft_3','ierr=',ierr)
#else
  call ZZFFT3D(isw,nr1,nr2,nr3,1.d0,f1,nrz1(igrid),nrz2(igrid),&
       &             f1,nrz1(igrid),nrz2(igrid),cw1,cw2,ierr)
#endif


  !     copy f1 back in f with odd dimensions
  !      call zcopy(nrx1*nrx2*nrx3,f1(1:nrx1,1:nrx2,1:nrx3),1,f,1)
  f(:,:,:)=f1(1:nrx1,1:nrx2,1:nrx3)
  deallocate(f1)
  deallocate(cw2)
  nullify(cw1)
#ifdef ASL
  nullify(iw)
#endif
  !
  return
#endif

#ifdef __AIX
#define PRESENT
!----------------------------------------------------------------------


subroutine cft_3 (f, n1, n2, n3, nx1, nx2, nx3, igrid, sign)
  !     ===============
  !
  !     ibm driver routine for 3d complex fft (essl library)
  !     nx1=n1+1 is allowed (in order to avoid memory conflicts)
  !     for compatibility: nx2=n2, nx3=n3. nx2 and nx3 are not used
  !----------------------------------------------------------------------
  use parameters, only : DP
  implicit none
  integer :: n1, n2, n3, nx1, nx2, nx3, sign, igrid

  complex(kind=DP) :: f (nx1 * nx2 * nx3)
  integer :: isign, naux
  parameter (naux = 60000)

  real(kind=DP) :: aux (naux), scale
  if (sign.ne. - 1.and.sign.ne.1) call errore ('cft_3', 'which fft ?' &
       &, 1)
  !
  !  ESSL sign convention for fft's is the opposite of the "usual" one
  !

  isign = - sign
  if (isign.gt.0) then
     scale = 1.0d0 / (n1 * n2 * n3)
  else
     scale = 1.0d0

  endif

  call dcft3 (f, nx1, nx1 * nx2, f, nx1, nx1 * nx2, n1, n2, n3, &
       isign, scale, aux, naux)
  return



end subroutine cft_3
#endif
#ifdef SUN
#define PRESENT
!----------------------------------------------------------------------


subroutine cft_3 (f, n1, n2, n3, nx1, nx2, nx3, igrid, sign)
  !     ===============
  !
  ! SUNperf library
  !
  !----------------------------------------------------------------------
  use parameters, only : DP

  implicit none

  integer :: n1, n2, n3, nx1, nx2, nx3, sign, igrid
  complex(kind=DP) :: f ( nx1 , nx2 , nx3 )

  external zfft3i, zfft3f, zfft3b, zdscal
!
! Local variables
!
  integer :: lwork
  complex(kind=DP) ::  work ( 4*(nx1 + nx2 + nx3) + 45 )
  real(kind=DP) :: scale

  lwork = 4 * ( nx1 + nx2 + nx3 ) + 45

  if (sign.ne. - 1.and.sign.ne.1) call errore ('cft_3', 'which fft ?', 1)

  call zfft3i ( nx1, nx2, nx3, work )

  if( sign.eq. 1 ) then
   call zfft3b ( n1, n2, n3, f, nx1, nx2, work, lwork )
  else
   call zfft3f ( n1, n2, n3, f, nx1, nx2, work, lwork )
   scale = 1.0d0 /dble( n1 * n2 * n3 )
   call zdscal ( nx1*nx2*nx3, scale, f, 1 )
  endif

  return

end subroutine cft_3
#endif
#ifdef DEC
#ifdef DXML
#define PRESENT
!----------------------------------------------------------------------

subroutine cft_3 (f, n1, n2, n3, nm1, nm2, nm3, igrid, sign)
  !     ===============
  !     driver routine for 3d complex fft using DXML/CXML libraries
  !
  !----------------------------------------------------------------------

  use parameters, only : DP
  implicit none
  integer :: n1, n2, n3, nm1, nm2, nm3, igrid, sign

  complex(kind=DP) :: f (nm1, nm2, nm3)
  STRUCTURE / DXML_Z_FFT_STRUCTURE /
  INTEGER :: N
  LOGICAL :: STRIDE_1_FLAG
  INTEGER :: N_TI (0:16)
  INTEGER :: N_K (0:16)
  INTEGER :: N_T (0:16)
  INTEGER :: TYPE (0:16)
  INTEGER :: NUM_STAGES
  INTEGER (8) :: ROTATION_VECTOR
  INTEGER :: ROTATION_VECTOR_SIZE
  INTEGER (8) :: TEMP_AREA
  INTEGER :: TEMP_AREA_SIZE
  INTEGER :: SET_BLOCK_SIZE
  INTEGER :: NUM_NON_SPECIAL_RADIX
  INTEGER :: NON_SPECIAL_RADIX (0:16)
  INTEGER :: NON_SPEC_RAD_TWIDDLE_SIZE
  INTEGER (8) :: NON_SPEC_RAD_TWIDDLE_VEC
  INTEGER (8) :: NON_SPECIAL_RADIX_COS (0:16)
  INTEGER (8) :: NON_SPECIAL_RADIX_SIN (0:16)
  INTEGER :: FUTURE_USE (20)
  INTEGER :: GK (0:11)
  ENDSTRUCTURE
  integer :: ngrid
  parameter (ngrid = 2)
  record / DXML_Z_FFT_STRUCTURE / fft_struct (ngrid)
  integer :: status, zfft_init_3d, zfft_exit_3d
  real(kind=DP) :: norm
  character (len=1) :: direction
  logical :: first (ngrid)
  data first / ngrid * .true. /
  save first, fft_struct

  norm = dfloat (n1 * n2 * n3)
  if (sign.eq.1) then
     call dscal (2 * nm1 * nm2 * nm3, norm, f, 1)
     direction = 'b'
  elseif (sign.eq. - 1) then
     call dscal (2 * nm1 * nm2 * nm3, 1.0d0 / norm, f, 1)
     direction = 'f'
  else
     call errore ('cft_3', 'sign unexpected',1)
  endif

  if (first (igrid) ) then
     status = zfft_exit_3d (fft_struct (igrid) )
     ! not sure whether the above call is useful or not
     status = zfft_init_3d (n1, n2, n3, fft_struct (igrid), .true.)
     first (igrid) = .false.
  endif

  call zfft_apply_3d ('C', 'C', direction, f, f, nm1, nm2, &
       fft_struct (igrid) , 1, 1, 1)
  return

end subroutine cft_3
#endif
#endif

#ifdef FUJ64
#define PRESENT
!
!---------------------------------------------------------------
      subroutine cft_3(ac,n1,n2,n3,nm1,nm2,nm3,igrid,isign)
!---------------------------------------------------------------
!
!      3d fft - FUJITSU, FFTVPLIB version (GRoma, 2001)
!
  use parameters, only: DP
  implicit none
  integer :: n1,n2,n3,nm1,nm2,nm3,igrid,isign
  integer, parameter :: ngrid=2, nmax=256
  real(kind=DP) ::  ac(2,nm1,nm2,nm3)
  real(kind=DP) :: workarray(2*nm1*nm2*nm3)
  real(kind=DP) :: trig(2*3*nmax+2*120,ngrid),norm(ngrid)

  integer :: iopt,ierr,trigdim(ngrid)
  integer :: idim(3,ngrid)
  character(len=2) :: init
  character(len=1), dimension(-1:1),      &
  parameter :: mode=(/'m','x','p'/),scale=(/'n','s','i'/)
  logical first(ngrid)
  data first/ngrid*.true./
  save first,trig,idim,trigdim,norm
!
!
  iopt=sign(1,isign)
  if(nm1/=n1)  &
    call errore('cft_3','not any more implemented',nm1)
  if (igrid.le.0.or.igrid.gt.ngrid)  &
    call errore('cft_3','which grid ?',1)
  if (isign.ne.-1 .and. isign.ne.1)  &
    call errore('cft_3','which fft ?', 2)
  if (n1.gt.nmax .or. n2.gt.nmax .or. n3.gt.nmax)  &
    call errore('cft_3','increase nmax',3)
!
  init(2:2)=scale(iopt)
  if(first(igrid)) then
    init(1:1)='i'
    idim(1,igrid)=nm1
    idim(2,igrid)=nm2
    idim(3,igrid)=nm3
    norm(igrid)=sqrt(dfloat(nm1*nm2*nm3))
    trigdim(igrid)= &
      2*(idim(1,igrid)+idim(2,igrid)+idim(3,igrid))+120
    if(n1.gt.nm1) call errore  &
      ('cft3','too large input dimension',n1)
    if( n2.gt.nm2) call errore &
      ('cft3','too large input dimension',n2)
    if( n3.gt.nm3) call errore &
      ('cft3','too large input dimension',n3)
!
!   The FFTVP library stores idim in a common (a single one)
!   so every time you change grid you have to reinitialise!!!
!   That's why the following line.
!   first = .true.
!   which in fact is not needed if one uses the modified version, fftvplib2

    first(igrid) = .false.
  else
    init(1:1)='r'
  end if
  call dftcbm(ac(1,:,:,:),ac(2,:,:,:),3,idim(:,igrid),workarray, &
    trig(:trigdim(igrid),igrid),mode(iopt),init,ierr)
  call errore('cft_3','problems in fft',ierr)
  if (iopt>0) then
    call DSCAL(2*nm1*nm2*nm3,norm(igrid),ac,1)
  else if(iopt<0) then
    call DSCAL(2*nm1*nm2*nm3,1.0_8/norm(igrid),ac,1)
  end if
!
  return
end subroutine cft_3
#endif

#ifndef PRESENT
subroutine cft_3 (f, n1, n2, n3, nm1, nm2, nm3, igrid, sign)
  use parameters, only : DP
  call errore ('cft_3', 'machine-specific routine not available', 1)
  return

end subroutine cft_3
#endif
#endif


