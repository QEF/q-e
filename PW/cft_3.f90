!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! ... This is a collection of serial fft drivers for several 
! ... machine-specific libraries (including some obsolete ones). 
! ... The performance of the code may depend heavily upon the
! ... performance of these routines.
!
! ... If a machine-specific routine is not available, a fake
! ... routine issuing an error message is compiled.
!
#include "f_defs.h"
!
#if defined (__SX6)
!
! ... SX6 case
!
MODULE afftnec
  !
  USE kinds
  !
  IMPLICIT NONE
  !
  INTEGER, PARAMETER :: ngrid = 2
  INTEGER, PARAMETER :: dim_iw = 60
  INTEGER            :: nrz1(ngrid), nrz2(ngrid), nrz3(ngrid)
  LOGICAL            :: first(ngrid)        ! is true at the first iteration
  !
  DATA first/ngrid*.TRUE./
  !
  REAL(DP), DIMENSION(ngrid) :: fact
  REAL(DP), ALLOCATABLE, TARGET, DIMENSION(:,:) :: auxp
  INTEGER, TARGET, DIMENSION(dim_iw,ngrid) :: iw0
  !
END MODULE afftnec
!
!----------------------------------------------------------------------
SUBROUTINE cft_3(f,nr1,nr2,nr3,nrx1,nrx2,nrx3,igrid,sign)
  !----------------------------------------------------------------------
  !
  !      3d fft for NEC SX6 - uses ASL library routines
  !      contributed by Guido Roma
  !
  USE kinds, ONLY : DP
  USE afftnec
  IMPLICIT NONE

  INTEGER :: &
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
  COMPLEX(DP) :: &       
       &       f(nrx1,nrx2,nrx3)    ! inp/out: the function to transform
  COMPLEX(DP), ALLOCATABLE, DIMENSION(:,:,:) :: f1 ! for ASL Library FFT routines 
#ifdef ASL
  INTEGER, POINTER, DIMENSION(:) :: iw
#if defined MICRO
  COMMON/NEC_ASL_PARA/nbtasks
  INTEGER :: nbtasks
#endif
#endif
  REAL(DP), POINTER, DIMENSION(:) :: cw1
  COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: cw2   

  !     allocate auxp at the first call (independently of the grid)
  IF (.NOT.ALLOCATED(auxp)) THEN
     ALLOCATE(auxp(195+2*(nr1+nr2+nr3),ngrid))
  END IF

  !
  !    test the sign and put the correct normalization on f
  !
  IF (first(igrid)) THEN
     nrz1(igrid)=nrx1
     nrz2(igrid)=nrx2
     nrz3(igrid)=nrx3
     IF (MOD(nrx1,2)==0) nrz1(igrid)=nrx1+1
     IF (MOD(nrx2,2)==0) nrz2(igrid)=nrx2+1
     IF (MOD(nrx3,2)==0) nrz3(igrid)=nrx3+1
  END IF
#ifdef ASL
  ALLOCATE(cw2(nrz1(igrid)*nrz2(igrid)*nrz3(igrid)))
#else
  ALLOCATE(cw2(6*nrz1(igrid)*nrz2(igrid)*nrz3(igrid)))
#endif
  ALLOCATE(f1(nrz1(igrid),nrz2(igrid),nrz3(igrid)))

  IF ( sign.EQ.-1 ) THEN
     fact(igrid)=1.0_8/DBLE(nr1*nr2*nr3)
     CALL DSCAL(2*nrx1*nrx2*nrx3,fact(igrid),f,1)
  ELSE IF ( sign.NE.1) THEN
     CALL errore('cft_3', 'wrong isign',1)
  ENDIF
  IF (igrid.LE.0.OR.igrid.GT.ngrid)&
       &  CALL errore('cft_3','which grid ?',1)

  !     copy f in the auxiliary f1 with odd dimensions
  !      call ZCOPY(nrx1*nrx2*nrx3,f,1,f1(1:nrx1,1:nrx2,1:nrx3),1)
  f1(1:nrx1,1:nrx2,1:nrx3)=f

#ifdef ASL
  CALL zfc3cl(f1,nr1,nr2,nr3,nrz1(igrid),nrz2(igrid),nrz3(igrid),ierr)
  CALL errore('cft_3', 'initialisation problem',ierr)
  iw=>iw0(:,igrid)
#endif
  cw1=>auxp(:,igrid)

  IF (first(igrid)) THEN
     isw=0
     first(igrid)=.FALSE.
     !         WRITE( stdout,*)'________________________________________________________'
     !         WRITE( stdout,*) 'igrid = ',igrid
     !         WRITE( stdout,*) '  nrxs => ',nrx1,nrx2,nrx3
     !         WRITE( stdout,*) '  nrzs => ',nrz1(igrid),nrz2(igrid),nrz3(igrid)
     !         WRITE( stdout,*) '  nrs => ',nr1,nr2,nr3
     !      WRITE( stdout,*)'size(auxp)',size(auxp,1),size(auxp,2)
     !      WRITE( stdout,*)'size(cw1)',size(cw1)
     !      WRITE( stdout,*)'size(iw)',size(iw)
     !         WRITE( stdout,*)'________________________________________________________'
#ifdef ASL
#if defined MICRO
     CALL hfc3fb(nr1,nr2,nr3,f1,nrz1(igrid),nrz2(igrid),nrz3(igrid),&
          &            isw,iw,cw1,cw2,nbtasks,ierr)
#else
     CALL zfc3fb(nr1,nr2,nr3,f1,nrz1(igrid),nrz2(igrid),nrz3(igrid),&
          &            isw,iw,cw1,cw2,ierr)
#endif
     IF (ierr.NE.0) CALL errore('cft_3','ierr=',ierr)
#else
     CALL ZZFFT3D(0,nr1,nr2,nr3,1.d0,f1,nrz1(igrid),nrz2(igrid),&
          &             f1,nrz1(igrid),nrz2(igrid),cw1,cw2,ierr)
#endif
  ENDIF

#ifdef ASL
  isw=-sign
#if defined MICRO
  CALL hfc3bf(nr1,nr2,nr3,f1,nrz1(igrid),nrz2(igrid),nrz3(igrid),&
       &            isw,iw,cw1,cw2,nbtasks,ierr)
#else
  CALL zfc3bf(nr1,nr2,nr3,f1,nrz1(igrid),nrz2(igrid),nrz3(igrid),&
       &            isw,iw,cw1,cw2,ierr)     
#endif
  IF (ierr.NE.0) CALL errore('cft_3','ierr=',ierr)
#else
  isw=sign
  CALL ZZFFT3D(isw,nr1,nr2,nr3,1.d0,f1,nrz1(igrid),nrz2(igrid),&
       &             f1,nrz1(igrid),nrz2(igrid),cw1,cw2,ierr)
#endif


  !     copy f1 back in f with odd dimensions
  !      call zcopy(nrx1*nrx2*nrx3,f1(1:nrx1,1:nrx2,1:nrx3),1,f,1)
  f(:,:,:)=f1(1:nrx1,1:nrx2,1:nrx3)
  DEALLOCATE(f1)
  DEALLOCATE(cw2)
  NULLIFY(cw1)
#ifdef ASL
  NULLIFY(iw)
#endif
  !
  RETURN
END SUBROUTINE cft_3
!
#elif defined(__SUN)
!
! ... SUN case
!
!----------------------------------------------------------------------
SUBROUTINE cft_3 (f, n1, n2, n3, nx1, nx2, nx3, igrid, sign)
  !----------------------------------------------------------------------
  !
  ! ... SUNperf library
  !
  USE kinds, ONLY : DP

  IMPLICIT NONE

  INTEGER :: n1, n2, n3, nx1, nx2, nx3, sign, igrid
  COMPLEX(DP) :: f ( nx1 , nx2 , nx3 )

  EXTERNAL zfft3i, zfft3f, zfft3b, zdscal
!
! Local variables
!
  INTEGER :: lwork
  COMPLEX(DP) ::  work ( 4*(nx1 + nx2 + nx3) + 45 )
  REAL(DP) :: scale

  lwork = 4 * ( nx1 + nx2 + nx3 ) + 45

  IF (sign.NE. - 1.AND.sign.NE.1) CALL errore ('cft_3', 'which fft ?', 1)

  CALL zfft3i ( nx1, nx2, nx3, work )

  IF( sign.EQ. 1 ) THEN
   CALL zfft3b ( n1, n2, n3, f, nx1, nx2, work, lwork )
  ELSE
   CALL zfft3f ( n1, n2, n3, f, nx1, nx2, work, lwork )
   scale = 1.0d0 /DBLE( n1 * n2 * n3 )
   CALL zdscal ( nx1*nx2*nx3, scale, f, 1 )
  ENDIF

  RETURN

END SUBROUTINE cft_3
!
#elif defined(DXML)
!
! ... DXML case
!
!----------------------------------------------------------------------
SUBROUTINE cft_3 (f, n1, n2, n3, nm1, nm2, nm3, igrid, sign)
  !----------------------------------------------------------------------
  !
  ! ... driver routine for 3d complex fft using DXML/CXML libraries
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  INTEGER :: n1, n2, n3, nm1, nm2, nm3, igrid, sign

  COMPLEX(DP) :: f (nm1, nm2, nm3)
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
  INTEGER :: ngrid
  PARAMETER (ngrid = 2)
  record / DXML_Z_FFT_STRUCTURE / fft_struct (ngrid)
  INTEGER :: status, zfft_init_3d, zfft_exit_3d
  REAL(DP) :: norm
  CHARACTER (len=1) :: direction
  LOGICAL :: first (ngrid)
  DATA first / ngrid * .TRUE. /
  SAVE first, fft_struct

  norm = DBLE (n1 * n2 * n3)
  IF (sign.EQ.1) THEN
     CALL dscal (2 * nm1 * nm2 * nm3, norm, f, 1)
     direction = 'b'
  ELSEIF (sign.EQ. - 1) THEN
     CALL dscal (2 * nm1 * nm2 * nm3, 1.0d0 / norm, f, 1)
     direction = 'f'
  ELSE
     CALL errore ('cft_3', 'sign unexpected',1)
  ENDIF

  IF (first (igrid) ) THEN
     status = zfft_exit_3d (fft_struct (igrid) )
     ! not sure whether the above call is useful or not
     status = zfft_init_3d (n1, n2, n3, fft_struct (igrid), .TRUE.)
     first (igrid) = .FALSE.
  ENDIF

  CALL zfft_apply_3d ('C', 'C', direction, f, f, nm1, nm2, &
       fft_struct (igrid) , 1, 1, 1)
  RETURN

END SUBROUTINE cft_3
!
#elif defined(FUJ64)
!
! ... FUJ64 case
!
!----------------------------------------------------------------------------
SUBROUTINE cft_3( ac, n1, n2, n3, nm1, nm2, nm3, igrid, isign )
  !----------------------------------------------------------------------------
  !
  ! ... 3d fft - FUJITSU, FFTVPLIB version (GRoma, 2001)
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER            :: n1,n2,n3,nm1,nm2,nm3,igrid,isign
  INTEGER, PARAMETER :: ngrid=2, nmax=256
  REAL(DP)      ::  ac(2,nm1,nm2,nm3)
  REAL(DP)      :: workarray(2*nm1*nm2*nm3)
  REAL(DP)      :: trig(2*3*nmax+2*120,ngrid),norm(ngrid)
  INTEGER            :: iopt,ierr,trigdim(ngrid)
  INTEGER            :: idim(3,ngrid)
  CHARACTER(LEN=2)   :: init
  LOGICAL            :: first(ngrid)
  !
  CHARACTER(LEN=1), DIMENSION(-1:1), PARAMETER :: &
    mode  = (/'m','x','p'/), &
    scale = (/'n','s','i'/)
  !
  DATA first/ngrid*.TRUE./
  !
  SAVE first, trig, idim, trigdim, norm
  !
  !
  iopt = SIGN( 1, isign )
  !
  IF ( nm1 /= n1 ) &
     CALL errore( 'cft_3', 'not any more implemented', nm1 )
  !
  IF ( igrid <= 0 .OR. igrid > ngrid ) &
     CALL errore( 'cft_3', 'which grid ?', 1 )
  !
  IF ( isign /= -1 .AND. isign /= 1 ) &
     CALL errore( 'cft_3', 'which fft ?', 2 )
  !
  IF ( n1 > nmax .OR. n2 > nmax .OR. n3 > nmax )  &
     CALL errore( 'cft_3', 'increase nmax', 3 )
  !
  init(2:2) = scale(iopt)
  !
  IF (first(igrid)) THEN
     !
     init(1:1)='i'
     !
     idim(1,igrid)=nm1
     idim(2,igrid)=nm2
     idim(3,igrid)=nm3
     !
     norm(igrid) = SQRT( REAL( nm1*nm2*nm3 ) )
     !
     trigdim(igrid)= 2*( idim(1,igrid) + idim(2,igrid) + idim(3,igrid) ) + 120
     !
     IF( n1 > nm1 ) CALL errore( 'cft3', 'too large input dimension', n1 )
     IF( n2 > nm2 ) CALL errore( 'cft3', 'too large input dimension', n2 )
     IF( n3 > nm3 ) CALL errore( 'cft3', 'too large input dimension', n3 )
     !
     ! ... The FFTVP library stores idim in a common (a single one)
     ! ... so every time you change grid you have to reinitialise!!!
     ! ... That's why the following line.
     ! ... first = .true.
     ! ... which in fact is not needed if one uses the modified version, 
     ! ... fftvplib2
     !
     first(igrid) = .FALSE.
     !
  ELSE
     !
     init(1:1) = 'r'
     !
  END IF
  !
  CALL dftcbm( ac(1,:,:,:), ac(2,:,:,:), 3, idim(:,igrid), workarray, &
               trig(:trigdim(igrid),igrid), mode(iopt), init, ierr )
  !
  CALL errore( 'cft_3', 'problems in fft', ierr )
  !
  IF ( iopt > 0 ) THEN
     !
     CALL DSCAL( 2*nm1*nm2*nm3, norm(igrid), ac, 1 )
     !
  ELSE IF ( iopt < 0 ) THEN
     !
     CALL DSCAL( 2*nm1*nm2*nm3, 1.D0 / norm(igrid), ac, 1 )
     !
  END IF
  !
  RETURN
  !
END SUBROUTINE cft_3
!
#else
!
! ... error routine
!
!----------------------------------------------------------------------------
SUBROUTINE cft_3( f, n1, n2, n3, nm1, nm2, nm3, igrid, sign )
  !---------------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  !
  !
  CALL errore( 'cft_3', 'machine-specific routine not available', 1 )
  !
  RETURN
  !
END SUBROUTINE cft_3
#endif
