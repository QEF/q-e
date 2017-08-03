!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=---------------------------------------------------------------------==!
!
!     Parallel 3D FFT high level Driver
!     ( Charge density and Wave Functions )
!
!     Written and maintained by Carlo Cavazzoni !     Last update Apr. 2009
!     Rewritten  by Stefano de Gironcoli        !     Sep-Nov 2016
!
!!=---------------------------------------------------------------------==!
!
MODULE fft_parallel
!
   USE fft_param
   IMPLICIT NONE
   SAVE
!
CONTAINS
!
!  General purpose driver
!
!----------------------------------------------------------------------------
SUBROUTINE tg_cft3s( f, dfft, isgn )
  !----------------------------------------------------------------------------
  !
  !! ... isgn = +-1 : parallel 3d fft for rho and for the potential
  !
  !! ... isgn = +-2 : parallel 3d fft for wavefunctions
  !
  !! ... isgn = +-3 : parallel 3d fft for wavefunctions with task group
  !
  !! ... isgn = +   : G-space to R-space, output = \sum_G f(G)exp(+iG*R)
  !! ...              fft along z using pencils        (cft_1z)
  !! ...              transpose across nodes           (fft_scatter_yz)
  !! ...              fft along y using pencils        (cft_1y)
  !! ...              transpose across nodes           (fft_scatter_xy)
  !! ...              fft along x using pencils        (cft_1x)
  !
  !! ... isgn = -   : R-space to G-space, output = \int_R f(R)exp(-iG*R)/Omega
  !! ...              fft along x using pencils        (cft_1x)
  !! ...              transpose across nodes           (fft_scatter_xy)
  !! ...              fft along y using pencils        (cft_1y)
  !! ...              transpose across nodes           (fft_scatter_yz)
  !! ...              fft along z using pencils        (cft_1z)
  !
  ! If task_group_fft_is_active the FFT acts on a number of wfcs equal to 
  ! dfft%nproc2, the number of Y-sections in which a plane is divided. 
  ! Data are reshuffled by the fft_scatter_tg routine so that each of the 
  ! dfft%nproc2 subgroups (made by dfft%nproc3 procs) deals with whole planes 
  ! of a single wavefunciton.
  !
  ! This driver is based on code written by Stefano de Gironcoli for PWSCF.
  !
  USE fft_scalar, ONLY : cft_1z
  USE scatter_mod,ONLY : fft_scatter_xy, fft_scatter_yz, fft_scatter_tg
  USE fft_types,  ONLY : fft_type_descriptor
  !
  IMPLICIT NONE
  !
  TYPE (fft_type_descriptor), INTENT(in) :: dfft   ! descriptor of fft data layout
  COMPLEX(DP), INTENT(inout)             :: f( : ) ! array containing data to be transformed
  INTEGER, INTENT(in)                    :: isgn   ! fft direction (potential: +/-1, wave: +/-2, wave_tg: +/-3)
  !
  INTEGER                  :: n1, n2, n3, nx1, nx2, nx3
  INTEGER                  :: nnr_
  INTEGER                  :: nsticks_x, nsticks_y, nsticks_z
  COMPLEX(DP), ALLOCATABLE :: aux (:)
  INTEGER                  :: i
  !
  !write (6,*) 'enter tg_cft3s ',isgn ; write(6,*) ; FLUSH(6)
  n1  = dfft%nr1  ; n2  = dfft%nr2  ; n3  = dfft%nr3
  nx1 = dfft%nr1x ; nx2 = dfft%nr2x ; nx3 = dfft%nr3x
  !
  if (abs(isgn) == 1 ) then       ! potential fft
     nnr_ = dfft%nnr
     nsticks_x = dfft%my_nr2p * dfft%my_nr3p
     nsticks_y = dfft%nr1p(dfft%mype2+1) * dfft%my_nr3p
     nsticks_z = dfft%nsp(dfft%mype+1)
  else if (abs(isgn) == 2 ) then  ! wave func fft
     nnr_ = dfft%nnr
     nsticks_x = dfft%my_nr2p * dfft%my_nr3p
     nsticks_y = dfft%nr1w(dfft%mype2+1) * dfft%my_nr3p
     nsticks_z = dfft%nsw(dfft%mype+1)
  else if (abs(isgn) == 3 ) then  ! wave func fft with task groups
     nnr_ = dfft%nnr_tg
     nsticks_x = dfft%nr2 * dfft%my_nr3p
     nsticks_y = dfft%nr1w_tg * dfft%my_nr3p
     nsticks_z = dfft%nsw_tg(dfft%mype+1)
  else
     CALL fftx_error__( ' tg_cft3s', ' wrong value of isgn ', 10+abs(isgn) )
  end if
  ALLOCATE( aux( nnr_ ) ) !;  aux = (0.d0,0.d0)
  !

  !write (6,99) f(1:400); write(6,*); FLUSH(6)
  IF ( isgn > 0 ) THEN  ! G -> R
  !write (6,*) 'enter scatter_tg' ; write(6,*); FLUSH(6)
     if (isgn==+3) call fft_scatter_tg ( dfft, f, aux, nnr_, isgn)
  !write (6,99) f(1:860); write(6,*); FLUSH(6)
     !
  !write (6,*) 'enter cft_1z' ; write(6,*); FLUSH(6)
     CALL cft_1z( f, nsticks_z, n3, nx3, isgn, aux )
  !write (6,99) aux(1:860); write(6,*); FLUSH(6)
     !
  !write (6,*) 'enter scatter_yz' ; write(6,*); FLUSH(6)
     CALL fft_scatter_yz ( dfft, aux, f, nnr_, isgn )
  !write (6,99) f(1:400); write(6,*); FLUSH(6)
     !
  !write (6,*) 'enter cft_1y' ; write(6,*); FLUSH(6)
     CALL cft_1z( f, nsticks_y, n2, nx2, isgn, aux )
  !write (6,99) aux(1:400); write(6,*); FLUSH(6)
     !
  !write (6,*) 'enter scatter_xy' ; write(6,*); FLUSH(6)
     CALL fft_scatter_xy ( dfft, aux, f, nnr_, isgn )
  !write (6,*) 'f after  scatter_xy' ; write(6,*); FLUSH(6)
  !write (6,99) f(1:400); write(6,*); FLUSH(6)
     !
  !write (6,*) 'enter cft_1x' ; write(6,*); FLUSH(6)
     CALL cft_1z( f, nsticks_x, n1, nx1, isgn, aux )
  !write (6,*) 'aux after  cft_1x' ; write(6,*); FLUSH(6)
  !write (6,99) aux(1:400); write(6,*); FLUSH(6)
     ! clean garbage beyond the intended dimension. should not be needed but apparently it is !
     if (nsticks_x*nx1 < nnr_) aux(nsticks_x*nx1+1:nnr_) = (0.0_DP,0.0_DP)
     !write (6,*) size(f), size(aux),nsticks_x*nx1,dfft%nnr ; FLUSH(6)
     f(1:nnr_)=aux(1:nnr_) ! not limiting the range to dfft%nnr may crash when size(f)>size(aux)
     !
  ELSE                  ! R -> G
     !
  !write (6,*) 'enter cft_1x' ; write(6,*); FLUSH(6)
     CALL cft_1z( f, nsticks_x, n1, nx1, isgn, aux )
  !write (6,99) aux(1:400); write(6,*); FLUSH(6)
     !
  !write (6,*) 'enter scatter_xy' ; write(6,*); FLUSH(6)
     CALL fft_scatter_xy ( dfft, f, aux, nnr_, isgn )
  !write (6,99) f(1:400); write(6,*); FLUSH(6)
     !
  !write (6,*) 'enter cft_1y' ; write(6,*); FLUSH(6)
     CALL cft_1z( f, nsticks_y, n2, nx2, isgn, aux )
  !write (6,99) aux(1:400); write(6,*); FLUSH(6)
     !
  !write (6,*) 'enter scatter_yz' ; write(6,*); FLUSH(6)
     CALL fft_scatter_yz ( dfft, f, aux, nnr_, isgn )
  !write (6,*) 'f after  scatter_yz' ; write(6,*); FLUSH(6)
  !write (6,99) f(1:400); write(6,*); FLUSH(6)
     !
  !write (6,*) 'enter cft_1z' ; write(6,*); FLUSH(6)
     CALL cft_1z( f, nsticks_z, n3, nx3, isgn, aux )
     ! clean garbage beyond the intended dimension. should not be needed but apparently it is !
     if (nsticks_z*nx3 < nnr_) aux(nsticks_z*nx3+1:nnr_) = (0.0_DP,0.0_DP)
     f(1:nnr_)=aux(1:nnr_) ! not limiting the range to dfft%nnr may crash when size(f)>size(aux)
     !
  !write (6,*) 'enter scatter_tg' ; write(6,*); FLUSH(6)
     if (isgn==-3) call fft_scatter_tg ( dfft, f, aux, nnr_, isgn)
  ENDIF
  !write (6,99) f(1:400); write(6,*); FLUSH(6)
  !
  DEALLOCATE( aux )
  !
  !if (.true.) stop
  RETURN
99 format ( 20 ('(',2f12.9,')') )

  !
END SUBROUTINE tg_cft3s

!--------------------------------------------------------------------------------
!   Auxiliary routines to read/write from/to a distributed array
!   NOT optimized for efficiency .... just to show how one can access the data
!--------------------------------------------------------------------------------
!
COMPLEX (DP) FUNCTION get_f_of_R (i,j,k,f,dfft)
!------  read from a distributed complex array f(:) in direct space
!
  USE fft_types,  ONLY : fft_type_descriptor
  IMPLICIT NONE
  TYPE (fft_type_descriptor), INTENT(IN) :: dfft
  INTEGER, INTENT (IN) :: i,j,k
  COMPLEX(DP), INTENT (IN) :: f(:)
  INTEGER :: kk, ii, jj, ip, ierr
  COMPLEX(DP) :: f_aux

  IF ( i <= 0 .OR. i > dfft%nr1 ) CALL fftx_error__( ' get_f_of_R', ' first  index out of range ', 1 )
  IF ( j <= 0 .OR. j > dfft%nr2 ) CALL fftx_error__( ' get_f_of_R', ' second index out of range ', 2 )
  IF ( k <= 0 .OR. k > dfft%nr3 ) CALL fftx_error__( ' get_f_of_R', ' third  index out of range ', 3 )

#if defined(__MPI)
  do ip = 1, dfft%nproc3
     if ( dfft%i0r3p(ip) < k ) kk = ip
  end do
  do ip = 1, dfft%nproc2
     if ( dfft%i0r2p(ip) < j ) jj = ip
  end do
  ii  = i + dfft%nr1x * ( j - dfft%i0r2p(jj) - 1 ) + dfft%nr1x * dfft%nr2p(jj) * ( k - dfft%i0r3p(kk) - 1 )
  f_aux = (0.d0,0.d0)
  if ( (jj == (dfft%mype2 + 1)) .and. (kk == (dfft%mype3 +1)) ) f_aux = f(ii)
  CALL MPI_ALLREDUCE( f_aux, get_f_of_R,   2, MPI_DOUBLE_PRECISION, MPI_SUM, dfft%comm, ierr )
#else
  ii = i + dfft%nr1x * (j-1) + dfft%nr1x*dfft%nr2x * (k-1)
  get_f_of_R = f(ii)
#endif
END FUNCTION get_f_of_R

SUBROUTINE put_f_of_R (f_in,i,j,k,f,dfft)
!------  write on a distributed complex array f(:) in direct space
!
  USE fft_types,  ONLY : fft_type_descriptor
  IMPLICIT NONE
  TYPE (fft_type_descriptor), INTENT(IN) :: dfft
  INTEGER, INTENT (IN) :: i,j,k
  COMPLEX(DP), INTENT (IN) :: f_in
  COMPLEX(DP), INTENT (INOUT) :: f(:)
  INTEGER :: kk, ii, jj, ip, ierr

  IF ( i <= 0 .OR. i > dfft%nr1 ) CALL fftx_error__( ' put_f_of_R', ' first  index out of range ', 1 )
  IF ( j <= 0 .OR. j > dfft%nr2 ) CALL fftx_error__( ' put_f_of_R', ' second index out of range ', 2 )
  IF ( k <= 0 .OR. k > dfft%nr3 ) CALL fftx_error__( ' put_f_of_R', ' third  index out of range ', 3 )

#if defined(__MPI)
  do ip = 1, dfft%nproc3
     if ( dfft%i0r3p(ip) < k ) kk = ip
  end do
  do ip = 1, dfft%nproc2
     if ( dfft%i0r2p(ip) < j ) jj = ip
  end do
  ii  = i + dfft%nr1x * ( j - dfft%i0r2p(jj) - 1 ) + dfft%nr1x * dfft%nr2p(jj) * ( k - dfft%i0r3p(kk) - 1 )
  if ( (jj == (dfft%mype2 + 1)) .and. (kk == (dfft%mype3 +1)) ) f(ii) = f_in
#else
  ii = i + dfft%nr1x * (j-1) + dfft%nr1x*dfft%nr2x * (k-1)
  f(ii) = f_in
#endif

END SUBROUTINE put_f_of_R

COMPLEX (DP) FUNCTION get_f_of_G (i,j,k,f,dfft)
!------  read from a distributed complex array f(:) in reciprocal space
!
  USE fft_types,  ONLY : fft_type_descriptor
  IMPLICIT NONE
  INTEGER, INTENT (IN) :: i,j,k
  COMPLEX(DP), INTENT (IN) :: f(:)
  TYPE (fft_type_descriptor), INTENT(IN) :: dfft
  INTEGER :: ii, jj, ip, ierr
  COMPLEX(DP) :: f_aux

  IF ( i <= 0 .OR. i > dfft%nr1 ) CALL fftx_error__( ' get_f_of_G', ' first  index out of range ', 1 )
  IF ( j <= 0 .OR. j > dfft%nr2 ) CALL fftx_error__( ' get_f_of_G', ' second index out of range ', 2 )
  IF ( k <= 0 .OR. k > dfft%nr3 ) CALL fftx_error__( ' get_f_of_G', ' third  index out of range ', 3 )

#if defined(__MPI)
  ii = i + dfft%nr1x * (j -1)
  jj = dfft%isind(ii)  ! if jj is zero this G vector does not belong to this processor
  f_aux = (0.d0,0.d0)
  if ( jj > 0 ) f_aux = f( k + dfft%nr3x * (jj -1))
  CALL MPI_ALLREDUCE( f_aux, get_f_of_G,   2, MPI_DOUBLE_PRECISION, MPI_SUM, dfft%comm, ierr )
#else
  ii = i + dfft%nr1 * (j-1) + dfft%nr1*dfft%nr2 * (k-1)
  get_f_of_G = f(ii)
#endif
END FUNCTION get_f_of_G

SUBROUTINE put_f_of_G (f_in,i,j,k,f,dfft)
!------  write on a distributed complex array f(:) in reciprocal space
!
  USE fft_types,  ONLY : fft_type_descriptor
  IMPLICIT NONE
  COMPLEX(DP), INTENT (IN) :: f_in
  INTEGER, INTENT (IN) :: i,j,k
  COMPLEX(DP), INTENT (INOUT) :: f(:)
  TYPE (fft_type_descriptor), INTENT(IN) :: dfft
  INTEGER :: ii, jj

  IF ( i <= 0 .OR. i > dfft%nr1 ) CALL fftx_error__( ' put_f_of_G', ' first  index out of range ', 1 )
  IF ( j <= 0 .OR. j > dfft%nr2 ) CALL fftx_error__( ' put_f_of_G', ' second index out of range ', 2 )
  IF ( k <= 0 .OR. k > dfft%nr3 ) CALL fftx_error__( ' put_f_of_G', ' third  index out of range ', 3 )

#if defined(__MPI)
  ii = i + dfft%nr1x * (j -1)
  jj = dfft%isind(ii)  ! if jj is zero this G vector does not belong to this processor
  if ( jj > 0 )   f( k + dfft%nr3x * (jj -1)) = f_in
#else
  ii = i + dfft%nr1 * (j-1) + dfft%nr1*dfft%nr2 * (k-1)
  f(ii) = f_in
#endif
END SUBROUTINE put_f_of_G

END MODULE fft_parallel
