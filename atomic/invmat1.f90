!----------------------------------------------------------------------
      subroutine invmat1( a, a_inv, lda, n )
!----------------------------------------------------------------------
! compute the inverse of a -> a_inv
!
! In all machines we use the lapack routines. Different names are to
! be used for single or double precision 
!
      implicit none

      integer, parameter :: lwork=100, dp=kind(1.d0)

      integer  :: &
            info, &      ! flag which says if the inversion was ok
            n,    &      ! the physical dimension
            lda,  &      ! the leading dimension 
            ipiv(lwork) ! a pivoting variable

      real(kind=dp) :: &
            a(lda,n),  &   ! the matrix to invert
            a_inv(lda,n), &! the inverse
            work(lwork) ! some work space

      if(n.gt.lwork)call errore('invmat','work space is too small ',n)
      call DCOPY(lda*n,a,1,a_inv,1)
      call DGETRF(n,n,a_inv,lda,ipiv,info)
      call errore('invmat','error in DGETRF',abs(info))
      call DGETRI(n,a_inv,lda,ipiv,work,lwork,info)
      call errore('invmat','error in DGETRI',abs(info))

      return
      end

