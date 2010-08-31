! FOR GWW
!
! Author: P. Umari
!
SUBROUTINE gram_schmidt_pwannier(npwan, o_mat, nrmat, a_real_mat, ncmat,a_comp_mat,nrside,b_real_mat)

!this subroutine orthonormalize the pruduct of (real) wavefunctions
! and update matrix elements of operators "on the fly"
!it updates matrix B just on the right side "on the fly"
! see notes
! #ifdef __GWW

  USE kinds,      ONLY : DP
  USE io_global,   ONLY : stdout
  USE wannier_gw, ONLY : real_matrix_pointer, complex_matrix_pointer
  USE mp_global,  ONLY : mpime, nproc, world_comm
  USE mp,         ONLY : mp_sum, mp_barrier, mp_bcast


  implicit none

  INTEGER, INTENT(in)          :: npwan!number of product functions
  REAL(kind=DP), INTENT(inout) :: o_mat(npwan,npwan)!overlap matrix (it's symmetric)
  INTEGER, INTENT(in)          :: nrmat!number of real matrices to be updated
  TYPE(real_matrix_pointer)    :: a_real_mat(nrmat)!pointers to real matrices
  INTEGER, INTENT(in)          :: ncmat!number of complex matrices to be updated
  TYPE(complex_matrix_pointer) :: a_comp_mat(ncmat)!pointers to complex matrices
  INTEGER, INTENT(in)          :: nrside!number of real matrices to be updated just on the right
  TYPE(real_matrix_pointer)    :: b_real_mat(nrside)!pointers to real matrices to be update just on the right


  INTEGER :: i,j,n, nr,nc
  REAL(kind=DP) :: alfa
  REAL(kind=DP), ALLOCATABLE :: r_row(:), r_column(:)!temporary arrays for updating
  COMPLEX(kind=DP), ALLOCATABLE :: c_row(:), c_column(:)!temporary arrays for updating
  REAL(kind=DP), EXTERNAL :: ddot


  ALLOCATE(r_row(npwan),r_column(npwan),c_row(npwan),c_column(npwan))

  do n=1,npwan!external loop
! calculates alfa

    alfa=0.d0
    write(stdout,*) 'alfa0',n,alfa
    do i=1,n-1
       if(mod(i,nproc)==mpime) then
          alfa=alfa-2.d0*o_mat(i,n)**2.
!          do j=1,n-1
!             alfa=alfa+o_mat(i,n)*o_mat(j,n)*o_mat(j,i)
!          enddo
          alfa=alfa+o_mat(i,n)*ddot(n-1,o_mat(:,n),1,o_mat(:,i),1)

        endif
    enddo
    call mp_sum(alfa, world_comm)
    alfa=alfa+o_mat(n,n)
    write(stdout,*) 'ALFA', alfa!ATTENZIONE

    alfa=1.d0/sqrt(alfa)

    write(stdout,*) '1/ALFA', alfa!ATTENZIONE

!updates real matrices

    do nr=1,nrmat!loop on matrices
      r_row(:)=0.d0
      r_column(:)=0.d0

!term A_{n,n}
      r_row(n)=0.d0

!      do i=1,n-1
!         if(mod(i,nproc)==mpime) then
!            r_row(n)=r_row(n)-o_mat(i,n)*a_real_mat(nr)%p(i,n) &
!                 & -o_mat(i,n)*a_real_mat(nr)%p(n,i)
!            do j=1,n-1
!               r_row(n)=r_row(n)+o_mat(i,n)*o_mat(j,n)*a_real_mat(nr)%p(i,j)
!            enddo
!         endif
!      enddo


      do j=1,n-1
         if(mod(j,nproc)==mpime) then
            do i=1,n-1
               r_row(n)=r_row(n)+o_mat(i,n)*o_mat(j,n)*a_real_mat(nr)%p(i,j)
            enddo
         endif
      enddo
      call mp_sum(r_row(n), world_comm)
      do i=1,n-1
          r_row(n)=r_row(n)-o_mat(i,n)*a_real_mat(nr)%p(i,n) &
                 & -o_mat(i,n)*a_real_mat(nr)%p(n,i)
       enddo


      r_row(n)=r_row(n)+a_real_mat(nr)%p(n,n)
      r_row(n)=r_row(n)*(alfa**2.d0)
      r_column(n)=r_row(n)

!terms A_{n,i}, A_{i,n}, i/=n

!      do i=1,npwan
!        if(i/=n) then
!           if(mod(i,nproc)==mpime) then
!              r_row(i)=a_real_mat(nr)%p(n,i)
!              r_column(i)=a_real_mat(nr)%p(i,n)
!              do j=1,n-1
!                 r_row(i)=r_row(i)-o_mat(j,n)*a_real_mat(nr)%p(j,i)
!                 r_column(i)=r_column(i)-o_mat(j,n)*a_real_mat(nr)%p(i,j)
!              enddo
!              r_row(i)=r_row(i)*alfa
!              r_column(i)=r_column(i)*alfa
!           endif
!        endif
!     enddo
!     call mp_barrier
!     do  i=1,npwan
!        if(i/=n) then
!           call mp_sum(r_row(i), group)!to distribute among processors
!           call mp_sum(r_column(i), group)!to distribute among processors
!        endif
!     enddo

      do i=1,npwan
         if(i/=n) then
            if(mod(i,nproc)==mpime) then
               r_row(i)=a_real_mat(nr)%p(n,i)
               do j=1,n-1
                  r_row(i)=r_row(i)-o_mat(j,n)*a_real_mat(nr)%p(j,i)
               enddo
               r_row(i)=r_row(i)*alfa
            endif
         endif
      enddo

      do i=1,npwan
         if(i/=n) then
            call mp_bcast(r_row(i),mod(i,nproc))
         endif
      enddo

      do j=1,n-1
         if(mod(j,nproc)==mpime) then
            do i=1,npwan
               if(i/=n) then
                  r_column(i)=r_column(i)-o_mat(j,n)*a_real_mat(nr)%p(i,j)
               endif
            enddo
         endif
      enddo
      do i=1,npwan
         if(i/=n) then
            call mp_sum(r_column(i))
            r_column(i)=r_column(i)+a_real_mat(nr)%p(i,n)
            r_column(i)=r_column(i)*alfa
         endif
      enddo



! updates A
     a_real_mat(nr)%p(n,:)=r_row(:)
     a_real_mat(nr)%p(:,n)=r_column(:)
  enddo

!updates complex matrices

    do nc=1,ncmat!loop on matrices
      c_row(:)=0.d0
      c_column(:)=0.d0

!term A_{n,n}
      c_row(n)=(0.d0,0.d0)
      do i=1,n-1
         if(mod(i,nproc)==mpime) then
            c_row(n)=c_row(n)-o_mat(i,n)*a_comp_mat(nc)%p(i,n) &
                 & -o_mat(i,n)*a_comp_mat(nc)%p(n,i)
            do j=1,n-1
               c_row(n)=c_row(n)+o_mat(i,n)*o_mat(j,n)*a_comp_mat(nc)%p(i,j)
            enddo
         endif
      enddo
      call mp_barrier
      call mp_sum(c_row(n), world_comm)
      c_row(n)=c_row(n)+a_comp_mat(nc)%p(n,n)
      call mp_sum(c_column(n), world_comm)
      c_row(n)=c_row(n)*(alfa**2.)
      c_column(n)=c_row(n)

!terms A_{n,i}, A_{i,n}, i/=n
      do i=1,npwan
        if(i/=n) then
           if(mod(i,nproc)==mpime) then
              c_row(i)=a_comp_mat(nc)%p(n,i)
              c_column(i)=a_comp_mat(nc)%p(i,n)
              do j=1,n-1
                 c_row(i)=c_row(i)-o_mat(j,n)*a_comp_mat(nc)%p(j,i)
                 c_column(i)=c_column(i)-o_mat(j,n)*a_comp_mat(nc)%p(i,j)
              enddo
              c_row(i)=c_row(i)*alfa
              c_column(i)=c_column(i)*alfa
           endif
        endif
     enddo
     call mp_barrier
     do i=1,npwan
        if(i/=n) then
           call mp_sum(c_row(i), world_comm)!for distributing among processors
           call mp_sum(c_column(i), world_comm)!for distributing among processors
        endif
     enddo
! updates A
     a_comp_mat(nc)%p(n,:)=c_row(:)
     a_comp_mat(nc)%p(:,n)=c_column(:)
  enddo

!now updates real B matrices on the right

   do nr=1,nrside!loop on matrices
      r_column(:)=0.d0
!terms b_{i,n}
!      do i=1,npwan
!         if(mod(i,nproc)==mpime) then
!            r_column(i)=b_real_mat(nr)%p(i,n)
!            do j=1,n-1
!               r_column(i)=r_column(i)-o_mat(j,n)*b_real_mat(nr)%p(i,j)
!            enddo
!            r_column(i)=r_column(i)*alfa
!         endif
!      enddo
!      call mp_barrier
!      do i=1,npwan
!         call mp_sum(r_column(i),  group)!for distributing among processors
!      enddo

      do j=1,n-1
         if((mod(j,nproc)==mpime)) then
            do i=1,npwan
               r_column(i)=r_column(i)-o_mat(j,n)*b_real_mat(nr)%p(i,j)
            enddo
         endif
      enddo
      do i=1,npwan
         call mp_sum(r_column(i))
         r_column(i)=r_column(i)+b_real_mat(nr)%p(i,n)
         r_column(i)=r_column(i)*alfa
      enddo



!updates B
      b_real_mat(nr)%p(:,n)=r_column(:)
   enddo




!updates o_mat
      r_row(:)=0.d0

!term o_{n,n}
      r_row(n)=1.d0

!terms o_{i<n,n}, o_{n,i<n}
      do i=1,n-1
         r_row(i)=0.d0
      enddo
!terms o_{n,i}, o_{i,n}, i/=n
      do i=n+1,npwan
        if(i/=n) then
           if(mod(i,nproc)==mpime) then
              r_row(i)=o_mat(n,i)
!              do j=1,n-1
!                 r_row(i)=r_row(i)-o_mat(j,n)*o_mat(j,i)
!              enddo
              r_row(i)=r_row(i)-ddot(n-1,o_mat(:,n),1,o_mat(:,i),1)
              r_row(i)=r_row(i)*alfa
           endif
        endif
     enddo
     do i=n+1,npwan
        call mp_sum(r_row(i), world_comm)
     enddo
     o_mat(n,:)=r_row(:)
     o_mat(:,n)=r_row(:)

  enddo!on n





  DEALLOCATE(r_row,r_column,c_row,c_column)
! #endif
END SUBROUTINE gram_schmidt_pwannier
