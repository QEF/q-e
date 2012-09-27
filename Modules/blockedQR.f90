!
module blockedQR
!
implicit none

PRIVATE


real*8 dnrm2,ddot
EXTERNAL dnrm2,dgemv,dscal,dger,daxpy,ddot
#ifdef __ELPA
public :: qr_rank2_real
contains

subroutine qr_rank2_real(a, lda, vmr, ldv, tmat, nbw, istep, cols, nblk, mpi_comm_rows, mpi_comm_cols, work, eps)
    use mpi

    implicit none

    integer lda, ldv, nbw, nblk, mpi_comm_rows, mpi_comm_cols, cols, istep
    real*8 :: a(lda,*), vmr(ldv,*), tmat(nbw,nbw,*), work(*)
    real*8 eps

    ! matrix organization
    integer ncol
    integer firstcol
    integer nrow
    integer local_rows 
    integer local_columns

    integer top11_row
    integer top11_col
 
    integer top22_row
    integer top22_col
 
    integer top11_column_comm_rank
    integer top11_row_comm_rank
    integer top22_column_comm_rank
    integer top22_row_comm_rank
    integer row
    integer n_cols
    
    real*8 top11
    real*8 top12
    real*8 top21
    real*8 top22
    real*8 tau1
    real*8 tau2
    real*8 beta1
    real*8 beta2
    real*8 alpha
    real*8 z2_div_tau1
    real*8 z2
    real*8 dot_second_second
    real*8 dot_first_first
    real*8 dot_first_second
 
    ! MPI stuff 
    integer column_comm_rank
    integer column_processes
    integer row_comm_rank
    integer row_processes
    integer mpierr

    ! BLAS stuff
    integer blas_rows
    integer blas_columns

    n_cols = cols

    ! in column communicator
    call mpi_comm_rank(mpi_comm_rows,column_comm_rank,mpierr)
    call mpi_comm_size(mpi_comm_rows,column_processes,mpierr)

    ! in row communicator 
    call mpi_comm_rank(mpi_comm_cols,row_comm_rank,mpierr)
    call mpi_comm_size(mpi_comm_cols,row_processes,mpierr)
 
    firstcol = istep*nbw + 1

    ! apply 2-blocked QR */
    do while (n_cols >= 2)
        ncol = istep*nbw + n_cols
        nrow = ncol - nbw

        if (nrow == 1) return ! nothing left to do

        ! figure out processes with "top" elements
        call mapEntrytoProcess(column_processes, 0, nrow-1, 0, nblk, top11_row, top11_column_comm_rank)
        call mapEntrytoProcess(row_processes, 0, ncol-1, 0, nblk, top11_col, top11_row_comm_rank)

        call mapEntrytoProcess(column_processes, 0, nrow-2, 0, nblk, top22_row, top22_column_comm_rank)
        call mapEntrytoProcess(row_processes, 0, ncol-2, 0, nblk, top22_col, top22_row_comm_rank)
 
        ! figure out size of sub matrix 
        call matrixGridInfo(local_rows, nrow, nblk, column_processes, column_comm_rank)
        call matrixGridInfo(local_columns, ncol, nblk, row_processes, row_comm_rank)
        call matrixGridInfo(blas_columns, firstcol-1, nblk, row_processes, row_comm_rank)
 
        local_columns = local_columns - blas_columns

        if ((top22_row_comm_rank .ne. top11_row_comm_rank) .or. (nrow == 2)) then
            ! needed vectors are spread across two different processor columns, 
            ! or meaningful data is only available for one vector
            ! => falling back to single householder vector generation
            if (top11_row_comm_rank == row_comm_rank) then
                ! calculate HV and tau 
                tau1 = houseGen_reversed(a(1,top11_col), nrow, mpi_comm_rows, nblk, work)
     
                ! store tau in seperate array
                work(1) = tau1
                work(2:local_rows+1) = a(1:local_rows,top11_col)

                ! update all following columns
                local_columns = local_columns - 1
            endif 

            ! broadcast HV and tau to other process columns
            call mpi_bcast(work, local_rows+1, mpi_real8, top11_row_comm_rank, mpi_comm_cols,mpierr)
  
            tau1 = work(1)
            tmat(n_cols,n_cols,istep)=tau1
            vmr(1:local_rows,n_cols) = work(2:local_rows+1)
  
            if (top11_column_comm_rank == column_comm_rank) then
                vmr(top11_row,n_cols)=1.0
            endif

            ! Apply HV to all other (local) columns
            if (local_columns > 0) then 
               call houseLeft_reversed(vmr(1,n_cols), ldv, nrow, tau1, a(1,blas_columns+1), lda, local_columns, mpi_comm_rows, nblk, work)
            endif

            ! go to next column
            n_cols = n_cols - 1
            
            cycle
        endif

        call matrixGridInfo(blas_rows, nrow-1, nblk, column_processes, column_comm_rank) ! matrix without top11/top12 elements 
 
        if (top11_row_comm_rank == row_comm_rank) then
            ! both vectors are available in current process column
            ! figure out size of sub matrix
    
            ! first collect needed data
            if (top11_column_comm_rank == column_comm_rank) then
                ! i got the top elements of the current column */
                top12 = a(top11_row,top22_col)
                top11 = a(top11_row,top22_col+1)

            else
                ! no special vector data
                top12 = 0.0
                top11 = 0.0
            endif

            if (top22_column_comm_rank == column_comm_rank) then
                top22 = a(top22_row,top22_col)
                top21 = a(top22_row,top22_col+1)
            else
                ! no special vector data
                top22 = 0.0
                top21 = 0.0
            endif

            if (blas_rows > 0) then ! only calculate if there is actual data, 
                                    ! but still take part in Allreduce later on */
                dot_second_second = dot_product(a(1:blas_rows,top22_col), a(1:blas_rows,top22_col))
                dot_first_first = dot_product(a(1:blas_rows,top22_col+1), a(1:blas_rows,top22_col+1))
                dot_first_second = dot_product(a(1:blas_rows,top22_col), a(1:blas_rows,top22_col+1))
            else
                dot_second_second = 0.0
                dot_first_first = 0.0
                dot_first_second = 0.0
            endif
 

            ! prepare allreduce buffer 
            work(1) = top11
            work(2) = top21
            work(3) = top12
            work(4) = top22
            work(5) = dot_first_first
            work(6) = dot_second_second
            work(7) = dot_first_second
            work(8) = 0.0 ! alignment

            ! step 1.5: allreduce data: dot_first_first, dot_second_second, dot_first_second, top11, top12, top21, top22, <null> (alignment) 
            call mpi_allreduce(work(1), work(9), 8, mpi_real8, MPI_SUM, mpi_comm_rows, mpierr)
            
            ! get reduced values
            top11 = work(9)
            top21 = work(10)
            top12 = work(11)
            top22 = work(12)
            dot_first_first = work(13)
            dot_second_second = work(14)
            dot_first_second = work(15)

            ! step 2: build householder vectors (completely independent)
            ! step 2.1: build beta1
            beta1 = dot_first_first + top11*top11
            if (top11 < 0.0) then 
                beta1 = -sqrt(beta1)
            else 
                beta1 = sqrt(beta1)
            endif

            if (abs((dot_first_second+top11*top12)*(dot_first_second+top11*top12)/((dot_first_first+top11*top11)*(dot_second_second+top12*top12))) > (eps/(1.0d0+eps))) then
                tau1 = top11 + beta1
 
                ! update columns and build householder vector 
                if (top11_column_comm_rank == column_comm_rank) then
                    a(top11_row,top22_col+1)=-beta1
                endif
     
                ! daxpy emulation, seems to be more efficient than actual daxpy
                ! call
                alpha = 1.0/tau1
                a(1:blas_rows,top22_col+1) = a(1:blas_rows,top22_col+1) * alpha

                tau1 = tau1/beta1
                ! apply houseleft operation to next column
                ! build y' * <second column> */
                if (top11_column_comm_rank == column_comm_rank) then
                    dot_first_second = a(top11_row,top22_col)
                else
                    dot_first_second = 0.0
                endif

                do row=1,blas_rows
                   dot_first_second =  dot_first_second + a(row,top22_col) * a(row,top22_col+1)
                enddo

                work(1) = dot_first_second
                
                call mpi_allreduce(work(1), work(3), 1, mpi_real8, MPI_SUM, mpi_comm_rows,mpierr)
                
                dot_first_second = work(3)
                dot_first_second = dot_first_second * tau1

                ! update column and build dot_second_second on the fly */
                if (top11_column_comm_rank == column_comm_rank) then
                    a(top11_row,top22_col) = a(top11_row,top22_col) - dot_first_second
                endif
                
                dot_second_second = 0.0
                do row=1,blas_rows
                     ! use top22 as temp variable
                     top22 = a(row,top22_col)  - a(row,top22_col+1) * dot_first_second
                     dot_second_second = dot_second_second + top22 * top22
                     a(row,top22_col) = top22
                enddo
 
                ! exchange top22 and dot_second_second
                if (top22_column_comm_rank == column_comm_rank) then
                    work(1) = a(top22_row,top22_col)
                else
                    work(1) = 0.0
                endif

                work(2) = dot_second_second
 
                call mpi_allreduce(work(1), work(3), 2, mpi_real8, MPI_SUM, mpi_comm_rows, mpierr)

                top22 = work(3)
                dot_second_second = work(4)
 
                ! build second householder vector 
                if (top22 < 0.0) then
                    beta2 = -sqrt(dot_second_second)
                else 
                    beta2 = sqrt(dot_second_second)
                endif
 
                tau2 = beta2 + top22
                if (top22_column_comm_rank == column_comm_rank) then
                    a(top22_row,top22_col) = -beta2
                    blas_rows = blas_rows - 1
                endif
 
                alpha = 1.0/tau2
                do row=1,blas_rows
                    a(row,top22_col) = a(row,top22_col) * alpha
                enddo
                tau2 = tau2/beta2
            else
                ! step 2.2: build temporary tau1 and tau2 (version 1)
                tau1 = top11 + beta1
                z2 = dot_first_second/beta1 + top12*(1+top11/beta1)
                z2_div_tau1 = z2 / tau1
                tau2 = top22 - z2_div_tau1*top21

                ! step 2.3: build beta2
                beta2 = dot_second_second
                beta2 = beta2 - (z2_div_tau1)*(2.0*dot_first_second - (z2_div_tau1)*(dot_first_first))

                if (tau2 < 0.0) then
                    beta2 = -sqrt(beta2)
                else 
                    beta2 = sqrt(beta2)
                endif

                ! step 2.4: build temporary tau2 (version 2)
                tau2 = tau2 + beta2
  
                ! step 2.5: update columns and build householder vectors
                if (top11_column_comm_rank == column_comm_rank) then
                    a(top11_row,top22_col+1)= -beta1
                endif
     
                ! daxpy emulation, seems to be more efficient than actual daxpy
                ! call
                alpha = 1.0/tau1
                a(1:blas_rows,top22_col+1) = a(1:blas_rows,top22_col+1) * alpha

                ! step 2.6.1: build real tau1
                tau1 = tau1/beta1
      
                if (top11_column_comm_rank == column_comm_rank) then
                    a(top11_row,top22_col) = top12 - z2
                endif

                ! handle case with zero elements in second vector
                if (top22_column_comm_rank == column_comm_rank) then
                    a(top22_row,top22_col) = -beta2
                    blas_rows = blas_rows - 1
                endif

                ! perform a special "daxpy" which uses 2 scalars and 2 vectors
                alpha = 1.0/tau2
                do row=1,blas_rows
                    a(row,top22_col) = (a(row,top22_col) - a(row,top22_col+1)*z2) * alpha
                enddo
         
                ! step 2.6.2: build real tau2, and add both tau values to their vector 
                tau2 = tau2/beta2
            endif ! fallback check
    
            work(1) = tau1
            work(2) = tau2
        endif ! END of householder vector generation */

        ! step 2.7: broadcast householder parts to other process columns
 
        if (local_rows > 0) then
            ! check if we are the column with the householder vectors,
            ! otherwise the processes copy their stuff just to lose it after broadcast */
            if (top11_row_comm_rank == row_comm_rank) then
                ! there is valuable householder vector data to broadcast
                work(3:local_rows+2) = a(1:local_rows,top22_col) ! second householder vector 
                work(local_rows+3:2*local_rows+2) = a(1:local_rows,top22_col+1) ! first householder vector
                
                ! sanitize vectors by inserting zeros and ones
                ! (but only on the process with the top elements)
                if (top11_column_comm_rank == column_comm_rank) then
                    work(2+top11_row) = 0.0
                    work(2+local_rows+top11_row) = 1.0
                endif
 
                if (top22_column_comm_rank == column_comm_rank) then
                    work(2+top22_row) = 1.0
                endif
            endif
        endif
  
      
        call mpi_bcast(work, 2+2*local_rows, mpi_real8, top11_row_comm_rank, mpi_comm_cols, mpierr)

        tau1 = work(1)
        tau2 = work(2)
  
        if (local_rows > 0) then
            vmr(1:local_rows,n_cols-1) = work(3:local_rows+2)
            vmr(1:local_rows,n_cols) = work(3+local_rows:2*local_rows+2)
        endif
 

        ! store tau values in future T matrix 
        tmat(n_cols-1,n_cols-1,istep)=tau2
        tmat(n_cols,n_cols,istep)=tau1

        ! if we produced the householder vectors, move to next 2 columns 
        if (top11_row_comm_rank == row_comm_rank) then
            local_columns=local_columns-2
        endif
 

        ! time for updating the remaining columns, if there are any 
        if (local_columns > 0) then
            ! fused approach
            call houseLeft_2rank_MPI_reversed(nrow, local_columns, nrow, a(1,blas_columns+1), lda, vmr(1,n_cols-1), ldv, work, mpi_comm_rows, nblk, work(3))

            ! normal approach
            !call houseLeft_reversed(vmr(1,n_cols), ldv, nrow, tau1, a(1,blas_columns+1), lda, local_columns, mpi_comm_rows, nblk, work(3))
            !call houseLeft_reversed(vmr(1,n_cols-1), ldv, nrow, tau2, a(1,blas_columns+1), lda, local_columns, mpi_comm_rows, nblk, work(3))
        endif
  
        ! advance to next two columns */
        n_cols = n_cols - 2
    end do
 
    ! handle last column the standard way
    if (n_cols .eq. 1) then
        ncol = istep*nbw + n_cols
        nrow = ncol - nbw
 
        if (nrow == 1) return ! nothing

        ! figure out size of sub matrix
        call matrixGridInfo(local_rows, nrow, nblk, column_processes, column_comm_rank)
        call matrixGridInfo(local_columns, ncol, nblk, row_processes, row_comm_rank)

        ! figure out process with "top" element
        call mapEntrytoProcess(column_processes, 0, nrow-1, 0, nblk, top11_row, top11_column_comm_rank)
        call mapEntrytoProcess(row_processes, 0, ncol-1, 0, nblk, top11_col, top11_row_comm_rank)

        if (top11_row_comm_rank == row_comm_rank) then
            ! calculate HV and tau 
            tau1 = houseGen_reversed(a(1,top11_col), nrow, mpi_comm_rows, nblk, work)
            ! store tau in seperate array
            work(1) = tau1
            work(2:local_rows+1) = a(1:local_rows,top11_col)
        endif

        ! broadcast HV and tau to other process columns
        call mpi_bcast(work, local_rows+1, mpi_real8, top11_row_comm_rank, mpi_comm_cols, mpierr)

        tau1 = work(1)
        tmat(n_cols,n_cols,istep)= tau1
        vmr(1:local_rows,n_cols) = work(2:local_rows+1)

        if (top11_column_comm_rank == column_comm_rank) then
            vmr(top11_row,n_cols) = 1.0
        endif
    endif
end subroutine qr_rank2_real

subroutine houseGen_scal(x,local_size,alpha)
    implicit none
    
    integer local_size
    real*8 alpha,x(*)

    call dscal(local_size, alpha, x, 1)
    !x(1:local_size) = x(1:local_size) * alpha
end subroutine houseGen_scal

real*8 function houseGen_dlapy2(x, y)
    implicit none

    real*8 x,y

    real*8 xabs
    real*8 yabs
    real*8 w 
    real*8 z

    xabs = abs(x)
    yabs = abs(y)

    w = max(xabs, yabs)
    z = min(xabs, yabs)

    if (z == 0.0) housegen_dlapy2 = w
    z = z / w
    houseGen_dlapy2 = w * sqrt(1.0+z*z)
end function houseGen_dlapy2

!
! Determines the Householder Vector which sets the 
! @param work: buffer, needs to be 2x process_count of size
! @param x pointer to part of build vector
! @param n size of total build vector
! @param b index from which all entries in the vector should be zeroed
! @param communicator MPI communicator over which the vector is spread
!
function houseGen_reversed(x, b, communicator, blk, work)
    use mpi

    implicit none
    real*8 houseGen_reversed

    real*8 x(*),work(*)
    integer b,blk,communicator

    real*8 beta
    real*8 local_norm 
    real*8 xnorm 
    real*8 alpha
    integer alpha_index
    integer alpha_rank

    real*8 tau
    integer my_grid_rank
    integer process_count
    real*8 temp
    real*8 safmin
    real*8 rsafmn
    integer knt
    integer b_entry
    integer i
    integer local_size
    integer mpierr
 
    ! what rank am i and what data parts do i have?
    call mpi_comm_rank(communicator, my_grid_rank,mpierr)

    ! how many processes are on this communicator 
    call mpi_comm_size(communicator, process_count,mpierr)
 
    ! max entries for this process
    call matrixGridInfo(local_size, b, blk, process_count, my_grid_rank)
 
    ! calculate norm
    local_norm = dnrm2(local_size, x, 1)
  
    alpha = 0.0
 
    ! determine alpha value (= first value in global vector part)
    call mapEntrytoProcess(process_count, 0, b-1, 0, blk, alpha_index, alpha_rank)
    if (my_grid_rank == alpha_rank) alpha = x(alpha_index)
 
    ! work buffer: first receive_bufffer, then send_buffer
    do i=1,process_count
        work(2*i-1) = local_norm
        work(2*i) = alpha
    enddo

    ! now exchange all local norms to all processes 
    call mpi_alltoall(work, 2, mpi_real8, work(2*process_count+1), 2, mpi_real8, communicator, mpierr)
 
    ! process with alpha rank provides alpha value
    alpha = work(2*process_count+1+2*alpha_rank+1)
 
    ! copy local norms into right position
    do i=1,process_count
        work(i) = work(2*process_count+2*i-1)
    enddo

    ! now calculate global norm
    xnorm = dnrm2(process_count, work, 1)

    ! use safe minimum as 0 indicator 
    !safmin = dlamch('S')

    ! General case 

    beta = xnorm
    if (alpha < 0.0) beta = -beta

    tau = (beta + alpha) / beta
    call houseGen_scal(x, local_size, 1.0/(alpha+beta))
    alpha = -beta

    if (alpha_rank == my_grid_rank) x(alpha_index) = alpha

    houseGen_reversed = tau
end function houseGen_reversed

! maps array slots of one process and old block size to a new process with new block size 
subroutine mapEntrytoProcess(processes, nb_old, entry_old, rank_old, nb_new, entry_new, rank_new)
    implicit none

    integer processes,nb_old,entry_old,rank_old,nb_new,entry_new,rank_new
    integer real_entry

    if ((nb_old == 0) .and. (nb_new == 0)) then
        entry_new = entry_old
        rank_new = 0
    endif

    if (nb_old == 0) then
        real_entry = entry_old
    else
        real_entry = rank_old * nb_old + (entry_old / nb_old) * (nb_old * processes) + modulo(entry_old,nb_old)
    endif

    if (nb_new == 0) then
        rank_new = 0
        entry_new = real_entry + 1 ! fortran adresses from 1 to n 
    else
        rank_new = modulo((real_entry / nb_new), processes)
        entry_new = (real_entry / (nb_new * processes)) * nb_new + modulo(real_entry,nb_new) + 1 ! fortran adresses from 1 to n
    endif
end subroutine mapEntrytoProcess

subroutine matrixGridInfo(my_entries, total_entries, nb, processes, rank)
    integer my_entries, total_entries, nb, processes, rank

    integer blocks
    integer entries
    integer blocks_extra

    if(nb == 0) then
        if(rank == 0) then
            my_entries = total_entries
        else
            my_entries = 0
        endif
    else
        blocks = total_entries / (nb * processes)
        entries = blocks * nb
        blocks_extra = modulo(total_entries, (nb * processes))

        if (blocks_extra > nb * rank) then
            blocks_extra = blocks_extra - nb * rank
            if (blocks_extra >= nb) then
                entries = entries + nb
            else
                entries = entries + blocks_extra
            endif
        endif

        my_entries = entries
    endif
end subroutine matrixGridInfo

subroutine houseLeft_reversed(hv, ldh, b, tau, A, lda, nb, communicator, blk, work)
    use mpi

    implicit none

    integer ldh, b, lda,nb, communicator, blk
    real*8 hv(*), tau, A(lda,*), work(*)

    integer column_size
    integer hv_size
    integer process_count
    integer my_rank
    integer x0_process_nr
    integer x0_index
    integer i
    integer incx, incy
    real*8 beta

    integer mpierr

    ! how many processes are there in my communicator
    call mpi_comm_size(communicator, process_count, mpierr)

    ! get my own rank in the communicator
    call mpi_comm_rank(communicator, my_rank, mpierr)

    ! determine minimal size of hv and A 
    call matrixGridInfo(column_size, b, blk, process_count, my_rank)
 
    if (column_size > 0) then
        call dgemv("T", column_size, nb, tau, A(1,1), lda, hv(1), 1, 0.0d0, work(1), 1)
    else
        work(1:nb) = 0
    endif

    !print *,'used hv:', b, my_rank, hv(1:column_size)


    call mpi_allreduce(work(1), work(nb+1), nb, mpi_real8, MPI_SUM, communicator, mpierr)
  
    !print *,'z real: ', b, my_rank, work(nb+1:2*nb)
 
    if (column_size > 0) then
        !do i=1,nb
        !    A(1:column_size,i) = A(1:column_size,i) - hv(1:column_size) * work(nb+i)
        !end do
        call dger(column_size, nb, -1.0d0, hv(1), 1, work(nb+1), 1, A(1,1), lda)
    endif
end subroutine houseLeft_reversed

! (there is a single zero gap between both householder vector parts )
! rows: number of rows of distributed matrix
! cols: columns per process 
! work: 2*cols+2 
subroutine houseLeft_2rank_MPI_reversed(rows, cols, b, a, lda, hv, ldh, tauvalues, communicator, blk, work)
    use mpi

    implicit none

    integer rows, cols, b, lda, ldh, communicator, blk
    real*8 a(lda,*), hv(ldh,*), tauvalues(*), work(*)

    integer local_rows, rank, processes, top11_process_rank, top11_row, top22_process_rank, top22_row
    real*8 hvdot, tau1, tau2, sum1, sum2, a_value, hv1, hv2
    integer coliter, row, rowoffset

    ! blas stuff 
    integer blas_rows
    real*8 beta

    integer mpierr

    tau1 = tauvalues(1)
    tau2 = tauvalues(2)

    call mpi_comm_size(communicator, processes, mpierr)
    call mpi_comm_rank(communicator, rank, mpierr)

    ! step 0: determine available data
    call mapEntrytoProcess(processes, 0, b-1, 0, blk, top11_row, top11_process_rank)
    call mapEntrytoProcess(processes, 0, b-2, 0, blk, top22_row, top22_process_rank) 
     
    call matrixGridInfo(blas_rows, b, blk, processes, rank)
 
    ! step 1: build partial sums 
    if (blas_rows == 0) then
        ! no more data for us, but initialize allreduce values to take part in collective operation */
        work(1:2*cols+2) = 0
    else
        ! build hvdot part
        hvdot = ddot(blas_rows,hv(1,1),1, hv(1,2),1)

        work(1) = hvdot
        work(2) = 0.0
 
        ! build first sum
        call dgemv("T", blas_rows, cols, tau1, a(1,1), lda, hv(1,2), 1, 0.0d0, work(3), 1)
  
        ! build second sum part - make use of zero entry in second householder vector 
        call dgemv("T", blas_rows, cols, tau2, a(1,1), lda, hv(1,1), 1, 0.0d0, work(3+cols), 1)
    endif

    ! step 2: allreduce sums and build final sum2
    call mpi_allreduce(work(1), work(2*cols+3), 2*cols+2, mpi_real8, MPI_SUM, communicator, mpierr)
  
    hvdot = work(2*cols+3)
    hvdot = hvdot * (-tau2)

    call daxpy(cols, hvdot, work(2*cols+5), 1, work(3*cols+5), 1)
 
    call matrixGridInfo(local_rows, b, blk, processes, rank)
    call matrixGridInfo(blas_rows, b-2, blk, processes, rank)

    ! step 3: update columns (if we have data to update) 
    if (local_rows > 0) then
        do coliter=1, cols
            ! handle first two elements if we got the top values 
            sum1 = work(2*cols+4+coliter)
            sum2 = work(3*cols+4+coliter)
 
            do row=1,blas_rows
                a_value = a(row,coliter)
                hv1 = hv(row,2)
                hv2 = hv(row,1)

                a_value = a_value - sum1 * hv1  
                a_value = a_value - sum2 * hv2
                a(row,coliter) = a_value
            enddo
 
            if (rank == top22_process_rank) then
                hv1 = hv(top22_row,2)
                a(top22_row,coliter) = a(top22_row,coliter) - sum1*hv1 - sum2
            endif
 
            if (rank == top11_process_rank) then
                a(top11_row,coliter) = a(top11_row,coliter) - sum1
            endif
        enddo
   endif
end subroutine houseLeft_2rank_MPI_reversed
#endif
end module blockedqr
