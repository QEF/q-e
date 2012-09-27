! ELPA2 -- 2-stage solver for ELPA
! 
! Copyright of the original code rests with the authors inside the ELPA
! consortium. The copyright of any additional modifications shall rest
! with their original authors, but shall adhere to the licensing terms
! distributed along with the original code in the file "COPYING".

module ELPA2

! Version 1.1.2, 2011-02-21

  USE ELPA1
  USE blockedQR

  implicit none

  PRIVATE ! By default, all routines contained are private

  ! The following routines are public:
#ifdef __ELPA

  public :: solve_evp_real_2stage
  public :: solve_evp_complex_2stage

  public :: bandred_real
  public :: tridiag_band_real
  public :: trans_ev_tridi_to_band_real
  public :: trans_ev_band_to_full_real

  public :: bandred_complex
  public :: tridiag_band_complex
  public :: trans_ev_tridi_to_band_complex
  public :: trans_ev_band_to_full_complex
  
  public :: band_band_real
  public :: divide_band
#endif
  
!-------------------------------------------------------------------------------  

  integer, public :: which_qr_decomposition = 1     ! defines, which QR-decomposition algorithm will be used
                                                    ! 0 for unblocked
                                                    ! 1 for rank-2

!-------------------------------------------------------------------------------

  ! The following array contains the Householder vectors of the
  ! transformation band -> tridiagonal.
  ! It is allocated and set in tridiag_band_real and used in
  ! trans_ev_tridi_to_band_real.
  ! It must be deallocated by the user after trans_ev_tridi_to_band_real!

  real*8, allocatable :: hh_trans_real(:,:)
  complex*16, allocatable :: hh_trans_complex(:,:)

!-------------------------------------------------------------------------------
#ifdef __ELPA
  include 'mpif.h'


!******
contains

subroutine solve_evp_real_2stage(na, nev, a, lda, ev, q, ldq, nblk, mpi_comm_rows, mpi_comm_cols, mpi_comm_all)

!-------------------------------------------------------------------------------
!  solve_evp_real_2stage: Solves the real eigenvalue problem with a 2 stage approach
!
!  Parameters
!
!  na          Order of matrix a
!
!  nev         Number of eigenvalues needed
!
!  a(lda,*)    Distributed matrix for which eigenvalues are to be computed.
!              Distribution is like in Scalapack.
!              The full matrix must be set (not only one half like in scalapack).
!              Destroyed on exit (upper and lower half).
!
!  lda         Leading dimension of a
!
!  ev(na)      On output: eigenvalues of a, every processor gets the complete set
!
!  q(ldq,*)    On output: Eigenvectors of a
!              Distribution is like in Scalapack.
!              Must be always dimensioned to the full size (corresponding to (na,na))
!              even if only a part of the eigenvalues is needed.
!
!  ldq         Leading dimension of q
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!  mpi_comm_all
!              MPI-Communicator for the total processor set
!
!-------------------------------------------------------------------------------

   implicit none

   integer, intent(in) :: na, nev, lda, ldq, nblk, mpi_comm_rows, mpi_comm_cols, mpi_comm_all
   real*8, intent(inout) :: a(lda,*), ev(na), q(ldq,*)

   integer my_pe, n_pes, my_prow, my_pcol, np_rows, np_cols, mpierr
   integer nbw, num_blocks
   real*8, allocatable :: tmat(:,:,:), e(:)
   real*8 ttt0, ttt1, ttts

   call mpi_comm_rank(mpi_comm_all,my_pe,mpierr)
   call mpi_comm_size(mpi_comm_all,n_pes,mpierr)

   call mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   call mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   call mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   call mpi_comm_size(mpi_comm_cols,np_cols,mpierr)

   ! Choose bandwidth, must be a multiple of nblk, set to a value >= 32

   nbw = (31/nblk+1)*nblk

   num_blocks = (na-1)/nbw + 1

   allocate(tmat(nbw,nbw,num_blocks))

   ! Reduction full -> band

   ttt0 = MPI_Wtime()
   ttts = ttt0
   call bandred_real(na, a, lda, nblk, nbw, mpi_comm_rows, mpi_comm_cols, tmat)
   ttt1 = MPI_Wtime()
   if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) &
      print 1,'Time bandred_real               :',ttt1-ttt0

   ! Reduction band -> tridiagonal

   allocate(e(na))

   ttt0 = MPI_Wtime()
   call tridiag_band_real(na, nbw, nblk, a, lda, ev, e, mpi_comm_rows, mpi_comm_cols, mpi_comm_all)
   ttt1 = MPI_Wtime()
   if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) &
      print 1,'Time tridiag_band_real          :',ttt1-ttt0

   call mpi_bcast(ev,na,MPI_REAL8,0,mpi_comm_all,mpierr)
   call mpi_bcast(e,na,MPI_REAL8,0,mpi_comm_all,mpierr)

   ttt1 = MPI_Wtime()
   time_evp_fwd = ttt1-ttts

   ! Solve tridiagonal system

   ttt0 = MPI_Wtime()
   call solve_tridi(na, nev, ev, e, q, ldq, nblk, mpi_comm_rows, mpi_comm_cols)
   ttt1 = MPI_Wtime()
   if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) &
      print 1,'Time solve_tridi                :',ttt1-ttt0
   time_evp_solve = ttt1-ttt0
   ttts = ttt1

   deallocate(e)

   ! Backtransform stage 1

   ttt0 = MPI_Wtime()
   call trans_ev_tridi_to_band_real(na, nev, nblk, nbw, q, ldq, mpi_comm_rows, mpi_comm_cols)
   ttt1 = MPI_Wtime()
   if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) &
      print 1,'Time trans_ev_tridi_to_band_real:',ttt1-ttt0

   ! We can now deallocate the stored householder vectors
   deallocate(hh_trans_real)

   ! Backtransform stage 2

   ttt0 = MPI_Wtime()
   call trans_ev_band_to_full_real(na, nev, nblk, nbw, a, lda, tmat, q, ldq, mpi_comm_rows, mpi_comm_cols)
   ttt1 = MPI_Wtime()
   if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) &
      print 1,'Time trans_ev_band_to_full_real :',ttt1-ttt0
   time_evp_back = ttt1-ttts

   deallocate(tmat)

1  format(a,f10.3)

end subroutine solve_evp_real_2stage

!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------

subroutine solve_evp_complex_2stage(na, nev, a, lda, ev, q, ldq, nblk, mpi_comm_rows, mpi_comm_cols, mpi_comm_all)

!-------------------------------------------------------------------------------
!  solve_evp_complex_2stage: Solves the complex eigenvalue problem with a 2 stage approach
!
!  Parameters
!
!  na          Order of matrix a
!
!  nev         Number of eigenvalues needed
!
!  a(lda,*)    Distributed matrix for which eigenvalues are to be computed.
!              Distribution is like in Scalapack.
!              The full matrix must be set (not only one half like in scalapack).
!              Destroyed on exit (upper and lower half).
!
!  lda         Leading dimension of a
!
!  ev(na)      On output: eigenvalues of a, every processor gets the complete set
!
!  q(ldq,*)    On output: Eigenvectors of a
!              Distribution is like in Scalapack.
!              Must be always dimensioned to the full size (corresponding to (na,na))
!              even if only a part of the eigenvalues is needed.
!
!  ldq         Leading dimension of q
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!  mpi_comm_all
!              MPI-Communicator for the total processor set
!
!-------------------------------------------------------------------------------

   implicit none

   integer, intent(in) :: na, nev, lda, ldq, nblk, mpi_comm_rows, mpi_comm_cols, mpi_comm_all
   complex*16, intent(inout) :: a(lda,*), q(ldq,*)
   real*8, intent(inout) :: ev(na)

   integer my_prow, my_pcol, np_rows, np_cols, mpierr
   integer l_cols, l_rows, l_cols_nev, nbw, num_blocks
   complex*16, allocatable :: tmat(:,:,:)
   real*8, allocatable :: q_real(:,:), e(:)
   real*8 ttt0, ttt1, ttts

   call mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   call mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   call mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   call mpi_comm_size(mpi_comm_cols,np_cols,mpierr)

   ! Choose bandwidth, must be a multiple of nblk, set to a value >= 32

   nbw = (31/nblk+1)*nblk

   num_blocks = (na-1)/nbw + 1

   allocate(tmat(nbw,nbw,num_blocks))

   ! Reduction full -> band

   ttt0 = MPI_Wtime()
   ttts = ttt0
   call bandred_complex(na, a, lda, nblk, nbw, mpi_comm_rows, mpi_comm_cols, tmat)
   ttt1 = MPI_Wtime()
   if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) &
      print 1,'Time bandred_complex               :',ttt1-ttt0

   ! Reduction band -> tridiagonal

   allocate(e(na))

   ttt0 = MPI_Wtime()
   call tridiag_band_complex(na, nbw, nblk, a, lda, ev, e, mpi_comm_rows, mpi_comm_cols, mpi_comm_all)
   ttt1 = MPI_Wtime()
   if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) &
      print 1,'Time tridiag_band_complex          :',ttt1-ttt0

   call mpi_bcast(ev,na,MPI_REAL8,0,mpi_comm_all,mpierr)
   call mpi_bcast(e,na,MPI_REAL8,0,mpi_comm_all,mpierr)

   ttt1 = MPI_Wtime()
   time_evp_fwd = ttt1-ttts

   l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a and q
   l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local columns of q
   l_cols_nev = local_index(nev, my_pcol, np_cols, nblk, -1) ! Local columns corresponding to nev

   allocate(q_real(l_rows,l_cols))

   ! Solve tridiagonal system

   ttt0 = MPI_Wtime()
   call solve_tridi(na, nev, ev, e, q_real, ubound(q_real,1), nblk, mpi_comm_rows, mpi_comm_cols)
   ttt1 = MPI_Wtime()
   if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times)  &
      print 1,'Time solve_tridi                   :',ttt1-ttt0
   time_evp_solve = ttt1-ttt0
   ttts = ttt1

   q(1:l_rows,1:l_cols_nev) = q_real(1:l_rows,1:l_cols_nev)

   deallocate(e, q_real)

   ! Backtransform stage 1

   ttt0 = MPI_Wtime()
   call trans_ev_tridi_to_band_complex(na, nev, nblk, nbw, q, ldq, mpi_comm_rows, mpi_comm_cols)
   ttt1 = MPI_Wtime()
   if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) &
      print 1,'Time trans_ev_tridi_to_band_complex:',ttt1-ttt0

   ! We can now deallocate the stored householder vectors
   deallocate(hh_trans_complex)

   ! Backtransform stage 2

   ttt0 = MPI_Wtime()
   call trans_ev_band_to_full_complex(na, nev, nblk, nbw, a, lda, tmat, q, ldq, mpi_comm_rows, mpi_comm_cols)
   ttt1 = MPI_Wtime()
   if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) &
      print 1,'Time trans_ev_band_to_full_complex :',ttt1-ttt0
   time_evp_back = ttt1-ttts

   deallocate(tmat)

1  format(a,f10.3)

end subroutine solve_evp_complex_2stage

!-------------------------------------------------------------------------------

subroutine bandred_real(na, a, lda, nblk, nbw, mpi_comm_rows, mpi_comm_cols, tmat)

!-------------------------------------------------------------------------------
!  bandred_real: Reduces a distributed symmetric matrix to band form
!
!  Parameters
!
!  na          Order of matrix
!
!  a(lda,*)    Distributed matrix which should be reduced.
!              Distribution is like in Scalapack.
!              Opposed to Scalapack, a(:,:) must be set completely (upper and lower half)
!              a(:,:) is overwritten on exit with the band and the Householder vectors
!              in the upper half.
!
!  lda         Leading dimension of a
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  nbw         semi bandwith of output matrix
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!
!  tmat(nbw,nbw,num_blocks)    where num_blocks = (na-1)/nbw + 1
!              Factors for the Householder vectors (returned), needed for back transformation
!
!-------------------------------------------------------------------------------

   implicit none

   integer na, lda, nblk, nbw, mpi_comm_rows, mpi_comm_cols
   real*8 a(lda,*), tmat(nbw,nbw,*)

   integer my_prow, my_pcol, np_rows, np_cols, mpierr
   integer l_cols, l_rows
   integer i, j, lcs, lce, lre, lc, lr, cur_pcol, n_cols, nrow
   integer istep, ncol, lch, lcx, nlc
   integer tile_size, l_rows_tile, l_cols_tile
   integer work_size
   
   real*8 eps

   real*8 vnorm2, xf, aux1(nbw), aux2(nbw), vrl, tau, vav(nbw,nbw)

   real*8, allocatable:: tmp(:,:), vr(:), vmr(:,:), umc(:,:)
   real*8, allocatable:: work(:) ! needed for blocked QR


   integer pcol, prow
   pcol(i) = MOD((i-1)/nblk,np_cols) !Processor col for global col number
   prow(i) = MOD((i-1)/nblk,np_rows) !Processor row for global row number


   call mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   call mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   call mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   call mpi_comm_size(mpi_comm_cols,np_cols,mpierr)

   ! Semibandwith nbw must be a multiple of blocksize nblk

   if(mod(nbw,nblk)/=0) then
      if(my_prow==0 .and. my_pcol==0) then
         print *,'ERROR: nbw=',nbw,', nblk=',nblk
         print *,'ELPA2 works only for nbw==n*nblk'
         call mpi_abort(mpi_comm_world,0,mpierr)
      endif
   endif

   ! Matrix is split into tiles; work is done only for tiles on the diagonal or above

   tile_size = nblk*least_common_multiple(np_rows,np_cols) ! minimum global tile size
   tile_size = ((128*max(np_rows,np_cols)-1)/tile_size+1)*tile_size ! make local tiles at least 128 wide

   l_rows_tile = tile_size/np_rows ! local rows of a tile
   l_cols_tile = tile_size/np_cols ! local cols of a tile
 
   if (which_qr_decomposition == 1) then
          l_rows = local_index(na, my_prow, np_rows, nblk, -1)
          work_size = max(4*np_rows,2*nbw)
          work_size = max(l_rows+1,work_size)
          work_size = max(16, work_size)
          work_size = max(2*(l_rows+1),work_size)
          work_size = max(2+4*(nbw+1),work_size)
          allocate(work(work_size))
          work = 0
          eps = 1.0d0
   endif

   do istep = (na-1)/nbw, 1, -1

      n_cols = MIN(na,(istep+1)*nbw) - istep*nbw ! Number of columns in current step

      ! Number of local columns/rows of remaining matrix
      l_cols = local_index(istep*nbw, my_pcol, np_cols, nblk, -1)
      l_rows = local_index(istep*nbw, my_prow, np_rows, nblk, -1)

      ! Allocate vmr and umc to their exact sizes so that they can be used in bcasts and reduces

      allocate(vmr(max(l_rows,1),2*n_cols))
      allocate(umc(max(l_cols,1),2*n_cols))

      allocate(vr(l_rows+1))

      vmr(1:l_rows,1:n_cols) = 0.
      vr(:) = 0
      tmat(:,:,istep) = 0

      ! Reduce current block to lower triangular form
      if (which_qr_decomposition == 1) then
          call qr_rank2_real(a, lda, vmr, max(l_rows,1), tmat, nbw, istep, n_cols, nblk, mpi_comm_rows, mpi_comm_cols, work, eps)
      else

      do lc = n_cols, 1, -1

         ncol = istep*nbw + lc ! absolute column number of householder vector
         nrow = ncol - nbw ! Absolute number of pivot row

         lr  = local_index(nrow, my_prow, np_rows, nblk, -1) ! current row length
         lch = local_index(ncol, my_pcol, np_cols, nblk, -1) ! HV local column number

         tau = 0

         if(nrow == 1) exit ! Nothing to do

         cur_pcol = pcol(ncol) ! Processor column owning current block

         if(my_pcol==cur_pcol) then

            ! Get vector to be transformed; distribute last element and norm of
            ! remaining elements to all procs in current column

            vr(1:lr) = a(1:lr,lch) ! vector to be transformed

            if(my_prow==prow(nrow)) then
               aux1(1) = dot_product(vr(1:lr-1),vr(1:lr-1))
               aux1(2) = vr(lr)
            else
               aux1(1) = dot_product(vr(1:lr),vr(1:lr))
               aux1(2) = 0.
            endif

            call mpi_allreduce(aux1,aux2,2,MPI_REAL8,MPI_SUM,mpi_comm_rows,mpierr)

            vnorm2 = aux2(1)
            vrl    = aux2(2)

            ! Householder transformation

            call hh_transform_real(vrl, vnorm2, xf, tau)

            ! Scale vr and store Householder vector for back transformation

            vr(1:lr) = vr(1:lr) * xf
            if(my_prow==prow(nrow)) then
               a(1:lr-1,lch) = vr(1:lr-1)
               a(lr,lch) = vrl
               vr(lr) = 1.
            else
               a(1:lr,lch) = vr(1:lr)
            endif

         endif

         ! Broadcast Householder vector and tau along columns

         vr(lr+1) = tau
         call MPI_Bcast(vr,lr+1,MPI_REAL8,cur_pcol,mpi_comm_cols,mpierr)
         vmr(1:lr,lc) = vr(1:lr)
         tau = vr(lr+1)
         tmat(lc,lc,istep) = tau ! Store tau in diagonal of tmat

         ! Transform remaining columns in current block with Householder vector

         ! Local dot product

         aux1 = 0

         nlc = 0 ! number of local columns
         do j=1,lc-1
            lcx = local_index(istep*nbw+j, my_pcol, np_cols, nblk, 0)
            if(lcx>0) then
               nlc = nlc+1
               if(lr>0) aux1(nlc) = dot_product(vr(1:lr),a(1:lr,lcx))
            endif
         enddo

         ! Get global dot products
         if(nlc>0) call mpi_allreduce(aux1,aux2,nlc,MPI_REAL8,MPI_SUM,mpi_comm_rows,mpierr)

         ! Transform

         nlc = 0
         do j=1,lc-1
            lcx = local_index(istep*nbw+j, my_pcol, np_cols, nblk, 0)
            if(lcx>0) then
               nlc = nlc+1
               a(1:lr,lcx) = a(1:lr,lcx) - tau*aux2(nlc)*vr(1:lr)
            endif
         enddo

      enddo

      endif

      ! Calculate scalar products of stored Householder vectors.
      ! This can be done in different ways, we use dsyrk

      vav = 0
      if(l_rows>0) &
         call dsyrk('U','T',n_cols,l_rows,1.d0,vmr,ubound(vmr,1),0.d0,vav,ubound(vav,1))
      call symm_matrix_allreduce(n_cols,vav,ubound(vav,1),mpi_comm_rows)

      ! Calculate triangular matrix T for block Householder Transformation

      do lc=n_cols,1,-1
         tau = tmat(lc,lc,istep)
         if(lc<n_cols) then
            call dtrmv('U','T','N',n_cols-lc,tmat(lc+1,lc+1,istep),ubound(tmat,1),vav(lc+1,lc),1)
            tmat(lc,lc+1:n_cols,istep) = -tau * vav(lc+1:n_cols,lc)
         endif
      enddo

      ! Transpose vmr -> vmc (stored in umc, second half)

      call elpa_transpose_vectors  (vmr, ubound(vmr,1), mpi_comm_rows, &
                                    umc(1,n_cols+1), ubound(umc,1), mpi_comm_cols, &
                                    1, istep*nbw, n_cols, nblk)

      ! Calculate umc = A**T * vmr
      ! Note that the distributed A has to be transposed
      ! Opposed to direct tridiagonalization there is no need to use the cache locality
      ! of the tiles, so we can use strips of the matrix

      umc(1:l_cols,1:n_cols) = 0.d0
      vmr(1:l_rows,n_cols+1:2*n_cols) = 0
      if(l_cols>0 .and. l_rows>0) then
         do i=0,(istep*nbw-1)/tile_size

            lcs = i*l_cols_tile+1
            lce = min(l_cols,(i+1)*l_cols_tile)
            if(lce<lcs) cycle

            lre = min(l_rows,(i+1)*l_rows_tile)
            call DGEMM('T','N',lce-lcs+1,n_cols,lre,1.d0,a(1,lcs),ubound(a,1), &
                       vmr,ubound(vmr,1),1.d0,umc(lcs,1),ubound(umc,1))

            if(i==0) cycle
            lre = min(l_rows,i*l_rows_tile)
            call DGEMM('N','N',lre,n_cols,lce-lcs+1,1.d0,a(1,lcs),lda, &
                       umc(lcs,n_cols+1),ubound(umc,1),1.d0,vmr(1,n_cols+1),ubound(vmr,1))
         enddo
      endif

      ! Sum up all ur(:) parts along rows and add them to the uc(:) parts
      ! on the processors containing the diagonal
      ! This is only necessary if ur has been calculated, i.e. if the
      ! global tile size is smaller than the global remaining matrix

      if(tile_size < istep*nbw) then
         call elpa_reduce_add_vectors  (vmr(1,n_cols+1),ubound(vmr,1),mpi_comm_rows, &
                                        umc, ubound(umc,1), mpi_comm_cols, &
                                        istep*nbw, n_cols, nblk)
      endif

      if(l_cols>0) then
         allocate(tmp(l_cols,n_cols))
         call mpi_allreduce(umc,tmp,l_cols*n_cols,MPI_REAL8,MPI_SUM,mpi_comm_rows,mpierr)
         umc(1:l_cols,1:n_cols) = tmp(1:l_cols,1:n_cols)
         deallocate(tmp)
      endif

      ! U = U * Tmat**T

      call dtrmm('Right','Upper','Trans','Nonunit',l_cols,n_cols,1.d0,tmat(1,1,istep),ubound(tmat,1),umc,ubound(umc,1))

      ! VAV = Tmat * V**T * A * V * Tmat**T = (U*Tmat**T)**T * V * Tmat**T

      call dgemm('T','N',n_cols,n_cols,l_cols,1.d0,umc,ubound(umc,1),umc(1,n_cols+1),ubound(umc,1),0.d0,vav,ubound(vav,1))
      call dtrmm('Right','Upper','Trans','Nonunit',n_cols,n_cols,1.d0,tmat(1,1,istep),ubound(tmat,1),vav,ubound(vav,1))

      call symm_matrix_allreduce(n_cols,vav,ubound(vav,1),mpi_comm_cols)

      ! U = U - 0.5 * V * VAV
      call dgemm('N','N',l_cols,n_cols,n_cols,-0.5d0,umc(1,n_cols+1),ubound(umc,1),vav,ubound(vav,1),1.d0,umc,ubound(umc,1))

      ! Transpose umc -> umr (stored in vmr, second half)

       call elpa_transpose_vectors  (umc, ubound(umc,1), mpi_comm_cols, &
                                     vmr(1,n_cols+1), ubound(vmr,1), mpi_comm_rows, &
                                     1, istep*nbw, n_cols, nblk)

      ! A = A - V*U**T - U*V**T

      do i=0,(istep*nbw-1)/tile_size
         lcs = i*l_cols_tile+1
         lce = min(l_cols,(i+1)*l_cols_tile)
         lre = min(l_rows,(i+1)*l_rows_tile)
         if(lce<lcs .or. lre<1) cycle
         call dgemm('N','T',lre,lce-lcs+1,2*n_cols,-1.d0, &
                    vmr,ubound(vmr,1),umc(lcs,1),ubound(umc,1), &
                    1.d0,a(1,lcs),lda)
      enddo

      deallocate(vmr, umc, vr)

   enddo
 
   if (which_qr_decomposition == 1) then
          deallocate(work)
   endif

end subroutine bandred_real

!-------------------------------------------------------------------------------

subroutine symm_matrix_allreduce(n,a,lda,comm)

!-------------------------------------------------------------------------------
!  symm_matrix_allreduce: Does an mpi_allreduce for a symmetric matrix A.
!  On entry, only the upper half of A needs to be set
!  On exit, the complete matrix is set
!-------------------------------------------------------------------------------

   implicit none
   integer n, lda, comm
   real*8 a(lda,*)

   integer i, nc, mpierr
   real*8 h1(n*n), h2(n*n)

   nc = 0
   do i=1,n
      h1(nc+1:nc+i) = a(1:i,i)
      nc = nc+i
   enddo

   call mpi_allreduce(h1,h2,nc,MPI_REAL8,MPI_SUM,comm,mpierr)

   nc = 0
   do i=1,n
      a(1:i,i) = h2(nc+1:nc+i)
      a(i,1:i-1) = a(1:i-1,i)
      nc = nc+i
   enddo

end subroutine symm_matrix_allreduce

!-------------------------------------------------------------------------------

subroutine trans_ev_band_to_full_real(na, nqc, nblk, nbw, a, lda, tmat, q, ldq, mpi_comm_rows, mpi_comm_cols)

!-------------------------------------------------------------------------------
!  trans_ev_band_to_full_real:
!  Transforms the eigenvectors of a band matrix back to the eigenvectors of the original matrix
!
!  Parameters
!
!  na          Order of matrix a, number of rows of matrix q
!
!  nqc         Number of columns of matrix q
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  nbw         semi bandwith
!
!  a(lda,*)    Matrix containing the Householder vectors (i.e. matrix a after bandred_real)
!              Distribution is like in Scalapack.
!
!  lda         Leading dimension of a
!
!  tmat(nbw,nbw,.) Factors returned by bandred_real
!
!  q           On input: Eigenvectors of band matrix
!              On output: Transformed eigenvectors
!              Distribution is like in Scalapack.
!
!  ldq         Leading dimension of q
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!
!-------------------------------------------------------------------------------

   implicit none

   integer na, nqc, lda, ldq, nblk, nbw, mpi_comm_rows, mpi_comm_cols
   real*8 a(lda,*), q(ldq,*), tmat(nbw, nbw, *)

   integer my_prow, my_pcol, np_rows, np_cols, mpierr
   integer max_blocks_row, max_blocks_col, max_local_rows, max_local_cols
   integer l_cols, l_rows, l_colh, n_cols
   integer istep, lc, ncol, nrow, nb, ns

   real*8, allocatable:: tmp1(:), tmp2(:), hvb(:), hvm(:,:)

   integer pcol, prow, i
   pcol(i) = MOD((i-1)/nblk,np_cols) !Processor col for global col number
   prow(i) = MOD((i-1)/nblk,np_rows) !Processor row for global row number


   call mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   call mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   call mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   call mpi_comm_size(mpi_comm_cols,np_cols,mpierr)


   max_blocks_row = ((na -1)/nblk)/np_rows + 1  ! Rows of A
   max_blocks_col = ((nqc-1)/nblk)/np_cols + 1  ! Columns of q!

   max_local_rows = max_blocks_row*nblk
   max_local_cols = max_blocks_col*nblk

   allocate(tmp1(max_local_cols*nbw))
   allocate(tmp2(max_local_cols*nbw))
   allocate(hvb(max_local_rows*nbw))
   allocate(hvm(max_local_rows,nbw))

   hvm = 0   ! Must be set to 0 !!!
   hvb = 0   ! Safety only

   l_cols = local_index(nqc, my_pcol, np_cols, nblk, -1) ! Local columns of q

   do istep=1,(na-1)/nbw

      n_cols = MIN(na,(istep+1)*nbw) - istep*nbw ! Number of columns in current step

      ! Broadcast all Householder vectors for current step compressed in hvb

      nb = 0
      ns = 0

      do lc = 1, n_cols
         ncol = istep*nbw + lc ! absolute column number of householder vector
         nrow = ncol - nbw ! absolute number of pivot row

         l_rows = local_index(nrow-1, my_prow, np_rows, nblk, -1) ! row length for bcast
         l_colh = local_index(ncol  , my_pcol, np_cols, nblk, -1) ! HV local column number

         if(my_pcol==pcol(ncol)) hvb(nb+1:nb+l_rows) = a(1:l_rows,l_colh)

         nb = nb+l_rows

         if(lc==n_cols .or. mod(ncol,nblk)==0) then
            call MPI_Bcast(hvb(ns+1),nb-ns,MPI_REAL8,pcol(ncol),mpi_comm_cols,mpierr)
            ns = nb
         endif
      enddo

      ! Expand compressed Householder vectors into matrix hvm

      nb = 0
      do lc = 1, n_cols
         nrow = (istep-1)*nbw+lc ! absolute number of pivot row
         l_rows = local_index(nrow-1, my_prow, np_rows, nblk, -1) ! row length for bcast

         hvm(1:l_rows,lc) = hvb(nb+1:nb+l_rows)
         if(my_prow==prow(nrow)) hvm(l_rows+1,lc) = 1.

         nb = nb+l_rows
      enddo

      l_rows = local_index(MIN(na,(istep+1)*nbw), my_prow, np_rows, nblk, -1)

      ! Q = Q - V * T**T * V**T * Q

      if(l_rows>0) then
         call dgemm('T','N',n_cols,l_cols,l_rows,1.d0,hvm,ubound(hvm,1), &
                    q,ldq,0.d0,tmp1,n_cols)
      else
         tmp1(1:l_cols*n_cols) = 0
      endif
      call mpi_allreduce(tmp1,tmp2,n_cols*l_cols,MPI_REAL8,MPI_SUM,mpi_comm_rows,mpierr)
      if(l_rows>0) then
         call dtrmm('L','U','T','N',n_cols,l_cols,1.0d0,tmat(1,1,istep),ubound(tmat,1),tmp2,n_cols)
         call dgemm('N','N',l_rows,l_cols,n_cols,-1.d0,hvm,ubound(hvm,1), &
                    tmp2,n_cols,1.d0,q,ldq)
      endif

   enddo

   deallocate(tmp1, tmp2, hvb, hvm)


end subroutine trans_ev_band_to_full_real

! --------------------------------------------------------------------------------------------------

subroutine tridiag_band_real(na, nb, nblk, a, lda, d, e, mpi_comm_rows, mpi_comm_cols, mpi_comm)

!-------------------------------------------------------------------------------
! tridiag_band_real:
! Reduces a real symmetric band matrix to tridiagonal form
!
!  na          Order of matrix a
!
!  nb          Semi bandwith
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  a(lda,*)    Distributed system matrix reduced to banded form in the upper diagonal
!
!  lda         Leading dimension of a
!
!  d(na)       Diagonal of tridiagonal matrix, set only on PE 0 (output)
!
!  e(na)       Subdiagonal of tridiagonal matrix, set only on PE 0 (output)
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!  mpi_comm
!              MPI-Communicator for the total processor set
!-------------------------------------------------------------------------------

   implicit none

   integer, intent(in) ::  na, nb, nblk, lda, mpi_comm_rows, mpi_comm_cols, mpi_comm
   real*8, intent(in)  :: a(lda,*)
   real*8, intent(out) :: d(na), e(na) ! set only on PE 0


   real*8 vnorm2, hv(nb), tau, x, h(nb), ab_s(1+nb), hv_s(nb), hv_new(nb), tau_new, hf
   real*8 hd(nb), hs(nb)

   integer i, j, n, nc, nr, ns, ne, istep, iblk, nblocks_total, nblocks, nt
   integer my_pe, n_pes, mpierr
   integer my_prow, np_rows, my_pcol, np_cols
   integer ireq_ab, ireq_hv
   integer na_s, nx, num_hh_vecs, num_chunks, local_size, max_blk_size, n_off
   integer, allocatable :: ireq_hhr(:), ireq_hhs(:), global_id(:,:), hh_cnt(:), hh_dst(:)
   integer, allocatable :: limits(:), snd_limits(:,:)
   integer, allocatable :: block_limits(:)
   real*8, allocatable :: ab(:,:), hh_gath(:,:,:), hh_send(:,:,:)
   ! dummies for calling redist_band
   complex*16 :: c_a(1,1), c_ab(1,1)


   call mpi_comm_rank(mpi_comm,my_pe,mpierr)
   call mpi_comm_size(mpi_comm,n_pes,mpierr)

   call mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   call mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   call mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   call mpi_comm_size(mpi_comm_cols,np_cols,mpierr)

   ! Get global_id mapping 2D procssor coordinates to global id

   allocate(global_id(0:np_rows-1,0:np_cols-1))
   global_id(:,:) = 0
   global_id(my_prow, my_pcol) = my_pe

   call mpi_allreduce(mpi_in_place, global_id, np_rows*np_cols, mpi_integer, mpi_sum, mpi_comm, mpierr)


   ! Total number of blocks in the band:

   nblocks_total = (na-1)/nb + 1

   ! Set work distribution

   allocate(block_limits(0:n_pes))
   call divide_band(nblocks_total, n_pes, block_limits)

   ! nblocks: the number of blocks for my task
   nblocks = block_limits(my_pe+1) - block_limits(my_pe)

   ! allocate the part of the band matrix which is needed by this PE
   ! The size is 1 block larger than needed to avoid extensive shifts
   allocate(ab(2*nb,(nblocks+1)*nb))
   ab = 0 ! needed for lower half, the extra block should also be set to 0 for safety

   ! n_off: Offset of ab within band
   n_off = block_limits(my_pe)*nb

   ! Redistribute band in a to ab
   call redist_band(.true., a, c_a, lda, na, nblk, nb, mpi_comm_rows, mpi_comm_cols, mpi_comm, ab, c_ab)

   ! Calculate the workload for each sweep in the back transformation
   ! and the space requirements to hold the HH vectors

   allocate(limits(0:np_rows))
   call determine_workload(na, nb, np_rows, limits)
   max_blk_size = maxval(limits(1:np_rows) - limits(0:np_rows-1))

   num_hh_vecs = 0
   num_chunks  = 0
   nx = na
   do n = 1, nblocks_total
      call determine_workload(nx, nb, np_rows, limits)
      local_size = limits(my_prow+1) - limits(my_prow)
      ! add to number of householder vectors
      ! please note: for nx==1 the one and only HH vector is 0 and is neither calculated nor send below!
      if(mod(n-1,np_cols) == my_pcol .and. local_size>0 .and. nx>1) then
         num_hh_vecs = num_hh_vecs + local_size
         num_chunks  = num_chunks+1
      endif
      nx = nx - nb
   enddo

   ! Allocate space for HH vectors

   allocate(hh_trans_real(nb,num_hh_vecs))

   ! Allocate and init MPI requests

   allocate(ireq_hhr(num_chunks)) ! Recv requests
   allocate(ireq_hhs(nblocks))    ! Send requests

   num_hh_vecs = 0
   num_chunks  = 0
   nx = na
   nt = 0
   do n = 1, nblocks_total
      call determine_workload(nx, nb, np_rows, limits)
      local_size = limits(my_prow+1) - limits(my_prow)
      if(mod(n-1,np_cols) == my_pcol .and. local_size>0 .and. nx>1) then
         num_chunks  = num_chunks+1
         call mpi_irecv(hh_trans_real(1,num_hh_vecs+1), nb*local_size, mpi_real8, nt, &
                        10+n-block_limits(nt), mpi_comm, ireq_hhr(num_chunks), mpierr)
         num_hh_vecs = num_hh_vecs + local_size
      endif
      nx = nx - nb
      if(n == block_limits(nt+1)) then
         nt = nt + 1
      endif
   enddo

   ireq_hhs(:) = MPI_REQUEST_NULL

   ! Buffers for gathering/sending the HH vectors

   allocate(hh_gath(nb,max_blk_size,nblocks)) ! gathers HH vectors
   allocate(hh_send(nb,max_blk_size,nblocks)) ! send buffer for HH vectors
   hh_gath(:,:,:) = 0
   hh_send(:,:,:) = 0

   ! Some counters

   allocate(hh_cnt(nblocks))
   allocate(hh_dst(nblocks))

   hh_cnt(:) = 1 ! The first transfomation vector is always 0 and not calculated at all
   hh_dst(:) = 0 ! PE number for receive

   ireq_ab = MPI_REQUEST_NULL
   ireq_hv = MPI_REQUEST_NULL

   ! Limits for sending

   allocate(snd_limits(0:np_rows,nblocks))

   do iblk=1,nblocks
      call determine_workload(na-(iblk+block_limits(my_pe)-1)*nb, nb, np_rows, snd_limits(:,iblk))
   enddo

   ! ---------------------------------------------------------------------------
   ! Start of calculations

   na_s = block_limits(my_pe)*nb + 1

   if(my_pe>0 .and. na_s<=na) then
      ! send first column to previous PE
      ! Only the PE owning the diagonal does that (sending 1 element of the subdiagonal block also)
      ab_s(1:nb+1) = ab(1:nb+1,na_s-n_off)
      call mpi_isend(ab_s,nb+1,mpi_real8,my_pe-1,1,mpi_comm,ireq_ab,mpierr)
   endif

   do istep=1,na-1

      if(my_pe==0) then
         n = MIN(na-na_s,nb) ! number of rows to be reduced
         hv(:) = 0
         tau = 0
         ! The last step (istep=na-1) is only needed for sending the last HH vectors.
         ! We don't want the sign of the last element flipped (analogous to the other sweeps)
         if(istep < na-1) then
            ! Transform first column of remaining matrix
            vnorm2 = sum(ab(3:n+1,na_s-n_off)**2)
            call hh_transform_real(ab(2,na_s-n_off),vnorm2,hf,tau)
            hv(1) = 1
            hv(2:n) = ab(3:n+1,na_s-n_off)*hf
         endif
         d(istep) = ab(1,na_s-n_off)
         e(istep) = ab(2,na_s-n_off)
         if(istep == na-1) then
            d(na) = ab(1,na_s+1-n_off)
            e(na) = 0
         endif
      else
         if(na>na_s) then
            ! Receive Householder vector from previous task, from PE owning subdiagonal
            call mpi_recv(hv,nb,mpi_real8,my_pe-1,2,mpi_comm,MPI_STATUS_IGNORE,mpierr)
            tau = hv(1)
            hv(1) = 1.
         endif
      endif

      na_s = na_s+1
      if(na_s-n_off > nb) then
         ab(:,1:nblocks*nb) = ab(:,nb+1:(nblocks+1)*nb)
         ab(:,nblocks*nb+1:(nblocks+1)*nb) = 0
         n_off = n_off + nb
      endif

      do iblk=1,nblocks

         ns = na_s + (iblk-1)*nb - n_off ! first column in block
         ne = ns+nb-1                    ! last column in block

         if(ns+n_off>na) exit

         ! Store Householder vector for back transformation

         hh_cnt(iblk) = hh_cnt(iblk) + 1

         hh_gath(1   ,hh_cnt(iblk),iblk) = tau
         hh_gath(2:nb,hh_cnt(iblk),iblk) = hv(2:nb)

         if(hh_cnt(iblk) == snd_limits(hh_dst(iblk)+1,iblk)-snd_limits(hh_dst(iblk),iblk)) then
            ! Wait for last transfer to finish
            call mpi_wait(ireq_hhs(iblk), MPI_STATUS_IGNORE, mpierr)
            ! Copy vectors into send buffer
            hh_send(:,1:hh_cnt(iblk),iblk) = hh_gath(:,1:hh_cnt(iblk),iblk)
            ! Send to destination
            call mpi_isend(hh_send(1,1,iblk), nb*hh_cnt(iblk), mpi_real8, &
                           global_id(hh_dst(iblk),mod(iblk+block_limits(my_pe)-1,np_cols)), &
                           10+iblk, mpi_comm, ireq_hhs(iblk), mpierr)
            ! Reset counter and increase destination row
            hh_cnt(iblk) = 0
            hh_dst(iblk) = hh_dst(iblk)+1
         endif

         ! The following code is structured in a way to keep waiting times for
         ! other PEs at a minimum, especially if there is only one block.
         ! For this reason, it requests the last column as late as possible
         ! and sends the Householder vector and the first column as early
         ! as possible.

         nc = MIN(na-ns-n_off+1,nb) ! number of columns in diagonal block
         nr = MIN(na-nb-ns-n_off+1,nb) ! rows in subdiagonal block (may be < 0!!!)
                                       ! Note that nr>=0 implies that diagonal block is full (nc==nb)!

         ! Multiply diagonal block and subdiagonal block with Householder vector

         if(iblk==nblocks .and. nc==nb) then

            ! We need the last column from the next PE.
            ! First do the matrix multiplications without last column ...

            ! Diagonal block, the contribution of the last element is added below!
            ab(1,ne) = 0
            call DSYMV('L',nc,tau,ab(1,ns),2*nb-1,hv,1,0.d0,hd,1)

            ! Subdiagonal block
            if(nr>0) call DGEMV('N',nr,nb-1,tau,ab(nb+1,ns),2*nb-1,hv,1,0.d0,hs,1)

            ! ... then request last column ...
            call mpi_recv(ab(1,ne),nb+1,mpi_real8,my_pe+1,1,mpi_comm,MPI_STATUS_IGNORE,mpierr)

            ! ... and complete the result
            hs(1:nr) = hs(1:nr) + ab(2:nr+1,ne)*tau*hv(nb)
            hd(nb) = hd(nb) + ab(1,ne)*hv(nb)*tau

         else

            ! Normal matrix multiply
            call DSYMV('L',nc,tau,ab(1,ns),2*nb-1,hv,1,0.d0,hd,1)
            if(nr>0) call DGEMV('N',nr,nb,tau,ab(nb+1,ns),2*nb-1,hv,1,0.d0,hs,1)

         endif

         ! Calculate first column of subdiagonal block and calculate new
         ! Householder transformation for this column

         hv_new(:) = 0 ! Needed, last rows must be 0 for nr < nb
         tau_new = 0

         if(nr>0) then

            ! complete (old) Householder transformation for first column

            ab(nb+1:nb+nr,ns) = ab(nb+1:nb+nr,ns) - hs(1:nr) ! Note: hv(1) == 1

            ! calculate new Householder transformation ...
            if(nr>1) then
               vnorm2 = sum(ab(nb+2:nb+nr,ns)**2)
               call hh_transform_real(ab(nb+1,ns),vnorm2,hf,tau_new)
               hv_new(1) = 1.
               hv_new(2:nr) = ab(nb+2:nb+nr,ns)*hf
               ab(nb+2:,ns) = 0
            endif

            ! ... and send it away immediatly if this is the last block

            if(iblk==nblocks) then
               call mpi_wait(ireq_hv,MPI_STATUS_IGNORE,mpierr)
               hv_s(1) = tau_new
               hv_s(2:) = hv_new(2:)
               call mpi_isend(hv_s,nb,mpi_real8,my_pe+1,2,mpi_comm,ireq_hv,mpierr)
            endif

         endif


         ! Transform diagonal block
         x = dot_product(hv(1:nc),hd(1:nc))*tau
         hd(1:nc) = hd(1:nc) - 0.5*x*hv(1:nc)

         if(my_pe>0 .and. iblk==1) then

            ! The first column of the diagonal block has to be send to the previous PE
            ! Calculate first column only ...

            ab(1:nc,ns) = ab(1:nc,ns) - hd(1:nc)*hv(1) - hv(1:nc)*hd(1)

            ! ... send it away ...

            call mpi_wait(ireq_ab,MPI_STATUS_IGNORE,mpierr)
            ab_s(1:nb+1) = ab(1:nb+1,ns)
            call mpi_isend(ab_s,nb+1,mpi_real8,my_pe-1,1,mpi_comm,ireq_ab,mpierr)

            ! ... and calculate remaining columns with rank-2 update
            if(nc>1) call DSYR2('L',nc-1,-1.d0,hd(2),1,hv(2),1,ab(1,ns+1),2*nb-1)
         else
            ! No need to  send, just a rank-2 update
            call DSYR2('L',nc,-1.d0,hd,1,hv,1,ab(1,ns),2*nb-1)
         endif

         ! Do the remaining double Householder transformation on the subdiagonal block cols 2 ... nb

         if(nr>0) then
            if(nr>1) then
               call DGEMV('T',nr,nb-1,tau_new,ab(nb,ns+1),2*nb-1,hv_new,1,0.d0,h(2),1)
               x = dot_product(hs(1:nr),hv_new(1:nr))*tau_new
               h(2:nb) = h(2:nb) - x*hv(2:nb)
               ! Unfortunately there is no BLAS routine like DSYR2 for a nonsymmetric rank 2 update ("DGER2")
               do i=2,nb
                  ab(2+nb-i:1+nb+nr-i,i+ns-1) = ab(2+nb-i:1+nb+nr-i,i+ns-1) - hv_new(1:nr)*h(i) - hs(1:nr)*hv(i)
               enddo
            else
               ! No double Householder transformation for nr=1, just complete the row
               do i=2,nb
                  ab(2+nb-i,i+ns-1) = ab(2+nb-i,i+ns-1) - hs(1)*hv(i)
               enddo
            endif
         endif

         ! Use new HH vector for the next block
         hv(:) = hv_new(:)
         tau = tau_new

      enddo

   enddo

   ! Finish the last outstanding requests
   call mpi_wait(ireq_ab,MPI_STATUS_IGNORE,mpierr)
   call mpi_wait(ireq_hv,MPI_STATUS_IGNORE,mpierr)

   call mpi_waitall(nblocks, ireq_hhs, MPI_STATUSES_IGNORE, mpierr)
   call mpi_waitall(num_chunks, ireq_hhr, MPI_STATUSES_IGNORE, mpierr)

   call mpi_barrier(mpi_comm,mpierr)

   deallocate(ab)
   deallocate(ireq_hhr, ireq_hhs)
   deallocate(hh_cnt, hh_dst)
   deallocate(hh_gath, hh_send)
   deallocate(limits, snd_limits)
   deallocate(block_limits)
   deallocate(global_id)

 end subroutine tridiag_band_real

! --------------------------------------------------------------------------------------------------

subroutine trans_ev_tridi_to_band_real(na, nev, nblk, nbw, q, ldq, mpi_comm_rows, mpi_comm_cols)

!-------------------------------------------------------------------------------
!  trans_ev_tridi_to_band_real:
!  Transforms the eigenvectors of a tridiagonal matrix back to the eigenvectors of the band matrix
!
!  Parameters
!
!  na          Order of matrix a, number of rows of matrix q
!
!  nev         Number eigenvectors to compute (= columns of matrix q)
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  nb          semi bandwith
!
!  q           On input: Eigenvectors of tridiagonal matrix
!              On output: Transformed eigenvectors
!              Distribution is like in Scalapack.
!
!  ldq         Leading dimension of q
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns/both
!
!-------------------------------------------------------------------------------

    implicit none

    integer, intent(in) :: na, nev, nblk, nbw, ldq, mpi_comm_rows, mpi_comm_cols
    real*8 q(ldq,*)

    integer np_rows, my_prow, np_cols, my_pcol

    integer i, j, ip, sweep, nbuf, l_nev, a_dim2
    integer current_n, current_local_n, current_n_start, current_n_end
    integer next_n, next_local_n, next_n_start, next_n_end
    integer bottom_msg_length, top_msg_length, next_top_msg_length
    integer thread_width, stripe_width, stripe_count, csw
    integer num_result_blocks, num_result_buffers, num_bufs_recvd
    integer a_off, current_tv_off, max_blk_size, b_off, b_len
    integer mpierr, src, src_offset, dst, offset, nfact, num_blk
    logical flag

    real*8, allocatable :: a(:,:,:,:), row(:)
    real*8, allocatable :: top_border_send_buffer(:,:), top_border_recv_buffer(:,:)
    real*8, allocatable :: bottom_border_send_buffer(:,:), bottom_border_recv_buffer(:,:)
    real*8, allocatable :: result_buffer(:,:,:)
    real*8, allocatable :: bcast_buffer(:,:)

    integer n_off
    integer, allocatable :: result_send_request(:), result_recv_request(:), limits(:)
    integer, allocatable :: top_send_request(:), bottom_send_request(:)
    integer, allocatable :: top_recv_request(:), bottom_recv_request(:)

    ! MPI send/recv tags, arbitrary

    integer, parameter :: bottom_recv_tag = 111
    integer, parameter :: top_recv_tag    = 222
    integer, parameter :: result_recv_tag = 333

    integer :: max_threads, my_thread
!$  integer :: omp_get_max_threads

    ! Just for measuring the kernel performance
    real*8 kernel_time
    integer*8 kernel_flops


    kernel_time = 1.d-100
    kernel_flops = 0

    max_threads = 1
!$  max_threads = omp_get_max_threads()

    call MPI_Comm_rank(mpi_comm_rows, my_prow, mpierr)
    call MPI_Comm_size(mpi_comm_rows, np_rows, mpierr)
    call MPI_Comm_rank(mpi_comm_cols, my_pcol, mpierr)
    call MPI_Comm_size(mpi_comm_cols, np_cols, mpierr)

    if(mod(nbw,nblk)/=0) then
      if(my_prow==0 .and. my_pcol==0) then
         print *,'ERROR: nbw=',nbw,', nblk=',nblk
         print *,'band backtransform works only for nbw==n*nblk'
         call mpi_abort(mpi_comm_world,0,mpierr)
      endif
    endif

    nfact = nbw / nblk


    ! local number of eigenvectors
    l_nev = local_index(nev, my_pcol, np_cols, nblk, -1)

    if(l_nev==0) then
        thread_width = 0
        stripe_width = 0
        stripe_count = 0
    else
        ! Suggested stripe width is 48 since 48*64 real*8 numbers should fit into
        ! every primary cache
        thread_width = (l_nev-1)/max_threads + 1 ! number of eigenvectors per OMP thread
        stripe_width = 48 ! Must be a multiple of 4
        stripe_count = (thread_width-1)/stripe_width + 1
        ! Adapt stripe width so that last one doesn't get too small
        stripe_width = (thread_width-1)/stripe_count + 1
        stripe_width = ((stripe_width+3)/4)*4 ! Must be a multiple of 4 !!!
    endif

    ! Determine the matrix distribution at the beginning

    allocate(limits(0:np_rows))

    call determine_workload(na, nbw, np_rows, limits)

    max_blk_size = maxval(limits(1:np_rows) - limits(0:np_rows-1))

    a_dim2 = max_blk_size + nbw

    allocate(a(stripe_width,a_dim2,stripe_count,max_threads))
    ! a(:,:,:,:) should be set to 0 in a parallel region, not here!

    allocate(row(l_nev))
    row(:) = 0

    ! Copy q from a block cyclic distribution into a distribution with contiguous rows,
    ! and transpose the matrix using stripes of given stripe_width for cache blocking.

    ! The peculiar way it is done below is due to the fact that the last row should be
    ! ready first since it is the first one to start below

    ! Please note about the OMP usage below:
    ! This is not for speed, but because we want the matrix a in the memory and
    ! in the cache of the correct thread (if possible)

!$omp parallel do private(my_thread), schedule(static, 1)
    do my_thread = 1, max_threads
        a(:,:,:,my_thread) = 0 ! if possible, do first touch allocation!
    enddo

    do ip = np_rows-1, 0, -1
        if(my_prow == ip) then
            ! Receive my rows which have not yet been received
            src_offset = local_index(limits(ip), my_prow, np_rows, nblk, -1)
            do i=limits(ip)+1,limits(ip+1)
                src = mod((i-1)/nblk, np_rows)
                if(src < my_prow) then
                    call MPI_Recv(row, l_nev, MPI_REAL8, src, 0, mpi_comm_rows, MPI_STATUS_IGNORE, mpierr)
!$omp parallel do private(my_thread), schedule(static, 1)
                    do my_thread = 1, max_threads
                        call unpack_row(row,i-limits(ip),my_thread)
                    enddo
                elseif(src==my_prow) then
                    src_offset = src_offset+1
                    row(:) = q(src_offset, 1:l_nev)
!$omp parallel do private(my_thread), schedule(static, 1)
                    do my_thread = 1, max_threads
                        call unpack_row(row,i-limits(ip),my_thread)
                    enddo
                endif
            enddo
            ! Send all rows which have not yet been send
            src_offset = 0
            do dst = 0, ip-1
              do i=limits(dst)+1,limits(dst+1)
                if(mod((i-1)/nblk, np_rows) == my_prow) then
                    src_offset = src_offset+1
                    row(:) = q(src_offset, 1:l_nev)
                    call MPI_Send(row, l_nev, MPI_REAL8, dst, 0, mpi_comm_rows, mpierr)
                endif
              enddo
            enddo
        else if(my_prow < ip) then
            ! Send all rows going to PE ip
            src_offset = local_index(limits(ip), my_prow, np_rows, nblk, -1)
            do i=limits(ip)+1,limits(ip+1)
                src = mod((i-1)/nblk, np_rows)
                if(src == my_prow) then
                    src_offset = src_offset+1
                    row(:) = q(src_offset, 1:l_nev)
                    call MPI_Send(row, l_nev, MPI_REAL8, ip, 0, mpi_comm_rows, mpierr)
                endif
            enddo
            ! Receive all rows from PE ip
            do i=limits(my_prow)+1,limits(my_prow+1)
                src = mod((i-1)/nblk, np_rows)
                if(src == ip) then
                    call MPI_Recv(row, l_nev, MPI_REAL8, src, 0, mpi_comm_rows, MPI_STATUS_IGNORE, mpierr)
!$omp parallel do private(my_thread), schedule(static, 1)
                    do my_thread = 1, max_threads
                        call unpack_row(row,i-limits(my_prow),my_thread)
                    enddo
                endif
            enddo
        endif
    enddo


    ! Set up result buffer queue

    num_result_blocks = ((na-1)/nblk + np_rows - my_prow) / np_rows

    num_result_buffers = 4*nfact
    allocate(result_buffer(l_nev,nblk,num_result_buffers))

    allocate(result_send_request(num_result_buffers))
    allocate(result_recv_request(num_result_buffers))
    result_send_request(:) = MPI_REQUEST_NULL
    result_recv_request(:) = MPI_REQUEST_NULL

    ! Queue up buffers

    if(my_prow > 0 .and. l_nev>0) then ! note: row 0 always sends
        do j = 1, min(num_result_buffers, num_result_blocks)
            call MPI_Irecv(result_buffer(1,1,j), l_nev*nblk, MPI_REAL8, 0, result_recv_tag, &
                           mpi_comm_rows, result_recv_request(j), mpierr)
        enddo
    endif

    num_bufs_recvd = 0 ! No buffers received yet

    ! Initialize top/bottom requests

    allocate(top_send_request(stripe_count))
    allocate(top_recv_request(stripe_count))
    allocate(bottom_send_request(stripe_count))
    allocate(bottom_recv_request(stripe_count))

    top_send_request(:) = MPI_REQUEST_NULL
    top_recv_request(:) = MPI_REQUEST_NULL
    bottom_send_request(:) = MPI_REQUEST_NULL
    bottom_recv_request(:) = MPI_REQUEST_NULL

    allocate(top_border_send_buffer(stripe_width*nbw*max_threads, stripe_count))
    allocate(top_border_recv_buffer(stripe_width*nbw*max_threads, stripe_count))
    allocate(bottom_border_send_buffer(stripe_width*nbw*max_threads, stripe_count))
    allocate(bottom_border_recv_buffer(stripe_width*nbw*max_threads, stripe_count))

    top_border_send_buffer(:,:) = 0
    top_border_recv_buffer(:,:) = 0
    bottom_border_send_buffer(:,:) = 0
    bottom_border_recv_buffer(:,:) = 0

    ! Initialize broadcast buffer

    allocate(bcast_buffer(nbw, max_blk_size))
    bcast_buffer = 0

    current_tv_off = 0 ! Offset of next row to be broadcast


    ! ------------------- start of work loop -------------------

    a_off = 0 ! offset in A (to avoid unnecessary shifts)

    top_msg_length = 0
    bottom_msg_length = 0

    do sweep = 0, (na-1)/nbw

        current_n = na - sweep*nbw
        call determine_workload(current_n, nbw, np_rows, limits)
        current_n_start = limits(my_prow)
        current_n_end   = limits(my_prow+1)
        current_local_n = current_n_end - current_n_start

        next_n = max(current_n - nbw, 0)
        call determine_workload(next_n, nbw, np_rows, limits)
        next_n_start = limits(my_prow)
        next_n_end   = limits(my_prow+1)
        next_local_n = next_n_end - next_n_start

        if(next_n_end < next_n) then
            bottom_msg_length = current_n_end - next_n_end
        else
            bottom_msg_length = 0
        endif

        if(next_local_n > 0) then
            next_top_msg_length = current_n_start - next_n_start
        else
            next_top_msg_length = 0
        endif

        if(sweep==0 .and. current_n_end < current_n .and. l_nev > 0) then
            do i = 1, stripe_count
                csw = min(stripe_width, thread_width-(i-1)*stripe_width) ! "current_stripe_width"
                b_len = csw*nbw*max_threads
                call MPI_Irecv(bottom_border_recv_buffer(1,i), b_len, MPI_REAL8, my_prow+1, bottom_recv_tag, &
                           mpi_comm_rows, bottom_recv_request(i), mpierr)
            enddo
        endif

        if(current_local_n > 1) then
            if(my_pcol == mod(sweep,np_cols)) then
                bcast_buffer(:,1:current_local_n) = hh_trans_real(:,current_tv_off+1:current_tv_off+current_local_n)
                current_tv_off = current_tv_off + current_local_n
            endif
            call mpi_bcast(bcast_buffer, nbw*current_local_n, MPI_REAL8, mod(sweep,np_cols), mpi_comm_cols, mpierr)
        else
            ! for current_local_n == 1 the one and only HH vector is 0 and not stored in hh_trans_real
            bcast_buffer(:,1) = 0
        endif

        if(l_nev == 0) cycle

        if(current_local_n > 0) then

          do i = 1, stripe_count

            ! Get real stripe width for strip i;
            ! The last OpenMP tasks may have an even smaller stripe with,
            ! but we don't care about this, i.e. we send/recv a bit too much in this case.
            ! csw: current_stripe_width

            csw = min(stripe_width, thread_width-(i-1)*stripe_width)

            !wait_b
            if(current_n_end < current_n) then
                call MPI_Wait(bottom_recv_request(i), MPI_STATUS_IGNORE, mpierr)
!$omp parallel do private(my_thread, n_off, b_len, b_off), schedule(static, 1)
                do my_thread = 1, max_threads
                    n_off = current_local_n+a_off
                    b_len = csw*nbw
                    b_off = (my_thread-1)*b_len
                    a(1:csw,n_off+1:n_off+nbw,i,my_thread) = &
                      reshape(bottom_border_recv_buffer(b_off+1:b_off+b_len,i), (/ csw, nbw /))
                enddo
                if(next_n_end < next_n) then
                    call MPI_Irecv(bottom_border_recv_buffer(1,i), csw*nbw*max_threads, &
                                   MPI_REAL8, my_prow+1, bottom_recv_tag, &
                                   mpi_comm_rows, bottom_recv_request(i), mpierr)
                endif
            endif

            if(current_local_n <= bottom_msg_length + top_msg_length) then

                !wait_t
                if(top_msg_length>0) then
                    call MPI_Wait(top_recv_request(i), MPI_STATUS_IGNORE, mpierr)
                endif

                !compute
!$omp parallel do private(my_thread, n_off, b_len, b_off), schedule(static, 1)
                do my_thread = 1, max_threads
                    if(top_msg_length>0) then
                        b_len = csw*top_msg_length
                        b_off = (my_thread-1)*b_len
                        a(1:csw,a_off+1:a_off+top_msg_length,i,my_thread) = &
                          reshape(top_border_recv_buffer(b_off+1:b_off+b_len,i), (/ csw, top_msg_length /))
                    endif
                    call compute_hh_trafo(0, current_local_n, i, my_thread)
                enddo

                !send_b
                call MPI_Wait(bottom_send_request(i), MPI_STATUS_IGNORE, mpierr)
                if(bottom_msg_length>0) then
                    n_off = current_local_n+nbw-bottom_msg_length+a_off
                    b_len = csw*bottom_msg_length*max_threads
                    bottom_border_send_buffer(1:b_len,i) = &
                        reshape(a(1:csw,n_off+1:n_off+bottom_msg_length,i,:), (/ b_len /))
                    call MPI_Isend(bottom_border_send_buffer(1,i), b_len, MPI_REAL8, my_prow+1, &
                                   top_recv_tag, mpi_comm_rows, bottom_send_request(i), mpierr)
                endif

            else

                !compute
!$omp parallel do private(my_thread, b_len, b_off), schedule(static, 1)
                do my_thread = 1, max_threads
                    call compute_hh_trafo(current_local_n - bottom_msg_length, bottom_msg_length, i, my_thread)
                enddo

                !send_b
                call MPI_Wait(bottom_send_request(i), MPI_STATUS_IGNORE, mpierr)
                if(bottom_msg_length > 0) then
                    n_off = current_local_n+nbw-bottom_msg_length+a_off
                    b_len = csw*bottom_msg_length*max_threads
                    bottom_border_send_buffer(1:b_len,i) = &
                      reshape(a(1:csw,n_off+1:n_off+bottom_msg_length,i,:), (/ b_len /))
                    call MPI_Isend(bottom_border_send_buffer(1,i), b_len, MPI_REAL8, my_prow+1, &
                                   top_recv_tag, mpi_comm_rows, bottom_send_request(i), mpierr)
                endif

                !compute
!$omp parallel do private(my_thread), schedule(static, 1)
                do my_thread = 1, max_threads
                    call compute_hh_trafo(top_msg_length, current_local_n-top_msg_length-bottom_msg_length, i, my_thread)
                enddo

                !wait_t
                if(top_msg_length>0) then
                    call MPI_Wait(top_recv_request(i), MPI_STATUS_IGNORE, mpierr)
                endif

                !compute
!$omp parallel do private(my_thread, b_len, b_off), schedule(static, 1)
                do my_thread = 1, max_threads
                    if(top_msg_length>0) then
                        b_len = csw*top_msg_length
                        b_off = (my_thread-1)*b_len
                        a(1:csw,a_off+1:a_off+top_msg_length,i,my_thread) = &
                          reshape(top_border_recv_buffer(b_off+1:b_off+b_len,i), (/ csw, top_msg_length /))
                    endif
                    call compute_hh_trafo(0, top_msg_length, i, my_thread)
                enddo
            endif

            if(next_top_msg_length > 0) then
                !request top_border data
                b_len = csw*next_top_msg_length*max_threads
                call MPI_Irecv(top_border_recv_buffer(1,i), b_len, MPI_REAL8, my_prow-1, &
                               top_recv_tag, mpi_comm_rows, top_recv_request(i), mpierr)
            endif

            !send_t
            if(my_prow > 0) then
                call MPI_Wait(top_send_request(i), MPI_STATUS_IGNORE, mpierr)
                b_len = csw*nbw*max_threads
                top_border_send_buffer(1:b_len,i) = reshape(a(1:csw,a_off+1:a_off+nbw,i,:), (/ b_len /))
                call MPI_Isend(top_border_send_buffer(1,i), b_len, MPI_REAL8, &
                               my_prow-1, bottom_recv_tag, &
                               mpi_comm_rows, top_send_request(i), mpierr)
            endif

            ! Care that there are not too many outstanding top_recv_request's
            if(stripe_count > 1) then
                if(i>1) then
                    call MPI_Wait(top_recv_request(i-1), MPI_STATUS_IGNORE, mpierr)
                else
                    call MPI_Wait(top_recv_request(stripe_count), MPI_STATUS_IGNORE, mpierr)
                endif
            endif

          enddo

          top_msg_length = next_top_msg_length

        else
            ! wait for last top_send_request
          do i = 1, stripe_count
            call MPI_Wait(top_send_request(i), MPI_STATUS_IGNORE, mpierr)
          enddo
        endif

        ! Care about the result

        if(my_prow == 0) then

            ! topmost process sends nbw rows to destination processes

            do j=0,nfact-1

                num_blk = sweep*nfact+j ! global number of destination block, 0 based
                if(num_blk*nblk >= na) exit

                nbuf = mod(num_blk, num_result_buffers) + 1 ! buffer number to get this block

                call MPI_Wait(result_send_request(nbuf), MPI_STATUS_IGNORE, mpierr)

                dst = mod(num_blk, np_rows)

                if(dst == 0) then
                    do i = 1, min(na - num_blk*nblk, nblk)
                        call pack_row(row, j*nblk+i+a_off)
                        q((num_blk/np_rows)*nblk+i,1:l_nev) = row(:)
                    enddo
                else
                    do i = 1, nblk
                        call pack_row(result_buffer(:,i,nbuf),j*nblk+i+a_off)
                    enddo
                    call MPI_Isend(result_buffer(1,1,nbuf), l_nev*nblk, MPI_REAL8, dst, &
                                   result_recv_tag, mpi_comm_rows, result_send_request(nbuf), mpierr)
                endif
            enddo

        else

           ! receive and store final result

            do j = num_bufs_recvd, num_result_blocks-1

                nbuf = mod(j, num_result_buffers) + 1 ! buffer number to get this block

                ! If there is still work to do, just test for the next result request
                ! and leave the loop if it is not ready, otherwise wait for all
                ! outstanding requests

                if(next_local_n > 0) then
                    call MPI_Test(result_recv_request(nbuf), flag, MPI_STATUS_IGNORE, mpierr)
                    if(.not.flag) exit
                else
                    call MPI_Wait(result_recv_request(nbuf), MPI_STATUS_IGNORE, mpierr)
                endif

                ! Fill result buffer into q
                num_blk = j*np_rows + my_prow ! global number of current block, 0 based
                do i = 1, min(na - num_blk*nblk, nblk)
                    q(j*nblk+i, 1:l_nev) = result_buffer(1:l_nev, i, nbuf)
                enddo

                ! Queue result buffer again if there are outstanding blocks left
                if(j+num_result_buffers < num_result_blocks) &
                    call MPI_Irecv(result_buffer(1,1,nbuf), l_nev*nblk, MPI_REAL8, 0, result_recv_tag, &
                                   mpi_comm_rows, result_recv_request(nbuf), mpierr)

            enddo
            num_bufs_recvd = j

        endif

        ! Shift the remaining rows to the front of A (if necessary)

        offset = nbw - top_msg_length
        if(offset<0) then
            print *,'internal error, offset for shifting = ',offset
            call MPI_Abort(MPI_COMM_WORLD, 1, mpierr)
        endif
        a_off = a_off + offset
        if(a_off + next_local_n + nbw > a_dim2) then
!$omp parallel do private(my_thread, i, j), schedule(static, 1)
            do my_thread = 1, max_threads
                do i = 1, stripe_count
                    do j = top_msg_length+1, top_msg_length+next_local_n
                       A(:,j,i,my_thread) = A(:,j+a_off,i,my_thread)
                    enddo
                enddo
            enddo
            a_off = 0
        endif

    enddo

    ! Just for safety:
    if(ANY(top_send_request    /= MPI_REQUEST_NULL)) print *,'*** ERROR top_send_request ***',my_prow,my_pcol
    if(ANY(bottom_send_request /= MPI_REQUEST_NULL)) print *,'*** ERROR bottom_send_request ***',my_prow,my_pcol
    if(ANY(top_recv_request    /= MPI_REQUEST_NULL)) print *,'*** ERROR top_recv_request ***',my_prow,my_pcol
    if(ANY(bottom_recv_request /= MPI_REQUEST_NULL)) print *,'*** ERROR bottom_recv_request ***',my_prow,my_pcol

    if(my_prow == 0) then
        call MPI_Waitall(num_result_buffers, result_send_request, MPI_STATUSES_IGNORE, mpierr)
    endif

    if(ANY(result_send_request /= MPI_REQUEST_NULL)) print *,'*** ERROR result_send_request ***',my_prow,my_pcol
    if(ANY(result_recv_request /= MPI_REQUEST_NULL)) print *,'*** ERROR result_recv_request ***',my_prow,my_pcol

    if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) &
        print '(" Kernel time:",f10.3," MFlops: ",f10.3)', kernel_time, kernel_flops/kernel_time*1.d-6

    ! deallocate all working space

    deallocate(a)
    deallocate(row)
    deallocate(limits)
    deallocate(result_send_request)
    deallocate(result_recv_request)
    deallocate(top_border_send_buffer)
    deallocate(top_border_recv_buffer)
    deallocate(bottom_border_send_buffer)
    deallocate(bottom_border_recv_buffer)
    deallocate(result_buffer)
    deallocate(bcast_buffer)
    deallocate(top_send_request)
    deallocate(top_recv_request)
    deallocate(bottom_send_request)
    deallocate(bottom_recv_request)

contains

    subroutine pack_row(row, n)
        real*8 row(:)
        integer n, i, noff, nl, nt

        do nt = 1, max_threads
            do i = 1, stripe_count
                noff = (nt-1)*thread_width + (i-1)*stripe_width
                nl   = min(stripe_width, nt*thread_width-noff, l_nev-noff)
                if(nl<=0) exit
                row(noff+1:noff+nl) = a(1:nl,n,i,nt)
            enddo
        enddo

    end subroutine

    subroutine unpack_row(row, n, my_thread)

        ! Private variables in OMP regions (my_thread) should better be in the argument list!
        integer, intent(in) :: n, my_thread
        real*8, intent(in)  :: row(:)
        integer i, noff, nl

        do i=1,stripe_count
            noff = (my_thread-1)*thread_width + (i-1)*stripe_width
            nl   = min(stripe_width, my_thread*thread_width-noff, l_nev-noff)
            if(nl<=0) exit
            a(1:nl,n,i,my_thread) = row(noff+1:noff+nl)
        enddo

    end subroutine

    subroutine compute_hh_trafo(off, ncols, istripe, my_thread)

        ! Private variables in OMP regions (my_thread) should better be in the argument list!
        integer, intent(in) :: off, ncols, istripe, my_thread
        integer j, nl, noff
        real*8 w(nbw,2), ttt

        ttt = mpi_wtime()
        if(istripe<stripe_count) then
          nl = stripe_width
        else
          noff = (my_thread-1)*thread_width + (istripe-1)*stripe_width
          nl = min(my_thread*thread_width-noff, l_nev-noff)
          if(nl<=0) return
        endif
        do j = ncols, 2, -2
            w(:,1) = bcast_buffer(1:nbw,j+off)
            w(:,2) = bcast_buffer(1:nbw,j+off-1)
            call double_hh_trafo(a(1,j+off+a_off-1,istripe,my_thread), w, nbw, nl, stripe_width, nbw)
        enddo
        if(j==1) call single_hh_trafo(a(1,1+off+a_off,istripe,my_thread),bcast_buffer(1,off+1), nbw, nl, stripe_width)
        if(my_thread==1) then
            kernel_flops = kernel_flops + 4*int(nl,8)*int(ncols,8)*int(nbw,8)
            kernel_time  = kernel_time + mpi_wtime()-ttt
        endif

    end subroutine

end subroutine

!-------------------------------------------------------------------------------

subroutine single_hh_trafo(q, hh, nb, nq, ldq)

    ! Perform single real Householder transformation.
    ! This routine is not performance critical and thus it is coded here in Fortran

    implicit none
    integer nb, nq, ldq
    real*8 q(ldq, *), hh(*)

    integer i
    real*8 v(nq)

    ! v = q * hh
    v(:) = q(1:nq,1)
    do i=2,nb
        v(:) = v(:) + q(1:nq,i) * hh(i)
    enddo

    ! v = v * tau
    v(:) = v(:) * hh(1)

    ! q = q - v * hh**T
    q(1:nq,1) = q(1:nq,1) - v(:)
    do i=2,nb
        q(1:nq,i) = q(1:nq,i) - v(:) * hh(i)
    enddo

end subroutine

!-------------------------------------------------------------------------------

subroutine determine_workload(na, nb, nprocs, limits)

    integer, intent(in) :: na, nb, nprocs
    integer, intent(out) :: limits(0:nprocs)

    integer i

    if(na <= 0) then
        limits(:) = 0
        return
    endif

    if(nb*nprocs > na) then
        ! there is not enough work for all
        do i = 0, nprocs
            limits(i) = min(na, i*nb)
        enddo
    else
        do i = 0, nprocs
            limits(i) = (i*na)/nprocs
        enddo
    endif

end subroutine

!-------------------------------------------------------------------------------

subroutine bandred_complex(na, a, lda, nblk, nbw, mpi_comm_rows, mpi_comm_cols, tmat)

!-------------------------------------------------------------------------------
!  bandred_complex: Reduces a distributed hermitian matrix to band form
!
!  Parameters
!
!  na          Order of matrix
!
!  a(lda,*)    Distributed matrix which should be reduced.
!              Distribution is like in Scalapack.
!              Opposed to Scalapack, a(:,:) must be set completely (upper and lower half)
!              a(:,:) is overwritten on exit with the band and the Householder vectors
!              in the upper half.
!
!  lda         Leading dimension of a
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  nbw         semi bandwith of output matrix
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!
!  tmat(nbw,nbw,num_blocks)    where num_blocks = (na-1)/nbw + 1
!              Factors for the Householder vectors (returned), needed for back transformation
!
!-------------------------------------------------------------------------------

   implicit none

   integer na, lda, nblk, nbw, mpi_comm_rows, mpi_comm_cols
   complex*16 a(lda,*), tmat(nbw,nbw,*)

   complex*16, parameter :: CZERO = (0.d0,0.d0), CONE = (1.d0,0.d0)

   integer my_prow, my_pcol, np_rows, np_cols, mpierr
   integer l_cols, l_rows
   integer i, j, lcs, lce, lre, lc, lr, cur_pcol, n_cols, nrow
   integer istep, ncol, lch, lcx, nlc
   integer tile_size, l_rows_tile, l_cols_tile

   real*8 vnorm2
   complex*16 xf, aux1(nbw), aux2(nbw), vrl, tau, vav(nbw,nbw)

   complex*16, allocatable:: tmp(:,:), vr(:), vmr(:,:), umc(:,:)

   integer pcol, prow
   pcol(i) = MOD((i-1)/nblk,np_cols) !Processor col for global col number
   prow(i) = MOD((i-1)/nblk,np_rows) !Processor row for global row number


   call mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   call mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   call mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   call mpi_comm_size(mpi_comm_cols,np_cols,mpierr)

   ! Semibandwith nbw must be a multiple of blocksize nblk

   if(mod(nbw,nblk)/=0) then
      if(my_prow==0 .and. my_pcol==0) then
         print *,'ERROR: nbw=',nbw,', nblk=',nblk
         print *,'ELPA2 works only for nbw==n*nblk'
         call mpi_abort(mpi_comm_world,0,mpierr)
      endif
   endif

   ! Matrix is split into tiles; work is done only for tiles on the diagonal or above

   tile_size = nblk*least_common_multiple(np_rows,np_cols) ! minimum global tile size
   tile_size = ((128*max(np_rows,np_cols)-1)/tile_size+1)*tile_size ! make local tiles at least 128 wide

   l_rows_tile = tile_size/np_rows ! local rows of a tile
   l_cols_tile = tile_size/np_cols ! local cols of a tile

   do istep = (na-1)/nbw, 1, -1

      n_cols = MIN(na,(istep+1)*nbw) - istep*nbw ! Number of columns in current step

      ! Number of local columns/rows of remaining matrix
      l_cols = local_index(istep*nbw, my_pcol, np_cols, nblk, -1)
      l_rows = local_index(istep*nbw, my_prow, np_rows, nblk, -1)

      ! Allocate vmr and umc to their exact sizes so that they can be used in bcasts and reduces

      allocate(vmr(max(l_rows,1),2*n_cols))
      allocate(umc(max(l_cols,1),2*n_cols))

      allocate(vr(l_rows+1))

      vmr(1:l_rows,1:n_cols) = 0.
      vr(:) = 0
      tmat(:,:,istep) = 0

      ! Reduce current block to lower triangular form

      do lc = n_cols, 1, -1

         ncol = istep*nbw + lc ! absolute column number of householder vector
         nrow = ncol - nbw ! Absolute number of pivot row

         lr  = local_index(nrow, my_prow, np_rows, nblk, -1) ! current row length
         lch = local_index(ncol, my_pcol, np_cols, nblk, -1) ! HV local column number

         tau = 0

         if(nrow == 1) exit ! Nothing to do

         cur_pcol = pcol(ncol) ! Processor column owning current block

         if(my_pcol==cur_pcol) then

            ! Get vector to be transformed; distribute last element and norm of
            ! remaining elements to all procs in current column

            vr(1:lr) = a(1:lr,lch) ! vector to be transformed

            if(my_prow==prow(nrow)) then
               aux1(1) = dot_product(vr(1:lr-1),vr(1:lr-1))
               aux1(2) = vr(lr)
            else
               aux1(1) = dot_product(vr(1:lr),vr(1:lr))
               aux1(2) = 0.
            endif

            call mpi_allreduce(aux1,aux2,2,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_comm_rows,mpierr)

            vnorm2 = aux2(1)
            vrl    = aux2(2)

            ! Householder transformation

            call hh_transform_complex(vrl, vnorm2, xf, tau)

            ! Scale vr and store Householder vector for back transformation

            vr(1:lr) = vr(1:lr) * xf
            if(my_prow==prow(nrow)) then
               a(1:lr-1,lch) = vr(1:lr-1)
               a(lr,lch) = vrl
               vr(lr) = 1.
            else
               a(1:lr,lch) = vr(1:lr)
            endif

         endif

         ! Broadcast Householder vector and tau along columns

         vr(lr+1) = tau
         call MPI_Bcast(vr,lr+1,MPI_DOUBLE_COMPLEX,cur_pcol,mpi_comm_cols,mpierr)
         vmr(1:lr,lc) = vr(1:lr)
         tau = vr(lr+1)
         tmat(lc,lc,istep) = conjg(tau) ! Store tau in diagonal of tmat

         ! Transform remaining columns in current block with Householder vector

         ! Local dot product

         aux1 = 0

         nlc = 0 ! number of local columns
         do j=1,lc-1
            lcx = local_index(istep*nbw+j, my_pcol, np_cols, nblk, 0)
            if(lcx>0) then
               nlc = nlc+1
               aux1(nlc) = dot_product(vr(1:lr),a(1:lr,lcx))
            endif
         enddo

         ! Get global dot products
         if(nlc>0) call mpi_allreduce(aux1,aux2,nlc,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_comm_rows,mpierr)

         ! Transform

         nlc = 0
         do j=1,lc-1
            lcx = local_index(istep*nbw+j, my_pcol, np_cols, nblk, 0)
            if(lcx>0) then
               nlc = nlc+1
               a(1:lr,lcx) = a(1:lr,lcx) - conjg(tau)*aux2(nlc)*vr(1:lr)
            endif
         enddo

      enddo

      ! Calculate scalar products of stored Householder vectors.
      ! This can be done in different ways, we use zherk

      vav = 0
      if(l_rows>0) &
         call zherk('U','C',n_cols,l_rows,CONE,vmr,ubound(vmr,1),CZERO,vav,ubound(vav,1))
      call herm_matrix_allreduce(n_cols,vav,ubound(vav,1),mpi_comm_rows)

      ! Calculate triangular matrix T for block Householder Transformation

      do lc=n_cols,1,-1
         tau = tmat(lc,lc,istep)
         if(lc<n_cols) then
            call ztrmv('U','C','N',n_cols-lc,tmat(lc+1,lc+1,istep),ubound(tmat,1),vav(lc+1,lc),1)
            tmat(lc,lc+1:n_cols,istep) = -tau * conjg(vav(lc+1:n_cols,lc))
         endif
      enddo

      ! Transpose vmr -> vmc (stored in umc, second half)

      call elpa_transpose_vectors  (vmr, 2*ubound(vmr,1), mpi_comm_rows, &
                                    umc(1,n_cols+1), 2*ubound(umc,1), mpi_comm_cols, &
                                    1, 2*istep*nbw, n_cols, 2*nblk)

      ! Calculate umc = A**T * vmr
      ! Note that the distributed A has to be transposed
      ! Opposed to direct tridiagonalization there is no need to use the cache locality
      ! of the tiles, so we can use strips of the matrix

      umc(1:l_cols,1:n_cols) = 0.d0
      vmr(1:l_rows,n_cols+1:2*n_cols) = 0
      if(l_cols>0 .and. l_rows>0) then
         do i=0,(istep*nbw-1)/tile_size

            lcs = i*l_cols_tile+1
            lce = min(l_cols,(i+1)*l_cols_tile)
            if(lce<lcs) cycle

            lre = min(l_rows,(i+1)*l_rows_tile)
            call ZGEMM('C','N',lce-lcs+1,n_cols,lre,CONE,a(1,lcs),ubound(a,1), &
                       vmr,ubound(vmr,1),CONE,umc(lcs,1),ubound(umc,1))

            if(i==0) cycle
            lre = min(l_rows,i*l_rows_tile)
            call ZGEMM('N','N',lre,n_cols,lce-lcs+1,CONE,a(1,lcs),lda, &
                       umc(lcs,n_cols+1),ubound(umc,1),CONE,vmr(1,n_cols+1),ubound(vmr,1))
         enddo
      endif

      ! Sum up all ur(:) parts along rows and add them to the uc(:) parts
      ! on the processors containing the diagonal
      ! This is only necessary if ur has been calculated, i.e. if the
      ! global tile size is smaller than the global remaining matrix

      if(tile_size < istep*nbw) then
         call elpa_reduce_add_vectors  (vmr(1,n_cols+1),2*ubound(vmr,1),mpi_comm_rows, &
                                        umc, 2*ubound(umc,1), mpi_comm_cols, &
                                        2*istep*nbw, n_cols, 2*nblk)
      endif

      if(l_cols>0) then
         allocate(tmp(l_cols,n_cols))
         call mpi_allreduce(umc,tmp,l_cols*n_cols,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_comm_rows,mpierr)
         umc(1:l_cols,1:n_cols) = tmp(1:l_cols,1:n_cols)
         deallocate(tmp)
      endif

      ! U = U * Tmat**T

      call ztrmm('Right','Upper','C','Nonunit',l_cols,n_cols,CONE,tmat(1,1,istep),ubound(tmat,1),umc,ubound(umc,1))

      ! VAV = Tmat * V**T * A * V * Tmat**T = (U*Tmat**T)**T * V * Tmat**T

      call zgemm('C','N',n_cols,n_cols,l_cols,CONE,umc,ubound(umc,1),umc(1,n_cols+1),ubound(umc,1),CZERO,vav,ubound(vav,1))
      call ztrmm('Right','Upper','C','Nonunit',n_cols,n_cols,CONE,tmat(1,1,istep),ubound(tmat,1),vav,ubound(vav,1))

      call herm_matrix_allreduce(n_cols,vav,ubound(vav,1),mpi_comm_cols)

      ! U = U - 0.5 * V * VAV
      call zgemm('N','N',l_cols,n_cols,n_cols,(-0.5d0,0.d0),umc(1,n_cols+1),ubound(umc,1),vav,ubound(vav,1),CONE,umc,ubound(umc,1))

      ! Transpose umc -> umr (stored in vmr, second half)

       call elpa_transpose_vectors  (umc, 2*ubound(umc,1), mpi_comm_cols, &
                                     vmr(1,n_cols+1), 2*ubound(vmr,1), mpi_comm_rows, &
                                     1, 2*istep*nbw, n_cols, 2*nblk)

      ! A = A - V*U**T - U*V**T

      do i=0,(istep*nbw-1)/tile_size
         lcs = i*l_cols_tile+1
         lce = min(l_cols,(i+1)*l_cols_tile)
         lre = min(l_rows,(i+1)*l_rows_tile)
         if(lce<lcs .or. lre<1) cycle
         call zgemm('N','C',lre,lce-lcs+1,2*n_cols,-CONE, &
                    vmr,ubound(vmr,1),umc(lcs,1),ubound(umc,1), &
                    CONE,a(1,lcs),lda)
      enddo

      deallocate(vmr, umc, vr)

   enddo

end subroutine bandred_complex

!-------------------------------------------------------------------------------

subroutine herm_matrix_allreduce(n,a,lda,comm)

!-------------------------------------------------------------------------------
!  herm_matrix_allreduce: Does an mpi_allreduce for a hermitian matrix A.
!  On entry, only the upper half of A needs to be set
!  On exit, the complete matrix is set
!-------------------------------------------------------------------------------

   implicit none
   integer n, lda, comm
   complex*16 a(lda,*)

   integer i, nc, mpierr
   complex*16 h1(n*n), h2(n*n)

   nc = 0
   do i=1,n
      h1(nc+1:nc+i) = a(1:i,i)
      nc = nc+i
   enddo

   call mpi_allreduce(h1,h2,nc,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,mpierr)

   nc = 0
   do i=1,n
      a(1:i,i) = h2(nc+1:nc+i)
      a(i,1:i-1) = conjg(a(1:i-1,i))
      nc = nc+i
   enddo

end subroutine herm_matrix_allreduce

!-------------------------------------------------------------------------------

subroutine trans_ev_band_to_full_complex(na, nqc, nblk, nbw, a, lda, tmat, q, ldq, mpi_comm_rows, mpi_comm_cols)

!-------------------------------------------------------------------------------
!  trans_ev_band_to_full_complex:
!  Transforms the eigenvectors of a band matrix back to the eigenvectors of the original matrix
!
!  Parameters
!
!  na          Order of matrix a, number of rows of matrix q
!
!  nqc         Number of columns of matrix q
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  nbw         semi bandwith
!
!  a(lda,*)    Matrix containing the Householder vectors (i.e. matrix a after bandred_complex)
!              Distribution is like in Scalapack.
!
!  lda         Leading dimension of a
!
!  tmat(nbw,nbw,.) Factors returned by bandred_complex
!
!  q           On input: Eigenvectors of band matrix
!              On output: Transformed eigenvectors
!              Distribution is like in Scalapack.
!
!  ldq         Leading dimension of q
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!
!-------------------------------------------------------------------------------

   implicit none

   integer na, nqc, lda, ldq, nblk, nbw, mpi_comm_rows, mpi_comm_cols
   complex*16 a(lda,*), q(ldq,*), tmat(nbw, nbw, *)

   complex*16, parameter :: CZERO = (0.d0,0.d0), CONE = (1.d0,0.d0)

   integer my_prow, my_pcol, np_rows, np_cols, mpierr
   integer max_blocks_row, max_blocks_col, max_local_rows, max_local_cols
   integer l_cols, l_rows, l_colh, n_cols
   integer istep, lc, ncol, nrow, nb, ns

   complex*16, allocatable:: tmp1(:), tmp2(:), hvb(:), hvm(:,:)

   integer pcol, prow, i
   pcol(i) = MOD((i-1)/nblk,np_cols) !Processor col for global col number
   prow(i) = MOD((i-1)/nblk,np_rows) !Processor row for global row number


   call mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   call mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   call mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   call mpi_comm_size(mpi_comm_cols,np_cols,mpierr)


   max_blocks_row = ((na -1)/nblk)/np_rows + 1  ! Rows of A
   max_blocks_col = ((nqc-1)/nblk)/np_cols + 1  ! Columns of q!

   max_local_rows = max_blocks_row*nblk
   max_local_cols = max_blocks_col*nblk

   allocate(tmp1(max_local_cols*nbw))
   allocate(tmp2(max_local_cols*nbw))
   allocate(hvb(max_local_rows*nbw))
   allocate(hvm(max_local_rows,nbw))

   hvm = 0   ! Must be set to 0 !!!
   hvb = 0   ! Safety only

   l_cols = local_index(nqc, my_pcol, np_cols, nblk, -1) ! Local columns of q

   do istep=1,(na-1)/nbw

      n_cols = MIN(na,(istep+1)*nbw) - istep*nbw ! Number of columns in current step

      ! Broadcast all Householder vectors for current step compressed in hvb

      nb = 0
      ns = 0

      do lc = 1, n_cols
         ncol = istep*nbw + lc ! absolute column number of householder vector
         nrow = ncol - nbw ! absolute number of pivot row

         l_rows = local_index(nrow-1, my_prow, np_rows, nblk, -1) ! row length for bcast
         l_colh = local_index(ncol  , my_pcol, np_cols, nblk, -1) ! HV local column number

         if(my_pcol==pcol(ncol)) hvb(nb+1:nb+l_rows) = a(1:l_rows,l_colh)

         nb = nb+l_rows

         if(lc==n_cols .or. mod(ncol,nblk)==0) then
            call MPI_Bcast(hvb(ns+1),nb-ns,MPI_DOUBLE_COMPLEX,pcol(ncol),mpi_comm_cols,mpierr)
            ns = nb
         endif
      enddo

      ! Expand compressed Householder vectors into matrix hvm

      nb = 0
      do lc = 1, n_cols
         nrow = (istep-1)*nbw+lc ! absolute number of pivot row
         l_rows = local_index(nrow-1, my_prow, np_rows, nblk, -1) ! row length for bcast

         hvm(1:l_rows,lc) = hvb(nb+1:nb+l_rows)
         if(my_prow==prow(nrow)) hvm(l_rows+1,lc) = 1.

         nb = nb+l_rows
      enddo

      l_rows = local_index(MIN(na,(istep+1)*nbw), my_prow, np_rows, nblk, -1)

      ! Q = Q - V * T**T * V**T * Q

      if(l_rows>0) then
         call zgemm('C','N',n_cols,l_cols,l_rows,CONE,hvm,ubound(hvm,1), &
                    q,ldq,CZERO,tmp1,n_cols)
      else
         tmp1(1:l_cols*n_cols) = 0
      endif
      call mpi_allreduce(tmp1,tmp2,n_cols*l_cols,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_comm_rows,mpierr)
      if(l_rows>0) then
         call ztrmm('L','U','C','N',n_cols,l_cols,CONE,tmat(1,1,istep),ubound(tmat,1),tmp2,n_cols)
         call zgemm('N','N',l_rows,l_cols,n_cols,-CONE,hvm,ubound(hvm,1), &
                    tmp2,n_cols,CONE,q,ldq)
      endif

   enddo

   deallocate(tmp1, tmp2, hvb, hvm)


end subroutine trans_ev_band_to_full_complex

!---------------------------------------------------------------------------------------------------

subroutine tridiag_band_complex(na, nb, nblk, a, lda, d, e, mpi_comm_rows, mpi_comm_cols, mpi_comm)

!-------------------------------------------------------------------------------
! tridiag_band_complex:
! Reduces a real symmetric band matrix to tridiagonal form
!
!  na          Order of matrix a
!
!  nb          Semi bandwith
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  a(lda,*)    Distributed system matrix reduced to banded form in the upper diagonal
!
!  lda         Leading dimension of a
!
!  d(na)       Diagonal of tridiagonal matrix, set only on PE 0 (output)
!
!  e(na)       Subdiagonal of tridiagonal matrix, set only on PE 0 (output)
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!  mpi_comm
!              MPI-Communicator for the total processor set
!-------------------------------------------------------------------------------

   implicit none

   integer, intent(in) ::  na, nb, nblk, lda, mpi_comm_rows, mpi_comm_cols, mpi_comm
   complex*16, intent(in) :: a(lda,*)
   real*8, intent(out) :: d(na), e(na) ! set only on PE 0


   real*8 vnorm2
   complex*16 hv(nb), tau, x, h(nb), ab_s(1+nb), hv_s(nb), hv_new(nb), tau_new, hf
   complex*16 hd(nb), hs(nb)

   integer i, j, n, nc, nr, ns, ne, istep, iblk, nblocks_total, nblocks, nt
   integer my_pe, n_pes, mpierr
   integer my_prow, np_rows, my_pcol, np_cols
   integer ireq_ab, ireq_hv
   integer na_s, nx, num_hh_vecs, num_chunks, local_size, max_blk_size, n_off
   integer, allocatable :: ireq_hhr(:), ireq_hhs(:), global_id(:,:), hh_cnt(:), hh_dst(:)
   integer, allocatable :: limits(:), snd_limits(:,:)
   integer, allocatable :: block_limits(:)
   complex*16, allocatable :: ab(:,:), hh_gath(:,:,:), hh_send(:,:,:)
   ! dummies for calling redist_band
   real*8 :: r_a(1,1), r_ab(1,1)


   call mpi_comm_rank(mpi_comm,my_pe,mpierr)
   call mpi_comm_size(mpi_comm,n_pes,mpierr)

   call mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   call mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   call mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   call mpi_comm_size(mpi_comm_cols,np_cols,mpierr)

   ! Get global_id mapping 2D procssor coordinates to global id

   allocate(global_id(0:np_rows-1,0:np_cols-1))
   global_id(:,:) = 0
   global_id(my_prow, my_pcol) = my_pe

   call mpi_allreduce(mpi_in_place, global_id, np_rows*np_cols, mpi_integer, mpi_sum, mpi_comm, mpierr)


   ! Total number of blocks in the band:

   nblocks_total = (na-1)/nb + 1

   ! Set work distribution

   allocate(block_limits(0:n_pes))
   call divide_band(nblocks_total, n_pes, block_limits)

   ! nblocks: the number of blocks for my task
   nblocks = block_limits(my_pe+1) - block_limits(my_pe)

   ! allocate the part of the band matrix which is needed by this PE
   ! The size is 1 block larger than needed to avoid extensive shifts
   allocate(ab(2*nb,(nblocks+1)*nb))
   ab = 0 ! needed for lower half, the extra block should also be set to 0 for safety

   ! n_off: Offset of ab within band
   n_off = block_limits(my_pe)*nb

   ! Redistribute band in a to ab
   call redist_band(.false., r_a, a, lda, na, nblk, nb, mpi_comm_rows, mpi_comm_cols, mpi_comm, r_ab, ab)

   ! Calculate the workload for each sweep in the back transformation
   ! and the space requirements to hold the HH vectors

   allocate(limits(0:np_rows))
   call determine_workload(na, nb, np_rows, limits)
   max_blk_size = maxval(limits(1:np_rows) - limits(0:np_rows-1))

   num_hh_vecs = 0
   num_chunks  = 0
   nx = na
   do n = 1, nblocks_total
      call determine_workload(nx, nb, np_rows, limits)
      local_size = limits(my_prow+1) - limits(my_prow)
      ! add to number of householder vectors
      ! please note: for nx==1 the one and only HH vector is 0 and is neither calculated nor send below!
      if(mod(n-1,np_cols) == my_pcol .and. local_size>0 .and. nx>1) then
         num_hh_vecs = num_hh_vecs + local_size
         num_chunks  = num_chunks+1
      endif
      nx = nx - nb
   enddo

   ! Allocate space for HH vectors

   allocate(hh_trans_complex(nb,num_hh_vecs))

   ! Allocate and init MPI requests

   allocate(ireq_hhr(num_chunks)) ! Recv requests
   allocate(ireq_hhs(nblocks))    ! Send requests

   num_hh_vecs = 0
   num_chunks  = 0
   nx = na
   nt = 0
   do n = 1, nblocks_total
      call determine_workload(nx, nb, np_rows, limits)
      local_size = limits(my_prow+1) - limits(my_prow)
      if(mod(n-1,np_cols) == my_pcol .and. local_size>0 .and. nx>1) then
         num_chunks  = num_chunks+1
         call mpi_irecv(hh_trans_complex(1,num_hh_vecs+1), nb*local_size, MPI_COMPLEX16, nt, &
                        10+n-block_limits(nt), mpi_comm, ireq_hhr(num_chunks), mpierr)
         num_hh_vecs = num_hh_vecs + local_size
      endif
      nx = nx - nb
      if(n == block_limits(nt+1)) then
         nt = nt + 1
      endif
   enddo

   ireq_hhs(:) = MPI_REQUEST_NULL

   ! Buffers for gathering/sending the HH vectors

   allocate(hh_gath(nb,max_blk_size,nblocks)) ! gathers HH vectors
   allocate(hh_send(nb,max_blk_size,nblocks)) ! send buffer for HH vectors
   hh_gath(:,:,:) = 0
   hh_send(:,:,:) = 0

   ! Some counters

   allocate(hh_cnt(nblocks))
   allocate(hh_dst(nblocks))

   hh_cnt(:) = 1 ! The first transfomation vector is always 0 and not calculated at all
   hh_dst(:) = 0 ! PE number for receive

   ireq_ab = MPI_REQUEST_NULL
   ireq_hv = MPI_REQUEST_NULL

   ! Limits for sending

   allocate(snd_limits(0:np_rows,nblocks))

   do iblk=1,nblocks
      call determine_workload(na-(iblk+block_limits(my_pe)-1)*nb, nb, np_rows, snd_limits(:,iblk))
   enddo

   ! ---------------------------------------------------------------------------
   ! Start of calculations

   na_s = block_limits(my_pe)*nb + 1

   if(my_pe>0 .and. na_s<=na) then
      ! send first column to previous PE
      ! Only the PE owning the diagonal does that (sending 1 element of the subdiagonal block also)
      ab_s(1:nb+1) = ab(1:nb+1,na_s-n_off)
      call mpi_isend(ab_s,nb+1,MPI_COMPLEX16,my_pe-1,1,mpi_comm,ireq_ab,mpierr)
   endif

   do istep=1,na-1

      if(my_pe==0) then
         n = MIN(na-na_s,nb) ! number of rows to be reduced
         hv(:) = 0
         tau = 0
         ! Transform first column of remaining matrix
         ! Opposed to the real case, the last step (istep=na-1) is needed here for making
         ! the last subdiagonal element a real number
         vnorm2 = sum(dble(ab(3:n+1,na_s-n_off))**2+dimag(ab(3:n+1,na_s-n_off))**2)
         if(n<2) vnorm2 = 0. ! Safety only
         call hh_transform_complex(ab(2,na_s-n_off),vnorm2,hf,tau)

         hv(1) = 1
         hv(2:n) = ab(3:n+1,na_s-n_off)*hf

         d(istep) = ab(1,na_s-n_off)
         e(istep) = ab(2,na_s-n_off)
         if(istep == na-1) then
            d(na) = ab(1,na_s+1-n_off)
            e(na) = 0
         endif
      else
         if(na>na_s) then
            ! Receive Householder vector from previous task, from PE owning subdiagonal
            call mpi_recv(hv,nb,MPI_COMPLEX16,my_pe-1,2,mpi_comm,MPI_STATUS_IGNORE,mpierr)
            tau = hv(1)
            hv(1) = 1.
         endif
      endif

      na_s = na_s+1
      if(na_s-n_off > nb) then
         ab(:,1:nblocks*nb) = ab(:,nb+1:(nblocks+1)*nb)
         ab(:,nblocks*nb+1:(nblocks+1)*nb) = 0
         n_off = n_off + nb
      endif

      do iblk=1,nblocks

         ns = na_s + (iblk-1)*nb - n_off ! first column in block
         ne = ns+nb-1                    ! last column in block

         if(ns+n_off>na) exit

         ! Store Householder vector for back transformation

         hh_cnt(iblk) = hh_cnt(iblk) + 1

         hh_gath(1   ,hh_cnt(iblk),iblk) = tau
         hh_gath(2:nb,hh_cnt(iblk),iblk) = hv(2:nb)

         if(hh_cnt(iblk) == snd_limits(hh_dst(iblk)+1,iblk)-snd_limits(hh_dst(iblk),iblk)) then
            ! Wait for last transfer to finish
            call mpi_wait(ireq_hhs(iblk), MPI_STATUS_IGNORE, mpierr)
            ! Copy vectors into send buffer
            hh_send(:,1:hh_cnt(iblk),iblk) = hh_gath(:,1:hh_cnt(iblk),iblk)
            ! Send to destination
            call mpi_isend(hh_send(1,1,iblk), nb*hh_cnt(iblk), MPI_COMPLEX16, &
                           global_id(hh_dst(iblk),mod(iblk+block_limits(my_pe)-1,np_cols)), &
                           10+iblk, mpi_comm, ireq_hhs(iblk), mpierr)
            ! Reset counter and increase destination row
            hh_cnt(iblk) = 0
            hh_dst(iblk) = hh_dst(iblk)+1
         endif

         ! The following code is structured in a way to keep waiting times for
         ! other PEs at a minimum, especially if there is only one block.
         ! For this reason, it requests the last column as late as possible
         ! and sends the Householder vector and the first column as early
         ! as possible.

         nc = MIN(na-ns-n_off+1,nb) ! number of columns in diagonal block
         nr = MIN(na-nb-ns-n_off+1,nb) ! rows in subdiagonal block (may be < 0!!!)
                                       ! Note that nr>=0 implies that diagonal block is full (nc==nb)!

         ! Multiply diagonal block and subdiagonal block with Householder vector

         if(iblk==nblocks .and. nc==nb) then

            ! We need the last column from the next PE.
            ! First do the matrix multiplications without last column ...

            ! Diagonal block, the contribution of the last element is added below!
            ab(1,ne) = 0
            call ZHEMV('L',nc,tau,ab(1,ns),2*nb-1,hv,1,(0.d0,0.d0),hd,1)

            ! Subdiagonal block
            if(nr>0) call ZGEMV('N',nr,nb-1,tau,ab(nb+1,ns),2*nb-1,hv,1,(0.d0,0.d0),hs,1)

            ! ... then request last column ...
            call mpi_recv(ab(1,ne),nb+1,MPI_COMPLEX16,my_pe+1,1,mpi_comm,MPI_STATUS_IGNORE,mpierr)

            ! ... and complete the result
            hs(1:nr) = hs(1:nr) + ab(2:nr+1,ne)*tau*hv(nb)
            hd(nb) = hd(nb) + ab(1,ne)*hv(nb)*tau

         else

            ! Normal matrix multiply
            call ZHEMV('L',nc,tau,ab(1,ns),2*nb-1,hv,1,(0.d0,0.d0),hd,1)
            if(nr>0) call ZGEMV('N',nr,nb,tau,ab(nb+1,ns),2*nb-1,hv,1,(0.d0,0.d0),hs,1)

         endif

         ! Calculate first column of subdiagonal block and calculate new
         ! Householder transformation for this column

         hv_new(:) = 0 ! Needed, last rows must be 0 for nr < nb
         tau_new = 0

         if(nr>0) then

            ! complete (old) Householder transformation for first column

            ab(nb+1:nb+nr,ns) = ab(nb+1:nb+nr,ns) - hs(1:nr) ! Note: hv(1) == 1

            ! calculate new Householder transformation ...
            if(nr>1) then
               vnorm2 = sum(dble(ab(nb+2:nb+nr,ns))**2+dimag(ab(nb+2:nb+nr,ns))**2)
               call hh_transform_complex(ab(nb+1,ns),vnorm2,hf,tau_new)
               hv_new(1) = 1.
               hv_new(2:nr) = ab(nb+2:nb+nr,ns)*hf
               ab(nb+2:,ns) = 0
            endif

            ! ... and send it away immediatly if this is the last block

            if(iblk==nblocks) then
               call mpi_wait(ireq_hv,MPI_STATUS_IGNORE,mpierr)
               hv_s(1) = tau_new
               hv_s(2:) = hv_new(2:)
               call mpi_isend(hv_s,nb,MPI_COMPLEX16,my_pe+1,2,mpi_comm,ireq_hv,mpierr)
            endif

         endif


         ! Transform diagonal block
         x = dot_product(hv(1:nc),hd(1:nc))*conjg(tau)
         hd(1:nc) = hd(1:nc) - 0.5*x*hv(1:nc)

         if(my_pe>0 .and. iblk==1) then

            ! The first column of the diagonal block has to be send to the previous PE
            ! Calculate first column only ...

            ab(1:nc,ns) = ab(1:nc,ns) - hd(1:nc)*conjg(hv(1)) - hv(1:nc)*conjg(hd(1))

            ! ... send it away ...

            call mpi_wait(ireq_ab,MPI_STATUS_IGNORE,mpierr)
            ab_s(1:nb+1) = ab(1:nb+1,ns)
            call mpi_isend(ab_s,nb+1,MPI_COMPLEX16,my_pe-1,1,mpi_comm,ireq_ab,mpierr)

            ! ... and calculate remaining columns with rank-2 update
            if(nc>1) call ZHER2('L',nc-1,(-1.d0,0.d0),hd(2),1,hv(2),1,ab(1,ns+1),2*nb-1)
         else
            ! No need to  send, just a rank-2 update
            call ZHER2('L',nc,(-1.d0,0.d0),hd,1,hv,1,ab(1,ns),2*nb-1)
         endif

         ! Do the remaining double Householder transformation on the subdiagonal block cols 2 ... nb

         if(nr>0) then
            if(nr>1) then
               call ZGEMV('C',nr,nb-1,tau_new,ab(nb,ns+1),2*nb-1,hv_new,1,(0.d0,0.d0),h(2),1)
               x = dot_product(hs(1:nr),hv_new(1:nr))*tau_new
               h(2:nb) = h(2:nb) - x*hv(2:nb)
               ! Unfortunately the is no BLAS routine like DGER2 for a nonsymmetric rank 2 update
               do i=2,nb
                  ab(2+nb-i:1+nb+nr-i,i+ns-1) = ab(2+nb-i:1+nb+nr-i,i+ns-1) - hv_new(1:nr)*conjg(h(i)) - hs(1:nr)*conjg(hv(i))
               enddo
            else
               ! No double Householder transformation for nr=1, just complete the row
               do i=2,nb
                  ab(2+nb-i,i+ns-1) = ab(2+nb-i,i+ns-1) - hs(1)*conjg(hv(i))
               enddo
            endif
         endif

         ! Use new HH vector for the next block
         hv(:) = hv_new(:)
         tau = tau_new

      enddo

   enddo

   ! Finish the last outstanding requests
   call mpi_wait(ireq_ab,MPI_STATUS_IGNORE,mpierr)
   call mpi_wait(ireq_hv,MPI_STATUS_IGNORE,mpierr)

   call mpi_waitall(nblocks, ireq_hhs, MPI_STATUSES_IGNORE, mpierr)
   call mpi_waitall(num_chunks, ireq_hhr, MPI_STATUSES_IGNORE, mpierr)

   call mpi_barrier(mpi_comm,mpierr)

   deallocate(ab)
   deallocate(ireq_hhr, ireq_hhs)
   deallocate(hh_cnt, hh_dst)
   deallocate(hh_gath, hh_send)
   deallocate(limits, snd_limits)
   deallocate(block_limits)
   deallocate(global_id)

end subroutine

!---------------------------------------------------------------------------------------------------

subroutine trans_ev_tridi_to_band_complex(na, nev, nblk, nbw, q, ldq, mpi_comm_rows, mpi_comm_cols)

!-------------------------------------------------------------------------------
!  trans_ev_tridi_to_band_complex:
!  Transforms the eigenvectors of a tridiagonal matrix back to the eigenvectors of the band matrix
!
!  Parameters
!
!  na          Order of matrix a, number of rows of matrix q
!
!  nev         Number eigenvectors to compute (= columns of matrix q)
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  nb          semi bandwith
!
!  q           On input: Eigenvectors of tridiagonal matrix
!              On output: Transformed eigenvectors
!              Distribution is like in Scalapack.
!
!  ldq         Leading dimension of q
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns/both
!
!-------------------------------------------------------------------------------

    implicit none

    integer, intent(in) :: na, nev, nblk, nbw, ldq, mpi_comm_rows, mpi_comm_cols
    complex*16 q(ldq,*)

    integer np_rows, my_prow, np_cols, my_pcol

    integer i, j, ip, sweep, nbuf, l_nev, a_dim2
    integer current_n, current_local_n, current_n_start, current_n_end
    integer next_n, next_local_n, next_n_start, next_n_end
    integer bottom_msg_length, top_msg_length, next_top_msg_length
    integer thread_width, stripe_width, stripe_count, csw
    integer num_result_blocks, num_result_buffers, num_bufs_recvd
    integer a_off, current_tv_off, max_blk_size, b_off, b_len
    integer mpierr, src, src_offset, dst, offset, nfact, num_blk
    logical flag

    complex*16, allocatable :: a(:,:,:,:), row(:)
    complex*16, allocatable :: top_border_send_buffer(:,:), top_border_recv_buffer(:,:)
    complex*16, allocatable :: bottom_border_send_buffer(:,:), bottom_border_recv_buffer(:,:)
    complex*16, allocatable :: result_buffer(:,:,:)
    complex*16, allocatable :: bcast_buffer(:,:)

    integer n_off
    integer, allocatable :: result_send_request(:), result_recv_request(:), limits(:)
    integer, allocatable :: top_send_request(:), bottom_send_request(:)
    integer, allocatable :: top_recv_request(:), bottom_recv_request(:)

    ! MPI send/recv tags, arbitrary

    integer, parameter :: bottom_recv_tag = 111
    integer, parameter :: top_recv_tag    = 222
    integer, parameter :: result_recv_tag = 333

    integer :: max_threads, my_thread
!$  integer :: omp_get_max_threads

    ! Just for measuring the kernel performance
    real*8 kernel_time
    integer*8 kernel_flops


    kernel_time = 1.d-100
    kernel_flops = 0

    max_threads = 1
!$  max_threads = omp_get_max_threads()

    call MPI_Comm_rank(mpi_comm_rows, my_prow, mpierr)
    call MPI_Comm_size(mpi_comm_rows, np_rows, mpierr)
    call MPI_Comm_rank(mpi_comm_cols, my_pcol, mpierr)
    call MPI_Comm_size(mpi_comm_cols, np_cols, mpierr)

    if(mod(nbw,nblk)/=0) then
      if(my_prow==0 .and. my_pcol==0) then
         print *,'ERROR: nbw=',nbw,', nblk=',nblk
         print *,'band backtransform works only for nbw==n*nblk'
         call mpi_abort(mpi_comm_world,0,mpierr)
      endif
    endif

    nfact = nbw / nblk


    ! local number of eigenvectors
    l_nev = local_index(nev, my_pcol, np_cols, nblk, -1)

    if(l_nev==0) then
        thread_width = 0
        stripe_width = 0
        stripe_count = 0
    else
        ! Suggested stripe width is 48 - should this be reduced for the complex case ???
        thread_width = (l_nev-1)/max_threads + 1 ! number of eigenvectors per OMP thread
        stripe_width = 48 ! Must be a multiple of 4
        stripe_count = (thread_width-1)/stripe_width + 1
        ! Adapt stripe width so that last one doesn't get too small
        stripe_width = (thread_width-1)/stripe_count + 1
        stripe_width = ((stripe_width+3)/4)*4 ! Must be a multiple of 4 !!!
    endif

    ! Determine the matrix distribution at the beginning

    allocate(limits(0:np_rows))

    call determine_workload(na, nbw, np_rows, limits)

    max_blk_size = maxval(limits(1:np_rows) - limits(0:np_rows-1))

    a_dim2 = max_blk_size + nbw

    allocate(a(stripe_width,a_dim2,stripe_count,max_threads))
    ! a(:,:,:,:) should be set to 0 in a parallel region, not here!

    allocate(row(l_nev))
    row(:) = 0

    ! Copy q from a block cyclic distribution into a distribution with contiguous rows,
    ! and transpose the matrix using stripes of given stripe_width for cache blocking.

    ! The peculiar way it is done below is due to the fact that the last row should be
    ! ready first since it is the first one to start below

    ! Please note about the OMP usage below:
    ! This is not for speed, but because we want the matrix a in the memory and
    ! in the cache of the correct thread (if possible)

!$omp parallel do private(my_thread), schedule(static, 1)
    do my_thread = 1, max_threads
        a(:,:,:,my_thread) = 0 ! if possible, do first touch allocation!
    enddo

    do ip = np_rows-1, 0, -1
        if(my_prow == ip) then
            ! Receive my rows which have not yet been received
            src_offset = local_index(limits(ip), my_prow, np_rows, nblk, -1)
            do i=limits(ip)+1,limits(ip+1)
                src = mod((i-1)/nblk, np_rows)
                if(src < my_prow) then
                    call MPI_Recv(row, l_nev, MPI_COMPLEX16, src, 0, mpi_comm_rows, MPI_STATUS_IGNORE, mpierr)
!$omp parallel do private(my_thread), schedule(static, 1)
                    do my_thread = 1, max_threads
                        call unpack_row(row,i-limits(ip),my_thread)
                    enddo
                elseif(src==my_prow) then
                    src_offset = src_offset+1
                    row(:) = q(src_offset, 1:l_nev)
!$omp parallel do private(my_thread), schedule(static, 1)
                    do my_thread = 1, max_threads
                        call unpack_row(row,i-limits(ip),my_thread)
                    enddo
                endif
            enddo
            ! Send all rows which have not yet been send
            src_offset = 0
            do dst = 0, ip-1
              do i=limits(dst)+1,limits(dst+1)
                if(mod((i-1)/nblk, np_rows) == my_prow) then
                    src_offset = src_offset+1
                    row(:) = q(src_offset, 1:l_nev)
                    call MPI_Send(row, l_nev, MPI_COMPLEX16, dst, 0, mpi_comm_rows, mpierr)
                endif
              enddo
            enddo
        else if(my_prow < ip) then
            ! Send all rows going to PE ip
            src_offset = local_index(limits(ip), my_prow, np_rows, nblk, -1)
            do i=limits(ip)+1,limits(ip+1)
                src = mod((i-1)/nblk, np_rows)
                if(src == my_prow) then
                    src_offset = src_offset+1
                    row(:) = q(src_offset, 1:l_nev)
                    call MPI_Send(row, l_nev, MPI_COMPLEX16, ip, 0, mpi_comm_rows, mpierr)
                endif
            enddo
            ! Receive all rows from PE ip
            do i=limits(my_prow)+1,limits(my_prow+1)
                src = mod((i-1)/nblk, np_rows)
                if(src == ip) then
                    call MPI_Recv(row, l_nev, MPI_COMPLEX16, src, 0, mpi_comm_rows, MPI_STATUS_IGNORE, mpierr)
!$omp parallel do private(my_thread), schedule(static, 1)
                    do my_thread = 1, max_threads
                        call unpack_row(row,i-limits(my_prow),my_thread)
                    enddo
                endif
            enddo
        endif
    enddo


    ! Set up result buffer queue

    num_result_blocks = ((na-1)/nblk + np_rows - my_prow) / np_rows

    num_result_buffers = 4*nfact
    allocate(result_buffer(l_nev,nblk,num_result_buffers))

    allocate(result_send_request(num_result_buffers))
    allocate(result_recv_request(num_result_buffers))
    result_send_request(:) = MPI_REQUEST_NULL
    result_recv_request(:) = MPI_REQUEST_NULL

    ! Queue up buffers

    if(my_prow > 0 .and. l_nev>0) then ! note: row 0 always sends
        do j = 1, min(num_result_buffers, num_result_blocks)
            call MPI_Irecv(result_buffer(1,1,j), l_nev*nblk, MPI_COMPLEX16, 0, result_recv_tag, &
                           mpi_comm_rows, result_recv_request(j), mpierr)
        enddo
    endif

    num_bufs_recvd = 0 ! No buffers received yet

    ! Initialize top/bottom requests

    allocate(top_send_request(stripe_count))
    allocate(top_recv_request(stripe_count))
    allocate(bottom_send_request(stripe_count))
    allocate(bottom_recv_request(stripe_count))

    top_send_request(:) = MPI_REQUEST_NULL
    top_recv_request(:) = MPI_REQUEST_NULL
    bottom_send_request(:) = MPI_REQUEST_NULL
    bottom_recv_request(:) = MPI_REQUEST_NULL

    allocate(top_border_send_buffer(stripe_width*nbw*max_threads, stripe_count))
    allocate(top_border_recv_buffer(stripe_width*nbw*max_threads, stripe_count))
    allocate(bottom_border_send_buffer(stripe_width*nbw*max_threads, stripe_count))
    allocate(bottom_border_recv_buffer(stripe_width*nbw*max_threads, stripe_count))

    top_border_send_buffer(:,:) = 0
    top_border_recv_buffer(:,:) = 0
    bottom_border_send_buffer(:,:) = 0
    bottom_border_recv_buffer(:,:) = 0

    ! Initialize broadcast buffer

    allocate(bcast_buffer(nbw, max_blk_size))
    bcast_buffer = 0

    current_tv_off = 0 ! Offset of next row to be broadcast


    ! ------------------- start of work loop -------------------

    a_off = 0 ! offset in A (to avoid unnecessary shifts)

    top_msg_length = 0
    bottom_msg_length = 0

    do sweep = 0, (na-1)/nbw

        current_n = na - sweep*nbw
        call determine_workload(current_n, nbw, np_rows, limits)
        current_n_start = limits(my_prow)
        current_n_end   = limits(my_prow+1)
        current_local_n = current_n_end - current_n_start

        next_n = max(current_n - nbw, 0)
        call determine_workload(next_n, nbw, np_rows, limits)
        next_n_start = limits(my_prow)
        next_n_end   = limits(my_prow+1)
        next_local_n = next_n_end - next_n_start

        if(next_n_end < next_n) then
            bottom_msg_length = current_n_end - next_n_end
        else
            bottom_msg_length = 0
        endif

        if(next_local_n > 0) then
            next_top_msg_length = current_n_start - next_n_start
        else
            next_top_msg_length = 0
        endif

        if(sweep==0 .and. current_n_end < current_n .and. l_nev > 0) then
            do i = 1, stripe_count
                csw = min(stripe_width, thread_width-(i-1)*stripe_width) ! "current_stripe_width"
                b_len = csw*nbw*max_threads
                call MPI_Irecv(bottom_border_recv_buffer(1,i), b_len, MPI_COMPLEX16, my_prow+1, bottom_recv_tag, &
                           mpi_comm_rows, bottom_recv_request(i), mpierr)
            enddo
        endif

        if(current_local_n > 1) then
            if(my_pcol == mod(sweep,np_cols)) then
                bcast_buffer(:,1:current_local_n) = hh_trans_complex(:,current_tv_off+1:current_tv_off+current_local_n)
                current_tv_off = current_tv_off + current_local_n
            endif
            call mpi_bcast(bcast_buffer, nbw*current_local_n, MPI_COMPLEX16, mod(sweep,np_cols), mpi_comm_cols, mpierr)
        else
            ! for current_local_n == 1 the one and only HH vector is 0 and not stored in hh_trans_complex
            bcast_buffer(:,1) = 0
        endif

        if(l_nev == 0) cycle

        if(current_local_n > 0) then

          do i = 1, stripe_count

            ! Get real stripe width for strip i;
            ! The last OpenMP tasks may have an even smaller stripe with,
            ! but we don't care about this, i.e. we send/recv a bit too much in this case.
            ! csw: current_stripe_width

            csw = min(stripe_width, thread_width-(i-1)*stripe_width)

            !wait_b
            if(current_n_end < current_n) then
                call MPI_Wait(bottom_recv_request(i), MPI_STATUS_IGNORE, mpierr)
!$omp parallel do private(my_thread, n_off, b_len, b_off), schedule(static, 1)
                do my_thread = 1, max_threads
                    n_off = current_local_n+a_off
                    b_len = csw*nbw
                    b_off = (my_thread-1)*b_len
                    a(1:csw,n_off+1:n_off+nbw,i,my_thread) = &
                      reshape(bottom_border_recv_buffer(b_off+1:b_off+b_len,i), (/ csw, nbw /))
                enddo
                if(next_n_end < next_n) then
                    call MPI_Irecv(bottom_border_recv_buffer(1,i), csw*nbw*max_threads, &
                                   MPI_COMPLEX16, my_prow+1, bottom_recv_tag, &
                                   mpi_comm_rows, bottom_recv_request(i), mpierr)
                endif
            endif

            if(current_local_n <= bottom_msg_length + top_msg_length) then

                !wait_t
                if(top_msg_length>0) then
                    call MPI_Wait(top_recv_request(i), MPI_STATUS_IGNORE, mpierr)
                endif

                !compute
!$omp parallel do private(my_thread, n_off, b_len, b_off), schedule(static, 1)
                do my_thread = 1, max_threads
                    if(top_msg_length>0) then
                        b_len = csw*top_msg_length
                        b_off = (my_thread-1)*b_len
                        a(1:csw,a_off+1:a_off+top_msg_length,i,my_thread) = &
                          reshape(top_border_recv_buffer(b_off+1:b_off+b_len,i), (/ csw, top_msg_length /))
                    endif
                    call compute_hh_trafo(0, current_local_n, i, my_thread)
                enddo

                !send_b
                call MPI_Wait(bottom_send_request(i), MPI_STATUS_IGNORE, mpierr)
                if(bottom_msg_length>0) then
                    n_off = current_local_n+nbw-bottom_msg_length+a_off
                    b_len = csw*bottom_msg_length*max_threads
                    bottom_border_send_buffer(1:b_len,i) = &
                        reshape(a(1:csw,n_off+1:n_off+bottom_msg_length,i,:), (/ b_len /))
                    call MPI_Isend(bottom_border_send_buffer(1,i), b_len, MPI_COMPLEX16, my_prow+1, &
                                   top_recv_tag, mpi_comm_rows, bottom_send_request(i), mpierr)
                endif

            else

                !compute
!$omp parallel do private(my_thread, b_len, b_off), schedule(static, 1)
                do my_thread = 1, max_threads
                    call compute_hh_trafo(current_local_n - bottom_msg_length, bottom_msg_length, i, my_thread)
                enddo

                !send_b
                call MPI_Wait(bottom_send_request(i), MPI_STATUS_IGNORE, mpierr)
                if(bottom_msg_length > 0) then
                    n_off = current_local_n+nbw-bottom_msg_length+a_off
                    b_len = csw*bottom_msg_length*max_threads
                    bottom_border_send_buffer(1:b_len,i) = &
                      reshape(a(1:csw,n_off+1:n_off+bottom_msg_length,i,:), (/ b_len /))
                    call MPI_Isend(bottom_border_send_buffer(1,i), b_len, MPI_COMPLEX16, my_prow+1, &
                                   top_recv_tag, mpi_comm_rows, bottom_send_request(i), mpierr)
                endif

                !compute
!$omp parallel do private(my_thread), schedule(static, 1)
                do my_thread = 1, max_threads
                    call compute_hh_trafo(top_msg_length, current_local_n-top_msg_length-bottom_msg_length, i, my_thread)
                enddo

                !wait_t
                if(top_msg_length>0) then
                    call MPI_Wait(top_recv_request(i), MPI_STATUS_IGNORE, mpierr)
                endif

                !compute
!$omp parallel do private(my_thread, b_len, b_off), schedule(static, 1)
                do my_thread = 1, max_threads
                    if(top_msg_length>0) then
                        b_len = csw*top_msg_length
                        b_off = (my_thread-1)*b_len
                        a(1:csw,a_off+1:a_off+top_msg_length,i,my_thread) = &
                          reshape(top_border_recv_buffer(b_off+1:b_off+b_len,i), (/ csw, top_msg_length /))
                    endif
                    call compute_hh_trafo(0, top_msg_length, i, my_thread)
                enddo
            endif

            if(next_top_msg_length > 0) then
                !request top_border data
                b_len = csw*next_top_msg_length*max_threads
                call MPI_Irecv(top_border_recv_buffer(1,i), b_len, MPI_COMPLEX16, my_prow-1, &
                               top_recv_tag, mpi_comm_rows, top_recv_request(i), mpierr)
            endif

            !send_t
            if(my_prow > 0) then
                call MPI_Wait(top_send_request(i), MPI_STATUS_IGNORE, mpierr)
                b_len = csw*nbw*max_threads
                top_border_send_buffer(1:b_len,i) = reshape(a(1:csw,a_off+1:a_off+nbw,i,:), (/ b_len /))
                call MPI_Isend(top_border_send_buffer(1,i), b_len, MPI_COMPLEX16, &
                               my_prow-1, bottom_recv_tag, &
                               mpi_comm_rows, top_send_request(i), mpierr)
            endif

            ! Care that there are not too many outstanding top_recv_request's
            if(stripe_count > 1) then
                if(i>1) then
                    call MPI_Wait(top_recv_request(i-1), MPI_STATUS_IGNORE, mpierr)
                else
                    call MPI_Wait(top_recv_request(stripe_count), MPI_STATUS_IGNORE, mpierr)
                endif
            endif

          enddo

          top_msg_length = next_top_msg_length

        else
            ! wait for last top_send_request
          do i = 1, stripe_count
            call MPI_Wait(top_send_request(i), MPI_STATUS_IGNORE, mpierr)
          enddo
        endif

        ! Care about the result

        if(my_prow == 0) then

            ! topmost process sends nbw rows to destination processes

            do j=0,nfact-1

                num_blk = sweep*nfact+j ! global number of destination block, 0 based
                if(num_blk*nblk >= na) exit

                nbuf = mod(num_blk, num_result_buffers) + 1 ! buffer number to get this block

                call MPI_Wait(result_send_request(nbuf), MPI_STATUS_IGNORE, mpierr)

                dst = mod(num_blk, np_rows)

                if(dst == 0) then
                    do i = 1, min(na - num_blk*nblk, nblk)
                        call pack_row(row, j*nblk+i+a_off)
                        q((num_blk/np_rows)*nblk+i,1:l_nev) = row(:)
                    enddo
                else
                    do i = 1, nblk
                        call pack_row(result_buffer(:,i,nbuf),j*nblk+i+a_off)
                    enddo
                    call MPI_Isend(result_buffer(1,1,nbuf), l_nev*nblk, MPI_COMPLEX16, dst, &
                                   result_recv_tag, mpi_comm_rows, result_send_request(nbuf), mpierr)
                endif
            enddo

        else

           ! receive and store final result

            do j = num_bufs_recvd, num_result_blocks-1

                nbuf = mod(j, num_result_buffers) + 1 ! buffer number to get this block

                ! If there is still work to do, just test for the next result request
                ! and leave the loop if it is not ready, otherwise wait for all
                ! outstanding requests

                if(next_local_n > 0) then
                    call MPI_Test(result_recv_request(nbuf), flag, MPI_STATUS_IGNORE, mpierr)
                    if(.not.flag) exit
                else
                    call MPI_Wait(result_recv_request(nbuf), MPI_STATUS_IGNORE, mpierr)
                endif

                ! Fill result buffer into q
                num_blk = j*np_rows + my_prow ! global number of current block, 0 based
                do i = 1, min(na - num_blk*nblk, nblk)
                    q(j*nblk+i, 1:l_nev) = result_buffer(1:l_nev, i, nbuf)
                enddo

                ! Queue result buffer again if there are outstanding blocks left
                if(j+num_result_buffers < num_result_blocks) &
                    call MPI_Irecv(result_buffer(1,1,nbuf), l_nev*nblk, MPI_COMPLEX16, 0, result_recv_tag, &
                                   mpi_comm_rows, result_recv_request(nbuf), mpierr)

            enddo
            num_bufs_recvd = j

        endif

        ! Shift the remaining rows to the front of A (if necessary)

        offset = nbw - top_msg_length
        if(offset<0) then
            print *,'internal error, offset for shifting = ',offset
            call MPI_Abort(MPI_COMM_WORLD, 1, mpierr)
        endif
        a_off = a_off + offset
        if(a_off + next_local_n + nbw > a_dim2) then
!$omp parallel do private(my_thread, i, j), schedule(static, 1)
            do my_thread = 1, max_threads
                do i = 1, stripe_count
                    do j = top_msg_length+1, top_msg_length+next_local_n
                       A(:,j,i,my_thread) = A(:,j+a_off,i,my_thread)
                    enddo
                enddo
            enddo
            a_off = 0
        endif

    enddo

    ! Just for safety:
    if(ANY(top_send_request    /= MPI_REQUEST_NULL)) print *,'*** ERROR top_send_request ***',my_prow,my_pcol
    if(ANY(bottom_send_request /= MPI_REQUEST_NULL)) print *,'*** ERROR bottom_send_request ***',my_prow,my_pcol
    if(ANY(top_recv_request    /= MPI_REQUEST_NULL)) print *,'*** ERROR top_recv_request ***',my_prow,my_pcol
    if(ANY(bottom_recv_request /= MPI_REQUEST_NULL)) print *,'*** ERROR bottom_recv_request ***',my_prow,my_pcol

    if(my_prow == 0) then
        call MPI_Waitall(num_result_buffers, result_send_request, MPI_STATUSES_IGNORE, mpierr)
    endif

    if(ANY(result_send_request /= MPI_REQUEST_NULL)) print *,'*** ERROR result_send_request ***',my_prow,my_pcol
    if(ANY(result_recv_request /= MPI_REQUEST_NULL)) print *,'*** ERROR result_recv_request ***',my_prow,my_pcol

    if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) &
        print '(" Kernel time:",f10.3," MFlops: ",f10.3)', kernel_time, kernel_flops/kernel_time*1.d-6

    ! deallocate all working space

    deallocate(a)
    deallocate(row)
    deallocate(limits)
    deallocate(result_send_request)
    deallocate(result_recv_request)
    deallocate(top_border_send_buffer)
    deallocate(top_border_recv_buffer)
    deallocate(bottom_border_send_buffer)
    deallocate(bottom_border_recv_buffer)
    deallocate(result_buffer)
    deallocate(bcast_buffer)
    deallocate(top_send_request)
    deallocate(top_recv_request)
    deallocate(bottom_send_request)
    deallocate(bottom_recv_request)

contains

    subroutine pack_row(row, n)
        complex*16 row(:)
        integer n, i, noff, nl, nt

        do nt = 1, max_threads
            do i = 1, stripe_count
                noff = (nt-1)*thread_width + (i-1)*stripe_width
                nl   = min(stripe_width, nt*thread_width-noff, l_nev-noff)
                if(nl<=0) exit
                row(noff+1:noff+nl) = a(1:nl,n,i,nt)
            enddo
        enddo

    end subroutine

    subroutine unpack_row(row, n, my_thread)

        ! Private variables in OMP regions (my_thread) should better be in the argument list!
        integer, intent(in) :: n, my_thread
        complex*16, intent(in)  :: row(:)
        integer i, noff, nl

        do i=1,stripe_count
            noff = (my_thread-1)*thread_width + (i-1)*stripe_width
            nl   = min(stripe_width, my_thread*thread_width-noff, l_nev-noff)
            if(nl<=0) exit
            a(1:nl,n,i,my_thread) = row(noff+1:noff+nl)
        enddo

    end subroutine

    subroutine compute_hh_trafo(off, ncols, istripe, my_thread)

        ! Private variables in OMP regions (my_thread) should better be in the argument list!
        integer, intent(in) :: off, ncols, istripe, my_thread
        integer j, nl, noff
        real*8 ttt

        ttt = mpi_wtime()
        if(istripe<stripe_count) then
          nl = stripe_width
        else
          noff = (my_thread-1)*thread_width + (istripe-1)*stripe_width
          nl = min(my_thread*thread_width-noff, l_nev-noff)
          if(nl<=0) return
        endif
        do j = ncols, 1, -1
          call single_hh_trafo_complex(a(1,j+off+a_off,istripe,my_thread),bcast_buffer(1,j+off),nbw,nl,stripe_width)
        enddo
        if(my_thread==1) then
          kernel_flops = kernel_flops + 4*4*int(nl,8)*int(ncols,8)*int(nbw,8)
          kernel_time  = kernel_time + mpi_wtime()-ttt
        endif

    end subroutine

end subroutine

! --------------------------------------------------------------------------------------------------
! redist_band: redistributes band from 2D block cyclic form to 1D band

subroutine redist_band(l_real, r_a, c_a, lda, na, nblk, nbw, mpi_comm_rows, mpi_comm_cols, mpi_comm, r_ab, c_ab)

   logical, intent(in)     :: l_real
   real*8, intent(in)      :: r_a(lda, *)
   complex*16, intent(in)  :: c_a(lda, *)
   integer, intent(in)     :: lda, na, nblk, nbw, mpi_comm_rows, mpi_comm_cols, mpi_comm
   real*8, intent(out)     :: r_ab(:,:)
   complex*16, intent(out) :: c_ab(:,:)

   integer, allocatable :: ncnt_s(:), nstart_s(:), ncnt_r(:), nstart_r(:), global_id(:,:), block_limits(:)
   real*8, allocatable :: r_sbuf(:,:,:), r_rbuf(:,:,:), r_buf(:,:)
   complex*16, allocatable :: c_sbuf(:,:,:), c_rbuf(:,:,:), c_buf(:,:)

   integer i, j, my_pe, n_pes, my_prow, np_rows, my_pcol, np_cols, nfact, np, npr, npc, mpierr, is, js
   integer nblocks_total, il, jl, l_rows, l_cols, n_off

   call mpi_comm_rank(mpi_comm,my_pe,mpierr)
   call mpi_comm_size(mpi_comm,n_pes,mpierr)

   call mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   call mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   call mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   call mpi_comm_size(mpi_comm_cols,np_cols,mpierr)

   ! Get global_id mapping 2D procssor coordinates to global id

   allocate(global_id(0:np_rows-1,0:np_cols-1))
   global_id(:,:) = 0
   global_id(my_prow, my_pcol) = my_pe

   call mpi_allreduce(mpi_in_place, global_id, np_rows*np_cols, mpi_integer, mpi_sum, mpi_comm, mpierr)


   ! Set work distribution

   nblocks_total = (na-1)/nbw + 1

   allocate(block_limits(0:n_pes))
   call divide_band(nblocks_total, n_pes, block_limits)


   allocate(ncnt_s(0:n_pes-1))
   allocate(nstart_s(0:n_pes-1))
   allocate(ncnt_r(0:n_pes-1))
   allocate(nstart_r(0:n_pes-1))


   nfact = nbw/nblk

   ! Count how many blocks go to which PE

   ncnt_s(:) = 0
   np = 0 ! receiver PE number
   do j=0,(na-1)/nblk ! loop over rows of blocks
      if(j/nfact==block_limits(np+1)) np = np+1
      if(mod(j,np_rows) == my_prow) then
         do i=0,nfact
            if(mod(i+j,np_cols) == my_pcol) then
               ncnt_s(np) = ncnt_s(np) + 1
            endif
         enddo
      endif
   enddo

   ! Allocate send buffer

   if(l_real) then
      allocate(r_sbuf(nblk,nblk,sum(ncnt_s)))
      r_sbuf(:,:,:) = 0.
   else
      allocate(c_sbuf(nblk,nblk,sum(ncnt_s)))
      c_sbuf(:,:,:) = 0.
   endif

   ! Determine start offsets in send buffer

   nstart_s(0) = 0
   do i=1,n_pes-1
      nstart_s(i) = nstart_s(i-1) + ncnt_s(i-1)
   enddo

   ! Fill send buffer

   l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a
   l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local columns of a

   np = 0
   do j=0,(na-1)/nblk ! loop over rows of blocks
      if(j/nfact==block_limits(np+1)) np = np+1
      if(mod(j,np_rows) == my_prow) then
         do i=0,nfact
            if(mod(i+j,np_cols) == my_pcol) then
               nstart_s(np) = nstart_s(np) + 1
               js = (j/np_rows)*nblk
               is = ((i+j)/np_cols)*nblk
               jl = MIN(nblk,l_rows-js)
               il = MIN(nblk,l_cols-is)
               if(l_real) then
                  r_sbuf(1:jl,1:il,nstart_s(np)) = r_a(js+1:js+jl,is+1:is+il)
               else
                  c_sbuf(1:jl,1:il,nstart_s(np)) = c_a(js+1:js+jl,is+1:is+il)
               endif
            endif
         enddo
      endif
   enddo

   ! Count how many blocks we get from which PE

   ncnt_r(:) = 0
   do j=block_limits(my_pe)*nfact,min(block_limits(my_pe+1)*nfact-1,(na-1)/nblk)
      npr = mod(j,np_rows)
      do i=0,nfact
         npc = mod(i+j,np_cols)
         np = global_id(npr,npc)
         ncnt_r(np) = ncnt_r(np) + 1
      enddo
   enddo

   ! Allocate receive buffer

   if(l_real) then
      allocate(r_rbuf(nblk,nblk,sum(ncnt_r)))
   else
      allocate(c_rbuf(nblk,nblk,sum(ncnt_r)))
   endif

   ! Set send counts/send offsets, receive counts/receive offsets
   ! now actually in variables, not in blocks

   ncnt_s(:) = ncnt_s(:)*nblk*nblk

   nstart_s(0) = 0
   do i=1,n_pes-1
      nstart_s(i) = nstart_s(i-1) + ncnt_s(i-1)
   enddo

   ncnt_r(:) = ncnt_r(:)*nblk*nblk

   nstart_r(0) = 0
   do i=1,n_pes-1
      nstart_r(i) = nstart_r(i-1) + ncnt_r(i-1)
   enddo

   ! Exchange all data with MPI_Alltoallv

   if(l_real) then
      call MPI_Alltoallv(r_sbuf,ncnt_s,nstart_s,MPI_REAL8,r_rbuf,ncnt_r,nstart_r,MPI_REAL8,mpi_comm,mpierr)
   else
      call MPI_Alltoallv(c_sbuf,ncnt_s,nstart_s,MPI_COMPLEX16,c_rbuf,ncnt_r,nstart_r,MPI_COMPLEX16,mpi_comm,mpierr)
   endif

   ! set band from receive buffer

   ncnt_r(:) = ncnt_r(:)/(nblk*nblk)

   nstart_r(0) = 0
   do i=1,n_pes-1
      nstart_r(i) = nstart_r(i-1) + ncnt_r(i-1)
   enddo

   if(l_real) then
      allocate(r_buf((nfact+1)*nblk,nblk))
   else
      allocate(c_buf((nfact+1)*nblk,nblk))
   endif

   ! n_off: Offset of ab within band
   n_off = block_limits(my_pe)*nbw

   do j=block_limits(my_pe)*nfact,min(block_limits(my_pe+1)*nfact-1,(na-1)/nblk)
      npr = mod(j,np_rows)
      do i=0,nfact
         npc = mod(i+j,np_cols)
         np = global_id(npr,npc)
         nstart_r(np) = nstart_r(np) + 1
         if(l_real) then
            r_buf(i*nblk+1:i*nblk+nblk,:) = transpose(r_rbuf(:,:,nstart_r(np)))
         else
            c_buf(i*nblk+1:i*nblk+nblk,:) = conjg(transpose(c_rbuf(:,:,nstart_r(np))))
         endif
      enddo
      do i=1,MIN(nblk,na-j*nblk)
         if(l_real) then
            r_ab(1:nbw+1,i+j*nblk-n_off) = r_buf(i:i+nbw,i)
         else
            c_ab(1:nbw+1,i+j*nblk-n_off) = c_buf(i:i+nbw,i)
         endif
      enddo
   enddo

   deallocate(ncnt_s, nstart_s)
   deallocate(ncnt_r, nstart_r)
   deallocate(global_id)
   deallocate(block_limits)
   if(l_real) then
      deallocate(r_sbuf, r_rbuf, r_buf)
   else
      deallocate(c_sbuf, c_rbuf, c_buf)
   endif

end subroutine

!---------------------------------------------------------------------------------------------------
! divide_band: sets the work distribution in band
! Proc n works on blocks block_limits(n)+1 .. block_limits(n+1)

subroutine divide_band(nblocks_total, n_pes, block_limits)

   integer, intent(in) :: nblocks_total ! total number of blocks in band
   integer, intent(in) :: n_pes         ! number of PEs for division
   integer, intent(out) :: block_limits(0:n_pes)

   integer :: n, nblocks, nblocks_left

   block_limits(0) = 0
   if(nblocks_total < n_pes) then
      ! Not enough work for all: The first tasks get exactly 1 block
      do n=1,n_pes
         block_limits(n) = min(nblocks_total,n)
      enddo
   else
      ! Enough work for all. If there is no exact loadbalance,
      ! the LAST tasks get more work since they are finishing earlier!
      nblocks = nblocks_total/n_pes
      nblocks_left = nblocks_total - n_pes*nblocks
      do n=1,n_pes
         if(n<=n_pes-nblocks_left) then
            block_limits(n) = block_limits(n-1) + nblocks
         else
            block_limits(n) = block_limits(n-1) + nblocks + 1
         endif
      enddo
   endif

end subroutine

!-------------------------------------------------------------------------------

subroutine band_band_real(na, nb, nb2, ab, ab2, d, e, mpi_comm)

!-------------------------------------------------------------------------------
! band_band_real:
! Reduces a real symmetric banded matrix to a real symmetric matrix with smaller bandwidth. Householder transformations are not stored.
! Matrix size na and original bandwidth nb have to be a multiple of the target bandwidth nb2. (Hint: expand your matrix with zero entries, if this 
! requirement doesn't hold)
!
!  na          Order of matrix
!
!  nb          Semi bandwidth of original matrix
!
!  nb2         Semi bandwidth of target matrix
!
!  ab          Input matrix with bandwidth nb. The leading dimension of the banded matrix has to be 2*nb. The parallel data layout 
!              has to be accordant to divide_band(), i.e. the matrix columns block_limits(n)*nb+1 to min(na, block_limits(n+1)*nb) 
!              are located on rank n.
!
!  ab2         Output matrix with bandwidth nb2. The leading dimension of the banded matrix is 2*nb2. The parallel data layout is
!              accordant to divide_band(), i.e. the matrix columns block_limits(n)*nb2+1 to min(na, block_limits(n+1)*nb2) are located
!              on rank n.
!
!  d(na)       Diagonal of tridiagonal matrix, set only on PE 0, set only if ab2 = 1 (output)
!
!  e(na)       Subdiagonal of tridiagonal matrix, set only on PE 0, set only if ab2 = 1 (output)
!
!  mpi_comm
!              MPI-Communicator for the total processor set
!-------------------------------------------------------------------------------

   implicit none

   integer, intent(in) ::  na, nb, nb2, mpi_comm
   real*8, intent(inout)  :: ab(2*nb,*)
   real*8, intent(inout)  :: ab2(2*nb2,*)
   real*8, intent(out) :: d(na), e(na) ! set only on PE 0

!----------------

   real*8 hv(nb,nb2), w(nb,nb2), w_new(nb,nb2), tau(nb2), hv_new(nb,nb2), tau_new(nb2), ab_s(1+nb,nb2), ab_r(1+nb,nb2), ab_s2(2*nb2,nb2), hv_s(nb,nb2)
   
   real*8 work(nb*nb2), work2(nb2*nb2)
   integer lwork, info
   
   integer istep, i, n, dest
   integer n_off, na_s
   integer my_pe, n_pes, mpierr
   integer nblocks_total, nblocks
   integer nblocks_total2, nblocks2
   integer ireq_ab, ireq_hv
   integer, allocatable :: block_limits(:), block_limits2(:), ireq_ab2(:)
   
!----------------

   integer j, nc, nr, ns, ne, iblk
   
   call mpi_comm_rank(mpi_comm,my_pe,mpierr)
   call mpi_comm_size(mpi_comm,n_pes,mpierr)

   ! Total number of blocks in the band:
   nblocks_total = (na-1)/nb + 1
   nblocks_total2 = (na-1)/nb2 + 1
   
   ! Set work distribution
   allocate(block_limits(0:n_pes))
   call divide_band(nblocks_total, n_pes, block_limits)
   
   allocate(block_limits2(0:n_pes))
   call divide_band(nblocks_total2, n_pes, block_limits2)

   ! nblocks: the number of blocks for my task
   nblocks = block_limits(my_pe+1) - block_limits(my_pe)
   nblocks2 = block_limits2(my_pe+1) - block_limits2(my_pe)
   
   allocate(ireq_ab2(1:nblocks2))
   ireq_ab2 = MPI_REQUEST_NULL
   if(nb2>1) then
       do i=0,nblocks2-1
           call mpi_irecv(ab2(1,i*nb2+1),2*nb2*nb2,mpi_real8,0,3,mpi_comm,ireq_ab2(i+1),mpierr)
       enddo
   endif

   ! n_off: Offset of ab within band
   n_off = block_limits(my_pe)*nb
   lwork = nb*nb2
   dest = 0

   ireq_ab = MPI_REQUEST_NULL
   ireq_hv = MPI_REQUEST_NULL
   
   ! ---------------------------------------------------------------------------
   ! Start of calculations

   na_s = block_limits(my_pe)*nb + 1

   if(my_pe>0 .and. na_s<=na) then
      ! send first nb2 columns to previous PE
      ! Only the PE owning the diagonal does that (sending 1 element of the subdiagonal block also)
      do i=1,nb2
      	ab_s(1:nb+1,i) = ab(1:nb+1,na_s-n_off+i-1)
      enddo
      call mpi_isend(ab_s,(nb+1)*nb2,mpi_real8,my_pe-1,1,mpi_comm,ireq_ab,mpierr)
   endif
   
   do istep=1,na/nb2
   
      if(my_pe==0) then
      
         n = MIN(na-na_s-nb2+1,nb) ! number of rows to be reduced
         hv(:,:) = 0
         tau(:) = 0
         
         ! The last step (istep=na-1) is only needed for sending the last HH vectors.
         ! We don't want the sign of the last element flipped (analogous to the other sweeps)
         if(istep < na/nb2) then
            
            ! Transform first block column of remaining matrix
            call dgeqrf(n, nb2, ab(1+nb2,na_s-n_off), 2*nb-1, tau, work, lwork, info);
                        
            do i=1,nb2
            	hv(i,i) = 1.0
            	hv(i+1:n,i) = ab(1+nb2+1:1+nb2+n-i,na_s-n_off+i-1)
            	ab(1+nb2+1:2*nb,na_s-n_off+i-1) = 0
            enddo
            
         endif
         
         if(nb2==1) then
            d(istep) = ab(1,na_s-n_off)
	    e(istep) = ab(2,na_s-n_off)
	    if(istep == na) then
	    	e(na) = 0
            endif
         else
            ab_s2 = 0
            ab_s2(:,:) = ab(1:nb2+1,na_s-n_off:na_s-n_off+nb2-1)
            if(block_limits2(dest+1)<istep) then
            	dest = dest+1
            endif
            call mpi_send(ab_s2,2*nb2*nb2,mpi_real8,dest,3,mpi_comm,mpierr)
         endif
         
      else
         if(na>na_s+nb2-1) then
            ! Receive Householder vectors from previous task, from PE owning subdiagonal
            call mpi_recv(hv,nb*nb2,mpi_real8,my_pe-1,2,mpi_comm,MPI_STATUS_IGNORE,mpierr)
            do i=1,nb2
	       	tau(i) = hv(i,i)
	       	hv(i,i) = 1.
            enddo
         endif
      endif
      
      na_s = na_s+nb2
      if(na_s-n_off > nb) then
         ab(:,1:nblocks*nb) = ab(:,nb+1:(nblocks+1)*nb)
         ab(:,nblocks*nb+1:(nblocks+1)*nb) = 0
         n_off = n_off + nb
      endif
      
      do iblk=1,nblocks
         ns = na_s + (iblk-1)*nb - n_off ! first column in block
         ne = ns+nb-nb2                    ! last column in block

         if(ns+n_off>na) exit
         
         nc = MIN(na-ns-n_off+1,nb) ! number of columns in diagonal block
         nr = MIN(na-nb-ns-n_off+1,nb) ! rows in subdiagonal block (may be < 0!!!)
                                       ! Note that nr>=0 implies that diagonal block is full (nc==nb)!
                                       
         call wy_gen(nc,nb2,w,hv,tau,work,nb)

         if(iblk==nblocks .and. nc==nb) then
             !request last nb2 columns
             call mpi_recv(ab_r,(nb+1)*nb2,mpi_real8,my_pe+1,1,mpi_comm,MPI_STATUS_IGNORE,mpierr)
             do i=1,nb2
	         ab(1:nb+1,ne+i-1) = ab_r(:,i)
             enddo
         endif
         
         hv_new(:,:) = 0 ! Needed, last rows must be 0 for nr < nb
         tau_new(:) = 0
         
         if(nr>0) then
             call wy_right(nr,nb,nb2,ab(nb+1,ns),2*nb-1,w,hv,work,nb)
             
             call dgeqrf(nr,nb2,ab(nb+1,ns),2*nb-1,tau_new,work,lwork,info);
                          
             do i=1,nb2
	     	 hv_new(i,i) = 1.0
	     	 hv_new(i+1:,i) = ab(nb+2:2*nb-i+1,ns+i-1)
	     	 ab(nb+2:,ns+i-1) = 0
	     enddo
	     
	     !send hh-vector
	     if(iblk==nblocks) then
	         call mpi_wait(ireq_hv,MPI_STATUS_IGNORE,mpierr)
	         hv_s = hv_new
	         do i=1,nb2
		     hv_s(i,i) = tau_new(i)
                 enddo
	         call mpi_isend(hv_s,nb*nb2,mpi_real8,my_pe+1,2,mpi_comm,ireq_hv,mpierr)
             endif
             
         endif
         
	 call wy_symm(nc,nb2,ab(1,ns),2*nb-1,w,hv,work,work2,nb)
         
         if(my_pe>0 .and. iblk==1) then
	     !send first nb2 columns to previous PE
	     call mpi_wait(ireq_ab,MPI_STATUS_IGNORE,mpierr)
	     do i=1,nb2
	         ab_s(1:nb+1,i) = ab(1:nb+1,ns+i-1)
	     enddo
	     call mpi_isend(ab_s,(nb+1)*nb2,mpi_real8,my_pe-1,1,mpi_comm,ireq_ab,mpierr)
         endif
         
         if(nr>0) then
             call wy_gen(nr,nb2,w_new,hv_new,tau_new,work,nb)
	     call wy_left(nb-nb2,nr,nb2,ab(nb+1-nb2,ns+nb2),2*nb-1,w_new,hv_new,work,nb)
         endif
         
         ! Use new HH vector for the next block
	 hv(:,:) = hv_new(:,:)
         tau = tau_new
         
     enddo

   enddo

   ! Finish the last outstanding requests
   call mpi_wait(ireq_ab,MPI_STATUS_IGNORE,mpierr)
   call mpi_wait(ireq_hv,MPI_STATUS_IGNORE,mpierr)
   call mpi_waitall(nblocks2,ireq_ab2,MPI_STATUSES_IGNORE,mpierr)

   call mpi_barrier(mpi_comm,mpierr)

   deallocate(block_limits)
   deallocate(block_limits2)
   deallocate(ireq_ab2)

end subroutine

! --------------------------------------------------------------------------------------------------

subroutine wy_gen(n, nb, W, Y, tau, mem, lda)
    
    integer, intent(in) :: n		!length of householder-vectors
    integer, intent(in) :: nb		!number of householder-vectors
    integer, intent(in) :: lda		!leading dimension of Y and W
    real*8, intent(in) :: Y(lda,nb)	!matrix containing nb householder-vectors of length b
    real*8, intent(in) :: tau(nb)	!tau values
    real*8, intent(out) :: W(lda,nb)	!output matrix W
    real*8, intent(in) :: mem(nb)	!memory for a temporary matrix of size nb
    
    integer i
    
    W(1:n,1) = tau(1)*Y(1:n,1)
    do i=2,nb
        W(1:n,i) = tau(i)*Y(1:n,i)
        call DGEMV('T',n,i-1,1.d0,Y,lda,W(1,i),1,0.d0,mem,1)
	call DGEMV('N',n,i-1,-1.d0,W,lda,mem,1,1.d0,W(1,i),1)
    enddo
    
end subroutine

! --------------------------------------------------------------------------------------------------

subroutine wy_left(n, m, nb, A, lda, W, Y, mem, lda2)

    integer, intent(in) :: n		!width of the matrix A
    integer, intent(in) :: m		!length of matrix W and Y
    integer, intent(in) :: nb		!width of matrix W and Y
    integer, intent(in) :: lda		!leading dimension of A
    integer, intent(in) :: lda2		!leading dimension of W and Y
    real*8, intent(inout) :: A(lda,*)	!matrix to be transformed
    real*8, intent(in) :: W(m,nb)	!blocked transformation matrix W
    real*8, intent(in) :: Y(m,nb)	!blocked transformation matrix Y
    real*8, intent(inout) :: mem(n,nb)	!memory for a temporary matrix of size n x nb
    
    call DGEMM('T', 'N', nb, n, m, 1.d0, W, lda2, A, lda, 0.d0, mem, nb)
    call DGEMM('N', 'N', m, n, nb, -1.d0, Y, lda2, mem, nb, 1.d0, A, lda)

end subroutine

! --------------------------------------------------------------------------------------------------

subroutine wy_right(n, m, nb, A, lda, W, Y, mem, lda2)

    integer, intent(in) :: n		!height of the matrix A
    integer, intent(in) :: m		!length of matrix W and Y
    integer, intent(in) :: nb		!width of matrix W and Y
    integer, intent(in) :: lda		!leading dimension of A
    integer, intent(in) :: lda2		!leading dimension of W and Y
    real*8, intent(inout) :: A(lda,*)	!matrix to be transformed
    real*8, intent(in) :: W(m,nb)	!blocked transformation matrix W
    real*8, intent(in) :: Y(m,nb)	!blocked transformation matrix Y
    real*8, intent(inout) :: mem(n,nb)	!memory for a temporary matrix of size n x nb
    
    call DGEMM('N', 'N', n, nb, m, 1.d0, A, lda, W, lda2, 0.d0, mem, n)
    call DGEMM('N', 'T', n, m, nb, -1.d0, mem, n, Y, lda2, 1.d0, A, lda)

end subroutine

! --------------------------------------------------------------------------------------------------

subroutine wy_symm(n, nb, A, lda, W, Y, mem, mem2, lda2)

    integer, intent(in) :: n		!width/heigth of the matrix A; length of matrix W and Y
    integer, intent(in) :: nb		!width of matrix W and Y
    integer, intent(in) :: lda		!leading dimension of A
    integer, intent(in) :: lda2		!leading dimension of W and Y
    real*8, intent(inout) :: A(lda,*)	!matrix to be transformed
    real*8, intent(in) :: W(n,nb)	!blocked transformation matrix W
    real*8, intent(in) :: Y(n,nb)	!blocked transformation matrix Y
    real*8 :: mem(n,nb)			!memory for a temporary matrix of size n x nb
    real*8 :: mem2(nb,nb)		!memory for a temporary matrix of size nb x nb
    
    call DSYMM('L', 'L', n, nb, 1.d0, A, lda, W, lda2, 0.d0, mem, n)
    call DGEMM('T', 'N', nb, nb, n, 1.d0, mem, n, W, lda2, 0.d0, mem2, nb)
    call DGEMM('N', 'N', n, nb, nb, -0.5d0, Y, lda2, mem2, nb, 1.d0, mem, n)
    call DSYR2K('L', 'N', n, nb, -1.d0, Y, lda2, mem, n, 1.d0, A, lda)

end subroutine

! --------------------------------------------------------------------------------------------------

#endif
end module ELPA2
