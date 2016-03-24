  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  ! 
  !---------------------------------------------------------------------------
  subroutine ephwan2blochp ( nmodes, xxq, irvec, ndegen, nrr_q, cuf, epmatf, nbnd, nrr_k )
  !---------------------------------------------------------------------------
  !
  ! even though this is for phonons, I use the same notations
  ! adopted for the electronic case (nmodes->nmodes etc)
  !
  !
#include "f_defs.h"
  USE kinds,         only : DP
  USE epwcom,        only : parallel_k, parallel_q, etf_mem
  USE io_epw,        only : iunepmatwp
  USE elph2,         only : epmatwp
  USE constants_epw, ONLY : twopi, ci, czero
  USE io_global,     ONLY : ionode
  USE io_files,      ONLY : prefix, diropn
#ifdef __PARA 
  USE mp_global,     ONLY : inter_pool_comm, intra_pool_comm, mp_sum
  USE mp_world,      ONLY : world_comm
  USE parallel_include
#endif
  implicit none
  !
  !  input variables
  !
  integer :: nmodes, nrr_q, irvec ( 3, nrr_q), ndegen (nrr_q), nbnd, nrr_k
  ! number of bands (possibly in tyhe optimal subspace)
  ! number of WS points
  ! coordinates of WS points
  ! degeneracy of WS points
  ! n of bands
  ! n of electronic WS points
  complex(kind=DP) :: epmatw ( nbnd, nbnd, nrr_k, nmodes), cuf (nmodes, nmodes)
  ! e-p matrix in Wanner representation
  ! rotation matrix U(k)
  real(kind=DP) :: xxq(3)
  ! kpoint for the interpolation (WARNING: this must be in crystal coord!)
  !
  !  output variables
  !
  complex(kind=DP) :: epmatf (nbnd, nbnd, nrr_k, nmodes)
  ! e-p matrix in Bloch representation, fine grid
  !
  ! work variables 
  !
  character (len=256) :: filint
  character (len=256) :: string 
  logical :: exst
  integer :: ibnd, jbnd, ir, ire, ir_start, ir_stop, imode,iunepmatwp2,ierr, i
  integer ::  ip , test
  integer(kind=8) ::  lrepmatw,  lrepmatw2
  real(kind=DP) :: rdotk
  complex(kind=DP) :: eptmp( nbnd, nbnd, nrr_k, nmodes)
  complex(kind=DP) :: cfac(nrr_q)
  complex(kind=DP):: aux( nbnd*nbnd*nrr_k*nmodes )
  !
  CALL start_clock('ephW2Bp')
  !----------------------------------------------------------
  !  STEP 3: inverse Fourier transform of g to fine k mesh
  !----------------------------------------------------------
  !
  !  g~ (k') = sum_R 1/ndegen(R) e^{-ik'R} g (R)
  !
  !  g~(k') is epmatf (nmodes, nmodes, ik )
  !  every pool works with its own subset of k points on the fine grid
  !
  IF (parallel_k) THEN
     CALL para_bounds(ir_start, ir_stop, nrr_q)
  ELSEIF (parallel_q) THEN
     ir_start = 1
     ir_stop  = nrr_q
  ELSE 
     CALL errore ('ephwan2blochp', 'Problem with parallel_k/q scheme', nrr_q)
  ENDIF
  !
#ifdef __PARA
  IF (.NOT. etf_mem) then
    filint = trim(prefix)//'.epmatwp1'
    CALL MPI_FILE_OPEN(intra_pool_comm,filint,MPI_MODE_RDONLY,MPI_INFO_NULL,iunepmatwp2,ierr)
    IF( ierr /= 0 ) CALL errore( 'ephwan2blochp', 'error in MPI_FILE_OPEN',1 )
  ENDIF
#endif
  ! CALL MPI_ERROR_STRING(ierr, string , i, ierr)
  ! inquire(FILE=filint,EXIST=exst)
  ! CALL MPI_FILE_GET_SIZE(iunepmatwp2,test,  ierr)
  !
  eptmp = czero
  cfac(:) = czero
  !
  DO ir = ir_start, ir_stop
     !   
     ! note xxq is assumed to be already in cryst coord
     !
     rdotk = twopi * dot_product ( xxq, dble(irvec(:, ir)) )
     cfac(ir) = exp( ci*rdotk ) / dble( ndegen(ir) )
  ENDDO
  ! 
  IF (etf_mem) then
    DO ir = ir_start, ir_stop
      eptmp(:,:,:,:) = eptmp(:,:,:,:) +&
        cfac(ir)*epmatwp( :, :, :, :, ir)
    ENDDO
  ELSE
    lrepmatw2   = 2 * nbnd * nbnd * nrr_k * nmodes
   ! IF( ierr /= 0 ) CALL errore( 'ephwan2blochp', 'error in mpi_set_view',ierr )
    DO ir = ir_start, ir_stop
#ifdef __PARA
      !  Direct read of epmatwp for this ir
      lrepmatw   = 2 * nbnd * nbnd * nrr_k * nmodes * 8 * (ir-1)
      ! SP: mpi view is used to set the position at which we should start
      ! reading the file. It is given in bits. 
      CALL MPI_FILE_SET_VIEW(iunepmatwp2,lrepmatw,MPI_DOUBLE_PRECISION,MPI_DOUBLE_PRECISION,'native',MPI_INFO_NULL,ierr)
      IF( ierr /= 0 ) CALL errore( 'ephwan2blochp', 'error in MPI_FILE_SET_VIEW',1 )
      CALL MPI_FILE_READ_ALL(iunepmatwp2, aux, lrepmatw2, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE,ierr)
      IF( ierr /= 0 ) CALL errore( 'ephwan2blochp', 'error in MPI_FILE_READ_ALL',1 )
      ! 
      i = 0
      DO imode = 1, nmodes
       DO ip = 1, nrr_k
        DO jbnd = 1, nbnd
         DO ibnd = 1, nbnd
           i = i + 1
           epmatw ( ibnd, jbnd, ip, imode ) = aux (i)
           ! 
         ENDDO
        ENDDO
       ENDDO
      ENDDO
#else      
      call rwepmatw ( epmatw, nbnd, nrr_k, nmodes, ir, iunepmatwp, -1)
#endif
      !
      !call rwepmatw ( epmatw, nbnd, nrr_k, nmodes, ir, iunepmatwp, -1)
      eptmp = eptmp + cfac(ir)*epmatw
    ENDDO
  ENDIF
  !
#ifdef __PARA
  IF (parallel_k) CALL mp_sum(eptmp, inter_pool_comm)
#endif
  !
  !----------------------------------------------------------
  !  STEP 4: un-rotate to Bloch space, fine grid
  !----------------------------------------------------------
  !
  ! epmatf(j) = sum_i eptmp(i) * uf(i,j)
  !
  Call zgemm( 'n', 'n', nbnd * nbnd * nrr_k, nmodes, nmodes, ( 1.d0, 0.d0 ),eptmp  , nbnd * nbnd * nrr_k, &
                                                                                    cuf, nmodes         , &
                                                             ( 0.d0, 0.d0 ),epmatf, nbnd * nbnd * nrr_k )
 ! DO ibnd = 1, nbnd
 !  DO jbnd = 1, nbnd
 !    !
 !    CALL zgemm ('n','n',nrr_k, nmodes, nmodes,(1.d0,0.d0),eptmp(ibnd,jbnd,:,:),nrr_k,&
 !        cuf,nmodes,(0.d0,0.d0), epmatf(ibnd,jbnd,:,:), nrr_k )
 !    !
 !  ENDDO
 ! ENDDO
  !
  CALL stop_clock('ephW2Bp')
  !
  end subroutine ephwan2blochp

