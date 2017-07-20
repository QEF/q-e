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
  subroutine ephwan2blochp_mem (imode, nmodes, xxq, irvec, ndegen, nrr_q, cuf, epmatf, nbnd, nrr_k )
  !---------------------------------------------------------------------------
  !!
  !! Even though this is for phonons, I use the same notations
  !! adopted for the electronic case (nmodes->nmodes etc)
  !!
  USE kinds,         only : DP
  USE epwcom,        only : parallel_k, parallel_q, etf_mem
  USE elph2,         only : epmatwp
  USE constants_epw, ONLY : twopi, ci, czero
  USE io_files,      ONLY : prefix, tmp_dir
  USE io_epw,        ONLY : iunepmatwp
  USE mp_global,     ONLY : mp_sum
  USE mp_world,      ONLY : world_comm
  USE parallel_include
  implicit none
  !
  !  input variables
  !
  INTEGER, INTENT (in) :: imode
  !! Current mode  
  INTEGER, INTENT (in) :: nmodes
  !! Total number of modes
  INTEGER, INTENT (in) :: nrr_q
  !! Number of WS points
  INTEGER, INTENT (in) :: irvec ( 3, nrr_q)
  !! Coordinates of WS points
  INTEGER, INTENT (in) :: ndegen (nrr_q)
  !! Number of degeneracy of WS points
  INTEGER, INTENT (in) :: nbnd
  !! Number of bands
  INTEGER, INTENT (in) ::  nrr_k
  !! Number of electronic WS points
  REAL(kind=DP) :: xxq(3)
  !! Kpoint for the interpolation (WARNING: this must be in crystal coord!)
  COMPLEX(kind=DP), INTENT (in) :: cuf (nmodes, nmodes)
  !! e-p matrix in Wanner representation
  COMPLEX(kind=DP), INTENT (out) :: epmatf (nbnd, nbnd, nrr_k)
  !! e-p matrix in Bloch representation, fine grid
  ! 
  ! Local variables 
  !
  CHARACTER (len=256) :: filint
  !! File name
  !
  INTEGER :: ir
  !! Real space WS index
  INTEGER :: ir_start
  !! Starting ir for this cores
  INTEGER :: ir_stop
  !! Ending ir for this pool
  INTEGER :: iunepmatwp2
  !! Return the file unit
  INTEGER :: ierr
  !! Return if there is an error
  INTEGER (kind=MPI_OFFSET_KIND) :: lrepmatw
  !! Offset to tell where to start reading the file
  INTEGER (kind=MPI_OFFSET_KIND) :: lrepmatw2
  !! Offset to tell where to start reading the file
  !
  REAL(kind=DP) :: rdotk
  !! Exponential for the FT
  !
  COMPLEX(kind=DP) :: cfac(nrr_q)
  !! Factor for the FT
  COMPLEX(kind=DP), ALLOCATABLE :: epmatw ( :,:,:)
  !! El-ph matrix elements
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
  CALL para_bounds(ir_start, ir_stop, nrr_q)
  !
  filint = trim(tmp_dir)//trim(prefix)//'.epmatwp1'
  CALL MPI_FILE_OPEN(world_comm,filint,MPI_MODE_RDONLY,MPI_INFO_NULL,iunepmatwp2,ierr)
  IF( ierr /= 0 ) CALL errore( 'ephwan2blochp_mem', 'error in MPI_FILE_OPEN',1 )
  !
  cfac(:) = czero
  !
  DO ir = ir_start, ir_stop
     !   
     ! note xxq is assumed to be already in cryst coord
     rdotk = twopi * dot_product ( xxq, dble(irvec(:, ir)) )
     cfac(ir) = exp( ci*rdotk ) / dble( ndegen(ir) )
  ENDDO
  ! 
  ALLOCATE(epmatw ( nbnd, nbnd, nrr_k))
  !
  lrepmatw2 = 2_MPI_OFFSET_KIND * INT( nbnd  , kind = MPI_OFFSET_KIND ) * &
                                  INT( nbnd  , kind = MPI_OFFSET_KIND ) * &
                                  INT( nrr_k , kind = MPI_OFFSET_KIND )
  ! 
  DO ir = ir_start, ir_stop
    !
    ! SP: The following needs a small explaination: although lrepmatw is correctly defined as kind 8 bits or 
    !     kind=MPI_OFFSET_KIND, the number "2" and "8" are default kind 4. The other as well. Therefore
    !     if the product is too large, this will crash. The solution (kind help recieved from Ian Bush) is below:
    lrepmatw = 2_MPI_OFFSET_KIND * INT( nbnd  , kind = MPI_OFFSET_KIND ) * &
                                   INT( nbnd  , kind = MPI_OFFSET_KIND ) * &
                                   INT( nrr_k , kind = MPI_OFFSET_KIND ) * &
                                   INT( nmodes, kind = MPI_OFFSET_KIND ) * &
             8_MPI_OFFSET_KIND * ( INT( ir    , kind = MPI_OFFSET_KIND ) - 1_MPI_OFFSET_KIND ) + &
             2_MPI_OFFSET_KIND *   INT( nbnd  , kind = MPI_OFFSET_KIND ) * &
                                   INT( nbnd  , kind = MPI_OFFSET_KIND ) * &
                                   INT( nrr_k , kind = MPI_OFFSET_KIND ) * &
             8_MPI_OFFSET_KIND * ( INT( imode , kind = MPI_OFFSET_KIND ) - 1_MPI_OFFSET_KIND )
    !
    ! SP: mpi seek is used to set the position at which we should start
    ! reading the file. It is given in bits. 
    ! Note : The process can be collective (=blocking) if using MPI_FILE_SET_VIEW & MPI_FILE_READ_ALL
    !        or noncollective (=non blocking) if using MPI_FILE_SEEK & MPI_FILE_READ. 
    !        Here we want non blocking because not all the process have the same nb of ir. 
    !
    CALL MPI_FILE_SEEK(iunepmatwp2,lrepmatw,MPI_SEEK_SET,ierr)
    IF( ierr /= 0 ) CALL errore( 'ephwan2blochp', 'error in MPI_FILE_SEEK',1 )
    CALL MPI_FILE_READ(iunepmatwp2, epmatw, lrepmatw2, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE,ierr)
    IF( ierr /= 0 ) CALL errore( 'ephwan2blochp', 'error in MPI_FILE_READ_ALL',1 )
    ! 
    !
    CALL ZAXPY(nbnd * nbnd * nrr_k, cfac(ir), epmatw, 1, epmatf, 1)
    ! 
  ENDDO
  DEALLOCATE(epmatw)
  !
  CALL mp_sum(epmatf, world_comm)
  !
  CALL MPI_FILE_CLOSE(iunepmatwp2,ierr)
  IF( ierr /= 0 ) CALL errore( 'ephwan2blochp_mem', 'error in MPI_FILE_CLOSE',1 )
  !
  CALL stop_clock('ephW2Bp')
  !
  end subroutine ephwan2blochp_mem
