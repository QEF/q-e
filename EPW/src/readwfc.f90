  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !--------------------------------------------------------
  subroutine readwfc ( ipool, recn, evc0 )
  !--------------------------------------------------------
  !
  !  open wfc files as direct access, read, and close again
  !
  ! RM - Nov/Dec 2014
  ! Imported the noncolinear case implemented by xlzhang
  !
  !-------------------------------------------------------------
  USE io_files, ONLY : prefix, tmp_dir
  USE units_ph, ONLY : lrwfc, iuwfc
  USE kinds,    ONLY : DP
  USE wvfct,    ONLY : npwx
  USE pwcom,    ONLY : nbnd
  USE noncollin_module,ONLY : npol
  USE mp_global,ONLY : nproc_pool, me_pool
  USE mp_global,ONLY : npool
  !
  implicit none
  integer :: recn, ipool
  !  kpoint number
  !  poolfile to be read (not used in serial case)
  complex(kind=DP) :: evc0 ( npwx*npol, nbnd )
  character (len=3) :: nd_nmbr0
  ! node number for shuffle
  !
  integer :: unf_recl, ios
  character (len=256) :: tempfile

  !
  !  open the wfc file, read and close
  !
  CALL set_ndnmbr ( ipool, me_pool, nproc_pool, npool, nd_nmbr0)
  !
#if defined (__ALPHA)
#  define DIRECT_IO_FACTOR 2    
# else
#  define DIRECT_IO_FACTOR 8
#endif

#if defined(__MPI)
  tempfile = trim(tmp_dir) // trim(prefix) // '.wfc' // nd_nmbr0
# else
  tempfile = trim(tmp_dir) // trim(prefix) // '.wfc'
#endif
  unf_recl = DIRECT_IO_FACTOR * lrwfc
  !
  open  (iuwfc, file = tempfile, form = 'unformatted', &
          access = 'direct', iostat=ios, recl = unf_recl)
  IF (ios /= 0) call errore ('readwfc', 'error opening wfc file', iuwfc)
  read  (iuwfc, rec = recn) evc0
  close (iuwfc, status = 'keep')
  !
  !
  end subroutine readwfc

