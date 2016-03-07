  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-------------------------------------------------------------
  subroutine readigk ( ipool, recn, npw0, igk0 )
  !-------------------------------------------------------------
  !
  !  open igk files as direct access, read, ad close again
  !
  !
  !-------------------------------------------------------------
  use io_files, only : prefix, iunigk, tmp_dir
  use wvfct,    only : npwx
#ifdef __PARA
  use mp_global,only : nproc_pool, me_pool
  use mp_global,only : npool,my_pool_id
  USE mp_world, ONLY : mpime
#endif
  !
  implicit none
  integer :: npw0, igk0 (npwx), recn, ipool
  !  number of planewaves
  !  igk's
  !  kpoint number
  !  poolfile to be read (not used in serial case)
#ifdef __PARA
  character (len=3) :: nd_nmbr0
#endif
  ! node number for shuffle
  !
  integer :: itmp, lrigk, unf_recl, itmp1
  character (len=256) :: tempfile
  !
  ! the following is for pgi on opteron, I should test the rigth numbers
  ! on other machines/compilers
  !
#if defined(__PGI)||defined(__AIX)
#  define INT_DIRECT_IO_FACTOR 4
#else
!this must be here to compile on a serial machine
!however the value is not tested
#  define INT_DIRECT_IO_FACTOR 4
!jn won't compile with this on civet:  call errore('INT_DIRECT_IO_FACTOR is not a tested quantity',-1)
#endif
  !
  ! filename
  !
  tempfile = trim(tmp_dir) // trim(prefix) // '.igk'
  !
#ifdef __PARA
  CALL set_ndnmbr (ipool, me_pool, nproc_pool, npool, nd_nmbr0)
  IF (ipool.ne.1.or.me_pool.ne.0) tempfile = trim(tempfile) // nd_nmbr0
#endif
  !
  ! record lenght: 1+1+npwx+1
  ! first and last byte is rec len, second is npw, 3dh to npwx-rh is igk
  !
  lrigk = npwx + 3
  unf_recl = INT_DIRECT_IO_FACTOR * lrigk
  !
  ! open, read, and close
  !
  open  (iunigk, file = tempfile, form = 'unformatted', &
         access = 'direct', recl = unf_recl)
  read  (iunigk, rec = recn) itmp, npw0, igk0, itmp1
  close (iunigk, status = 'keep')
  !
  end subroutine readigk

