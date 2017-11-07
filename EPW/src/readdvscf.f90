  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !--------------------------------------------------------
  subroutine readdvscf ( dvscf, recn, iq, nqc )
  !--------------------------------------------------------
  !!
  !!  open dvscf files as direct access, read, ad close again
  !!
  !! SP - Nov 2017
  !! Replaced fstat by Fortran instric function inquire. 
  !! 
  !! RM - Nov/Dec 2014
  !! Imported the noncolinear case implemented by xlzhang
  !!
  !-------------------------------------------------------------
#if defined(__ALPHA)
#  define DIRECT_IO_FACTOR 2
#else
#  define DIRECT_IO_FACTOR 8
#endif
  !
  USE io_files,  ONLY : prefix
  USE units_ph,  ONLY : lrdrho
  USE kinds,     ONLY : DP
  USE fft_base,  ONLY : dfftp
  USE pwcom
  USE epwcom,    ONLY : dvscf_dir
  USE noncollin_module, ONLY : nspin_mag
  USE io_epw,    ONLY : iudvscf
  !
  implicit none
  ! 
  INTEGER, INTENT (in) :: recn
  !! perturbation number
  INTEGER, INTENT (in) :: iq
  !! the current q point
  INTEGER, INTENT (in) :: nqc
  !! the total number of qpoints in the list
  !
  COMPLEX(kind=DP), INTENT (out) :: dvscf ( dfftp%nnr , nspin_mag) 
  !! dVscf potential is read from file
  !
  ! Local variables
  INTEGER :: unf_recl,ios
  CHARACTER (len=256) :: tempfile
  CHARACTER (len=3) :: filelab
  INTEGER(kind=8) :: mult_unit, file_size
  !
  !  the call to set_ndnmbr is just a trick to get quickly
  !  a file label by exploiting an existing subroutine
  !  (if you look at the sub you will find that the original 
  !  purpose was for pools and nodes)

  ! DBSP:
  !  Iotemp is a output variable and it does not matter whether it is integer or any other type. 
  !   
  CALL set_ndnmbr ( 0, iq, 1, nqc, filelab)
  tempfile = trim(dvscf_dir) // trim(prefix) // '.dvscf_q' // filelab
  unf_recl = DIRECT_IO_FACTOR * lrdrho
  !unf_recl = iofactor * lrdrho
  !DBSP
  !print*,'iofactor ',iofactor
  mult_unit = unf_recl
  mult_unit = recn * mult_unit
  !
  !
  !  open the dvscf file, read and close
  !
  open  (iudvscf, file = tempfile, form = 'unformatted', &
          access = 'direct', iostat=ios,recl = unf_recl,status='old')
  IF (ios /= 0) call errore ('readdvscf','error opening ' // tempfile, iudvscf)
  !
  ! check that the binary file is long enough
  INQUIRE(FILE=tempfile, SIZE=file_size)
  IF (mult_unit .gt. file_size) call errore('readdvscf', &
       trim(tempfile)//' too short, check ecut', iudvscf)
  !
  read  (iudvscf, rec = recn) dvscf
  close (iudvscf, status = 'keep')
  !
  !
  end subroutine readdvscf
