!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

   subroutine write_vpot_matrix( vmat, ort)
!this subroutine writes the coulomb potential on the basis
! of orthonormalized products of wanniers
!to be read by GWW code


  USE kinds, ONLY : DP
  USE wannier_gw, ONLY : numw_prod
  USE io_global, ONLY : stdout
  USE io_files, ONLY : prefix, tmp_dir


  implicit none

  INTEGER, EXTERNAL :: find_free_unit

  REAL(kind=DP) :: vmat(numw_prod,numw_prod)
  INTEGER :: ort!if ort==0 writes nonorthogonal file, if ort==1 writes orthogonal file, 
                !if ort==2 writes nonorthogonal ^1/2 file
  INTEGER :: iunu, iw

  iunu = find_free_unit()

  if(ort == 1) then
     open(unit=iunu,file=trim(tmp_dir)//trim(prefix)//'.vpot',status='unknown',form='unformatted')
     !     open(unit=iunu,file=trim(tmp_dir)//trim(prefix)//'.vpot',status='unknown',form='formatted')
  else if(ort==0) then 
     open(unit=iunu,file=trim(tmp_dir)//trim(prefix)//'.vpot_no',status='unknown',form='unformatted')
!     open(unit=iunu,file=trim(tmp_dir)//trim(prefix)//'.vpot_no',status='unknown',form='formatted')
  else if(ort==2) then
     open(unit=iunu,file=trim(tmp_dir)//trim(prefix)//'.vpot_no_sym',status='unknown',form='unformatted')
  else if(ort==3) then
     open(unit=iunu,file=trim(tmp_dir)//trim(prefix)//'.vpot_no_zero',status='unknown',form='unformatted')
  else if(ort==4) then
     open(unit=iunu,file=trim(tmp_dir)//trim(prefix)//'.vpot_no_sym_zero',status='unknown',form='unformatted')
  endif

  write(iunu) numw_prod
!  write(iunu,*) numw_prod

  do iw=1,numw_prod
    write(iunu) vmat(1:numw_prod,iw)
!     write(iunu,*) vmat(1:numw_prod,iw)
  enddo

  close(iunu)

  return
  end subroutine

