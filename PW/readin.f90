!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!-----------------------------------------------------------------------
subroutine readpp
  !-----------------------------------------------------------------------
  !
  !    Read pseudopotentials
  !
  USE kinds,      ONLY : DP
  USE pseudo_types
  USE read_pseudo_module
  USE upf_to_internal
  USE atom,       ONLY : chi, nchi, oc, mesh, rab, numeric, xmin, dx
  USE uspp_param, ONLY : iver, tvanp, newpseudo
  USE ions_base,  ONLY : ntyp => nsp
  USE funct,      ONLY : get_iexch, get_icorr, get_igcx, get_igcc
  USE io_files,   ONLY : pseudo_dir, psfile
  USE io_global,  ONLY : stdout
  use ions_base,  only: zv
  use pseud,      only: zp, lmax, lloc
  USE uspp_param, ONLY : lll, nbeta
  implicit none
  !
  TYPE (pseudo_upf) :: upf
  !
  character(len=256) :: file_pseudo
  ! file name complete with path
  real(DP), allocatable :: chi2r(:)
  real(DP):: norm, eps = 1.0D-08
  integer :: iunps, isupf, l, nt, nb, ios
  integer :: iexch_, icorr_, igcx_, igcc_
  integer, external :: pseudo_type
  !
  iunps = 4
  l = len_trim (pseudo_dir)
  do nt = 1, ntyp
     !
     ! iver, xmin, dx are not read from UPF format
     !
     iver(:,nt) = 0
     xmin(nt) = 0.d0
     dx(nt) = 0.d0
     !
     ! obsolescent variables
     !
     lmax(nt) = -1
     lloc(nt) = -1
     numeric (nt) = .true.
     !
     ! add / if needed
     !
     if (pseudo_dir (l:l) .ne.'/') then
        file_pseudo = pseudo_dir (1:l) //'/'//psfile (nt)
     else
        file_pseudo = pseudo_dir (1:l) //psfile (nt)
     endif
     !
     open (unit = iunps, file = file_pseudo, status = 'old', form = &
          'formatted', iostat = ios)
     call errore ('readpp', 'file '//trim(file_pseudo)//' not found', ios)
     !
     ! read UPF  pseudopotentials - the UPF format is detected via the
     ! presence of the keyword '<PP_HEADER>' at the beginning of the file
     !
     call read_pseudo_upf(iunps, upf, isupf)
     !
     if (isupf == 0) then
        call set_pseudo_upf (nt, upf)
        CALL deallocate_pseudo_upf( upf )
        ! UPF is RRKJ3-like
        newpseudo (nt) = .true.
     else
        rewind (unit = iunps)
        !
        !     The type of the pseudopotential is determined by the file name:
        !    *.vdb or *.van  Vanderbilt US pseudopotential code  pseudo_type=1
        !    *.RRKJ3         Andrea's   US new code              pseudo_type=2
        !    none of the above: PWSCF norm-conserving format     pseudo_type=0
        !
        if ( pseudo_type (psfile (nt) ) == 1 .or. &
             pseudo_type (psfile (nt) ) == 2 ) then
           !
           !    newpseudo distinguishes beteween US pseudopotentials
           !    produced by Vanderbilt code and those produced
           !    by Andrea's atomic code.
           !
           newpseudo (nt) = ( pseudo_type (psfile (nt) ) == 2 )
           !
           if ( newpseudo (nt) ) then
              call readrrkj (nt, iunps)
           else
              call readvan (nt, iunps)
           endif
           !
           zp (nt) = zv (nt)
           lmax(nt) = max ( lmax(nt), MAXVAL( lll( 1:nbeta(nt), nt) ) )
           !
        else
           tvanp (nt) = .false.
           newpseudo (nt) = .false.
           ! numeric is read inside read_ncpp
           call read_ncpp (nt, iunps)
           !
        endif
     endif
     close (iunps)
     !
     if (nt == 1) then
        iexch_ = get_iexch()
        icorr_ = get_icorr()
        igcx_  = get_igcx()
        igcc_  = get_igcc()
     else
        if ( iexch_ /= get_iexch() .or. icorr_ /= get_icorr() .or. &
             igcx_  /= get_igcx()  .or. igcc_ /= get_igcc() ) then
           CALL errore( 'readpp','inconsistent DFT read',nt)
        end if
     end if
     !
     ! Check that there are no zero wavefunctions
     !
     allocate ( chi2r (mesh(nt)) )
     do nb = 1, nchi (nt)
        chi2r(:) = chi ( :mesh(nt), nb, nt ) **2
        call simpson (mesh(nt), chi2r(1), rab(1,nt), norm)
        !
        if ( norm < eps ) then
           WRITE( stdout,'(5X,"WARNING: atomic wfc # ",i2, &
                & " for atom type",i2," has zero norm")') nb, nt
           !
           ! set occupancy to a small negative number so that this wfc
           ! is not going to be used for starting wavefunctions
           !
           oc (nb, nt) = -eps
        end if
     enddo
     deallocate ( chi2r )
     !
  enddo
  !
  return
end subroutine readpp
!-----------------------------------------------------------------------
integer function pseudo_type (psfile)
  !-----------------------------------------------------------------------
  implicit none
  character (len=*) :: psfile
  integer :: l
  !
  l = len_trim (psfile)
  pseudo_type = 0
  if (psfile (l - 3:l) .eq.'.vdb'.or.psfile (l - 3:l) .eq.'.van') &
       pseudo_type = 1
  if (l > 5) then
     if (psfile (l - 5:l) .eq.'.RRKJ3') pseudo_type = 2
  end if
  !
  return

end function pseudo_type

