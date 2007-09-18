!
! Copyright (C) 2001-2007 PWSCF group
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
  USE read_upf_module , ONLY : read_pseudo_upf
  USE read_uspp_module, ONLY : readvan, readrrkj
  USE upf_to_internal,  ONLY : set_pseudo_upf
  USE paw,              ONLY : set_paw_upf
  USE atom,       ONLY : chi, nchi, oc, msh, numeric, rgrid
  USE uspp_param, ONLY : zp, iver, tvanp, newpseudo
  USE ions_base,  ONLY : ntyp => nsp
  USE funct,      ONLY : get_iexch, get_icorr, get_igcx, get_igcc
  USE io_files,   ONLY : pseudo_dir, psfile
  USE io_global,  ONLY : stdout
  USE ions_base,  ONLY : zv
  USE pseud,      ONLY : lmax, lloc
  USE uspp_param, ONLY : lll, nbeta
  USE parameters, ONLY : nchix !PAW
  USE grid_paw_variables, ONLY : tpawp
  USE read_paw_module,    ONLY : paw_io, allocate_pseudo_paw, deallocate_pseudo_paw
  USE paw_to_internal,    ONLY : set_pseudo_paw
  implicit none
  !
  real(DP), parameter :: rcut = 10.d0, eps = 1.0D-08
  !
  TYPE (pseudo_upf) :: upf
  TYPE(paw_t) :: pawset
  !
  character(len=256) :: file_pseudo
  ! file name complete with path
  real(DP), allocatable :: chi2r(:)
  real(DP):: norm
  integer :: iunps, isupf, l, nt, nb, ir, ios
  integer :: iexch_, icorr_, igcx_, igcc_
  integer, external :: pseudo_type
  !
  iunps = 4
  l = len_trim (pseudo_dir)
  do nt = 1, ntyp
     tpawp(nt) = .false.
     !
     ! obsolescent variables, not read from UPF format, no longer used
     !
     iver(:,nt) = 0
     rgrid(nt)%xmin = 0.d0
     rgrid(nt)%dx = 0.d0
     lmax(nt) = -1
     lloc(nt) = -1
     ! 
     ! for compatibility - numeric is read by read_ncpp
     !
     numeric (nt) = .true.
     !
     ! add slash at the end if needed
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
        call set_paw_upf (nt, upf)
        CALL deallocate_pseudo_upf( upf )
        ! for compatibility with old formats
        newpseudo (nt) = .true.
        lmax(nt) = max ( lmax(nt), MAXVAL( lll( 1:nbeta(nt), nt) ) )
        !
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
           lmax(nt) = max ( lmax(nt), MAXVAL( lll( 1:nbeta(nt), nt) ) )
           !
        else if (pseudo_type (psfile (nt) ) ==3) then
           !
           !    PSEUDO PAW in temporary format. Use with care
           !
           !tpaw(nt)=.true.
           numeric (nt) = .true.
           newpseudo (nt) = .true.
           tvanp (nt) = .true.
           open (unit = iunps, file = file_pseudo, status = 'old', &
                 form='formatted', iostat = ios)
           call paw_io (pawset, iunps, "INP") !,ndmx,nchix,lmaxx)
           close (iunps)
           call set_pseudo_paw (nt, pawset)
           call deallocate_pseudo_paw (pawset)
           !
        else
           tvanp (nt) = .false.
           newpseudo (nt) = .false.
           ! 
           call read_ncpp (nt, iunps)
           !
        endif
     endif
     close (iunps)
     !
     ! ... Zv = valence charge of the (pseudo-)atom, read from PP files,
     ! ... is set equal to Zp = pseudo-charge of the pseudopotential
     !
     zv(nt) = zp(nt)
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
     ! the radial grid is defined up to r(mesh) but we introduce 
     ! an auxiliary variable msh to limit the grid up to rcut=10 a.u. 
     ! This is used to cut off the numerical noise arising from the
     ! large-r tail in cases like the integration of V_loc-Z/r
     !
     do ir = 1, rgrid(nt)%mesh
        if (rgrid(nt)%r(ir) > rcut) then
           msh (nt) = ir
           goto 5
        endif
     enddo
     msh (nt) = rgrid(nt)%mesh 
     !
     ! force msh to be odd for simpson integration (maybe obsolete)
     !
5    msh (nt) = 2 * ( (msh (nt) + 1) / 2) - 1
     !
     ! Check that there are no zero wavefunctions
     !
     allocate ( chi2r (rgrid(nt)%mesh) )
     do nb = 1, nchi (nt)
        chi2r(:) = chi ( :rgrid(nt)%mesh, nb, nt ) **2
        call simpson (rgrid(nt)%mesh, chi2r(1), rgrid(nt)%rab, norm)
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
  if (l > 3) then
     if (psfile (l - 3:l) .eq.'.PAW') pseudo_type = 3
  end if
  !
  return

end function pseudo_type

