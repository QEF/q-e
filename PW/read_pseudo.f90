!
! Copyright (C) 2001-2007 Quantum-Espresso group
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
  USE pseudo_types,     ONLY : pseudo_upf, nullify_pseudo_upf, deallocate_pseudo_upf
  USE read_uspp_module, ONLY : readvan, readrrkj
  USE upf_to_internal,  ONLY : set_pseudo_upf
  USE atom,             ONLY :  msh, rgrid
  USE uspp_param, ONLY : newpseudo
  USE ions_base,  ONLY : ntyp => nsp
  USE funct,      ONLY : get_iexch, get_icorr, get_igcx, get_igcc
  USE io_files,   ONLY : pseudo_dir, psfile
  USE io_global,  ONLY : stdout
  USE ions_base,  ONLY : zv
  USE uspp_param, ONLY : upf
  use upf_module, ONLY : read_upf
  use radial_grids, ONLY : deallocate_radial_grid, nullify_radial_grid

  implicit none
  !
  real(DP), parameter :: rcut = 10.d0, eps = 1.0D-08
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

  IF( ALLOCATED( rgrid ) ) THEN
     DO nt = 1, SIZE( rgrid )
        CALL deallocate_radial_grid( rgrid( nt ) )
        CALL nullify_radial_grid( rgrid( nt ) )
     END DO
     DEALLOCATE( rgrid )
     DEALLOCATE( msh )
  END IF

  ALLOCATE( rgrid( ntyp ), msh( ntyp ) )

  DO nt = 1, ntyp
     CALL nullify_radial_grid( rgrid( nt ) )
  END DO

  IF( ALLOCATED( upf ) ) THEN
     DO nt = 1, SIZE( upf )
        CALL deallocate_pseudo_upf( upf( nt ) )
        CALL nullify_pseudo_upf( upf( nt ) )
     END DO
     DEALLOCATE( upf )
  END IF
  !
  ALLOCATE ( upf( ntyp ) )
  !
  !  nullify upf objects as soon as they are instantiated
  !
  do nt = 1, ntyp 
     CALL nullify_pseudo_upf( upf( nt ) )
  end do

  do nt = 1, ntyp
     !
     ! variables not necessary for USPP, but necessary for PAW, 
     ! they will be read from file if it is a PAW dataset.
     !
     rgrid(nt)%xmin = 0.d0
     rgrid(nt)%dx = 0.d0
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
     upf(nt)%grid => rgrid(nt)
     call read_upf(upf(nt), rgrid(nt), isupf, unit=iunps)
     !
     if (isupf == 0) then
        call set_pseudo_upf (nt, upf(nt))
        ! 
        ! UPF is assumed to be multi-projector
        !
        newpseudo (nt) = .true.
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
           ! PPs produced by Andrea Dal Corso's atomic code are assumed to
           ! be multiprojector; NCPP produced by Vanderbilt's core are not
           !    
           newpseudo (nt) = ( pseudo_type (psfile (nt) ) == 2 )
           !
           IF ( newpseudo (nt) ) THEN
              call readrrkj (iunps, nt, upf(nt))
           ELSE
              CALL readvan (iunps, nt, upf(nt))
           ENDIF
           CALL set_pseudo_upf (nt, upf(nt), rgrid(nt))
           !
        else
           newpseudo (nt) = .false.
           ! 
           call read_ncpp (iunps, nt, upf(nt))
           !
           CALL set_pseudo_upf (nt, upf(nt), rgrid(nt)) 
           !
        endif
        !
     endif
        !
     close (iunps)
     ! ... Zv = valence charge of the (pseudo-)atom, read from PP files,
     ! ... is set equal to Zp = pseudo-charge of the pseudopotential
     !
     zv(nt) = upf(nt)%zp
     !
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
     ! force msh to be odd for simpson integration
     !
5    msh (nt) = 2 * ( (msh (nt) + 1) / 2) - 1
     !
     ! Check that there are no zero wavefunctions
     !
     allocate ( chi2r (rgrid(nt)%mesh) )
     do nb = 1, upf(nt)%nwfc
        chi2r(:) = upf(nt)%chi (1:rgrid(nt)%mesh, nb ) **2
        call simpson (rgrid(nt)%mesh, chi2r(1), rgrid(nt)%rab, norm)
        !
        if ( norm < eps ) then
           WRITE( stdout,'(5X,"WARNING: atomic wfc # ",i2, &
                & " for atom type",i2," has zero norm")') nb, nt
           !
           ! set occupancy to a small negative number so that this wfc
           ! is not going to be used for starting wavefunctions
           !
           upf(nt)%oc (nb) = -eps
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

