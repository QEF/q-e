!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!--------------------------------------------------------------------
subroutine openfil
  !--------------------------------------------------------------------
  !
  !    This routine opens all files needed to the self consistent run,
  !    sets various file names, units, record lengths
  !
  USE io_global,  ONLY : stdout
  use pwcom
  use io_files, only: prefix
  use restart_module, only: readfile_new
#ifdef __PARA
  use para
#endif
  implicit none
  logical :: exst
  integer :: ndr, kunittmp, ierr
  real(kind=DP) :: edum(1,1), wdum(1,1)
  !
  !     iunwfc contains the wavefunctions
  !
  iunwfc = 10
  !
  !     iunoldwfc contains the old wavefunctions, used in molecular dynamics
  !
  iunoldwfc = 11
  iunoldwfc2= 12
  !
  !     nwordwfc is the record length for the direct-access file
  !     containing wavefunctions
  !
  nwordwfc = 2 * nbnd * npwx
  !
  call diropn (iunwfc, trim(prefix)//'.wfc', nwordwfc, exst)

  if (startingwfc.eq.'file'.and..not.exst) then

#if defined __NEW_PUNCH

     ndr      = 4
     kunittmp = 1
#  ifdef __PARA
     kunittmp = kunit
#  endif

     call readfile_new( 'wave', ndr, edum, wdum, kunittmp, nwordwfc, iunwfc, ierr )
     if ( ierr > 0 ) then

#else

       WRITE( stdout, '(5x,"Cannot read wfc file: not found")')
       startingwfc='atomic'

#endif

#if defined __NEW_PUNCH

     end if

#endif

  endif
  !
  !
  ! Needed for LDA+U
  !
  !     iunat contains the orthogonalized wfcs
  !
  iunat = 13
  nwordatwfc = 2 * npwx * natomwfc
  if (lda_plus_u) &
       call diropn (iunat, trim(prefix)//'.atwfc', nwordatwfc, exst)
  !
  !     iunocc contains the atomic occupations computed in new_ns
  !     it is opened and closed for each reading-writing operation
  !
  iunocc = 14
  !
  !     if extrapolation of wfc's is requested (order=2)
  !     another file is needed to store the "old" wfc's
  !
  if (order > 1) &
       call diropn (iunoldwfc, trim(prefix)//'.oldwfc',nwordwfc, exst)
  !
  if (order > 2) &
       call diropn (iunoldwfc2,trim(prefix)//'.oldwfc2', nwordwfc, exst)
  !
  !    iunigk contains the number of PW and the indices igk
  !    Note that unit 15 is reserved for error messages 
  !
  iunigk = 16
  call seqopn (iunigk, trim(prefix)//'.igk', 'unformatted', exst)
  !
  return

end subroutine openfil

