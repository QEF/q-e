!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine openfilq
  !-----------------------------------------------------------------------
  !
  !     This subroutine opens all the files necessary for the phononq
  !     calculation.
  !

  use pwcom
  use mp, only: mp_end
  use io_files, only: prefix, iunigk
  USE kinds, only : DP
  use phcom
  use control_flags, ONLY : twfcollect
#ifdef __PARA
  use para
#endif
  use restart_module, only: readfile_new
  implicit none
  integer :: ios
  ! integer variable for I/O control
  character (len=256) :: filint
  ! the name of the file
  logical :: exst
  ! logical variable to check file existe
  !
  real(kind=DP) :: edum(1,1), wdum(1,1)
  integer :: ndr, ierr, kunittmp

  if (len_trim(prefix) == 0) call errore ('openfilq', 'wrong prefix', 1)

  !
  !     There are six direct access files to be opened in the tmp area
  !
  !     The file with the wavefunctions
  !
  iuwfc = 20

  lrwfc = 2 * nbnd * npwx
  filint = trim(prefix)//'.wfc'
  call diropn (iuwfc, filint, lrwfc, exst)

  if (.not.exst) then
    ndr      = 4
    kunittmp = 1
#  ifdef __PARA
    kunittmp = kunit
#  endif
    call readfile_new( 'wave', ndr, edum, wdum, kunittmp, lrwfc, iuwfc, ierr )
    if( ierr > 0 ) &
      call errore ('openfilq', 'file '//filint//' not found', 1)
    twfcollect=.not.exst
  end if
  !
  !    The file with deltaV_{bare} * psi
  !
  iubar = 21
  lrbar = 2 * nbnd * npwx
  filint = trim(prefix) //'.bar'
  call diropn (iubar, filint, lrbar, exst)
  if (recover.and..not.exst) call errore ('openfilq','file bar not found', 1)
  !
  !    The file with the solution delta psi
  !
  iudwf = 22
  lrdwf = 2 * nbnd * npwx
  filint = trim(prefix) //'.dwf'
  call diropn (iudwf, filint, lrdwf, exst)
  if (recover.and..not.exst) call errore ('openfilq','file dwf not found', 1)
  !
  !   open a file with the static change of the charge
  !
  if (okvan) then
     iudrhous = 25
     lrdrhous = 2 * nrxx * nspin
     filint = trim(prefix) //'.prd'
     call diropn (iudrhous, filint, lrdrhous, exst)
     if (recover.and..not.exst) call errore ('openfilq','file prod not found', 1)
  endif
  !
  !   An optional file for testing purposes containing the deltarho
  !
  if (fildrho.ne.' ') then
     iudrho = 23
     lrdrho = 2 * nrx1 * nrx2 * nrx3 * nspin
#ifdef __PARA
     if (me.ne.1) goto 300
#endif
     if(epsil) then
        filint = trim(fildrho)//".E"
     else
        filint = trim(fildrho)//".u"
     end if
     call diropn (iudrho, filint, lrdrho, exst)
#ifdef __PARA
300  continue
#endif
  endif
  !
  !   Here the sequential files
  !
  !   The igk at a given k (and k+q if q!=0)
  !
  iunigk = 24
  filint = trim(prefix) //'.igk'
  call seqopn (iunigk, filint, 'unformatted', exst)
  !
  !   a formatted file which contains the dynamical matrix in cartesian
  !   coordinates is opened in the current directory
#ifdef __PARA
  !   ... by the first node only, other nodes write on unit 6 (i.e./dev/nu
  !   exception: electron-phonon calculation from saved data
  !  (iudyn is read, not written, by all nodes)
  !
  if ( (me.ne.1.or.mypool.ne.1) .and. (.not.elph.or.trans) ) then
     iudyn = 6
     goto 400
  endif
#endif
  if (trans.or.elph) then
     iudyn = 26
     open (unit = iudyn, file = fildyn, status = 'unknown', err = &
          100, iostat = ios)
100  call errore ('openfilq', 'opening file'//fildyn, abs (ios) )
     rewind (iudyn)
  endif
  !
  !   An optional file for electron-phonon calculations containing deltaVs
  !
400 if (fildvscf.ne.' ') then
     iudvscf = 27
     call seqopn (iudvscf, fildvscf, 'unformatted', exst)
     rewind (iudvscf)
  end if
  !
  !    In the USPP case we need two files for the Commutator, the first is
  !    given by filbar and a second which just contains P_c x |psi>,
  !    which is required for the calculation of the Born effective carges
  !
  if (okvan .and. (epsil .or. zue)) then
     iucom = 28
     lrcom = 2 * nbnd * npwx
     filint = trim(prefix) //'.com'
     call diropn (iucom, filint, lrcom, exst)
     if (recover.and..not.exst) &
         call errore ('openfilq', 'file com not found', 1)
  !
  !    In the USPP case we also need a file in  order to store derivatives 
  !    ok kb projectors
  !  
     iudvkb3 = 29
     lrdvkb3 = 2 * npwx * nkb * 3
     filint = trim(prefix) //'.dvkb3'
     call diropn (iudvkb3, filint, lrdvkb3, exst)
     if (recover.and..not.exst) &
         call errore ('openfilq', 'file dvkb3 not found', 1)
  endif
  if (epsil .or. zue) then
     iuebar = 30
     lrebar = 2 * nbnd * npwx
     filint = trim(prefix) //'.ebar'
     call diropn (iuebar, filint, lrebar, exst)
     if (recover.and..not.exst) &
        call errore ('openfilq','file ebar not found', 1)
  endif
  
  return
end subroutine openfilq
