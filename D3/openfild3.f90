!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------

subroutine openfild3
  !-----------------------------------------------------------------------
  !
  !     This subroutine opens all the files necessary for the
  !     third derivative calculation.
  !
  use pwcom
  use phcom
  use d3com
#ifdef __PARA
  use para
#endif

  implicit none
  integer :: ios
  ! integer variable for I/O control
  character (len=256) :: filint
  ! the name of the file
  logical :: exst
  ! logical variable to check file existe


  if (len_trim(filpun).eq.0) call errore ('openfild3', 'wrong filpun name', 1)
  !
  !     The file with the wavefunctions
  !
  iuwfc = 20

  lrwfc = 2 * nbnd * npwx
  filint = trim(filpun) //'.wfc'
  call diropn (iuwfc, filint, lrwfc, exst)
  if (.not.exst) call errore ('openfild3', 'file '//filint//' not found', 1)
  !
  !    The file with deltaV_{bare} * psi
  !
  iubar = 21
  lrbar = 2 * nbnd * npwx
  filint = trim(filpun) //'.bar'
  call diropn (iubar, filint, lrbar, exst)
  if (recover.and..not.exst) call errore ('openfild3', 'file bar not &
       &found', 1)
  !
  !    The file with the solution delta psi
  !
  iudwf = 22
  lrdwf = 2 * nbnd * npwx
  filint = trim(filpun) //'.dwf'
  call diropn (iudwf, filint, lrdwf, exst)
  if (recover.and..not.exst) call errore ('openfild3', 'file dwf not &
       &found', 1)
  !
  !   Here the sequential files
  !
  !   The igk at a given k (and k+q if q!=0)
  !
  iunigk = 24
  filint = trim(filpun) //'.igk'
  call seqopn (iunigk, filint, 'unformatted', exst)
  !
  !   a formatted file which contains the dynamical matrix in cartesian
  !   coordinates is opened in the current directory
#ifdef __PARA
  ! ... by the first node only, other nodes write on unit 6 (i.e. /dev/null)
  !
  if (me.ne.1.or.mypool.ne.1) then
     iudyn = 6
     goto 100
     return
  endif
#endif
  iudyn = 26
  open (unit = iudyn, file = fildyn, status = 'unknown', err = 110, &
       iostat = ios)
110 call errore ('openfild3', 'opening file'//fildyn, abs (ios) )
  rewind (iudyn)
#ifdef __PARA



100 continue
#endif
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  ! Variation of the charge density with respect to a perturbation
  ! with a generic q
  !
  iudrho = 25
  iud0rho = 33
  if (lgamma) iud0rho = iudrho
  lrdrho = 2 * nrx1 * nrx2 * nrx3 * nspin
#ifdef __PARA
  !
  !   is opened only by the first task of each pool
  !
  if (me.ne.1) goto 120
#endif
  filint = trim(fildrho)
  call diropn (iudrho, filint, lrdrho, exst)
  !
  ! Variation of the charge density with respect to a perturbation with q=
  ! Not needed if q=0
  !
  if (.not.lgamma) then
     filint = trim(fild0rho)
     call diropn (iud0rho, filint, lrdrho, exst)
  endif
#ifdef __PARA
120 continue
#endif
  !
  ! If q=0, we need only one file with the variation of the wavefunctions
  !
  iud0qwf = iudwf
  iudqwf = iudwf
  if (.not.lgamma) then
     !
     !    Open the file with the solution q=0 delta psi
     !
     iud0qwf = 34
     filint = trim(filpun) //'.d0wf'
     call diropn (iud0qwf, filint, lrdwf, exst)
     !
     !    Open the file with the solution q=0 delta psi
     !
     iudqwf = 35
     filint = trim(filpun) //'.dqwf'
     call diropn (iudqwf, filint, lrdwf, exst)
  endif
  !
  !    The file with   <psi| dqV |psi>
  !
  iupdqvp = 36
  lrpdqvp = 2 * nbnd * nbnd
  filint = trim(filpun) //'.pdp'
  call diropn (iupdqvp, filint, lrpdqvp, exst)
  !
  !    The file with   <psi| d0V |psi>
  !
  iupd0vp = iupdqvp
  if (.not.lgamma) then
     iupd0vp = 37
     filint = trim(filpun) //'.p0p'
     call diropn (iupd0vp, filint, lrpdqvp, exst)

  endif
  if (degauss.ne.0.d0) then
     !
     !    The file with   <dqpsi| dqV |psi> (only in the metallic case)
     !
     iudpdvp_1 = 38
     lrdpdvp = 2 * nbnd * nbnd
     filint = trim(filpun) //'.pv1'
     call diropn (iudpdvp_1, filint, lrdpdvp, exst)
     !
     !    The file with   <dqpsi| d0V |psi>
     !
     iudpdvp_2 = iudpdvp_1
     iudpdvp_3 = iudpdvp_1
     if (.not.lgamma) then
        iudpdvp_2 = 39
        filint = trim(filpun) //'.pv2'
        call diropn (iudpdvp_2, filint, lrdpdvp, exst)
        !
        !    The file with   <d0psi| dqV |psi>
        !
        iudpdvp_3 = 40
        filint = trim(filpun) //'.pv3'
        call diropn (iudpdvp_3, filint, lrdpdvp, exst)
     endif
     !
     ! The file containing the variation of the FermiEnergy ef_sh
#ifdef __PARA
     ! opened only by the first task of the first pool
     !

     if (me.ne.1.or.mypool.ne.1) goto 130
#endif
     iuef = 41
     filint = trim(filpun) //'.efs'
     call seqopn (iuef, filint, 'unformatted', exst)
#ifdef __PARA

130  continue
#endif

  endif
  return
end subroutine openfild3
