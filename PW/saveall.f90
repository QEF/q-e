!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine saveall (iun, iflag)
  !-----------------------------------------------------------------------
  !
  !      This subroutine writes to (iflag=1) or reads from (iflag!=1) unit
  !      iun all quantities needed in subsequent calculations that are not
  !      i) dynamically allocated, ii) distributed in parallel execution
  !
  use pwcom
  use io
  use funct
  implicit none

  integer :: iun, iflag
  ! units where reads or writes
  ! if 1 writes otherwise reads
  integer :: ios
  ! integer variable for I/O control
  character (len=55) :: dummy_tmp_dir
  !
  !
  if (iflag.eq.1) then
     write (iun) celldm, at, bg, alat, omega, tpiba, tpiba2, ibrav, symm_type
     write (iun) iswitch, istep, modenum
     write (iun) nat, ntyp, nbnd, npwx, nbndx, natomwfc
     write (iun) nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, gcutm, ecutwfc, dual
     write (iun) nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, nrxxs, &
          doublegrid, gcutms
     write (iun) ef, degauss, nelec, nks, nkstot, ngauss, lgauss
     write (iun) k1, k2, k3, nk1, nk2, nk3, ntetra, ltetra
     write (iun) s, ftau, nsym, invsym, noinv
     write (iun) zmesh, xmin, dx, r, rab, vnl, chi, oc, rho_at, &
          rho_atc, mesh, msh, nchi, lchi, numeric
     write (iun) cc, alpc, zp, aps, alps, zv, nlc, nnl, lmax, lloc, bhstype
     write (iun) dion, betar, qqq, qfunc, qfcoef, rinner, nbeta, &
          kkbeta, nqf, nqlc, ifqopt, lll, iver, tvanp, okvan
     write (iun) newpseudo
     write (iun) iexch, icorr, igcx, igcc, lsda
     write (iun) a_nlcc, b_nlcc, alpha_nlcc, nlcc
     write (iun) lforce
     write (iun) tefield 
     write (iun) edir
     write (iun) emaxpos, eopreg, eamp
     write (iun) title, crystal, atm, psd, sname, tmp_dir
  else
     read (iun, err = 100, iostat = ios) celldm, at, bg, alat, &
          omega, tpiba, tpiba2, ibrav, symm_type
     read (iun, err = 100, iostat = ios) iswitch, istep, modenum
     read (iun, err = 100, iostat = ios) nat, ntyp, nbnd, npwx, nbndx, natomwfc
     read (iun, err = 100, iostat = ios) nr1, nr2, nr3, nrx1, nrx2, &
          nrx3, nrxx, gcutm, ecutwfc, dual
     read (iun, err = 100, iostat = ios) nr1s, nr2s, nr3s, nrx1s, &
          nrx2s, nrx3s, nrxxs, doublegrid, gcutms
     read (iun, err = 100, iostat = ios) ef, degauss, nelec, nks, nkstot,&
          ngauss, lgauss
     read (iun, err = 100, iostat = ios) k1, k2, k3, nk1, nk2, nk3, &
          ntetra, ltetra
     read (iun, err = 100, iostat = ios) s, ftau, nsym, invsym, noinv
     read (iun, err = 100, iostat = ios) zmesh, xmin, dx, r, rab, &
          vnl, chi, oc, rho_at, rho_atc, mesh, msh, nchi, lchi, numeric
     read (iun, err = 100, iostat = ios) cc, alpc, zp, aps, alps, &
          zv, nlc, nnl, lmax, lloc, bhstype
     read (iun, err = 100, iostat = ios) dion, betar, qqq, qfunc, &
          qfcoef, rinner, nbeta, kkbeta, nqf, nqlc, ifqopt, lll, iver, &
          tvanp, okvan
     read (iun, err = 100, iostat = ios) newpseudo
     read (iun, err = 100, iostat = ios) iexch, icorr, igcx, igcc, lsda
     read (iun, err = 100, iostat = ios) a_nlcc, b_nlcc, alpha_nlcc, nlcc
     read (iun, err = 100, iostat = ios) lforce
     read (iun, err = 100, iostat = ios) tefield 
     read (iun, err = 100, iostat = ios) edir
     read (iun, err = 100, iostat = ios) emaxpos, eopreg, eamp
     read (iun, err = 100, iostat = ios) title, crystal, atm, psd, &
          sname, dummy_tmp_dir
     !
100  call errore ('saveall', 'reading file', abs (ios) )
  endif
  return
end subroutine saveall

