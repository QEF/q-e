!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine elphon
  !-----------------------------------------------------------------------
  !
  ! Electron-phonon calculation from data saved in fildvscf
  !
#include "f_defs.h"
  !
  USE ions_base, ONLY : nat, ntyp => nsp, ityp, tau
  use pwcom
  USE kinds, only : DP
  use phcom
  use el_phon
  implicit none
  integer :: irr, imode0, ipert, is
  ! counter on the representations
  ! counter on the modes
  ! the change of Vscf due to perturbations
  complex(kind=DP), pointer :: dvscfin(:,:,:), dvscfins (:,:,:)

  call start_clock ('elphon')

  !
  ! read Delta Vscf and calculate electron-phonon coefficients
  !
  rewind (iudvscf)
  imode0 = 0
  do irr = 1, nirr
     allocate (dvscfin ( nrxx , nspin , npert(irr)) )
     read (iudvscf) dvscfin
     if (doublegrid) then
        allocate (dvscfins ( nrxxs , nspin , npert(irr)) )
        do is = 1, nspin
           do ipert = 1, npert(irr)
              call cinterpolate (dvscfin(1,is,ipert),dvscfins(1,is,ipert),-1)
           enddo 
        enddo
     else
        dvscfins => dvscfin
     endif
     call newdq (dvscfin, npert(irr))
     call elphel (npert (irr), imode0, dvscfins)
     !
     imode0 = imode0 + npert (irr)
     if (doublegrid) deallocate (dvscfins)
     deallocate (dvscfin)
  enddo
  !
  ! now read the eigenvalues and eigenvectors of the dynamical matrix
  ! calculated in a previous run
  !
  if (.not.trans) call readmat (iudyn, ibrav, celldm, nat, ntyp, &
       ityp, omega, amass, tau, xq, w2, dyn)
  !
  call stop_clock ('elphon')
  return
end subroutine elphon
!
!-----------------------------------------------------------------------
subroutine readmat (iudyn, ibrav, celldm, nat, ntyp, ityp, omega, &
     amass, tau, q, w2, dyn)
  !-----------------------------------------------------------------------
  !
#include "f_defs.h"
  USE kinds, only : DP
  implicit none
  ! Input
  integer :: iudyn, ibrav, nat, ntyp, ityp (nat)
  real(kind=DP) :: celldm (6), amass (ntyp), tau (3, nat), q (3), &
       omega
  ! output
  real(kind=DP) :: w2 (3 * nat)
  complex(kind=DP) :: dyn (3 * nat, 3 * nat)
  ! local (control variables)
  integer :: ntyp_, nat_, ibrav_, ityp_
  real(kind=DP) :: celldm_ (6), amass_, tau_ (3), q_ (3)
  ! local
  real(kind=DP) :: dynr (2, 3, nat, 3, nat)
  character(len=80) :: line
  character(len=3)  :: atm
  integer :: nt, na, nb, naa, nbb, nu, mu, i, j
  !
  !
  rewind (iudyn)
  read (iudyn, '(a)') line
  read (iudyn, '(a)') line
  read (iudyn, * ) ntyp_, nat_, ibrav_, celldm_
  if (ntyp.ne.ntyp_.or.nat.ne.nat_.or.ibrav_.ne.ibrav.or.abs ( &
       celldm_ (1) - celldm (1) ) .gt.1.0d-5) call errore ('readmat', &
       'inconsistent data', 1)
  do nt = 1, ntyp
     read (iudyn, * ) i, atm, amass_
     if (nt.ne.i.or.abs (amass_ - amass (nt) ) .gt.1.0d-5) call errore ( &
          'readmat', 'inconsistent data', 1 + nt)
  enddo
  do na = 1, nat
     read (iudyn, * ) i, ityp_, tau_
     if (na.ne.i.or.ityp_.ne.ityp (na) ) call errore ('readmat', &
          'inconsistent data', 10 + na)
  enddo
  read (iudyn, '(a)') line
  read (iudyn, '(a)') line
  read (iudyn, '(a)') line
  read (iudyn, '(a)') line
  read (line (11:80), * ) (q_ (i), i = 1, 3)
  read (iudyn, '(a)') line
  do na = 1, nat
     do nb = 1, nat
        read (iudyn, * ) naa, nbb
        if (na.ne.naa.or.nb.ne.nbb) call errore ('readmat', 'error reading &
             &file', nb)
        read (iudyn, * ) ( (dynr (1, i, na, j, nb), dynr (2, i, na, j, nb) &
             , j = 1, 3), i = 1, 3)
     enddo
  enddo
  !
  ! divide the dynamical matrix by the masses
  !
  do nb = 1, nat
     do j = 1, 3
        do na = 1, nat
           do i = 1, 3
              dynr (1, i, na, j, nb) = dynr (1, i, na, j, nb) / sqrt (amass ( &
                   ityp (na) ) * amass (ityp (nb) ) )
              dynr (2, i, na, j, nb) = dynr (2, i, na, j, nb) / sqrt (amass ( &
                   ityp (na) ) * amass (ityp (nb) ) )
           enddo
        enddo
     enddo
  enddo
  !
  ! solve the eigenvalue problem.
  ! NOTA BENE: eigenvectors are overwritten on dyn
  !
  call cdiagh (3 * nat, dynr, 3 * nat, w2, dyn)
  !
  ! divide by sqrt(mass) to get displacements
  !
  do nu = 1, 3 * nat
     do mu = 1, 3 * nat
        na = (mu - 1) / 3 + 1
        dyn (mu, nu) = dyn (mu, nu) / sqrt (amass (ityp (na) ) )
     enddo
  enddo
  !
  !
  return
end subroutine readmat
!
!-----------------------------------------------------------------------
subroutine elphel (npe, imode0, dvscfins)
  !-----------------------------------------------------------------------
  !
  !      Calculation of the electron-phonon matrix elements el_ph_mat
  !         <\psi(k+q)|dV_{SCF}/du^q_{i a}|\psi(k)>
  !      Original routine written by Francesco Mauri
  !
#include "f_defs.h"
  use pwcom
  USE wavefunctions_module,  ONLY: evc
  USE kinds, only : DP
  USE io_files, ONLY: iunigk
  use phcom
  use el_phon
  implicit none
  !
  integer :: npe, imode0
  complex(kind=DP) :: dvscfins (nrxxs, nspin, npe)
  ! LOCAL variables
  integer :: ik, ikk, ikq, ipert, mode, nrec, ibnd, jbnd, ir, ig, &
       ios
  complex(kind=DP) , allocatable :: aux1 (:), elphmat (:,:,:)
  complex(kind=DP) :: ZDOTC
  !
  allocate (aux1    ( nrxxs))    
  allocate (elphmat ( nbnd , nbnd , npe))    
  !
  !  Start the loops over the k-points
  !
  if (nksq.gt.1) rewind (unit = iunigk)
  do ik = 1, nksq
     if (nksq.gt.1) then
        read (iunigk, err = 100, iostat = ios) npw, igk
100     call errore ('elphel', 'reading igk', abs (ios) )
     endif
     !
     !  ik = counter of k-points with vector k
     !  ikk= index of k-point with vector k
     !  ikq= index of k-point with vector k+q
     !       k and k+q are alternated if q!=0, are the same if q=0
     !
     if (lgamma) then
        ikk = ik
        ikq = ik
        npwq = npw
     else
        ikk = 2 * ik - 1
        ikq = ikk + 1
     endif
     if (lsda) current_spin = isk (ikk)
     if (.not.lgamma.and.nksq.gt.1) then
        read (iunigk, err = 200, iostat = ios) npwq, igkq
200     call errore ('elphel', 'reading igkq', abs (ios) )
     endif
     !
     call init_us_2 (npwq, igkq, xk (1, ikq), vkb)
     !
     ! read unperturbed wavefuctions psi(k) and psi(k+q)
     !
     if (nksq.gt.1) then
        if (lgamma) then
           call davcio (evc, lrwfc, iuwfc, ikk, - 1)
        else
           call davcio (evc, lrwfc, iuwfc, ikk, - 1)
           call davcio (evq, lrwfc, iuwfc, ikq, - 1)
        endif
     endif
     !
     do ipert = 1, npe
        nrec = (ipert - 1) * nksq + ik
        !
        !  dvbare_q*psi_kpoint is read from file (if available) or recalculated
        !
        if (trans) then
           call davcio (dvpsi, lrbar, iubar, nrec, - 1)
        else
           mode = imode0 + ipert
           ! TODO : .false. or .true. ???
           call dvqpsi_us (ik, mode, u (1, mode), .false. )
        endif
        !
        ! calculate dvscf_q*psi_k
        !
        do ibnd = 1, nbnd
           aux1(:) = (0.d0, 0.d0)
           do ig = 1, npw
              aux1 (nls (igk (ig) ) ) = evc (ig, ibnd)
           enddo
           call cft3s (aux1, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, + 2)
           do ir = 1, nrxxs
              aux1 (ir) = aux1 (ir) * dvscfins (ir, current_spin, ipert)
           enddo
           call cft3s (aux1, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, - 2)
           do ig = 1, npwq
              dvpsi (ig, ibnd) = dvpsi (ig, ibnd) + aux1 (nls (igkq (ig) ) )
           enddo
        end do
        call adddvscf (ipert, ik)

        !
        ! calculate elphmat(j,i)=<psi_{k+q,j}|dvscf_q*psi_{k,i}> for this pertur
        !
        do ibnd =1, nbnd
           do jbnd = 1, nbnd
              elphmat (jbnd, ibnd, ipert) = ZDOTC (npwq, evq (1, jbnd), 1, &
                   dvpsi (1, ibnd), 1)
           enddo
           !
        enddo
     enddo
#ifdef __PARA
     call reduce (2 * nbnd * nbnd * npe, elphmat)
#endif
     !
     !  save all e-ph matrix elements into el_ph_mat
     !
     do ipert = 1, npe
        do jbnd = 1, nbnd
           do ibnd = 1, nbnd
              el_ph_mat (ibnd, jbnd, ik, ipert + imode0) = elphmat (ibnd, jbnd, ipert)
           enddo
        enddo
     enddo
  enddo
  !
  deallocate (elphmat)
  deallocate (aux1)
  !
  return
end subroutine elphel
!
!-----------------------------------------------------------------------
subroutine elphsum
  !-----------------------------------------------------------------------
  !
  !      Sum over BZ of the electron-phonon matrix elements el_ph_mat
  !      Original routine written by Francesco Mauri
  !
#include "f_defs.h"
  USE ions_base, ONLY : nat
  use pwcom
  USE kinds, only : DP
  use phcom
  use el_phon
#ifdef __PARA
  use para
#endif
  implicit none
  ! eps = 20 cm^-1, in Ry
  real(kind=DP) :: eps
  parameter (eps = 20.d0 / 13.6058d0 / 8065.5d0)
  !
  integer :: ik, ikk, ikq, isig, ibnd, jbnd, ipert, jpert, nu, mu, &
       vu, ngauss1, nsig, iuelph, ios
  real(kind=DP) :: weight, w0g1, w0g2, w0gauss, degauss1, dosef, dos_ef, &
       ef1, phase_space, lambda, gamma
  external dos_ef
  !
  complex(kind=DP) :: el_ph_sum (3*nat,3*nat)

  !
  write (6, '(5x,"electron-phonon interaction  ..."/)')
  ngauss1 = 1
  nsig = 10
  if (filelph.ne.' ') then
#ifdef __PARA
     ! parallel case: only first node writes
     if (me.ne.1) then
        iuelph = 0
     else
#endif
     iuelph = 4
     open (unit = iuelph, file = filelph, status = 'unknown', err = &
          100, iostat = ios)
100  call errore ('elphon', 'opening file'//filelph, abs (ios) )
     rewind (iuelph)
     write (iuelph, '(3f15.8,2i8)') xq, nsig, 3 * nat
     write (iuelph, '(6e14.6)') (w2 (nu) , nu = 1, nmodes)
#ifdef __PARA
     end if
#endif
  else
     iuelph = 0
  endif
  !
  !
  do isig = 1, nsig
     degauss1 = 0.01 * isig
     el_ph_sum(:,:) = (0.d0, 0.d0)
     phase_space = 0.d0
     !
     ! Recalculate the Fermi energy Ef=ef1 and the DOS at Ef, dosef = N(Ef)
     ! for this gaussian broadening
     !
     ! Note that the weights of k+q points must be set to zero for the
     ! following call to yield correct results
     !
     call efermig (et, nbnd, nks, nelec, wk, degauss1, ngauss1, ef1)
     dosef = dos_ef (ngauss1, degauss1, ef1, et, wk, nks, nbnd)
     ! N(Ef) is the DOS per spin, not summed over spin
     dosef = dosef / 2.d0
     !
     ! Sum over bands with gaussian weights
     !
     do ik = 1, nksq

        !
        ! see subroutine elphel for the logic of indices
        !
        if (lgamma) then
           ikk = ik
           ikq = ik
        else
           ikk = 2 * ik - 1
           ikq = ikk + 1
        endif
        do ibnd = 1, nbnd
           w0g1 = w0gauss ( (ef1 - et (ibnd, ikk) ) / degauss1, ngauss1) &
                / degauss1
           do jbnd = 1, nbnd
              w0g2 = w0gauss ( (ef1 - et (jbnd, ikq) ) / degauss1, ngauss1) &
                   / degauss1
              ! note that wk(ikq)=wk(ikk)
              weight = wk (ikk) * w0g1 * w0g2
              do jpert = 1, 3 * nat
                 do ipert = 1, 3 * nat
                    el_ph_sum (ipert, jpert) = el_ph_sum (ipert, jpert)  +  weight * &
                                        conjg (el_ph_mat (jbnd, ibnd, ik, ipert) ) * &
                                               el_ph_mat (jbnd, ibnd, ik, jpert)
                 enddo
              enddo
              phase_space = phase_space+weight
           enddo
        enddo

     enddo
     !
     ! el_ph_sum(mu,nu)=\sum_k\sum_{i,j}[ <psi_{k+q,j}|dvscf_q(mu)*psi_{k,i}>
     !                                  x <psi_{k+q,j}|dvscf_q(nu)*psi_{k,i}>
     !                                  x \delta(e_{k,i}-Ef) \delta(e_{k+q,j}
#ifdef __PARA
     !
     ! collect contributions from all pools (sum over k-points)
     !
     call poolreduce (2 * 3 * nat * 3 * nat, el_ph_sum)
     call poolreduce (1, phase_space)
#endif
     !
     ! symmetrize el_ph_sum(mu,nu) : it transforms as the dynamical matrix
     !
     call symdyn_munu (el_ph_sum, u, xq, s, invs, rtau, irt, irgq, at, &
          bg, nsymq, nat, irotmq, minus_q)
     !
     write (6, 9000) degauss1, ngauss1
     write (6, 9005) dosef, ef1 * 13.6058
     write (6, 9006) phase_space
     if (iuelph.ne.0) then
        write (iuelph, 9000) degauss1, ngauss1
        write (iuelph, 9005) dosef, ef1 * 13.6058
     endif
     !
     do nu = 1, nmodes
        gamma = 0.0
        do mu = 1, 3 * nat
           do vu = 1, 3 * nat
              gamma = gamma + real (conjg (dyn (mu, nu) ) * el_ph_sum (mu, vu) &
                   * dyn (vu, nu) )
           enddo
        enddo
        gamma = 3.1415926 * gamma / 2.d0
        !
        ! the factor 2 comes from the factor sqrt(hbar/2/M/omega) that appears
        ! in the definition of the electron-phonon matrix element g
        ! The sqrt(1/M) factor is actually hidden into the normal modes
        !
        ! gamma = \pi \sum_k\sum_{i,j} \delta(e_{k,i}-Ef) \delta(e_{k+q,j}-Ef)
        !         | \sum_mu z(mu,nu) <psi_{k+q,j}|dvscf_q(mu)*psi_{k,i}> |^2
        ! where z(mu,nu) is the mu component of normal mode nu (z = dyn)
        ! gamma(nu) is the phonon linewidth of mode nu
        !
        ! The factor N(Ef)^2 that appears in most formulations of el-ph interact
        ! is absent because we sum, not average, over the Fermi surface.
        ! The factor 2 is provided by the sum over spins
        !
        if (sqrt (abs (w2 (nu) ) ) .gt.eps) then
           ! lambda is the adimensional el-ph coupling for mode nu:
           ! lambda(nu)= gamma(nu)/(pi N(Ef) \omega_{q,nu}^2)
           lambda = gamma / 3.1415926 / w2 (nu) / dosef
        else
           lambda = 0.0
        endif
        ! 3.289828x10^6 is the conversion factor from Ry to GHz
        write (6, 9010) nu, lambda, gamma * 3.289828d6
        if (iuelph.ne.0) write (iuelph, 9010) nu, lambda, gamma * &
             3.289828d6
     enddo
  enddo
9000 format(5x,'Gaussian Broadening: ',f7.3,' Ry, ngauss=',i4)
9005 format(5x,'DOS =',f10.6,' states/spin/Ry/Unit Cell at Ef=', &
       &       f10.6,' eV')
9006 format(5x,'double delta at Ef =',f10.6)
9010 format(5x,'lambda(',i2,')=',f10.6,'   gamma=',f10.6,' GHz')
  !
  !
  if (iuelph.ne.0) close (unit = iuelph)
  return
end subroutine elphsum
!
!-----------------------------------------------------------------------
function dos_ef (ngauss, degauss, ef, et, wk, nks, nbnd)
  !-----------------------------------------------------------------------
  !
  USE kinds, only : DP
  implicit none
  real(kind=DP) :: dos_ef
  integer :: ngauss, nbnd, nks
  real(kind=DP) :: et (nbnd, nks), wk (nks), ef, degauss
  !
  integer :: ik, ibnd
  real(kind=DP) :: w0gauss
  !
  !     Compute DOS at E_F (states per Ry per unit cell)
  !
  dos_ef = 0.0
  do ik = 1, nks
     do ibnd = 1, nbnd
        dos_ef = dos_ef + wk (ik) * w0gauss ( (et (ibnd, ik) - ef) &
             / degauss, ngauss) / degauss
     enddo
  enddo
#ifdef __PARA
  !
  !    Collects partial sums on k-points from all pools
  !
  call poolreduce (1, dos_ef)
#endif
  !
  return
end function dos_ef
