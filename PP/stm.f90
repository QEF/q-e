!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!--------------------------------------------------------------------
subroutine stm (wf, sample_bias, z, dz, stm_wfc_matching, stmdos)
  !--------------------------------------------------------------------
  !
  !     This routine calculates an stm image defined as the local density
  !     of states at the fermi energy.
  !     To this end uses a matched exponential behaved wavefunction, in th
  !     spirit of the Tersoff and Hamann approximation (PRB 31, 805 (1985)
  !     The matching with the true wavefunction is decided by the variable
  !     in celldm(1) units, then stm images are calculateted every dz.
  !     The bias of the sample is decided by sample_bias, states between
  !     ef and ef + sample_bias are taken into account.
  !     It needs the workfunction wf. On output wf contains the number of
  !     states used to compute the image.
  !     The slab must be orientated with the main axis along celldm(3).
  !     It may not properly work if the slab has two symmetric surfaces.
  !
#include "machine.h"
  use pwcom
!
  implicit none
  logical :: stm_wfc_matching
  real(kind=DP) :: sample_bias, z, dz, stmdos (nrx1 * nrx2 * nrx3)
  ! the stm density of states
  !
  !    And here the local variables
  !

  logical :: uguale

  integer :: istates, igs, npws, ir, ir1, irx, iry, irz, ig, ibnd, &
       ik, nbnd_ocp
  ! the number of states to compute the im
  ! counter on surface G vectors
  ! number of surfac g-vectors
  ! counters on 3D r points
  ! counter on g vectors
  ! counter on bands
  ! counter on k points

  real(kind=DP) :: emin, emax, fac, wf, wf1, x, y, zz, &
       w1 , w0gauss, up, up1, down, down1, t0, scnds
  complex(kind=DP), parameter :: i= (0.d0, 1.d0)

  real(kind=DP), allocatable :: gs (:,:)
  complex(kind=DP), allocatable :: a (:), psi (:,:)
  ! the coefficients of the matching wfc
  ! plane stm wfc

  external w0gauss

  t0 = scnds ()
  allocate (gs( 2, npwx))    
  allocate (a ( npwx))    
  allocate (psi(nrx1, nrx2))    
  !
  write ( * ,  * ) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  write ( * ,  * ) '!this routine has not been tested since last changes!'
  write ( * ,  * ) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  !
  !     if matching is .true. then matches the wfc's and uses their
  !     exponential behaviour, otherwise uses the true wfc's on the fft gr
  !
  call setv (nrx1 * nrx2 * nrx3, 0.d0, stmdos, 1)

  if (.not.stm_wfc_matching) call setv (nrxx, 0.d0, rho, 1)
  if (stm_wfc_matching) then
     z = z * alat

     dz = dz * alat
     write (6, '(5x,"Matching plane at z  =",f6.2, &
          &       " alat units")') z / alat
     write (6, '(5x,"Next planes every dz =",f6.2, &
          &       " atomic units")') dz
     write (6, '(5x,"Sample bias          =",f8.4, &
          &       " eV")') sample_bias * rytoev
  else
     write (6, '(5x,"Use the true wfcs")')
     write (6, '(5x,"Sample bias          =",f8.4, &
          &       " eV")') sample_bias * rytoev

  endif
  !
  if (.not.lgauss) then
     !
     !  for semiconductors, add small broadening
     !
     nbnd_ocp = nint (nelec) / degspin

     if (nbnd.le.nbnd_ocp + 1) call errore ('stm', 'not enough bands', 1)
     emin = et (nbnd_ocp + 1, 1)
     do ik = 2, nks
        emin = min (emin, et (nbnd_ocp + 1, ik) )
     enddo
#ifdef __PARA
     ! find the minimum across pools
     call poolextreme (emin, - 1)
#endif
     emax = et (nbnd_ocp, 1)
     do ik = 2, nks
        emax = max (emax, et (nbnd_ocp, ik) )
     enddo
#ifdef __PARA
     ! find the maximum across pools
     call poolextreme (emax, 1)
#endif
     ef = (emin + emax) * 0.5d0
     degauss = 0.00001

     ngauss = 0
     write (6, '(/5x,"Occupied bands: ",i6)') nbnd_ocp
     write (6, '(/5x,"  Fermi energy: ",f10.2," eV")') ef * rytoev
     write (6, '(/5x,"    Gap energy: ",f10.2," eV")')  (emax - emin)  * rytoev
  endif
  !
  !     take only the states in the energy window above or below the fermi
  !     as determined by the bias of the sample
  !
  if (sample_bias.gt.0) then
     up = ef + sample_bias
     up1 = ef + 3 * sample_bias
     down = ef
     down1 = ef - 2 * sample_bias
  else
     up = ef
     up1 = ef - 2 * sample_bias
     down = ef + sample_bias
     down1 = ef + 3 * sample_bias

  endif
  do ik = 1, nks
     do ibnd = 1, nbnd
        if (et (ibnd, ik) .gt.down.and.et (ibnd, ik) .lt.up) then
           wg (ibnd, ik) = wk (ik)
        elseif (et (ibnd, ik) .lt.down) then
           wg (ibnd, ik) = wk (ik) * w0gauss ( (down - et (ibnd, ik) ) &
                / degauss, ngauss)
        elseif (et (ibnd, ik) .gt.up) then
           wg (ibnd, ik) = wk (ik) * w0gauss ( (up - et (ibnd, ik) ) &
                / degauss, ngauss)
        endif
     enddo


  enddo
  istates = 0
  !
  !     here we sum for each k point the contribution
  !     of the wavefunctions to the stm dos
  !
  do ik = 1, nks
     call gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
     call davcio (evc, nwordwfc, iunwfc, ik, - 1)
     !
     !     find the surface G vectors, the first is the first G of the list
     !
     if (stm_wfc_matching) then
        npws = 1
        gs (1, npws) = g (1, igk (1) )
        gs (2, npws) = g (2, igk (1) )
        do ig = 1, npw
           uguale = .false.
           do igs = 1, npws
              !
              !     if the surface part of G is equal to at least one of the surface v
              !     already found then uguale = .true.
              !
              uguale = uguale.or. (g (1, igk (ig) ) .eq.gs (1, igs) .and.g ( &
                   2, igk (ig) ) .eq.gs (2, igs) )
           enddo
           !
           !     if G is not equal to any surface vector then G is a new surface ve
           !
           if (.not.uguale) then
              npws = npws + 1
              gs (1, npws) = g (1, igk (ig) )
              gs (2, npws) = g (2, igk (ig) )
           endif
        enddo
     endif
     !
     do ibnd = 1, nbnd
        w1 = wg (ibnd, ik) / omega
        if (et (ibnd, ik) .gt.up1.or.et (ibnd, ik) .lt.down1) goto 10
        write (6, * ) w1, ibnd, ik
        !
        !     istates is a counter on the states used to compute the image
        !
        istates = istates + 1
        !
        !     for this state the work function is modified accordingly to its en
        !
        !
        !     find the coefficients of the matching wfcs
        !
        if (stm_wfc_matching) then
           wf1 = wf - (et (ibnd, ik) - ef - sample_bias)
           do igs = 1, npws
              a (igs) = (0.d0, 0.d0)
              fac = exp (z * sqrt (wf1 + ( (xk (1, ik) + gs (1, igs) ) **2 + &
                   (xk (2, ik) + gs (2, igs) ) **2) * tpiba2) )
              do ig = 1, npw
                 !
                 !     sum over the z-component of the G vector
                 !
                 if (g (1, igk (ig) ) .eq.gs (1, igs) .and.g (2, igk (ig) ) &
                      .eq.gs (2, igs) ) then
                    a (igs) = a (igs) + evc (ig, ibnd) * fac * exp (i * z * &
                         (xk (3, ik) + g (3, igk (ig) ) ) * tpiba)
                 endif
              enddo
           enddo
           !
           !     reconstruct the wfc for the z in interest for this k-point and thi
           !     uses nr3/2 planes only, the other nr3/2 are empty -> only the uppe
           !     surface is used. N.B. it may not properly work if the upper surfac
           !     is connected to the lower surface by a symmetry operation (one
           !     should take the average of the two surfaces...).
           !
           do irz = 2, nr3 / 2
              !
              !     zz is the new z
              !
              zz = z + dz * (irz - 2)
              call setv (2 * nrx1 * nrx2, 0.d0, psi, 1)
              do igs = 1, npws
                 fac = exp ( - sqrt (wf1 + ( (xk (1, ik) + gs (1, igs) ) **2 + &
                      (xk (2, ik) + gs (2, igs) ) **2) * tpiba2) * zz)
                 do iry = 1, nr2
                    do irx = 1, nr1
                       !
                       !     works for the z axis orthogonal to the xy plane
                       !
                       x = at (1, 1) * real (irx - 1) / nr1 + at (1, 2) * real (iry - &
                            1) / nr2
                       y = at (2, 1) * real (irx - 1) / nr1 + at (2, 2) * real (iry - &
                            1) / nr2
                       !
                       !     psi is the wfc in the plane xy at height zz
                       !
                       psi (irx, iry) = psi (irx, iry) + a (igs) * fac * exp (tpi * i &
                            * ( (xk (1, ik) + gs (1, igs) ) * x + (xk (2, ik) + gs (2, igs) &
                            ) * y) )
                    enddo
                 enddo
              enddo
#ifdef __PARA
              call reduce (2 * nrx1 * nrx2, psi)
#endif
              !
              !     now sum for each k-point and for each band the square modulus of t
              !     wfc times the weighting factor
              !
              do iry = 1, nr2
                 do irx = 1, nr1
                    ir = irx + (iry - 1) * nrx1 + (irz - 1) * nrx1 * nrx2
                    stmdos (ir) = stmdos (ir) + w1 * psi (irx, iry) * conjg (psi ( &
                         irx, iry) )
                 enddo
              enddo
           enddo
        else
           !
           !     do not match
           !
           call setv (2 * nrxx, 0.d0, psic, 1)
           do ig = 1, npw
              psic (nl (igk (ig) ) ) = evc (ig, ibnd)

           enddo

           call cft3 (psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1)
           do ir = 1, nrxx
              rho (ir, 1) = rho (ir, 1) + w1 * (real (psic (ir) ) **2 + &
                   DIMAG (psic (ir) ) **2)
           enddo
        endif
10      continue
     enddo
  enddo
  !
  !     symmetrization of the stm dos
  !
#ifdef __PARA
  if (stm_wfc_matching) then
     call poolreduce (nrx1 * nrx2 * nrx3, stmdos)
     call symrho (stmdos, nrx1, nrx2, nrx3, nr1, nr2, nr3, nsym, s, &
          ftau)
  else
     call poolreduce (nrxx, rho)
     call psymrho (rho, nrx1, nrx2, nrx3, nr1, nr2, nr3, nsym, s, &
          ftau)
     call gather (rho, stmdos)
  endif
#else
  if (.not.stm_wfc_matching) call DCOPY (nrxx, rho, 1, stmdos, 1)
  call symrho (stmdos, nrx1, nrx2, nrx3, nr1, nr2, nr3, nsym, s, ftau)
#endif
  deallocate(psi)
  deallocate(a)
  deallocate(gs)
  write (6, '(/5x,"stm took ",f10.2," cpu secs")') scnds ()-t0
  !
  !     use wf to store istates
  !
  wf = istates
#ifdef __PARA
  call poolreduce (1, wf)
#endif
  z = z / alat

  dz = dz / alat
  return
end subroutine stm
