! Copyright (C) 2004 PWSCF group 
! This file is distributed under the terms of the 
! GNU General Public License. See the file `License' 
! in the root directory of the present distribution, 
! or http://www.gnu.org/copyleft/gpl.txt . 
! 
!----------------------------------------------------------------------- 
program pw2casino
  !----------------------------------------------------------------------- 

  ! This subroutine writes the file pwfn.data containing the plane wave
  ! coefficients and other stuff needed by the QMC code CASINO. 

  ! #include "machine.h"

  use io_files, only: nd_nmbr, prefix, outdir
#ifdef __PARA 
  use para,       only : me 
  use mp, only: mp_bcast
#endif 
  implicit none
  integer :: ios, ionode_id = 0

  namelist / inputpp / prefix, outdir

  call start_postproc(nd_nmbr)
  ! 
  !   set default values for variables in namelist 
  ! 
  prefix = 'pwscf'
  outdir = './'
#ifdef __PARA 
  if (me == 1)  then 
#endif 
     read (5, inputpp, err=200, iostat=ios)
200  call errore('pw2casino', 'reading inputpp namelist', abs(ios))
#ifdef __PARA 
  end if
  ! 
  ! ... Broadcast variables 
  ! 
  CALL mp_bcast( prefix, ionode_id ) 
#endif 
  !
  call read_file
  call openfil
  !
  call compute_casino
  !
  call stop_pp
  stop

end program pw2casino


subroutine compute_casino

  use kinds, ONLY: DP
  use atom, only: zmesh
  use basis, only: nat, ntyp, ityp, tau
  use brilz, only: omega, alat, tpiba2, at, bg
  use char, only: title
  use constants, only: tpi
  use ener, only: ewld, ehart, etxc, vtxc, etot, etxcc
  use gvect, only: ngm, gstart, nr1, nr2, nr3, nrx1, nrx2, nrx3, &
       nrxx, g, gg, ecutwfc, gcutm, nl, igtongl
  use klist , only: nks, nelec, xk
  use lsda_mod, only: lsda, nspin
  use pseud, only: zv
  use scf, only: rho, rho_core
  use vlocal, only: vloc, vnew, strf
  use wvfct, only: npw, npwx, nbnd, gamma_only, igk, g2kin, wg, et
  use us, only: nkb, vkb, nh, dvan
  use becmod,   only: becp 
  use io_global, only: stdout
  use io_files, only: nd_nmbr, nwordwfc, iunwfc
  USE wavefunctions_module, ONLY : evc
  implicit none
  integer :: ig, ibnd, ik, io, na, j, ispin, nbndup, nbnddown, &
       nk, ngtot, ig7, ikk, nt, ijkb0, ikb, ih, jh, jkb 
  integer, allocatable :: index(:), igtog(:)
  logical :: exst, found
  real(kind=DP) :: ek, eloc, enl,charge
  complex(kind=DP), allocatable :: aux(:), hpsi(:,:)
  integer :: ios
  REAL (KIND=DP), EXTERNAL :: ewald

  call init_us_1
  call newd
  io = 77

  write (6,'(/,5x,''Writing file pwfn.data for program CASINO'')')

  call seqopn( 77, 'pwfn.data', 'formatted',exst)  

  allocate (hpsi(npwx, nbnd))
  allocate (aux(nrxx))
  ! four times npwx should be enough
  allocate (index (4*npwx) )
  allocate (igtog (4*npwx) )

  hpsi (:,:) = (0.d0, 0.d0)
  index(:) = 0
  igtog(:) = 0

  if( lsda ) then
     nbndup = nbnd
     nbnddown = nbnd
     nk = nks/2
     !     nspin = 2
  else
     nbndup = nbnd
     nbnddown = 0
     nk = nks
     !     nspin = 1
  endif

  !  if(nks > 1) rewind(iunigk)
  !  do ik=1,nks
  !     if(nks > 1) read(iunigk) npw, igk
  !     
  !  if(nks > 1) rewind(iunigk)
  ek  = 0.d0
  eloc= 0.d0
  enl = 0.d0
  do ispin = 1, nspin 
     !
     !     calculate the local contribution to the total energy
     !
     !      bring rho to G-space
     !
     aux(:) = DCMPLX ( rho(:,ispin), 0.d0)
     call cft3(aux,nr1,nr2,nr3,nrx1,nrx2,nrx3,-1)
     !
     do nt=1,ntyp
        do ig = gstart, ngm
           eloc = eloc + vloc(igtongl(ig),nt) * strf(ig,nt) &
                * conjg(aux(nl(ig)))
        enddo
     enddo

     do ik = 1, nk
        ikk = ik + nk*(ispin-1)
        call gk_sort (xk (1, ikk), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
        call davcio (evc, nwordwfc, iunwfc, ikk, - 1)
        call init_us_2 (npw, igk, xk (1, ikk), vkb)
        call ccalbec (nkb, npwx, npw, nbnd, becp, vkb, evc)

        do ig =1, npw
           if( igk(ig) > 4*npwx ) & 
                call errore ('pw2casino','increase allocation of index', ig)
           index( igk(ig) ) = 1
        enddo

        !
        ! calculate the kinetic energy
        !
        do ibnd = 1, nbnd
           do j = 1, npw
              hpsi(j,ibnd) =  g2kin(j) * evc(j,ibnd)
              ek = ek +  conjg(evc(j,ibnd))*hpsi(j,ibnd) * wg(ibnd,ikk)
           end do

           !
           ! Calculate Non-local energy
           !
           ijkb0 = 0
           do nt = 1, ntyp
              do na = 1, nat
                 if (ityp (na) .eq.nt) then
                    do ih = 1, nh (nt)
                       ikb = ijkb0 + ih
                       enl=enl+conjg(becp(ikb,ibnd))*becp(ikb,ibnd) &
                            *wg(ibnd,ikk)* dvan(ih,ih,nt)
                       DO jh = ( ih + 1 ), nh(nt)
                          jkb = ijkb0 + jh
                          enl=enl + &
                               (conjg(becp(ikb,ibnd))*becp(jkb,ibnd)+&
                               conjg(becp(jkb,ibnd))*becp(ikb,ibnd))&
                               * wg(ibnd,ikk) * dvan(ih,jh,nt)

                       END DO

                    enddo
                    ijkb0 = ijkb0 + nh (nt)
                 endif
              enddo
           enddo
        enddo
     enddo
  enddo

#ifdef __PARA
  call reduce(1,eloc)
  call reduce(1,ek)
  call poolreduce(1,ek)
  call poolreduce(1,enl)
#endif
  eloc = eloc * omega 
  ek = ek * tpiba2

  ngtot = 0
  do ig = 1, 4*npwx
     if( index(ig) == 1 ) then
        ngtot = ngtot + 1
        igtog(ngtot) = ig
     endif
  enddo
  !
  ! compute ewald contribution
  !
  ewld = ewald( alat, nat, ntyp, ityp, zv, at, bg, tau, omega, &
       g, gg, ngm, gcutm, gstart, gamma_only, strf )
  !
  ! compute hartree and xc contribution
  !
  CALL v_of_rho( rho, rho_core, nr1, nr2, nr3, nrx1, nrx2, nrx3, &
       nrxx, nl, ngm, gstart, nspin, g, gg, alat, omega, &
       ehart, etxc, vtxc, charge, vnew )
  !
  etot=(ek + (etxc-etxcc)+ehart+eloc+enl+ewld)
  !
  write(io,'(a)') title
  write(io,'(a)')
  write(io,'(a)') ' BASIC INFO'
  write(io,'(a)') ' ----------'
  write(io,'(a)') ' Generated by:'
  write(io,'(a)') ' PWSCF'                    
  write(io,'(a)') ' Method:'
  write(io,'(a)') ' DFT'
  write(io,'(a)') ' DFT Functional:'         
  write(io,'(a)') ' unknown'
  write(io,'(a)') ' Pseudopotential'
  write(io,'(a)') ' unknown'
  write(io,'(a)') ' Plane wave cutoff (au)'              
  write(io,*) ecutwfc/2
  write(io,'(a)') ' Spin polarized:'
  write(io,*)lsda 
  write(io,'(a)') ' Total energy (au per primitive cell)' 
  write(io,*)etot/2                
  write(io,'(a)') ' Kinetic energy (au per primitive cell)' 
  write(io,*)ek/2              
  write(io,'(a)') ' Local potential energy (au per primitive cell)' 
  write(io,*)eloc/2 
  write(io,'(a)') ' Non local potential energy(au per primitive cel)'
  write(io,*)enl/2
  write(io,'(a)') ' Electron electron energy (au per primitive cell)' 
  write(io,*)ehart/2    
  write(io,'(a)') ' Ion ion energy (au per primitive cell)' 
  write(io,*)ewld/2
  write(io,'(a)') ' Number of electrons per primitive cell'                 
  write(io,*)nint(nelec)
  write(io,'(a)') ' '                 
  write(io,'(a)') ' GEOMETRY'
  write(io,'(a)') ' -------- '
  write(io,'(a)') ' Number of atoms per primitive cell '
  write(io,*) nat
  write(io,'(a)')' Atomic number and position of the atoms(au) '
  do na = 1, nat
     write(io,'(i6,3f20.14)') int(zmesh(ityp(na))), (alat*tau(j,na),j=1,3)
  enddo
  write(io,'(a)') ' Primitive lattice vectors (au) '
  write(io,100) alat*at(1,1), alat*at(2,1), alat*at(3,1)
  write(io,100) alat*at(1,2), alat*at(2,2), alat*at(3,2)
  write(io,100) alat*at(1,3), alat*at(2,3), alat*at(3,3)
  write(io,'(a)') ' '
  write(io,'(a)') ' G VECTORS'
  write(io,'(a)') ' ---------'
  write(io,'(a)') ' Number of G-vectors'
  write(io,*) ngtot
  write(io,'(a)') ' Gx Gy Gz (au)'
  do ig = 1, ngtot
     write(io,100) tpi/alat*g(1,igtog(ig)), tpi/alat*g(2,igtog(ig)), &
          tpi/alat* g(3,igtog(ig))
  enddo

100 FORMAT (3(1x,f20.15))

  write(io,'(a)') ' '
  write(io,'(a)') ' WAVE FUNCTION'
  write(io,'(a)') ' -------------'
  write(io,'(a)') ' Number of k-points'
  write(io,*) nks
  !  if(nks > 1) rewind(iunigk)
  do ispin = 1, nspin 
     do ik = 1, nks
        ikk = ik + nks*(ispin-1)
        if( nks > 1 ) then
           !           read(iunigk) npw, igk
           call gk_sort (xk (1, ikk), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
           call davcio(evc,nwordwfc,iunwfc,ikk,-1)
        endif

        write(io,'(a)') ' k-point # ; # of bands (up spin/down spin); &
             &           k-point coords (au)'
        write(io,'(3i4,3f20.16)') ik, nbndup, nbnddown, &
             (tpi/alat*xk(j,ikk),j=1,3)
        do ibnd = 1, nbnd
           write(io,'(a)') ' Band, spin, eigenvalue (au)'
           write(io,*) ibnd, ispin, et(ibnd,ikk)/2 
           write(io,'(a)') ' Eigenvectors coefficients'
           do ig=1, ngtot
              ! now for all G vectors find the PW coefficient for this k-point
              found = .false.
              do ig7 = 1, npw
                 if( igk(ig7) == igtog(ig) )then
                    write(io,*) evc(ig7,ibnd)
                    found = .true.
                    goto 17
                 endif
              enddo
              ! if can't find the coefficient this is zero
17            if( .not. found ) write(io,*) (0.d0, 0.d0)
           enddo
        enddo
     enddo
  enddo

  write (stdout,*) 'Kinetic energy   '  , ek/2
  write (stdout,*) 'Local energy     ', eloc/2
  write (stdout,*) 'Non-Local energy ', enl/2
  write (stdout,*) 'Ewald energy     ', ewld/2
  write (stdout,*) 'xc contribution  ',(etxc-etxcc)/2
  write (stdout,*) 'hartree energy   ', ehart/2
  write (stdout,*) 'Total energy     ', (ek + (etxc-etxcc)+ehart+eloc+enl+ewld)/2


end subroutine compute_casino


