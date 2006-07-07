!
! Copyright (C) 2004 PWSCF group 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
! 
!----------------------------------------------------------------------- 
PROGRAM pw2casino
  !----------------------------------------------------------------------- 

  ! This subroutine writes the file "prefix".pwfn.data containing the 
  ! plane wave coefficients and other stuff needed by the QMC code CASINO. 

#include "f_defs.h"

  USE io_files,  ONLY : nd_nmbr, prefix, outdir, tmp_dir, trimcheck
  USE io_global, ONLY : ionode, ionode_id
  USE mp,        ONLY : mp_bcast
  !
  IMPLICIT NONE
  INTEGER :: ios

  NAMELIST / inputpp / prefix, outdir

  CALL start_postproc(nd_nmbr)
  ! 
  !   set default values for variables in namelist 
  ! 
  prefix = 'pwscf'
  outdir = './'

  IF ( ionode )  THEN 
     !
     READ (5, inputpp, err=200, iostat=ios)
200  CALL errore('pw2casino', 'reading inputpp namelist', ABS(ios))
     tmp_dir = trimcheck (outdir)
     !
  END IF
  ! 
  ! ... Broadcast variables 
  ! 
  CALL mp_bcast( prefix, ionode_id ) 
  CALL mp_bcast(tmp_dir, ionode_id ) 
  !
  CALL read_file
  CALL openfil_pp
  !
  CALL compute_casino
  !
  CALL stop_pp
  STOP

END PROGRAM pw2casino


SUBROUTINE compute_casino

  USE kinds, ONLY: DP
  USE atom, ONLY: zmesh
  USE ions_base, ONLY : nat, ntyp => nsp, ityp, tau, zv
  USE cell_base, ONLY: omega, alat, tpiba2, at, bg
  USE char, ONLY: title
  USE constants, ONLY: tpi
  USE ener, ONLY: ewld, ehart, etxc, vtxc, etot, etxcc
  USE gvect, ONLY: ngm, gstart, nr1, nr2, nr3, nrx1, nrx2, nrx3, &
       nrxx, g, gg, ecutwfc, gcutm, nl, igtongl
  USE klist , ONLY: nks, nelec, xk
  USE lsda_mod, ONLY: lsda, nspin
  USE scf, ONLY: rho, rho_core, rhog, rhog_core
  USE vlocal, ONLY: vloc, vnew, strf
  USE wvfct, ONLY: npw, npwx, nbnd, gamma_only, igk, g2kin, wg, et
  USE uspp, ONLY: nkb, vkb, dvan
  USE uspp_param, ONLY: nh
  USE becmod,   ONLY: becp 
  USE io_global, ONLY: stdout
  USE io_files, ONLY: nd_nmbr, nwordwfc, iunwfc
  USE wavefunctions_module, ONLY : evc
  IMPLICIT NONE
  INTEGER :: ig, ibnd, ik, io, na, j, ispin, nbndup, nbnddown, &
       nk, ngtot, ig7, ikk, nt, ijkb0, ikb, ih, jh, jkb 
  INTEGER, ALLOCATABLE :: INDEX(:), igtog(:)
  LOGICAL :: exst, found
  REAL(DP) :: ek, eloc, enl, charge, etotefield
  COMPLEX(DP), ALLOCATABLE :: aux(:), hpsi(:,:)
  INTEGER :: ios
  REAL (DP), EXTERNAL :: ewald

  CALL init_us_1
  CALL newd
  io = 77

  WRITE (6,'(/,5x,''Writing file pwfn.data for program CASINO'')')

  CALL seqopn( 77, 'pwfn.data', 'formatted',exst)  

  ALLOCATE (hpsi(npwx, nbnd))
  ALLOCATE (aux(nrxx))
  ALLOCATE (becp (nkb,nbnd))
  ! four times npwx should be enough
  ALLOCATE (INDEX (4*npwx) )
  ALLOCATE (igtog (4*npwx) )

  hpsi (:,:) = (0.d0, 0.d0)
  INDEX(:) = 0
  igtog(:) = 0

  IF( lsda ) THEN
     nbndup = nbnd
     nbnddown = nbnd
     nk = nks/2
     !     nspin = 2
  ELSE
     nbndup = nbnd
     nbnddown = 0
     nk = nks
     !     nspin = 1
  ENDIF

  !  if(nks > 1) rewind(iunigk)
  !  do ik=1,nks
  !     if(nks > 1) read(iunigk) npw, igk
  !     
  !  if(nks > 1) rewind(iunigk)
  ek  = 0.d0
  eloc= 0.d0
  enl = 0.d0
  !
  DO ispin = 1, nspin 
     !
     !     calculate the local contribution to the total energy
     !
     !      bring rho to G-space
     !
     aux(:) = CMPLX ( rho(:,ispin), 0.d0)
     CALL cft3(aux,nr1,nr2,nr3,nrx1,nrx2,nrx3,-1)
     !
     DO nt=1,ntyp
        DO ig = 1, ngm
           eloc = eloc + vloc(igtongl(ig),nt) * strf(ig,nt) &
                * CONJG(aux(nl(ig)))
        ENDDO
     ENDDO

     DO ik = 1, nk
        ikk = ik + nk*(ispin-1)
        CALL gk_sort (xk (1, ikk), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
        CALL davcio (evc, nwordwfc, iunwfc, ikk, - 1)
        CALL init_us_2 (npw, igk, xk (1, ikk), vkb)
        CALL ccalbec (nkb, npwx, npw, nbnd, becp, vkb, evc)

        DO ig =1, npw
           IF( igk(ig) > 4*npwx ) & 
                CALL errore ('pw2casino','increase allocation of index', ig)
           INDEX( igk(ig) ) = 1
        ENDDO
        !
        ! calculate the kinetic energy
        !
        DO ibnd = 1, nbnd
           DO j = 1, npw
              hpsi(j,ibnd) =  g2kin(j) * evc(j,ibnd)
              ek = ek +  CONJG(evc(j,ibnd))*hpsi(j,ibnd) * wg(ibnd,ikk)
           END DO

           !
           ! Calculate Non-local energy
           !
           ijkb0 = 0
           DO nt = 1, ntyp
              DO na = 1, nat
                 IF (ityp (na) .EQ.nt) THEN
                    DO ih = 1, nh (nt)
                       ikb = ijkb0 + ih
                       enl=enl+CONJG(becp(ikb,ibnd))*becp(ikb,ibnd) &
                            *wg(ibnd,ikk)* dvan(ih,ih,nt)
                       DO jh = ( ih + 1 ), nh(nt)
                          jkb = ijkb0 + jh
                          enl=enl + &
                               (CONJG(becp(ikb,ibnd))*becp(jkb,ibnd)+&
                               CONJG(becp(jkb,ibnd))*becp(ikb,ibnd))&
                               * wg(ibnd,ikk) * dvan(ih,jh,nt)

                       END DO

                    ENDDO
                    ijkb0 = ijkb0 + nh (nt)
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

#ifdef __PARA
  CALL reduce(1,eloc)
  CALL reduce(1,ek)
  CALL poolreduce(1,ek)
  CALL poolreduce(1,enl)
#endif
  eloc = eloc * omega 
  ek = ek * tpiba2

  ngtot = 0
  DO ig = 1, 4*npwx
     IF( INDEX(ig) == 1 ) THEN
        ngtot = ngtot + 1
        igtog(ngtot) = ig
     ENDIF
  ENDDO
  !
  ! compute ewald contribution
  !
  ewld = ewald( alat, nat, ntyp, ityp, zv, at, bg, tau, omega, &
       g, gg, ngm, gcutm, gstart, gamma_only, strf )
  !
  ! compute hartree and xc contribution
  !
  CALL v_of_rho( rho, rhog, rho_core, rhog_core, &
                 ehart, etxc, vtxc, etotefield, charge, vnew )
  !
  etot=(ek + (etxc-etxcc)+ehart+eloc+enl+ewld)
  !
  WRITE(io,'(a)') title
  WRITE(io,'(a)')
  WRITE(io,'(a)') ' BASIC INFO'
  WRITE(io,'(a)') ' ----------'
  WRITE(io,'(a)') ' Generated by:'
  WRITE(io,'(a)') ' PWSCF'                    
  WRITE(io,'(a)') ' Method:'
  WRITE(io,'(a)') ' DFT'
  WRITE(io,'(a)') ' DFT Functional:'         
  WRITE(io,'(a)') ' unknown'
  WRITE(io,'(a)') ' Pseudopotential'
  WRITE(io,'(a)') ' unknown'
  WRITE(io,'(a)') ' Plane wave cutoff (au)'              
  WRITE(io,*) ecutwfc/2
  WRITE(io,'(a)') ' Spin polarized:'
  WRITE(io,*)lsda 
  WRITE(io,'(a)') ' Total energy (au per primitive cell)' 
  WRITE(io,*)etot/2                
  WRITE(io,'(a)') ' Kinetic energy (au per primitive cell)' 
  WRITE(io,*)ek/2              
  WRITE(io,'(a)') ' Local potential energy (au per primitive cell)' 
  WRITE(io,*)eloc/2 
  WRITE(io,'(a)') ' Non local potential energy(au per primitive cel)'
  WRITE(io,*)enl/2
  WRITE(io,'(a)') ' Electron electron energy (au per primitive cell)' 
  WRITE(io,*)ehart/2    
  WRITE(io,'(a)') ' Ion ion energy (au per primitive cell)' 
  WRITE(io,*)ewld/2
  WRITE(io,'(a)') ' Number of electrons per primitive cell'                 
  WRITE(io,*)NINT(nelec)
  WRITE(io,'(a)') ' '                 
  WRITE(io,'(a)') ' GEOMETRY'
  WRITE(io,'(a)') ' -------- '
  WRITE(io,'(a)') ' Number of atoms per primitive cell '
  WRITE(io,*) nat
  WRITE(io,'(a)')' Atomic number and position of the atoms(au) '
  DO na = 1, nat
     WRITE(io,'(i6,3f20.14)') INT(zmesh(ityp(na))), (alat*tau(j,na),j=1,3)
  ENDDO
  WRITE(io,'(a)') ' Primitive lattice vectors (au) '
  WRITE(io,100) alat*at(1,1), alat*at(2,1), alat*at(3,1)
  WRITE(io,100) alat*at(1,2), alat*at(2,2), alat*at(3,2)
  WRITE(io,100) alat*at(1,3), alat*at(2,3), alat*at(3,3)
  WRITE(io,'(a)') ' '
  WRITE(io,'(a)') ' G VECTORS'
  WRITE(io,'(a)') ' ---------'
  WRITE(io,'(a)') ' Number of G-vectors'
  WRITE(io,*) ngtot
  WRITE(io,'(a)') ' Gx Gy Gz (au)'
  DO ig = 1, ngtot
     WRITE(io,100) tpi/alat*g(1,igtog(ig)), tpi/alat*g(2,igtog(ig)), &
          tpi/alat* g(3,igtog(ig))
  ENDDO

100 FORMAT (3(1x,f20.15))

  WRITE(io,'(a)') ' '
  WRITE(io,'(a)') ' WAVE FUNCTION'
  WRITE(io,'(a)') ' -------------'
  WRITE(io,'(a)') ' Number of k-points'
  WRITE(io,*) nk
  !  if(nks > 1) rewind(iunigk)
  
  DO ik = 1, nk
     WRITE(io,'(a)') ' k-point # ; # of bands (up spin/down spin); &
          &           k-point coords (au)'
     WRITE(io,'(3i4,3f20.16)') ik, nbndup, nbnddown, &
          (tpi/alat*xk(j,ik),j=1,3)
     DO ispin = 1, nspin 
        ikk = ik + nk*(ispin-1)
        IF( nks > 1 ) THEN
           CALL gk_sort (xk (1, ikk), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
           CALL davcio(evc,nwordwfc,iunwfc,ikk,-1)
        ENDIF
        DO ibnd = 1, nbnd
           WRITE(io,'(a)') ' Band, spin, eigenvalue (au)'
           WRITE(io,*) ibnd, ispin, et(ibnd,ikk)/2 
           WRITE(io,'(a)') ' Eigenvectors coefficients'
           DO ig=1, ngtot
              ! now for all G vectors find the PW coefficient for this k-point
              found = .FALSE.
              DO ig7 = 1, npw
                 IF( igk(ig7) == igtog(ig) )THEN
                    WRITE(io,*) evc(ig7,ibnd)
                    found = .TRUE.
                    GOTO 17
                 ENDIF
              ENDDO
              ! if can't find the coefficient this is zero
17            IF( .NOT. found ) WRITE(io,*) (0.d0, 0.d0)
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  CLOSE(io)

  WRITE (stdout,*) 'Kinetic energy   '  , ek/2
  WRITE (stdout,*) 'Local energy     ', eloc/2
  WRITE (stdout,*) 'Non-Local energy ', enl/2
  WRITE (stdout,*) 'Ewald energy     ', ewld/2
  WRITE (stdout,*) 'xc contribution  ',(etxc-etxcc)/2
  WRITE (stdout,*) 'hartree energy   ', ehart/2
  WRITE (stdout,*) 'Total energy     ', (ek + (etxc-etxcc)+ehart+eloc+enl+ewld)/2

  DEALLOCATE (igtog)
  DEALLOCATE (index)
  DEALLOCATE (becp)
  DEALLOCATE (aux)
  DEALLOCATE (hpsi)

END SUBROUTINE compute_casino


