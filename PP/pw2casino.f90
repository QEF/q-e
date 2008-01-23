!
! Copyright (C) 2004-2007 Quantum-Espresso group
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
  CALL get_env( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM( outdir ) == ' ' ) outdir = './'
  ios = 0
  IF ( ionode )  THEN 
     !
     READ (5, inputpp, iostat=ios)
     tmp_dir = trimcheck (outdir)
     !
  END IF
  ! 
  CALL mp_bcast( ios, ionode_id ) 
  IF ( ios/=0 ) CALL errore('pw2casino', 'reading inputpp namelist', ABS(ios))
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
  USE ions_base, ONLY : nat, ntyp => nsp, ityp, tau, zv, atm
  USE cell_base, ONLY: omega, alat, tpiba2, at, bg
  USE char, ONLY: title
  USE constants, ONLY: tpi, e2
  USE ener, ONLY: ewld, ehart, etxc, vtxc, etot, etxcc, demet, ef
  USE gvect, ONLY: ngm, gstart, nr1, nr2, nr3, nrx1, nrx2, nrx3, &
       nrxx, g, gg, ecutwfc, gcutm, nl, igtongl
  USE klist , ONLY: nks, nelec, xk, wk, degauss, ngauss
  USE lsda_mod, ONLY: lsda, nspin
  USE scf, ONLY: rho, rho_core, rhog_core, vnew
  USE ldaU, ONLY : lda_plus_u, eth, Hubbard_lmax
  USE vlocal, ONLY: vloc, strf
  USE wvfct, ONLY: npw, npwx, nbnd, igk, g2kin, wg, et
  USE control_flags, ONLY : gamma_only
  USE uspp, ONLY: nkb, vkb, dvan
  USE uspp_param, ONLY: nh
  USE becmod,   ONLY: becp, calbec
  USE io_global, ONLY: stdout
  USE io_files, ONLY: nd_nmbr, nwordwfc, iunwfc
  USE wavefunctions_module, ONLY : evc
  USE funct, ONLY : dft_is_meta
  USE mp_global,            ONLY: inter_pool_comm, intra_pool_comm
  USE mp,                   ONLY: mp_sum

  IMPLICIT NONE
  INTEGER :: ig, ibnd, ik, io, na, j, ispin, nbndup, nbnddown, &
       nk, ngtot, ig7, ikk, nt, ijkb0, ikb, ih, jh, jkb, at_num 
  INTEGER, ALLOCATABLE :: idx(:), igtog(:)
  LOGICAL :: exst, found
  REAL(DP) :: ek, eloc, enl, charge, etotefield
  COMPLEX(DP), ALLOCATABLE :: aux(:)
  REAL(DP), allocatable :: v_h_new(:,:,:,:), kedtaur_new(:,:)
  INTEGER :: ios
  INTEGER, EXTERNAL :: atomic_number
  REAL (DP), EXTERNAL :: ewald, w1gauss

  CALL init_us_1
  CALL newd
  io = 77

  WRITE (6,'(/,5x,''Writing file pwfn.data for program CASINO'')')

  CALL seqopn( 77, 'pwfn.data', 'formatted',exst)  

  ALLOCATE (aux(nrxx))
  ALLOCATE (becp (nkb,nbnd))
  ! four times npwx should be enough
  ALLOCATE (idx (4*npwx) )
  ALLOCATE (igtog (4*npwx) )
  if (lda_plus_u) allocate(v_h_new(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat))
  if (dft_is_meta()) then
     allocate (kedtaur_new(nrxx,nspin))
  else
     allocate (kedtaur_new(1,nspin))
  endif
  idx(:) = 0
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

  ek  = 0.d0
  eloc= 0.d0
  enl = 0.d0
  demet=0.d0
  !
  DO ispin = 1, nspin 
     !
     !     calculate the local contribution to the total energy
     !
     !      bring rho to G-space
     !
     aux(:) = CMPLX ( rho%of_r(:,ispin), 0.d0)
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
        CALL calbec ( npw, vkb, evc, becp )
        !
        ! -TS term for metals (if any)
        !
        IF ( degauss > 0.0_dp) THEN
           DO ibnd = 1, nbnd
              demet = demet + wk (ik) * &
                 degauss * w1gauss ( (ef-et(ibnd,ik)) / degauss, ngauss)
           END DO
        END IF
        !
        DO ig =1, npw
           IF( igk(ig) > 4*npwx ) & 
                CALL errore ('pw2casino','increase allocation of index', ig)
           idx( igk(ig) ) = 1
        ENDDO
        !
        ! calculate the kinetic energy
        !
        DO ibnd = 1, nbnd
           DO j = 1, npw
              ek = ek +  CONJG(evc(j,ibnd)) * evc(j,ibnd) * &
                              g2kin(j) * wg(ibnd,ikk)
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
  call mp_sum( eloc,  intra_pool_comm )
  call mp_sum( ek,    intra_pool_comm )
  call mp_sum( ek,    inter_pool_comm )
  call mp_sum( enl,   inter_pool_comm )
  call mp_sum( demet, inter_pool_comm )
#endif
  eloc = eloc * omega 
  ek = ek * tpiba2

  ngtot = 0
  DO ig = 1, 4*npwx
     IF( idx(ig) == 1 ) THEN
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
  CALL v_of_rho( rho, rho_core, rhog_core, &
                 ehart, etxc, vtxc, eth, etotefield, charge, vnew )
  !
  etot=(ek + (etxc-etxcc)+ehart+eloc+enl+ewld)+demet
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
  IF ( degauss > 0.0_dp ) THEN
     WRITE(io,'(a)') ' Total energy (au per primitive cell; includes -TS term)'
     WRITE(io,*)etot/e2, demet/e2
  ELSE
     WRITE(io,'(a)') ' Total energy (au per primitive cell)' 
     WRITE(io,*)etot/e2                
  END IF
  WRITE(io,'(a)') ' Kinetic energy (au per primitive cell)' 
  WRITE(io,*)ek/e2              
  WRITE(io,'(a)') ' Local potential energy (au per primitive cell)' 
  WRITE(io,*)eloc/e2 
  WRITE(io,'(a)') ' Non local potential energy(au per primitive cel)'
  WRITE(io,*)enl/e2
  WRITE(io,'(a)') ' Electron electron energy (au per primitive cell)' 
  WRITE(io,*)ehart/e2    
  WRITE(io,'(a)') ' Ion ion energy (au per primitive cell)' 
  WRITE(io,*)ewld/e2
  WRITE(io,'(a)') ' Number of electrons per primitive cell'                 
  WRITE(io,*)NINT(nelec)
  WRITE(io,'(a)') ' '                 
  WRITE(io,'(a)') ' GEOMETRY'
  WRITE(io,'(a)') ' -------- '
  WRITE(io,'(a)') ' Number of atoms per primitive cell '
  WRITE(io,*) nat
  WRITE(io,'(a)')' Atomic number and position of the atoms(au) '
  DO na = 1, nat
     nt = ityp(na)
     at_num = atomic_number(TRIM(atm(nt)))
     WRITE(io,'(i6,3f20.14)') at_num, (alat*tau(j,na),j=1,3)
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
           WRITE(io,*) ibnd, ispin, et(ibnd,ikk)/e2 
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

  WRITE (stdout,*) 'Kinetic energy   '  , ek/e2
  WRITE (stdout,*) 'Local energy     ', eloc/e2
  WRITE (stdout,*) 'Non-Local energy ', enl/e2
  WRITE (stdout,*) 'Ewald energy     ', ewld/e2
  WRITE (stdout,*) 'xc contribution  ',(etxc-etxcc)/e2
  WRITE (stdout,*) 'hartree energy   ', ehart/e2
  IF ( degauss > 0.0_dp ) &
  WRITE (stdout,*) 'Smearing (-TS)   ', demet/e2
  WRITE (stdout,*) 'Total energy     ', etot/e2

  DEALLOCATE (igtog)
  DEALLOCATE (idx)
  DEALLOCATE (becp)
  DEALLOCATE (aux)

END SUBROUTINE compute_casino


