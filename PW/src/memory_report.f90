!
! Copyright (C) 2007-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE memory_report()
  !----------------------------------------------------------------------------
  !
  ! Rough estimate of the dynamical memory allocated by the pw.x code
  ! Should be called after the first steps of initialization are done,
  ! but before large arrays are actually allocated. 
  ! Not guaranteed to be accurate for all cases (especially exotic ones).
  ! Originally written by PG, much improved by Pietro Bonfa' with:
  ! * Detailed memory report in verbose mode
  ! * Memory buffers for LDA+U projectors now included.
  ! * Local potential and structure factors now included
  !  (small but sometimes not negligible).
  ! * Q functions (qrad) now included.
  ! * Initial estimate of memory used in runs with real space augmentation.
  ! * Added psi and vc (or vr) in diagonalization routines.
  ! * Added memory allocated during initialization (wfc_init + rotate_wfc).
  ! * Added memory used for force evaluation.
  ! * Added a few comments.
  !
  USE io_global, ONLY : stdout
  USE kinds,     ONLY : dp 
  USE constants, ONLY : tpi, fpi, pi, eps16
  USE wvfct,     ONLY : nbnd, nbndx
  USE basis,     ONLY : natomwfc, starting_wfc
  USE cell_base, ONLY : omega, bg, alat
  USE exx,       ONLY : ecutfock, nkqs, use_ace
  USE fft_base,  ONLY : dffts, dfftp
  USE gvect,     ONLY : ngm, ngl, ngm_g, g, gcutm
  USE gvecs,     ONLY : ngms, doublegrid
  USE gvecw,     ONLY : ecutwfc, gcutw
  USE klist,     ONLY : nks, nkstot, xk, qnorm
  USE cellmd,    ONLY : cell_factor
  USE uspp,      ONLY : nkb, okvan
  USE atom,      ONLY : rgrid
  USE funct,     ONLY : dft_is_meta, dft_is_hybrid
  USE ldaU,      ONLY : lda_plus_u, U_projection, nwfcU
  USE fixed_occ, ONLY : one_atom_occupations
  USE wannier_new,ONLY: use_wannier
  USE lsda_mod,  ONLY : nspin
  USE uspp_param,ONLY : lmaxkb, upf, nh, nbetam
  USE us,        ONLY : dq
  USE noncollin_module, ONLY : npol, nspin_mag
  USE control_flags,    ONLY: isolve, nmix, imix, gamma_only, lscf, io_level, &
       lxdm, smallmem, tqr, iverbosity
  USE force_mod, ONLY : lforce, lstres
  USE ions_base, ONLY : nat, ntyp => nsp, ityp
  USE mp_diag,   ONLY : np_ortho
  USE mp_bands,  ONLY : nproc_bgrp, nbgrp
  USE mp_pools,  ONLY : npool
  USE mp_images, ONLY : nproc_image  
  !
  IMPLICIT NONE
  !
  INTEGER, PARAMETER :: MB=1024*1024
  INTEGER, PARAMETER :: GB=1024*MB
  INTEGER :: g_fact, mix_type_size, scf_type_size
  INTEGER :: nk, nbnd_l, npwx_g, npwx_l, ngxx_g, nexx_l, ngm_l
  INTEGER :: maxbnd, maxnab, maxnij, nab, na, nij, nt, lmaxq, nqxq
  INTEGER :: indm, ijv, roughestimate
  REAL(DP):: mbr, mbx, mby, mbz, dmbx, dmby, dmbz
  !
  INTEGER, EXTERNAL :: n_plane_waves
  !
  ! these quantities are real in order to prevent integer overflow
  !
  REAL(dp), PARAMETER :: complex_size=16_dp, real_size=8_dp, int_size=4_dp
  REAL(dp) :: ram, ram_, ram1, ram2, maxram, totram, add
  !
  IF ( gamma_only) THEN
     g_fact = 2  ! use half plane waves or G-vectors
  ELSE
     g_fact = 1  ! use  all plane waves or G-vectors
  END IF
  !
  ! npwx (max number of plane waves) is roughly estimated from the volume
  ! of the G-vector sphere, divided by the volume of the Brillouin Zone
  ! (the exact value can be computed only after G-vectors are computed)
  ! npwx_l is npwx on this processor (local)
  ! npwx_g is npwx summed over all processors (global)
  !
  npwx_g = NINT ( fpi/3.0_dp * SQRT(ecutwfc)**3 / (tpi**3/omega) / g_fact )
  npwx_l = npwx_g/nproc_bgrp
  !
  ! ram   = dynamically (permanently) allocated RAM, per process
  ! maxram= "high watermark": max ram needed during a run
  ! totram= max ram needed summed over all processors
  !
  !=====================================================================
  ! Wavefunctions (including buffer space)
  !
  !                files: allocate_wfc.f90:28,30,32
  !                       orthoatwfc.f90:102
  !=====================================================================
  IF ( io_level > 0 .OR. nks == 1) THEN
     nk = 1      ! store just one wavefunction array at the time
  ELSE
     nk = nks+1  ! store nks wavefunctions in memory (buffers)
  END IF
  ram =  complex_size * nbnd * npol * npwx_l * nk
  IF ( iverbosity > 0 ) WRITE( stdout, 1013 ) 'wfc', ram/nk/MB
  IF ( iverbosity > 0 ) WRITE( stdout, 1013 ) 'wfc (w. buffer)', ram/MB
  ! atomic wavefunctions 
  IF ( one_atom_occupations .OR. use_wannier ) THEN
     add = complex_size * natomwfc * npol * npwx_l 
     IF ( iverbosity > 0 ) WRITE( stdout, 1013 ) 'atomic wfc', add/MB
     ram = ram + add
  END IF
  ! Hubbard wavefunctions
  IF ( lda_plus_u .AND. U_projection .NE. 'pseudo' ) THEN
     add = complex_size * nwfcU * npol * npwx_l * nk ! also buffer 
     IF ( iverbosity > 0 ) WRITE( stdout, 1013 ) 'U proj.', add/nk/MB
     IF ( iverbosity > 0 ) WRITE( stdout, 1013 ) 'U proj. (w. buff.)', add/MB
     ram = ram + add
  END IF
  !
  ! hybrid functionals
  IF ( dft_is_hybrid () ) THEN
     ! ngxx_g = estimated global number of G-vectors used in V_x\psi products
     ! nexx_l = estimated local size of the FFT grid used in V_x\psi products
     ngxx_g = NINT ( fpi/3.0_dp * SQRT(ecutfock)**3 / (tpi**3/omega) )
     nexx_l = 16*ngxx_g/nproc_bgrp
     ! nbnd_l : estimated number of bands per proc with band parallelization
     nbnd_l = NINT( DBLE(nbnd) / nbgrp )
     ! Stored wavefunctions in real space 
     add = complex_size/g_fact * nexx_l * npol * nbnd_l * nkqs
     ! ACE Projectors
     IF (use_ace) add = add + complex_size * npwx_l * npol * nbnd * nks
     IF ( iverbosity > 0 ) WRITE( stdout, 1013 ) 'EXX', add/MB
     ram = ram + add
  END IF
  !=====================================================================
  !
  !=====================================================================
  ! Local pseudopotential, allocate_locpot.f90
  !=====================================================================
  add = complex_size * ngm * ntyp ! structure factor
  IF ( iverbosity > 0 ) WRITE( stdout, 1013 ) 'str. fact', add/MB
  ram = ram + add
  add = real_size * ngl * ntyp    ! local 
  IF ( iverbosity > 0 ) WRITE( stdout, 1013 ) 'local pot', add/MB
  ram = ram + add
  !=====================================================================
  !
  !=====================================================================
  ! Nonlocal pseudopotentials V_NL (beta functions), reciprocal space
  !=====================================================================
  add = complex_size * nkb * npwx_l ! allocate_nlpot.f90:88 vkb
  IF ( iverbosity > 0 ) WRITE( stdout, 1013 ) 'nlocal pot', add/MB
  ram = ram + add
  ! other (possibly minor) data loads
  lmaxq = 2*lmaxkb+1
  IF (lmaxq > 0) THEN
     ! not accurate if spline_ps .and. cell_factor <= 1.1d0
     nqxq = int( ( (sqrt(gcutm) + qnorm) / dq + 4) * cell_factor )
     ! allocate_nlpot.f90:87 qrad
     add = real_size * nqxq * nbetam*(nbetam+1)/2 * lmaxq * ntyp
     IF ( iverbosity > 0 ) WRITE( stdout, 1013 ) 'qrad', add/MB
     ram = ram + add
  END IF
  !=====================================================================
  !
  !=====================================================================
  ! Charge density and potentials - see scf_type in scf_mod
  !=====================================================================
  scf_type_size =  (complex_size * ngm + real_size * dfftp%nnr ) * nspin ! scf_mod.f90:94-95
  IF ( dft_is_meta() .or. lxdm ) scf_type_size =  2 * scf_type_size
  ! rho, v, vnew (allocate_fft.f90:56) 
  add = 3 * scf_type_size
  IF ( iverbosity > 0 ) WRITE( stdout, 1013 ) 'rho,v,vnew', add/MB
  ram =  ram + add
  
  ! vltot, vrs, rho_core, rhog_core, psic, strf, kedtau if needed
  ram =  ram + complex_size * ( dfftp%nnr + ngm *( 1 + ntyp ) ) + &
       real_size * dfftp%nnr*(2+nspin)
  IF ( dft_is_meta() ) ram = ram + real_size * dfftp%nnr*nspin
  ! arrays for rho mixing
  IF ( lscf ) THEN
     ! rhoin (electrons.f90:439)
     ram =  ram + scf_type_size
     IF ( iverbosity > 0 ) WRITE( stdout, 1013 ) 'rhoin', DBLE(scf_type_size)/MB
     ! see mix_type in scf_mod
     mix_type_size =  complex_size * ngm * nspin
     IF ( dft_is_meta() .or. lxdm ) mix_type_size =  2 * mix_type_size
     ! df, dv (if kept in memory)
     IF ( io_level < 2 ) THEN
        add = mix_type_size * 2 * nmix
        IF ( iverbosity > 0 ) WRITE( stdout, 1013 ) 'rho*nmix', add/MB
        ram = ram + add
     END IF
  END IF
  !=====================================================================
  !
  !=====================================================================
  ! G-vectors: g, gg, mill, nl, nlm, ig_l2g, igtongl
  !=====================================================================
  add = real_size * ngm * 4 + int_size * ngm * 7
  ! double grid: nls, nlsm
  add = add + int_size * ngms * 2
  IF ( iverbosity > 0 ) WRITE( stdout, 1013 ) 'G-vectors', add/MB
  ram = ram + add
  !=====================================================================
  !  
  !=====================================================================
  ! Real treatment of augmentation charge.
  !=====================================================================
  IF ( okvan .and. tqr ) THEN
    ! realus.f90:422
    mbr = 0.d0
    DO nt = 1, ntyp
        IF ( .not. upf(nt)%tvanp ) CYCLE
        DO ijv = 1, upf(nt)%nbeta*(upf(nt)%nbeta+1)/2
          DO indm = upf(nt)%mesh,1,-1
              !
              IF ( maxval(abs( upf(nt)%qfuncl(indm,ijv,:) )) > eps16 ) THEN
                mbr = max( mbr,  rgrid(nt)%r(indm) )
                exit
              ENDIF
              !
          ENDDO
        ENDDO
    END DO
    mbr = mbr / alat

    mbx = mbr*sqrt( bg(1,1)**2 + bg(1,2)**2 + bg(1,3)**2 )
    mby = mbr*sqrt( bg(2,1)**2 + bg(2,2)**2 + bg(2,3)**2 )
    mbz = mbr*sqrt( bg(3,1)**2 + bg(3,2)**2 + bg(3,3)**2 )
    !
    dmbx = 2*anint( mbx*dfftp%nr1x ) + 2
    dmby = 2*anint( mby*dfftp%my_nr2p ) + 2  ! approximation!
    dmbz = 2*anint( mbz*dfftp%my_nr3p ) + 2  ! approximation!
    !
    roughestimate = anint( dble( dmbx*dmby*dmbz ) * pi / 6.D0 )
    add  = 0.d0
    DO na = 1, nat
      nt = ityp(na)
      IF ( .not. upf(nt)%tvanp ) CYCLE
      !                          mbia             nfuncs
      add = add + real_size*roughestimate*(nh(nt)*(nh(nt)+1)/2)
    END DO
    IF ( iverbosity > 0 ) WRITE( stdout, 1013 ) 'qr (very rough)', add/MB
    ram = ram + add
  END IF 
  !=====================================================================
  !
  ! compute ram_: scratch space that raises the "high watermark"
  !
  !=====================================================================
  ! ram1:  scratch space allocated in iterative diagonalization 
  !        hpsi, spsi, hr and sr matrices, scalar products
  !        nbnd_l is the estimated dimension of distributed matrices
  !
  nbnd_l = nbndx/np_ortho(1)
  ram1 = complex_size/g_fact * ( 3*nbnd_l**2 ) ! hr,sr,vr/hc,sc,vc 
  IF ( iverbosity > 0 ) WRITE( stdout, 1013 ) 'h,s,v(r/c)', ram1/MB
  add = complex_size/g_fact * ( nkb*npol*nbnd ) ! <psi|beta>, becmod.f90:353
  IF ( iverbosity > 0 ) WRITE( stdout, 1013 ) '<psi|beta>', add/MB
  ram1 = ram1 + add
  !
  ! 
  IF ( isolve == 0 ) THEN
     add = complex_size * nbndx * npol * npwx_l              ! hpsi
     ram1 = ram1 + add
     IF ( iverbosity > 0 ) WRITE( stdout, 1013 ) 'psi', add/MB
     ram1 = ram1 + add
     IF ( iverbosity > 0 ) WRITE( stdout, 1013 ) 'hpsi', add/MB
     IF ( okvan ) THEN
        add = complex_size * nbndx * npol * npwx_l ! spsi
        IF ( iverbosity > 0 ) WRITE( stdout, 1013 ) 'spsi', add/MB 
        ram1 = ram1 + add
     END IF
  END IF
  ram_ = ram1
  !=====================================================================
  !
  !=====================================================================
  ! arrays allocated in approx_screening2 during charge mixing 
  !=====================================================================
  IF ( lscf .AND. imix > 1 ) &
     ram_ = MAX( ram_, complex_size * ngm * 27 + real_size * dffts%nnr )
  !=====================================================================
  !
  !=====================================================================
  ! ram1: scratch space allocated in initialization
  !       wfcinit + rotatewfcgamma can be larger than regterg cegterg if
  !       starting with many atomic wavefunctions
  !       files: wfcinit.f90:246 
  !              rotate_wfc_gamma.f90:45
  !              rotate_wfc_k.f90:55
  !=====================================================================
  IF ( starting_wfc(1:6) == 'atomic' ) THEN
     maxbnd = MAX( natomwfc, nbnd )
  ELSE
     maxbnd = nbnd
  END IF
  ram1 = complex_size * npwx_l * maxbnd  * ( npol + 1) + &
         real_size * maxbnd*(3*maxbnd + 1)
  !
  IF ( iverbosity > 0 ) WRITE( stdout, 1013 ) 'wfcinit/wfcrot', ram1/MB 
  ram_ = MAX( ram_,  ram1 )
  !=====================================================================
  !
  !=====================================================================
  ! ram1: arrays allocated in addusdens_g and addusforce & addusstress
  !       for ultrasoft/paw pp.
  !
  !       N.B: newq is always smaller than addusdens_g
  !
  !       files: addusdens.f90  : 92,93,114
  !              addusforce.f90 : 78,102,105,117,132-133
  !              addusstress.f90
  !=====================================================================
  IF ( okvan ) THEN
     IF ( .not. tqr ) THEN
        ngm_l  = ngm/npool
        maxnab = 0
        maxnij = 0
        DO nt = 1, ntyp
            IF ( upf(nt)%tvanp ) THEN
              !
              ! nij = max number of (ih,jh) pairs per atom type nt
              !
              nij = nh(nt)*(nh(nt)+1)/2
              IF ( nij > maxnij ) maxnij = nij
              !
              ! count max number of atoms of type nt
              !
              nab = 0
              DO na = 1, nat
                  IF ( ityp(na) == nt ) nab = nab + 1
              ENDDO
              IF ( nab > maxnab ) maxnab = nab
            END IF
        END DO
        !                               ylmk0      qmod
        ram1 = real_size * ngm_l * ( lmaxq * lmaxq + 1 )
        !                                aux                   skk   aux2   qgm
        ram1 = ram1 + complex_size * ( ngm*nspin_mag + ngm_l*(maxnab+maxnij+1) )
        !
        IF ( iverbosity > 0 ) WRITE( stdout, 1013 ) 'addusdens', ram1/MB
        ram_ = MAX ( ram_, ram1 )
        !
        ! forces
        !
        IF (lforce) THEN
           !                      vg                       ylmk0     qmod
           ram1 = real_size * (ngm*nspin_mag + ngm_l*( lmaxq*lmaxq + 1 ) )
           !                                    qgm      aux1
           ram1 = ram1 + complex_size * ngm_l * ( maxnij + nat*3 )
           !                           ddeeq
           ram1 = ram1 + real_size * ( maxnij * nat * 3 * nspin_mag )
           IF ( iverbosity > 0 ) WRITE( stdout, 1013 ) 'addusforce', ram1/MB
           !
           ram_ = MAX ( ram_, ram1 )
        END IF
        !
        ! stress
        !
        IF (lstres) THEN
           !                      vg                      ylmk0,dylmk0  qmod
           ram1 = real_size *  (ngm*nspin_mag + ngm_l*( 2*lmaxq*lmaxq + 1 ) )
           !                                    qgm      aux1  aux2
           ram1 = ram1 + complex_size * ngm_l * ( maxnij + 3 + nspin )
           !                           ddeeq
           IF ( iverbosity > 0 ) WRITE( stdout, 1013 ) 'addusstress', ram1/MB
           !
           ram_ = MAX ( ram_, ram1 )
        END IF
     ELSE
        ! nothing allocated in addusdens_r, hurray!
        IF (lforce) THEN
           ram1 = real_size*roughestimate*(nh(nt)*(nh(nt)+1)/2)*3
           IF ( iverbosity > 0 ) WRITE( stdout, 1013 ) 'addusforce_r', ram1/MB
           ram_ = MAX ( ram_, ram1 )
        END IF
     END IF
  END IF
  !
  maxram = ram + ram_
  !
  ! arrays used for global sorting in ggen:
  !    igsrt, g2l, g2sort_g, total dimensions:
  !
  IF ( .NOT. smallmem ) maxram = MAX ( maxram, &
       int_size * 2 * ngm_g + real_size * ngm_g )
  !
  totram = maxram * nproc_image
  IF ( iverbosity > 0 ) THEN
     IF ( ram .lt. GB ) WRITE( stdout, 1010 ) ram/MB, ' MB'
     IF ( ram .ge. GB ) WRITE( stdout, 1010 ) ram/GB, ' GB'
  END IF

  IF ( maxram .lt. GB ) WRITE( stdout, 1011 ) maxram/MB, ' MB'
  IF ( maxram .ge. GB ) WRITE( stdout, 1011 ) maxram/GB, ' GB'
  
  IF ( nproc_image > 1) THEN
     IF ( totram .lt. GB ) WRITE( stdout, 1012 ) totram/MB, ' MB'
     IF ( totram .ge. GB ) WRITE( stdout, 1012 ) totram/GB, ' GB'
  END IF
  !
 1010 format (/5x,'Estimated static dynamical RAM per process > ', F10.2, A3)
 1011 format (/5x,'Estimated max dynamical RAM per process > ', F10.2, A3)
 1012 format (/5x,'Estimated total dynamical RAM > ', F10.2, A3)
 1013 format (/5x,'Dynamical RAM for ', A19, ': ', F10.2, ' MB')
 1014 format (/5x,'Dynamical RAM for ', A19, ': ', '?? MB (not yet implemented')
END subroutine memory_report
