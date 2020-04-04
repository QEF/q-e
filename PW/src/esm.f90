!
! Copyright (C) 2007-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE esm
  !--------------------------------------------------------------------------
  !! This module contains the variables and subroutines needed for the 
  !! EFFECTIVE SCREENING MEDIUM (ESM) METHOD.
  !
  !! Original version by Minoru Otani (AIST), Yoshio Miura (Tohoku U.),
  !! Nicephore Bonet (MIT), Nicola Marzari (MIT), Brandon Wood (LLNL), 
  !! Tadashi Ogitsu (LLNL).  
  !! Constant bias potential (constant-mu) method by Minoru Otani (AIST) and
  !! Nicephore Bonnet (AIST).
  !
  !! Contains subroutines for implementation of:  
  !! 1) ESM (Effective Screening Medium Method) developed by M. Otani and 
  !!    O. Sugino (see PRB 73, 115407 [2006]);  
  !! 2) Constant-mu method developed by N. Bonnet, T. Morishita, O. Sugino, 
  !!    and M. Otani (see PRL 109, 266101 [2012]).
  !
  !! ESM enables description of a surface slab sandwiched between two 
  !! semi-infinite media, making it possible to deal with polarized surfaces 
  !! without using dipole corrections. It is useful for simulating interfaces 
  !! with vacuum, one or more electrodes, or an electrolyte.
  !
  !! Constant-mu scheme with the boundary condition 'bc2' and 'bc3' enables
  !! description of the system is connected to a potentiostat which preserves
  !! the Fermi energy of the system as the target Fermi energy (mu).
  !
  !! Modified subroutines for calculating the Hartree potential, the local 
  !! potential, and the Ewald sum are contained here, along with subroutines for
  !! calculating force contributions based on the modified local potential and 
  !! Ewald term. Constant-mu parts are contained in the fcp.f90.
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE
  !
  PUBLIC :: do_comp_esm, esm_nfit, esm_efield, esm_w, esm_a, esm_bc,   &
            mill_2d, imill_2d, ngm_2d,                                 &
            esm_init, esm_hartree, esm_local, esm_ewald, esm_force_lc, &
            esm_force_ew, esm_printpot, esm_summary
  PUBLIC :: esm_stres_har, esm_stres_ewa, esm_stres_loclong
  !
  !
  LOGICAL :: do_comp_esm=.FALSE.
  !! if TRUE the ESM is active
  INTEGER :: esm_nfit
  !! number of grid points for fit at edges
  REAL(DP) :: esm_efield
  !! field strength
  REAL(DP) :: esm_w
  !! ESM offset from cell edge
  REAL(DP) :: esm_a
  !! smoothness parameter 
  CHARACTER(LEN=3) :: esm_bc
  !! pbc, bc1, ..., bc4
  INTEGER, ALLOCATABLE :: mill_2d(:,:)
  !! miller index on vectors
  INTEGER, ALLOCATABLE :: imill_2d(:,:)
  INTEGER :: ngm_2d = 0
  !! ngm_2d = total number of vectors (h,k) on this proc, excluding 
  !! duplicates with different l values
  !
  REAL(DP), EXTERNAL   :: qe_erf, qe_erfc
  !
  !
  CONTAINS
  !
  !
  !  ... ESM ENERGY AND POTENTIAL SUBROUTINES
  !
  !--------------------------------------------------------------
  SUBROUTINE esm_hartree( rhog, ehart, aux )
    !------------------------------------------------------------
    !! Wrapper for the calculation of the Hartree energy.
    !
    USE kinds,      ONLY : DP    
    USE gvect,      ONLY : ngm    
    USE fft_base,   ONLY : dfftp    
    !
    IMPLICIT NONE
    !
    REAL(DP) :: ehart
    !! Hartree energy    
    COMPLEX(DP) :: rhog(ngm)
    !! n(G)    
    COMPLEX(DP) :: aux(dfftp%nnr)
    !! v_h(G)    
    !
    IF ( esm_bc == 'pbc' ) THEN
       CALL esm_hartree_pbc( rhog, ehart, aux )    
    ELSEIF ( esm_bc == 'bc1' ) THEN
       CALL esm_hartree_bc1( rhog, ehart, aux )    
    ELSEIF ( esm_bc == 'bc2' ) THEN
       CALL esm_hartree_bc2( rhog, ehart, aux )    
    ELSEIF ( esm_bc == 'bc3' ) THEN    
       CALL esm_hartree_bc3( rhog, ehart, aux )    
    ELSEIF ( esm_bc == 'bc4' ) THEN    
       CALL esm_hartree_bc4( rhog, ehart, aux )    
    ENDIF    
    !
  END SUBROUTINE esm_hartree
  !
  !
  !--------------------------------------------------------------
  SUBROUTINE esm_ewaldr( alpha_g, ewr )
    !------------------------------------------------------------
    !! Wrapper for the calculation of the Ewald R-space terms of 
    !! energy.
    !
    USE kinds,    ONLY : DP
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN)  :: alpha_g
    !! alpha term in Ewald sum
    REAL(DP), INTENT(OUT) :: ewr
    !! Ewald R-space term
    !    
    IF ( esm_bc == 'pbc' ) THEN    
       CALL esm_ewaldr_pbc( alpha_g, ewr )    
    ELSEIF ( esm_bc == 'bc1' ) THEN    
       CALL esm_ewaldr_pbc( alpha_g, ewr )    
    ELSEIF ( esm_bc == 'bc2' ) THEN    
       CALL esm_ewaldr_pbc( alpha_g, ewr )    
    ELSEIF ( esm_bc == 'bc3' ) THEN    
       CALL esm_ewaldr_pbc( alpha_g, ewr )    
    ELSEIF ( esm_bc == 'bc4' ) THEN    
       CALL esm_ewaldr_bc4( alpha_g, ewr )    
    ENDIF    
    !
  END SUBROUTINE esm_ewaldr
  !
  !
  !------------------------------------------------------------
  SUBROUTINE esm_ewaldg( alpha_g, ewg )
    !----------------------------------------------------------
    !! Wrapper for the calculation of the Ewald G-space terms of
    !! energy.
    !
    USE kinds,    ONLY : DP
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN)  :: alpha_g
    !! alpha term in Ewald sum
    REAL(DP), INTENT(OUT) :: ewg
    !! Ewald G-space term
    !
    IF ( esm_bc == 'pbc' ) THEN    
       CALL esm_ewaldg_pbc( alpha_g, ewg )    
    ELSEIF ( esm_bc == 'bc1' ) THEN    
       CALL esm_ewaldg_bc1( alpha_g, ewg )    
    ELSEIF ( esm_bc == 'bc2' ) THEN    
       CALL esm_ewaldg_bc2( alpha_g, ewg )    
    ELSEIF ( esm_bc == 'bc3' ) THEN    
       CALL esm_ewaldg_bc3( alpha_g, ewg )    
    ELSEIF ( esm_bc == 'bc4' ) THEN    
       CALL esm_ewaldg_bc4( alpha_g, ewg )    
    ENDIF    
    !
  END SUBROUTINE esm_ewaldg
  !
  !
  !---------------------------------------------------------------
  SUBROUTINE esm_local( aux )
    !-------------------------------------------------------------
    !! Wrapper for the calucaltion of ESM local potential.
    !
    USE kinds,    ONLY : DP    
    USE fft_base, ONLY : dfftp
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(INOUT) :: aux(dfftp%nnr)
    !! The local potential
    !
    IF ( esm_bc == 'pbc' ) THEN    
       CALL esm_local_pbc( aux )    
    ELSEIF ( esm_bc == 'bc1' ) THEN    
       CALL esm_local_bc1( aux )    
    ELSEIF ( esm_bc == 'bc2' ) THEN    
       CALL esm_local_bc2( aux )    
    ELSEIF ( esm_bc == 'bc3' ) THEN    
       CALL esm_local_bc3( aux )    
    ELSEIF ( esm_bc == 'bc4' ) THEN    
       CALL esm_local_bc4( aux )    
    ENDIF
    !
  END SUBROUTINE esm_local
  !
  !
  !
  !  ...  ESM FORCE SUBROUTINES
  !
  !------------------------------------------------------------------
  SUBROUTINE esm_force_ewr( alpha_g, forceion )
    !----------------------------------------------------------------- 
    !! Wrapper for the calculation of the Ewald R-space terms of the 
    !! force.
    !
    USE kinds,      ONLY : DP
    USE ions_base,  ONLY : nat
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN) :: alpha_g
    !! alpha term in Ewald sum
    REAL(DP), INTENT(INOUT) :: forceion(3,nat)
    !! Ewald R-space term of the force
    !
    IF ( esm_bc == 'pbc' ) THEN
       CALL esm_force_ewr_pbc( alpha_g, forceion )
    ELSEIF ( esm_bc == 'bc1' ) THEN
       CALL esm_force_ewr_pbc( alpha_g, forceion )
    ELSEIF ( esm_bc == 'bc2' ) THEN
       CALL esm_force_ewr_pbc( alpha_g, forceion )
    ELSEIF ( esm_bc == 'bc3' ) THEN
       CALL esm_force_ewr_pbc( alpha_g, forceion )
    ELSEIF ( esm_bc == 'bc4' ) THEN
       CALL esm_force_ewr_bc4( alpha_g, forceion )
    ENDIF    
    !
  END SUBROUTINE esm_force_ewr
  !
  !
  !------------------------------------------------------------------
  SUBROUTINE esm_force_ewg( alpha_g, forceion )
    !----------------------------------------------------------------
    !! Wrapper for the calculation of the Ewald G-space terms of the 
    !! force.
    !
    USE kinds,       ONLY : DP
    USE ions_base,   ONLY : nat
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN) :: alpha_g
    !! alpha term in Ewald sum
    REAL(DP), INTENT(OUT) :: forceion(3,nat) 
    !! Ewald G-space term of the force
    !
    IF ( esm_bc == 'pbc' ) THEN
       CALL esm_force_ewg_pbc( alpha_g, forceion )
    ELSEIF ( esm_bc == 'bc1' ) THEN
       CALL esm_force_ewg_bc1( alpha_g, forceion )
    ELSEIF ( esm_bc == 'bc2' ) THEN
       CALL esm_force_ewg_bc2( alpha_g, forceion )
    ELSEIF ( esm_bc == 'bc3' ) THEN
       CALL esm_force_ewg_bc3( alpha_g, forceion )
    ELSEIF ( esm_bc == 'bc4' ) THEN
       CALL esm_force_ewg_bc4( alpha_g, forceion )
    ENDIF
    !
  END SUBROUTINE esm_force_ewg
  !
  !
  !------------------------------------------------------------------
  SUBROUTINE esm_force_lc( aux, forcelc )
    !----------------------------------------------------------------
    !! Wrapper for the calculation of the local ESM force.
    !
    USE kinds,     ONLY : DP
    USE ions_base, ONLY : nat
    USE fft_base,  ONLY : dfftp
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(IN) :: aux(dfftp%nnr)
    !! aux contains n(G) (input)
    REAL(DP), INTENT(INOUT) :: forcelc(3,nat)
    !! ESM local force
    !
    IF ( esm_bc == 'pbc' ) THEN
       CALL esm_force_lc_pbc( aux, forcelc )
    ELSEIF ( esm_bc == 'bc1' ) THEN
       CALL esm_force_lc_bc1( aux, forcelc )
    ELSEIF ( esm_bc == 'bc2' ) THEN
       CALL esm_force_lc_bc2( aux, forcelc )
    ELSEIF ( esm_bc == 'bc3' ) THEN
       CALL esm_force_lc_bc3( aux, forcelc )
    ELSEIF ( esm_bc == 'bc4' ) THEN
       CALL esm_force_lc_bc4( aux, forcelc )
    ENDIF
    !
  END SUBROUTINE esm_force_lc
  !
  !
  !
  !  ...  ESM STRESS SUBROUTINES
  !
  !------------------------------------------------------------------
  SUBROUTINE esm_stres_har( sigmahar, rhog )
    !----------------------------------------------------------------
    !! Wrapper for the calculation of the Hartree term of the stress.
    !
    USE kinds,    ONLY : DP
    USE gvect,    ONLY : ngm
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(OUT) :: sigmahar(3,3)
    !! Hartree term of the stress
    COMPLEX(DP), INTENT(IN) :: rhog(ngm)
    !! n(G)
    !
    SELECT CASE( esm_bc )
    CASE( 'pbc' )
       STOP 'esm_stres_har must not be called for esm_bc = pbc'
    CASE( 'bc1' )
       CALL esm_stres_har_bc1( sigmahar, rhog )
    CASE( 'bc2' )
       CALL esm_stres_har_bc2( sigmahar, rhog )
    CASE( 'bc3' )
       CALL esm_stres_har_bc3( sigmahar, rhog )
    CASE( 'bc4' )
       STOP 'esm_stres_har has not yet implemented for esm_bc = bc4'
    END SELECT
    !
    RETURN
    !
  END SUBROUTINE esm_stres_har
  !
  !
  !------------------------------------------------------------------------
  SUBROUTINE esm_stres_ewa( sigmaewa )
    !-----------------------------------------------------------------------
    !! Calculates Ewald stresswith both G- and R-space terms. 
    !! Determines optimal alpha. Should hopefully work for any structure.
    !
    USE kinds,     ONLY : DP
    USE constants, ONLY : tpi
    USE cell_base, ONLY : tpiba2
    USE ions_base, ONLY : zv, nat, ityp
    USE gvect,     ONLY : gcutm
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(OUT) :: sigmaewa(3,3)
    !! output: the Ewald stress
    !
    ! ... local variables
    !
    INTEGER :: ia
    ! counter on atoms
    REAL(DP) :: charge, alpha, upperbound
    ! total ionic charge in the cell
    ! alpha term in ewald sum
    ! the maximum radius to consider real space sum
    REAL(DP) :: sigmaewg(3,3), sigmaewr(3,3)
    ! ewald stress computed in reciprocal space
    ! ewald stress computed in real space
    !
    charge = SUM(zv(ityp(:)))
    !
    ! choose alpha in order to have convergence in the sum over G
    ! upperbound is a safe upper bound for the error in the sum over G
    alpha = 2.9d0
    !
    DO
       alpha = alpha - 0.1d0
       IF (alpha <= 0.d0) CALL errore( 'esm_stres_ewa', 'optimal alpha not found', 1 )
       upperbound = 2.d0 * charge**2 * SQRT(2.d0 * alpha / tpi) * &
            qe_erfc ( SQRT(tpiba2 * gcutm / 4.d0 / alpha) )
       IF ( upperbound < 1.0d-7 ) EXIT
    ENDDO
    !
    ! G-space sum here.
    ! Determine IF this processor contains G=0 and set the constant term
    CALL esm_stres_ewg( alpha, sigmaewg )
    !
    ! R-space sum here (only for the processor that contains G=0)
    CALL esm_stres_ewr( alpha, sigmaewr )
    !
    sigmaewa(:,:) = sigmaewg(:,:) + sigmaewr(:,:)
    !
    RETURN
    !
  END SUBROUTINE esm_stres_ewa
  !
  !
  !------------------------------------------------------------------------------
  SUBROUTINE esm_stres_loclong( sigmaloclong, rhog )
    !---------------------------------------------------------------------------
    !! Wrapper for the calculation of the local term of the stress.
    !
    USE kinds,    ONLY : DP
    USE gvect,    ONLY : ngm
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(OUT)   :: sigmaloclong(3,3)
    !! the local term of the stress
    COMPLEX(DP), INTENT(IN) :: rhog(ngm)
    !! n(G)
    !
    SELECT CASE( esm_bc )
    CASE( 'pbc' )
       STOP 'esm_stres_loclong must not be called for esm_bc = pbc'
    CASE( 'bc1' )
       CALL esm_stres_loclong_bc1( sigmaloclong, rhog )
    CASE( 'bc2' )
       CALL esm_stres_loclong_bc2( sigmaloclong, rhog )
    CASE( 'bc3' )
       CALL esm_stres_loclong_bc3( sigmaloclong, rhog )
    CASE( 'bc4' )
       STOP 'esm_stres_loclong has not yet implemented for esm_bc = bc4'
    END SELECT
    !
    RETURN
    !
  END SUBROUTINE esm_stres_loclong
  !
  !
  !------------------------------------------------------------------------
  SUBROUTINE esm_rgen_2d( dtau, rmax, mxr, at, bg, r, r2, nrm )
    !-----------------------------------------------------------------------
    !! Generates neighbours shells (cartesian, in units of lattice parameter)
    !! with length < rmax, and returns them in order of increasing length.
    !
    USE kinds, ONLY : DP
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: mxr
    !! maximum number of vectors
    INTEGER, INTENT(OUT) :: nrm
    !! the number of vectors with r^2 < rmax^2
    REAL(DP), INTENT(IN) :: at(3,3)
    !! lattice vectors ( a1=at(:,1), a2=at(:,2), a3=at(:,3) )
    REAL(DP), INTENT(IN) :: bg(3,3)
    !! reciprocal lattice vectors ( b1=bg(:,1), b2=bg(:,2), b3=bg(:,3) )
    REAL(DP), INTENT(IN) :: dtau(3)
    !! increase of vector r
    REAL(DP), INTENT(IN) :: rmax
    !! length of neighbours shells < rmax
    REAL(DP), INTENT(OUT) :: r(3,mxr)
    !! r(:) = i*a1(:) + j*a2(:) + k*a3(:) - dtau(:)
    REAL(DP), INTENT(OUT) :: r2(mxr)
    !! r2 = r^2
    !
    ! ... local variables
    !
    INTEGER, ALLOCATABLE :: irr(:)
    INTEGER  ::  nm1, nm2, i, j, ipol, ir, indsw, iswap
    REAL(DP) :: ds(3), dtau0(3)
    REAL(DP) :: t(3), tt, swap
    REAL(DP), EXTERNAL :: dnrm2
    !
    !
    nrm = 0
    IF (rmax==0.d0) RETURN
    !
    ! bring dtau into the unit cell centered on the origin - prevents trouble
    ! IF atomic positions are not centered around the origin but displaced
    ! far away (remember that translational invariance allows this!)
    !
    ds(:) = MATMUL( dtau(:), bg(:,:) )
    ds(:) = ds(:) - ANINT(ds(:))
    dtau0(:) = MATMUL( at(:,:), ds(:) )
    !
    ALLOCATE( irr(mxr) )
    !
    ! these are estimates of the maximum values of needed integer indices
    !
    nm1 = INT(dnrm2(3,bg(1,1),1) * rmax) + 2
    nm2 = INT(dnrm2(3,bg(1,2),1) * rmax) + 2
    !
    DO i = -nm1, nm1
       DO j = -nm2, nm2
          tt = 0.d0
          DO ipol = 1, 3
             t(ipol) = i*at(ipol,1) + j*at(ipol,2) - dtau0(ipol)
             tt = tt + t(ipol) * t(ipol)
          ENDDO
          IF (tt<=rmax**2 .AND. ABS(tt) >1.d-10) THEN
             nrm = nrm + 1
             IF (nrm>mxr) CALL errore( 'esm_rgen_2d', 'too many r-vectors', nrm )
             DO ipol = 1, 3
                r(ipol,nrm) = t(ipol)
             ENDDO
             r2(nrm) = tt
          ENDIF
       ENDDO
    ENDDO
    !
    !   reorder the vectors in order of increasing magnitude
    !
    !   initialize the index inside sorting routine
    !
    irr (1) = 0
    IF (nrm>1) CALL hpsort (nrm, r2, irr)
    DO ir = 1, nrm - 1
20     indsw = irr (ir)
       IF (indsw/=ir) THEN
          DO ipol = 1, 3
             swap = r(ipol,indsw)
             r(ipol,indsw) = r(ipol,irr(indsw))
             r(ipol,irr(indsw)) = swap
          ENDDO
          iswap = irr(ir)
          irr(ir) = irr(indsw)
          irr(indsw) = iswap
          GOTO 20
       ENDIF
    ENDDO
    DEALLOCATE( irr )
    !
    RETURN
    !
  END SUBROUTINE esm_rgen_2d
  !
  !
  !-----------------------------------------------------------
  SUBROUTINE esm_init()
    !--------------------------------------------------------
    !! Wrapper to ESM initialization.
    !
    USE fft_base, ONLY : dfftp
    !
    IMPLICIT NONE
    !
    CALL esm_ggen_2d()
    !
  END SUBROUTINE esm_init
  !
  !
  !------------------------------------------------------------
  SUBROUTINE esm_ggen_2d()
    !-----------------------------------------------------------
    !! ESM initialization routine.
    !
    USE fft_base,         ONLY : dfftp
    USE gvect,            ONLY : ngm, mill
    !
    IMPLICIT NONE
    !
    INTEGER :: n1xh, n2xh, ng, n1, n2, ng_2d
    LOGICAL, ALLOCATABLE  :: do_mill_2d(:,:)
    !
    ! Make g parallel array
    !
    n1xh = dfftp%nr1x/2
    n2xh = dfftp%nr2x/2
    ALLOCATE( do_mill_2d(-n1xh:n1xh,-n2xh:n2xh) )
    do_mill_2d(:,:) = .FALSE.
    !
    DO ng = 1, ngm
       n1 = mill(1,ng)
       n2 = mill(2,ng)
       do_mill_2d(n1,n2) = .TRUE.
    ENDDO
    !
    ngm_2d = COUNT( do_mill_2d )
    !
    !*** do_mill_2d(h,k) = .true. means there is an h,k vector on this proc
    !*** ngm_2d = total number of vectors (h,k) on this proc, excluding 
    !*** duplicates with different l values
    !
    IF( .NOT. ALLOCATED(mill_2d ) ) ALLOCATE( mill_2d(2,ngm_2d) )
    IF( .NOT. ALLOCATED(imill_2d) ) ALLOCATE( imill_2d(-n1xh:n1xh,-n2xh:n2xh) )
    !
    mill_2d(:,:) = 0
    imill_2d(:,:) = 0
    ng_2d = 1
    !
    DO n1 = -n1xh, n1xh
      DO n2 = -n2xh, n2xh
        IF( do_mill_2d(n1,n2) ) THEN
           mill_2d(1,ng_2d) = n1
           mill_2d(2,ng_2d) = n2
           imill_2d(n1,n2) = ng_2d
           ng_2d = ng_2d + 1
        ENDIF
      ENDDO
    ENDDO
    !
    DEALLOCATE( do_mill_2d )
    !
    !**** mill_2d(:,ig) = h,k indices of vector ig
    !**** imill_2d(h,k) = 2d index of vector with h,k indices
    !**** ng_2d = total number of 2d g vectors on this proc
    !
    RETURN
    !
  END SUBROUTINE esm_ggen_2d
  !
  !
  !
  !  ...  ESM HARTREE SUBROUTINE
  !
  !-----------------------------------------------------------------------------------
  SUBROUTINE esm_hartree_pbc( rhog, ehart, aux )
    !-------------------------------------------------------------------------------
    !! Error routine: esm_hartree must not be called for esm_bc = pbc.
    !
    USE gvect,      ONLY : ngm
    USE fft_base,   ONLY : dfftp
    !
    IMPLICIT NONE
    !
    REAL(DP) :: ehart
    !! Hartree energy
    COMPLEX(DP) :: rhog(ngm)
    !! n(G)
    COMPLEX(DP) :: aux(dfftp%nnr)
    !! v_h(G)
    !
    STOP 'esm_hartree must not be called for esm_bc = pbc'
    !
  END SUBROUTINE esm_hartree_pbc
  !
  !
  !-------------------------------------------------------------------------------
  SUBROUTINE esm_hartree_bc1( rhog, ehart, aux )
    !-----------------------------------------------------------------------------
    !! Calculation of Hartree energy and potential - bc1.
    !
    USE constants,        ONLY : tpi, fpi, e2
    USE gvect,            ONLY : ngm, mill
    USE cell_base,        ONLY : omega, alat, tpiba2, at, bg
    USE control_flags,    ONLY : gamma_only
    USE mp_bands,         ONLY : intra_bgrp_comm
    USE mp,               ONLY : mp_sum
    USE fft_base,         ONLY : dfftp
    USE fft_scalar,       ONLY : cft_1z
    !
    IMPLICIT NONE
    !
    REAL(DP) :: ehart
    !! Hartree energy
    COMPLEX(DP) :: rhog(ngm)
    !! n(G)
    COMPLEX(DP) :: aux(dfftp%nnr)
    !! v_h(G)   
    !
    ! ... local variables
    !
    INTEGER :: k1, k2, k3, iz, ng, n1, n2, n3, nz_r, nz_l, ng_2d
    REAL(DP) :: t(2), z, z0, gp, gp2, kn, cc0, ss0, L, z_l, z_r, eh,  &
                arg1, arg2
    COMPLEX(DP) :: ci, f1, f2, f3, f4, a0, a1, a2, a3, c_r, c_l, s_r, &
                   s_l, rg3, tmp1, tmp2, tmp3
    COMPLEX(DP), ALLOCATABLE  :: rhog3(:,:), vg(:), vg_r(:), vg3(:,:)
    !
    ALLOCATE( rhog3(dfftp%nr3,ngm_2d), vg3(dfftp%nr3,ngm_2d) )
    !
    ! Map to FFT mesh (dfftp%nr3,ngm_2d)
    rhog3(:,:) = (0.d0,0.d0)
    !
    DO ng = 1, ngm
       n1 = mill(1,ng)
       n2 = mill(2,ng)
       ng_2d = imill_2d(n1,n2)
       n3 = mill(3,ng)+1
       IF (n3<1) n3 = n3 + dfftp%nr3
       rg3 = rhog(ng)
       rhog3(n3,ng_2d)=rg3
       IF ( gamma_only .AND. n1==0 .AND. n2==0 ) THEN
          n3 = -mill(3,ng)+1
          IF (n3<1) n3 = n3 + dfftp%nr3
          rhog3(n3,ng_2d) = CONJG(rg3)
       ENDIF
    ENDDO
    ! End mapping
    !
    vg3(:,:) = (0.d0,0.d0)
    L = at(3,3)*alat
    z0 = L/2.d0
    ci = (0.d0,1.d0)
    !
    !****For gp!=0 case ********************
    ALLOCATE( vg(dfftp%nr3), vg_r(dfftp%nr3) )
    !
    DO ng_2d = 1, ngm_2d
       k1 = mill_2d(1,ng_2d)
       k2 = mill_2d(2,ng_2d)
       IF (k1==0 .AND. k2==0) CYCLE
       t(1:2) = k1 * bg(1:2,1) + k2 * bg(1:2,2)
       gp2 = SUM(t(:) * t(:)) * tpiba2
       gp = SQRT(gp2)
       tmp1 = (0.d0,0.d0); tmp2 = (0.d0,0.d0);
       vg(:) = (0.d0,0.d0)
       !
       DO iz=1, dfftp%nr3
          IF (iz<=dfftp%nr3/2) kn = DBLE(iz-1)*tpi/L
          IF (iz> dfftp%nr3/2) kn = DBLE(iz-1-dfftp%nr3)*tpi/L
          cc0 = COS(kn*z0)
          ss0 = SIN(kn*z0)
          rg3 = rhog3(iz,ng_2d)
          ! bc1
          vg(iz) = fpi*rg3/(gp**2+kn**2)
          tmp1 = tmp1+rg3*(cc0+ci*ss0)/(gp-ci*kn)
          tmp2 = tmp2+rg3*(cc0-ci*ss0)/(gp+ci*kn)
       ENDDO
       !
       vg3(:,ng_2d)=vg(:)
       !
       ! real part
       vg_r(:) = (0.d0,0.d0)
       !
       DO iz = 1, dfftp%nr3
          k3 = iz-1
          IF (k3 > dfftp%nr3/2) k3 = iz-dfftp%nr3-1
          z = DBLE(k3)/DBLE(dfftp%nr3)*L
          ! bc1
          arg1 =  gp*(z-z0)
          arg2 = -gp*(z+z0)
          vg_r(iz) = -tpi/gp*(EXP(arg1)*tmp1+EXP(arg2)*tmp2)
       ENDDO
       !
       CALL cft_1z( vg_r, 1, dfftp%nr3, dfftp%nr3, -1, vg )
       !
       vg3(:,ng_2d) = (vg3(:,ng_2d)+vg(:))*e2 ! factor e2: hartree -> Ry.
       !
    ENDDO
    !
    DEALLOCATE( vg, vg_r )
    !
    !****For gp=0 case ****
    ng_2d = imill_2d(0,0)
    !
    IF ( ng_2d > 0 ) THEN
       !
       ALLOCATE( vg(dfftp%nr3), vg_r(dfftp%nr3) )
       tmp1=(0.d0,0.d0); tmp2=(0.d0,0.d0); tmp3=(0.d0,0.d0);
       vg(:)=(0.d0,0.d0);
       rg3 = rhog3(1,ng_2d)
       vg(1) = -tpi*z0**2*rg3
       !
       DO iz = 2, dfftp%nr3
          IF (iz<=dfftp%nr3/2) kn = DBLE(iz-1)*tpi/L
          IF (iz> dfftp%nr3/2) kn = DBLE(iz-1-dfftp%nr3)*tpi/L
          rg3 = rhog3(iz,ng_2d)
          cc0 = COS(kn*z0)
          ss0 = SIN(kn*z0)
          tmp1 = tmp1 + rg3*ci*(cc0+ci*ss0)/kn
          tmp2 = tmp2 + rg3*ci*(cc0-ci*ss0)/kn
          tmp3 = tmp3 + rg3*cc0/kn**2
          vg(iz) = fpi*rg3/(kn**2)
       ENDDO
       !
       vg3(:,ng_2d) = vg(:)
       !
       ! real part
       vg_r(:) = (0.d0,0.d0)
       rg3=rhog3(1,ng_2d)
       DO iz = 1, dfftp%nr3
          k3 = iz-1
          IF (k3 > dfftp%nr3/2) k3 = iz-dfftp%nr3-1
          z = DBLE(k3)/DBLE(dfftp%nr3)*L
          ! bc1
          vg_r(iz) = -tpi*z**2*rg3    &
                     -tpi*(z-z0)*tmp1 &
                     -tpi*(z+z0)*tmp2 &
                     -fpi*tmp3
       ENDDO
       !
       ! start smoothing
       nz_l = dfftp%nr3/2+1+esm_nfit
       nz_r = dfftp%nr3/2+1-esm_nfit
       !
       z_l = DBLE(nz_l-1)*L/DBLE(dfftp%nr3)-L
       z_r = DBLE(nz_r-1)*L/DBLE(dfftp%nr3)
       !
       f1 = -tpi*z_r**2*rg3    &
            -tpi*(z_r-z0)*tmp1 &
            -tpi*(z_r+z0)*tmp2 &
            -fpi*tmp3
       f2 = -tpi*z_l**2*rg3    &
            -tpi*(z_l-z0)*tmp1 &
            -tpi*(z_l+z0)*tmp2 &
            -fpi*tmp3
       f3 = -fpi*z_r*rg3 &
            -tpi*tmp1    &
            -tpi*tmp2
       f4 = -fpi*z_l*rg3 &
            -tpi*tmp1    &
            -tpi*tmp2
       !
       z_r = z_r
       z_l = z_l + L
       !
       a0 = (f1*z_l**2*(z_l-3.d0*z_r)+z_r*(f3*z_l**2*(-z_l+z_r)    &
            +z_r*(f2*(3.d0*z_l-z_r)+f4*z_l*(-z_l+z_r))))/(z_l-z_r)**3
       a1 = (f3*z_l**3+z_l*(6.d0*f1-6.d0*f2+(f3+2.d0*f4)*z_l)*z_r  &
            -(2*f3+f4)*z_l*z_r**2-f4*z_r**3)/(z_l-z_r)**3
       a2 = (-3*f1*(z_l+z_r)+3.d0*f2*(z_l+z_r)-(z_l-z_r)*(2*f3*z_l &
            +f4*z_l+f3*z_r+2*f4*z_r))/(z_l-z_r)**3
       a3 = (2.d0*f1-2.d0*f2+(f3+f4)*(z_l-z_r))/(z_l-z_r)**3
       !
       DO iz = nz_r,nz_l
          z = DBLE(iz-1)/DBLE(dfftp%nr3)*L
          vg_r(iz) = (a0+a1*z+a2*z**2+a3*z**3)
       ENDDO
       ! end smoothing
       !
       CALL cft_1z( vg_r, 1, dfftp%nr3, dfftp%nr3, -1, vg )
       vg3(:,ng_2d) = (vg3(:,ng_2d)+vg(:))*e2 ! factor e2: hartree -> Ry.
       !
       DEALLOCATE( vg, vg_r )
       !
    ENDIF ! if( ng_2d > 0 )
    !
    ! Hartree Energy
    ehart = 0.d0
    eh = 0d0
    !
    DO ng_2d = 1, ngm_2d
       k1 = mill_2d(1,ng_2d)
       k2 = mill_2d(2,ng_2d)
       eh = eh + SUM( vg3(:,ng_2d)*CONJG(rhog3(:,ng_2d)) )
    ENDDO
    !
    ehart = ehart+eh
    IF ( gamma_only ) THEN
       ehart = ehart * 2d0
       ng_2d = imill_2d(0,0)
       IF ( ng_2d > 0 ) THEN
          ehart = ehart - SUM( vg3(:,ng_2d)*CONJG(rhog3(:,ng_2d)) )
       ENDIF
    ENDIF
    ehart = ehart*omega*0.5d0
    !
    CALL mp_sum( ehart, intra_bgrp_comm )
    !
    ! Map to FFT mesh (dfftp%nrx)
    aux = 0.0d0
    !
    DO ng = 1, ngm
       n1 = mill(1,ng)
       n2 = mill(2,ng)
       ng_2d = imill_2d(n1,n2)
       n3 = mill(3,ng) + 1
       IF (n3<1) n3 = n3 + dfftp%nr3
       aux(dfftp%nl(ng)) = aux(dfftp%nl(ng)) + vg3(n3,ng_2d)
       IF (gamma_only) THEN
          aux(dfftp%nlm(ng)) = CONJG(aux(dfftp%nl(ng)))
       ENDIF
    ENDDO
    !
    DEALLOCATE( rhog3, vg3 )
    !
    RETURN
    !
  END SUBROUTINE esm_hartree_bc1
  !
  !
  !--------------------------------------------------------------------
  SUBROUTINE esm_hartree_bc2( rhog, ehart, aux )
    !------------------------------------------------------------------
    !! Calculation of Hartree energy and potential - bc2
    !
    USE constants,        ONLY : tpi, fpi, e2
    USE gvect,            ONLY : ngm, mill
    USE cell_base,        ONLY : omega, alat, tpiba2, at, bg
    USE control_flags,    ONLY : gamma_only
    USE mp_bands,         ONLY : intra_bgrp_comm
    USE mp,               ONLY : mp_sum
    USE fft_base,         ONLY : dfftp
    USE fft_scalar,       ONLY : cft_1z
    !
    IMPLICIT NONE
    !
    REAL(DP) :: ehart
    !! Hartree energy
    COMPLEX(DP) :: rhog(ngm)
    !! n(G)
    COMPLEX(DP) :: aux(dfftp%nnr)
    !! v_h(G)
    !
    ! ... local variables
    !
    INTEGER :: k1, k2, k3, iz, ng, n1, n2, n3, nz_r, nz_l, ng_2d
    REAL(DP):: t(2), z, z0, z1, gp, gp2, kn, cc0, ss0, L, &
               z_l, z_r, eh, arg1, arg2, arg3, arg4, arg5
    COMPLEX(DP) :: ci, f1, f2, f3, f4, a0, a1, a2, a3, c_r, c_l, &
                   s_r, s_l, rg3, tmp, tmp1, tmp2, tmp3, tmp4
    COMPLEX(DP), ALLOCATABLE :: rhog3(:,:), vg(:), vg_r(:), vg3(:,:)
    !
    ALLOCATE( rhog3(dfftp%nr3,ngm_2d), vg3(dfftp%nr3,ngm_2d) )
    !
    ! Map to FFT mesh (dfftp%nr3,ngm_2d)
    rhog3(:,:) = (0.d0,0.d0)
    !
    DO ng = 1, ngm
       n1 = mill(1,ng)
       n2 = mill(2,ng)
       ng_2d = imill_2d(n1,n2)
       n3 = mill(3,ng)+1
       IF (n3<1) n3 = n3 + dfftp%nr3    
       rg3 = rhog(ng)
       rhog3(n3,ng_2d) = rg3
       IF ( gamma_only .AND. n1==0 .AND. n2==0 ) THEN
          n3 = -mill(3,ng)+1
          IF (n3<1) n3 = n3 + dfftp%nr3
          rhog3(n3,ng_2d) = CONJG(rg3)
       ENDIF
    ENDDO
    ! End mapping
    !
    vg3(:,:) = (0.d0,0.d0)
    L  = at(3,3)*alat
    z0 = L/2.d0
    z1 = z0 + esm_w
    ci = (0.d0,1.d0)
    !
    !****For gp!=0 case ****
    ALLOCATE( vg(dfftp%nr3), vg_r(dfftp%nr3) )
    !
    DO ng_2d = 1, ngm_2d
       k1 = mill_2d(1,ng_2d)
       k2 = mill_2d(2,ng_2d)
       IF (k1==0 .AND. k2==0) CYCLE
       t(1:2) = k1 * bg(1:2,1) + k2 * bg(1:2,2)
       gp2 = SUM(t(:) * t(:)) * tpiba2
       gp = SQRT(gp2)
       tmp1 = (0.d0,0.d0); tmp2 = (0.d0,0.d0);
       vg(:) = (0.d0,0.d0)
       !
       DO iz = 1, dfftp%nr3
          IF (iz<=dfftp%nr3/2) kn = DBLE(iz-1)*tpi/L
          IF (iz> dfftp%nr3/2) kn = DBLE(iz-1-dfftp%nr3)*tpi/L
          cc0 = COS(kn*z0)
          ss0 = SIN(kn*z0)
          rg3 = rhog3(iz,ng_2d)
          ! bc2
          arg1 =  gp*(z1-z0)
          arg2 = -gp*(z1-z0)
          vg(iz) = fpi*rg3/(gp**2+kn**2)
          tmp  = ((gp+ci*kn)*EXP(arg1)+(gp-ci*kn)*EXP(arg2))/(2.d0*gp)
          tmp1 = tmp1+rg3*(cc0+ci*ss0)/(gp**2+kn**2)*tmp
          tmp  = ((gp-ci*kn)*EXP(arg1)+(gp+ci*kn)*EXP(arg2))/(2.d0*gp)
          tmp2 = tmp2+rg3*(cc0-ci*ss0)/(gp**2+kn**2)*tmp
       ENDDO
       !
       vg3(:,ng_2d) = vg(:)
       !
       ! real part
       vg_r(:) = (0.d0,0.d0)
       DO iz = 1, dfftp%nr3
          k3 = iz-1
          IF (k3 > dfftp%nr3/2) k3 = iz-dfftp%nr3-1
          z = DBLE(k3)/DBLE(dfftp%nr3)*L
          ! bc2
          arg1 =  gp*(z-z1)
          arg2 = -gp*(z+z1)
          arg3 =  gp*(z-3.d0*z1)
          arg4 = -gp*(z+3.d0*z1)
          arg5 = -4.d0*gp*z1
          vg_r(iz) = - fpi*(EXP(arg1)-EXP(arg4))*tmp1/(1.d0-EXP(arg5)) &
                     + fpi*(EXP(arg3)-EXP(arg2))*tmp2/(1.d0-EXP(arg5))
       ENDDO
       !
       CALL cft_1z( vg_r, 1, dfftp%nr3, dfftp%nr3, -1, vg )
       !
       vg3(:,ng_2d) = (vg3(:,ng_2d)+vg(:))*e2 ! factor e2: hartree -> Ry.
       !
    ENDDO
    !
    DEALLOCATE( vg, vg_r )
    !
    !****For gp=0 case ****
    ng_2d = imill_2d(0,0)
    !
    IF ( ng_2d > 0 ) THEN
       !
       ALLOCATE( vg(dfftp%nr3), vg_r(dfftp%nr3) )
       tmp1=(0.d0,0.d0); tmp2=(0.d0,0.d0)
       tmp3=(0.d0,0.d0); tmp4=(0.d0,0.d0)
       vg(:) = (0.d0,0.d0)
       rg3 = rhog3(1,ng_2d)
       vg(1) = tpi*(2.d0*z1-z0)*z0*rg3
       !
       DO iz = 2, dfftp%nr3
          IF (iz<=dfftp%nr3/2) kn = DBLE(iz-1)*tpi/L
          IF (iz> dfftp%nr3/2) kn = DBLE(iz-1-dfftp%nr3)*tpi/L
          rg3 = rhog3(iz,ng_2d)
          cc0 = COS(kn*z0)
          ss0 = SIN(kn*z0)
          tmp1 = tmp1+rg3*(cc0+ci*ss0)/kn**2
          tmp2 = tmp2+rg3*(cc0-ci*ss0)/kn**2
          tmp3 = tmp3+rg3*ci*cc0/kn
          tmp4 = tmp4+rg3*ss0/kn
          vg(iz) = fpi*rg3/(kn**2)
       ENDDO
       !
       vg3(:,ng_2d)=vg(:)
       !  
       ! real part
       vg_r(:) = (0.d0,0.d0)
       rg3 = rhog3(1,ng_2d)
       !
       DO iz = 1, dfftp%nr3
          k3 = iz-1
          IF (k3 > dfftp%nr3/2) k3 = iz-dfftp%nr3-1
          z = DBLE(k3)/DBLE(dfftp%nr3)*L
          vg_r(iz) = -tpi*z**2*rg3 &
                     -tpi*(z+z1)*tmp1/z1 &
                     +tpi*(z-z1)*tmp2/z1 &
                     -fpi*z*(z1-z0)/z1*tmp3 &
                     +fpi*(z1-z0)*tmp4
       ENDDO
       !
       !  start smoothing
       nz_l = dfftp%nr3/2+1+esm_nfit
       nz_r = dfftp%nr3/2+1-esm_nfit
       !
       z_l = DBLE(nz_l-1)*L/DBLE(dfftp%nr3)-L
       z_r = DBLE(nz_r-1)*L/DBLE(dfftp%nr3)
       !
       f1 = -tpi*z_r**2*rg3 &
            -tpi*(z_r+z1)*tmp1/z1 &
            +tpi*(z_r-z1)*tmp2/z1 &
            -fpi*z_r*(z1-z0)/z1*tmp3 &
            +fpi*(z1-z0)*tmp4
       f2 = -tpi*z_l**2*rg3 &
            -tpi*(z_l+z1)*tmp1/z1 &
            +tpi*(z_l-z1)*tmp2/z1 &
            -fpi*z_l*(z1-z0)/z1*tmp3 &
            +fpi*(z1-z0)*tmp4
       f3 = -fpi*z_r*rg3 &
            -tpi*tmp1/z1 &
            +tpi*tmp2/z1 &
            -fpi*(z1-z0)/z1*tmp3
       f4 = -fpi*z_l*rg3 &
            -tpi*tmp1/z1 &
            +tpi*tmp2/z1 &
            -fpi*(z1-z0)/z1*tmp3
       !
       z_r = z_r
       z_l = z_l+L
       !
       a0 = (f1*z_l**2*(z_l-3.d0*z_r)+z_r*(f3*z_l**2*(-z_l+z_r)    &
            +z_r*(f2*(3.d0*z_l-z_r)+f4*z_l*(-z_l+z_r))))/(z_l-z_r)**3
       a1 = (f3*z_l**3+z_l*(6.d0*f1-6.d0*f2+(f3+2.d0*f4)*z_l)*z_r  &
            -(2*f3+f4)*z_l*z_r**2-f4*z_r**3)/(z_l-z_r)**3
       a2 = (-3*f1*(z_l+z_r)+3.d0*f2*(z_l+z_r)-(z_l-z_r)*(2*f3*z_l &
            +f4*z_l+f3*z_r+2*f4*z_r))/(z_l-z_r)**3
       a3 = (2.d0*f1-2.d0*f2+(f3+f4)*(z_l-z_r))/(z_l-z_r)**3
       !
       DO iz = nz_r, nz_l
          z = DBLE(iz-1)/DBLE(dfftp%nr3)*L
          vg_r(iz) = (a0+a1*z+a2*z**2+a3*z**3)
       ENDDO
       ! end smoothing
       !
       CALL cft_1z( vg_r, 1, dfftp%nr3, dfftp%nr3, -1, vg )
       vg3(:,ng_2d) = (vg3(:,ng_2d)+vg(:))*e2 ! factor e2: hartree -> Ry.
       !
       DEALLOCATE( vg, vg_r )
       !
    ENDIF ! if( ng_2d > 0 )
    !
    ! Hartree Energy
    ehart = 0.d0
    eh = 0d0
    !
    DO ng_2d = 1, ngm_2d
       k1 = mill_2d(1,ng_2d)
       k2 = mill_2d(2,ng_2d)
       eh = eh + SUM( vg3(:,ng_2d)*CONJG(rhog3(:,ng_2d)) )
    ENDDO
    !
    ehart = ehart + eh
    !
    IF ( gamma_only ) THEN
       ehart = ehart * 2d0
       ng_2d = imill_2d(0,0)
       IF ( ng_2d > 0 ) THEN
          ehart = ehart - SUM( vg3(:,ng_2d)*CONJG(rhog3(:,ng_2d)) )
       ENDIF
    ENDIF
    !
    ehart = ehart*omega*0.5d0
    !
    CALL mp_sum( ehart, intra_bgrp_comm )
    !
    ! Map to FFT mesh (dfftp%nrx)
    aux = 0.0d0
    !
    DO ng = 1, ngm
       n1 = mill(1,ng)
       n2 = mill(2,ng)
       ng_2d = imill_2d(n1,n2)
       n3 = mill(3,ng) + 1
       IF (n3<1) n3 = n3 + dfftp%nr3
       aux(dfftp%nl(ng)) = aux(dfftp%nl(ng)) + vg3(n3,ng_2d)
       IF (gamma_only) THEN
          aux(dfftp%nlm(ng)) = CONJG(aux(dfftp%nl(ng)))
       ENDIF
    ENDDO
    !
    DEALLOCATE( rhog3, vg3 )
    !
    RETURN
    !
  END SUBROUTINE esm_hartree_bc2
  !
  !
  !----------------------------------------------------------------
  SUBROUTINE esm_hartree_bc3( rhog, ehart, aux )
    !---------------------------------------------------------------
    !! Calculation of Hartree energy and potential - bc3
    !
    USE constants,        ONLY : tpi, fpi, e2
    USE gvect,            ONLY : ngm, mill
    USE cell_base,        ONLY : omega, alat, tpiba2, at, bg
    USE control_flags,    ONLY : gamma_only
    USE mp_bands,         ONLY : intra_bgrp_comm
    USE mp,               ONLY : mp_sum
    USE fft_base,         ONLY : dfftp
    USE fft_scalar,       ONLY : cft_1z
    !
    IMPLICIT NONE
    !
    REAL(DP) :: ehart
    !! Hartree energy
    COMPLEX(DP) :: rhog(ngm)
    !! n(G)
    COMPLEX(DP) :: aux(dfftp%nnr)
    !! v_h(G)
    !
    ! ... local variables
    !
    INTEGER :: k1, k2, k3, iz, ng, n1, n2, n3, nz_r, nz_l, ng_2d
    REAL(DP) :: t(2), z, z0, z1, gp, gp2, kn, cc0, ss0, L, z_l,  &
                z_r, eh, arg1, arg2, arg3
    COMPLEX(DP) :: ci, f1, f2, f3, f4, a0, a1, a2, a3, c_r, c_l, &
                   s_r, s_l, rg3, tmp, tmp1, tmp2, tmp3
    COMPLEX(DP), ALLOCATABLE :: rhog3(:,:), vg(:), vg_r(:), vg3(:,:)
    !
    ALLOCATE( rhog3(dfftp%nr3,ngm_2d), vg3(dfftp%nr3,ngm_2d) )
    !
    ! Map to FFT mesh (dfftp%nr3,ngm_2d)
    rhog3(:,:)=(0.d0,0.d0)
    !
    DO ng = 1, ngm
       n1 = mill(1,ng)
       n2 = mill(2,ng)
       ng_2d = imill_2d(n1,n2)
       n3 = mill(3,ng)+1
       IF (n3<1) n3 = n3 + dfftp%nr3    
       rg3 = rhog(ng)
       rhog3(n3,ng_2d)=rg3
       IF ( gamma_only .AND. n1==0 .AND. n2==0 ) THEN
          n3 = -mill(3,ng)+1
          IF (n3<1) n3 = n3 + dfftp%nr3
          rhog3(n3,ng_2d) = CONJG(rg3)
       ENDIF
    ENDDO
    ! End mapping
    !
    vg3(:,:) = (0.d0,0.d0)
    L  = at(3,3)*alat
    z0 = L/2.d0
    z1 = z0+esm_w
    ci = (0.d0,1.d0)
    !
    !****For gp!=0 case ****
    ALLOCATE( vg(dfftp%nr3), vg_r(dfftp%nr3) )
    !
    DO ng_2d = 1, ngm_2d
       k1 = mill_2d(1,ng_2d)
       k2 = mill_2d(2,ng_2d)
       IF (k1==0 .AND. k2==0) CYCLE
       t(1:2) = k1 * bg (1:2, 1) + k2 * bg (1:2, 2)
       gp2 = SUM( t(:) * t(:) ) * tpiba2
       gp = SQRT(gp2)
       tmp1=(0.d0,0.d0); tmp2=(0.d0,0.d0);
       vg(:) = (0.d0,0.d0)
       DO iz = 1, dfftp%nr3
          IF (iz<=dfftp%nr3/2) kn = DBLE(iz-1)*tpi/L
          IF (iz> dfftp%nr3/2) kn = DBLE(iz-1-dfftp%nr3)*tpi/L
          cc0 = COS(kn*z0)
          ss0 = SIN(kn*z0)
          rg3 = rhog3(iz,ng_2d)
          ! bc3
          arg1 =  gp*(z1-z0)
          arg2 = -gp*(z1-z0)
          vg(iz) = fpi*rg3/(gp**2+kn**2)
          tmp = ((gp+ci*kn)*EXP(arg1)+(gp-ci*kn)*EXP(arg2))/(2.d0*gp)
          tmp1 = tmp1+rg3*(cc0+ci*ss0)/(gp**2+kn**2)*tmp
          tmp = (gp-ci*kn)/gp
          tmp2 = tmp2+rg3*(cc0-ci*ss0)/(gp**2+kn**2)*tmp
       ENDDO
       vg3(:,ng_2d)=vg(:)
       !
       ! real part
       vg_r(:)=(0.d0,0.d0)
       DO iz = 1, dfftp%nr3
          k3 = iz-1
          IF (k3 > dfftp%nr3/2) k3 = iz-dfftp%nr3-1
          z = DBLE(k3)/DBLE(dfftp%nr3)*L
          ! bc3
          arg1 =  gp*(z-z1)
          arg2 = -gp*(z+z0)
          arg3 =  gp*(z-z0-2.d0*z1)
          vg_r(iz) = -fpi*EXP(arg1)*tmp1+tpi*(EXP(arg3)-EXP(arg2))*tmp2
       ENDDO
       !
       CALL cft_1z(vg_r,1,dfftp%nr3,dfftp%nr3,-1,vg)
       !
       vg3(:,ng_2d)=(vg3(:,ng_2d)+vg(:))*e2 ! factor e2: hartree -> Ry.
       !
    ENDDO
    !
    DEALLOCATE( vg, vg_r )
    !
    !****For gp=0 case ****
    ng_2d = imill_2d(0,0)
    !
    IF ( ng_2d > 0 ) THEN
       ALLOCATE( vg(dfftp%nr3), vg_r(dfftp%nr3) )
       tmp1=(0.d0,0.d0); tmp2=(0.d0,0.d0); tmp3=(0.d0,0.d0);
       vg(:) = (0.d0,0.d0)
       rg3 = rhog3(1,ng_2d)
       vg(1) = tpi*(4.d0*z1-z0)*z0*rg3
       !
       DO iz=2,dfftp%nr3
          IF (iz<=dfftp%nr3/2) kn = DBLE(iz-1)*tpi/L
          IF (iz> dfftp%nr3/2) kn = DBLE(iz-1-dfftp%nr3)*tpi/L
          rg3 = rhog3(iz,ng_2d)
          cc0 = COS(kn*z0)
          ss0 = SIN(kn*z0)
          tmp1 = tmp1+rg3*(cc0+ci*ss0)/kn**2
          tmp2 = tmp2+rg3*(cc0-ci*ss0)/kn
          tmp3 = tmp3+rg3*(cc0+ci*ss0)/kn
          vg(iz) = fpi*rg3/(kn**2)
       ENDDO
       !
       vg3(:,ng_2d)=vg(:)
       ! 
       ! real part
       vg_r(:) = (0.d0,0.d0) 
       rg3 = rhog3(1,ng_2d)
       !
       DO iz = 1, dfftp%nr3
          k3 = iz-1
          IF (k3 > dfftp%nr3/2) k3=iz-dfftp%nr3-1
          z = DBLE(k3)/DBLE(dfftp%nr3)*L
          vg_r(iz) = -tpi*(z**2+2.d0*z*z0)*rg3 &
                     -fpi*tmp1                 &
                     -fpi*ci*(z-z1)*tmp2       &
                     -fpi*ci*(z1-z0)*tmp3
       ENDDO
       !
       ! start smoothing
       nz_l = dfftp%nr3/2+1+esm_nfit
       nz_r = dfftp%nr3/2+1-esm_nfit
       !
       z_l = DBLE(nz_l-1)*L/DBLE(dfftp%nr3)-L
       z_r = DBLE(nz_r-1)*L/DBLE(dfftp%nr3)
       !
       f1 = -tpi*(z_r**2+2.d0*z_r*z0)*rg3 &
            -fpi*tmp1 &
            -fpi*ci*(z_r-z1)*tmp2 &
            -fpi*ci*(z1 -z0)*tmp3
       f2 = -tpi*(z_l**2+2.d0*z_l*z0)*rg3 &
            -fpi*tmp1 &
            -fpi*ci*(z_l-z1)*tmp2 &
            -fpi*ci*(z1 -z0)*tmp3
       f3 = -tpi*(2.d0*z_r+2.d0*z0)*rg3 &
            -fpi*ci*tmp2
       f4 = -tpi*(2.d0*z_l+2.d0*z0)*rg3 &
            -fpi*ci*tmp2
       !
       z_r = z_r
       z_l = z_l+L
       !
       a0 = (f1*z_l**2*(z_l-3.d0*z_r)+z_r*(f3*z_l**2*(-z_l+z_r) &
            +z_r*(f2*(3.d0*z_l-z_r)+f4*z_l*(-z_l+z_r))))/(z_l-z_r)**3
       a1 = (f3*z_l**3+z_l*(6.d0*f1-6.d0*f2+(f3+2.d0*f4)*z_l)*z_r &
            -(2*f3+f4)*z_l*z_r**2-f4*z_r**3)/(z_l-z_r)**3
       a2 = (-3*f1*(z_l+z_r)+3.d0*f2*(z_l+z_r)-(z_l-z_r)*(2*f3*z_l &
            +f4*z_l+f3*z_r+2*f4*z_r))/(z_l-z_r)**3
       a3 = (2.d0*f1-2.d0*f2+(f3+f4)*(z_l-z_r))/(z_l-z_r)**3
       !
       DO iz = nz_r, nz_l
          z = DBLE(iz-1)/DBLE(dfftp%nr3)*L
          vg_r(iz) = (a0+a1*z+a2*z**2+a3*z**3)
       ENDDO
       ! end smoothing
       !
       CALL cft_1z( vg_r, 1, dfftp%nr3, dfftp%nr3, -1, vg )
       vg3(:,ng_2d) = (vg3(:,ng_2d)+vg(:))*e2 ! factor e2: hartree -> Ry.
       !
       DEALLOCATE( vg, vg_r )
    ENDIF ! if( ng_2d > 0 )
    !
    ! Hartree Energy
    ehart = 0.d0
    eh = 0d0
    !
    DO ng_2d = 1, ngm_2d
       k1 = mill_2d(1,ng_2d)
       k2 = mill_2d(2,ng_2d)
       eh = eh + SUM( vg3(:,ng_2d)*CONJG(rhog3(:,ng_2d)) )
    ENDDO
    !
    ehart = ehart+eh
    !
    IF ( gamma_only ) THEN
       ehart = ehart * 2d0
       ng_2d = imill_2d(0,0)
       IF ( ng_2d > 0 ) THEN
          ehart = ehart - SUM( vg3(:,ng_2d)*CONJG(rhog3(:,ng_2d)) )
       ENDIF
    ENDIF
    !
    ehart = ehart*omega*0.5d0
    !
    CALL mp_sum( ehart, intra_bgrp_comm )
    !
    ! Map to FFT mesh (dfftp%nrx)
    aux = 0.0d0
    !
    DO ng = 1, ngm
       n1 = mill(1,ng)
       n2 = mill(2,ng)
       ng_2d = imill_2d(n1,n2)
       n3 = mill(3,ng) + 1
       IF (n3<1) n3 = n3 + dfftp%nr3
       aux(dfftp%nl(ng)) = aux(dfftp%nl(ng)) + vg3(n3,ng_2d)
       IF (gamma_only) THEN
          aux(dfftp%nlm(ng)) = CONJG(aux(dfftp%nl(ng)))
       ENDIF
    ENDDO
    !
    DEALLOCATE( rhog3, vg3 )
    !
    RETURN
    !
  END SUBROUTINE esm_hartree_bc3
  !
  !
  !-------------------------------------------------------------------------
  SUBROUTINE esm_hartree_bc4( rhog, ehart, aux )
    !------------------------------------------------------------------------
    !! Calculation of Hartree energy and potential - bc4
    !
    USE constants,        ONLY : pi, tpi, fpi, e2
    USE gvect,            ONLY : ngm, mill
    USE cell_base,        ONLY : omega, alat, tpiba2, at, bg
    USE control_flags,    ONLY : gamma_only
    USE mp_bands,         ONLY : intra_bgrp_comm
    USE mp,               ONLY : mp_sum
    USE fft_base,         ONLY : dfftp
    USE fft_scalar,       ONLY : cft_1z
    !
    IMPLICIT NONE
    !
    REAL(DP) :: ehart
    !! Hartree energy
    COMPLEX(DP) :: rhog(ngm)
    !! n(G)
    COMPLEX(DP) :: aux(dfftp%nnr)
    !! v_h(G)
    !
    ! ... local variables
    !
    INTEGER :: k1, k2, k3, iz, ng, n1, n2, n3, nz_r, nz_l, ng_2d
    REAL(DP) :: t(2), z, z0, z1, gp, gp2, kn, cc0, ss0, L, &
                z_l, z_r, eh, aaa, cc1, ss1, alpha, beta,  &
                chi, xi, kappa, lambda, arg1, arg2, arg3,  &
                arg4, argr1, argr2, argr3, argr4, argr5
    COMPLEX(DP) :: ci, f1, f2, f3, f4, a0, a1, a2, a3, c_r, c_l, &
                   s_r, s_l, rg3, tmp, tmp1, tmp2, tmp3, tmp4,   &
                   tmpr1, tmpr2, tmpr3, tmpr4
    COMPLEX(DP), ALLOCATABLE :: rhog3(:,:), vg(:), vg_r(:), vg3(:,:), &
                                vr(:), vr_r(:)
    !
    ALLOCATE( rhog3(dfftp%nr3,ngm_2d), vg3(dfftp%nr3,ngm_2d) )
    !
    ! Map to FFT mesh (dfftp%nr3,ngm_2d)
    rhog3(:,:) = (0.d0,0.d0)
    !
    DO ng=1,ngm
       n1 = mill(1,ng)
       n2 = mill(2,ng)
       ng_2d = imill_2d(n1,n2)
       n3 = mill(3,ng)+1
       IF (n3<1) n3 = n3 + dfftp%nr3    
       rg3 = rhog(ng)
       rhog3(n3,ng_2d)=rg3
       IF ( gamma_only .AND. n1==0 .AND. n2==0 ) THEN
          n3 = -mill(3,ng)+1
          IF (n3<1) n3 = n3 + dfftp%nr3
          rhog3(n3,ng_2d) = CONJG(rg3)
       ENDIF
    ENDDO
    ! End mapping
    !
    vg3(:,:) = (0.d0,0.d0)
    L = at(3,3)*alat
    z0 = L/2.d0
    z1 = z0+esm_w
    aaa = esm_a
    ci = (0.d0,1.d0)
    !
    !****For gp!=0 case ****
    ALLOCATE( vg(dfftp%nr3), vg_r(dfftp%nr3), vr(dfftp%nr3), &
              vr_r(dfftp%nr3) )
    !
    DO ng_2d = 1, ngm_2d
       !
       k1 = mill_2d(1,ng_2d)
       k2 = mill_2d(2,ng_2d)
       !
       IF (k1==0 .AND. k2==0) CYCLE
       !
       t(1:2) = k1*bg(1:2,1) + k2 * bg(1:2,2)
       gp2 = SUM(t(:) * t(:)) * tpiba2
       gp = SQRT(gp2)
       !
       tmp1 =(0.d0,0.d0); tmp2 =(0.d0,0.d0); tmp3 =(0.d0,0.d0)
       tmp4 =(0.d0,0.d0); tmpr1=(0.d0,0.d0); tmpr2=(0.d0,0.d0)
       tmpr3=(0.d0,0.d0); tmpr4=(0.d0,0.d0)
       vr(:)=(0.d0,0.d0); vg(:)=(0.d0,0.d0)
       !
       DO iz = 1, dfftp%nr3
          !
          IF (iz<=dfftp%nr3/2) kn = DBLE(iz-1)*tpi/L
          IF (iz> dfftp%nr3/2) kn = DBLE(iz-1-dfftp%nr3)*tpi/L
          !
          cc0 = COS(kn*z0)
          ss0 = SIN(kn*z0)
          rg3 = rhog3(iz,ng_2d)
          !
          ! bc4
          vg(iz) = fpi*rg3/(gp**2+kn**2)
          vr(iz) = fpi*rg3/(gp**2+kn**2+ci*aaa*kn)
          !
          cc1 = COS(kn*z1)
          ss1 = SIN(kn*z1)
          !
          alpha = aaa+gp+SQRT(aaa**2+gp**2)
          beta  = aaa+gp-SQRT(aaa**2+gp**2)
          kappa = aaa-gp+SQRT(aaa**2+gp**2)
          xi    = aaa   +SQRT(aaa**2+gp**2)
          chi   = aaa   -SQRT(aaa**2+gp**2)
          lambda=        SQRT(aaa**2+gp**2)
          !
          tmp1 = tmp1 + rg3*(cc0+ci*ss0)/(xi-ci*kn)/alpha
          tmp2 = tmp2 + rg3*(cc0-ci*ss0)/(gp+ci*kn)/gp
          tmp3 = tmp3 + rg3*kappa/alpha*(cc0-ci*ss0)/(gp+ci*kn)/gp
          tmp4 = tmp4 + rg3*kappa*(cc1+ci*ss1)/(xi-ci*kn)/(gp**2+kn**2)
          !
          tmpr1 = tmpr1 + rg3*(cc0-ci*ss0)/(gp+ci*kn)/alpha
          tmpr2 = tmpr2 + rg3*(cc0+ci*ss0)/(xi-ci*kn)/lambda
          tmpr3 = tmpr3 + rg3*beta/alpha*(cc0+ci*ss0)/(xi-ci*kn)/lambda
          tmpr4 = tmpr4 + rg3*beta*(cc1+ci*ss1)/(gp+ci*kn) &
                          /(gp**2+kn**2+ci*2.d0*aaa*kn)
          !
       ENDDO
       !
       CALL cft_1z( vg, 1, dfftp%nr3, dfftp%nr3, 1, vg_r )
       ! bc4
       CALL cft_1z( vr, 1, dfftp%nr3, dfftp%nr3, 1, vr_r )
       !
       DO iz = 1, dfftp%nr3
          k3 = iz-1
          IF (k3 > dfftp%nr3/2) k3 = iz-dfftp%nr3-1
          z = DBLE(k3)/DBLE(dfftp%nr3)*L
          ! bc4
          arg1 =  gp*(z-z1)-xi*(z0-z1)
          arg2 = -gp*(z+z0)
          arg3 = -gp*(z0+z1)+gp*(z-z1)
          arg4 =  gp*(z-z1)
          !
          argr1 = -gp*(z0+z1)-xi*(z-z1)
          argr2 = -xi*(z0-z1)-chi*(z-z1)
          argr3 = -xi*(z-z1)-xi*(z0-z1)
          argr4 = -xi*(z-z1)
          argr5 = -2.d0*aaa*(z-z1)
          !
          IF (z < z1) THEN
             vg_r(iz) = vg_r(iz) - fpi*EXP(arg1)*tmp1-tpi*EXP(arg2)*tmp2 &
                                 + tpi*EXP(arg3)*tmp3-fpi*EXP(arg4)*tmp4
          ELSE
             vg_r(iz) = vr_r(iz)*EXP(argr5) &
                       - fpi*EXP(argr1)*tmpr1-tpi*EXP(argr2)*tmpr2 &
                       + tpi*EXP(argr3)*tmpr3-fpi*EXP(argr4)*tmpr4
          ENDIF
       ENDDO
       !
       CALL cft_1z( vg_r, 1, dfftp%nr3, dfftp%nr3, -1, vg )
       !
       vg3(:,ng_2d)=vg(:)*e2 ! factor e2: hartree -> Ry.
       !
    ENDDO
    !
    DEALLOCATE( vg, vg_r, vr, vr_r )
    !
    !****For gp=0 case ****
    ng_2d = imill_2d(0,0)
    !
    IF ( ng_2d > 0 ) THEN
       !
       ALLOCATE( vg(dfftp%nr3), vg_r(dfftp%nr3), vr(dfftp%nr3), vr_r(dfftp%nr3) )
       tmp1=(0.d0,0.d0); tmp2=(0.d0,0.d0); tmp3=(0.d0,0.d0); tmp4=(0.d0,0.d0)
       vg(:)=(0.d0,0.d0); vr(:)=(0.d0,0.d0)
       !
       !for smoothing
       f1=(0.d0,0.d0); f2=(0.d0,0.d0); f3=(0.d0,0.d0); f4=(0.d0,0.d0)
       !
       nz_l = dfftp%nr3/2+1+esm_nfit
       nz_r = dfftp%nr3/2+1-esm_nfit
       !
       z_l = DBLE(nz_l-1)*L/DBLE(dfftp%nr3)-L
       z_r = DBLE(nz_r-1)*L/DBLE(dfftp%nr3)
       !
       rg3 = rhog3(1,ng_2d)
       ! bc4
       arg1  = -2.d0*aaa*(z0-z1)
       vg(1) = tpi*((z0+z1)/aaa+2.d0*z0*z1+z1**2)*rg3 &
              - pi *(EXP(arg1)-1.d0)/aaa**2*rg3
       vr(1) = tpi*(z0+0.5d0/aaa)/aaa*rg3
       !
       DO iz = 2, dfftp%nr3
          !
          IF (iz<=dfftp%nr3/2) kn = DBLE(iz-1)*tpi/L
          IF (iz> dfftp%nr3/2) kn = DBLE(iz-1-dfftp%nr3)*tpi/L
          rg3 = rhog3(iz,ng_2d)
          !
          ! bc4
          cc0 = COS(kn*z0)
          ss0 = SIN(kn*z0)
          cc1 = COS(kn*z1)
          ss1 = SIN(kn*z1)
          tmp1 = tmp1+rg3*(cc1+ci*ss1)/(2.d0*aaa-ci*kn)/kn**2
          tmp2 = tmp2+rg3*(cc0-ci*ss0)/kn
          tmp3 = tmp3+rg3*(cc0+ci*ss0)/(2.d0*aaa-ci*kn)
          tmp4 = tmp4+(0.d0,0.d0)
          !
          vg(iz)=fpi*rg3/(kn**2)
          ! bc4
          vr(iz)=fpi*rg3/(kn**2+ci*2.d0*aaa*kn)
          !
          !for smoothing
          c_r = COS(kn*z_r)
          s_r = SIN(kn*z_r)
          c_l = COS(kn*z_l)
          s_l = SIN(kn*z_l)
          ! bc4
          f1 = f1+fpi*   rg3*(c_r+ci*s_r)/(kn**2+ci*2.d0*aaa*kn)
          f2 = f2+fpi*   rg3*(c_l+ci*s_l)/kn**2
          f3 = f3+fpi*ci*rg3*(c_r+ci*s_r)/kn
          f4 = f4+fpi*ci*rg3*(c_l+ci*s_l)/kn
          !
       ENDDO
       !
       CALL cft_1z( vg, 1, dfftp%nr3, dfftp%nr3, 1, vg_r )
       ! bc4
       CALL cft_1z( vr, 1, dfftp%nr3, dfftp%nr3, 1, vr_r )
       !
       rg3 = rhog3(1,ng_2d)
       !
       DO iz = 1, dfftp%nr3
          k3 = iz-1
          IF (k3 > dfftp%nr3/2) k3 = iz-dfftp%nr3-1
          z = DBLE(k3)/DBLE(dfftp%nr3)*L
          ! bc4
          arg1 = -2.d0*aaa*(z0-z1)
          arg2 = -2.d0*aaa*(z-z1)
          !
          IF (z < z1) THEN
             vg_r(iz) = vg_r(iz)   &
                -fpi*2.d0*aaa*tmp1 &
                -tpi*ci*(2.d0*(z -z1)-1.d0/aaa)*tmp2 &
                -tpi*EXP(arg1)/aaa*tmp3 &
                -tpi*z*(z+2.d0*z0)*rg3
          ELSE
             vg_r(iz) = vr_r(iz)*EXP(arg2) &
                +tpi*ci*EXP(arg2)/aaa*tmp2 &
                -tpi*EXP(arg1)/aaa*tmp3    &
                +tpi*EXP(arg2)*z/aaa*rg3   &
                -pi *EXP(arg1)/aaa**2*rg3
          ENDIF
       ENDDO
       !
       !for smoothing
       ! bc4
       arg1 = -2.d0*aaa*(z0-z1)
       arg2 = -2.d0*aaa*(z_r-z1)
       !
       f1 = f1 + tpi*(z0+0.5d0/aaa)/aaa*rg3
       f1 = f1 * EXP(arg2) &
            +tpi*ci*EXP(arg2)/aaa*tmp2 &
            -tpi*EXP(arg1)/aaa*tmp3    &
            +tpi*EXP(arg2)*z_r/aaa*rg3 &
            -pi *EXP(arg1)/aaa**2*rg3
       f2 = f2 + tpi*((z0+z1)/aaa+2.d0*z0*z1+z1**2)*rg3 &
            -pi *(EXP(arg1)-1.d0)/aaa**2*rg3
       f2 = f2 &
            -fpi*2.d0*aaa*tmp1 &
            -tpi*ci*(2.d0*(z_l-z1)-1.d0/aaa)*tmp2 &
            -tpi*EXP(arg1)/aaa*tmp3 &
            -tpi*z_l*(z_l+2.d0*z0)*rg3
       f3 = f3*EXP(arg2) &
            -fpi*ci*EXP(arg2)*tmp2 &
            -fpi*(z_r+z0)*EXP(arg2)*rg3 
       f4 = f4-fpi*ci*tmp2-fpi*(z_l+z0)*rg3
       !
       ! for smoothing
       ! factor e2 will be multiplied later (at vg3 <= vg)
       ! f1=f1*e2; f2=f2*e2; f3=f3*e2; f4=f4*e2
       !
       z_r = z_r
       z_l = z_l+L
       !
       a0 = (f1*z_l**2*(z_l-3.d0*z_r)+z_r*(f3*z_l**2*(-z_l+z_r)    &
            +z_r*(f2*(3.d0*z_l-z_r)+f4*z_l*(-z_l+z_r))))/(z_l-z_r)**3
       a1 = (f3*z_l**3+z_l*(6.d0*f1-6.d0*f2+(f3+2.d0*f4)*z_l)*z_r  &
            -(2*f3+f4)*z_l*z_r**2-f4*z_r**3)/(z_l-z_r)**3
       a2 = (-3*f1*(z_l+z_r)+3.d0*f2*(z_l+z_r)-(z_l-z_r)*(2*f3*z_l &
            +f4*z_l+f3*z_r+2*f4*z_r))/(z_l-z_r)**3
       a3 = (2.d0*f1-2.d0*f2+(f3+f4)*(z_l-z_r))/(z_l-z_r)**3
       !
       DO iz = nz_r, nz_l
          z = DBLE(iz-1)/DBLE(dfftp%nr3)*L
          vg_r(iz) = (a0+a1*z+a2*z**2+a3*z**3)
       ENDDO
       !
       CALL cft_1z(vg_r,1,dfftp%nr3,dfftp%nr3,-1,vg)
       !
       vg3(:,ng_2d) = vg(:)*e2 ! factor e2: hartree -> Ry.
       !
       DEALLOCATE( vg, vg_r, vr, vr_r )
    ENDIF ! if( ng_2d > 0 )
    !
    ! Hartree Energy
    ehart=0.d0
    eh = 0d0
    DO ng_2d = 1, ngm_2d
       k1 = mill_2d(1,ng_2d)
       k2 = mill_2d(2,ng_2d)
       eh = eh + SUM( vg3(:,ng_2d)*CONJG(rhog3(:,ng_2d)) )
    ENDDO
    !
    ehart = ehart+eh
    !
    IF ( gamma_only ) THEN
       ehart = ehart * 2d0
       ng_2d = imill_2d(0,0)
       IF ( ng_2d > 0 ) THEN
          ehart = ehart - SUM( vg3(:,ng_2d)*CONJG(rhog3(:,ng_2d)) )
       ENDIF
    ENDIF
    ehart = ehart*omega*0.5d0
    !
    CALL mp_sum( ehart, intra_bgrp_comm )
    !
    ! Map to FFT mesh (dfftp%nrx)
    aux = 0.0d0
    !
    DO ng = 1, ngm
       n1 = mill(1,ng)
       n2 = mill(2,ng)
       ng_2d = imill_2d(n1,n2)
       n3 = mill(3,ng) + 1
       IF (n3<1) n3 = n3 + dfftp%nr3
       aux(dfftp%nl(ng)) = aux(dfftp%nl(ng)) + vg3(n3,ng_2d)
       IF (gamma_only) THEN
          aux(dfftp%nlm(ng)) = CONJG(aux(dfftp%nl(ng)))
       ENDIF
    ENDDO
    !
    DEALLOCATE( rhog3, vg3 )
    !
    RETURN
    !
  END SUBROUTINE esm_hartree_bc4
  !
  !
  !----------------------------------------------------------------------
  COMPLEX(DP) FUNCTION qe_exp( x )
    !--------------------------------------------------------------------
    !
    COMPLEX(DP), INTENT(IN) :: x
    REAL(DP) :: r, i, c, s
    !
    r = DREAL(x)
    i = DIMAG(x)
    c = COS(i)
    s = SIN(i)
    !
    qe_exp = EXP(r)*CMPLX(c,s,KIND=DP)
    !
  END FUNCTION qe_exp
  !
  !
  !-------------------------------------------------------------------------
  COMPLEX(DP) FUNCTION qe_sinh( x )
    !-----------------------------------------------------------------------
    !
    COMPLEX(DP), INTENT(IN) :: x
    REAL(DP) :: r, i, c, s
    !
    r = DREAL(x)
    i = DIMAG(x)
    c = COS(i)
    s = SIN(i)
    !
    qe_sinh = 0.5d0*(EXP(r)*CMPLX(c,s,KIND=DP) - EXP(-r)*CMPLX(c,-s,KIND=DP) )
    !
  END FUNCTION qe_sinh
  !
  !
  !-------------------------------------------------------------------------
  COMPLEX(DP) FUNCTION qe_cosh( x )
    !-----------------------------------------------------------------------
    !
    COMPLEX(DP), INTENT(IN) :: x
    REAL(DP) :: r, i, c, s
    !
    r = DREAL(x)
    i = DIMAG(x)
    c = COS(i)
    s = SIN(i)
    !
    qe_cosh = 0.5d0*( EXP(r)*CMPLX(c,s,KIND=DP) + EXP(-r)*CMPLX(c,-s,KIND=DP) )
    !
  END FUNCTION qe_cosh
  !
  !
  !---------------------------------------------------------------------------
  SUBROUTINE esm_stres_har_bc1( sigmahar, rhog )
    !-------------------------------------------------------------------------
    !!  Calculation of Hartree stress tensor - bc1.
    !
    USE kinds,           ONLY : DP
    USE gvect,           ONLY : ngm, mill
    USE constants,       ONLY : tpi, fpi, e2
    USE cell_base,       ONLY : omega, alat, at, tpiba, bg
    USE control_flags,   ONLY : gamma_only
    USE fft_base,        ONLY : dfftp
    USE fft_scalar,      ONLY : cft_1z
    USE mp_bands,        ONLY : intra_bgrp_comm
    USE mp,              ONLY : mp_sum
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(OUT) :: sigmahar(3,3)
    !! Hartree stress
    COMPLEX(DP), INTENT(IN) :: rhog(ngm)
    !! n(G)
    !
    ! ... local variables
    !
    INTEGER :: ig, iga, igb, igz, igp, la, mu, iz
    REAL(DP) :: L, S, z0, z
    REAL(DP) :: g(2), gp, gz
    COMPLEX(DP), PARAMETER :: ci = DCMPLX(0.0d0, 1.0d0)
    COMPLEX(DP) :: rg3
    COMPLEX(DP) :: sum1p, sum1m, sum1c, sum2c, sum2p, sum2m
    REAL(DP)    :: z_l, z_r
    COMPLEX(DP) :: f1, f2, f3, f4, a0, a1, a2, a3
    COMPLEX(DP) :: poly_fr, poly_fl, poly_dfr, poly_dfl
    COMPLEX(DP) :: poly_a, poly_b, poly_c, poly_d
    REAL(DP), PARAMETER :: delta(2,2) = RESHAPE( (/ 1.0d0, 0.0d0, 0.0d0, 1.0d0 /), (/2,2/) )
    REAL(DP) :: dgp_deps(2,2)     ! dgp/deps
    REAL(DP) :: dgp2_deps(2,2)    ! dgp^2/deps
    REAL(DP) :: dinvgp_deps(2,2)  ! dgp^-1/deps
    !
    COMPLEX(DP), ALLOCATABLE :: rhog3(:,:)
    COMPLEX(DP), ALLOCATABLE :: dVr_deps(:,:,:)
    COMPLEX(DP), ALLOCATABLE :: dVg_deps(:,:,:)
    COMPLEX(DP), ALLOCATABLE :: Vr(:)
    COMPLEX(DP), ALLOCATABLE :: Vg(:)
    !
    ! cell settings
    L  = at(3,3)*alat
    S  = omega/L
    z0 = L/2.d0
    !
    ! initialize
    sigmahar(:,:) = 0.0d0
    !
    ALLOCATE( rhog3(dfftp%nr3,ngm_2d) )
    ALLOCATE( dVr_deps(dfftp%nr3,2,2) )
    ALLOCATE( dVg_deps(dfftp%nr3,2,2) )
    ALLOCATE( Vr(dfftp%nr3) )
    ALLOCATE( Vg(dfftp%nr3) )
    !
    ! reconstruct rho(gz,gp)
    rhog3(:,:)=(0.d0,0.d0)
    !
    DO ig = 1, ngm
       !
       iga = mill(1,ig)
       igb = mill(2,ig)
       igz = mill(3,ig)+1
       igp = imill_2d(iga,igb)
       !
       IF ( igz<1 ) THEN
          igz = igz + dfftp%nr3
       ENDIF
       !
       rg3 = rhog(ig)
       rhog3(igz,igp) = rg3
       !
       ! expand function symmetrically to gz<0
       IF ( gamma_only .AND. iga==0 .AND. igb==0 ) THEN
          igz = 1-mill(3,ig)
          if( igz<1 ) THEN
             igz = igz + dfftp%nr3
          ENDIF
          rhog3(igz,igp) = CONJG(rg3)
       ENDIF
    ENDDO ! ig
    !
    !
    !****For gp!=0 case ****
    DO igp=1, ngm_2d
       iga = mill_2d(1,igp)
       igb = mill_2d(2,igp)
       g(1:2) = (iga*bg(1:2,1) + igb*bg(1:2,2))*tpiba
       gp = SQRT(g(1)*g(1) + g(2)*g(2))
       !
       IF ( gp==0.0d0 ) CYCLE ! skip gp=0
       !
       ! derivatives by strain tensor
       DO la = 1, 2
          DO mu = 1, 2
             dgp_deps(la,mu)    = -g(la)*g(mu)/gp
             dgp2_deps(la,mu)   = -g(la)*g(mu)*2.0d0
             dinvgp_deps(la,mu) = +g(la)*g(mu)/gp**3
          ENDDO
       ENDDO
       !
       ! summations over gz
       sum1p = (0.d0,0.d0)
       sum1m = (0.d0,0.d0)
       sum2p = (0.d0,0.d0)
       sum2m = (0.d0,0.d0)
       !
       DO igz = 1, dfftp%nr3
          IF ( igz<=dfftp%nr3/2 )THEN
             gz = DBLE(igz-1)*tpi/L
          ELSE
             gz = DBLE(igz-1-dfftp%nr3)*tpi/L
          ENDIF
          !
          rg3 = rhog3(igz,igp)
          sum1p = sum1p + rg3*qe_EXP(+ci*gz*z0)/(gp-ci*gz)
          sum1m = sum1m + rg3*qe_EXP(-ci*gz*z0)/(gp+ci*gz)
          sum2p = sum2p + rg3*qe_EXP(+ci*gz*z0)/(gp-ci*gz)**2
          sum2m = sum2m + rg3*qe_EXP(-ci*gz*z0)/(gp+ci*gz)**2
          !
       ENDDO ! igz
       !
       ! calculate dV(z)/deps
       DO iz = 1, dfftp%nr3
          IF ( iz-1 > dfftp%nr3/2 ) THEN
             z = DBLE(iz-1-dfftp%nr3)/DBLE(dfftp%nr3)*L
          ELSE
             z = DBLE(iz-1)/DBLE(dfftp%nr3)*L
          ENDIF
          !
          dVr_deps(iz,:,:) = &
               -( dgp_deps(:,:)*tpi/gp**2*(gp*(z-z0)-1.0d0) &
               -  delta(:,:)*tpi/gp )    &
               * EXP(+gp*(z-z0)) * sum1p &
               + dgp_deps(:,:)*tpi/gp * EXP(+gp*(z-z0)) * sum2p &
               +( dgp_deps(:,:)*tpi/gp**2*(gp*(z+z0)+1.0d0)     &
               +  delta(:,:)*tpi/gp )    &
               * EXP(-gp*(z+z0)) * sum1m &
               + dgp_deps(:,:)*tpi/gp * EXP(-gp*(z+z0)) * sum2m
       ENDDO ! iz
       !
       ! convert dV(z)/deps to dV(gz)/deps
       DO la = 1, 2
          DO mu = 1, 2
             CALL cft_1z( dVr_deps(:,la,mu), 1, dfftp%nr3, dfftp%nr3, -1, &
                          dVg_deps(:,la,mu) )
          ENDDO
       ENDDO
       !
       ! add bare coulomn terms to dV(gz)/deps
       DO igz = 1, dfftp%nr3
          !
          IF ( igz<=dfftp%nr3/2 ) THEN
             gz = DBLE(igz-1)*tpi/L
          ELSE
             gz = DBLE(igz-1-dfftp%nr3)*tpi/L
          ENDIF
          !
          rg3 = rhog3(igz,igp)
          !
          dVg_deps(igz,:,:) = dVg_deps(igz,:,:) &
               - delta(:,:)* fpi*rg3/(gp**2+gz**2) &
               - dgp2_deps(:,:)* fpi*rg3/(gp**2+gz**2)**2
       ENDDO ! igz
       !
       ! modifications
       IF ( gamma_only ) THEN
          dVg_deps(:,:,:) = dVg_deps(:,:,:)*2.0d0
       ENDIF
       !
       ! calculate stress tensor
       DO igz=1, dfftp%nr3
          rg3 = rhog3(igz,igp)
          sigmahar(1:2,1:2) = sigmahar(1:2,1:2) &
               + DBLE( CONJG(rg3) * dVg_deps(igz,:,:) )
       ENDDO ! igz
    ENDDO ! igp
    !
    !****For gp=0 case ****
    IF ( imill_2d(0,0) > 0 ) THEN
       ! summations over gz
       sum1c = (0.d0,0.d0)
       sum2c = (0.d0,0.d0)
       DO igz = 2, dfftp%nr3
          IF ( igz<=dfftp%nr3/2 ) THEN
             gz = DBLE(igz-1)*tpi/L
          ELSE
             gz = DBLE(igz-1-dfftp%nr3)*tpi/L
          ENDIF
          !
          rg3 = rhog3(igz,imill_2d(0,0))
          sum1c = sum1c + rg3*ci*COS(gz*z0)/gz
          sum2c = sum2c + rg3*COS(gz*z0)/gz**2
       ENDDO ! igz
       !
       ! calculate V(z)
       DO iz = 1, dfftp%nr3
          IF ( iz-1 > dfftp%nr3/2 ) THEN
             z = DBLE(iz-1-dfftp%nr3)/DBLE(dfftp%nr3)*L
          ELSE
             z = DBLE(iz-1)/DBLE(dfftp%nr3)*L
          ENDIF
          !
          rg3 = rhog3(1,imill_2d(0,0))
          Vr(iz) = &
               - tpi*z**2*rg3 &
               - tpi*z0**2*rg3 &
               - fpi*z*sum1c &
               - fpi*sum2c
       ENDDO ! iz
       !
       ! separation by polynomial
       z_l = -z0
       z_r = +z0
       !
       f1 = -tpi*z_r**2*rg3 &
            -tpi*z0**2*rg3 &
            -fpi*z_r*sum1c &
            -fpi*sum2c
       f2 = -tpi*z_l**2*rg3 &
            -tpi*z0**2*rg3 &
            -fpi*z_l*sum1c &
            -fpi*sum2c
       f3 = -fpi*z_r*rg3 &
            -fpi*sum1c
       f4 = -fpi*z_l*rg3 &
            -fpi*sum1c
       !
       a0=(f1*z_l**2*(z_l-3.d0*z_r)+z_r*(f3*z_l**2*(-z_l+z_r) &
            +z_r*(f2*(3.d0*z_l-z_r)+f4*z_l*(-z_l+z_r))))/(z_l-z_r)**3
       a1=(f3*z_l**3+z_l*(6.d0*f1-6.d0*f2+(f3+2.d0*f4)*z_l)*z_r &
            -(2*f3+f4)*z_l*z_r**2-f4*z_r**3)/(z_l-z_r)**3
       a2=(-3*f1*(z_l+z_r)+3.d0*f2*(z_l+z_r)-(z_l-z_r)*(2*f3*z_l &
            +f4*z_l+f3*z_r+2*f4*z_r))/(z_l-z_r)**3
       a3=(2.d0*f1-2.d0*f2+(f3+f4)*(z_l-z_r))/(z_l-z_r)**3
       !
       ! remove polynomial from V(z)
       DO iz = 1, dfftp%nr3
          IF ( iz-1 > dfftp%nr3/2 ) THEN
             z = DBLE(iz-1-dfftp%nr3)/DBLE(dfftp%nr3)*L
          ELSE
             z = DBLE(iz-1)/DBLE(dfftp%nr3)*L
          ENDIF
          Vr(iz) = Vr(iz) - (a0+a1*z+a2*z**2+a3*z**3)
       ENDDO
       !
       ! convert V(z) to V(gz) without polynomial
       CALL cft_1z( Vr, 1, dfftp%nr3, dfftp%nr3, -1, Vg )
       !
       ! add polynomial to V(gz)
       DO igz = 2, dfftp%nr3
          IF ( igz<=dfftp%nr3/2 ) THEN
             gz = DBLE(igz-1)*tpi/L
          ELSE
             gz = DBLE(igz-1-dfftp%nr3)*tpi/L
          ENDIF
          !
          Vg(igz) = Vg(igz) &
               + a1 * ci * COS(gz*z0)/gz &
               + a2 * 2.0d0 * COS(gz*z0)/gz**2 &
               + a3 * ci * z0**2 * COS(gz*z0)/gz &
               - a3 * ci * 6.0d0 * COS(gz*z0)/gz**3
       ENDDO
       !
       Vg(1) = Vg(1) + a0*1.0d0 + a2*z0**2/3.0d0
       !
       ! add bare coulomn terms to V(gz)
       DO igz = 2, dfftp%nr3
          IF ( igz<=dfftp%nr3/2 ) THEN
             gz = DBLE(igz-1)*tpi/L
          ELSE
             gz = DBLE(igz-1-dfftp%nr3)*tpi/L
          ENDIF
          !
          rg3 = rhog3(igz,imill_2d(0,0))
          Vg(igz) = Vg(igz) + fpi*rg3/gz**2
       ENDDO ! igz
       !
       ! calculate dV/deps(gz)
       DO igz=1, dfftp%nr3
          dVg_deps(igz,:,:) = -delta(:,:) * Vg(igz)
       ENDDO ! igz
       !
       ! calculate stress tensor
       DO igz=1, dfftp%nr3
          rg3 = rhog3(igz,imill_2d(0,0))
          sigmahar(1:2,1:2) = sigmahar(1:2,1:2) &
               + DBLE( CONJG(rg3) * dVg_deps(igz,1:2,1:2) )
       ENDDO ! igz
       !
    ENDIF ! imill_2d(0,0) > 0
    !
    ! half means removing duplications.
    ! e2 means hartree -> Ry.
    sigmahar(:,:) = sigmahar(:,:) * (-0.5d0*e2)
    !
    CALL mp_sum( sigmahar, intra_bgrp_comm )
    !
    DEALLOCATE( rhog3 )
    DEALLOCATE( dVr_deps )
    DEALLOCATE( dVg_deps )
    DEALLOCATE( Vr )
    DEALLOCATE( Vg )
    !
    RETURN
    !
  END SUBROUTINE esm_stres_har_bc1
  !
  !
  !------------------------------------------------------------------
  SUBROUTINE esm_stres_har_bc2( sigmahar, rhog )
    !----------------------------------------------------------------
    !! Calculation of Hartree stress tensor - bc2.
    !
    USE kinds,         ONLY : DP
    USE gvect,         ONLY : ngm, mill
    USE constants,     ONLY : tpi, fpi, e2
    USE cell_base,     ONLY : omega, alat, at, tpiba, bg
    USE control_flags, ONLY : gamma_only
    USE fft_base,      ONLY : dfftp
    USE fft_scalar,    ONLY : cft_1z
    USE mp_bands,      ONLY : intra_bgrp_comm
    USE mp,            ONLY : mp_sum
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(OUT)   :: sigmahar(3,3)
    !! the Hartree term of the stress tensor
    COMPLEX(DP), INTENT(IN) :: rhog(ngm)
    !! n(G)
    !
    ! ... local variables
    !
    INTEGER :: ig, iga, igb, igz, igp, la, mu, iz
    REAL(DP) :: L, S, z0, z1, z
    REAL(DP) :: g(2), gp, gz
    COMPLEX(DP), PARAMETER :: ci = DCMPLX(0.0d0, 1.0d0)
    COMPLEX(DP) :: rg3
    COMPLEX(DP) :: sum1p, sum1m, sum1c, sum2c, sum2p, sum2m
    COMPLEX(DP) :: sum1sp, sum1sm, sum1cp, sum1cm, sum2sp, sum2sm
    REAL(DP)    :: z_l, z_r
    COMPLEX(DP) :: f1, f2, f3, f4, a0, a1, a2, a3
    REAL(DP), PARAMETER :: delta(2,2) = RESHAPE( (/ 1.0d0, 0.0d0, 0.0d0, 1.0d0 /), (/2,2/) )
    REAL(DP) :: dgp_deps(2,2)     !! dgp/deps
    REAL(DP) :: dgp2_deps(2,2)    !! dgp^2/deps
    REAL(DP) :: dinvgp_deps(2,2)  !! dgp^-1/deps
    !
    COMPLEX(DP), ALLOCATABLE :: rhog3(:,:)
    COMPLEX(DP), ALLOCATABLE :: dVr_deps(:,:,:)
    COMPLEX(DP), ALLOCATABLE :: dVg_deps(:,:,:)
    COMPLEX(DP), ALLOCATABLE :: Vr(:)
    COMPLEX(DP), ALLOCATABLE :: Vg(:)
    !
    ! cell settings
    L  = at(3,3)*alat
    S  = omega/L
    z0 = L/2.d0
    z1 = z0+esm_w
    !
    ! initialize
    sigmahar(:,:) = 0.0d0
    !
    ALLOCATE( rhog3(dfftp%nr3,ngm_2d) )
    ALLOCATE( dVr_deps(dfftp%nr3,2,2) )
    ALLOCATE( dVg_deps(dfftp%nr3,2,2) )
    ALLOCATE( Vr(dfftp%nr3) )
    ALLOCATE( Vg(dfftp%nr3) )
    !
    ! reconstruct rho(gz,gp)
    rhog3(:,:) = (0.d0,0.d0)
    !
    DO ig=1, ngm
       iga = mill(1,ig)
       igb = mill(2,ig)
       igz = mill(3,ig)+1
       igp = imill_2d(iga,igb)
       IF ( igz<1 ) THEN
          igz = igz + dfftp%nr3
       ENDIF
       rg3 = rhog(ig)
       rhog3(igz,igp) = rg3
       !
       ! expand function symmetrically to gz<0
       IF ( gamma_only .AND. iga==0 .AND. igb==0 ) THEN
          igz = 1-mill(3,ig)
          IF ( igz<1 ) THEN
             igz = igz + dfftp%nr3
          ENDIF
          rhog3(igz,igp) = CONJG(rg3)
       ENDIF
       !
    ENDDO ! ig
    !
    !****For gp!=0 case ****
    DO igp = 1, ngm_2d
       !
       iga = mill_2d(1,igp)
       igb = mill_2d(2,igp)
       g(1:2) = (iga*bg(1:2,1) + igb*bg(1:2,2))*tpiba
       gp = SQRT(g(1)*g(1) + g(2)*g(2))
       !
       IF ( gp==0.0d0 ) CYCLE ! skip gp=0
       !
       ! derivatives by strain tensor
       DO la = 1, 2
          DO mu = 1, 2
             dgp_deps(la,mu)    = -g(la)*g(mu)/gp
             dgp2_deps(la,mu)   = -g(la)*g(mu)*2.0d0
             dinvgp_deps(la,mu) = +g(la)*g(mu)/gp**3
          ENDDO
       ENDDO
       !
       ! summations over gz
       sum1p = (0.d0,0.d0)
       sum1m = (0.d0,0.d0)
       sum2p = (0.d0,0.d0)
       sum2m = (0.d0,0.d0)
       sum1sp = (0.d0,0.d0)
       sum1sm = (0.d0,0.d0)
       sum1cp = (0.d0,0.d0)
       sum1cm = (0.d0,0.d0)
       sum2sp = (0.d0,0.d0)
       sum2sm = (0.d0,0.d0)
       !
       DO igz = 1, dfftp%nr3
          IF ( igz<=dfftp%nr3/2 )THEN
             gz = DBLE(igz-1)*tpi/L
          ELSE
             gz = DBLE(igz-1-dfftp%nr3)*tpi/L
          ENDIF
          !
          rg3 = rhog3(igz,igp)
          !
          sum1p = sum1p + rg3*qe_EXP(+ci*gz*z0)/(gp-ci*gz)
          sum1m = sum1m + rg3*qe_EXP(-ci*gz*z0)/(gp+ci*gz)
          sum2p = sum2p + rg3*qe_EXP(+ci*gz*z0)/(gp-ci*gz)**2
          sum2m = sum2m + rg3*qe_EXP(-ci*gz*z0)/(gp+ci*gz)**2
          !
          sum1sp = sum1sp + rg3*qe_sinh(gp*z0+ci*gz*z0)/(gp+ci*gz)
          sum1sm = sum1sm + rg3*qe_sinh(gp*z0-ci*gz*z0)/(gp-ci*gz)
          !
          sum1cp = sum1cp + rg3*qe_cosh(gp*z0+ci*gz*z0)/(gp+ci*gz)*z0
          sum1cm = sum1cm + rg3*qe_cosh(gp*z0-ci*gz*z0)/(gp-ci*gz)*z0
          !
          sum2sp = sum2sp + rg3*qe_sinh(gp*z0+ci*gz*z0)/(gp+ci*gz)**2
          sum2sm = sum2sm + rg3*qe_sinh(gp*z0-ci*gz*z0)/(gp-ci*gz)**2
       ENDDO ! igz
       !
       ! calculate dV(z)/deps
       DO iz = 1, dfftp%nr3
          !
          IF ( iz-1 > dfftp%nr3/2 ) THEN
             z = DBLE(iz-1-dfftp%nr3)/DBLE(dfftp%nr3)*L
          ELSE
             z = DBLE(iz-1)/DBLE(dfftp%nr3)*L
          ENDIF
          !
          ! BC1 terms
          dVr_deps(iz,:,:) = &
               -( dgp_deps(:,:)*tpi/gp**2*(gp*(z-z0)-1.0d0) &
               -  delta(:,:)*tpi/gp ) &
               * EXP(+gp*(z-z0)) * sum1p &
               + dgp_deps(:,:)*tpi/gp * EXP(+gp*(z-z0)) * sum2p &
               +( dgp_deps(:,:)*tpi/gp**2*(gp*(z+z0)+1.0d0) &
               +  delta(:,:)*tpi/gp ) &
               * EXP(-gp*(z+z0)) * sum1m &
               + dgp_deps(:,:)*tpi/gp * EXP(-gp*(z+z0)) * sum2m
          !
          ! BC2 terms
          dVr_deps(iz,:,:) = dVr_deps(iz,:,:) &
               + dgp_deps(:,:) * ( &
               -   tpi/gp**2 * ( EXP(-gp*(z+2*z1))-EXP(+gp*z) )/SINH(2*gp*z1) &
               +   tpi/gp * ( (-z-2*z1)*EXP(-gp*(z+2*z1)) - z*EXP(+gp*z) )/SINH(2*gp*z1) &
               -   tpi/gp * ( EXP(-gp*(z+2*z1))-EXP(+gp*z) )/SINH(2*gp*z1)**2 * 2*z1     &
               *   COSH(2*gp*z1) ) * sum1sp   &
               + tpi/gp * ( EXP(-gp*(z+2*z1))-EXP(+gp*z) )/SINH(2*gp*z1)    &
               * ( -delta(:,:)*sum1sp + dgp_deps(:,:)*( sum1cp - sum2sp ) ) &
               + dgp_deps(:,:) * ( &
               -   tpi/gp**2 * (EXP(+gp*(z-2*z1))-EXP(-gp*z) )/SINH(2*gp*z1) &
               +   tpi/gp * ( (+z-2*z1)*EXP(+gp*(z-2*z1)) + z*EXP(-gp*z) )/SINH(2*gp*z1) &
               -   tpi/gp * ( EXP(+gp*(z-2*z1))-EXP(-gp*z) )/SINH(2*gp*z1)**2 * 2*z1     &
               *   COSH(2*gp*z1) ) * sum1sm   &
               + tpi/gp * ( EXP(+gp*(z-2*z1))-EXP(-gp*z) )/SINH(2*gp*z1) &
               * ( -delta(:,:)*sum1sm + dgp_deps(:,:)*( sum1cm - sum2sm ) )
       ENDDO ! iz
       !
       ! convert dV(z)/deps to dV(gz)/deps
       DO la=1, 2
          DO mu=1, 2
             CALL cft_1z( dVr_deps(:,la,mu), 1, dfftp%nr3, dfftp%nr3, -1, dVg_deps(:,la,mu) )
          ENDDO
       ENDDO
       !
       ! add bare couloum terms to dV(gz)/deps
       DO igz = 1, dfftp%nr3
          IF ( igz<=dfftp%nr3/2 )THEN
             gz = DBLE(igz-1)*tpi/L
          ELSE
             gz = DBLE(igz-1-dfftp%nr3)*tpi/L
          ENDIF
          rg3 = rhog3(igz,igp)
          !
          dVg_deps(igz,:,:) = dVg_deps(igz,:,:)    &
               - delta(:,:)* fpi*rg3/(gp**2+gz**2) &
               - dgp2_deps(:,:)* fpi*rg3/(gp**2+gz**2)**2
       ENDDO ! igz
       !
       ! modifications
       IF ( gamma_only ) THEN
          dVg_deps(:,:,:) = dVg_deps(:,:,:)*2.0d0
       ENDIF
       !
       ! calculate stress tensor
       DO igz = 1, dfftp%nr3
          rg3 = rhog3(igz,igp)
          sigmahar(1:2,1:2) = sigmahar(1:2,1:2) &
               + DBLE( CONJG(rg3) * dVg_deps(igz,:,:) )
       ENDDO ! igz
       !
    ENDDO ! igp
    !
    !****For gp=0 case ****
    IF ( imill_2d(0,0) > 0 ) THEN
       ! summations over gz
       sum1c = (0.d0,0.d0)
       sum2c = (0.d0,0.d0)
       !
       DO igz = 2, dfftp%nr3
          IF ( igz<=dfftp%nr3/2 ) THEN
             gz = DBLE(igz-1)*tpi/L
          ELSE
             gz = DBLE(igz-1-dfftp%nr3)*tpi/L
          ENDIF
          !
          rg3 = rhog3(igz,imill_2d(0,0))
          sum1c = sum1c + rg3*ci*COS(gz*z0)/gz
          sum2c = sum2c + rg3*COS(gz*z0)/gz**2
       ENDDO ! igz
       !
       ! calculate V(z)
       DO iz = 1, dfftp%nr3
          IF ( iz-1 > dfftp%nr3/2 ) THEN
             z = DBLE(iz-1-dfftp%nr3)/DBLE(dfftp%nr3)*L
          ELSE
             z = DBLE(iz-1)/DBLE(dfftp%nr3)*L
          ENDIF
          !
          rg3 = rhog3(1,imill_2d(0,0))
          !
          ! BC1 terms
          Vr(iz) = &
               - tpi*z**2*rg3  &
               - tpi*z0**2*rg3 &
               - fpi*z*sum1c   &
               - fpi*sum2c
          !
          ! BC2 terms
          Vr(iz) = Vr(iz) &
               + tpi*z1*2*z0 * rg3 - tpi*(-z/z1)*2*z0*sum1c
       ENDDO ! iz
       !
       ! separation by polynomial
       z_l=-z0
       z_r=+z0
       !
       f1 = -tpi*z_r**2*rg3 &
            -tpi*z0**2*rg3  &
            -fpi*z_r*sum1c  &
            -fpi*sum2c
       f1 = f1 &
            + tpi*z1*2*z0 * rg3 - tpi*(-z_r/z1)*2*z0*sum1c

       f2 = -tpi*z_l**2*rg3 &
            -tpi*z0**2*rg3  &
            -fpi*z_l*sum1c  &
            -fpi*sum2c
       f2 = f2 &
            + tpi*z1*2*z0 * rg3 - tpi*(-z_l/z1)*2*z0*sum1c

       f3 = -fpi*z_r*rg3 &
            -fpi*sum1c
       f3 = f3 &
            - tpi*(-1.0d0/z1)*2*z0*sum1c

       f4 = -fpi*z_l*rg3 &
            -fpi*sum1c
       f4 = f4 &
            - tpi*(-1.0d0/z1)*2*z0*sum1c
       !
       a0 = (f1*z_l**2*(z_l-3.d0*z_r)+z_r*(f3*z_l**2*(-z_l+z_r)    &
            +z_r*(f2*(3.d0*z_l-z_r)+f4*z_l*(-z_l+z_r))))/(z_l-z_r)**3
       a1 = (f3*z_l**3+z_l*(6.d0*f1-6.d0*f2+(f3+2.d0*f4)*z_l)*z_r  &
            -(2*f3+f4)*z_l*z_r**2-f4*z_r**3)/(z_l-z_r)**3
       a2 = (-3*f1*(z_l+z_r)+3.d0*f2*(z_l+z_r)-(z_l-z_r)*(2*f3*z_l &
            +f4*z_l+f3*z_r+2*f4*z_r))/(z_l-z_r)**3
       a3 = (2.d0*f1-2.d0*f2+(f3+f4)*(z_l-z_r))/(z_l-z_r)**3
       !
       ! remove polynomial from V(z)
       DO iz=1, dfftp%nr3
          IF ( iz-1 > dfftp%nr3/2 ) THEN
             z = DBLE(iz-1-dfftp%nr3)/DBLE(dfftp%nr3)*L
          ELSE
             z = DBLE(iz-1)/DBLE(dfftp%nr3)*L
          ENDIF
          Vr(iz) = Vr(iz) - (a0+a1*z+a2*z**2+a3*z**3)
       ENDDO
       !
       ! convert V(z) to V(gz) without polynomial
       CALL cft_1z( Vr, 1, dfftp%nr3, dfftp%nr3, -1, Vg )
       !
       ! add polynomial to V(gz)
       DO igz = 2, dfftp%nr3
          IF ( igz<=dfftp%nr3/2 )THEN
             gz = DBLE(igz-1)*tpi/L
          ELSE
             gz = DBLE(igz-1-dfftp%nr3)*tpi/L
          ENDIF
          !
          Vg(igz) = Vg(igz) &
               + a1 * ci * COS(gz*z0)/gz &
               + a2 * 2.0d0 * COS(gz*z0)/gz**2   &
               + a3 * ci * z0**2 * COS(gz*z0)/gz &
               - a3 * ci * 6.0d0 * COS(gz*z0)/gz**3
       ENDDO
       !
       Vg(1) = Vg(1) + a0*1.0d0 + a2*z0**2/3.0d0
       !
       ! add bare coulomn terms to V(gz)
       DO igz=2, dfftp%nr3
          IF ( igz<=dfftp%nr3/2 ) THEN
             gz = DBLE(igz-1)*tpi/L
          ELSE
             gz = DBLE(igz-1-dfftp%nr3)*tpi/L
          ENDIF
          !
          rg3 = rhog3(igz,imill_2d(0,0))
          Vg(igz) = Vg(igz) + fpi*rg3/gz**2
       ENDDO ! igz
       !
       ! calculate dV/deps(gz)
       DO igz = 1, dfftp%nr3
          dVg_deps(igz,:,:) = -delta(:,:) * Vg(igz)
       ENDDO ! igz
       !
       ! calculate stress tensor
       DO igz = 1, dfftp%nr3
          rg3 = rhog3(igz,imill_2d(0,0))
          sigmahar(1:2,1:2) = sigmahar(1:2,1:2) &
               + DBLE( CONJG(rg3) * dVg_deps(igz,1:2,1:2) )
       ENDDO ! igz
    ENDIF ! imill_2d(0,0) > 0
    !
    ! half means removing duplications.
    ! e2 means hartree -> Ry.
    sigmahar(:,:) = sigmahar(:,:) * (-0.5d0*e2)
    !
    CALL mp_sum( sigmahar, intra_bgrp_comm )
    !
    DEALLOCATE( rhog3 )
    DEALLOCATE( dVr_deps )
    DEALLOCATE( dVg_deps )
    DEALLOCATE( Vr )
    DEALLOCATE( Vg )
    !
    RETURN
    !
  END SUBROUTINE esm_stres_har_bc2
  !
  !
  !------------------------------------------------------------------
  SUBROUTINE esm_stres_har_bc3( sigmahar, rhog )
    !----------------------------------------------------------------
    !! Calculation of Hartree stress tensor - bc2.
    !
    USE kinds,         ONLY : DP
    USE gvect,         ONLY : ngm, mill
    USE constants,     ONLY : tpi, fpi, e2
    USE cell_base,     ONLY : omega, alat, at, tpiba, bg
    USE control_flags, ONLY : gamma_only
    USE fft_base,      ONLY : dfftp
    USE fft_scalar,    ONLY : cft_1z
    USE mp_bands,      ONLY : intra_bgrp_comm
    USE mp,            ONLY : mp_sum
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(OUT) :: sigmahar(3,3)
    !!  the Hartree term of the stress tensor
    COMPLEX(DP), INTENT(IN) :: rhog(ngm)
    !! n(G)
    !
    ! ... local variables
    !
    INTEGER :: ig, iga, igb, igz, igp, la, mu, iz
    REAL(DP) :: L, S, z0, z1, z
    REAL(DP) :: g(2), gp, gz
    COMPLEX(DP), PARAMETER :: ci = DCMPLX(0.0d0, 1.0d0)
    COMPLEX(DP) :: rg3
    COMPLEX(DP) :: sum1p, sum1m, sum2p, sum2m, sum1c, sum2c
    COMPLEX(DP) :: sum1sh, sum1ch, sum2sh
    REAL(DP)    :: z_l, z_r
    COMPLEX(DP) :: f1, f2, f3, f4, a0, a1, a2, a3
    COMPLEX(DP) :: poly_fr, poly_fl, poly_dfr, poly_dfl
    COMPLEX(DP) :: poly_a, poly_b, poly_c, poly_d
    REAL(DP), PARAMETER :: delta(2,2) = RESHAPE( (/ 1.0d0, 0.0d0, 0.0d0, 1.0d0 /), (/2,2/) )
    REAL(DP) :: dgp_deps(2,2)     ! dgp/deps
    REAL(DP) :: dgp2_deps(2,2)    ! dgp^2/deps
    REAL(DP) :: dinvgp_deps(2,2)  ! dgp^-1/deps
    !
    COMPLEX(DP), ALLOCATABLE :: rhog3(:,:)
    COMPLEX(DP), ALLOCATABLE :: dVr_deps(:,:,:)
    COMPLEX(DP), ALLOCATABLE :: dVg_deps(:,:,:)
    COMPLEX(DP), ALLOCATABLE :: Vr(:)
    COMPLEX(DP), ALLOCATABLE :: Vg(:)
    !
    REAL(DP) :: sigmahar_bc1(3,3)
    !
    ! cell settings
    L  = at(3,3)*alat
    S  = omega/L
    z0 = L/2.d0
    z1 = z0+esm_w
    !
    ! initialize
    sigmahar(:,:) = 0.0d0
    !
    ALLOCATE( rhog3(dfftp%nr3,ngm_2d) )
    ALLOCATE( dVr_deps(dfftp%nr3,2,2) )
    ALLOCATE( dVg_deps(dfftp%nr3,2,2) )
    ALLOCATE( Vr(dfftp%nr3) )
    ALLOCATE( Vg(dfftp%nr3) )
    !
    ! reconstruct rho(gz,gp)
    rhog3(:,:) = (0.d0,0.d0)
    !
    DO ig = 1, ngm
       iga = mill(1,ig)
       igb = mill(2,ig)
       igz = mill(3,ig)+1
       igp = imill_2d(iga,igb)
       IF ( igz<1 ) THEN
          igz = igz + dfftp%nr3
       ENDIF
       rg3 = rhog(ig)
       rhog3(igz,igp) = rg3
       !
       ! expand function symmetrically to gz<0
       IF ( gamma_only .AND. iga==0 .AND. igb==0 ) THEN
          igz = 1-mill(3,ig)
          IF ( igz<1 ) THEN
             igz = igz + dfftp%nr3
          ENDIF
          rhog3(igz,igp) = CONJG(rg3)
       ENDIF
    ENDDO ! ig
    !
    !****For gp!=0 case ****
    DO igp = 1, ngm_2d
       iga = mill_2d(1,igp)
       igb = mill_2d(2,igp)
       g(1:2) = (iga*bg(1:2,1) + igb*bg(1:2,2))*tpiba
       gp = SQRT(g(1)*g(1) + g(2)*g(2))
       !
       IF ( gp==0.0d0 ) CYCLE ! skip gp=0

       ! derivatives by strain tensor
       DO la=1, 2
          DO mu=1, 2
             dgp_deps(la,mu)    = -g(la)*g(mu)/gp
             dgp2_deps(la,mu)   = -g(la)*g(mu)*2.0d0
             dinvgp_deps(la,mu) = +g(la)*g(mu)/gp**3
          ENDDO
       ENDDO
       !
       ! summations over gz
       sum1p = (0.d0,0.d0)
       sum1m = (0.d0,0.d0)
       sum2p = (0.d0,0.d0)
       sum2m = (0.d0,0.d0)
       sum1sh = (0.d0,0.d0)
       sum1ch = (0.d0,0.d0)
       sum2sh = (0.d0,0.d0)
       !
       DO igz = 1, dfftp%nr3
          IF ( igz<=dfftp%nr3/2 ) THEN
             gz = DBLE(igz-1)*tpi/L
          ELSE
             gz = DBLE(igz-1-dfftp%nr3)*tpi/L
          ENDIF
          !
          rg3 = rhog3(igz,igp)
          sum1p = sum1p + rg3*qe_EXP(+ci*gz*z0)/(gp-ci*gz)
          sum1m = sum1m + rg3*qe_EXP(-ci*gz*z0)/(gp+ci*gz)
          sum2p = sum2p + rg3*qe_EXP(+ci*gz*z0)/(gp-ci*gz)**2
          sum2m = sum2m + rg3*qe_EXP(-ci*gz*z0)/(gp+ci*gz)**2
          sum1sh = sum1sh + rg3*qe_sinh(gp*z0+ci*gz*z0)/(gp+ci*gz)
          sum1ch = sum1ch + rg3*qe_cosh(gp*z0+ci*gz*z0)/(gp+ci*gz)*z0
          sum2sh = sum2sh + rg3*qe_sinh(gp*z0+ci*gz*z0)/(gp+ci*gz)**2
       ENDDO ! igz
       !
       ! calculate dV(z)/deps
       DO iz = 1, dfftp%nr3
          IF ( iz-1 > dfftp%nr3/2 ) THEN
             z = DBLE(iz-1-dfftp%nr3)/DBLE(dfftp%nr3)*L
          ELSE
             z = DBLE(iz-1)/DBLE(dfftp%nr3)*L
          ENDIF
          !
          ! BC1 terms
          dVr_deps(iz,:,:) = &
               -( dgp_deps(:,:)*tpi/gp**2*(gp*(z-z0)-1.0d0) &
               -  delta(:,:)*tpi/gp ) &
               * EXP(+gp*(z-z0)) * sum1p &
               + dgp_deps(:,:)*tpi/gp * EXP(+gp*(z-z0)) * sum2p &
               +( dgp_deps(:,:)*tpi/gp**2*(gp*(z+z0)+1.0d0) &
               +  delta(:,:)*tpi/gp ) &
               * EXP(-gp*(z+z0)) * sum1m &
               + dgp_deps(:,:)*tpi/gp * EXP(-gp*(z+z0)) * sum2m
          !
          ! BC3 termn
          dVr_deps(iz,:,:) = dVr_deps(iz,:,:) &
               - dgp_deps(:,:) * ( &
               -   fpi/gp**2 *EXP(-gp*(-z+2*z1)) &
               -   fpi/gp*(-z+2*z1)*EXP(-gp*(-z+2*z1)) &
               ) * sum1sh &
               - fpi/gp*EXP(-gp*(-z+2*z1))* ( &
               -   delta(:,:)* sum1sh &
               +   dgp_deps(:,:) * (sum1ch - sum2sh) )
       ENDDO ! iz
       !
       ! convert dV(z)/deps to dV(gz)/deps
       DO la = 1, 2
          DO mu = 1, 2
             CALL cft_1z( dVr_deps(:,la,mu), 1, dfftp%nr3, dfftp%nr3, -1, dVg_deps(:,la,mu) )
          ENDDO
       ENDDO
       !
       ! add bare coulomn terms to dV(gz)/deps
       DO igz = 1, dfftp%nr3
          IF ( igz<=dfftp%nr3/2 ) THEN
             gz = DBLE(igz-1)*tpi/L
          ELSE
             gz = DBLE(igz-1-dfftp%nr3)*tpi/L
          ENDIF
          rg3 = rhog3(igz,igp)
          !
          dVg_deps(igz,:,:) = dVg_deps(igz,:,:) &
               - delta(:,:)* fpi*rg3/(gp**2+gz**2) &
               - dgp2_deps(:,:)* fpi*rg3/(gp**2+gz**2)**2
       ENDDO ! igz
       !
       ! modifications
       IF ( gamma_only ) THEN
          dVg_deps(:,:,:) = dVg_deps(:,:,:)*2.0d0
       ENDIF
       !
       ! calculate stress tensor
       DO igz = 1, dfftp%nr3
          rg3 = rhog3(igz,igp)
          sigmahar(1:2,1:2) = sigmahar(1:2,1:2) &
               + DBLE( CONJG(rg3) * dVg_deps(igz,:,:) )
       ENDDO ! igz
       !
    ENDDO ! igp
    !
    !****For gp=0 case ****
    IF ( imill_2d(0,0) > 0 ) THEN
       ! summations over gz
       sum1c = (0.d0,0.d0)
       sum2c = (0.d0,0.d0)
       !
       DO igz = 2, dfftp%nr3
          IF ( igz<=dfftp%nr3/2 ) THEN
             gz = DBLE(igz-1)*tpi/L
          ELSE
             gz = DBLE(igz-1-dfftp%nr3)*tpi/L
          ENDIF
          !
          rg3 = rhog3(igz,imill_2d(0,0))
          sum1c = sum1c + rg3*ci*COS(gz*z0)/gz
          sum2c = sum2c + rg3*COS(gz*z0)/gz**2
       ENDDO ! igz
       !
       ! calculate V(z)
       DO iz = 1, dfftp%nr3
          IF ( iz-1 > dfftp%nr3/2 ) THEN
             z = DBLE(iz-1-dfftp%nr3)/DBLE(dfftp%nr3)*L
          ELSE
             z = DBLE(iz-1)/DBLE(dfftp%nr3)*L
          ENDIF
          !
          rg3 = rhog3(1,imill_2d(0,0))
          !! BC1 terms
          Vr(iz) = &
               - tpi*z**2*rg3  &
               - tpi*z0**2*rg3 &
               - fpi*z*sum1c   &
               - fpi*sum2c
          !
          ! BC3 terms
          Vr(iz) = Vr(iz) - tpi*(z-2*z1)*2*z0 * rg3 + fpi*z0*sum1c
       ENDDO ! iz
       !
       ! separation by polynomial
       z_l = -z0
       z_r = +z0
       !
       f1  = -tpi*z_r**2*rg3 &
             -tpi*z0**2*rg3  &
             -fpi*z_r*sum1c  &
             -fpi*sum2c
       f1 =  f1 &
             - tpi*(z_r-2*z1)*2*z0 * rg3 + fpi*z0*sum1c
       f2 = -tpi*z_l**2*rg3 &
            -tpi*z0**2*rg3  &
            -fpi*z_l*sum1c  &
            -fpi*sum2c
       f2 = f2 &
            - tpi*(z_l-2*z1)*2*z0 * rg3 + fpi*z0*sum1c
       f3 = -fpi*z_r*rg3 &
            -fpi*sum1c
       f3 = f3 &
            - tpi*(1.0d0)*2*z0 * rg3
       f4 = -fpi*z_l*rg3 &
            -fpi*sum1c
       f4 = f4 &
            - tpi*(1.0d0)*2*z0 * rg3
       !
       a0 = (f1*z_l**2*(z_l-3.d0*z_r)+z_r*(f3*z_l**2*(-z_l+z_r)    &
             +z_r*(f2*(3.d0*z_l-z_r)+f4*z_l*(-z_l+z_r))))/(z_l-z_r)**3
       a1 = (f3*z_l**3+z_l*(6.d0*f1-6.d0*f2+(f3+2.d0*f4)*z_l)*z_r  &
             -(2*f3+f4)*z_l*z_r**2-f4*z_r**3)/(z_l-z_r)**3
       a2 = (-3*f1*(z_l+z_r)+3.d0*f2*(z_l+z_r)-(z_l-z_r)*(2*f3*z_l &
             +f4*z_l+f3*z_r+2*f4*z_r))/(z_l-z_r)**3
       a3 = (2.d0*f1-2.d0*f2+(f3+f4)*(z_l-z_r))/(z_l-z_r)**3
       !
       ! remove polynomial from V(z)
       DO iz = 1, dfftp%nr3
          IF ( iz-1 > dfftp%nr3/2 ) THEN
             z = DBLE(iz-1-dfftp%nr3)/DBLE(dfftp%nr3)*L
          ELSE
             z = DBLE(iz-1)/DBLE(dfftp%nr3)*L
          ENDIF
          Vr(iz) = Vr(iz) - (a0+a1*z+a2*z**2+a3*z**3)
       ENDDO
       !
       ! convert V(z) to V(gz) without polynomial
       CALL cft_1z( Vr, 1, dfftp%nr3, dfftp%nr3, -1, Vg )
       !
       ! add polynomial to V(gz)
       DO igz = 2, dfftp%nr3
          IF ( igz<=dfftp%nr3/2 ) THEN
             gz = DBLE(igz-1)*tpi/L
          ELSE
             gz = DBLE(igz-1-dfftp%nr3)*tpi/L
          ENDIF
          !
          Vg(igz) = Vg(igz) &
               + a1 * ci * COS(gz*z0)/gz &
               + a2 * 2.0d0 * COS(gz*z0)/gz**2 &
               + a3 * ci * z0**2 * COS(gz*z0)/gz &
               - a3 * ci * 6.0d0 * COS(gz*z0)/gz**3
       ENDDO
       !
       Vg(1) = Vg(1) + a0*1.0d0 + a2*z0**2/3.0d0
       !
       ! add bare coulomn terms to V(gz)
       DO igz = 2, dfftp%nr3
          IF ( igz<=dfftp%nr3/2 ) THEN
             gz = DBLE(igz-1)*tpi/L
          ELSE
             gz = DBLE(igz-1-dfftp%nr3)*tpi/L
          ENDIF
          !
          rg3 = rhog3(igz,imill_2d(0,0))
          Vg(igz) = Vg(igz) + fpi*rg3/gz**2
       ENDDO ! igz
       !
       ! calculate dV/deps(gz)
       DO igz = 1, dfftp%nr3
          dVg_deps(igz,:,:) = -delta(:,:) * Vg(igz)
       ENDDO ! igz
       !
       ! calculate stress tensor
       DO igz = 1, dfftp%nr3
          rg3 = rhog3(igz,imill_2d(0,0))
          !
          sigmahar(1:2,1:2) = sigmahar(1:2,1:2) &
               + DBLE( CONJG(rg3) * dVg_deps(igz,1:2,1:2) )
       ENDDO ! igz
       !
    ENDIF ! imill_2d(0,0) > 0
    !
    ! half means removing duplications.
    ! e2 means hartree -> Ry.
    sigmahar(:,:) = sigmahar(:,:) * (-0.5d0*e2)
    !
    CALL mp_sum( sigmahar, intra_bgrp_comm )
    !
    DEALLOCATE( rhog3 )
    DEALLOCATE( dVr_deps )
    DEALLOCATE( dVg_deps )
    DEALLOCATE( Vr )
    DEALLOCATE( Vg )
    !
    RETURN
    !
  END SUBROUTINE esm_stres_har_bc3
  !
  !
  !----------------------------------------------------------------------
  FUNCTION esm_ewald()
    !---------------------------------------------------------------------
    !! Calculates Ewald energy with both G- and R-space terms.  
    !! Determines optimal alpha. Should hopefully work for any structure.
    !
    USE kinds,       ONLY : DP
    USE constants,   ONLY : tpi, e2
    USE mp_bands,    ONLY : intra_bgrp_comm
    USE mp,          ONLY : mp_sum
    USE cell_base,   ONLY : tpiba2
    USE ions_base,   ONLY : zv, nat, ityp
    USE gvect,       ONLY : gcutm
    !
    IMPLICIT NONE
    !
    REAL(DP) :: esm_ewald
    !! output: the ewald energy
    !
    ! ... local variables
    !
    INTEGER :: na
    ! counter on atoms
    !
    REAL(DP) :: charge, ewaldg, ewaldr, alpha, upperbound
    ! total ionic charge in the cell
    ! ewald energy computed in reciprocal space
    ! ewald energy computed in real space
    ! alpha term in ewald sum
    ! the maximum radius to consider real space sum
    !
    charge = 0.d0
    DO na = 1, nat
       charge = charge + zv( ityp(na) )
    ENDDO
    !
    ! choose alpha in order to have convergence in the sum over G
    ! upperbound is a safe upper bound for the error in the sum over G
    alpha = 2.9d0
    DO
       alpha = alpha - 0.1d0
       IF (alpha <= 0.d0) CALL errore( 'esm_ewald', 'optimal alpha not found', 1 )
       upperbound = 2.d0 * charge**2 * SQRT(2.d0 * alpha / tpi) * &
                    qe_erfc ( SQRT(tpiba2 * gcutm / 4.d0 / alpha) )
       IF ( upperbound < 1.0d-7 ) EXIT
    ENDDO
    !
    ! G-space sum here.
    ! Determine IF this processor contains G=0 and set the constant term
    CALL esm_ewaldg( alpha, ewaldg )
    !
    ! R-space sum here (only for the processor that contains G=0)
    CALL esm_ewaldr( alpha, ewaldr )
    !
    esm_ewald = 0.5d0 * e2 * ( ewaldg + ewaldr )
    !
    CALL mp_sum( esm_ewald, intra_bgrp_comm )
    !write( *,'(5x,"alpha used in ewald term: ",f5.2 )')alpha
    !
    RETURN
    !
  END FUNCTION esm_ewald
  !
  !
  !  ... ESM EWALD RSUM SUBROUTINE
  !
  !-----------------------------------------------------------------------
  SUBROUTINE esm_ewaldr_pbc( alpha_g, ewr )
    !----------------------------------------------------------------------
    !! Ewald R-space sum - pbc.
    !
    USE io_global,        ONLY : stdout
    USE constants,        ONLY : pi, tpi, fpi, e2
    USE gvect,            ONLY : gstart
    USE cell_base,        ONLY : alat, tpiba2, at, bg
    USE ions_base,        ONLY : zv, nat, tau, ityp
    USE control_flags,    ONLY : iverbosity
    USE mp,               ONLY : mp_rank, mp_size
    USE mp_bands,         ONLY : intra_bgrp_comm
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN)  :: alpha_g
    !! alpha term in ewald sum
    REAL(DP), INTENT(OUT) :: ewr
    !! R-space term for Ewald energy
    !
    ! ... local variables
    !
    INTEGER :: na, nb, nr, nrm, np, ip
    ! counter on atoms
    ! counter on atoms
    ! counter over direct vectors
    ! number of R vectors included in r sum
    INTEGER, PARAMETER :: mxr = 500
    ! the maximum number of R vectors included in r
    REAL(DP) :: dtau(3), r(3,mxr), r2(mxr)
    ! the difference tau_s - tau_s'
    ! neighbering shell vector
    ! the square modulus of R_j-tau_s-tau_s'
    ! buffer variable
    ! buffer variable
    !
    ! ESM variables
    REAL(DP) :: tmp, fac, ss, ew, rmax0, rr
    !
    ewr = 0.d0
    !
    tmp = SQRT(alpha_g)
    rmax0 = 4.d0/tmp/alat
    !
    ip = mp_rank( intra_bgrp_comm )
    np = mp_size( intra_bgrp_comm )
    !
    ew = 0.d0
    !
    DO na = ip+1, nat, np
       DO nb = 1, nat
          !
          dtau(:) = tau(:,na)-tau(:,nb)
          fac = zv(ityp(nb))*zv(ityp(na))
          !
          ! generates nearest-neighbors shells
          CALL rgen( dtau, rmax0, mxr, at, bg, r, r2, nrm )
          !
          ! and sum to the real space part
          DO nr = 1, nrm
             rr = SQRT(r2(nr))*alat
             ew = ew + fac*qe_erfc(tmp*rr)/rr
          ENDDO
          !
       ENDDO
       ! Here add the other constant term
       ew = ew-zv(ityp(na))**2*tmp/SQRT(pi)*2.d0 ! 2.d0: fit to original code
    ENDDO
    !
    ewr = ewr + ew
    !
  END SUBROUTINE esm_ewaldr_pbc
  !
  !
  !------------------------------------------------------------------------------
  SUBROUTINE esm_ewaldr_bc4( alpha_g, ewr )
    !----------------------------------------------------------------------------
    !! Ewald R-space sum - bc4.
    !
    USE io_global,        ONLY : stdout
    USE constants,        ONLY : pi, tpi, fpi, e2
    USE gvect,            ONLY : gstart
    USE cell_base,        ONLY : alat, tpiba2, at, bg
    USE ions_base,        ONLY : zv, nat, tau, ityp
    USE control_flags,    ONLY : iverbosity
    USE mp,               ONLY : mp_rank, mp_size
    USE mp_bands,         ONLY : intra_bgrp_comm
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN)  :: alpha_g
    !! alpha term in ewald sum
    REAL(DP), INTENT(OUT) :: ewr
    !! R-space term for Ewald energy
    !
    ! ... local variables
    !
    INTEGER :: na, nb, nr, nrm, np, ip
    ! counter on atoms
    ! counter on atoms
    ! counter over direct vectors
    ! number of R vectors included in r sum
    INTEGER, PARAMETER :: mxr = 500
    ! the maximum number of R vectors included in r
    REAL(DP) :: dtau(3), r(3,mxr), r2(mxr), rxy, rxyz
    ! the difference tau_s - tau_s'
    ! neighbering shell vector
    ! the square modulus of R_j-tau_s-tau_s'
    ! buffer variable
    ! buffer variable
    !
    ! ESM variables
    REAL(DP) :: L, z, zp, z0, z1, aaa, tmp, fac, ss, ew, err, ss0, &
                gpmax, rmax0, rmax, zbuff, znrm, rr
    ! gpmax: upper bound of g_parallel integral
    ! rmax: the maximum radius to consider real space sum
    ! zbuff: smearing width to avoid the singularity of the Force
    ! znrm: threashold value for normal RSUM and Smooth-ESM's RSUM
    REAL(DP), PARAMETER :: eps=1.d-11, epsneib=1.d-6
    !
    ewr = 0.d0
    L = at(3,3)*alat
    z0 = L/2.d0
    z1 = z0 + esm_w
    aaa = esm_a
    tmp = SQRT(alpha_g)
    zbuff = 1.d0
    !
    ! Define upperbound for g_parallel integral
    err=1.d0; ss0=0.d0; gpmax=1.d0
    !
    DO
       gpmax = gpmax+1.d0
       IF (gpmax > 1000.d0) &
          CALL errore( 'esm_ewaldr', 'optimal gpmax not found', 1 )
       CALL qromb( vl11, aaa, tmp, z1, z1-zbuff, z1-zbuff, 0.0_DP, gpmax, ss )
       err = ABS(ss-ss0); ss0=ss
       IF (err < eps) EXIT
    ENDDO
    !
    ! Define znrm using the deviation from the constant term in RSUM
    znrm = z1
    DO
       znrm = znrm-0.01d0
       IF ( znrm <= -z0 ) &
          CALL errore( 'esm_ewaldr', 'optimal znrm not found', 1 )
       CALL qromb( vl11, aaa, tmp, z1, znrm, znrm, 0.0_DP, gpmax, ss )
       err = -2.d0*tmp/SQRT(pi)-ss*2.d0
       IF (ABS(err) < eps) EXIT
    ENDDO
    ! Define rmax for real space sum
    rmax = 1.d0
    !
    DO
       rmax = rmax+1.d0
       IF (rmax > 200.d0) &
          CALL errore ('esm_ewaldr', 'optimal rmax not found', 1)
       CALL qromb( vl11j0, aaa, tmp, z1, z1-zbuff, z1-zbuff, rmax, gpmax, ss )
       err = 1.d0/rmax+ss*2.d0
       IF (ABS(err) < epsneib) EXIT
    ENDDO
    !
    rmax = rmax/alat
    !
    IF (iverbosity > 0) THEN
       WRITE( stdout, '(5x,"=== Smooth-ESM RSUM parameters (Energy) ===")')
       WRITE( stdout, '(5x,A,F10.2,A)') &
          'Upper bound of g_parallel integral:      ',gpmax,' (1/a.u.)'
       WRITE( stdout, '(5x,A,F10.2,A)') &
          'Boundary for normal RSUM|Smooth-ESM RSUM:',z1-znrm,' (a.u.)'
       WRITE( stdout, '(5x,A,F10.2,A)') &
          'Upper bound of real-space summation:     ',rmax*alat,' (a.u.)'
       WRITE( stdout, '(5x,"===========================================")')
    ENDIF
    !
    ip = mp_rank( intra_bgrp_comm )
    np = mp_size( intra_bgrp_comm )
    !
    ew = 0.d0
    !
    DO na = ip+1, nat, np
       z = tau(3,na)
       IF (z > at(3,3)*0.5) z = z-at(3,3)
       z = z*alat
       !
       DO nb = 1, nat
          zp = tau(3,nb)
          IF (zp>at(3,3)*0.5) zp = zp-at(3,3)
          zp = zp*alat
          dtau(1:2) = tau(1:2,na)-tau(1:2,nb)
          dtau(3) = (z-zp)/alat
          fac = zv(ityp(nb))*zv(ityp(na))
          !
          IF ( z < znrm ) THEN
             !
             IF ( zp < znrm ) THEN ! z in I, zp in I (normal RSUM)
                rmax0 = 4.d0/tmp/alat
                !
                ! generates nearest-neighbors shells
                !
                CALL rgen( dtau, rmax0, mxr, at, bg, r, r2, nrm )
                !
                ! and sum to the real space part
                !
                DO nr = 1, nrm
                   rr = SQRT(r2(nr))*alat
                   ew = ew + fac*qe_erfc(tmp*rr)/rr
                ENDDO
                !
             ELSEIF ( zp < z1 ) THEN ! z in I, zp in I
                !
                CALL esm_rgen_2d( dtau, rmax, mxr, at, bg, r, r2, nrm )
                !
                DO nr = 1, nrm
                   rxy = SQRT(r2(nr))*alat
                   rxyz = SQRT(r2(nr)+dtau(3)**2)*alat
                   CALL qromb( vl11j0, aaa, tmp, z1, z, zp, rxy, gpmax, ss )
                   ew = ew + fac*(1.d0/rxyz+ss*2.d0)
                ENDDO
             ELSE ! z in I, zp in II
                CALL esm_rgen_2d( dtau, rmax, mxr, at, bg, r, r2, nrm )
                DO nr = 1, nrm
                   rxy = SQRT(r2(nr))*alat
                   CALL qromb( vl12j0, aaa, tmp, z1, z, zp, rxy, gpmax, ss )
                   ew = ew + fac*ss*2.d0
                ENDDO
             ENDIF ! IF for zp
             !
          ELSEIF ( z < z1 ) THEN ! znrm < z < z1
             !
             CALL esm_rgen_2d( dtau, rmax, mxr, at, bg, r, r2, nrm )
             !
             IF ( zp < z1 ) THEN ! z in I, zp in I
                DO nr = 1, nrm
                   rxy = SQRT(r2(nr))*alat
                   rxyz = SQRT(r2(nr)+dtau(3)**2)*alat
                   CALL qromb( vl11j0, aaa, tmp, z1, z, zp, rxy, gpmax, ss )
                   ew = ew + fac*(1.d0/rxyz+ss*2.d0)
                ENDDO
             ELSE ! z in I, zp in II
                DO nr = 1, nrm
                   rxy = SQRT(r2(nr))*alat
                   CALL qromb(vl12j0,aaa,tmp,z1,z,zp,rxy,gpmax,ss)
                   ew=ew+fac*ss*2.d0
                ENDDO
             ENDIF ! IF for zp
             !
          ELSE ! z1 < z
             !
             CALL esm_rgen_2d( dtau, rmax, mxr, at, bg, r, r2, nrm )
             !
             IF ( zp < z1 ) THEN ! z in II, zp in I
                DO nr = 1, nrm
                   rxy = SQRT(r2(nr))*alat
                   CALL qromb( vl21j0, aaa, tmp, z1, z, zp, rxy, gpmax, ss )
                   ew = ew + fac*ss*2.d0
                ENDDO
             ELSE ! z in II, zp in II
                DO nr = 1, nrm
                   rxy = SQRT(r2(nr))*alat
                   rxyz = SQRT(r2(nr)+dtau(3)**2)*alat
                   CALL qromb( vl22j0, aaa, tmp, z1, z, zp, rxy, gpmax, ss )
                   ew = ew + fac*(EXP(-aaa*(rxyz+z+zp-2.d0*z1))/rxyz+ss*2.d0)
                ENDDO
             ENDIF
             !
          ENDIF ! IF for z
          !
       ENDDO
       !
       IF (z < znrm ) THEN
          ss = -tmp/SQRT(pi)
       ELSEIF (z < z1) THEN
          CALL qromb( vl11, aaa, tmp, z1, z, z, 0.0_DP, gpmax, ss )
       ELSE
          CALL qromb( vl22, aaa, tmp, z1, z, z, 0.0_DP, gpmax, ss )
       ENDIF
       !
       ew = ew + zv(ityp(na))**2*ss*2.d0 ! 2.0: fit to original code
    ENDDO
    !
    ewr = ewr + ew
    !
  END SUBROUTINE esm_ewaldr_bc4
  !
  !
  !------------------------------------------------------------------------------
  SUBROUTINE esm_stres_ewr( alpha, sigmaewa )
    !----------------------------------------------------------------------------
    !! Wrapper for the calculation of R-space Ewald term - pbc.
    !
    USE kinds,     ONLY : DP
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN)  :: alpha
    !! alpha term in ewald sum
    REAL(DP), INTENT(OUT) :: sigmaewa(3,3)
    !! the Ewald term of the stress tensor
    !
    SELECT CASE( esm_bc )
    CASE( 'pbc')
       STOP 'esm_stres_ewa must not be called for esm_bc = pbc'
    CASE( 'bc1' )
       CALL esm_stres_ewr_pbc( alpha, sigmaewa )
    CASE( 'bc2' )
       CALL esm_stres_ewr_pbc( alpha, sigmaewa )
    CASE( 'bc3' )
       CALL esm_stres_ewr_pbc( alpha, sigmaewa )
    CASE( 'bc4' )
       STOP 'esm_stres_ewa has not yet implemented for esm_bc = bc4'
    END SELECT
    !
    RETURN
    !
  END SUBROUTINE esm_stres_ewr
  !
  !
  !---------------------------------------------------------------------------------
  SUBROUTINE esm_stres_ewr_pbc( alpha, sigmaewa )
    !-------------------------------------------------------------------------------
    !! Calculation of R-space Ewald term - pbc.
    !
    USE kinds,     ONLY : DP
    USE constants, ONLY : pi, sqrtpm1, tpi, fpi, e2
    USE cell_base, ONLY : omega, alat, at, tpiba, bg
    USE ions_base, ONLY : zv, nat, tau, ityp
    USE gvect,     ONLY : gstart
    USE mp_bands,  ONLY : intra_bgrp_comm
    USE mp,        ONLY : mp_sum
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN)  :: alpha
    !! alpha term in ewald sum
    REAL(DP), INTENT(OUT) :: sigmaewa(3,3)
    !! the Ewald term of the stress tensor
    !
    ! ... local variables
    !
    INTEGER, PARAMETER :: mxr = 50
    ! the maximum number of R vectors included in r sum
    INTEGER :: ia, ib, nr, nrm, la, mu
    REAL(DP) :: Qa, Qb, dtau(3), rmax
    REAL(DP) :: salp, r(3,mxr), r2(mxr), rr, fac
    !
    salp = SQRT(alpha)
    !
    ! initialize
    sigmaewa(:,:) = 0.d0
    !
    ! R-space sum here (only for the processor that contains G=0)
    !
    IF ( gstart == 2 ) THEN
       rmax = 4.0d0/salp/alat
       !
       ! with this choice terms up to ZiZj*erfc(5) are counted (erfc(5)=2x10^-1
       !
       DO ib = 1, nat
          Qb = (-1.0d0)*zv(ityp(ib))
          DO ia = 1, nat
             Qa = (-1.0d0)*zv(ityp(ia))
             !
             !     generates nearest-neighbors shells r(i)=R(i)-dtau(i)
             !
             dtau(:) = tau(:,ib) - tau(:,ia)
             CALL rgen( dtau, rmax, mxr, at, bg, r, r2, nrm )
             !
             DO nr = 1, nrm
                rr = SQRT(r2(nr))*alat
                r(:,nr) = r(:,nr)*alat

                fac = Qb*Qa/rr**3 &
                     * ( qe_erfc(salp*rr) &
                     + rr*2.0d0*salp*sqrtpm1 * EXP(-alpha*rr**2) )
                DO la=1, 3
                   DO mu=1, 3
                      sigmaewa(la,mu) = sigmaewa(la,mu) + fac*r(la,nr)*r(mu,nr)
                   ENDDO ! mu
                ENDDO ! la
                !
             ENDDO ! nr
             !
          ENDDO ! ia
       ENDDO ! ib
       !
    ENDIF
    !
    sigmaewa(:,:) = sigmaewa(:,:)*(e2/2.0d0/omega)
    !
    CALL mp_sum( sigmaewa, intra_bgrp_comm )
    !
    RETURN
    !
  END SUBROUTINE esm_stres_ewr_pbc
  !
  !
  ! ... ESM EWALD GSUM SUBROUTINE
  !
  !-----------------------------------------------------------------------
  SUBROUTINE esm_ewaldg_pbc( alpha_g, ewg )
    !---------------------------------------------------------------------
    !! Ewald G-space sum - pbc.
    !
    USE constants,        ONLY : tpi
    USE gvect,            ONLY : gstart
    USE cell_base,        ONLY : omega, tpiba2
    USE ions_base,        ONLY : zv, nat, nsp, ityp
    USE control_flags,    ONLY : gamma_only
    USE gvect,            ONLY : ngm, gg
    USE vlocal,           ONLY : strf
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN)  :: alpha_g
    !! alpha term in ewald sum
    REAL(DP), INTENT(OUT) :: ewg
    !! Ewald G-space sum
    !
    ! ... local variables
    !
    INTEGER :: ng
    REAL(DP) :: charge, fact
    COMPLEX(DP) :: rhon
    !
    charge = SUM( zv( ityp(1:nat) ) )
    !
    ! same of the GSUM part in ewald.f90
    IF ( gstart == 2 ) THEN
       ewg = - charge**2 / alpha_g / 4.0d0
    ELSE
       ewg = 0.0d0
    ENDIF
    !
    IF ( gamma_only ) THEN
       fact = 2.d0
    ELSE
       fact = 1.d0
    ENDIF
    !
    DO ng = gstart, ngm
       rhon = SUM( zv(1:nsp) * CONJG( strf (ng, 1:nsp) ) )
       ewg = ewg + fact * ABS(rhon)**2 * &
             EXP( -gg(ng)*tpiba2/alpha_g/4.d0 )/gg(ng)/tpiba2
    ENDDO
    !
    ewg = 2.d0 * tpi / omega * ewg
    !
  END SUBROUTINE esm_ewaldg_pbc
  !
  !
  !-----------------------------------------------------------------------
  SUBROUTINE esm_ewaldg_bc1( alpha_g, ewg )
    !-----------------------------------------------------------------------
    !! Ewald G-space sum - bc1.
    !
    USE constants,        ONLY : pi, tpi, fpi
    USE gvect,            ONLY : gstart
    USE cell_base,        ONLY : omega, alat, tpiba2, at, bg
    USE ions_base,        ONLY : zv, nat, nsp, ityp, tau
    USE control_flags,    ONLY : gamma_only
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN)  :: alpha_g
    !! alpha term in ewald sum
    REAL(DP), INTENT(OUT) :: ewg
    !! Ewald G-space sum
    !
    ! ... local variables
    !
    INTEGER :: k1, k2, it1, it2, ng_2d
    REAL(DP) :: gp2, t(2), gp, sa, z, zp, z0, L, t1, t2, tt,    &
                tmp, cc1, cc2, kk1, kk2, ff, ew,arg001, arg002, &
                arg101, arg102
    !
    ewg = 0.d0
    L = at(3,3)*alat
    z0 = L/2.d0
    tmp = SQRT(alpha_g)
    sa = omega/L
    ew = 0d0
    !
    DO it1 = 1, nat
      DO it2 = 1, nat
         !
         z = tau(3,it1)
         IF (z>at(3,3)*0.5) z = z-at(3,3)
         z = z*alat
         zp = tau(3,it2)
         IF (zp>at(3,3)*0.5) zp = zp-at(3,3)
         zp = zp*alat
         tt = zv(ityp(it1))*zv(ityp(it2))*fpi/sa
         ! bc1 
         arg001=-tmp**2*(z-zp)**2
         arg101= tmp*(z-zp)
         kk1=0.5d0*(-(z-zp)*qe_erf(arg101)-EXP(arg001)/tmp/SQRT(pi))
         kk2=0.d0
         !
         cc1=0.d0
         cc2=0.d0
         !
         DO ng_2d = 1, ngm_2d
            k1 = mill_2d(1,ng_2d)
            k2 = mill_2d(2,ng_2d)
            !
            IF ( k1==0 .AND. k2==0 ) CYCLE
            t(1:2) = k1*bg(1:2, 1)+k2*bg(1:2, 2)
            !
            gp2 = SUM(t(:)*t(:))*tpiba2
            gp = SQRT(gp2)
            !
            ff = ((k1*bg(1,1)+k2*bg(1,2))*(tau(1,it1)-tau(1,it2))  &
                 +(k1*bg(2,1)+k2*bg(2,2))*(tau(2,it1)-tau(2,it2)))*tpi
            ! bc1
            arg001 =-gp*(z-zp)
            arg002 = gp*(z-zp)
            arg101 = gp/2.d0/tmp-tmp*(z-zp)
            arg102 = gp/2.d0/tmp+tmp*(z-zp)
            !
            t1 = exp_erfc(arg001,arg101)
            t2 = exp_erfc(arg002,arg102)
            !
            cc1 = cc1 + COS(ff)*(t1+t2)/4.d0/gp
            cc2 = 0.d0
         ENDDO
         !
         IF ( gamma_only ) THEN
            cc1=cc1*2d0
            cc2=cc2*2d0
         ENDIF
         !
         ew = ew + tt*(cc1+cc2)
         IF (gstart==2) ew = ew + tt*(kk1+kk2)
         !
      ENDDO
    ENDDO
    !
    ewg = ewg + ew
    !
    RETURN
    !
  END SUBROUTINE esm_ewaldg_bc1
  !
  !
  !-----------------------------------------------------------------------
  SUBROUTINE esm_ewaldg_bc2( alpha_g, ewg )
    !---------------------------------------------------------------------
    !! Ewald G-space sum - bc2.
    !
    USE constants,        ONLY : pi, tpi, fpi, e2
    USE gvect,            ONLY : gstart
    USE cell_base,        ONLY : omega, alat, tpiba2, at, bg
    USE ions_base,        ONLY : zv, nat, nsp, ityp, tau
    USE control_flags,    ONLY : gamma_only
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN)  :: alpha_g
    !! alpha term in ewald sum
    REAL(DP), INTENT(OUT) :: ewg
    !! Ewald G-space sum
    !
    ! ... local variables
    !
    INTEGER :: k1, k2, it1, it2, ng_2d
    REAL(DP) :: gp2, t(2), gp, sa, z, zp, z1, z0, L, t1, t2, tt, &
                tmp, cc1, cc2, kk1, kk2, ff, ew, arg001, arg002, &
                arg003, arg004, arg005, arg006, arg007, arg101,  &
                arg102
    !
    ewg = 0.d0
    L = at(3,3)*alat
    z0 = L/2.d0
    z1 = z0 + esm_w
    tmp = SQRT(alpha_g)
    sa = omega/L
    ew = 0d0
    !
    DO it1 = 1, nat
      DO it2 = 1, nat
         !
         z = tau(3,it1)
         IF (z > at(3,3)*0.5) z = z-at(3,3)
         z = z*alat
         zp = tau(3,it2)
         IF (zp > at(3,3)*0.5) zp = zp-at(3,3)
         zp = zp*alat
         tt = zv(ityp(it1))*zv(ityp(it2))*fpi/sa
         !
         IF (gstart==2) THEN
            IF (it1==it2) THEN
               ! add coulomb energy of ions under efield
               ew = ew - zv(ityp(it1))*(z1-z)*esm_efield/e2*2.0
            ENDIF
         ENDIF
         !
         ! bc2
         arg001 = -tmp**2*(z-zp)**2
         arg101 =  tmp*(z-zp)
         !
         kk1 = 0.5d0*(-(z-zp)*qe_erf(arg101)-EXP(arg001)/tmp/SQRT(pi))
         kk2 = 0.5d0*(z1-z*zp/z1)
         !
         cc1 = 0.d0
         cc2 = 0.d0
         !
         DO ng_2d = 1, ngm_2d
            k1 = mill_2d(1,ng_2d)
            k2 = mill_2d(2,ng_2d)
            IF ( k1==0 .AND. k2==0 ) CYCLE
            t(1:2) = k1*bg(1:2, 1)+k2*bg(1:2, 2)
            gp2 = SUM(t(:)*t(:))*tpiba2
            gp = SQRT(gp2)
            ff = ((k1*bg(1,1)+k2*bg(1,2))*(tau(1,it1)-tau(1,it2))  &
                 +(k1*bg(2,1)+k2*bg(2,2))*(tau(2,it1)-tau(2,it2)))*tpi
            ! bc2
            arg001 = -gp*(z-zp)
            arg002 =  gp*(z-zp)
            arg003 = -gp*(z+zp+2.d0*z1)
            arg004 =  gp*(z+zp-2.d0*z1)
            arg005 = -gp*(z-zp+4.d0*z1)
            arg006 =  gp*(z-zp-4.d0*z1)
            arg007 = -4.d0*gp*z1
            arg101 = gp/2.d0/tmp-tmp*(z-zp)
            arg102 = gp/2.d0/tmp+tmp*(z-zp)
            !
            t1=exp_erfc(arg001,arg101)
            t2=exp_erfc(arg002,arg102)
            !
            cc1 = cc1 + COS(ff)*(t1+t2)/4.d0/gp
            cc2 = cc2 + COS(ff)*(EXP(arg006)+EXP(arg005) &
                        -EXP(arg004)-EXP(arg003) ) &
                        /(1.d0-EXP(arg007))/2.d0/gp
         ENDDO
         !
         IF ( gamma_only ) THEN
            cc1=cc1*2d0
            cc2=cc2*2d0
         ENDIF
         !
         ew = ew + tt*(cc1+cc2)
         IF (gstart==2) ew = ew + tt*(kk1+kk2)
         !
      ENDDO
    ENDDO
    !
    ewg = ewg + ew
    !
    RETURN
    !
  END SUBROUTINE esm_ewaldg_bc2
  !
  !
  !-----------------------------------------------------------------------------
  SUBROUTINE esm_ewaldg_bc3( alpha_g, ewg )
    !---------------------------------------------------------------------------
    !! Ewald G-space sum - bc3.
    !
    USE constants,        ONLY : pi, tpi, fpi
    USE gvect,            ONLY : gstart
    USE cell_base,        ONLY : omega, alat, tpiba2, at, bg
    USE ions_base,        ONLY : zv, nat, nsp, ityp, tau
    USE control_flags,    ONLY : gamma_only
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN)  :: alpha_g
    !! alpha term in ewald sum
    REAL(DP), INTENT(OUT) :: ewg
    !! Ewald G-space sum
    !
    ! ... local variables
    !
    INTEGER :: k1, k2, it1, it2, ng_2d
    REAL(DP) :: gp2, t(2), gp, sa, z, zp, z1, z0, L, t1, t2, tt, &
                tmp, cc1, cc2, kk1, kk2, ff, ew, arg001, arg002, &
                arg003, arg101, arg102
    !
    ewg = 0.d0
    L = at(3,3)*alat
    z0 = L/2.d0
    z1 = z0+esm_w
    tmp = SQRT(alpha_g)
    sa = omega/L
    ew = 0d0
    !
    DO it1 = 1, nat
      DO it2 = 1, nat
         !
         z = tau(3,it1)
         IF (z > at(3,3)*0.5) z = z-at(3,3)
         z = z*alat
         zp = tau(3,it2)
         IF (zp > at(3,3)*0.5) zp = zp-at(3,3)
         zp = zp*alat
         tt = zv(ityp(it1))*zv(ityp(it2))*fpi/sa
         ! bc3
         arg001 = -tmp**2*(z-zp)**2
         arg101 =  tmp*(z-zp)
         !
         kk1 = 0.5d0*(-(z-zp)*qe_erf(arg101)-EXP(arg001)/tmp/SQRT(pi))
         kk2 = 0.5d0*(2.d0*z1-z-zp)
         !
         cc1 = 0.d0
         cc2 = 0.d0
         !
         DO ng_2d = 1, ngm_2d
            k1 = mill_2d(1,ng_2d)
            k2 = mill_2d(2,ng_2d)
            IF ( k1==0 .AND. k2==0 ) CYCLE
            t(1:2) = k1*bg(1:2, 1)+k2*bg(1:2,2)
            !
            gp2 = SUM(t(:)*t(:))*tpiba2
            gp = SQRT(gp2)
            !
            ff = ((k1*bg(1,1)+k2*bg(1,2))*(tau(1,it1)-tau(1,it2))  &
                 +(k1*bg(2,1)+k2*bg(2,2))*(tau(2,it1)-tau(2,it2)))*tpi
            ! bc3
            arg001 =-gp*(z-zp)
            arg002 = gp*(z-zp)
            arg003 = gp*(z+zp-2.d0*z1)
            arg101 = gp/2.d0/tmp-tmp*(z-zp)
            arg102 = gp/2.d0/tmp+tmp*(z-zp)
            !
            t1 = exp_erfc(arg001,arg101)
            t2 = exp_erfc(arg002,arg102)
            !
            cc1 = cc1 + COS(ff)*(t1+t2)/4.d0/gp
            cc2 = cc2 + COS(ff)*(-EXP(arg003))/2.d0/gp
         ENDDO
         !
         IF ( gamma_only ) THEN
            cc1 = cc1*2d0
            cc2 = cc2*2d0
         ENDIF
         !
         ew = ew + tt*(cc1+cc2)
         IF (gstart==2) ew = ew + tt*(kk1+kk2)
         !
      ENDDO
    ENDDO
    !
    ewg = ewg + ew
    !
    RETURN
    !
  END SUBROUTINE esm_ewaldg_bc3
  !
  !
  !------------------------------------------------------------------------------
  SUBROUTINE esm_ewaldg_bc4( alpha_g, ewg )
    !----------------------------------------------------------------------------
    !! Ewald G-space sum - bc3.
    !
    USE constants,        ONLY : pi, tpi, fpi
    USE gvect,            ONLY : gstart
    USE cell_base,        ONLY : omega, alat, tpiba2, at, bg
    USE ions_base,        ONLY : zv, nat, nsp, ityp, tau
    USE control_flags,    ONLY : gamma_only
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN)  :: alpha_g
    !! alpha term in ewald sum
    REAL(DP), INTENT(OUT) :: ewg
    !! Ewald G-space sum
    !
    ! ... local variables
    !
    INTEGER :: k1, k2, it1, it2, ng_2d
    REAL(DP) :: gp2, t(2), gp, sa, z, zp, z1, z0, L, t1, t2, tt, &
                tmp, cc1, cc2, kk1, kk2, ff, ew,arg001, arg002,  &
                arg003, arg005, arg006, arg007, arg008,          &
                arg009, arg011, arg101, arg102, arg103, arg104,  &
                arg106, arg107, arg109, arg111, arg113, aaa, t3, &
                alpha, beta, kappa, lambda, xi, chi
  
    ewg = 0.d0
    L = at(3,3)*alat
    z0 = L/2.d0
    z1 = z0 + esm_w
    aaa = esm_a
    tmp = SQRT(alpha_g)
    sa = omega/L
    ew = 0d0
    !
    DO it1 = 1, nat
      DO it2 = 1, nat
         !
         z = tau(3,it1)
         IF (z > at(3,3)*0.5) z = z-at(3,3)
         z = z*alat
         zp = tau(3,it2)
         IF (zp > at(3,3)*0.5) zp = zp-at(3,3)
         zp = zp*alat
         tt = zv(ityp(it1))*zv(ityp(it2))*fpi/sa
         !
         ! bc4
         arg001 = -tmp**2*(z-zp)**2
         arg002 = -tmp**2*(z1-zp)**2
         arg005 = -2.d0*aaa*(z-z1)
         arg006 =  aaa**2/tmp**2+2.d0*aaa*(z1-zp)
         arg101 =  tmp*(z-zp)
         arg102 =  tmp*(z1-zp)
         arg104 =  aaa/tmp+tmp*(z-zp)
         arg106 =  aaa/tmp+tmp*(z1-zp)
         !
         IF (z < z1) THEN
            t1 = -(z-zp)*qe_erf(arg101)+(0.5d0/aaa+z1-zp)*qe_erf(arg102)
            t2 = 0.5d0/aaa*exp_erfc(arg006,arg106)
            t3 = 0.5d0/aaa-(z-z1)+EXP(arg002)/tmp/SQRT(pi) &
                 -EXP(arg001)/tmp/SQRT(pi)
            kk1 = (t1+t2)/2.d0
            kk2 = t3/2.d0
         ELSE
            t1 = -exp_erfc(arg005,arg101)/aaa 
            t2 =  exp_erfc(arg006,arg104)/aaa
            t3 =  EXP(arg005)/aaa
            kk1 = (t1+t2)/4.d0
            kk2 = t3/2.d0
         ENDIF
         !
         cc1 = 0.d0
         cc2 = 0.d0
         !
         DO ng_2d = 1, ngm_2d
            !
            k1 = mill_2d(1,ng_2d)
            k2 = mill_2d(2,ng_2d)
            !
            IF( k1==0 .AND. k2==0 ) CYCLE
            t(1:2) = k1*bg(1:2, 1)+k2*bg(1:2,2)
            gp2 = SUM(t(:)*t(:))*tpiba2
            gp = SQRT(gp2)
            ff = ((k1*bg(1,1)+k2*bg(1,2))*(tau(1,it1)-tau(1,it2))  &
                 +(k1*bg(2,1)+k2*bg(2,2))*(tau(2,it1)-tau(2,it2)))*tpi
            !
            ! bc4
            alpha = aaa+gp+SQRT(aaa**2+gp**2)
            beta  = aaa+gp-SQRT(aaa**2+gp**2)
            kappa = aaa-gp+SQRT(aaa**2+gp**2)
            xi    = aaa   +SQRT(aaa**2+gp**2)
            chi   = aaa   -SQRT(aaa**2+gp**2)
            lambda=        SQRT(aaa**2+gp**2)
            !
            arg001 =  gp*(z-zp)
            arg002 = -gp*(z-zp)
            arg003 =  gp*(z+zp-2.d0*z1)
            arg005 = -gp*(z1-zp)-xi*(z-z1)
            arg006 = aaa/2.d0/tmp**2*xi+gp*(z-z1)+xi*(z1-zp)
            arg008 = aaa/2.d0/tmp**2*xi+xi*(z1-zp)-chi*(z-z1)
            arg009 = aaa/2.d0/tmp**2*xi+xi*(z1-zp)-xi*(z-z1)
            arg011 = aaa/2.d0/tmp**2*chi+chi*(z1-zp)-xi*(z-z1)
            arg101 =  gp/2.d0/tmp+tmp*(z-zp)
            arg102 =  gp/2.d0/tmp-tmp*(z-zp)
            arg103 =  gp/2.d0/tmp+tmp*(z1-zp)
            arg104 =  gp/2.d0/tmp-tmp*(z1-zp)
            arg107 =  xi/2.d0/tmp+tmp*(z-zp)
            arg109 =  xi/2.d0/tmp+tmp*(z1-zp)
            arg111 = chi/2.d0/tmp+tmp*(z-zp)
            arg113 = chi/2.d0/tmp+tmp*(z1-zp)
            !
            IF (z < z1) THEN
               t1 = exp_erfc(arg001,arg101)-exp_erfc(arg001,arg103)
               t2 = exp_erfc(arg002,arg102) &
                    -kappa/alpha*exp_erfc(arg003,arg104)
               t3 = exp_erfc(arg006,arg109)/alpha
               !
               cc1 = cc1 + COS(ff)*(t1+t2)/4.d0/gp
               cc2 = cc2 + COS(ff)*t3/2.d0
            ELSE
               t1 = exp_erfc(arg011,arg113)-exp_erfc(arg011,arg111)
               t2 = exp_erfc(arg008,arg107) &
                   -beta/alpha*exp_erfc(arg009,arg109)
               t3 = exp_erfc(arg005,arg104)/alpha
               cc1 = cc1 + COS(ff)*(t1+t2)/4.d0/lambda
               cc2 = cc2 + COS(ff)*t3/2.d0
            ENDIF
            !
         ENDDO
         !
         IF ( gamma_only ) THEN
            cc1 = cc1*2d0
            cc2 = cc2*2d0
         ENDIF
         !
         ew = ew + tt*(cc1+cc2)
         IF (gstart==2) ew = ew + tt*(kk1+kk2)
         !
      ENDDO
    ENDDO
    !
    ewg = ewg + ew
    !
    RETURN
    !
  END SUBROUTINE esm_ewaldg_bc4
  !
  !
  !----------------------------------------------------------------------------
  SUBROUTINE esm_stres_ewg( alpha, sigmaewa )
    !--------------------------------------------------------------------------
    !! Wrapper to \(\textrm{esm_stres_ewg_bc..}\) routines.
    !
    USE kinds,     ONLY : DP
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN) :: alpha
    !! alpha term in ewald sum
    REAL(DP), INTENT(OUT) :: sigmaewa(3,3)
    !! Ewald contribution to the stress
    !
    SELECT CASE( esm_bc )
    CASE( 'pbc' )
       STOP 'esm_stres_ewa must not be called for esm_bc = pbc'
    CASE( 'bc1' )
       CALL esm_stres_ewg_bc1( alpha, sigmaewa )
    CASE( 'bc2' )
       CALL esm_stres_ewg_bc2( alpha, sigmaewa )
    CASE( 'bc3' )
       CALL esm_stres_ewg_bc3( alpha, sigmaewa )
    CASE( 'bc4' )
       STOP 'esm_stres_ewa must not be called for esm_bc = bc4'
    END SELECT
    !
    RETURN
    !
  END SUBROUTINE esm_stres_ewg
  !
  !
  !-----------------------------------------------------------------
  FUNCTION qe_gauss(x) RESULT(gauss)
    !---------------------------------------------------------------
    !
    USE kinds,      ONLY : DP
    USE constants,  ONLY : sqrtpm1  ! 1/SQRT(pi)
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN) :: x
    REAL(DP) :: gauss
    !
    gauss = 2.0d0*sqrtpm1*EXP(-x*x)
    !
  END FUNCTION qe_gauss
  !
  !
  !-------------------------------------------------------------------
  SUBROUTINE esm_stres_ewg_bc1( alpha, sigmaewa )
    !-----------------------------------------------------------------
    !! Ewald G-space contribution to the stress - bc1.
    !
    USE kinds,           ONLY : DP
    USE constants,       ONLY : pi, sqrtpm1, tpi, fpi, e2
    USE cell_base,       ONLY : omega, alat, at, tpiba, bg
    USE ions_base,       ONLY : zv, nat, tau, ityp
    USE control_flags,   ONLY : gamma_only
    USE gvect,           ONLY : gstart
    USE mp_bands,        ONLY : intra_bgrp_comm
    USE mp,              ONLY : mp_sum
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN)  :: alpha
    !! alpha term in ewald sum
    REAL(DP), INTENT(OUT) :: sigmaewa(3,3)
    !! contribution to the stress tensor
    !
    ! ... local variables
    !
    INTEGER :: ia, ib, igp, iga, igb, la, mu, iz
    REAL(DP) :: L, S, salp
    REAL(DP) :: Qa, Qb, ra(2), rb(2), za, zb
    REAL(DP) :: g(2), gp, Vr
    REAL(DP) :: cosgpr, experfcm, experfcp, dexperfcm_dgp, dexperfcp_dgp
    !
    REAL(DP) :: dE_deps(2,2)
    REAL(DP), PARAMETER :: delta(2,2) = RESHAPE( (/ 1.0d0, 0.0d0, 0.0d0, 1.0d0 /), (/2,2/) )
    REAL(DP) :: dgp_deps(2,2)  !! dgp/deps
    REAL(DP) :: dinvgp_deps(2,2)  !! dgp^-1/deps
    !
    ! cell settings
    L  = at(3,3)*alat
    S  = omega/L
    salp = SQRT(alpha)
    !
    ! initialize
    sigmaewa(:,:) = 0.0d0
    !  
    !****For gp!=0 case ****
    DO ib = 1, nat
       Qb = (-1.0d0)*zv(ityp(ib))
       rb(1:2) = tau(1:2,ib)*alat
       zb = tau(3,ib)*alat
       !
       IF ( zb > L*0.5d0 ) THEN
          zb = zb - L
       ENDIF
       !
       DO ia = 1, nat
          Qa = (-1.0d0)*zv(ityp(ia))
          ra(1:2) = tau(1:2,ia)*alat
          za = tau(3,ia)*alat
          IF ( za > L*0.5d0 ) THEN
             za = za - L
          ENDIF
          !
          ! summations over gp
          dE_deps(:,:) = 0.0d0
          DO igp=1, ngm_2d
             iga = mill_2d(1,igp)
             igb = mill_2d(2,igp)
             g(1:2) = (iga*bg(1:2,1) + igb*bg(1:2,2))*tpiba
             gp = SQRT(g(1)*g(1) + g(2)*g(2))
             !
             IF (gp == 0.0d0) CYCLE ! skip gp=0
             !
             ! derivatives by strain tensor
             DO la = 1, 2
                DO mu = 1, 2
                   dgp_deps(la,mu) = -g(la)*g(mu)/gp
                   dinvgp_deps(la,mu) = +g(la)*g(mu)/gp**3
                ENDDO
             ENDDO
             !
             ! coefficients
             cosgpr = COS(g(1)*(rb(1)-ra(1)) + g(2)*(rb(2)-ra(2)))
             experfcm = exp_erfc( -gp*(zb-za), gp/2.d0/salp-salp*(zb-za) )
             experfcp = exp_erfc( +gp*(zb-za), gp/2.d0/salp+salp*(zb-za) )
             !
             dexperfcm_dgp = -(zb-za)*exp_erfc( -gp*(zb-za), gp/2.d0/salp-salp*(zb-za) ) &
                  - EXP( -gp*(zb-za) ) * qe_gauss( gp/2.d0/salp-salp*(zb-za) )/2.d0/salp
             dexperfcp_dgp = +(zb-za)*exp_erfc( +gp*(zb-za), gp/2.d0/salp+salp*(zb-za) ) &
                  - EXP( +gp*(zb-za) ) * qe_gauss( gp/2.d0/salp+salp*(zb-za) )/2.d0/salp
             !
             dE_deps(:,:) = dE_deps(:,:) &
                  + gp*dinvgp_deps(:,:) * pi/gp * Qb*Qa/S * cosgpr * experfcm &
                  - pi/gp * delta(:,:) * Qb*Qa/S * cosgpr * experfcm &
                  + pi/gp * Qb*Qa/S * cosgpr * dgp_deps(:,:) * dexperfcm_dgp &
                  + gp*dinvgp_deps(:,:) * pi/gp * Qb*Qa/S * cosgpr * experfcp &
                  - pi/gp * delta(:,:) * Qb*Qa/S * cosgpr * experfcp &
                  + pi/gp * Qb*Qa/S * cosgpr * dgp_deps(:,:) * dexperfcp_dgp
          ENDDO ! igp
          !
          ! modifications
          IF ( gamma_only ) THEN
             dE_deps(:,:) = dE_deps(:,:)*2.0d0
          ENDIF
          !
          ! calculate stress tensor
          sigmaewa(1:2,1:2) = sigmaewa(1:2,1:2) - dE_deps(1:2,1:2)/omega
          !
       ENDDO ! ia
    ENDDO ! ib
    !
    !****For gp=0 case ****
    IF ( gstart==2 ) THEN
       DO ib = 1, nat
          Qb = (-1.0d0)*zv(ityp(ib))
          rb(1:2) = tau(1:2,ib)*alat
          zb = tau(3,ib)*alat
          IF ( zb > L*0.5d0 ) THEN
             zb = zb - L
          ENDIF
          !
          Vr = 0.0d0
          DO ia = 1, nat
             Qa = (-1.0d0)*zv(ityp(ia))
             ra(1:2) = tau(1:2,ia)*alat
             za = tau(3,ia)*alat
             IF ( za > L*0.5d0 ) THEN
                za = za - L
             ENDIF
             !
             Vr = Vr - tpi * Qa/S &
                  * ( (zb-za)*qe_erf(salp*(zb-za)) &
                  + EXP(-alpha*(zb-za)**2)*sqrtpm1/salp )
          ENDDO ! ia
          !
          dE_deps(1:2,1:2) = - delta(1:2,1:2) * Vr*Qb
          !
          ! calculate stress tensor
          sigmaewa(1:2,1:2) = sigmaewa(1:2,1:2) - dE_deps(1:2,1:2)/omega
       ENDDO ! ib
    ENDIF
    !
    ! half means removing duplications.
    ! e2 means hartree -> Ry.
    sigmaewa(:,:) = sigmaewa(:,:) * (0.5d0*e2)
    !
    CALL mp_sum( sigmaewa, intra_bgrp_comm )
    !
    RETURN
    !
  END SUBROUTINE esm_stres_ewg_bc1
  !
  !
  !---------------------------------------------------------------------------
  SUBROUTINE esm_stres_ewg_bc2( alpha, sigmaewa )
    !-------------------------------------------------------------------------
    !! Ewald G-space contribution to the stress - bc2.
    !
    USE kinds,            ONLY : DP
    USE constants,        ONLY : pi, sqrtpm1, tpi, fpi, e2
    USE cell_base,        ONLY : omega, alat, at, tpiba, bg
    USE ions_base,        ONLY : zv, nat, tau, ityp
    USE control_flags,    ONLY : gamma_only
    USE gvect,            ONLY : gstart
    USE mp_bands,         ONLY : intra_bgrp_comm
    USE mp,               ONLY : mp_sum
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN)  :: alpha
    !! alpha term in ewald sum
    REAL(DP), INTENT(OUT) :: sigmaewa(3,3)
    !! contribution to the stress tensor
    !
    ! ... local variables
    !
    INTEGER :: ia, ib, igp, iga, igb, la, mu
    REAL(DP) :: L, S, salp, z0, z1
    REAL(DP) :: Qa, Qb, ra(2), rb(2), za, zb
    REAL(DP) :: g(2), gp, Vr
    REAL(DP) :: cosgpr, experfcm, experfcp, dexperfcm_dgp, dexperfcp_dgp
    REAL(DP) :: exph1, exph2, exph3
    !
    REAL(DP) :: dE_deps(2,2)
    REAL(DP), PARAMETER :: delta(2,2) = RESHAPE( (/ 1.0d0, 0.0d0, 0.0d0, 1.0d0 /), (/2,2/) )
    REAL(DP) :: dgp_deps(2,2)  !! dgp/deps
    REAL(DP) :: dinvgp_deps(2,2)  !! dgp^-1/deps
    !
    ! cell settings
    L  = at(3,3)*alat
    S  = omega/L
    z0 = L/2.d0
    z1 = z0 + esm_w
    salp = SQRT(alpha)
    !
    ! initialize
    sigmaewa(:,:) = 0.0d0
    !
    !****For gp!=0 case ****
    DO ib = 1, nat
       Qb = (-1.0d0)*zv(ityp(ib))
       rb(1:2) = tau(1:2,ib)*alat
       zb = tau(3,ib)*alat
       IF ( zb > L*0.5d0 ) THEN
          zb = zb - L
       ENDIF
       !
       DO ia = 1, nat
          Qa = (-1.0d0)*zv(ityp(ia))
          ra(1:2) = tau(1:2,ia)*alat
          za = tau(3,ia)*alat
          IF ( za > L*0.5d0 ) THEN
             za = za - L
          ENDIF
          !
          ! summations over gp
          dE_deps(:,:) = 0.0d0
          DO igp = 1, ngm_2d
             iga = mill_2d(1,igp)
             igb = mill_2d(2,igp)
             g(1:2) = (iga*bg(1:2,1) + igb*bg(1:2,2))*tpiba
             gp = SQRT(g(1)*g(1) + g(2)*g(2))
             !
             IF (gp == 0.0d0) CYCLE ! skip gp=0
             !
             ! derivatives by strain tensor
             DO la = 1, 2
                DO mu = 1, 2
                   dgp_deps(la,mu) = -g(la)*g(mu)/gp
                   dinvgp_deps(la,mu) = +g(la)*g(mu)/gp**3
                ENDDO
             ENDDO
             !
             ! coefficients
             cosgpr = COS(g(1)*(rb(1)-ra(1)) + g(2)*(rb(2)-ra(2)))
             !
             experfcm = exp_erfc( -gp*(zb-za), gp/2.d0/salp-salp*(zb-za) )
             experfcp = exp_erfc( +gp*(zb-za), gp/2.d0/salp+salp*(zb-za) )
             !
             dexperfcm_dgp = -(zb-za)*exp_erfc( -gp*(zb-za), gp/2.d0/salp-salp*(zb-za) ) &
                  - EXP( -gp*(zb-za) ) * qe_gauss( gp/2.d0/salp-salp*(zb-za) )/2.d0/salp
             dexperfcp_dgp = +(zb-za)*exp_erfc( +gp*(zb-za), gp/2.d0/salp+salp*(zb-za) ) &
                  - EXP( +gp*(zb-za) ) * qe_gauss( gp/2.d0/salp+salp*(zb-za) )/2.d0/salp
             !
             exph1  = (COSH(gp*(zb-za))*EXP(-2*gp*z1) - COSH(gp*(zb+za)) )/SINH(2*gp*z1)
             exph2  = ( (zb-za)*SINH(gp*(zb-za))*EXP(-2*gp*z1) &
                      - 2*z1*COSH(gp*(zb-za))*EXP(-2*gp*z1)    &
                      - (zb+za)*SINH(gp*(zb+za)) )/SINH(2*gp*z1)
             exph3  = - (COSH(gp*(zb-za))*EXP(-2*gp*z1) - COSH(gp*(zb+za)) )/ &
                        SINH(2*gp*z1)**2 * 2*z1*COSH(2*gp*z1)
             !
             ! BC1 terms
             dE_deps(:,:) = dE_deps(:,:) &
                  + gp*dinvgp_deps(:,:) * pi/gp * Qb*Qa/S * cosgpr * experfcm &
                  - pi/gp * delta(:,:) * Qb*Qa/S * cosgpr * experfcm          &
                  + pi/gp * Qb*Qa/S * cosgpr * dgp_deps(:,:) * dexperfcm_dgp  &
                  + gp*dinvgp_deps(:,:) * pi/gp * Qb*Qa/S * cosgpr * experfcp &
                  - pi/gp * delta(:,:) * Qb*Qa/S * cosgpr * experfcp          &
                  + pi/gp * Qb*Qa/S * cosgpr * dgp_deps(:,:) * dexperfcp_dgp
             !
             ! BC2 terms
             dE_deps(:,:) = dE_deps(:,:) &
                  + gp*dinvgp_deps(:,:) * tpi/gp * Qb*Qa/S * cosgpr * exph1 &
                  - tpi/gp * delta(:,:) * Qb*Qa/S * cosgpr * exph1 &
                  + tpi/gp * Qb*Qa/S * cosgpr * dgp_deps(:,:) * (exph2+exph3)
          ENDDO ! igp
          !
          ! modifications
          IF ( gamma_only ) THEN
             dE_deps(:,:) = dE_deps(:,:)*2.0d0
          ENDIF
          !
          ! calculate stress tensor
          sigmaewa(1:2,1:2) = sigmaewa(1:2,1:2) - dE_deps(1:2,1:2)/omega
          !
       ENDDO ! ia
    ENDDO ! ib
    !
    !****For gp=0 case ****
    IF (gstart == 2) THEN
       DO ib = 1, nat
          !
          Qb = (-1.0d0)*zv(ityp(ib))
          rb(1:2) = tau(1:2,ib)*alat
          zb = tau(3,ib)*alat
          IF ( zb > L*0.5d0 ) THEN
             zb = zb - L
          ENDIF
          !
          ! [note] this Vr does not contain a term due to efield z*efield
          ! because it vanishes in the differentiation with respect to strain.
          Vr = 0.0d0
          DO ia = 1, nat
             Qa = (-1.0d0)*zv(ityp(ia))
             ra(1:2) = tau(1:2,ia)*alat
             za = tau(3,ia)*alat
             IF ( za > L*0.5d0 ) THEN
                za = za - L
             ENDIF
             !
             ! BC1 terms
             Vr = Vr - tpi * Qa/S &
                  * ( (zb-za)*qe_erf(salp*(zb-za)) &
                  + EXP(-alpha*(zb-za)**2)*sqrtpm1/salp )

             !! BC2 terms
             Vr = Vr + tpi * Qa/S * ( -zb*za + z1*z1 )/z1
          ENDDO ! ia
          !
          dE_deps(1:2,1:2) = - delta(1:2,1:2) * Vr*Qb

          ! calculate stress tensor
          sigmaewa(1:2,1:2) = sigmaewa(1:2,1:2) - dE_deps(1:2,1:2)/omega
          !
       ENDDO ! ib
    ENDIF
    !
    ! half means removing duplications.
    ! e2 means hartree -> Ry.
    sigmaewa(:,:) = sigmaewa(:,:) * (0.5d0*e2)
    !
    CALL mp_sum( sigmaewa, intra_bgrp_comm )
    !
    RETURN
    !
  END SUBROUTINE esm_stres_ewg_bc2
  !
  !
  !-----------------------------------------------------------------------------
  SUBROUTINE esm_stres_ewg_bc3( alpha, sigmaewa )
    !--------------------------------------------------------------------------
    !! Ewald G-space contribution to the stress - bc3.
    !
    USE kinds,         ONLY : DP
    USE constants,     ONLY : pi, sqrtpm1, tpi, fpi, e2
    USE cell_base,     ONLY : omega, alat, at, tpiba, bg
    USE ions_base,     ONLY : zv, nat, tau, ityp
    USE control_flags, ONLY : gamma_only
    USE gvect,         ONLY : gstart
    USE mp_bands,      ONLY : intra_bgrp_comm
    USE mp,            ONLY : mp_sum
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN)  :: alpha
    !! alpha term in ewald sum
    REAL(DP), INTENT(OUT) :: sigmaewa(3,3)
    !! contribution to the stress tensor
    !
    ! ... local variables
    !
    INTEGER  :: ia, ib, igp, iga, igb, la, mu
    REAL(DP) :: L, S, salp, z0, z1
    REAL(DP) :: Qa, Qb, ra(2), rb(2), za, zb
    REAL(DP) :: g(2), gp, Vr
    REAL(DP) :: cosgpr, experfcm, experfcp, dexperfcm_dgp, dexperfcp_dgp, expm
    !
    REAL(DP) :: dE_deps(2,2)
    REAL(DP), PARAMETER :: delta(2,2) = RESHAPE( (/ 1.0d0, 0.0d0, 0.0d0, 1.0d0 /), (/2,2/) )
    REAL(DP) :: dgp_deps(2,2)  !! dgp/deps
    REAL(DP) :: dinvgp_deps(2,2)  !! dgp^-1/deps
    !
    ! cell settings
    L  = at(3,3)*alat
    S  = omega/L
    z0 = L/2.d0
    z1 = z0+esm_w
    salp = SQRT(alpha)
    !
    ! initialize
    sigmaewa(:,:) = 0.0d0
    ! 
    !****For gp!=0 case ****
    DO ib = 1, nat
       Qb = (-1.0d0)*zv(ityp(ib))
       rb(1:2) = tau(1:2,ib)*alat
       zb = tau(3,ib)*alat
       IF ( zb > L*0.5d0 ) THEN
          zb = zb - L
       ENDIF
       !
       DO ia = 1, nat
          Qa = (-1.0d0)*zv(ityp(ia))
          ra(1:2) = tau(1:2,ia)*alat
          za = tau(3,ia)*alat
          IF ( za > L*0.5d0 ) THEN
             za = za - L
          ENDIF
          !
          ! summations over gp
          dE_deps(:,:) = 0.0d0
          !
          DO igp = 1, ngm_2d
             iga = mill_2d(1,igp)
             igb = mill_2d(2,igp)
             g(1:2) = (iga*bg(1:2,1) + igb*bg(1:2,2))*tpiba
             gp = SQRT(g(1)*g(1) + g(2)*g(2))
             !
             IF ( gp==0.0d0 ) CYCLE ! skip gp=0
             !
             ! derivatives by strain tensor
             DO la = 1, 2
                DO mu = 1, 2
                   dgp_deps(la,mu) = -g(la)*g(mu)/gp
                   dinvgp_deps(la,mu) = +g(la)*g(mu)/gp**3
                ENDDO
             ENDDO
             !
             ! coefficients
             cosgpr = COS(g(1)*(rb(1)-ra(1)) + g(2)*(rb(2)-ra(2)))
             experfcm = exp_erfc( -gp*(zb-za), gp/2.d0/salp-salp*(zb-za) )
             experfcp = exp_erfc( +gp*(zb-za), gp/2.d0/salp+salp*(zb-za) )
             dexperfcm_dgp = -(zb-za)*exp_erfc( -gp*(zb-za), gp/2.d0/salp-salp*(zb-za) ) &
                  - EXP( -gp*(zb-za) ) * qe_gauss( gp/2.d0/salp-salp*(zb-za) )/2.d0/salp
             dexperfcp_dgp = +(zb-za)*exp_erfc( +gp*(zb-za), gp/2.d0/salp+salp*(zb-za) ) &
                  - EXP( +gp*(zb-za) ) * qe_gauss( gp/2.d0/salp+salp*(zb-za) )/2.d0/salp
             expm = EXP( -gp*(-zb+2*z1-za) )
             !
             ! BC1 terms
             dE_deps(:,:) = dE_deps(:,:) &
                  + gp*dinvgp_deps(:,:) * pi/gp * Qb*Qa/S * cosgpr * experfcm &
                  - pi/gp * delta(:,:) * Qb*Qa/S * cosgpr * experfcm          &
                  + pi/gp * Qb*Qa/S * cosgpr * dgp_deps(:,:) * dexperfcm_dgp  &
                  + gp*dinvgp_deps(:,:) * pi/gp * Qb*Qa/S * cosgpr * experfcp &
                  - pi/gp * delta(:,:) * Qb*Qa/S * cosgpr * experfcp          &
                  + pi/gp * Qb*Qa/S * cosgpr * dgp_deps(:,:) * dexperfcp_dgp
             !
             ! BC3 terms
             dE_deps(:,:) = dE_deps(:,:) &
                  - gp*dinvgp_deps(:,:) * tpi/gp * Qb*Qa/S * cosgpr * expm &
                  + tpi/gp * delta(:,:) * Qb*Qa/S * cosgpr * expm          &
                  + tpi/gp * Qb*Qa/S * cosgpr * dgp_deps(:,:) * (-zb+2*z1-za) * expm
          ENDDO ! igp
          !
          ! modifications
          IF ( gamma_only ) THEN
             dE_deps(:,:) = dE_deps(:,:)*2.0d0
          ENDIF
          !
          ! calculate stress tensor
          sigmaewa(1:2,1:2) = sigmaewa(1:2,1:2) - dE_deps(1:2,1:2)/omega
          !
       ENDDO ! ia
    ENDDO ! ib
    !
    !****For gp=0 case ****
    IF ( gstart==2 ) THEN
       DO ib = 1, nat
          Qb = (-1.0d0)*zv(ityp(ib))
          rb(1:2) = tau(1:2,ib)*alat
          zb = tau(3,ib)*alat
          IF ( zb > L*0.5d0 ) THEN
             zb = zb - L
          ENDIF
          !
          Vr = 0.0d0
          DO ia = 1, nat
             Qa = (-1.0d0)*zv(ityp(ia))
             ra(1:2) = tau(1:2,ia)*alat
             za = tau(3,ia)*alat
             IF ( za > L*0.5d0 ) THEN
                za = za - L
             ENDIF
             !
             ! BC1 terms
             Vr = Vr - tpi * Qa/S &
                  * ( (zb-za)*qe_erf(salp*(zb-za)) &
                  + EXP(-alpha*(zb-za)**2)*sqrtpm1/salp )
             !
             ! BC3 terms
             Vr = Vr + tpi * Qa/S * ( -zb+2*z1-za )
          ENDDO ! ia
          !
          dE_deps(1:2,1:2) = - delta(1:2,1:2) * Vr*Qb
          !
          ! calculate stress tensor
          sigmaewa(1:2,1:2) = sigmaewa(1:2,1:2) - dE_deps(1:2,1:2)/omega
       ENDDO ! ib
    ENDIF
    !
    ! half means removing duplications.
    ! e2 means hartree -> Ry.
    sigmaewa(:,:) = sigmaewa(:,:) * (0.5d0*e2)
    !
    CALL mp_sum( sigmaewa, intra_bgrp_comm )
    !
    RETURN
    !
  END SUBROUTINE esm_stres_ewg_bc3
  !
  !
  ! ... ESM LOCAL POTENTIAL SUBROUTINES
  !
  !-----------------------------------------------------------------------
  SUBROUTINE esm_local_pbc( aux )
    !----------------------------------------------------------------------
    !
    USE fft_base, ONLY : dfftp
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(INOUT) :: aux(dfftp%nnr)
    !
    STOP 'esm_local must not be called for esm_bc = pbc'
    !
  END SUBROUTINE esm_local_pbc
  !
  !
  !--------------------------------------------------------------------------
  SUBROUTINE esm_local_bc1( aux )
    !-----------------------------------------------------------------------
    !! Calculates ESM local potential - bc1.
    !
    USE constants,        ONLY : pi, tpi, fpi, e2
    USE gvect,            ONLY : ngm, mill
    USE control_flags,    ONLY : gamma_only
    USE cell_base,        ONLY : at, bg, alat, tpiba2, omega
    USE ions_base,        ONLY : zv, nat, tau, ityp
    USE fft_base,         ONLY : dfftp
    USE fft_scalar,       ONLY : cft_1z
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(INOUT) :: aux(dfftp%nnr)
    !! aux contains v_loc_short(G) (input) and v_loc(G) (output)
    !
    ! ... local variables
    !
    REAL(DP) :: t(2), tt, gp, gp2, sa, z0, pp, cc, ss, &
                t1, t2, z, zp, tmp, L, z_l, z_r,       &
                arg001, arg002, arg101, arg102
    INTEGER :: iz, it, k1, k2, k3, ng, n1, n2, n3, nz_l, &
               nz_r, ng_2d
    COMPLEX(DP) :: cs, cc1, cc2, a0, a1, a2, a3, f1, f2, &
                   f3, f4
    COMPLEX(DP), ALLOCATABLE :: vloc3(:,:), vg(:), vg_r(:)
    !
    L = at(3,3)*alat
    sa = omega/L
    z0 = L/2.d0
    tmp = 1.d0   ! Gaussian width
    !
    ALLOCATE( vloc3(dfftp%nr3,ngm_2d) )
    !
    !FOR gp!=0
    ALLOCATE( vg(dfftp%nr3), vg_r(dfftp%nr3) )
    !
    DO ng_2d = 1, ngm_2d
       k1 = mill_2d(1,ng_2d)
       k2 = mill_2d(2,ng_2d)
       IF (k1==0 .AND. k2==0) CYCLE
       t(1:2) = k1*bg(1:2,1)+k2*bg(1:2,2)
       gp2 = SUM(t(:)*t(:))*tpiba2
       gp  = SQRT(gp2)
       vg_r(:) = (0.d0,0.d0)
       !
       DO it = 1, nat
          tt = -fpi*zv(ityp(it))/sa
          pp = -tpi*(tau(1,it)*(k1*bg(1,1)+k2*bg(1,2)) &
                    +tau(2,it)*(k1*bg(2,1)+k2*bg(2,2)))
          cc = COS(pp)
          ss = SIN(pp)
          cs = CMPLX(cc,ss,KIND=DP)
          !
          zp = tau(3,it)
          IF (zp > at(3,3)*0.5) zp = zp-at(3,3)
          zp = zp*alat
          !
          DO iz = 1, dfftp%nr3
             k3 = iz-1
             IF (k3>dfftp%nr3/2) k3 = iz-dfftp%nr3-1
             z = DBLE(k3)/DBLE(dfftp%nr3)*L
             ! bc1
             arg001 =  gp*(z-zp)
             arg002 = -gp*(z-zp)
             arg101 =  gp/2.d0/tmp+tmp*(z-zp)
             arg102 =  gp/2.d0/tmp-tmp*(z-zp)
             !
             t1 = exp_erfc(arg002,arg102)
             t2 = exp_erfc(arg001,arg101)
             !
             cc1 = cs*(t1+t2)/4.d0/gp
             cc2 = (0.d0,0.d0)
             !
             vg_r(iz) = vg_r(iz)+tt*(cc1+cc2)*e2 ! factor e2: hartree -> Ry.
          ENDDO
          !
       ENDDO
       !
       CALL cft_1z( vg_r, 1, dfftp%nr3, dfftp%nr3, -1, vg )
       DO iz = 1, dfftp%nr3
          vloc3(iz,ng_2d) = vg(iz)
       ENDDO
    ENDDO
    !
    DEALLOCATE( vg, vg_r )
    !
    ng_2d = imill_2d(0,0)
    !
    IF ( ng_2d > 0 ) THEN
       ALLOCATE( vg(dfftp%nr3), vg_r(dfftp%nr3) )
       vg_r(:) = (0.d0,0.d0)
       !for smoothing
       f1=(0.d0,0.d0); f2=(0.d0,0.d0); f3=(0.d0,0.d0); f4=(0.d0,0.d0)
       nz_l = dfftp%nr3/2+1+esm_nfit
       nz_r = dfftp%nr3/2+1-esm_nfit
       z_l = DBLE(nz_l-1)*L/DBLE(dfftp%nr3)-L
       z_r = DBLE(nz_r-1)*L/DBLE(dfftp%nr3)
       !for gp = 0
       DO it = 1, nat
          tt = -fpi*zv(ityp(it))/sa
          zp = tau(3,it)
          IF (zp > at(3,3)*0.5) zp = zp-at(3,3)
          zp = zp*alat
          DO iz = 1, dfftp%nr3
             k3 = iz-1
             IF (k3>dfftp%nr3/2) k3 = iz-dfftp%nr3-1
             z = DBLE(k3)/DBLE(dfftp%nr3)*L
             ! bc1
             arg001 = -tmp**2*(z-zp)**2
             arg101 =  tmp*(z-zp)
             cc1 = 0.5d0*(-(z-zp)*qe_erf(arg101)-EXP(arg001)/tmp/SQRT(pi))
             cc2 = (0.d0,0.d0)
             !
             vg_r(iz) = vg_r(iz)+tt*(cc1+cc2)*e2 ! factor e2: hartree -> Ry.
          ENDDO
          ! smoothing cell edge potential (avoiding unphysical oscillation)
          ! bc1
          f1 = f1+tt*0.5d0*(-(z_r-zp)*qe_erf(tmp*(z_r-zp)) &
               -EXP(-tmp**2*(z_r-zp)**2)/tmp/SQRT(pi))
          f2 = f2+tt*0.5d0*(-(z_l-zp)*qe_erf(tmp*(z_l-zp)) &
               -EXP(-tmp**2*(z_l-zp)**2)/tmp/SQRT(pi))
          f3 = f3-tt*0.5d0*qe_erf(tmp*(z_r-zp))
          f4 = f4-tt*0.5d0*qe_erf(tmp*(z_l-zp))
       ENDDO
       ! for smoothing
       ! factor e2: hartree -> Ry.
       f1=f1*e2; f2=f2*e2; f3=f3*e2; f4=f4*e2
       !
       z_r = z_r
       z_l = z_l+L
       !
       a0 = (f1*z_l**2*(z_l-3.d0*z_r)+z_r*(f3*z_l**2*(-z_l+z_r)    &
             +z_r*(f2*(3.d0*z_l-z_r)+f4*z_l*(-z_l+z_r))))/(z_l-z_r)**3
       a1 = (f3*z_l**3+z_l*(6.d0*f1-6.d0*f2+(f3+2.d0*f4)*z_l)*z_r  &
             -(2*f3+f4)*z_l*z_r**2-f4*z_r**3)/(z_l-z_r)**3
       a2 = (-3*f1*(z_l+z_r)+3.d0*f2*(z_l+z_r)-(z_l-z_r)*(2*f3*z_l &
             +f4*z_l+f3*z_r+2*f4*z_r))/(z_l-z_r)**3
       a3 = (2.d0*f1-2.d0*f2+(f3+f4)*(z_l-z_r))/(z_l-z_r)**3
       !
       DO iz = nz_r, nz_l
          z = DBLE(iz-1)/DBLE(dfftp%nr3)*L
          vg_r(iz) = (a0+a1*z+a2*z**2+a3*z**3)
       ENDDO
       !
       CALL cft_1z( vg_r, 1, dfftp%nr3, dfftp%nr3, -1, vg )
       !
       DO iz = 1, dfftp%nr3
          vloc3(iz,ng_2d) = vg(iz)
       ENDDO
       !
       DEALLOCATE( vg, vg_r )
    ENDIF ! if( ng_2d > 0 )
    !
    ! Map to FFT mesh (dfftp%nrx)
    DO ng = 1, ngm
       n1 = mill(1,ng)
       n2 = mill(2,ng)
       ng_2d = imill_2d(n1,n2)
       n3 = mill(3,ng) + 1 
       IF (n3<1) n3 = n3 + dfftp%nr3
       aux(dfftp%nl(ng)) = aux(dfftp%nl(ng)) + vloc3(n3,ng_2d)
       IF (gamma_only) THEN
          aux(dfftp%nlm(ng)) = CONJG(aux(dfftp%nl(ng)))
       ENDIF
    ENDDO
    !
    DEALLOCATE( vloc3 )
    !
    !
    RETURN
    !
  END SUBROUTINE esm_local_bc1
  !
  !
  !-------------------------------------------------------------------------
  SUBROUTINE esm_local_bc2( aux )
    !-----------------------------------------------------------------------
    !! Calculates ESM local potential - bc2.
    !
    USE constants,        ONLY : pi, tpi, fpi, e2
    USE gvect,            ONLY : ngm, mill
    USE control_flags,    ONLY : gamma_only
    USE cell_base,        ONLY : at, bg, alat, tpiba2, omega
    USE ions_base,        ONLY : zv, nat, tau, ityp
    USE fft_base,         ONLY : dfftp
    USE fft_scalar,       ONLY : cft_1z
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(INOUT) :: aux( dfftp%nnr )
    !! aux contains v_loc_short(G) (input) and v_loc(G) (output)
    !
    ! ... local variables
    !
    REAL(DP) :: t(2), tt, gp, gp2, sa, z1, z0, pp, cc, ss, &
                t1, t2, z, zp, v0, tmp, L, z_l, z_r,       &
                arg001, arg002, arg003, arg004, arg005,    &
                arg006, arg007, arg101, arg102
    INTEGER :: iz, it, k1, k2, k3, ng, n1, n2, n3, nz_l,   &
               nz_r, ng_2d
    COMPLEX(DP) :: cs, cc1, cc2, a0, a1, a2, a3, f1, f2,   &
                   f3, f4
    COMPLEX(DP), ALLOCATABLE :: vloc3(:,:), vg(:), vg_r(:)
    !
    L = at(3,3)*alat
    sa = omega/L
    z0 = L/2.d0
    tmp = 1.d0   ! Gaussian width
    z1 = z0+esm_w
    v0 = esm_efield*z1*2.d0/e2 ! factor 1/e2: unit Ry. -> hartree
    ALLOCATE( vloc3(dfftp%nr3,ngm_2d) )
    !
    ! for gp!=0
    ALLOCATE( vg(dfftp%nr3), vg_r(dfftp%nr3) )
    !
    DO ng_2d = 1, ngm_2d
       k1 = mill_2d(1,ng_2d)
       k2 = mill_2d(2,ng_2d)
       IF (k1==0 .AND. k2==0) CYCLE
       t(1:2) = k1*bg(1:2, 1)+k2*bg (1:2, 2)
       gp2 = SUM(t(:)*t(:))*tpiba2
       gp  = SQRT(gp2)
       vg_r(:) = (0.d0,0.d0)
       !
       DO it = 1, nat
          tt = -fpi*zv(ityp(it))/sa
          pp = -tpi*(tau(1,it)*(k1*bg(1,1)+k2*bg(1,2)) &
                    +tau(2,it)*(k1*bg(2,1)+k2*bg(2,2)))
          cc = COS(pp)
          ss = SIN(pp)
          cs = CMPLX(cc,ss,KIND=DP)
          zp = tau(3,it)
          IF (zp > at(3,3)*0.5) zp = zp-at(3,3)
          zp = zp*alat
          !
          DO iz = 1, dfftp%nr3
             k3 = iz-1
             IF (k3>dfftp%nr3/2) k3 = iz-dfftp%nr3-1
             z = DBLE(k3)/DBLE(dfftp%nr3)*L
             ! bc2
             arg001 =  gp*(z-zp)
             arg002 = -gp*(z-zp)
             arg003 = -gp*(z+zp+2.d0*z1)
             arg004 =  gp*(z+zp-2.d0*z1)
             arg005 = -gp*(z-zp+4.d0*z1)
             arg006 =  gp*(z-zp-4.d0*z1)
             arg007 = -4.d0*gp*z1
             arg101 =  gp/2.d0/tmp+tmp*(z-zp)
             arg102 =  gp/2.d0/tmp-tmp*(z-zp)
             !
             t1 = exp_erfc(arg002,arg102)
             t2 = exp_erfc(arg001,arg101)
             !
             cc1 = cs*(t1+t2)/4.d0/gp
             cc2 = cs*(EXP(arg006)+EXP(arg005)-EXP(arg004)-EXP(arg003)) &
                   /(1.d0-EXP(arg007))/2.d0/gp 
             !
             vg_r(iz) = vg_r(iz)+tt*(cc1+cc2)*e2 ! factor e2: hartree -> Ry.
          ENDDO
          !
       ENDDO
       !
       CALL cft_1z( vg_r, 1, dfftp%nr3, dfftp%nr3, -1, vg )
       !
       DO iz = 1,dfftp%nr3
          vloc3(iz,ng_2d) = vg(iz)
       ENDDO
       !
    ENDDO
    !
    DEALLOCATE( vg, vg_r )
    !
    ng_2d = imill_2d(0,0)
    !
    IF ( ng_2d > 0 ) THEN
       ALLOCATE( vg(dfftp%nr3), vg_r(dfftp%nr3) )
       vg_r(:)=(0.d0,0.d0)
       ! for smoothing
       f1=(0.d0,0.d0); f2=(0.d0,0.d0); f3=(0.d0,0.d0); f4=(0.d0,0.d0)
       nz_l = dfftp%nr3/2+1+esm_nfit
       nz_r = dfftp%nr3/2+1-esm_nfit
       !
       z_l = DBLE(nz_l-1)*L/DBLE(dfftp%nr3)-L
       z_r = DBLE(nz_r-1)*L/DBLE(dfftp%nr3)
       ! add constant potential (capacitor term)
       ! bc2
       DO iz = 1, dfftp%nr3
          k3 = iz-1
          IF (k3>dfftp%nr3/2) k3 = iz-dfftp%nr3-1
          z = DBLE(k3)/DBLE(dfftp%nr3)*L
          vg_r(iz) = -0.5d0*v0*(z-z1)/z1*e2 ! factor e2: hartree -> Ry.
       ENDDO
       !
       f1 = -0.5d0*v0*(z_r-z1)/z1 ! unit: hartree
       f2 = -0.5d0*v0*(z_l-z1)/z1 ! unit: hartree
       f3 = -0.5d0*v0/z1 ! unit: hartree/a.u.
       f4 = -0.5d0*v0/z1 ! unit: harteee/a.u.
       !
       ! for gp=0
       DO it = 1, nat
          tt = -fpi*zv(ityp(it))/sa
          zp = tau(3,it)
          IF (zp>at(3,3)*0.5) zp = zp-at(3,3)
          zp = zp*alat
          DO iz = 1, dfftp%nr3
             k3=iz-1
             IF (k3 > dfftp%nr3/2) k3 = iz-dfftp%nr3-1
             z = DBLE(k3)/DBLE(dfftp%nr3)*L
             ! bc2
             arg001 = -tmp**2*(z-zp)**2
             arg101 =  tmp*(z-zp)
             !
             cc1 = 0.5d0*(-(z-zp)*qe_erf(arg101)-EXP(arg001)/tmp/SQRT(pi))
             cc2 = 0.5d0*(z1-z*zp/z1)
             !
             vg_r(iz) = vg_r(iz)+tt*(cc1+cc2)*e2 ! factor e2: hartree -> Ry.
          ENDDO
          ! smoothing cell edge potential (avoiding unphysical oscillation)
          ! bc2
          f1 = f1 + tt*0.5d0*(-(z_r-zp)*qe_erf(tmp*(z_r-zp)) &
                    -EXP(-tmp**2*(z_r-zp)**2)/tmp/SQRT(pi))
          f2 = f2 + tt*0.5d0*(-(z_l-zp)*qe_erf(tmp*(z_l-zp)) &
                    -EXP(-tmp**2*(z_l-zp)**2)/tmp/SQRT(pi))
          f3 = f3 - tt*0.5d0*qe_erf(tmp*(z_r-zp))
          f4 = f4 - tt*0.5d0*qe_erf(tmp*(z_l-zp))
          !
          f1 = f1 + tt*0.5d0*(z1-z_r*zp/z1)
          f2 = f2 + tt*0.5d0*(z1-z_l*zp/z1)
          f3 = f3 + tt*(-0.5d0*(zp/z1))
          f4 = f4 + tt*(-0.5d0*(zp/z1))
       ENDDO
       ! for smoothing
       ! factor e2: hartree -> Ry.
       f1=f1*e2; f2=f2*e2; f3=f3*e2; f4=f4*e2
       !
       z_r = z_r
       z_l = z_l+L
       !
       a0 = (f1*z_l**2*(z_l-3.d0*z_r)+z_r*(f3*z_l**2*(-z_l+z_r) &
             +z_r*(f2*(3.d0*z_l-z_r)+f4*z_l*(-z_l+z_r))))/(z_l-z_r)**3
       a1 = (f3*z_l**3+z_l*(6.d0*f1-6.d0*f2+(f3+2.d0*f4)*z_l)*z_r &
             -(2*f3+f4)*z_l*z_r**2-f4*z_r**3)/(z_l-z_r)**3
       a2 = (-3*f1*(z_l+z_r)+3.d0*f2*(z_l+z_r)-(z_l-z_r)*(2*f3*z_l &
             +f4*z_l+f3*z_r+2*f4*z_r))/(z_l-z_r)**3
       a3 = (2.d0*f1-2.d0*f2+(f3+f4)*(z_l-z_r))/(z_l-z_r)**3
       !
       DO iz = nz_r, nz_l
          z = DBLE(iz-1)/DBLE(dfftp%nr3)*L
          vg_r(iz) = (a0+a1*z+a2*z**2+a3*z**3)
       ENDDO
       !
       CALL cft_1z( vg_r, 1, dfftp%nr3, dfftp%nr3, -1, vg )
       !
       DO iz = 1, dfftp%nr3
          vloc3(iz,ng_2d) = vg(iz)
       ENDDO
       !
       DEALLOCATE( vg, vg_r )
       !
    ENDIF ! if( ng_2d > 0 )
    !
    !Map to FFT mesh (dfftp%nrx)
    DO ng = 1, ngm
       n1 = mill(1,ng)
       n2 = mill(2,ng)
       ng_2d = imill_2d(n1,n2)
       n3 = mill(3,ng) + 1 
       IF (n3<1) n3 = n3 + dfftp%nr3
       aux(dfftp%nl(ng)) = aux(dfftp%nl(ng)) + vloc3(n3,ng_2d)
       IF (gamma_only) THEN
          aux(dfftp%nlm(ng)) = CONJG(aux(dfftp%nl(ng)))
       ENDIF
    ENDDO
    !
    DEALLOCATE(vloc3)
    !
    RETURN
    !
  END SUBROUTINE esm_local_bc2
  !
  !
  !--------------------------------------------------------------------------
  SUBROUTINE esm_local_bc3( aux )
    !------------------------------------------------------------------------
    !! Calculates ESM local potential - bc3.
    !
    USE constants,        ONLY : pi, tpi, fpi, e2
    USE gvect,            ONLY : ngm, mill
    USE control_flags,    ONLY : gamma_only
    USE cell_base,        ONLY : at, bg, alat, tpiba2, omega
    USE ions_base,        ONLY : zv, nat, tau, ityp
    USE fft_base,         ONLY : dfftp
    USE fft_scalar,       ONLY : cft_1z
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(INOUT) :: aux(dfftp%nnr)
    !! aux contains v_loc_short(G) (input) and v_loc(G) (output)
    !
    ! ... local variables
    !
    REAL(DP) :: t(2), tt, gp, gp2, sa, z1, z0, pp, cc, ss, &
                t1, t2, z, zp, tmp, L, z_l, z_r,           &
                arg001, arg002, arg003, arg101, arg102
    INTEGER :: iz, it, k1, k2, k3, ng, n1, n2, n3, nz_l,   &
               nz_r, ng_2d
    COMPLEX(DP) :: cs, cc1, cc2, a0, a1, a2, a3, f1, f2,   &
                   f3, f4
    COMPLEX(DP), ALLOCATABLE :: vloc3(:,:), vg(:), vg_r(:)
    !
    L = at(3,3)*alat
    sa = omega/L
    z0 = L/2.d0
    tmp = 1.d0   ! Gaussian width
    z1 = z0+esm_w
    !
    ALLOCATE( vloc3(dfftp%nr3,ngm_2d) )
    !
    ! for gp!=0
    ALLOCATE( vg(dfftp%nr3), vg_r(dfftp%nr3) )
    !
    DO ng_2d = 1, ngm_2d
       k1 = mill_2d(1,ng_2d)
       k2 = mill_2d(2,ng_2d)
       IF (k1==0 .AND. k2==0) CYCLE
       t(1:2) = k1*bg(1:2,1)+k2*bg(1:2,2)
       gp2 = SUM(t(:)*t(:))*tpiba2
       gp  = SQRT(gp2)
       vg_r(:) = (0.d0,0.d0)
       !
       DO it = 1, nat
          tt = -fpi*zv(ityp(it))/sa
          pp = -tpi*(tau(1,it)*(k1*bg(1,1)+k2*bg(1,2)) &
                    +tau(2,it)*(k1*bg(2,1)+k2*bg(2,2)))
          cc = COS(pp)
          ss = SIN(pp)
          cs = CMPLX(cc,ss,KIND=DP)
          zp = tau(3,it)
          IF (zp > at(3,3)*0.5) zp = zp-at(3,3)
          zp = zp*alat
          DO iz = 1, dfftp%nr3
             k3 = iz-1
             IF (k3 > dfftp%nr3/2) k3 = iz-dfftp%nr3-1
             z = DBLE(k3)/DBLE(dfftp%nr3)*L
             ! bc3
             arg001 =  gp*(z-zp)
             arg002 = -gp*(z-zp)
             arg003 =  gp*(z+zp-2.d0*z1)
             arg101 =  gp/2.d0/tmp+tmp*(z-zp)
             arg102 =  gp/2.d0/tmp-tmp*(z-zp)
             !
             t1 = exp_erfc(arg002,arg102)
             t2 = exp_erfc(arg001,arg101)
             !
             cc1 = cs*(t1+t2)/4.d0/gp
             cc2 = cs*(-EXP(arg003))/2.d0/gp
             !
             vg_r(iz) = vg_r(iz)+tt*(cc1+cc2)*e2 ! factor e2: hartree -> Ry.
          ENDDO
       ENDDO
       !
       CALL cft_1z( vg_r, 1, dfftp%nr3, dfftp%nr3, -1, vg )
       !
       DO iz = 1, dfftp%nr3
          vloc3(iz,ng_2d) = vg(iz)
       ENDDO
       !
    ENDDO
    !
    DEALLOCATE(vg,vg_r)
    !
    ng_2d = imill_2d(0,0)
    !
    IF ( ng_2d > 0 ) THEN
       ALLOCATE( vg(dfftp%nr3), vg_r(dfftp%nr3) )
       vg_r(:) = (0.d0,0.d0)
       ! for smoothing
       f1=(0.d0,0.d0); f2=(0.d0,0.d0); f3=(0.d0,0.d0); f4=(0.d0,0.d0)
       nz_l = dfftp%nr3/2+1+esm_nfit
       nz_r = dfftp%nr3/2+1-esm_nfit
       z_l = DBLE(nz_l-1)*L/DBLE(dfftp%nr3)-L
       z_r = DBLE(nz_r-1)*L/DBLE(dfftp%nr3)
       !
       ! for gp=0
       DO it = 1, nat
          tt = -fpi*zv(ityp(it))/sa
          zp = tau(3,it)
          IF (zp>at(3,3)*0.5) zp = zp-at(3,3)
          zp = zp*alat
          DO iz = 1, dfftp%nr3
             k3 = iz-1
             IF (k3 > dfftp%nr3/2) k3 = iz-dfftp%nr3-1
             z = DBLE(k3)/DBLE(dfftp%nr3)*L
             ! bc3
             arg001 = -tmp**2*(z-zp)**2
             arg101 =  tmp*(z-zp)
             !
             cc1 = 0.5d0*(-(z-zp)*qe_erf(arg101)-EXP(arg001)/tmp/SQRT(pi))
             cc2 = 0.5d0*(2.d0*z1-z-zp)
             !
             vg_r(iz) = vg_r(iz)+tt*(cc1+cc2)*e2 ! factor e2: hartree -> Ry.
          ENDDO
          ! smoothing cell edge potential (avoiding unphysical oscillation)
          ! bc3
          f1 = f1 + tt*0.5d0*(-(z_r-zp)*qe_erf(tmp*(z_r-zp)) &
                    -EXP(-tmp**2*(z_r-zp)**2)/tmp/SQRT(pi))
          f2 = f2 + tt*0.5d0*(-(z_l-zp)*qe_erf(tmp*(z_l-zp)) &
                    -EXP(-tmp**2*(z_l-zp)**2)/tmp/SQRT(pi))
          f3 = f3 - tt*0.5d0*qe_erf(tmp*(z_r-zp))
          f4 = f4 - tt*0.5d0*qe_erf(tmp*(z_l-zp))
          !
          f1 = f1+tt*0.5d0*(2.d0*z1-z_r-zp)
          f2 = f2+tt*0.5d0*(2.d0*z1-z_l-zp)
          f3 = f3-tt*0.5d0
          f4 = f4-tt*0.5d0
       ENDDO
       ! for smoothing
       ! factor e2: hartree -> Ry.
       f1=f1*e2; f2=f2*e2; f3=f3*e2; f4=f4*e2
       z_r = z_r
       z_l = z_l+L
       !
       a0 = (f1*z_l**2*(z_l-3.d0*z_r)+z_r*(f3*z_l**2*(-z_l+z_r) &
              +z_r*(f2*(3.d0*z_l-z_r)+f4*z_l*(-z_l+z_r))))/(z_l-z_r)**3
       a1 = (f3*z_l**3+z_l*(6.d0*f1-6.d0*f2+(f3+2.d0*f4)*z_l)*z_r &
              -(2*f3+f4)*z_l*z_r**2-f4*z_r**3)/(z_l-z_r)**3
       a2 = (-3*f1*(z_l+z_r)+3.d0*f2*(z_l+z_r)-(z_l-z_r)*(2*f3*z_l &
              +f4*z_l+f3*z_r+2*f4*z_r))/(z_l-z_r)**3
       a3 = (2.d0*f1-2.d0*f2+(f3+f4)*(z_l-z_r))/(z_l-z_r)**3
       !
       DO iz = nz_r, nz_l
          z = DBLE(iz-1)/DBLE(dfftp%nr3)*L
          vg_r(iz) = a0+a1*z+a2*z**2+a3*z**3
       ENDDO
       !
       CALL cft_1z( vg_r, 1, dfftp%nr3, dfftp%nr3, -1, vg )
       !
       DO iz = 1, dfftp%nr3
          vloc3(iz,ng_2d) = vg(iz)
       ENDDO
       !
       DEALLOCATE( vg, vg_r )
       !
    ENDIF ! if( ng_2d > 0 )
    !
    ! Map to FFT mesh (dfftp%nrx)
    DO ng = 1, ngm
       n1 = mill(1,ng)
       n2 = mill(2,ng)
       ng_2d = imill_2d(n1,n2)
       n3 = mill(3,ng) + 1 
       IF (n3 < 1) n3 = n3 + dfftp%nr3
       aux(dfftp%nl(ng)) = aux(dfftp%nl(ng)) + vloc3(n3,ng_2d)
       IF (gamma_only) THEN
          aux(dfftp%nlm(ng)) = CONJG(aux(dfftp%nl(ng)))
       ENDIF
    ENDDO
    !
    DEALLOCATE( vloc3 )
    !
    RETURN
    !
  END SUBROUTINE esm_local_bc3
  !
  !
  !--------------------------------------------------------------------
  SUBROUTINE esm_local_bc4( aux )
    !------------------------------------------------------------------
    !! Calculates ESM local potential - bc4.
    !
    USE constants,        ONLY : pi, tpi, fpi, e2
    USE gvect,            ONLY : ngm, mill
    USE control_flags,    ONLY : gamma_only
    USE cell_base,        ONLY : at, bg, alat, tpiba2, omega
    USE ions_base,        ONLY : zv, nat, tau, ityp
    USE fft_base,         ONLY : dfftp
    USE fft_scalar,       ONLY : cft_1z
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(INOUT) :: aux(dfftp%nnr)
    !! aux contains v_loc_short(G) (input) and v_loc(G) (output)
    !
    ! ... local variables
    !
    REAL(DP) :: t(2), tt, gp, gp2, sa, z1, z0, pp, cc, ss, &
                t1, t2, z, zp, tmp, L, z_l, z_r,           &
                arg001, arg002, arg003, arg005,            &
                arg006, arg008, arg009, arg011,            &
                arg101, arg102, arg103, arg104, arg106,    &
                arg107, arg109, arg111, arg113, aaa, t3,   &
                alpha, beta, kappa, lambda, xi, chi
    INTEGER :: iz, it, k1, k2, k3, ng, n1, n2, n3, nz_l,   &
               nz_r, ng_2d
    COMPLEX(DP) :: cs, cc1, cc2, a0, a1, a2, a3, f1, f2,   &
                   f3, f4
    COMPLEX(DP), ALLOCATABLE :: vloc3(:,:), vg(:), vg_r(:)
    !
    L = at(3,3)*alat
    sa = omega/L
    z0 = L/2.d0
    tmp = 1.d0   ! Gaussian width
    z1 = z0+esm_w
    aaa = esm_a
    !
    ALLOCATE( vloc3(dfftp%nr3,ngm_2d) )
    !
    ! for gp!=0
    ALLOCATE( vg(dfftp%nr3), vg_r(dfftp%nr3) )
    !
    DO ng_2d = 1, ngm_2d
       k1 = mill_2d(1,ng_2d)
       k2 = mill_2d(2,ng_2d)
       IF (k1==0 .AND. k2==0) CYCLE
       t(1:2) = k1*bg(1:2,1)+k2*bg(1:2,2)
       gp2 = SUM(t(:)*t(:))*tpiba2
       gp  = SQRT(gp2)
       vg_r(:) = (0.d0,0.d0)
       !
       DO it = 1, nat
          tt = -fpi*zv(ityp(it))/sa
          pp = -tpi*(tau(1,it)*(k1*bg(1,1)+k2*bg(1,2)) &
                    +tau(2,it)*(k1*bg(2,1)+k2*bg(2,2)))
          cc = COS(pp)
          ss = SIN(pp)
          cs = CMPLX(cc,ss,KIND=DP)
          zp = tau(3,it)
          IF (zp > at(3,3)*0.5d0) zp = zp-at(3,3)
          zp = zp*alat
          !
          DO iz = 1, dfftp%nr3
             k3 = iz-1
             IF (k3 > dfftp%nr3/2) k3 = iz-dfftp%nr3-1
             z = DBLE(k3)/DBLE(dfftp%nr3)*L
             ! bc4
             alpha = aaa+gp+SQRT(aaa**2+gp**2)
             beta  = aaa+gp-SQRT(aaa**2+gp**2)
             kappa = aaa-gp+SQRT(aaa**2+gp**2)
             xi    = aaa   +SQRT(aaa**2+gp**2)
             chi   = aaa   -SQRT(aaa**2+gp**2)
             lambda=        SQRT(aaa**2+gp**2)
             !
             arg001= gp*(z-zp)
             arg002=-gp*(z-zp)
             arg003= gp*(z+zp-2.d0*z1)
             arg005=-gp*(z1-zp)-xi*(z-z1)
             arg006= aaa/2.d0/tmp**2*xi+gp*(z-z1)+xi*(z1-zp)
             arg008= aaa/2.d0/tmp**2*xi+xi*(z1-zp)-chi*(z-z1)
             arg009= aaa/2.d0/tmp**2*xi+xi*(z1-zp)-xi*(z-z1)
             arg011= aaa/2.d0/tmp**2*chi+chi*(z1-zp)-xi*(z-z1)
             arg101=  gp/2.d0/tmp+tmp*(z-zp)
             arg102=  gp/2.d0/tmp-tmp*(z-zp)
             arg103=  gp/2.d0/tmp+tmp*(z1-zp)
             arg104=  gp/2.d0/tmp-tmp*(z1-zp)
             arg107=  xi/2.d0/tmp+tmp*(z-zp)
             arg109=  xi/2.d0/tmp+tmp*(z1-zp)
             arg111= chi/2.d0/tmp+tmp*(z-zp)
             arg113= chi/2.d0/tmp+tmp*(z1-zp)
             !
             IF (z < z1) THEN
                t1 = exp_erfc(arg001,arg101)-exp_erfc(arg001,arg103)
                t2 = exp_erfc(arg002,arg102) &
                   -kappa/alpha*exp_erfc(arg003,arg104)
                t3 = exp_erfc(arg006,arg109)/alpha
                !
                cc1 = cs*(t1+t2)/4.d0/gp
                cc2 = cs*t3/2.d0
             ELSE
                t1 =  exp_erfc(arg011,arg113)-exp_erfc(arg011,arg111)
                t2 =  exp_erfc(arg008,arg107) &
                      -beta/alpha*exp_erfc(arg009,arg109)
                t3 =  exp_erfc(arg005,arg104)/alpha
                !
                cc1 = cs*(t1+t2)/4.d0/lambda
                cc2 = cs*t3/2.d0
             ENDIF
             !
             vg_r(iz) = vg_r(iz)+tt*(cc1+cc2)*e2 ! factor e2: hartree -> Ry.
             !
          ENDDO
          !
       ENDDO
       !
       CALL cft_1z( vg_r, 1, dfftp%nr3, dfftp%nr3, -1, vg )
       !
       DO iz = 1, dfftp%nr3
          vloc3(iz,ng_2d) = vg(iz)
       ENDDO
    ENDDO
    !
    DEALLOCATE( vg, vg_r )
    !
    ng_2d = imill_2d(0,0)
    !
    IF ( ng_2d > 0 ) THEN
       ALLOCATE( vg(dfftp%nr3), vg_r(dfftp%nr3) )
       vg_r(:) = (0.d0,0.d0)
       ! for smoothing
       f1=(0.d0,0.d0); f2=(0.d0,0.d0); f3=(0.d0,0.d0); f4=(0.d0,0.d0)
       nz_l = dfftp%nr3/2+1+esm_nfit
       nz_r = dfftp%nr3/2+1-esm_nfit
       z_l = DBLE(nz_l-1)*L/DBLE(dfftp%nr3)-L
       z_r = DBLE(nz_r-1)*L/DBLE(dfftp%nr3)
       ! for gp=0
       DO it = 1,nat
          tt = -fpi*zv(ityp(it))/sa
          zp = tau(3,it)
          IF (zp > at(3,3)*0.5) zp = zp-at(3,3)
          zp = zp*alat
          DO iz = 1, dfftp%nr3
             k3 = iz-1
             IF (k3 > dfftp%nr3/2) k3 = iz-dfftp%nr3-1
             z = DBLE(k3)/DBLE(dfftp%nr3)*L
             ! bc4
             arg001 = -tmp**2*(z-zp)**2
             arg002 = -tmp**2*(z1-zp)**2
             arg005 = -2.d0*aaa*(z-z1)
             arg006 = aaa**2/tmp**2+2.d0*aaa*(z1-zp)
             arg101 = tmp*(z-zp)
             arg102 = tmp*(z1-zp)
             arg104 = aaa/tmp+tmp*(z-zp)
             arg106 = aaa/tmp+tmp*(z1-zp)
             IF (z < z1) THEN
                t1 = -(z-zp)*qe_erf(arg101)+(0.5d0/aaa+z1-zp)*qe_erf(arg102)
                t2 = 0.5d0/aaa*exp_erfc(arg006,arg106)
                t3 = 0.5d0/aaa-(z-z1)+EXP(arg002)/tmp/SQRT(pi) &
                     -EXP(arg001)/tmp/SQRT(pi)
                !
                cc1 = (t1+t2)/2.d0
                cc2 = t3/2.d0
             ELSE
                t1 =-exp_erfc(arg005,arg101)/aaa
                t2 = exp_erfc(arg006,arg104)/aaa
                t3 = EXP(arg005)/aaa
                !
                cc1 = (t1+t2)/4.d0
                cc2 = t3/2.d0
             ENDIF
             !
             vg_r(iz) = vg_r(iz)+tt*(cc1+cc2)*e2 ! factor e2: hartree -> Ry.
             !
          ENDDO
          ! smoothing cell edge potential (avoiding unphysical oscillation)
          ! bc4
          arg002 = -tmp**2*(z1-zp)**2
          arg006 =  aaa**2/tmp**2+2.d0*aaa*(z1-zp)
          arg102 =  tmp*(z1-zp)
          arg106 =  aaa/tmp+tmp*(z1-zp)
          !-right only
          arg005 = -2.d0*aaa*(z_r-z1)
          arg101 =  tmp*(z_r-zp)
          arg104 =  aaa/tmp+tmp*(z_r-zp)
          !--
          t1 = -exp_erfc(arg005,arg101)/aaa
          t2 =  exp_erfc(arg006,arg104)/aaa
          t3 =  EXP(arg005)/aaa
          !
          f1 = f1 + tt*((t1+t2)/2.d0+t3)/2.d0
          f3 = f3 - tt*0.5d0*EXP(arg005)*(1.d0+qe_erf(arg101))
          !-left only
          arg001 = -tmp**2*(z_l-zp)**2
          arg101 =  tmp*(z_l-zp)
          !--
          t1 = -(z_l-zp)*qe_erf(arg101)+(0.5d0/aaa+z1-zp)*qe_erf(arg102)
          t2 = 0.5d0/aaa*exp_erfc(arg006,arg106)
          t3 = 0.5d0/aaa-(z_l-z1)+EXP(arg002)/tmp/SQRT(pi) &
               -EXP(arg001)/tmp/SQRT(pi)
          f2 = f2 + tt*(t1+t2+t3)/2.d0
          f4 = f4 - tt*0.5d0*(1.d0+qe_erf(arg101))
       ENDDO
       ! for smoothing
       ! factor e2: hartree -> Ry.
       f1=f1*e2; f2=f2*e2; f3=f3*e2; f4=f4*e2
       !
       z_r = z_r
       z_l = z_l+L
       !
       a0 = (f1*z_l**2*(z_l-3.d0*z_r)+z_r*(f3*z_l**2*(-z_l+z_r) &
            +z_r*(f2*(3.d0*z_l-z_r)+f4*z_l*(-z_l+z_r))))/(z_l-z_r)**3
       a1 = (f3*z_l**3+z_l*(6.d0*f1-6.d0*f2+(f3+2.d0*f4)*z_l)*z_r &
            -(2*f3+f4)*z_l*z_r**2-f4*z_r**3)/(z_l-z_r)**3
       a2 = (-3*f1*(z_l+z_r)+3.d0*f2*(z_l+z_r)-(z_l-z_r)*(2*f3*z_l &
            +f4*z_l+f3*z_r+2*f4*z_r))/(z_l-z_r)**3
       a3 = (2.d0*f1-2.d0*f2+(f3+f4)*(z_l-z_r))/(z_l-z_r)**3
       !
       DO iz = nz_r, nz_l
          z = DBLE(iz-1)/DBLE(dfftp%nr3)*L
          vg_r(iz) = (a0+a1*z+a2*z**2+a3*z**3)
       ENDDO
       !
       CALL cft_1z( vg_r, 1, dfftp%nr3, dfftp%nr3, -1, vg )
       !
       DO iz = 1, dfftp%nr3
          vloc3(iz,ng_2d) = vg(iz)
       ENDDO
       !
       DEALLOCATE( vg, vg_r )
       !
    ENDIF ! if( ng_2d > 0 )
    !
    ! Map to FFT mesh (dfftp%nrx)
    DO ng = 1, ngm
       n1 = mill(1,ng)
       n2 = mill(2,ng)
       ng_2d = imill_2d(n1,n2)
       n3 = mill(3,ng) + 1 
       IF (n3 < 1) n3 = n3 + dfftp%nr3
       aux(dfftp%nl(ng)) = aux(dfftp%nl(ng)) + vloc3(n3,ng_2d)
       IF (gamma_only) THEN
          aux(dfftp%nlm(ng)) = CONJG(aux(dfftp%nl(ng)))
       ENDIF
    ENDDO
    !
    DEALLOCATE( vloc3 )
    !
    RETURN
    !
  END SUBROUTINE esm_local_bc4
  !
  !
  !-----------------------------------------------------------------------
  SUBROUTINE esm_stres_loclong_bc1( sigmaloclong, rhog )
    !---------------------------------------------------------------------
    !! ESM local stress tensor - bc1.
    !
    USE kinds,         ONLY : DP
    USE gvect,         ONLY : ngm, mill
    USE constants,     ONLY : pi, sqrtpm1, tpi, fpi, e2
    USE cell_base,     ONLY : omega, alat, at, tpiba, bg
    USE ions_base,     ONLY : zv, nat, tau, ityp
    USE control_flags, ONLY : gamma_only
    USE fft_base,      ONLY : dfftp
    USE fft_scalar,    ONLY : cft_1z
    USE mp_bands,      ONLY : intra_bgrp_comm
    USE mp,            ONLY : mp_sum
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(OUT) :: sigmaloclong(3,3)
    !! the ESM local contribution to the stress tensor - bc1.
    COMPLEX(DP) :: rhog(ngm)
    !! n(G)
    !
    ! ... local variables
    !
    INTEGER :: ig, iga, igb, igz, igp, la, mu, iz, ia
    REAL(DP) :: L, S, z0, alpha, salp, z
    REAL(DP) :: Qa, ra(2), za
    REAL(DP) :: g(2), gp, gz
    COMPLEX(DP), PARAMETER :: ci = DCMPLX(0.0d0, 1.0d0)
    COMPLEX(DP) :: rg3
    COMPLEX(DP) :: expimgpr, experfcm, experfcp, dexperfcm_dgp, dexperfcp_dgp
    REAL(DP) :: z_r, z_l
    COMPLEX(DP) :: a0, a1, a2, a3, f1, f2, f3, f4
    COMPLEX(DP) :: poly_fr, poly_fl, poly_dfr, poly_dfl
    COMPLEX(DP) :: poly_a, poly_b, poly_c, poly_d
    REAL(DP), PARAMETER :: delta(2,2) = RESHAPE( (/ 1.0d0, 0.0d0, 0.0d0, 1.0d0 /), (/2,2/) )
    REAL(DP) :: dgp_deps(2,2)     ! dgp/deps
    REAL(DP) :: dinvgp_deps(2,2)  ! dgp^-1/deps
    !
    COMPLEX(DP), ALLOCATABLE :: rhog3(:,:)
    COMPLEX(DP), ALLOCATABLE :: dVr_deps(:,:,:)
    COMPLEX(DP), ALLOCATABLE :: dVg_deps(:,:,:)
    COMPLEX(DP), ALLOCATABLE :: Vr(:)
    COMPLEX(DP), ALLOCATABLE :: Vg(:)
    !
    ALLOCATE( rhog3(dfftp%nr3,ngm_2d) )
    ALLOCATE( dVr_deps(dfftp%nr3,2,2) )
    ALLOCATE( dVg_deps(dfftp%nr3,2,2) )
    ALLOCATE( Vr(dfftp%nr3) )
    ALLOCATE( Vg(dfftp%nr3) )
    !
    ! reconstruct rho(gz,gp)
    rhog3(:,:) = (0.d0,0.d0)
    !
    DO ig = 1, ngm
       iga = mill(1,ig)
       igb = mill(2,ig)
       igz = mill(3,ig)+1
       igp = imill_2d(iga,igb)
       !
       IF (igz < 1) THEN
          igz = igz + dfftp%nr3
       ENDIF
       rg3 = rhog(ig)
       rhog3(igz,igp) = rg3
       !
       IF ( gamma_only .AND. iga==0 .AND. igb==0 ) THEN
          igz = 1-mill(3,ig)
          if( igz<1 ) THEN
             igz = igz + dfftp%nr3
          ENDIF
          rhog3(igz,igp) = CONJG(rg3)
       ENDIF
    ENDDO ! ig
    !
    ! cell settings
    L  = at(3,3)*alat
    S  = omega/L
    z0 = L/2.d0
    alpha = 1.0d0
    salp = SQRT(alpha)
    !
    ! initialize
    sigmaloclong(:,:) = 0.0d0
    !
    !****For gp!=0 case ****
    DO igp = 1, ngm_2d
       !
       iga = mill_2d(1,igp)
       igb = mill_2d(2,igp)
       g(1:2) = (iga*bg(1:2,1) + igb*bg(1:2,2))*tpiba
       gp = SQRT(g(1)*g(1) + g(2)*g(2))
       !
       IF ( gp==0.0d0 ) CYCLE ! skip gp=0
       !
       ! derivatives by strain tensor
       DO la = 1, 2
          DO mu = 1, 2
             dgp_deps(la,mu)    = -g(la)*g(mu)/gp
             dinvgp_deps(la,mu) = +g(la)*g(mu)/gp**3
          ENDDO
       ENDDO
       !
       ! calculate dV(z)/deps
       DO iz = 1, dfftp%nr3
          !
          IF ( iz-1 > dfftp%nr3/2 ) THEN
             z = DBLE(iz-1-dfftp%nr3)/DBLE(dfftp%nr3)*L
          ELSE
             z = DBLE(iz-1)/DBLE(dfftp%nr3)*L
          ENDIF
          !
          ! summations over all atoms
          dVr_deps(iz,:,:) = ( 0.0d0, 0.0d0 )
          !
          DO ia = 1, nat
             Qa = (-1.0d0)*zv(ityp(ia))
             ra(1:2) = tau(1:2,ia)*alat
             za = tau(3,ia)*alat
             if( za > L*0.5d0 ) THEN
                za = za - L
             ENDIF
             !
             expimgpr = qe_exp( - ci*(g(1)*ra(1) + g(2)*ra(2)) )
             experfcm = exp_erfc( -gp*(z-za), gp/2.d0/salp-salp*(z-za) )
             experfcp = exp_erfc( +gp*(z-za), gp/2.d0/salp+salp*(z-za) )
             !
             dexperfcm_dgp = -(z-za)*exp_erfc( -gp*(z-za), gp/2.d0/salp-salp*(z-za) ) &
                  - EXP( -gp*(z-za) ) * qe_gauss( gp/2.d0/salp-salp*(z-za) )/2.d0/salp
             dexperfcp_dgp = +(z-za)*exp_erfc( +gp*(z-za), gp/2.d0/salp+salp*(z-za) ) &
                  - EXP( +gp*(z-za) ) * qe_gauss( gp/2.d0/salp+salp*(z-za) )/2.d0/salp
             !
             dVr_deps(iz,:,:) = dVr_deps(iz,:,:) &
                  + gp*dinvgp_deps(:,:) * pi/gp * Qa/S * expimgpr * experfcm &
                  - pi/gp * delta(:,:) * Qa/S * expimgpr * experfcm &
                  + pi/gp * Qa/S * expimgpr * dgp_deps(:,:) * dexperfcm_dgp
             !
             dVr_deps(iz,:,:) = dVr_deps(iz,:,:) &
                  + gp*dinvgp_deps(:,:) * pi/gp * Qa/S * expimgpr * experfcp &
                  - pi/gp * delta(:,:) * Qa/S * expimgpr * experfcp &
                  + pi/gp * Qa/S * expimgpr * dgp_deps(:,:) * dexperfcp_dgp
          ENDDO ! ia
       ENDDO ! iz
       !
       ! convert dV(z)/deps to dV(gz)/deps
       DO la = 1, 2
          DO mu = 1, 2
             CALL cft_1z( dVr_deps(:,la,mu), 1, dfftp%nr3, dfftp%nr3, -1, dVg_deps(:,la,mu) )
          ENDDO
       ENDDO
       !
       ! modifications
       IF ( gamma_only ) THEN
          dVg_deps(:,:,:) = dVg_deps(:,:,:)*2.0d0
       ENDIF
       !
       ! calculate stress tensor
       DO igz = 1, dfftp%nr3
          rg3 = rhog3(igz,igp)
          sigmaloclong(1:2,1:2) = sigmaloclong(1:2,1:2) &
               - DBLE( CONJG(rg3) * dVg_deps(igz,1:2,1:2) )
       ENDDO ! igz
       !
    ENDDO ! igp
    !
    !****For gp=0 case ****
    IF (imill_2d(0,0) > 0) THEN
       ! calculate V(z)
       Vr(:) = 0.0d0
       ! separation by polynomial
       f1=(0.d0,0.d0); f2=(0.d0,0.d0); f3=(0.d0,0.d0); f4=(0.d0,0.d0)
       z_l = -z0
       z_r = +z0
       !
       DO ia = 1, nat
          Qa = (-1.0d0)*zv(ityp(ia))
          ra(1:2) = tau(1:2,ia)*alat
          za = tau(3,ia)*alat
          !
          IF (za > L*0.5d0) THEN
             za = za - L
          ENDIF
          !
          DO iz = 1, dfftp%nr3
             IF ( iz-1 > dfftp%nr3/2 ) THEN
                z = DBLE(iz-1-dfftp%nr3)/DBLE(dfftp%nr3)*L
             ELSE
                z = DBLE(iz-1)/DBLE(dfftp%nr3)*L
             ENDIF
             !
             Vr(iz) = Vr(iz) - tpi * Qa/S &
                  * ( (z-za)*qe_erf(salp*(z-za)) &
                  + EXP(-alpha*(z-za)**2)*sqrtpm1/salp )
          ENDDO ! iz
          !
          f1 = f1 - tpi * Qa/S &
               * ( (z_r-za)*qe_erf(salp*(z_r-za)) &
               + EXP(-alpha*(z_r-za)**2)*sqrtpm1/salp )
          f2 = f2 - tpi * Qa/S &
               * ( (z_l-za)*qe_erf(salp*(z_l-za)) &
               + EXP(-alpha*(z_l-za)**2)*sqrtpm1/salp )
          f3 = f3 - tpi * Qa/S &
               * qe_erf(salp*(z_r-za))
          f4 = f4 - tpi * Qa/S &
               * qe_erf(salp*(z_l-za))
       ENDDO ! ia
       !
       a0 = (f1*z_l**2*(z_l-3.d0*z_r)+z_r*(f3*z_l**2*(-z_l+z_r)    &
             +z_r*(f2*(3.d0*z_l-z_r)+f4*z_l*(-z_l+z_r))))/(z_l-z_r)**3
       a1 = (f3*z_l**3+z_l*(6.d0*f1-6.d0*f2+(f3+2.d0*f4)*z_l)*z_r  &
             -(2*f3+f4)*z_l*z_r**2-f4*z_r**3)/(z_l-z_r)**3
       a2 = (-3*f1*(z_l+z_r)+3.d0*f2*(z_l+z_r)-(z_l-z_r)*(2*f3*z_l &
             +f4*z_l+f3*z_r+2*f4*z_r))/(z_l-z_r)**3
       a3 = (2.d0*f1-2.d0*f2+(f3+f4)*(z_l-z_r))/(z_l-z_r)**3
       !
       ! remove polynomial from V(z)
       DO iz = 1, dfftp%nr3
          IF ( iz-1 > dfftp%nr3/2 ) THEN
             z = DBLE(iz-1-dfftp%nr3)/DBLE(dfftp%nr3)*L
          ELSE
             z = DBLE(iz-1)/DBLE(dfftp%nr3)*L
          ENDIF
          Vr(iz) = Vr(iz) - (a0+a1*z+a2*z**2+a3*z**3)
       ENDDO
       !
       ! convert V(z) to V(gz) without polynomial
       CALL cft_1z( Vr, 1, dfftp%nr3, dfftp%nr3, -1, Vg )
       !
       ! add polynomial to V(gz)
       DO igz = 2, dfftp%nr3
          IF ( igz<=dfftp%nr3/2 )THEN
             gz = DBLE(igz-1)*tpi/L
          ELSE
             gz = DBLE(igz-1-dfftp%nr3)*tpi/L
          ENDIF
          !
          Vg(igz) = Vg(igz) &
               + a1 * ci * COS(gz*z0)/gz &
               + a2 * 2.0d0 * COS(gz*z0)/gz**2 &
               + a3 * ci * z0**2 * COS(gz*z0)/gz &
               - a3 * ci * 6.0d0 * COS(gz*z0)/gz**3
       ENDDO
       !
       Vg(1) = Vg(1) + a0*1.0d0 + a2*z0**2/3.0d0
       !
       ! calculate dV/deps(gz)
       DO igz = 1, dfftp%nr3
          dVg_deps(igz,:,:) = -delta(:,:) * Vg(igz)
       ENDDO ! igz
       !
       ! calculate stress tensor
       DO igz = 1, dfftp%nr3
          rg3 = rhog3(igz,imill_2d(0,0))
          sigmaloclong(1:2,1:2) = sigmaloclong(1:2,1:2) &
               - REAL( CONJG(rg3) * dVg_deps(igz,1:2,1:2) )
       ENDDO ! igz
       !
    ENDIF ! imill_2d(0,0) > 0
    !
    ! e2 means hartree -> Ry.
    sigmaloclong(:,:) = sigmaloclong(:,:)*(e2)
    !
    CALL mp_sum( sigmaloclong, intra_bgrp_comm )
    !
    DEALLOCATE( rhog3 )
    DEALLOCATE( dVr_deps )
    DEALLOCATE( dVg_deps )
    DEALLOCATE( Vr )
    DEALLOCATE( Vg )
    !
    RETURN
    !
  END SUBROUTINE esm_stres_loclong_bc1
  !
  !
  !---------------------------------------------------------------------------
  SUBROUTINE esm_stres_loclong_bc2( sigmaloclong, rhog )
    !-------------------------------------------------------------------------
    !! The ESM local contribution to the stress tensor - bc2.
    !
    USE kinds,         ONLY : DP
    USE gvect,         ONLY : ngm, mill
    USE constants,     ONLY : pi, sqrtpm1, tpi, fpi, e2
    USE cell_base,     ONLY : omega, alat, at, tpiba, bg
    USE ions_base,     ONLY : zv, nat, tau, ityp
    USE control_flags, ONLY : gamma_only
    USE fft_base,      ONLY : dfftp
    USE fft_scalar,    ONLY : cft_1z
    USE mp_bands,      ONLY : intra_bgrp_comm
    USE mp,            ONLY : mp_sum
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(OUT) :: sigmaloclong(3,3)
    !! contribution to the stress tensor
    COMPLEX(DP) :: rhog(ngm)   
    !! n(G)
    !
    ! ... local variables
    !
    INTEGER :: ig, iga, igb, igz, igp, la, mu, iz, ia
    REAL(DP) :: L, S, z0, z1, alpha, salp, z
    REAL(DP) :: Qa, ra(2), za
    REAL(DP) :: g(2), gp, gz
    COMPLEX(DP), PARAMETER :: ci = DCMPLX(0.0d0, 1.0d0)
    COMPLEX(DP) :: rg3
    COMPLEX(DP) :: expimgpr, experfcm, experfcp, dexperfcm_dgp, dexperfcp_dgp
    COMPLEX(DP) :: exph1, exph2, exph3
    REAL(DP) :: z_r, z_l
    COMPLEX(DP) :: a0, a1, a2, a3, f1, f2, f3, f4
    REAL(DP), PARAMETER :: delta(2,2) = RESHAPE( (/ 1.0d0, 0.0d0, 0.0d0, 1.0d0 /), (/2,2/) )
    REAL(DP) :: dgp_deps(2,2)     ! dgp/deps
    REAL(DP) :: dinvgp_deps(2,2)  ! dgp^-1/deps
    !
    COMPLEX(DP), ALLOCATABLE :: rhog3(:,:)
    COMPLEX(DP), ALLOCATABLE :: dVr_deps(:,:,:)
    COMPLEX(DP), ALLOCATABLE :: dVg_deps(:,:,:)
    COMPLEX(DP), ALLOCATABLE :: Vr(:)
    COMPLEX(DP), ALLOCATABLE :: Vg(:)
    !
    ALLOCATE( rhog3(dfftp%nr3,ngm_2d) )
    ALLOCATE( dVr_deps(dfftp%nr3,2,2) )
    ALLOCATE( dVg_deps(dfftp%nr3,2,2) )
    ALLOCATE( Vr(dfftp%nr3) )
    ALLOCATE( Vg(dfftp%nr3) )
    !
    ! reconstruct rho(gz,gp)
    rhog3(:,:) = (0.d0,0.d0)
    !
    DO ig = 1, ngm
       iga = mill(1,ig)
       igb = mill(2,ig)
       igz = mill(3,ig)+1
       igp = imill_2d(iga,igb)
       IF (igz < 1) THEN
          igz = igz + dfftp%nr3
       ENDIF
       !
       rg3 = rhog(ig)
       rhog3(igz,igp) = rg3
       !
       IF ( gamma_only .AND. iga==0 .AND. igb==0 ) THEN
          igz = 1-mill(3,ig)
          IF ( igz<1 ) THEN
             igz = igz + dfftp%nr3
          ENDIF
          rhog3(igz,igp) = CONJG(rg3)
       ENDIF
    ENDDO ! ig
    !
    ! cell settings
    L  = at(3,3)*alat
    S  = omega/L
    z0 = L/2.d0
    z1 = z0+esm_w
    alpha = 1.0d0
    salp = SQRT(alpha)
    !
    ! initialize
    sigmaloclong(:,:) = 0.0d0
    !
    !****For gp!=0 case ****
    DO igp = 1, ngm_2d
       iga = mill_2d(1,igp)
       igb = mill_2d(2,igp)
       g(1:2) = (iga*bg(1:2,1) + igb*bg(1:2,2))*tpiba
       gp = SQRT(g(1)*g(1) + g(2)*g(2))
       !
       IF ( gp==0.0d0 ) CYCLE ! skip gp=0
       !
       ! derivatives by strain tensor
       DO la = 1, 2
          DO mu = 1, 2
             dgp_deps(la,mu) = -g(la)*g(mu)/gp
             dinvgp_deps(la,mu) = +g(la)*g(mu)/gp**3
          ENDDO
       ENDDO
       !
       ! calculate dV(z)/deps
       DO iz = 1, dfftp%nr3
          IF ( iz-1 > dfftp%nr3/2 ) THEN
             z = DBLE(iz-1-dfftp%nr3)/DBLE(dfftp%nr3)*L
          ELSE
             z = DBLE(iz-1)/DBLE(dfftp%nr3)*L
          ENDIF
          !
          ! summations over all atoms
          dVr_deps(iz,:,:) = ( 0.0d0, 0.0d0 )
          !
          DO ia = 1, nat
             Qa = (-1.0d0)*zv(ityp(ia))
             ra(1:2) = tau(1:2,ia)*alat
             za = tau(3,ia)*alat
             IF ( za > L*0.5d0 ) THEN
                za = za - L
             ENDIF
             !
             expimgpr = qe_EXP( - ci*(g(1)*ra(1)+g(2)*ra(2)) )
             experfcm = exp_erfc( -gp*(z-za), gp/2.d0/salp-salp*(z-za) )
             experfcp = exp_erfc( +gp*(z-za), gp/2.d0/salp+salp*(z-za) )
             !
             dexperfcm_dgp = -(z-za)*exp_erfc( -gp*(z-za), gp/2.d0/salp-salp*(z-za) ) &
                  - EXP( -gp*(z-za) ) * qe_gauss( gp/2.d0/salp-salp*(z-za) )/2.d0/salp
             dexperfcp_dgp = +(z-za)*exp_erfc( +gp*(z-za), gp/2.d0/salp+salp*(z-za) ) &
                  - EXP( +gp*(z-za) ) * qe_gauss( gp/2.d0/salp+salp*(z-za) )/2.d0/salp
             !
             exph1  = (COSH(gp*(z-za))*EXP(-2*gp*z1) - COSH(gp*(z+za)) )/SINH(2*gp*z1)
             exph2  = ( (z-za)*SINH(gp*(z-za))*EXP(-2*gp*z1) &
                        - 2*z1*COSH(gp*(z-za))*EXP(-2*gp*z1) &
                        -(z+za)*SINH(gp*(z+za)) )/SINH(2*gp*z1)
             exph3  = - (COSH(gp*(z-za))*EXP(-2*gp*z1) - COSH(gp*(z+za)) )/SINH(2*gp*z1)**2 &
                        *2*z1*COSH(2*gp*z1)
             !
             ! BC1 terms
             dVr_deps(iz,:,:) = dVr_deps(iz,:,:) &
                  + gp*dinvgp_deps(:,:) * pi/gp * Qa/S * expimgpr * experfcm &
                  - pi/gp * delta(:,:) * Qa/S * expimgpr * experfcm &
                  + pi/gp * Qa/S * expimgpr * dgp_deps(:,:) * dexperfcm_dgp
             !
             ! BC1 terms
             dVr_deps(iz,:,:) = dVr_deps(iz,:,:) &
                  + gp*dinvgp_deps(:,:) * pi/gp * Qa/S * expimgpr * experfcp &
                  - pi/gp * delta(:,:) * Qa/S * expimgpr * experfcp &
                  + pi/gp * Qa/S * expimgpr * dgp_deps(:,:) * dexperfcp_dgp
             !
             ! BC2 terms
             dVr_deps(iz,:,:) = dVr_deps(iz,:,:) &
                  + gp*dinvgp_deps(:,:) * tpi/gp * Qa/S * expimgpr * exph1 &
                  - tpi/gp * delta(:,:) * Qa/S * expimgpr * exph1 &
                  + tpi/gp * Qa/S * expimgpr * dgp_deps(:,:) * (exph2+exph3)
          ENDDO ! ia
          !
       ENDDO ! iz
       !
       ! convert dV(z)/deps to dV(gz)/deps
       DO la = 1, 2
          DO mu = 1, 2
             CALL cft_1z( dVr_deps(:,la,mu), 1, dfftp%nr3, dfftp%nr3, -1, dVg_deps(:,la,mu) )
          ENDDO
       ENDDO
       !
       ! modifications
       IF ( gamma_only ) THEN
          dVg_deps(:,:,:) = dVg_deps(:,:,:)*2.0d0
       ENDIF
       !
       ! calculate stress tensor
       DO igz = 1, dfftp%nr3
          rg3 = rhog3(igz,igp)
          sigmaloclong(1:2,1:2) = sigmaloclong(1:2,1:2) &
               - DBLE( CONJG(rg3) * dVg_deps(igz,1:2,1:2) )
       ENDDO ! igz
       !
    ENDDO ! igp
    !
    !****For gp=0 case ****
    IF ( imill_2d(0,0) > 0 ) THEN
       ! calculate V(z)
       ! [note] this Vr does not contain a term due to efield z*efield
       ! because it vanishes in the differentiation with respect to strain.
       Vr(:) = 0.0d0
       ! separation by polynomial
       f1=(0.d0,0.d0); f2=(0.d0,0.d0); f3=(0.d0,0.d0); f4=(0.d0,0.d0)
       !
       z_l = -z0
       z_r = +z0
       !
       DO ia = 1, nat
          !
          Qa = (-1.0d0)*zv(ityp(ia))
          ra(1:2) = tau(1:2,ia)*alat
          za = tau(3,ia)*alat
          IF ( za > L*0.5d0 ) THEN
             za = za - L
          ENDIF
          !
          DO iz = 1, dfftp%nr3
             IF (iz-1 > dfftp%nr3/2) THEN
                z = DBLE(iz-1-dfftp%nr3)/DBLE(dfftp%nr3)*L
             ELSE
                z = DBLE(iz-1)/DBLE(dfftp%nr3)*L
             ENDIF
             !
             ! BC1 terms
             Vr(iz) = Vr(iz) - tpi * Qa/S &
                  * ( (z-za)*qe_erf(salp*(z-za)) &
                  + EXP(-alpha*(z-za)**2)*sqrtpm1/salp )
             !
             ! BC2 terms
             Vr(iz) = Vr(iz) + tpi * Qa/S * ( -z*za + z1*z1 )/z1
          ENDDO ! iz
          !
          f1 = f1 - tpi * Qa/S &
               * ( (z_r-za)*qe_erf(salp*(z_r-za)) &
               + EXP(-alpha*(z_r-za)**2)*sqrtpm1/salp )
          f1 = f1 + tpi * Qa/S * ( -z_r*za + z1*z1 )/z1
          !
          f2 = f2 - tpi * Qa/S &
               * ( (z_l-za)*qe_erf(salp*(z_l-za)) &
               + EXP(-alpha*(z_l-za)**2)*sqrtpm1/salp )
          f2 = f2 + tpi * Qa/S * ( -z_l*za + z1*z1 )/z1
          !
          f3 = f3 - tpi * Qa/S &
               * qe_erf(salp*(z_r-za))
          f3 = f3 + tpi * Qa/S * ( -za )/z1
          !
          f4 = f4 - tpi * Qa/S &
               * qe_erf(salp*(z_l-za))
          f4 = f4 + tpi * Qa/S * ( -za )/z1
          !
       ENDDO ! ia
       !
       a0 = (f1*z_l**2*(z_l-3.d0*z_r)+z_r*(f3*z_l**2*(-z_l+z_r)    &
            +z_r*(f2*(3.d0*z_l-z_r)+f4*z_l*(-z_l+z_r))))/(z_l-z_r)**3
       a1 = (f3*z_l**3+z_l*(6.d0*f1-6.d0*f2+(f3+2.d0*f4)*z_l)*z_r  &
            -(2*f3+f4)*z_l*z_r**2-f4*z_r**3)/(z_l-z_r)**3
       a2 = (-3*f1*(z_l+z_r)+3.d0*f2*(z_l+z_r)-(z_l-z_r)*(2*f3*z_l &
            +f4*z_l+f3*z_r+2*f4*z_r))/(z_l-z_r)**3
       a3 = (2.d0*f1-2.d0*f2+(f3+f4)*(z_l-z_r))/(z_l-z_r)**3
       !
       ! remove polynomial from V(z) 
       DO iz=1, dfftp%nr3
          IF ( iz-1 > dfftp%nr3/2 ) THEN
             z = DBLE(iz-1-dfftp%nr3)/DBLE(dfftp%nr3)*L
          ELSE
             z = DBLE(iz-1)/DBLE(dfftp%nr3)*L
          ENDIF
          Vr(iz) = Vr(iz) - (a0+a1*z+a2*z**2+a3*z**3)
       ENDDO
       !
       ! convert V(z) to V(gz) without polynomial
       CALL cft_1z( Vr, 1, dfftp%nr3, dfftp%nr3, -1, Vg )
       !
       ! add polynomial to V(gz)
       DO igz = 2, dfftp%nr3
          IF ( igz <= dfftp%nr3/2 )THEN
             gz = DBLE(igz-1)*tpi/L
          ELSE
             gz = DBLE(igz-1-dfftp%nr3)*tpi/L
          ENDIF
          !
          Vg(igz) = Vg(igz) &
               + a1 * ci * COS(gz*z0)/gz &
               + a2 * 2.0d0 * COS(gz*z0)/gz**2   &
               + a3 * ci * z0**2 * COS(gz*z0)/gz &
               - a3 * ci * 6.0d0 * COS(gz*z0)/gz**3
       ENDDO
       Vg(1) = Vg(1) + a0*1.0d0 + a2*z0**2/3.0d0
       !
       ! calculate dV/deps(gz)
       DO igz = 1, dfftp%nr3
          dVg_deps(igz,:,:) = -delta(:,:) * Vg(igz)
       ENDDO ! igz
       !
       ! calculate stress tensor
       DO igz = 1, dfftp%nr3
          rg3 = rhog3(igz,imill_2d(0,0))
          sigmaloclong(1:2,1:2) = sigmaloclong(1:2,1:2) &
               - DBLE( CONJG(rg3) * dVg_deps(igz,1:2,1:2) )
       ENDDO ! igz
    ENDIF ! imill_2d(0,0) > 0
    !
    ! e2 means hartree -> Ry.
    sigmaloclong(:,:) = sigmaloclong(:,:)*(e2)
    !
    CALL mp_sum( sigmaloclong, intra_bgrp_comm )
    !
    DEALLOCATE( rhog3 )
    DEALLOCATE( dVr_deps )
    DEALLOCATE( dVg_deps )
    DEALLOCATE( Vr )
    DEALLOCATE( Vg )
    !
    RETURN
    !
  END SUBROUTINE esm_stres_loclong_bc2
  !
  !
  !----------------------------------------------------------------------------------
  SUBROUTINE esm_stres_loclong_bc3( sigmaloclong, rhog )
    !--------------------------------------------------------------------------------
    !! The ESM local contribution to the stress tensor - bc3.
    !
    USE kinds,           ONLY : DP
    USE gvect,           ONLY : ngm, mill
    USE constants,       ONLY : pi, sqrtpm1, tpi, fpi, e2
    USE cell_base,       ONLY : omega, alat, at, tpiba, bg
    USE ions_base,       ONLY : zv, nat, tau, ityp
    USE control_flags,   ONLY : gamma_only
    USE fft_base,        ONLY : dfftp
    USE fft_scalar,      ONLY : cft_1z
    USE mp_bands,        ONLY : intra_bgrp_comm
    USE mp,              ONLY : mp_sum
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(OUT) :: sigmaloclong(3,3)
    !! contribution to the stress tensor
    COMPLEX(DP) :: rhog(ngm)
    !! n(G)
    !
    ! ... local variables
    !
    INTEGER :: ig, iga, igb, igz, igp, la, mu, iz, ia
    REAL(DP) :: L, S, z0, z1, alpha, salp, z
    REAL(DP) :: Qa, ra(2), za
    REAL(DP) :: g(2), gp, gz
    COMPLEX(DP), PARAMETER :: ci = DCMPLX(0.0d0, 1.0d0)
    COMPLEX(DP) :: rg3
    COMPLEX(DP) :: expimgpr, experfcm, experfcp, dexperfcm_dgp, dexperfcp_dgp
    COMPLEX(DP) :: expm
    REAL(DP) :: z_r, z_l
    COMPLEX(DP) :: a0, a1, a2, a3, f1, f2, f3, f4
    REAL(DP), PARAMETER :: delta(2,2) = RESHAPE( (/ 1.0d0, 0.0d0, 0.0d0, 1.0d0 /), (/2,2/) )
    REAL(DP) :: dgp_deps(2,2)     ! dgp/deps
    REAL(DP) :: dinvgp_deps(2,2)  ! dgp^-1/deps
    !
    COMPLEX(DP), ALLOCATABLE :: rhog3(:,:)
    COMPLEX(DP), ALLOCATABLE :: dVr_deps(:,:,:)
    COMPLEX(DP), ALLOCATABLE :: dVg_deps(:,:,:)
    COMPLEX(DP), ALLOCATABLE :: Vr(:)
    COMPLEX(DP), ALLOCATABLE :: Vg(:)
    !
    REAL(DP) :: sigmaloclong_bc1(3,3)
    !
    ALLOCATE( rhog3(dfftp%nr3,ngm_2d) )
    ALLOCATE( dVr_deps(dfftp%nr3,2,2) )
    ALLOCATE( dVg_deps(dfftp%nr3,2,2) )
    ALLOCATE( Vr(dfftp%nr3) )
    ALLOCATE( Vg(dfftp%nr3) )
    !
    ! reconstruct rho(gz,gp)
    rhog3(:,:) = (0.d0,0.d0)
    !
    DO ig = 1, ngm
       iga = mill(1,ig)
       igb = mill(2,ig)
       igz = mill(3,ig)+1
       igp = imill_2d(iga,igb)
       IF (igz < 1) THEN
          igz = igz + dfftp%nr3
       ENDIF
       rg3 = rhog(ig)
       rhog3(igz,igp) = rg3
       !
       IF ( gamma_only .AND. iga==0 .AND. igb==0 ) THEN
          igz = 1-mill(3,ig)
          IF ( igz<1 ) THEN
             igz = igz + dfftp%nr3
          ENDIF
          rhog3(igz,igp) = CONJG(rg3)
       ENDIF
    ENDDO ! ig
    !
    ! cell settings
    L  = at(3,3)*alat
    S  = omega/L
    z0 = L/2.d0
    z1 = z0+esm_w
    alpha = 1.0d0
    salp = SQRT(alpha)
    !
    ! initialize
    sigmaloclong(:,:) = 0.0d0
    !
    !****For gp!=0 case ****
    DO igp = 1, ngm_2d
       iga = mill_2d(1,igp)
       igb = mill_2d(2,igp)
       g(1:2) = (iga*bg(1:2,1) + igb*bg(1:2,2))*tpiba
       gp = SQRT(g(1)*g(1) + g(2)*g(2))
       !
       IF (gp == 0.0d0) CYCLE ! skip gp=0
       !
       ! derivatives by strain tensor
       DO la = 1, 2
          DO mu = 1, 2
             dgp_deps(la,mu) = -g(la)*g(mu)/gp
             dinvgp_deps(la,mu) = +g(la)*g(mu)/gp**3
          ENDDO
       ENDDO
       !
       ! calculate dV(z)/deps
       DO iz = 1, dfftp%nr3
          IF ( iz-1 > dfftp%nr3/2 ) THEN
             z = DBLE(iz-1-dfftp%nr3)/DBLE(dfftp%nr3)*L
          ELSE
             z = DBLE(iz-1)/DBLE(dfftp%nr3)*L
          ENDIF
          !
          ! summations over all atoms
          dVr_deps(iz,:,:) = ( 0.0d0, 0.0d0 )
          !
          DO ia = 1, nat
             Qa = (-1.0d0)*zv(ityp(ia))
             ra(1:2) = tau(1:2,ia)*alat
             za = tau(3,ia)*alat
             IF ( za > L*0.5d0 ) THEN
                za = za - L
             ENDIF
             !
             expimgpr = qe_exp( - ci*(g(1)*ra(1)+g(2)*ra(2)) )
             experfcm = exp_erfc( -gp*(z-za), gp/2.d0/salp-salp*(z-za) )
             experfcp = exp_erfc( +gp*(z-za), gp/2.d0/salp+salp*(z-za) )
             !
             dexperfcm_dgp = -(z-za)*exp_erfc( -gp*(z-za), gp/2.d0/salp-salp*(z-za) ) &
                  - EXP( -gp*(z-za) ) * qe_gauss( gp/2.d0/salp-salp*(z-za) )/2.d0/salp
             dexperfcp_dgp = +(z-za)*exp_erfc( +gp*(z-za), gp/2.d0/salp+salp*(z-za) ) &
                  - EXP( +gp*(z-za) ) * qe_gauss( gp/2.d0/salp+salp*(z-za) )/2.d0/salp
             !
             expm = EXP( -gp*(-z+2*z1-za) )
             !
             ! BC1 terms
             dVr_deps(iz,:,:) = dVr_deps(iz,:,:) &
                  + gp*dinvgp_deps(:,:) * pi/gp * Qa/S * expimgpr * experfcm &
                  - pi/gp * delta(:,:) * Qa/S * expimgpr * experfcm &
                  + pi/gp * Qa/S * expimgpr * dgp_deps(:,:) * dexperfcm_dgp
             !
             ! BC1 terms
             dVr_deps(iz,:,:) = dVr_deps(iz,:,:) &
                  + gp*dinvgp_deps(:,:) * pi/gp * Qa/S * expimgpr * experfcp &
                  - pi/gp * delta(:,:) * Qa/S * expimgpr * experfcp &
                  + pi/gp * Qa/S * expimgpr * dgp_deps(:,:) * dexperfcp_dgp

             ! BC3 terms
             dVr_deps(iz,:,:) = dVr_deps(iz,:,:) &
                  - gp*dinvgp_deps(:,:) * tpi/gp * Qa/S * expimgpr * expm &
                  + tpi/gp * delta(:,:) * Qa/S * expimgpr * expm &
                  + tpi/gp * Qa/S * expimgpr * dgp_deps(:,:) * (-z+2*z1-za) * expm
          ENDDO ! ia
          !
       ENDDO ! iz
       !
       ! convert dV(z)/deps to dV(gz)/deps
       DO la = 1, 2
          DO mu = 1, 2
             CALL cft_1z( dVr_deps(:,la,mu), 1, dfftp%nr3, dfftp%nr3, -1, dVg_deps(:,la,mu) )
          ENDDO
       ENDDO
       !
       ! modifications
       IF ( gamma_only ) THEN
          dVg_deps(:,:,:) = dVg_deps(:,:,:)*2.0d0
       ENDIF
       !
       ! calculate stress tensor
       DO igz = 1, dfftp%nr3
          rg3 = rhog3(igz,igp)
          sigmaloclong(1:2,1:2) = sigmaloclong(1:2,1:2) &
             - DBLE( CONJG(rg3) * dVg_deps(igz,1:2,1:2) )
       ENDDO ! igz
       !
    ENDDO ! igp
    !
    !****For gp=0 case ****
    IF (imill_2d(0,0) > 0) THEN
       ! calculate V(z)
       Vr(:) = 0.0d0
       ! separation by polynomial
       f1=(0.d0,0.d0); f2=(0.d0,0.d0); f3=(0.d0,0.d0); f4=(0.d0,0.d0)
       !
       z_l = -z0
       z_r = +z0
       !
       DO ia = 1, nat
          Qa = (-1.0d0)*zv(ityp(ia))
          ra(1:2) = tau(1:2,ia)*alat
          za = tau(3,ia)*alat
          IF (za > L*0.5d0) THEN
             za = za - L
          ENDIF
          !
          DO iz = 1, dfftp%nr3
             !
             IF ( iz-1 > dfftp%nr3/2 ) THEN
                z = DBLE(iz-1-dfftp%nr3)/DBLE(dfftp%nr3)*L
             ELSE
                z = DBLE(iz-1)/DBLE(dfftp%nr3)*L
             ENDIF
             !
             ! BC1 terms
             Vr(iz) = Vr(iz) - tpi * Qa/S &
                  * ( (z-za)*qe_erf(salp*(z-za)) &
                  + EXP(-alpha*(z-za)**2)*sqrtpm1/salp )
             !
             ! BC3 terms
             Vr(iz) = Vr(iz) + tpi * Qa/S * (-z+2*z1-za)
          ENDDO ! iz
          !
          f1 = f1 - tpi * Qa/S &
               * ( (z_r-za)*qe_erf(salp*(z_r-za)) &
               + EXP(-alpha*(z_r-za)**2)*sqrtpm1/salp )
          f1 = f1 + tpi * Qa/S * (-z_r+2*z1-za)
          !
          f2 = f2 - tpi * Qa/S &
               * ( (z_l-za)*qe_erf(salp*(z_l-za)) &
               + EXP(-alpha*(z_l-za)**2)*sqrtpm1/salp )
          f2 = f2 + tpi * Qa/S * (-z_l+2*z1-za)
          !
          f3 = f3 - tpi * Qa/S &
               * qe_erf(salp*(z_r-za))
          f3 = f3 + tpi * Qa/S * (-1.0d0)
          !
          f4 = f4 - tpi * Qa/S &
               * qe_erf(salp*(z_l-za))
          f4 = f4 + tpi * Qa/S * (-1.0d0)
          !
       ENDDO ! ia
       !
       a0 = (f1*z_l**2*(z_l-3.d0*z_r)+z_r*(f3*z_l**2*(-z_l+z_r) &
             +z_r*(f2*(3.d0*z_l-z_r)+f4*z_l*(-z_l+z_r))))/(z_l-z_r)**3
       a1 = (f3*z_l**3+z_l*(6.d0*f1-6.d0*f2+(f3+2.d0*f4)*z_l)*z_r &
             -(2*f3+f4)*z_l*z_r**2-f4*z_r**3)/(z_l-z_r)**3
       a2 = (-3*f1*(z_l+z_r)+3.d0*f2*(z_l+z_r)-(z_l-z_r)*(2*f3*z_l &
             +f4*z_l+f3*z_r+2*f4*z_r))/(z_l-z_r)**3
       a3 = (2.d0*f1-2.d0*f2+(f3+f4)*(z_l-z_r))/(z_l-z_r)**3
       !
       ! remove polynomial from V(z)
       DO iz = 1, dfftp%nr3
          IF ( iz-1 > dfftp%nr3/2 ) THEN
             z = DBLE(iz-1-dfftp%nr3)/DBLE(dfftp%nr3)*L
          ELSE
             z = DBLE(iz-1)/DBLE(dfftp%nr3)*L
          ENDIF
          Vr(iz) = Vr(iz) - (a0+a1*z+a2*z**2+a3*z**3)
       ENDDO
       !
       ! convert V(z) to V(gz) without polynomial
       CALL cft_1z( Vr, 1, dfftp%nr3, dfftp%nr3, -1, Vg )
       !
       ! add polynomial to V(gz)
       DO igz = 2, dfftp%nr3
          IF (igz <= dfftp%nr3/2)THEN
             gz = DBLE(igz-1)*tpi/L
          ELSE
             gz = DBLE(igz-1-dfftp%nr3)*tpi/L
          ENDIF
          !
          Vg(igz) = Vg(igz) &
               + a1 * ci * COS(gz*z0)/gz &
               + a2 * 2.0d0 * COS(gz*z0)/gz**2 &
               + a3 * ci * z0**2 * COS(gz*z0)/gz &
               - a3 * ci * 6.0d0 * COS(gz*z0)/gz**3
       ENDDO
       !
       Vg(1) = Vg(1) + a0*1.0d0 + a2*z0**2/3.0d0
       !
       ! calculate dV/deps(gz)
       DO igz = 1, dfftp%nr3
          dVg_deps(igz,:,:) = -delta(:,:) * Vg(igz)
       ENDDO ! igz
       !
       ! calculate stress tensor
       DO igz = 1, dfftp%nr3
          rg3 = rhog3(igz,imill_2d(0,0))
          sigmaloclong(1:2,1:2) = sigmaloclong(1:2,1:2) &
               - DBLE( CONJG(rg3) * dVg_deps(igz,1:2,1:2) )
       ENDDO ! igz
       !
    ENDIF ! imill_2d(0,0) > 0
    !
    ! e2 means hartree -> Ry.
    sigmaloclong(:,:) = sigmaloclong(:,:)*(e2)
    !
    CALL mp_sum( sigmaloclong, intra_bgrp_comm )
    !
    DEALLOCATE( rhog3 )
    DEALLOCATE( dVr_deps )
    DEALLOCATE( dVg_deps )
    DEALLOCATE( Vr )
    DEALLOCATE( Vg )
    !
    RETURN
    !
  END SUBROUTINE esm_stres_loclong_bc3
  !
  !
  !----------------------------------------------------------------------
  SUBROUTINE esm_force_ew( forceion )
    !-----------------------------------------------------------------------
    !! This routine computes the Ewald contribution to the forces,
    !! both the real- and reciprocal-space terms are present.
    !
    USE kinds
    USE constants,   ONLY : tpi, e2
    USE mp_bands,    ONLY : intra_bgrp_comm
    USE mp,          ONLY : mp_sum
    USE ions_base,   ONLY : zv, nat, ityp
    USE gvect,       ONLY : gcutm
    USE cell_base,   ONLY : tpiba2
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(OUT) :: forceion( 3, nat )
    !! the ewald part of the forces
    !
    ! ... local variables
    !
    REAL(DP) :: alpha, charge, upperbound
    ! the alpha parameter
    ! the total charge
    ! used to determine alpha
    !
    !
    forceion(:,:) = 0.d0
    charge = SUM( zv( ityp(1:nat) ) )
    !
    ! choose alpha in order to have convergence in the sum over G
    ! upperbound is a safe upper bound for the error ON THE ENERGY
    !
    alpha = 2.9d0
    !
    DO
       alpha = alpha - 0.1d0
       IF (alpha == 0.d0) THEN
          CALL errore( 'esm_force_ew', 'optimal alpha not found', 1 )
       ENDIF
       upperbound = e2 * charge**2 * SQRT(2.d0 * alpha / tpi) * &
          qe_erfc ( SQRT(tpiba2 * gcutm / 4.d0 / alpha) )
       IF (upperbound < 1.0d-7) EXIT
    ENDDO
    !
    !write(*,'(5X,A,F5.2)')'alpha used in esm ewald force :',alpha
    !
    CALL esm_force_ewg( alpha, forceion )
    !
    CALL esm_force_ewr( alpha, forceion )
    !
    CALL mp_sum( forceion, intra_bgrp_comm )
    !
    RETURN
    !
  END SUBROUTINE esm_force_ew
  !
  !
  !  ...  ESM EWALD-DERIVED FORCE (RSUM) SUBROUTINE 
  !
  !-----------------------------------------------------------------------
  SUBROUTINE esm_force_ewr_pbc( alpha_g, forceion )
    !---------------------------------------------------------------------
    !! ESM Ewald-derived force (R-sum) - pbc.
    !
    USE constants,        ONLY : pi, e2
    USE cell_base,        ONLY : alat, at, bg
    USE ions_base,        ONLY : zv, nat, tau, ityp
    USE mp,               ONLY : mp_rank, mp_size
    USE mp_bands,         ONLY : intra_bgrp_comm
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN) :: alpha_g
    !! alpha term in Ewald sum
    REAL(DP), INTENT(INOUT) :: forceion(3,nat)
    !! Ewald derived force
    !
    ! ... local variables
    !
    INTEGER :: na, nb, nr, nrm, ip, np
    ! counter on atoms
    ! counter on atoms
    ! counter over direct vectors
    ! number of R vectors included in r sum
    INTEGER, PARAMETER :: mxr = 1000
    ! the maximum number of R vectors included in r
    REAL(DP) :: dtau(3), r(3,mxr), r2(mxr)
    ! the difference tau_s - tau_s'
    ! neighbering shell vector
    ! the square modulus of R_j-tau_s-tau_s'
    !
    ! ESM variables
    REAL(DP) :: tmp, fac, rmax0, rr
    ! rmax0: the maximum radius to consider real space sum
    REAL(DP), ALLOCATABLE :: force(:,:)
    !
    tmp = SQRT(alpha_g)
    rmax0 = 5.d0 / tmp / alat
    !
    ip = mp_rank( intra_bgrp_comm )
    np = mp_size( intra_bgrp_comm )
    !
    ALLOCATE( force(3,nat) )
    !
    force(:,:) = 0.d0
    !
    DO na = ip+1, nat, np
       DO nb = 1, nat
          IF (nb == na) CYCLE
          dtau(:) = tau(:,na)-tau(:,nb)
          fac = zv(ityp(na))*zv(ityp(nb))*e2
          !
          ! generates nearest-neighbors shells r(i)=R(i)-dtau(i)
          !
          CALL rgen( dtau, rmax0, mxr, at, bg, r, r2, nrm )
          !
          ! and sum to the real space part
          !
          DO nr = 1, nrm
             rr = SQRT(r2(nr))*alat
             force(:,na) = force(:,na) &
                -fac/rr**2*(qe_erfc(tmp*rr)/rr+2.d0*tmp/SQRT(pi) &
                *EXP(-tmp**2*rr**2))*r(:,nr)*alat
          ENDDO
          !
       ENDDO
    ENDDO
    !
    forceion(:,:) = forceion(:,:) + force(:,:)
    !
    DEALLOCATE( force )
    !
  END SUBROUTINE esm_force_ewr_pbc
  !
  !
  !-----------------------------------------------------------------------
  SUBROUTINE esm_force_ewr_bc4( alpha_g, forceion )
    !---------------------------------------------------------------------
    !! ESM Ewald-derived force (R-sum) - bc4.
    !
    USE io_global,        ONLY : stdout
    USE constants,        ONLY : pi, tpi, fpi, e2
    USE gvect,            ONLY : gstart
    USE cell_base,        ONLY : alat, tpiba2, at, bg
    USE ions_base,        ONLY : zv, nat, tau, ityp
    USE control_flags,    ONLY : iverbosity
    USE mp,               ONLY : mp_rank, mp_size
    USE mp_bands,         ONLY : intra_bgrp_comm
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN)   :: alpha_g
    !! alpha term in ewald sum
    REAL(DP), INTENT(INOUT) :: forceion(3,nat)
    !! ESM Ewald-derived force
    !
    ! ... local variables
    !
    INTEGER :: na, nb, nr, nrm, ipol, ip, np
    ! counter on atoms
    ! counter on atoms
    ! counter over direct vectors
    ! number of R vectors included in r sum
    INTEGER, PARAMETER :: mxr = 1000
    ! the maximum number of R vectors included in r
    REAL(DP) :: dtau(3), r(3,mxr), r2(mxr), rxy, rxyz
    ! the difference tau_s - tau_s'
    ! neighbering shell vector
    ! the square modulus of R_j-tau_s-tau_s'
    ! buffer variable
    ! buffer variable
    !
    ! ESM variables
    REAL(DP) :: L, z, zp, z0, z1, aaa, tmp, ss, fac, err, ss0, &
                gpmax, rmax0, rmax, zbuff, znrm, rr
    ! gpmax: upper bound of g_parallel integral
    ! rmax: the maximum radius to consider real space sum
    ! zbuff: smearing width to avoid the singularity of the Force
    ! znrm: threashold value for normal RSUM and Smooth-ESM's RSUM
    REAL(DP), PARAMETER :: eps=1.d-11, epsneib=1.d-6
    REAL(DP), ALLOCATABLE :: force(:,:)
    !
    L = at(3,3)*alat
    z0 = L/2.d0
    z1 = z0+esm_w
    aaa = esm_a
    tmp = SQRT(alpha_g)
    zbuff = 1.d0
    !
    ! Define upperbound for g_parallel integral
    err=1.d0; ss0=0.d0; gpmax=1.d0
    !
    DO
       gpmax = gpmax+1.d0
       IF (gpmax > 1000.d0) &
          CALL errore( 'esm_force_ewr', 'optimal gpmax not found', 1 )
       CALL qromb( vl11, aaa, tmp, z1, z1-zbuff, z1-zbuff, 0.0_DP, gpmax, ss )
       err = ABS(ss-ss0); ss0 = ss
       IF (err < eps) EXIT
    ENDDO
    !
    ! Define znrm using the deviation from the constant term in RSUM
    znrm = z1
    !
    DO
       znrm = znrm-0.01d0
       IF (znrm <= -z0) &
          CALL errore( 'esm_force_ewr', 'optimal znrm not found', 1 )
       CALL qromb( vl11, aaa, tmp, z1, znrm, znrm, 0.0_DP, gpmax, ss )
       err = -2.d0*tmp/SQRT(pi)-ss*2.d0
       IF (ABS(err) < eps) EXIT
    ENDDO
    ! Define rmax for real space sum
    rmax = 1.d0
    !
    DO
       rmax = rmax+1.d0
       IF (rmax > 200.d0) &
          CALL errore( 'esm_force_ewr', 'optimal rmax not found', 1 )
       CALL qromb( dvl11j0, aaa, tmp, z1, z1-zbuff, z1-zbuff, rmax, gpmax, ss )
       err = ss
       IF (ABS(err) < epsneib) EXIT
    ENDDO
    !
    rmax = rmax/alat
    !
    IF (iverbosity > 0) THEN
       WRITE( stdout, '(5x,"=== Smooth-ESM RSUM parameters (Force) ===")')
       WRITE( stdout, '(5x,A,F10.2,A)') &
          'Upper bound of g_parallel integral:      ',gpmax,' (1/a.u.)'
       WRITE( stdout, '(5x,A,F10.2,A)') &
          'Boundary for normal RSUM|Smooth-ESM RSUM:',z1-znrm,' (a.u.)'
       WRITE( stdout, '(5x,A,F10.2,A)') &
          'Upper bound of real-space summation:     ',rmax*alat,' (a.u.)'
       WRITE( stdout, '(5x,"==========================================")')
    ENDIF
    !
    ip = mp_rank( intra_bgrp_comm )
    np = mp_size( intra_bgrp_comm )
    !
    ALLOCATE( force(3,nat) )
    !
    force(:,:) = 0.d0
    !
    DO na = ip+1, nat, np
       z = tau(3,na)
       IF (z > at(3,3)*0.5) z = z-at(3,3)
       z = z*alat
       !
       DO nb = 1, nat
          IF (nb == na) CYCLE
          zp = tau(3,nb)
          IF (zp > at(3,3)*0.5) zp = zp-at(3,3)
          zp = zp*alat
          dtau(1:2) = tau(1:2,na)-tau(1:2,nb)
          dtau(3) = (z-zp)/alat
          fac = zv(ityp(na))*zv(ityp(nb))*e2
          !
          IF (z < znrm) THEN
             IF (zp < znrm) THEN ! z in I, zp in I (normal RSUM)
                !
                rmax0 = 5.d0 / tmp / alat
                !
                ! generates nearest-neighbors shells r(i)=R(i)-dtau(i)
                CALL rgen( dtau, rmax0, mxr, at, bg, r, r2, nrm )
                !
                ! and sum to the real space part
                DO nr = 1, nrm
                   rr = SQRT(r2(nr))*alat
                   DO ipol = 1, 3
                      force(ipol,na) = force(ipol,na) &
                         -fac/rr**2*(qe_erfc(tmp*rr)/rr+2.d0*tmp/SQRT(pi) &
                         *EXP(-tmp**2*rr**2))*r(ipol,nr)*alat
                   ENDDO
                ENDDO
                !
             ELSEIF (zp < z1) THEN ! z in I, zp in I
                !
                CALL esm_rgen_2d( dtau, rmax, mxr, at, bg, r, r2, nrm )
                !
                DO nr = 1, nrm
                   rxy = SQRT(r2(nr))*alat
                   rxyz = SQRT(r2(nr)+dtau(3)**2)*alat
                   CALL qromb( vl11j1, aaa, tmp, z1, z, zp, rxy, gpmax, ss )
                   force(1:2,nb) = force(1:2,nb) - fac*(1.d0/rxyz**3 + &
                                                   1.d0/rxy*ss)*r(1:2,nr)*alat
                   CALL qromb( dvl11j0, aaa, tmp, z1, z, zp, rxy, gpmax, ss )
                   force(3,nb)=force(3,nb)-fac*((z-zp)/rxyz**3+ss)
                ENDDO
                !
             ELSE ! z in I, zp in II
                !
                CALL esm_rgen_2d( dtau, rmax, mxr, at, bg, r, r2, nrm )
                !
                DO nr = 1, nrm
                   rxy = SQRT(r2(nr))*alat
                   CALL qromb( vl12j1, aaa, tmp, z1, z, zp, rxy, gpmax, ss )
                   force(1:2,nb) = force(1:2,nb)-fac*ss/rxy*r(1:2,nr)*alat
                   CALL qromb( dvl12j0, aaa, tmp, z1, z, zp, rxy, gpmax, ss )
                   force(3,nb) = force(3,nb)-fac*ss
                ENDDO
                !
             ENDIF ! IF for zp
             !
          ELSEIF (z < z1) THEN ! znrm < z < z1
             !
             CALL esm_rgen_2d( dtau, rmax, mxr, at, bg, r, r2, nrm )
             !
             IF (zp < z1) THEN ! z in I, zp in I
                !
                DO nr = 1, nrm
                   rxy = SQRT(r2(nr))*alat
                   rxyz = SQRT(r2(nr)+dtau(3)**2)*alat
                   CALL qromb( vl11j1, aaa, tmp, z1, z, zp, rxy, gpmax, ss )
                   force(1:2,nb) = force(1:2,nb)-fac*(1.d0/rxyz**3 + &
                                      1.d0/rxy*ss)*r(1:2,nr)*alat
                   CALL qromb( dvl11j0, aaa, tmp, z1, z, zp, rxy, gpmax, ss )
                   force(3,nb) = force(3,nb)-fac*((z-zp)/rxyz**3+ss)
                ENDDO
                !
             ELSE ! z in I, zp in II
                !
                DO nr = 1, nrm
                   rxy = SQRT(r2(nr))*alat
                   CALL qromb( vl12j1, aaa, tmp, z1, z, zp, rxy, gpmax, ss )
                   force(1:2,nb) = force(1:2,nb) &
                      -fac*ss/rxy*r(1:2,nr)*alat
                   CALL qromb( dvl12j0, aaa, tmp, z1, z, zp, rxy, gpmax, ss )
                   force(3,nb) = force(3,nb) - fac*ss
                ENDDO
                !
             ENDIF ! IF for zp
             !
          ELSE ! z1 < z
             !
             CALL esm_rgen_2d( dtau, rmax, mxr, at, bg, r, r2, nrm )
             !
             IF (zp < z1 ) THEN ! z in II, zp in I
                !
                DO nr = 1, nrm
                   rxy = SQRT(r2(nr))*alat
                   CALL qromb( vl21j1, aaa, tmp, z1, z, zp, rxy, gpmax, ss )
                   force(1:2,nb) = force(1:2,nb)-fac*ss/rxy*r(1:2,nr)*alat
                   CALL qromb( dvl21j0, aaa, tmp, z1, z, zp, rxy, gpmax, ss )
                   force(3,nb) = force(3,nb) - fac*ss
                ENDDO
                !
             ELSE ! z in II, zp in II
                !
                DO nr = 1, nrm
                   rxy = SQRT(r2(nr))*alat
                   rxyz = SQRT(r2(nr)+dtau(3)**2)*alat
                   CALL qromb( vl22j1, aaa, tmp, z1, z, zp, rxy, gpmax, ss )
                   force(1:2,nb) = force(1:2,nb) &
                      -(EXP(-aaa*(rxyz+z+zp-2.d0*z1))*(aaa+1.d0/rxyz)/rxyz**2 &
                      +ss/rxy)*fac*r(1:2,nr)*alat
                   CALL qromb( dvl22j0, aaa, tmp, z1, z, zp, rxy, gpmax, ss )
                   force(3,nb) = force(3,nb) &
                      -(EXP(-aaa*(rxyz+z+zp-2.d0*z1))*(aaa+1.d0/rxyz)/rxyz**2 &
                      *(z-zp)-aaa*EXP(-aaa*(rxyz+z+zp-2.d0*z1))/rxyz+ss)*fac
                ENDDO
                !
             ENDIF ! IF for zp
             !
          ENDIF
          !
       ENDDO
       !
       IF (z < znrm) THEN
          ss = 0.d0
       ELSEIF (z < z1) THEN
          CALL qromb( dvl11, aaa, tmp, z1, z, z, 0.0_DP, gpmax, ss )
       ELSE
          CALL qromb( dvl22, aaa, tmp, z1, z, z, 0.0_DP, gpmax, ss )
       ENDIF
       ! factor e2: hartree -> Ry.
       force(3,na) = force(3,na)-zv(ityp(na))**2*e2*ss
    ENDDO
    !
    forceion(:,:) = forceion(:,:) + force(:,:)
    !
    DEALLOCATE( force )
    !
  END SUBROUTINE esm_force_ewr_bc4
  !
  !
  ! ...  ESM EWALD-DERIVED FORCE (GSUM) SUBROUTINE
  !
  !-----------------------------------------------------------------------
  SUBROUTINE esm_force_ewg_pbc( alpha_g, forceion )
    !---------------------------------------------------------------------
    !! ESM Ewald-derived force (Gsum) - pbc.
    !
    USE constants,        ONLY : tpi, e2
    USE cell_base,        ONLY : omega, alat, tpiba2, at, bg
    USE control_flags,    ONLY : gamma_only
    USE ions_base,        ONLY : zv, nat, nsp, tau, ityp
    USE gvect,            ONLY : gstart, ngm, gg, g
    USE vlocal,           ONLY : strf
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN) :: alpha_g
    !! alpha term in ewald sum
    REAL(DP), INTENT(OUT) :: forceion(3,nat)
    !! ESM Ewald-derived force
    !
    ! ... local variables
    !
    INTEGER :: nt, ig, na, ipol
    REAL(DP) :: fact, arg, sumnb
    COMPLEX(DP), ALLOCATABLE :: aux(:)
    !
    forceion(:,:) = 0.d0
    !
    ! same of the GSUM part in force_ew.f90
    ALLOCATE( aux(ngm) )
    !
    aux(:) = (0.d0, 0.d0)
    !
    DO nt = 1, nsp
       DO ig = gstart, ngm
          aux(ig) = aux(ig) + zv(nt) * CONJG(strf(ig,nt))
       ENDDO
    ENDDO
    !
    DO ig = gstart, ngm
       aux(ig) = aux(ig) * EXP(-gg(ig) * tpiba2 / alpha_g / 4.d0) &
                 / (gg(ig) * tpiba2)
    ENDDO
    !
    IF (gamma_only) THEN
       fact = 4.d0
    ELSE
       fact = 2.d0
    ENDIF
    !
    DO na = 1, nat
       DO ig = gstart, ngm
          arg = tpi * ( g(1,ig) * tau(1,na) + g(2,ig) * tau(2,na) &
                + g(3,ig) * tau(3,na) )
          sumnb = COS(arg) * AIMAG(aux(ig)) - SIN(arg) * DBLE(aux(ig))
          forceion(1,na) = forceion(1,na) + g(1,ig) * sumnb
          forceion(2,na) = forceion(2,na) + g(2,ig) * sumnb
          forceion(3,na) = forceion(3,na) + g(3,ig) * sumnb 
       ENDDO
       !
       DO ipol = 1, 3
          forceion (ipol, na) = - zv (ityp (na) ) * fact * e2 * tpi**2 / &
             omega / alat * forceion (ipol, na)
       ENDDO
    ENDDO
    !
    DEALLOCATE( aux )
    !
    RETURN
    !
  END SUBROUTINE esm_force_ewg_pbc
  !
  !
  !-----------------------------------------------------------------------------
  SUBROUTINE esm_force_ewg_bc1( alpha_g, forceion )
    !---------------------------------------------------------------------------
    !! ESM Ewald-derived force (Gsum) - bc1.
    !
    USE constants,        ONLY : tpi, fpi, e2
    USE cell_base,        ONLY : omega, alat, tpiba2, at, bg
    USE control_flags,    ONLY : gamma_only
    USE ions_base,        ONLY : zv, nat, nsp, tau, ityp
    USE gvect,            ONLY : gstart, ngm, g
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN) :: alpha_g
    !! alpha term in Ewald sum
    REAL(DP), INTENT(OUT) :: forceion(3,nat)
    !! ESM Ewald-derived force
    !
    ! ... local variables
    !
    INTEGER :: it1, it2, k1, k2, ng_2d
    REAL(DP) :: for(3, nat), for_g(3, nat), t1_for, t2_for,     &
                c1_for(3), c2_for(3), kk1_for, kk2_for, t1, t2, &
                ff, z0, z1, z, zp, tmp, gp2, gp, t(2), L, sa,   &
                arg001, arg002, arg101, arg102
    !
    forceion(:,:) = 0.d0
    for_g(:,:) = 0.d0
    L = at(3,3)*alat
    sa = omega/L
    z0 = L/2.d0
    z1 = z0+esm_w
    tmp = SQRT(alpha_g)
    !
    for = 0.d0
    !
    DO it1 = 1, nat
      DO it2 = 1, nat
         !
         z = tau(3,it1)
         IF (z > at(3,3)*0.5) z = z-at(3,3)
         z = z*alat
         zp = tau(3,it2)
         IF (zp > at(3,3)*0.5) zp = zp-at(3,3)
         zp = zp*alat
         !
         IF (gamma_only) THEN 
            t1_for = zv(ityp(it1))*zv(ityp(it2))*fpi/sa*2.d0
         ELSE
            t1_for = zv(ityp(it1))*zv(ityp(it2))*fpi/sa
         ENDIF
         !
         t2_for = zv(ityp(it1))*zv(ityp(it2))*fpi/sa
         ! bc1
         kk1_for = 0.5d0*qe_erf(tmp*(z-zp))
         kk2_for = 0.d0
         !
         c1_for(:)=0.d0; c2_for(:)=0.d0
         !
         DO ng_2d = 1, ngm_2d
            k1 = mill_2d(1,ng_2d)
            k2 = mill_2d(2,ng_2d)
            IF (k1==0 .AND. k2==0) CYCLE
            t(1:2) = k1*bg(1:2,1)+k2*bg(1:2,2)
            !
            gp2 = SUM(t(:)*t(:))*tpiba2
            gp = SQRT(gp2)
            !
            ff = ((k1*bg(1,1)+k2*bg(1,2))*(tau(1,it1)-tau(1,it2))  &
                 +(k1*bg(2,1)+k2*bg(2,2))*(tau(2,it1)-tau(2,it2)))*tpi
            ! bc1
            arg001 = -gp*(z-zp)
            arg002 =  gp*(z-zp)
            arg101 =  gp/2.d0/tmp-tmp*(z-zp)
            arg102 =  gp/2.d0/tmp+tmp*(z-zp)
            !
            t1 = exp_erfc(arg001,arg101)
            t2 = exp_erfc(arg002,arg102)
            !
            c1_for(1) = c1_for(1) + SIN(ff)*(t1+t2)/4.d0/gp*k1
            c1_for(2) = c1_for(2) + SIN(ff)*(t1+t2)/4.d0/gp*k2
            c1_for(3) = c1_for(3) + COS(ff)*(t1-t2)/4.d0
         ENDDO
         !
         for(:,it2) = for(:,it2)+t1_for*(c1_for(:)+c2_for(:))
         !
         IF (gstart == 2) THEN
            for(3,it2) = for(3,it2)+t2_for*(kk1_for+kk2_for)
         ENDIF
         !
      ENDDO
    ENDDO
    !
    for_g(:,:) = for_g(:,:) + for(:,:)
    !
    for_g(:,:) = for_g(:,:)*e2 ! factor e2: hartree -> Ry.
    !
    DO it1 = 1, nat
       forceion(1,it1) = -SUM( for_g(1:2,it1)*bg(1,1:2) )*SQRT(tpiba2)
       forceion(2,it1) = -SUM( for_g(1:2,it1)*bg(2,1:2) )*SQRT(tpiba2)
       forceion(3,it1) = -for_g(3,it1)
    ENDDO
    !
    RETURN
    !
  END SUBROUTINE esm_force_ewg_bc1
  !
  !
  !------------------------------------------------------------------------
  SUBROUTINE esm_force_ewg_bc2( alpha_g, forceion )
    !----------------------------------------------------------------------
    !! ESM Ewald-derived force (Gsum) - bc2.
    !
    USE constants,        ONLY : tpi, fpi, e2
    USE cell_base,        ONLY : omega, alat, tpiba2, at, bg
    USE control_flags,    ONLY : gamma_only
    USE ions_base,        ONLY : zv, nat, nsp, tau, ityp
    USE gvect,            ONLY : gstart, ngm, g
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN) :: alpha_g
    !! alpha term in Ewald sum
    REAL(DP), INTENT(OUT) :: forceion(3,nat)
    !! ESM Ewald-derived force
    !
    ! ... local variables
    !
    INTEGER :: it1, it2, k1, k2, ng_2d
    REAL(DP) :: for(3, nat), for_g(3, nat), t1_for, t2_for,     &
                c1_for(3), c2_for(3), kk1_for, kk2_for, t1, t2, &
                ff, z0, z1, z, zp, tmp, gp2, gp, t(2), L, sa,   &
                arg001, arg002, arg003, arg004, arg005,         &
                arg006, arg007, arg101, arg102
    !
    forceion(:,:) = 0.d0
    for_g(:,:) = 0.d0
    L = at(3,3)*alat
    sa = omega/L
    z0 = L/2.d0
    z1 = z0+esm_w
    tmp = SQRT(alpha_g)
    !
    for = 0.d0
    !
    DO it1 = 1, nat
      DO it2 = 1, nat
         !
         z = tau(3,it1)
         IF (z > at(3,3)*0.5) z = z-at(3,3)
         z = z*alat
         zp = tau(3,it2)
         IF (zp > at(3,3)*0.5) zp = zp-at(3,3)
         zp = zp*alat
         IF (gamma_only) THEN 
            t1_for = zv(ityp(it1))*zv(ityp(it2))*fpi/sa*2.d0
         ELSE
            t1_for = zv(ityp(it1))*zv(ityp(it2))*fpi/sa
         ENDIF
         t2_for = zv(ityp(it1))*zv(ityp(it2))*fpi/sa
         ! bc2
         kk1_for =  0.5d0*qe_erf(tmp*(z-zp))
         kk2_for = -0.5d0*(z/z1)
         !
         c1_for(:)=0.d0; c2_for(:)=0.d0
         !
         DO ng_2d = 1, ngm_2d
            !
            k1 = mill_2d(1,ng_2d)
            k2 = mill_2d(2,ng_2d)
            IF (k1==0 .AND. k2==0) CYCLE
            t(1:2) = k1*bg(1:2,1)+k2*bg(1:2,2)
            !
            gp2 = SUM(t(:)*t(:))*tpiba2
            gp = SQRT(gp2)
            !
            ff = ((k1*bg(1,1)+k2*bg(1,2))*(tau(1,it1)-tau(1,it2))  &
                 +(k1*bg(2,1)+k2*bg(2,2))*(tau(2,it1)-tau(2,it2)))*tpi
            ! bc2
            arg001 = -gp*(z-zp)
            arg002 =  gp*(z-zp)
            arg003 = -gp*(z+zp+2.d0*z1)
            arg004 =  gp*(z+zp-2.d0*z1)
            arg005 = -gp*(z-zp+4.d0*z1)
            arg006 =  gp*(z-zp-4.d0*z1)
            arg007 = -4.d0*gp*z1
            arg101 =  gp/2.d0/tmp-tmp*(z-zp)
            arg102 =  gp/2.d0/tmp+tmp*(z-zp)
            !
            t1 = exp_erfc(arg001,arg101)
            t2 = exp_erfc(arg002,arg102)
            !
            c1_for(1) = c1_for(1) + SIN(ff)*(t1+t2)/4.d0/gp*k1
            c1_for(2) = c1_for(2) + SIN(ff)*(t1+t2)/4.d0/gp*k2
            c1_for(3) = c1_for(3) + COS(ff)*(t1-t2)/4.d0
            c2_for(1) = c2_for(1) + SIN(ff)*(EXP(arg006)+EXP(arg005) &
                                  - EXP(arg004)-EXP(arg003))/(1.d0   &
                                  - EXP(arg007))/2.d0/gp*k1
            c2_for(2) = c2_for(2) + SIN(ff)*(EXP(arg006)+EXP(arg005) &
                                  - EXP(arg004)-EXP(arg003))/(1.d0   &
                                  - EXP(arg007))/2.d0/gp*k2
            c2_for(3) = c2_for(3) - COS(ff)*(EXP(arg006)-EXP(arg005) &
                                  + EXP(arg004)-EXP(arg003))/(1.d0   &
                                  - EXP(arg007))/2.d0
            !
         ENDDO
         !
         for(:,it2) = for(:,it2) + t1_for*(c1_for(:) + c2_for(:))
         !
         IF (gstart == 2) THEN
            for(3,it2) = for(3,it2) + t2_for*(kk1_for + kk2_for)
         ENDIF
         !
      ENDDO
    ENDDO
    !
    for_g(:,:) = for_g(:,:) + for(:,:)
    for_g(:,:) = for_g(:,:)*e2 ! factor e2: hartree -> Ry.
    !
    DO it1 = 1, nat
       forceion(1,it1) = -SUM( for_g(1:2,it1)*bg(1,1:2) )*SQRT(tpiba2)
       forceion(2,it1) = -SUM( for_g(1:2,it1)*bg(2,1:2) )*SQRT(tpiba2)
       forceion(3,it1) = -for_g(3,it1)
       IF (gstart == 2) THEN
          !! add coulomb fource of ions under efield
          forceion(3,it1) = forceion(3,it1) - zv(ityp(it1))*esm_efield
       ENDIF
    ENDDO
    !
    RETURN
    !
  END SUBROUTINE esm_force_ewg_bc2
  !
  !
  !---------------------------------------------------------------------------
  SUBROUTINE esm_force_ewg_bc3( alpha_g, forceion )
    !--------------------------------------------------------------------------
    !! ESM Ewald-derived force (Gsum) - bc3.
    !
    USE constants,        ONLY : tpi, fpi, e2
    USE cell_base,        ONLY : omega, alat, tpiba2, at, bg
    USE control_flags,    ONLY : gamma_only
    USE ions_base,        ONLY : zv, nat, nsp, tau, ityp
    USE gvect,            ONLY : gstart, ngm, g
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN) :: alpha_g
    !! alpha term in Ewald sum
    REAL(DP), INTENT(OUT) :: forceion(3,nat)
    !! ESM Ewald-derived force 
    !
    ! ... local variables
    !
    INTEGER :: it1, it2, k1, k2, ng_2d
    REAL(DP) :: for(3, nat), for_g(3, nat), t1_for, t2_for,     &
                c1_for(3), c2_for(3), kk1_for, kk2_for, t1, t2, &
                ff, z0, z1, z, zp, tmp, gp2, gp, t(2), L, sa,   &
                arg001, arg002, arg003, arg101, arg102
    !
    forceion(:,:) = 0.d0
    for_g(:,:) = 0.d0
    L = at(3,3)*alat
    sa = omega/L
    z0 = L/2.d0
    z1 = z0+esm_w
    tmp = SQRT(alpha_g)
    !
    for = 0.d0
    !
    DO it1 = 1, nat
      DO it2 = 1, nat
         !
         z = tau(3,it1)
         IF (z > at(3,3)*0.5) z = z-at(3,3)
         z = z*alat
         zp = tau(3,it2)
         IF (zp > at(3,3)*0.5) zp = zp-at(3,3)
         zp = zp*alat
         IF (gamma_only) THEN 
            t1_for = zv(ityp(it1))*zv(ityp(it2))*fpi/sa*2.d0
         ELSE
            t1_for = zv(ityp(it1))*zv(ityp(it2))*fpi/sa
         ENDIF
         t2_for = zv(ityp(it1))*zv(ityp(it2))*fpi/sa
         ! bc3
         kk1_for = 0.5d0*qe_erf(tmp*(z-zp))
         kk2_for = -0.5d0
         !
         c1_for(:)=0.d0; c2_for(:)=0.d0
         !
         DO ng_2d = 1, ngm_2d
            k1 = mill_2d(1,ng_2d)
            k2 = mill_2d(2,ng_2d)
            IF (k1==0 .AND. k2==0) CYCLE
            !
            t(1:2)=k1*bg(1:2,1)+k2*bg(1:2,2)
            gp2 = SUM(t(:)*t(:))*tpiba2
            gp = SQRT(gp2)
            ff = ((k1*bg(1,1)+k2*bg(1,2))*(tau(1,it1)-tau(1,it2))  &
                 +(k1*bg(2,1)+k2*bg(2,2))*(tau(2,it1)-tau(2,it2)))*tpi
            !
            ! bc3
            arg001 = -gp*(z-zp)
            arg002 =  gp*(z-zp)
            arg003 =  gp*(z+zp-2.d0*z1)
            arg101 =  gp/2.d0/tmp-tmp*(z-zp)
            arg102 =  gp/2.d0/tmp+tmp*(z-zp)
            !
            t1 = exp_erfc(arg001,arg101)
            t2 = exp_erfc(arg002,arg102)
            !
            c1_for(1) = c1_for(1) + SIN(ff)*(t1+t2)/4.d0/gp*k1
            c1_for(2) = c1_for(2) + SIN(ff)*(t1+t2)/4.d0/gp*k2
            c1_for(3) = c1_for(3) + COS(ff)*(t1-t2)/4.d0
            c2_for(1) = c2_for(1) + SIN(ff)*(-EXP(arg003))/2.d0/gp*k1
            c2_for(2) = c2_for(2) + SIN(ff)*(-EXP(arg003))/2.d0/gp*k2
            c2_for(3) = c2_for(3) + COS(ff)*(-EXP(arg003))/2.d0
         ENDDO
         !
         for(:,it2) = for(:,it2)+t1_for*(c1_for(:)+c2_for(:))
         !
         IF (gstart == 2) THEN
            for(3,it2) = for(3,it2)+t2_for*(kk1_for+kk2_for)
         ENDIF
         !
      ENDDO
    ENDDO
    !
    for_g(:,:) = for_g(:,:) + for(:,:)
    for_g(:,:) = for_g(:,:)*e2 ! factor e2: hartree -> Ry.
    !
    DO it1 = 1, nat
       forceion(1,it1) = -SUM( for_g(1:2,it1)*bg(1,1:2) )*SQRT(tpiba2)
       forceion(2,it1) = -SUM( for_g(1:2,it1)*bg(2,1:2) )*SQRT(tpiba2)
       forceion(3,it1) = -for_g(3,it1)
    ENDDO
    !
    RETURN
    !
  END SUBROUTINE esm_force_ewg_bc3
  !
  !
  !-------------------------------------------------------------------------
  SUBROUTINE esm_force_ewg_bc4( alpha_g, forceion )
    !-----------------------------------------------------------------------
    !! ESM Ewald-derived force (Gsum) - bc4.
    !
    USE constants,        ONLY : tpi, fpi, e2
    USE cell_base,        ONLY : omega, alat, tpiba2, at, bg
    USE control_flags,    ONLY : gamma_only
    USE ions_base,        ONLY : zv, nat, nsp, tau, ityp
    USE gvect,            ONLY : gstart, ngm, g
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN) :: alpha_g
    !! alpha term in Ewald sum
    REAL(DP), INTENT(OUT) :: forceion(3,nat)
    !! ESM Ewald-derived force 
    !
    ! ... local variables
    !
    INTEGER :: it1, it2, k1, k2, ng_2d
    REAL(DP) :: for(3, nat), for_g(3, nat), t1_for, t2_for,     &
                c1_for(3), c2_for(3), kk1_for, kk2_for, t1, t2, &
                ff, z0, z1, z, zp, tmp, gp2, gp, t(2), L, sa,   &
                arg001, arg002, arg003, arg004, arg005, &
                arg006, arg007, arg008, arg009, arg010, &
                arg011, arg012, arg101, arg102, arg103, &
                arg104, arg105, arg106, arg107, arg108, &
                arg109, arg110, arg111, arg112, arg113, &
                arg114, aaa, t3, alpha, beta, kappa, lambda, &
                xi, chi
    ! auxiliary space
    !
    forceion(:,:) = 0.d0
    for_g(:,:) = 0.d0
    L  = at(3,3)*alat
    sa = omega/L
    z0 = L/2.d0
    aaa = esm_a
    z1  = z0+esm_w
    tmp = SQRT(alpha_g)
    !
    for = 0.d0
    !
    DO it1 = 1, nat
      DO it2 = 1, nat
         z = tau(3,it1)
         IF (z > at(3,3)*0.5) z = z-at(3,3)
         z = z*alat
         zp = tau(3,it2)
         IF (zp > at(3,3)*0.5) zp = zp-at(3,3)
         zp = zp*alat
         !
         IF (gamma_only) THEN 
            t1_for = zv(ityp(it1))*zv(ityp(it2))*fpi/sa*2.d0
         ELSE
            t1_for = zv(ityp(it1))*zv(ityp(it2))*fpi/sa
         ENDIF
         t2_for = zv(ityp(it1))*zv(ityp(it2))*fpi/sa
         !
         ! bc4
         arg004 = -2.d0*aaa*(zp-z1)
         arg006 =  aaa**2/tmp**2+2.d0*aaa*(z1-zp)
         arg101 =  tmp*(z-zp)
         arg102 =  tmp*(z1-zp)
         arg104 =  aaa/tmp+tmp*(z-zp)
         arg106 =  aaa/tmp+tmp*(z1-zp)
         !
         IF (z < z1)THEN  ! factor 1/2 <- non-reciprocality
            IF (zp < z1)THEN
               kk1_for = 0.5d0*(qe_erf(arg101)-qe_erf(arg102))/2.d0 &
                        -0.5d0*exp_erfc(arg006,arg106)/2.d0
               kk2_for =-0.5d0*qe_erfc(arg101)/2.d0
            ELSE
               kk1_for = 0.5d0*(qe_erf(arg101)-qe_erf(arg102))/2.d0 &
                        -0.5d0*exp_erfc(arg006,arg106)/2.d0
               kk2_for =-0.5d0*exp_erfc(arg004,arg101)/2.d0
            ENDIF
         ELSE
            IF ( zp < z1 )THEN
               kk1_for = -0.5d0*exp_erfc(arg006,arg104)/2.d0
               kk2_for = -0.5d0*qe_erfc(arg101)/2.d0
            ELSE
               kk1_for = -0.5d0*exp_erfc(arg006,arg104)/2.d0
               kk2_for = -0.5d0*exp_erfc(arg004,arg101)/2.d0
            ENDIF
         ENDIF
         !
         c1_for(:)=0.d0; c2_for(:)=0.d0
         !
         DO ng_2d = 1, ngm_2d
            k1 = mill_2d(1,ng_2d)
            k2 = mill_2d(2,ng_2d)
            IF (k1==0 .AND. k2==0) CYCLE
            !
            t(1:2) = k1*bg(1:2,1)+k2*bg(1:2,2)
            gp2 = SUM(t(:)*t(:))*tpiba2
            gp = SQRT(gp2)
            !
            ff = ((k1*bg(1,1)+k2*bg(1,2))*(tau(1,it1)-tau(1,it2))  &
                 +(k1*bg(2,1)+k2*bg(2,2))*(tau(2,it1)-tau(2,it2)))*tpi
            !
            ! bc4
            alpha = aaa+gp+SQRT(aaa**2+gp**2)
            beta  = aaa+gp-SQRT(aaa**2+gp**2)
            kappa = aaa-gp+SQRT(aaa**2+gp**2)
            xi    = aaa   +SQRT(aaa**2+gp**2)
            chi   = aaa   -SQRT(aaa**2+gp**2)
            lambda=        SQRT(aaa**2+gp**2)
            !
            arg001 = gp*(z-zp)
            arg002 =-gp*(z-zp)
            arg003 = gp*(z+zp-2.d0*z1)
            arg004 = gp*(z-z1)+xi*(z1-zp)
            arg005 =-gp*(z1-zp)-xi*(z-z1)
            arg006 = aaa/2.d0/tmp**2*xi+gp*(z-z1)+xi*(z1-zp)
            arg007 = aaa/2.d0/tmp**2*xi-gp*(z1-zp)-xi*(z-z1)
            arg008 = aaa/2.d0/tmp**2*xi+xi*(z1-zp)-chi*(z-z1)
            arg009 = aaa/2.d0/tmp**2*xi+xi*(z1-zp)-xi*(z-z1)
            arg010 = aaa/2.d0/tmp**2*xi+chi*(z1-zp)-xi*(z-z1)
            arg011 = aaa/2.d0/tmp**2*chi+chi*(z1-zp)-xi*(z-z1)
            arg012 = aaa/2.d0/tmp**2*chi+xi*(z1-zp)-chi*(z-z1)
            arg101 =  gp/2.d0/tmp+tmp*(z-zp)
            arg102 =  gp/2.d0/tmp-tmp*(z-zp)
            arg103 =  gp/2.d0/tmp+tmp*(z1-zp)
            arg104 =  gp/2.d0/tmp-tmp*(z1-zp)
            arg105 =  gp/2.d0/tmp+tmp*(z-z1)
            arg106 =  gp/2.d0/tmp-tmp*(z-z1)
            arg107 =  xi/2.d0/tmp+tmp*(z-zp)
            arg108 =  xi/2.d0/tmp-tmp*(z-zp)
            arg109 =  xi/2.d0/tmp+tmp*(z1-zp)
            arg110 =  xi/2.d0/tmp-tmp*(z-z1)
            arg111 = chi/2.d0/tmp+tmp*(z-zp)
            arg112 = chi/2.d0/tmp-tmp*(z-zp)
            arg113 = chi/2.d0/tmp+tmp*(z1-zp)
            arg114 = chi/2.d0/tmp-tmp*(z-z1)
            !
            IF (z < z1) THEN ! factor 1/2 <- non-reciprocality
               !
               IF (zp < z1) THEN
                  !
                  t1 = exp_erfc(arg001,arg101)-exp_erfc(arg001,arg103)
                  t2 = exp_erfc(arg002,arg102)     &
                      -kappa/alpha*exp_erfc(arg003,arg104)
                  t3 = exp_erfc(arg006,arg109)/alpha
                  !
                  c1_for(1) = c1_for(1)+SIN(ff)*((t1+t2)/4.d0/gp+t3/2.d0)*k1/2.d0
                  c1_for(2) = c1_for(2)+SIN(ff)*((t1+t2)/4.d0/gp+t3/2.d0)*k2/2.d0
                  !
                  t1 = exp_erfc(arg002,arg102)-exp_erfc(arg002,arg106)
                  t2 = exp_erfc(arg001,arg101)     &
                      -kappa/alpha*exp_erfc(arg003,arg105)
                  t3 = exp_erfc(arg007,arg110)/alpha
                  !
                  c2_for(1) = c2_for(1)+SIN(ff)*((t1+t2)/4.d0/gp+t3/2.d0)*k1/2.d0
                  c2_for(2) = c2_for(2)+SIN(ff)*((t1+t2)/4.d0/gp+t3/2.d0)*k2/2.d0
                  !
                  t1 = exp_erfc(arg001,arg103)-exp_erfc(arg001,arg101)
                  t2 = exp_erfc(arg002,arg102)     &
                      -kappa/alpha*exp_erfc(arg003,arg104)
                  t3 =-xi/alpha*exp_erfc(arg006,arg109)
                  !
                  c1_for(3) = c1_for(3)+COS(ff)*((t1+t2)/4.d0+t3/2.d0)/2.d0 
                  !
                  t1 = exp_erfc(arg002,arg102)-exp_erfc(arg002,arg106)
                  t2 =-exp_erfc(arg001,arg101)     &
                      -kappa/alpha*exp_erfc(arg003,arg105)
                  t3 = gp/alpha*exp_erfc(arg007,arg110)
                  !
                  c2_for(3) = c2_for(3)+COS(ff)*((t1+t2)/4.d0+t3/2.d0)/2.d0
                  !
               ELSE
                  !
                  t1 = exp_erfc(arg001,arg101)-exp_erfc(arg001,arg103)
                  t2 = exp_erfc(arg002,arg102)     &
                      -kappa/alpha*exp_erfc(arg003,arg104)
                  t3 = exp_erfc(arg006,arg109)/alpha
                  !
                  c1_for(1) = c1_for(1)+SIN(ff)*((t1+t2)/4.d0/gp+t3/2.d0)*k1/2.d0
                  c1_for(2) = c1_for(2)+SIN(ff)*((t1+t2)/4.d0/gp+t3/2.d0)*k2/2.d0
                  !
                  t1 = exp_erfc(arg012,arg114)-exp_erfc(arg012,arg112)
                  t2 = exp_erfc(arg010,arg108)     &
                      -beta/alpha*exp_erfc(arg009,arg110)
                  t3 = exp_erfc(arg004,arg105)/alpha
                  !
                  c2_for(1) = c2_for(1)+SIN(ff)*((t1+t2)/4.d0/gp+t3/2.d0)*k1/2.d0
                  c2_for(2) = c2_for(2)+SIN(ff)*((t1+t2)/4.d0/gp+t3/2.d0)*k2/2.d0
                  !
                  t1 = exp_erfc(arg001,arg103)-exp_erfc(arg001,arg101)
                  t2 = exp_erfc(arg002,arg102)     &
                      -kappa/alpha*exp_erfc(arg003,arg104)
                  t3 =-xi/alpha*exp_erfc(arg006,arg109)
                  !
                  c1_for(3) = c1_for(3)+COS(ff)*((t1+t2)/4.d0+t3/2.d0)/2.d0 
                  !
                  t1 = xi*(exp_erfc(arg012,arg112)-exp_erfc(arg012,arg114))
                  t2 =-chi*exp_erfc(arg010,arg108) &
                      +xi*beta/alpha*exp_erfc(arg009,arg110)
                  t3 =-xi/alpha*exp_erfc(arg004,arg105)
                  !
                  c2_for(3) = c2_for(3)+COS(ff)*((t1+t2)/4.d0+t3/2.d0)/2.d0
                  !
               ENDIF
               !
            ELSE
               !
               IF (zp < z1) THEN
                  !
                  t1 = exp_erfc(arg011,arg113)-exp_erfc(arg011,arg111)
                  t2 = exp_erfc(arg008,arg107)     &
                      -beta/alpha*exp_erfc(arg009,arg109)
                  t3 = exp_erfc(arg005,arg104)/alpha
                  !
                  c1_for(1) = c1_for(1)+SIN(ff)*((t1+t2)/4.d0/lambda+t3/2.d0)*k1/2.d0
                  c1_for(2) = c1_for(2)+SIN(ff)*((t1+t2)/4.d0/lambda+t3/2.d0)*k2/2.d0
                  !
                  t1 =  exp_erfc(arg002,arg102)-exp_erfc(arg002,arg106)
                  t2 =  exp_erfc(arg001,arg101)    &
                       -kappa/alpha*exp_erfc(arg003,arg105)
                  t3 = exp_erfc(arg007,arg110)/alpha
                  !
                  c2_for(1) = c2_for(1)+SIN(ff)*((t1+t2)/4.d0/lambda+t3/2.d0)*k1/2.d0
                  c2_for(2) = c2_for(2)+SIN(ff)*((t1+t2)/4.d0/lambda+t3/2.d0)*k2/2.d0
                  !
                  t1 = chi*(exp_erfc(arg011,arg111)-exp_erfc(arg011,arg113))
                  t2 =-xi*exp_erfc(arg008,arg107)  &
                      +xi*beta/alpha*exp_erfc(arg009,arg109)
                  t3 = gp/alpha*exp_erfc(arg005,arg104)
                  !
                  c1_for(3) = c1_for(3)+COS(ff)*((t1+t2)/4.d0/lambda+t3/2.d0)/2.d0
                  !
                  t1 = exp_erfc(arg002,arg102)-exp_erfc(arg002,arg106)
                  t2 =-exp_erfc(arg001,arg101)     &
                      -kappa/alpha*exp_erfc(arg003,arg105)
                  t3 = gp/alpha*exp_erfc(arg007,arg110)
                  !
                  c2_for(3) = c2_for(3)+COS(ff)*((t1+t2)/4.d0/lambda+t3/2.d0)/2.d0
                  !
               ELSE
                  !
                  t1 = exp_erfc(arg011,arg113)-exp_erfc(arg011,arg111)
                  t2 = exp_erfc(arg008,arg107)     &
                      -beta/alpha*exp_erfc(arg009,arg109)
                  t3 = exp_erfc(arg005,arg104)/alpha
                  !
                  c1_for(1) = c1_for(1)+SIN(ff)*((t1+t2)/4.d0/lambda+t3/2.d0)*k1/2.d0
                  c1_for(2) = c1_for(2)+SIN(ff)*((t1+t2)/4.d0/lambda+t3/2.d0)*k2/2.d0
                  !
                  t1 = exp_erfc(arg012,arg114)-exp_erfc(arg012,arg112)
                  t2 = exp_erfc(arg010,arg108)     &
                      -beta/alpha*exp_erfc(arg009,arg110)
                  t3 = exp_erfc(arg004,arg105)/alpha
                  !
                  c2_for(1) = c2_for(1)+SIN(ff)*((t1+t2)/4.d0/lambda+t3/2.d0)*k1/2.d0
                  c2_for(2) = c2_for(2)+SIN(ff)*((t1+t2)/4.d0/lambda+t3/2.d0)*k2/2.d0
                  !
                  t1 = chi*(exp_erfc(arg011,arg111)-exp_erfc(arg011,arg113))
                  t2 =-xi*exp_erfc(arg008,arg107)  &
                      +xi*beta/alpha*exp_erfc(arg009,arg109)
                  t3 = gp/alpha*exp_erfc(arg005,arg104)
                  !
                  c1_for(3) = c1_for(3)+COS(ff)*((t1+t2)/4.d0/lambda+t3/2.d0)/2.d0
                  !
                  t1 = xi*(exp_erfc(arg012,arg112)-exp_erfc(arg012,arg114))
                  t2 =-chi*exp_erfc(arg010,arg108) &
                      +xi*beta/alpha*exp_erfc(arg009,arg110)
                  t3 =-xi/alpha*exp_erfc(arg004,arg105)
                  !
                  c2_for(3) = c2_for(3)+COS(ff)*((t1+t2)/4.d0/lambda+t3/2.d0)/2.d0
                  !
               ENDIF
               !
            ENDIF
            !
         ENDDO
         !
         for(:,it2)=for(:,it2)+t1_for*(c1_for(:)+c2_for(:))
         !
         IF (gstart == 2) THEN
            for(3,it2) = for(3,it2)+t2_for*(kk1_for+kk2_for)
         ENDIF
         !
      ENDDO
      !
    ENDDO
    !
    for_g(:,:) = for_g(:,:) + for(:,:)
    !
    for_g(:,:) = for_g(:,:)*e2 ! factor e2: hartree -> Ry.
    !
    DO it1 = 1, nat
       forceion(1,it1) = -SUM( for_g(1:2,it1)*bg(1,1:2) )*SQRT(tpiba2)
       forceion(2,it1) = -SUM( for_g(1:2,it1)*bg(2,1:2) )*SQRT(tpiba2)
       forceion(3,it1) = -for_g(3,it1)
    ENDDO
    !
    RETURN
    !
  END SUBROUTINE esm_force_ewg_bc4
  !
  !
  ! ... ESM LOCAL POTENTIAL-DERIVED FORCE SUBROUTINE
  !
  !-----------------------------------------------------------------------
  SUBROUTINE esm_force_lc_pbc( aux, forcelc )
    !---------------------------------------------------------------------
    !! ERROR routine: esm_force_lc must not be called for esm_bc = pbc.
    !
    USE ions_base, ONLY : nat
    USE fft_base,  ONLY : dfftp
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(IN) :: aux(dfftp%nnr)
    !! aux contains n(G) (input)
    REAL(DP), INTENT(INOUT) :: forcelc(3,nat)
    !! lc force
    !
    STOP 'esm_force_lc must not be called for esm_bc = pbc'
    !
  END SUBROUTINE esm_force_lc_pbc
  !
  !
  !------------------------------------------------------------------------
  SUBROUTINE esm_force_lc_bc1( aux, forcelc )
    !----------------------------------------------------------------------
    !! ESM local-potential derived force - bc1.
    !
    USE constants,        ONLY : tpi, fpi, e2
    USE gvect,            ONLY : ngm, mill
    USE cell_base,        ONLY : omega, alat, tpiba2, at, bg
    USE control_flags,    ONLY : gamma_only
    USE ions_base,        ONLY : zv, nat, tau, ityp
    USE fft_base,         ONLY : dfftp
    USE fft_scalar,       ONLY : cft_1z
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(IN) :: aux(dfftp%nnr)
    !! aux contains n(G) (input)   
    REAL(DP), INTENT(INOUT) :: forcelc(3,nat)
    !! lc force
    !
    ! ... local variables
    !
    INTEGER :: iz, it, k1, k2, k3, ng, n1, n2, n3, ng_2d
    REAL(DP), ALLOCATABLE :: for(:,:), for_g(:,:)
    REAL(DP) :: t(2), tt, gp, gp2, sa, z1, z0, pp, cc, ss, t1, &
                t2, z, zp, L, tmp, r1, r2, f1(3), f2(3),       &
                arg001, arg002, arg101, arg102
    COMPLEX(DP), ALLOCATABLE :: vg_f(:,:), vg_f_r(:,:), rhog3(:,:)
    COMPLEX(DP) :: c1(3), c2(3), cc1, cc2
    !
    ! Map to FFT mesh 
    ALLOCATE( rhog3(dfftp%nr3,ngm_2d) )
    !
    rhog3(:,:) = (0.d0,0.d0)
    !
    DO ng = 1, ngm
        n1 = mill(1,ng)
        n2 = mill(2,ng)
        ng_2d = imill_2d(n1,n2)
        n3 = mill(3,ng) + 1
        IF (n3 < 1) n3 = n3 + dfftp%nr3
        rhog3(n3,ng_2d) = aux(dfftp%nl(ng))
        IF (gamma_only .AND. n1==0 .AND. n2==0) THEN
           n3 = -mill(3,ng)+1
           IF (n3 < 1)n3 = n3+dfftp%nr3
           rhog3(n3,ng_2d) = aux(dfftp%nlm(ng))
        ENDIF  
    ENDDO
    !
    L = at(3,3)*alat
    sa = omega/L
    z0 = L/2.d0
    tmp = 1.d0
    z1 = z0+esm_w
    !
    ALLOCATE( for_g(3,nat) )
    for_g(:,:)=0.d0
    !
    !**** for gp!=0 ****
    ALLOCATE( for(3,nat), vg_f(dfftp%nr3x,3), vg_f_r(dfftp%nr3x,3) )
    for(:,:) = 0.d0
    vg_f_r(:,:) = (0.d0,0.d0)
    !
    DO ng_2d = 1, ngm_2d
       k1 = mill_2d(1,ng_2d)
       k2 = mill_2d(2,ng_2d)
       IF (k1==0 .AND. k2==0) CYCLE
       !
       t(1:2) = k1*bg(1:2,1)+k2*bg(1:2,2)
       gp2 = SUM(t(:)*t(:))*tpiba2
       gp  = SQRT(gp2)
       !
       DO it = 1, nat
          IF (gamma_only) THEN
             tt = -fpi*zv(ityp(it))/sa*2.d0
          ELSE 
             tt = -fpi*zv(ityp(it))/sa
          ENDIF
          pp = -tpi*(tau(1,it)*(k1*bg(1,1)+k2*bg(1,2)) &
                    +tau(2,it)*(k1*bg(2,1)+k2*bg(2,2)))
          cc = COS(pp)
          ss = SIN(pp)
          zp = tau(3,it)
          IF (zp > at(3,3)*0.5) zp = zp-at(3,3)
          zp = zp*alat
          !
          DO iz = 1, dfftp%nr3
             k3 = iz-1
             IF (k3 > dfftp%nr3/2) k3 = iz-dfftp%nr3-1
             z = DBLE(k3)/DBLE(dfftp%nr3)*L
             ! bc1
             arg001 = gp*(z-zp)
             arg002 =-gp*(z-zp)
             arg101 = gp/2.d0/tmp+tmp*(z-zp)
             arg102 = gp/2.d0/tmp-tmp*(z-zp)
             !
             t1 = exp_erfc(arg002,arg102)
             t2 = exp_erfc(arg001,arg101)
             !
             c1(1) = CMPLX(ss,-cc,KIND=DP)*(t1+t2)/4.d0/gp*k1
             c1(2) = CMPLX(ss,-cc,KIND=DP)*(t1+t2)/4.d0/gp*k2
             c1(3) = CMPLX(cc, ss,KIND=DP)*(t1-t2)/4.d0
             c2(:) = (0.d0,0.d0)
             !
             vg_f_r(iz,:) = tt*(c1(:)+c2(:))
             !
          ENDDO
          !
          CALL cft_1z( vg_f_r(:,1), 1, dfftp%nr3, dfftp%nr3, -1, vg_f(:,1) )
          CALL cft_1z( vg_f_r(:,2), 1, dfftp%nr3, dfftp%nr3, -1, vg_f(:,2) )
          CALL cft_1z( vg_f_r(:,3), 1, dfftp%nr3, dfftp%nr3, -1, vg_f(:,3) )
          !
          DO iz = 1, dfftp%nr3
             r1 = DBLE(rhog3(iz,ng_2d))
             r2 = AIMAG(rhog3(iz,ng_2d))
             !
             f1(:) =  DBLE(vg_f(iz,:))
             f2(:) = AIMAG(vg_f(iz,:))
             !
             for(:,it) = for(:,it)-r1*f1(:)-r2*f2(:)
          ENDDO
          !
       ENDDO
       !
    ENDDO
    !
    for_g(:,:) = for_g(:,:) + for(:,:)
    !
    DEALLOCATE( for, vg_f, vg_f_r )
    !
    !**** for gp==0****
    ng_2d = imill_2d(0,0)
    !
    IF ( ng_2d > 0 ) THEN
       !
       ALLOCATE( vg_f(dfftp%nr3x,1), vg_f_r(dfftp%nr3x,1) )
       vg_f_r(:,1) = (0.d0,0.d0)
       !
       DO it = 1, nat
          tt = -fpi*zv(ityp(it))/sa
          zp = tau(3,it)
          IF (zp>at(3,3)*0.5) zp = zp-at(3,3)
          zp = zp*alat
          !
          DO iz = 1, dfftp%nr3
             k3 = iz-1
             IF (k3 > dfftp%nr3/2) k3 = iz-dfftp%nr3-1
             z = DBLE(k3)/DBLE(dfftp%nr3)*L
             ! bc1
             cc1 = 0.5d0*qe_erf(tmp*(z-zp))
             cc2 = (0.d0,0.d0)
             !
             vg_f_r(iz,1) = tt*(cc1+cc2)
          ENDDO
          !
          CALL cft_1z( vg_f_r(:,1), 1, dfftp%nr3, dfftp%nr3, -1, vg_f(:,1) )
          !
          DO iz = 1, dfftp%nr3
             r1 =  DBLE(rhog3(iz,ng_2d))
             r2 = AIMAG(rhog3(iz,ng_2d))
             !
             f1(3) =  DBLE(vg_f(iz,1))
             f2(3) = AIMAG(vg_f(iz,1))
             !
             for_g(3,it) = for_g(3,it)-r1*f1(3)-r2*f2(3)
          ENDDO
          !
       ENDDO
       !
       DEALLOCATE( vg_f, vg_f_r )
       !
    ENDIF ! if( ng_2d > 0 )
    !
    !**** sum short_range part and long_range part in local potential force 
    !**** at cartecian coordinate
    !
    DO it = 1, nat
       ! factor e2: hartree -> Ry.
       forcelc(1,it) = forcelc(1,it) &
                      +SUM(for_g(1:2,it)*bg(1,1:2))*SQRT(tpiba2)*omega*e2
       forcelc(2,it) = forcelc(2,it) &
                      +SUM(for_g(1:2,it)*bg(2,1:2))*SQRT(tpiba2)*omega*e2
       forcelc(3,it) = forcelc(3,it)+for_g(3,it)*omega*e2
    ENDDO
    !
    DEALLOCATE( for_g )
    DEALLOCATE( rhog3 )
    !
    RETURN
    !
  END SUBROUTINE esm_force_lc_bc1
  !
  !
  !-------------------------------------------------------------------------
  SUBROUTINE esm_force_lc_bc2( aux, forcelc )
    !------------------------------------------------------------------------
    !! ESM local-potential derived force - bc2.
    !
    USE constants,        ONLY : tpi, fpi, e2
    USE gvect,            ONLY : ngm, mill
    USE cell_base,        ONLY : omega, alat, tpiba2, at, bg
    USE control_flags,    ONLY : gamma_only
    USE ions_base,        ONLY : zv, nat, tau, ityp
    USE fft_base,         ONLY : dfftp
    USE fft_scalar,       ONLY : cft_1z
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(IN) :: aux(dfftp%nnr)
    !! aux contains n(G) (input)   
    REAL(DP), INTENT(INOUT) :: forcelc(3,nat)
    !! lc force
    !
    ! ... local variables
    !
    INTEGER :: iz, it, k1, k2, k3, ng, n1, n2, n3, ng_2d
    REAL(DP), ALLOCATABLE :: for(:,:), for_g(:,:)
    REAL(DP) :: t(2), tt, gp, gp2, sa, z1, z0, pp, cc, ss, t1, &
                t2, z, zp, L, tmp, r1, r2, f1(3), f2(3),       &
                arg001, arg002, arg003, arg005,                &
                arg006, arg008, arg009, arg101, arg102
    COMPLEX(DP), ALLOCATABLE :: vg_f(:,:), vg_f_r(:,:), rhog3(:,:)
    COMPLEX(DP) :: c1(3), c2(3), cc1, cc2
    !
    ! Map to FFT mesh
    ALLOCATE( rhog3(dfftp%nr3,ngm_2d) )
    !
    rhog3(:,:) = (0.d0,0.d0)
    !
    DO ng = 1, ngm
        n1 = mill(1,ng)
        n2 = mill(2,ng)
        ng_2d = imill_2d(n1,n2)
        n3 = mill(3,ng) + 1
        IF (n3 < 1) n3 = n3 + dfftp%nr3
        rhog3(n3,ng_2d)=aux(dfftp%nl(ng))
        IF (gamma_only .AND. n1==0 .AND. n2==0) THEN
           n3 = -mill(3,ng)+1
           IF(n3 < 1)n3 = n3+dfftp%nr3
           rhog3(n3,ng_2d) = aux(dfftp%nlm(ng))
        ENDIF
    ENDDO
    !
    L = at(3,3)*alat
    sa = omega/L
    z0 = L/2.d0
    tmp = 1.d0
    z1 = z0+esm_w
    !
    ALLOCATE( for_g(3,nat) )
    for_g(:,:) = 0.d0
    !
    !**** for gp!=0 ****
    ALLOCATE( for(3,nat), vg_f(dfftp%nr3x,3), vg_f_r(dfftp%nr3x,3) )
    for(:,:) = 0.d0
    vg_f_r(:,:) = (0.d0,0.d0)
    !
    DO ng_2d = 1, ngm_2d
       !
       k1 = mill_2d(1,ng_2d)
       k2 = mill_2d(2,ng_2d)
       IF (k1==0 .AND. k2==0) CYCLE
       !
       t(1:2) = k1*bg(1:2,1) + k2*bg(1:2,2)
       gp2 = SUM(t(:)*t(:))*tpiba2
       gp = SQRT(gp2)
       !
       DO it = 1, nat
          IF (gamma_only) THEN
             tt = -fpi*zv(ityp(it))/sa*2.d0
          ELSE
             tt = -fpi*zv(ityp(it))/sa
          ENDIF
          pp = -tpi*(tau(1,it)*(k1*bg(1,1)+k2*bg(1,2)) &
                    +tau(2,it)*(k1*bg(2,1)+k2*bg(2,2)))
          cc = COS(pp)
          ss = SIN(pp)
          !
          zp = tau(3,it)
          IF (zp > at(3,3)*0.5) zp = zp-at(3,3)
          zp = zp*alat
          !
          DO iz = 1, dfftp%nr3
             k3 = iz-1
             IF (k3 > dfftp%nr3/2) k3 = iz-dfftp%nr3-1
             z = DBLE(k3)/DBLE(dfftp%nr3)*L
             ! bc2
             arg001 = gp*(z-zp)
             arg002 =-gp*(z-zp)
             arg003 =-gp*(z+zp+2.d0*z1)
             arg005 = gp*(z+zp-2.d0*z1)
             arg006 =-gp*(z-zp+4.d0*z1)
             arg008 = gp*(z-zp-4.d0*z1)
             arg009 =-4.d0*gp*z1
             arg101 = gp/2.d0/tmp+tmp*(z-zp)
             arg102 = gp/2.d0/tmp-tmp*(z-zp)
             !
             t1 = exp_erfc(arg002,arg102)
             t2 = exp_erfc(arg001,arg101)
             !
             c1(1) = CMPLX(ss,-cc,KIND=DP)*(t1+t2)/4.d0/gp*k1
             c1(2) = CMPLX(ss,-cc,KIND=DP)*(t1+t2)/4.d0/gp*k2
             c1(3) = CMPLX(cc, ss,KIND=DP)*(t1-t2)/4.d0
             c2(1) = CMPLX(ss,-cc,KIND=DP)*(EXP(arg008)+EXP(arg006)  &
                     -EXP(arg005)-EXP(arg003))/(1.d0-EXP(arg009))/2.d0/gp*k1
             c2(2) = CMPLX(ss,-cc,KIND=DP)*(EXP(arg008)+EXP(arg006)  &
                     -EXP(arg005)-EXP(arg003))/(1.d0-EXP(arg009))/2.d0/gp*k2
             c2(3) = CMPLX(cc, ss,KIND=DP)*(-EXP(arg008)+EXP(arg006) &
                     -EXP(arg005)+EXP(arg003))/(1.d0-EXP(arg009))/2.d0
             !
             vg_f_r(iz,:) = tt*(c1(:)+c2(:))
          ENDDO
          !
          CALL cft_1z( vg_f_r(:,1), 1, dfftp%nr3, dfftp%nr3, -1, vg_f(:,1) )
          CALL cft_1z( vg_f_r(:,2), 1, dfftp%nr3, dfftp%nr3, -1, vg_f(:,2) )
          CALL cft_1z( vg_f_r(:,3), 1, dfftp%nr3, dfftp%nr3, -1, vg_f(:,3) )
          !
          DO iz = 1, dfftp%nr3
             r1 =  DBLE(rhog3(iz,ng_2d))
             r2 = AIMAG(rhog3(iz,ng_2d))
             f1(:) =  DBLE(vg_f(iz,:))
             f2(:) = AIMAG(vg_f(iz,:))
             for(:,it) = for(:,it)-r1*f1(:)-r2*f2(:)
          ENDDO
          !
       ENDDO
       !
    ENDDO
    !
    for_g(:,:) = for_g(:,:) + for(:,:)
    !
    DEALLOCATE( for, vg_f, vg_f_r )
    !
    !**** for gp==0****
    ng_2d = imill_2d(0,0)
    !
    IF ( ng_2d > 0 ) THEN
       ALLOCATE( vg_f(dfftp%nr3x,1), vg_f_r(dfftp%nr3x,1) )
       vg_f_r(:,1) = (0.d0,0.d0)
       !
       DO it = 1, nat
          tt = -fpi*zv(ityp(it))/sa
          zp = tau(3,it)
          IF (zp > at(3,3)*0.5) zp = zp-at(3,3)
          zp = zp*alat
          !
          DO iz = 1, dfftp%nr3
             k3 = iz-1
             IF (k3 > dfftp%nr3/2) k3 = iz-dfftp%nr3-1
             z = DBLE(k3)/DBLE(dfftp%nr3)*L
             ! bc2
             cc1 = 0.5d0*qe_erf(tmp*(z-zp))
             cc2 = -0.5d0*(z/z1)
             vg_f_r(iz,1) = tt*(cc1+cc2)
          ENDDO
          !
          CALL cft_1z( vg_f_r(:,1), 1, dfftp%nr3, dfftp%nr3, -1, vg_f(:,1) )
          !
          DO iz = 1, dfftp%nr3
             r1 =  DBLE(rhog3(iz,ng_2d))
             r2 = AIMAG(rhog3(iz,ng_2d))
             f1(3) =  DBLE(vg_f(iz,1))
             f2(3) = AIMAG(vg_f(iz,1))
             for_g(3,it) = for_g(3,it)-r1*f1(3)-r2*f2(3)
          ENDDO
       ENDDO
       !
       DEALLOCATE( vg_f, vg_f_r )
       !
    ENDIF ! if( ng_2d > 0 )
    !
    !
    !**** sum short_range part and long_range part in local potential force 
    !**** at cartecian coordinate
    !
    DO it = 1, nat
       ! factor e2: hartree -> Ry.
       forcelc(1,it) = forcelc(1,it) &
                       +SUM(for_g(1:2,it)*bg(1,1:2))*SQRT(tpiba2)*omega*e2
       forcelc(2,it) = forcelc(2,it) &
                       +SUM(for_g(1:2,it)*bg(2,1:2))*SQRT(tpiba2)*omega*e2
       forcelc(3,it) = forcelc(3,it)+for_g(3,it)*omega*e2
    ENDDO
    !
    DEALLOCATE( for_g )
    DEALLOCATE( rhog3 )
    !
    RETURN
    !
  END SUBROUTINE esm_force_lc_bc2
  !
  !
  !---------------------------------------------------------------------------------
  SUBROUTINE esm_force_lc_bc3( aux, forcelc )
    !-------------------------------------------------------------------------------
    !! ESM local-potential derived force - bc3.
    !
    USE constants,        ONLY : tpi, fpi, e2
    USE gvect,            ONLY : ngm, mill
    USE cell_base,        ONLY : omega, alat, tpiba2, at, bg
    USE control_flags,    ONLY : gamma_only
    USE ions_base,        ONLY : zv, nat, tau, ityp
    USE fft_base,         ONLY : dfftp
    USE fft_scalar,       ONLY : cft_1z
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(IN) :: aux(dfftp%nnr)
    !! aux contains n(G) (input)   
    REAL(DP), INTENT(INOUT) :: forcelc(3,nat)
    !! lc force
    !
    ! ... local variables
    !
    INTEGER :: iz, it, k1, k2, k3, ng, n1, n2, n3, ng_2d
    REAL(DP), ALLOCATABLE :: for(:,:), for_g(:,:)
    REAL(DP) :: t(2), tt, gp, gp2, sa, z1, z0, pp, cc, ss, t1, &
                t2, z, zp, L, tmp, r1, r2, f1(3), f2(3),       &
                arg001, arg002, arg003, arg101, arg102
    COMPLEX(DP), ALLOCATABLE :: vg_f(:,:), vg_f_r(:,:), rhog3(:,:)
    COMPLEX(DP) :: c1(3), c2(3), cc1, cc2
    !
    ! Map to FFT mesh
    ALLOCATE( rhog3(dfftp%nr3,ngm_2d) )
    rhog3(:,:)=(0.d0,0.d0)
    !
    DO ng=1,ngm
        n1 = mill(1,ng)
        n2 = mill(2,ng)
        ng_2d = imill_2d(n1,n2)
        n3 = mill(3,ng) + 1
        IF (n3<1) n3 = n3 + dfftp%nr3
        rhog3(n3,ng_2d)=aux(dfftp%nl(ng))
        IF (gamma_only .AND. n1==0 .AND. n2==0) THEN
           n3 = -mill(3,ng)+1
           IF(n3<1)n3=n3+dfftp%nr3
           rhog3(n3,ng_2d)=aux(dfftp%nlm(ng))
        ENDIF  
    ENDDO
    !
    L = at(3,3)*alat
    sa = omega/L
    z0 = L/2.d0
    tmp = 1.d0
    z1 = z0 + esm_w
    !
    ALLOCATE( for_g(3,nat) )
    !
    for_g(:,:) = 0.d0
    !
    !**** for gp!=0 ****
    ALLOCATE( for(3,nat), vg_f(dfftp%nr3x,3), vg_f_r(dfftp%nr3x,3) )
    for(:,:) = 0.d0
    vg_f_r(:,:) = (0.d0,0.d0)
    !
    DO ng_2d = 1, ngm_2d
       k1 = mill_2d(1,ng_2d)
       k2 = mill_2d(2,ng_2d)
       IF (k1==0 .AND. k2==0) CYCLE
       !
       t(1:2) = k1*bg(1:2,1)+k2*bg(1:2,2)
       gp2 = SUM(t(:)*t(:))*tpiba2
       gp = SQRT(gp2)
       !
       DO it = 1, nat
          IF (gamma_only) THEN
             tt = -fpi*zv(ityp(it))/sa*2.d0
          ELSE
             tt = -fpi*zv(ityp(it))/sa
          ENDIF 
          pp = -tpi*(tau(1,it)*(k1*bg(1,1)+k2*bg(1,2)) &
                  +tau(2,it)*(k1*bg(2,1)+k2*bg(2,2)))
          cc = COS(pp)
          ss = SIN(pp)
          zp = tau(3,it)
          IF (zp > at(3,3)*0.5) zp = zp-at(3,3)
          zp = zp*alat
          !
          DO iz = 1, dfftp%nr3
             k3 = iz-1
             IF (k3 > dfftp%nr3/2) k3 = iz-dfftp%nr3-1
             z = DBLE(k3)/DBLE(dfftp%nr3)*L
             ! bc3
             arg001 = gp*(z-zp)
             arg002 =-gp*(z-zp)
             arg003 = gp*(z+zp-2.d0*z1)
             arg101 = gp/2.d0/tmp+tmp*(z-zp)
             arg102 = gp/2.d0/tmp-tmp*(z-zp)
             !
             t1 = exp_erfc(arg002,arg102)
             t2 = exp_erfc(arg001,arg101)
             !
             c1(1) = CMPLX(ss,-cc,KIND=DP)*(t1+t2)/4.d0/gp*k1
             c1(2) = CMPLX(ss,-cc,KIND=DP)*(t1+t2)/4.d0/gp*k2
             c1(3) = CMPLX(cc, ss,KIND=DP)*(t1-t2)/4.d0
             c2(1) = CMPLX(ss,-cc,KIND=DP)*(-EXP(arg003))/2.d0/gp*k1
             c2(2) = CMPLX(ss,-cc,KIND=DP)*(-EXP(arg003))/2.d0/gp*k2
             c2(3) = CMPLX(cc, ss,KIND=DP)*(-EXP(arg003))/2.d0
             !
             vg_f_r(iz,:) = tt*(c1(:)+c2(:))
          ENDDO
          !
          CALL cft_1z( vg_f_r(:,1), 1, dfftp%nr3, dfftp%nr3, -1, vg_f(:,1) )
          CALL cft_1z( vg_f_r(:,2), 1, dfftp%nr3, dfftp%nr3, -1, vg_f(:,2) )
          CALL cft_1z( vg_f_r(:,3), 1, dfftp%nr3, dfftp%nr3, -1, vg_f(:,3) )
          !
          DO iz = 1, dfftp%nr3
             r1 =  DBLE(rhog3(iz,ng_2d))
             r2 = AIMAG(rhog3(iz,ng_2d))
             !
             f1(:) =  DBLE(vg_f(iz,:))
             f2(:) = AIMAG(vg_f(iz,:))
             !
             for(:,it) = for(:,it)-r1*f1(:)-r2*f2(:)
          ENDDO
          !
       ENDDO
    ENDDO
    !
    for_g(:,:) = for_g(:,:) + for(:,:)
    DEALLOCATE( for, vg_f, vg_f_r )
    !
    !**** for gp==0****
    ng_2d = imill_2d(0,0)
    !
    IF ( ng_2d > 0 ) THEN
       ALLOCATE( vg_f(dfftp%nr3x,1), vg_f_r(dfftp%nr3x,1) )
       vg_f_r(:,1) = (0.d0,0.d0)
       !
       DO it = 1, nat
          tt = -fpi*zv(ityp(it))/sa
          zp = tau(3,it)
          IF (zp>at(3,3)*0.5) zp=zp-at(3,3)
          zp = zp*alat
          !
          DO iz = 1, dfftp%nr3
             k3 = iz-1
             IF (k3>dfftp%nr3/2) k3=iz-dfftp%nr3-1
             z=DBLE(k3)/DBLE(dfftp%nr3)*L
             ! bc3
             cc1=0.5d0*qe_erf(tmp*(z-zp))
             cc2=-0.5d0
             vg_f_r(iz,1) = tt*(cc1+cc2)
          ENDDO
          !
          CALL cft_1z( vg_f_r(:,1), 1, dfftp%nr3, dfftp%nr3, -1, vg_f(:,1) )
          !
          DO iz = 1, dfftp%nr3
             r1 =  DBLE(rhog3(iz,ng_2d))
             r2 = AIMAG(rhog3(iz,ng_2d))
             f1(3) =  DBLE(vg_f(iz,1))
             f2(3) = AIMAG(vg_f(iz,1))
             for_g(3,it)=for_g(3,it)-r1*f1(3)-r2*f2(3)
          ENDDO
       ENDDO
       !
       DEALLOCATE( vg_f, vg_f_r )
       !
    ENDIF ! if( ng_2d > 0 )
    !
    !**** sum short_range part and long_range part in local potential force 
    !**** at cartecian coordinate
    !
    DO it = 1, nat
       ! factor e2: hartree -> Ry.
       forcelc(1,it) = forcelc(1,it) &
                      +SUM(for_g(1:2,it)*bg(1,1:2))*SQRT(tpiba2)*omega*e2
       forcelc(2,it) = forcelc(2,it) &
                      +SUM(for_g(1:2,it)*bg(2,1:2))*SQRT(tpiba2)*omega*e2
       forcelc(3,it) = forcelc(3,it)+for_g(3,it)*omega*e2
    ENDDO
    !
    DEALLOCATE( for_g )
    DEALLOCATE( rhog3 )
    !
    RETURN
    !
  END SUBROUTINE esm_force_lc_bc3
  !
  !
  !----------------------------------------------------------------------------------
  SUBROUTINE esm_force_lc_bc4( aux, forcelc )
    !--------------------------------------------------------------------------------
    !! ESM local-potential derived force - bc4.
    !
    USE constants,        ONLY : tpi, fpi, e2
    USE gvect,            ONLY : ngm, mill
    USE cell_base,        ONLY : omega, alat, tpiba2, at, bg
    USE control_flags,    ONLY : gamma_only
    USE ions_base,        ONLY : zv, nat, tau, ityp
    USE fft_base,         ONLY : dfftp
    USE fft_scalar,       ONLY : cft_1z
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(IN) :: aux(dfftp%nnr)
    !! aux contains n(G) (input)   
    REAL(DP), INTENT(INOUT) :: forcelc(3,nat)
    !! lc force
    !
    ! ... local variables
    !
    INTEGER :: iz, it, k1, k2, k3, ng, n1, n2, n3, ng_2d
    REAL(DP), ALLOCATABLE :: for(:,:), for_g(:,:)
    REAL(DP) :: t(2), tt, gp, gp2, sa, z1, z0, pp, cc, ss, t1, &
                t2, z, zp, L, tmp, r1, r2, f1(3),       &
                f2(3), arg001, arg002, arg003, arg005,  &
                arg006, arg008, arg009, arg011, arg101, &
                arg102, arg103, arg104, arg106, arg107, &
                arg109, arg111, arg113, aaa, t3, alpha, beta,  &
                kappa, lambda, xi, chi
    COMPLEX(DP), ALLOCATABLE :: vg_f(:,:), vg_f_r(:,:), rhog3(:,:)
    COMPLEX(DP) :: c1(3), c2(3), cc1, cc2
    !
    ! Map to FFT mesh
    ALLOCATE( rhog3(dfftp%nr3,ngm_2d) )
    rhog3(:,:) = (0.d0,0.d0)
    !
    DO ng = 1, ngm
        n1 = mill(1,ng)
        n2 = mill(2,ng)
        ng_2d = imill_2d(n1,n2)
        n3 = mill(3,ng) + 1
        IF (n3 < 1) n3 = n3 + dfftp%nr3
        rhog3(n3,ng_2d) = aux(dfftp%nl(ng))
        IF (gamma_only .AND. n1==0 .AND. n2==0) THEN
           n3 = -mill(3,ng)+1
           IF (n3 < 1) n3 = n3+dfftp%nr3
           rhog3(n3,ng_2d) = aux(dfftp%nlm(ng))
        ENDIF  
    ENDDO
    !
    L = at(3,3)*alat
    sa = omega/L
    z0 = L/2.d0
    tmp = 1.d0
    z1 = z0+esm_w
    aaa = esm_a
    !
    ALLOCATE( for_g(3,nat) )
    for_g(:,:) = 0.d0
    !
    !**** for gp!=0 ****
    ALLOCATE( for(3,nat), vg_f(dfftp%nr3x,3), vg_f_r(dfftp%nr3x,3) )
    for(:,:) = 0.d0
    vg_f_r(:,:) = (0.d0,0.d0)
    !
    DO ng_2d = 1, ngm_2d
       k1 = mill_2d(1,ng_2d)
       k2 = mill_2d(2,ng_2d)
       IF (k1==0 .AND. k2==0) CYCLE
       !
       t(1:2) = k1*bg(1:2,1)+k2*bg(1:2,2)
       gp2 = SUM(t(:)*t(:))*tpiba2
       gp = SQRT(gp2)
       !
       DO it = 1, nat
          IF (gamma_only) THEN
             tt = -fpi*zv(ityp(it))/sa*2.d0
          ELSE
             tt = -fpi*zv(ityp(it))/sa
          ENDIF
          !
          pp = -tpi*(tau(1,it)*(k1*bg(1,1)+k2*bg(1,2)) &
                    +tau(2,it)*(k1*bg(2,1)+k2*bg(2,2)))
          cc = COS(pp)
          ss = SIN(pp)
          zp = tau(3,it)
          IF (zp>at(3,3)*0.5) zp=zp-at(3,3)
          zp = zp*alat
          !
          DO iz = 1, dfftp%nr3
             !
             k3 = iz-1
             IF (k3 > dfftp%nr3/2) k3 = iz-dfftp%nr3-1
             z = DBLE(k3)/DBLE(dfftp%nr3)*L
             ! bc4
             alpha = aaa+gp+SQRT(aaa**2+gp**2)
             beta  = aaa+gp-SQRT(aaa**2+gp**2)
             kappa = aaa-gp+SQRT(aaa**2+gp**2)
             xi    = aaa   +SQRT(aaa**2+gp**2)
             chi   = aaa   -SQRT(aaa**2+gp**2)
             lambda=        SQRT(aaa**2+gp**2)
             !
             arg001 = gp*(z-zp)
             arg002 =-gp*(z-zp)
             arg003 = gp*(z+zp-2.d0*z1)
             arg005 =-gp*(z1-zp)-xi*(z-z1)
             arg006 = aaa/2.d0/tmp**2*xi+gp*(z-z1)+xi*(z1-zp)
             arg008 = aaa/2.d0/tmp**2*xi+xi*(z1-zp)-chi*(z-z1)
             arg009 = aaa/2.d0/tmp**2*xi+xi*(z1-zp)-xi*(z-z1)
             arg011 = aaa/2.d0/tmp**2*chi+chi*(z1-zp)-xi*(z-z1)
             arg101 =  gp/2.d0/tmp+tmp*(z-zp)
             arg102 =  gp/2.d0/tmp-tmp*(z-zp)
             arg103 =  gp/2.d0/tmp+tmp*(z1-zp)
             arg104 =  gp/2.d0/tmp-tmp*(z1-zp)
             arg107 =  xi/2.d0/tmp+tmp*(z-zp)
             arg109 =  xi/2.d0/tmp+tmp*(z1-zp)
             arg111 = chi/2.d0/tmp+tmp*(z-zp)
             arg113 = chi/2.d0/tmp+tmp*(z1-zp)
             !
             IF (z < z1) THEN
                t1 = exp_erfc(arg001,arg101)-exp_erfc(arg001,arg103)
                t2 = exp_erfc(arg002,arg102) &
                    -kappa/alpha*exp_erfc(arg003,arg104)
                t3 = exp_erfc(arg006,arg109)/alpha
                !
                c1(1) = CMPLX(ss,-cc,KIND=DP)*(t1+t2)/4.d0/gp*k1
                c1(2) = CMPLX(ss,-cc,KIND=DP)*(t1+t2)/4.d0/gp*k2
                c2(1) = CMPLX(ss,-cc,KIND=DP)*t3/2.d0*k1
                c2(2) = CMPLX(ss,-cc,KIND=DP)*t3/2.d0*k2
                !
                t1 = exp_erfc(arg001,arg103)-exp_erfc(arg001,arg101)
                t2 = exp_erfc(arg002,arg102) &
                    -kappa/alpha*exp_erfc(arg003,arg104)
                t3 =-xi/alpha*exp_erfc(arg006,arg109)
                !
                c1(3) = CMPLX(cc,ss,KIND=DP)*(t1+t2)/4.d0
                c2(3) = CMPLX(cc,ss,KIND=DP)*t3/2.d0
             ELSE
                t1 = exp_erfc(arg011,arg113)-exp_erfc(arg011,arg111)
                t2 = exp_erfc(arg008,arg107) &
                    -beta/alpha*exp_erfc(arg009,arg109)
                t3 = exp_erfc(arg005,arg104)/alpha
                !
                c1(1) = CMPLX(ss,-cc,KIND=DP)*(t1+t2)/4.d0/lambda*k1
                c1(2) = CMPLX(ss,-cc,KIND=DP)*(t1+t2)/4.d0/lambda*k2
                c2(1) = CMPLX(ss,-cc,KIND=DP)*t3/2.d0*k1
                c2(2) = CMPLX(ss,-cc,KIND=DP)*t3/2.d0*k2
                !
                t1= chi*(exp_erfc(arg011,arg111)-exp_erfc(arg011,arg113))
                t2=-xi*(exp_erfc(arg008,arg107) &
                   +beta/alpha*exp_erfc(arg009,arg109))
                t3= gp/alpha*exp_erfc(arg005,arg104)
                !
                c1(3) = CMPLX(cc, ss,KIND=DP)*(t1+t2)/4.d0/lambda
                c2(3) = CMPLX(cc, ss,KIND=DP)*t3/2.d0
             ENDIF
             vg_f_r(iz,:) = tt*(c1(:)+c2(:))
          ENDDO
          !
          CALL cft_1z( vg_f_r(:,1), 1, dfftp%nr3, dfftp%nr3, -1, vg_f(:,1) )
          CALL cft_1z( vg_f_r(:,2), 1, dfftp%nr3, dfftp%nr3, -1, vg_f(:,2) )
          CALL cft_1z( vg_f_r(:,3), 1, dfftp%nr3, dfftp%nr3, -1, vg_f(:,3) )
          !
          DO iz = 1, dfftp%nr3
             r1 =  DBLE(rhog3(iz,ng_2d))
             r2 = AIMAG(rhog3(iz,ng_2d))
             f1(:) =  DBLE(vg_f(iz,:))
             f2(:) = AIMAG(vg_f(iz,:))
             for(:,it) = for(:,it)-r1*f1(:)-r2*f2(:)
          ENDDO
       ENDDO
    ENDDO
    !
    for_g(:,:) = for_g(:,:) + for(:,:)
    !
    DEALLOCATE( for, vg_f, vg_f_r )
    !
    !**** for gp==0****
    ng_2d = imill_2d(0,0)
    !
    IF ( ng_2d > 0 ) THEN
       !
       ALLOCATE( vg_f(dfftp%nr3x,1), vg_f_r(dfftp%nr3x,1) )
       vg_f_r(:,1)=(0.d0,0.d0)
       !
       DO it = 1, nat
          tt = -fpi*zv(ityp(it))/sa
          zp = tau(3,it)
          IF (zp > at(3,3)*0.5) zp = zp-at(3,3)
          zp = zp*alat
          !
          DO iz = 1, dfftp%nr3
             k3 = iz-1
             IF (k3 > dfftp%nr3/2) k3 = iz-dfftp%nr3-1
             z = DBLE(k3)/DBLE(dfftp%nr3)*L
             ! bc4
             arg006 = aaa**2/tmp**2+2.d0*aaa*(z1-zp)
             arg101 = tmp*(z-zp)
             arg102 = tmp*(z1-zp)
             arg104 = aaa/tmp+tmp*(z-zp)
             arg106 = aaa/tmp+tmp*(z1-zp)
             IF (z < z1)THEN
                cc1 = 0.5d0*(qe_erf(arg101)-qe_erf(arg102))
                cc2 =-0.5d0*exp_erfc(arg006,arg106)
             ELSE
                cc1 = 0.d0
                cc2 =-0.5d0*exp_erfc(arg006,arg104)
             ENDIF
             vg_f_r(iz,1) = tt*(cc1+cc2)
          ENDDO
          !
          CALL cft_1z( vg_f_r(:,1), 1, dfftp%nr3, dfftp%nr3, -1, vg_f(:,1) )
          !
          DO iz = 1, dfftp%nr3
             r1 =  DBLE(rhog3(iz,ng_2d))
             r2 = AIMAG(rhog3(iz,ng_2d))
             f1(3) =  DBLE(vg_f(iz,1))
             f2(3) = AIMAG(vg_f(iz,1))
             for_g(3,it) = for_g(3,it)-r1*f1(3)-r2*f2(3)
          ENDDO
          !
       ENDDO
       !
       DEALLOCATE( vg_f, vg_f_r )
       !
    ENDIF ! if( ng_2d > 0 )
    !
    !
    !**** sum short_range part and long_range part in local potential force 
    !**** at cartecian coordinate
    !
    DO it = 1, nat
       ! factor e2: hartree -> Ry.
       forcelc(1,it) = forcelc(1,it) &
                      +SUM(for_g(1:2,it)*bg(1,1:2))*SQRT(tpiba2)*omega*e2
       forcelc(2,it) = forcelc(2,it) &
                      +SUM(for_g(1:2,it)*bg(2,1:2))*SQRT(tpiba2)*omega*e2
       forcelc(3,it) = forcelc(3,it)+for_g(3,it)*omega*e2
    ENDDO
    !
    DEALLOCATE( for_g )
    DEALLOCATE( rhog3 )
    !
    RETURN
    !
  END SUBROUTINE esm_force_lc_bc4
  !
  !
  !  ...  ESM FINAL PRINTOUT SUBROUTINE
  !
  !--------------------------------------------------------------------------
  SUBROUTINE esm_printpot( rhog )
    !-------------------------------------------------------------------------
    !! Prints out vlocal and vhartree to stdout once electrons are converged.  
    !! Format: z, rho(r), v_hartree, v_local, (v_hartree + v_local).
    !
    USE kinds,         ONLY : DP
    USE gvect,         ONLY : ngm, mill, igtongl
    USE constants,     ONLY : pi, sqrtpm1, tpi, fpi, e2
    USE constants,     ONLY : AUTOEV, BOHR_RADIUS_ANGS
    USE cell_base,     ONLY : omega, alat, at, tpiba, bg
    USE ions_base,     ONLY : zv, nat, tau, ityp, ntyp => nsp
    USE vlocal,        ONLY : strf, vloc
    USE control_flags, ONLY : gamma_only
    USE fft_base,      ONLY : dfftp
    USE fft_scalar,    ONLY : cft_1z
    USE io_files,      ONLY : prefix, tmp_dir
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(IN) :: rhog(ngm)
    !! n(G)
    !
    ! ... local variables
    !
    INTEGER :: ig, iga, igb, igz, iz, ia
    REAL(DP) :: L, S, z0, z1, z, gz, alpha, salp
    REAL(DP) :: Qa, ra(2), za
    COMPLEX(DP), PARAMETER :: ci = DCMPLX(0.0d0, 1.0d0)
    COMPLEX(DP) :: rg3, vg3, sum1c, sum2c
    COMPLEX(DP) :: expimgpr, experfcm, experfcp, dexperfcm_dgp, dexperfcp_dgp
    !
    COMPLEX(DP), ALLOCATABLE :: rho0r(:), rho0g(:)
    COMPLEX(DP), ALLOCATABLE :: Vhar0r(:), Vhar0g(:)
    COMPLEX(DP), ALLOCATABLE :: Vloc0r(:), Vloc0g(:)
    CHARACTER(LEN=256) :: esm1_file = 'os.esm1'
    !
    IF ( imill_2d(0,0) == 0 ) RETURN
    !
    !****For gp=0 case ****
    !
    ! cell settings
    L  = at(3,3)*alat
    S  = omega/L
    z0 = L/2.d0
    z1 = z0 + esm_w
    alpha = 1.0d0
    salp = SQRT(alpha)
    !
    ALLOCATE( rho0r(dfftp%nr3),  rho0g(dfftp%nr3)  )
    ALLOCATE( Vhar0r(dfftp%nr3), Vhar0g(dfftp%nr3) )
    ALLOCATE( Vloc0r(dfftp%nr3), Vloc0g(dfftp%nr3) )
    !
    rho0g(:) = (0.d0,0.d0)
    !
    !---- calculate density potential
    DO ig = 1, ngm
       iga = mill(1,ig)
       igb = mill(2,ig)
       igz = mill(3,ig)+1
       !
       IF ( .NOT. (iga==0 .AND. igb==0) ) CYCLE
       !
       IF (igz < 1) THEN
          igz = igz + dfftp%nr3
       ENDIF
       !
       rg3 = rhog(ig)
       rho0g(igz) = rg3
       !
       IF ( gamma_only .AND. iga==0 .AND. igb==0 ) THEN
          igz = 1-mill(3,ig)
          IF (igz < 1) THEN
             igz = igz + dfftp%nr3
          ENDIF
          rho0g(igz) = CONJG(rg3)
       ENDIF
    ENDDO ! ig
    !
    CALL cft_1z( rho0g, 1, dfftp%nr3, dfftp%nr3, +1, rho0r )
    !
    !--- calculate Hartree potential
    Vhar0g(:) = 0.0d0
    !
    DO igz=2, dfftp%nr3
       IF (igz <= dfftp%nr3/2) THEN
          gz = DBLE(igz-1)*tpi/L
       ELSE
          gz = DBLE(igz-1-dfftp%nr3)*tpi/L
       ENDIF
       !
       rg3 = rho0g(igz)
       Vhar0g(igz) = fpi*rg3/gz**2
    ENDDO ! igz
    !
    CALL cft_1z( Vhar0g, 1, dfftp%nr3, dfftp%nr3, +1, Vhar0r )
    !
    ! summations over gz
    sum1c = (0.d0,0.d0)
    sum2c = (0.d0,0.d0)
    !
    DO igz = 2, dfftp%nr3
       IF (igz <= dfftp%nr3/2) THEN
          gz = DBLE(igz-1)*tpi/L
       ELSE
          gz = DBLE(igz-1-dfftp%nr3)*tpi/L
       ENDIF
       !
       rg3 = rho0g(igz)
       sum1c = sum1c + rg3*ci*COS(gz*z0)/gz
       sum2c = sum2c + rg3*COS(gz*z0)/gz**2
    ENDDO ! igz
    !
    rg3 = rho0g(1)
    !
    DO iz = 1, dfftp%nr3
       !
       IF( iz-1 > dfftp%nr3/2 ) THEN
          z = DBLE(iz-1-dfftp%nr3)/DBLE(dfftp%nr3)*L
       ELSE
          z = DBLE(iz-1)/DBLE(dfftp%nr3)*L
       ENDIF
       !
       ! BC1 terms
       Vhar0r(iz) = Vhar0r(iz) &
            - tpi*z**2*rg3     &
            - tpi*z0**2*rg3    &
            - fpi*z*sum1c      &
            - fpi*sum2c
       !
       IF ( esm_bc == 'bc2' ) THEN
          !! BC2 terms
          Vhar0r(iz) = Vhar0r(iz) &
               + tpi*z1*2*z0 * rg3 - tpi*(-z/z1)*2*z0*sum1c
       ELSEIF ( esm_bc == 'bc3' ) THEN
          !! BC3 terms
          Vhar0r(iz) = Vhar0r(iz) &
               - tpi*(z-2*z1)*2*z0 * rg3 + fpi*z0*sum1c
       ENDIF
       !
    ENDDO ! iz
    !
    !--- calculate local potential
    ! short range
    Vloc0g(:) = 0.0d0
    !
    DO ig = 1, ngm
       iga = mill(1,ig)
       igb = mill(2,ig)
       !
       IF ( .NOT. (iga==0 .AND. igb==0) ) CYCLE
       !
       igz = mill(3,ig)+1
       IF (igz < 1) THEN
          igz = igz + dfftp%nr3
       ENDIF
       !
       vg3 = 0.0d0
       !
       DO ia = 1, ntyp
          vg3 = vg3 + vloc(igtongl(ig),ia)*strf(ig,ia)/e2
       ENDDO
       !
       Vloc0g(igz) = vg3
       !
       IF ( gamma_only .AND. iga==0 .AND. igb==0 ) THEN
          igz = 1-mill(3,ig)
          IF (igz < 1) THEN
             igz = igz + dfftp%nr3
          ENDIF
          !
          Vloc0g(igz) = CONJG(vg3)
       ENDIF
    ENDDO
    !
    CALL cft_1z( Vloc0g, 1, dfftp%nr3, dfftp%nr3, +1, Vloc0r )
    !
    ! long range
    DO iz = 1, dfftp%nr3
       !
       IF ( iz-1 > dfftp%nr3/2 ) THEN
          z = DBLE(iz-1-dfftp%nr3)/DBLE(dfftp%nr3)*L
       ELSE
          z = DBLE(iz-1)/DBLE(dfftp%nr3)*L
       ENDIF
       !
       IF ( esm_bc == 'bc2' ) THEN
          Vloc0r(iz) = Vloc0r(iz) + (z1-z)*esm_efield/e2
       ENDIF
       !
       DO ia = 1, nat
          Qa = (-1.0d0)*zv(ityp(ia))
          ra(1:2) = tau(1:2,ia)*alat
          za = tau(3,ia)*alat
          IF ( za > L*0.5d0 ) THEN
             za = za - L
          ENDIF
          !
          ! BC1 terms
          Vloc0r(iz) = Vloc0r(iz) - tpi * Qa/S &
               * ( (z-za)*qe_erf(salp*(z-za)) &
               + EXP(-alpha*(z-za)**2)*sqrtpm1/salp )
          !
          IF ( esm_bc == 'bc2' ) THEN
             !! BC2 terms
             Vloc0r(iz) = Vloc0r(iz) &
                  + tpi * Qa/S * ( -z*za + z1*z1 )/z1
          ELSEIF( esm_bc == 'bc3' ) THEN
             !! BC3 terms
             Vloc0r(iz) = Vloc0r(iz) &
                  + tpi * Qa/S * (-z+2*z1-za)
          ENDIF

       ENDDO ! ia
    ENDDO ! iz
    !
    !---- output potentials
    esm1_file = TRIM( tmp_dir ) // TRIM( prefix ) // ".esm1"
    OPEN( UNIT = 4, FILE = esm1_file, STATUS = "UNKNOWN", &
          ACTION = "WRITE" )
    WRITE( UNIT = 4, FMT = 9050 )
    !
    DO iz = dfftp%nr3/2+2, dfftp%nr3
       z = DBLE(iz-1-dfftp%nr3)/DBLE(dfftp%nr3)*L
       !
       WRITE( UNIT = 4, FMT = 9051 ) z*BOHR_RADIUS_ANGS, &
                     DBLE(rho0r(iz))*S/BOHR_RADIUS_ANGS, &
                     DBLE(Vhar0r(iz))*AUTOEV, &
                     DBLE(Vloc0r(iz))*AUTOEV, &
                     DBLE(Vhar0r(iz))*AUTOEV + DBLE(Vloc0r(iz))*AUTOEV
    ENDDO
    !
    DO iz = 1, dfftp%nr3/2+1
       z = DBLE(iz-1)/DBLE(dfftp%nr3)*L
       !
       WRITE( UNIT = 4, FMT = 9051 ) z*BOHR_RADIUS_ANGS, &
                     DBLE(rho0r(iz))*S/BOHR_RADIUS_ANGS, &
                     DBLE(Vhar0r(iz))*AUTOEV, &
                     DBLE(Vloc0r(iz))*AUTOEV, &
                     DBLE(Vhar0r(iz))*AUTOEV + DBLE(Vloc0r(iz))*AUTOEV
    ENDDO
    !
    CLOSE( UNIT = 4 )
    !
    DEALLOCATE( rho0r,  rho0g  )
    DEALLOCATE( Vhar0r, Vhar0g )
    DEALLOCATE( Vloc0r, Vloc0g )
    !
    RETURN
    !
9050 FORMAT( '#z (A)',2X,'Tot chg (e/A)',2X,'Avg v_hartree (eV)',2X,&
         &'Avg v_local (eV)',2x,'Avg v_hart+v_loc (eV)' )
9051 FORMAT( F6.2,F20.7,F20.7,F18.7,F18.7 )
    !
  END SUBROUTINE esm_printpot
  !
  !
  !  ...  ESM SUMMARY PRINTOUT SUBROUTINE
  !
  !-----------------------------------------------------------------------
  SUBROUTINE esm_summary()
     !--------------------------------------------------------------------
     !! Prints summary of ESM parameters to stdout.
     !  
     USE io_global,        ONLY : stdout, ionode  
     USE constants,        ONLY : rytoev, BOHR_RADIUS_ANGS  
     !  
     IMPLICIT NONE  
     !  
     IF ( .NOT. ionode ) RETURN  
     !  
     WRITE( UNIT = stdout,                                          &  
            FMT  = '(/,5x, "Effective Screening Medium Method",     &  
                    &/,5x, "=================================")' )  
     !  
     SELECT CASE( TRIM(esm_bc) )  
     CASE( 'pbc' )  
        WRITE( UNIT = stdout,                                     &  
             FMT  = '(5x, "Ordinary Periodic Boundary Conditions")' )  
     CASE( 'bc1' )  
        WRITE( UNIT = stdout,                                     &  
             FMT  = '(5x, "Boundary Conditions: Vacuum-Slab-Vacuum")' )  
     CASE( 'bc2' )  
        WRITE( UNIT = stdout,                                     &  
             FMT  = '(5x, "Boundary Conditions: Metal-Slab-Metal")' )  
     CASE( 'bc3' )  
        WRITE( UNIT = stdout,                                     &  
             FMT  = '(5x, "Boundary Conditions: Vacuum-Slab-Metal")' )  
     CASE( 'bc4' )  
        WRITE( UNIT = stdout,                                     &  
             FMT  = '(5x, "Boundary Conditions: Vacuum-Slab-smooth ESM)")' )  
     END SELECT  
     !  
     IF ( esm_efield /= 0.0_DP ) THEN  
        WRITE( UNIT = stdout, FMT = 9051 ) esm_efield  
     ENDIF  
     !  
     IF ( esm_w /= 0.0_DP ) THEN  
        WRITE( UNIT = stdout, FMT = 9052 ) esm_w*BOHR_RADIUS_ANGS, esm_w  
     ENDIF  
     !  
     IF (esm_bc .EQ. 'bc4') THEN  
       WRITE( UNIT = stdout, FMT = 9054 ) esm_a  
     ENDIF  
     !  
     WRITE( UNIT = stdout, FMT = 9053 ) esm_nfit  
     !  
     WRITE( stdout, * )  
     !  
9051 FORMAT( '     field strength                   = ', F8.2,' Ry/a.u.')  
9052 FORMAT( '     ESM offset from cell edge        = ', F8.2,' A' &  
            /'                                      = ', F8.2,' a.u.')  
9053 FORMAT( '     grid points for fit at edges     = ', I8,' ')  
9054 FORMAT( '     smoothness parameter             = ', F8.2,' 1/a.u.' )  
     !
  END SUBROUTINE esm_summary
  !
  !
  !---------------------------------------------------------------------------
  FUNCTION vl11j0( gp, aaa, tmp, z1, z, zp, rxy )
    !-------------------------------------------------------------------------
    !
    USE kinds,  ONLY : DP
    !
    IMPLICIT NONE
    !
    REAL(DP) :: vl11j0
    REAL(DP), INTENT(IN) :: gp, aaa, tmp, z1, z, zp, rxy
    !
    ! ... local variables
    !
    REAL(DP) :: alpha, kappa, xi, t1, t2, t3, arg001, arg002,   &
                arg003, arg006, arg101, arg102, arg103, arg104, &
                arg109
    !
    alpha = aaa+gp+SQRT(aaa**2+gp**2)
    kappa = aaa-gp+SQRT(aaa**2+gp**2)
    xi    = aaa   +SQRT(aaa**2+gp**2)
    !
    arg001 =  gp*(z-zp)
    arg002 = -gp*(z-zp)
    arg003 =  gp*(z+zp-2.d0*z1)
    arg006 =  aaa/2.d0/tmp**2*xi+gp*(z-z1)+xi*(z1-zp)
    arg101 =   gp/2.d0/tmp+tmp*( z-zp)
    arg102 =   gp/2.d0/tmp-tmp*( z-zp)
    arg103 =   gp/2.d0/tmp+tmp*(z1-zp)
    arg104 =   gp/2.d0/tmp-tmp*(z1-zp)
    arg109 =   xi/2.d0/tmp+tmp*(z1-zp)
    !
    t1 = -EXP(arg003)*kappa/alpha
    t2 =  exp_erfc(arg001,arg101)-exp_erfc(arg001,arg103) &
         +exp_erfc(arg002,arg102)-kappa/alpha*exp_erfc(arg003,arg104)
    t3 =  exp_erfc(arg006,arg109)/alpha
    !
    vl11j0 = (t1/2.d0-(t2/4.d0+gp*t3/2.d0))*dbesj0(gp*rxy)
    !
    RETURN
    !
  END FUNCTION vl11j0
  !
  !
  !------------------------------------------------------------------------
  FUNCTION vl11j1( gp, aaa, tmp, z1, z, zp, rxy )
    !-----------------------------------------------------------------------
    !
    USE kinds,  ONLY : DP
    !
    IMPLICIT NONE
    !
    REAL(DP) :: vl11j1
    REAL(DP), INTENT(IN) :: gp, aaa, tmp, z1, z, zp, rxy
    !
    ! ... local variables
    !
    REAL(DP) :: alpha, kappa, xi, t1, t2, t3, arg001, arg002,   &
                arg003, arg006, arg007, arg101, arg102, arg103, &
                arg104, arg105, arg106, arg109, arg110
    !
    alpha = aaa+gp+SQRT(aaa**2+gp**2)
    kappa = aaa-gp+SQRT(aaa**2+gp**2)
    xi    = aaa   +SQRT(aaa**2+gp**2)
    !
    arg001 =  gp*(z-zp)
    arg002 = -gp*(z-zp)
    arg003 =  gp*(z+zp-2.d0*z1)
    arg006 =  aaa/2.d0/tmp**2*xi+gp*( z-z1)+xi*(z1-zp)
    arg007 =  aaa/2.d0/tmp**2*xi-gp*(z1-zp)-xi*( z-z1)
    arg101 =   gp/2.d0/tmp+tmp*( z-zp)
    arg102 =   gp/2.d0/tmp-tmp*( z-zp)
    arg103 =   gp/2.d0/tmp+tmp*(z1-zp)
    arg104 =   gp/2.d0/tmp-tmp*(z1-zp)
    arg105 =   gp/2.d0/tmp+tmp*( z-z1)
    arg106 =   gp/2.d0/tmp-tmp*( z-z1)
    arg109 =   xi/2.d0/tmp+tmp*(z1-zp)
    arg110 =   xi/2.d0/tmp-tmp*( z-z1)
    !
    t1 = -EXP(arg003)*kappa/alpha
    t2 =  exp_erfc(arg001,arg101)-exp_erfc(arg001,arg103) &
         +exp_erfc(arg002,arg102)-kappa/alpha*exp_erfc(arg003,arg104) &
         +exp_erfc(arg006,arg109)*2.d0*gp/alpha
    t3 =  exp_erfc(arg002,arg102)-exp_erfc(arg002,arg106) &
         +exp_erfc(arg001,arg101)-kappa/alpha*exp_erfc(arg003,arg105) &
         +exp_erfc(arg007,arg110)*2.d0*gp/alpha
    !
    vl11j1 = gp*(t1-(t2+t3)/4.d0)*dbesj1(gp*rxy)
    !
    RETURN
    !
  END FUNCTION vl11j1
  !
  !
  !----------------------------------------------------------------------
  FUNCTION vl12j0( gp, aaa, tmp, z1, z, zp, rxy )
    !--------------------------------------------------------------------- 
    !
    USE kinds, ONLY : DP
    !
    IMPLICIT NONE
    !
    REAL(DP) :: vl12j0
    REAL(DP), INTENT(IN) :: gp, aaa, tmp, z1, z, zp, rxy
    !
    ! ... local variables
    !
    REAL(DP) :: alpha, kappa, xi, t1, t2, t3, arg001, arg002,   &
                arg003, arg004, arg006, arg101, arg102, arg103, &
                arg104, arg109
    !
    alpha = aaa+gp+SQRT(aaa**2+gp**2)
    kappa = aaa-gp+SQRT(aaa**2+gp**2)
    xi    = aaa   +SQRT(aaa**2+gp**2)
    !
    arg001 =  gp*(z-zp)
    arg002 = -gp*(z-zp)
    arg003 =  gp*(z+zp-2.d0*z1)
    arg004 =  gp*(z-z1)+xi*(z1-zp)
    arg006 =  aaa/2.d0/tmp**2*xi+gp*(z-z1)+xi*(z1-zp)
    arg101 =   gp/2.d0/tmp+tmp*(z-zp)
    arg102 =   gp/2.d0/tmp-tmp*(z-zp)
    arg103 =   gp/2.d0/tmp+tmp*(z1-zp)
    arg104 =   gp/2.d0/tmp-tmp*(z1-zp)
    arg109 =   xi/2.d0/tmp+tmp*(z1-zp)
    !
    t1 =  EXP(arg004)/alpha
    t2 =  exp_erfc(arg001,arg101)-exp_erfc(arg001,arg103) &
         +exp_erfc(arg002,arg102)-kappa/alpha*exp_erfc(arg003,arg104)
    t3 =  exp_erfc(arg006,arg109)/alpha
    !
    vl12j0 = (gp*t1-(t2/4.d0+gp*t3/2.d0))*dbesj0(gp*rxy)
    !
  END FUNCTION vl12j0
  !
  !
  !----------------------------------------------------------------------
  FUNCTION vl12j1( gp, aaa, tmp, z1, z, zp, rxy )
    !---------------------------------------------------------------------
    !
    USE kinds,  ONLY : DP
    !
    IMPLICIT NONE
    !
    REAL(DP) :: vl12j1
    REAL(DP), INTENT(IN) :: gp, aaa, tmp, z1, z, zp, rxy
    !
    ! ... local variables
    !
    REAL(DP) :: alpha, beta, kappa, xi, chi, lambda, t1, t2, t3, &
                arg001, arg002, arg003, arg004, arg006, arg009,  &
                arg010, arg012, arg101, arg102, arg103, arg104,  &
                arg105, arg108, arg109, arg110, arg112, arg114
    !
    alpha  = aaa+gp+SQRT(aaa**2+gp**2)
    beta   = aaa+gp-SQRT(aaa**2+gp**2)
    kappa  = aaa-gp+SQRT(aaa**2+gp**2)
    xi     = aaa   +SQRT(aaa**2+gp**2)
    chi    = aaa   -SQRT(aaa**2+gp**2)
    lambda =        SQRT(aaa**2+gp**2)
    !
    arg001 =  gp*(z-zp)
    arg002 = -gp*(z-zp)
    arg003 =  gp*(z+zp-2.d0*z1)
    arg004 =  gp*(z-z1)+xi*(z1-zp)
    arg006 =  aaa/2.d0/tmp**2*xi+gp*(z-z1)+xi*(z1-zp)
    arg009 =  aaa/2.d0/tmp**2*xi+xi*(z1-zp)-xi*(z-z1)
    arg010 =  aaa/2.d0/tmp**2*xi+chi*(z1-zp)-xi*(z-z1)
    arg012 =  aaa/2.d0/tmp**2*chi+xi*(z1-zp)-chi*(z-z1)
    arg101 =   gp/2.d0/tmp+tmp*(z-zp)
    arg102 =   gp/2.d0/tmp-tmp*(z-zp)
    arg103 =   gp/2.d0/tmp+tmp*(z1-zp)
    arg104 =   gp/2.d0/tmp-tmp*(z1-zp)
    arg105 =   gp/2.d0/tmp+tmp*(z-z1)
    arg108 =   xi/2.d0/tmp-tmp*(z-zp)
    arg109 =   xi/2.d0/tmp+tmp*(z1-zp)
    arg110 =   xi/2.d0/tmp-tmp*(z-z1)
    arg112 =  chi/2.d0/tmp-tmp*(z-zp)
    arg114 =  chi/2.d0/tmp-tmp*(z-z1)
    !
    t1 =  EXP(arg004)/alpha
    t2 =  exp_erfc(arg001,arg101)-exp_erfc(arg001,arg103) &
         +exp_erfc(arg002,arg102)-kappa/alpha*exp_erfc(arg003,arg104) &
         +exp_erfc(arg006,arg109)*2.d0*gp/alpha
    t3 =  exp_erfc(arg012,arg114)-exp_erfc(arg012,arg112) &
         +exp_erfc(arg010,arg108)-beta/alpha*exp_erfc(arg009,arg110) &
         +exp_erfc(arg004,arg105)*2.d0*lambda/alpha
    !
    vl12j1 = (2.d0*t1*gp**2-(gp*t2+gp**2*t3/lambda)/4.d0)*dbesj1(gp*rxy)
    !
  END FUNCTION vl12j1
  !
  !
  !-------------------------------------------------------------------------
  FUNCTION vl21j0( gp, aaa, tmp, z1, z, zp, rxy )
    !-----------------------------------------------------------------------
    !
    USE kinds,  ONLY : DP
    !
    IMPLICIT NONE
    !
    REAL(DP) :: vl21j0
    REAL(DP),INTENT(IN) :: gp, aaa, tmp, z1, z, zp, rxy
    !
    ! ... local variables
    !
    REAL(DP) :: alpha, beta, xi, chi, lambda, t1, t2, t3, arg005, &
                arg008, arg009, arg011, arg104, arg107, arg109,   &
                arg111, arg113
    !
    alpha = aaa+gp+SQRT(aaa**2+gp**2)
    beta  = aaa+gp-SQRT(aaa**2+gp**2)
    xi    = aaa   +SQRT(aaa**2+gp**2)
    chi   = aaa   -SQRT(aaa**2+gp**2)
    lambda=        SQRT(aaa**2+gp**2)
    !
    arg005 = -gp*(z1-zp)-xi*(z-z1)
    arg008 =  aaa/2.d0/tmp**2*xi + xi*(z1-zp)-chi*(z-z1)
    arg009 =  aaa/2.d0/tmp**2*xi + xi*(z1-zp)- xi*(z-z1)
    arg011 =  aaa/2.d0/tmp**2*chi+chi*(z1-zp)- xi*(z-z1)
    arg104 =   gp/2.d0/tmp-tmp*(z1-zp)
    arg107 =   xi/2.d0/tmp+tmp*( z-zp)
    arg109 =   xi/2.d0/tmp+tmp*(z1-zp)
    arg111 =  chi/2.d0/tmp+tmp*( z-zp)
    arg113 =  chi/2.d0/tmp+tmp*(z1-zp)
    !
    t1 =  EXP(arg005)/alpha
    t2 =  exp_erfc(arg011,arg113)-exp_erfc(arg011,arg111) &
         +exp_erfc(arg008,arg107)-beta/alpha*exp_erfc(arg009,arg109)
    t3 =  exp_erfc(arg005,arg104)/alpha
    !
    vl21j0 = gp*(t1-(t2/4.d0/lambda+t3/2.d0))*dbesj0(gp*rxy)
    !
  END FUNCTION vl21j0
  !
  !
  !-----------------------------------------------------------------------
  FUNCTION vl21j1( gp, aaa, tmp, z1, z, zp, rxy )
    !---------------------------------------------------------------------
    USE kinds,  ONLY : DP
    !
    IMPLICIT NONE
    !
    REAL(DP) :: vl21j1
    REAL(DP), INTENT(IN) :: gp, aaa, tmp, z1, z, zp, rxy
    !
    ! ... local variables
    !
    REAL(DP) :: alpha, beta, kappa, xi, chi, lambda, t1, t2, t3, &
                arg001, arg002, arg003, arg005, arg007,arg008,   &
                arg009, arg011, arg101, arg102, arg104, arg105,  &
                arg106, arg107, arg109, arg110, arg111, arg113
    !
    alpha = aaa+gp+SQRT(aaa**2+gp**2)
    beta  = aaa+gp-SQRT(aaa**2+gp**2)
    kappa = aaa-gp+SQRT(aaa**2+gp**2)
    xi    = aaa   +SQRT(aaa**2+gp**2)
    chi   = aaa   -SQRT(aaa**2+gp**2)
    lambda=        SQRT(aaa**2+gp**2)
    !
    arg001 =  gp*(z-zp)
    arg002 = -gp*(z-zp)
    arg003 =  gp*(z+zp-2.d0*z1)
    arg005 = -gp*(z1-zp)-xi*(z-z1)
    arg007 =  aaa/2.d0/tmp**2*xi - gp*(z1-zp)- xi*(z-z1)
    arg008 =  aaa/2.d0/tmp**2*xi + xi*(z1-zp)-chi*(z-z1)
    arg009 =  aaa/2.d0/tmp**2*xi + xi*(z1-zp)- xi*(z-z1)
    arg011 =  aaa/2.d0/tmp**2*chi+chi*(z1-zp)- xi*(z-z1)
    arg101 =   gp/2.d0/tmp+tmp*( z-zp)
    arg102 =   gp/2.d0/tmp-tmp*( z-zp)
    arg104 =   gp/2.d0/tmp-tmp*(z1-zp)
    arg105 =   gp/2.d0/tmp+tmp*( z-z1)
    arg106 =   gp/2.d0/tmp-tmp*( z-z1)
    arg107 =   xi/2.d0/tmp+tmp*( z-zp)
    arg109 =   xi/2.d0/tmp+tmp*(z1-zp)
    arg110 =   xi/2.d0/tmp-tmp*( z-z1)
    arg111 =  chi/2.d0/tmp+tmp*( z-zp)
    arg113 =  chi/2.d0/tmp+tmp*(z1-zp)
    !
    t1 =  EXP(arg005)/alpha
    t2 =  exp_erfc(arg011,arg113)-exp_erfc(arg011,arg111) &
         +exp_erfc(arg008,arg107)-beta/alpha*exp_erfc(arg009,arg109)  &
         +exp_erfc(arg005,arg104)*2.d0*lambda/alpha
    t3 =  exp_erfc(arg002,arg102)-exp_erfc(arg002,arg106) &
         +exp_erfc(arg001,arg101)-kappa/alpha*exp_erfc(arg003,arg105) &
         +exp_erfc(arg007,arg110)*2.d0*gp/alpha
    !
    vl21j1 = (2.d0*t1*gp**2-(gp**2*t2/lambda+gp*t3)/4.d0)*dbesj1(gp*rxy)
    !
  END FUNCTION vl21j1
  !
  !
  !--------------------------------------------------------------------------
  FUNCTION vl22j0( gp, aaa, tmp, z1, z, zp, rxy )
    !------------------------------------------------------------------------
    !
    USE kinds,  ONLY : DP
    !
    IMPLICIT NONE
    !
    REAL(DP) :: vl22j0
    REAL(DP), INTENT(IN) :: gp, aaa, tmp, z1, z, zp, rxy
    !
    ! ... local variables
    !
    REAL(DP) :: alpha, beta, xi, chi, lambda, t1, t2, t3, arg000, &
                arg005, arg008, arg009, arg011, arg104, arg107,   &
                arg109, arg111, arg113
    !
    alpha = aaa+gp+SQRT(aaa**2+gp**2)
    beta  = aaa+gp-SQRT(aaa**2+gp**2)
    xi    = aaa   +SQRT(aaa**2+gp**2)
    chi   = aaa   -SQRT(aaa**2+gp**2)
    lambda=        SQRT(aaa**2+gp**2)
    !
    arg000 = -xi*(z+zp-2.d0*z1)
    arg005 = -gp*(z1-zp)-xi*(z-z1)
    arg008 =  aaa/2.d0/tmp**2*xi+xi*(z1-zp)-chi*(z-z1)
    arg009 =  aaa/2.d0/tmp**2*xi+xi*(z1-zp)-xi*(z-z1)
    arg011 =  aaa/2.d0/tmp**2*chi+chi*(z1-zp)-xi*(z-z1)
    arg104 =   gp/2.d0/tmp-tmp*(z1-zp)
    arg107 =   xi/2.d0/tmp+tmp*( z-zp)
    arg109 =   xi/2.d0/tmp+tmp*(z1-zp)
    arg111 =  chi/2.d0/tmp+tmp*( z-zp)
    arg113 =  chi/2.d0/tmp+tmp*(z1-zp)
    !
    t1 = -EXP(arg000)*beta/alpha
    t2 =  exp_erfc(arg011,arg113)-exp_erfc(arg011,arg111) &
         +exp_erfc(arg008,arg107)-beta/alpha*exp_erfc(arg009,arg109)
    t3 =  exp_erfc(arg005,arg104)/alpha
    !
    vl22j0 = gp*(t1/2.d0/lambda-(t2/4.d0/lambda+t3/2.d0))*dbesj0(gp*rxy)
    !
    RETURN
    !
  END FUNCTION vl22j0
  !
  !
  !--------------------------------------------------------------------------
  FUNCTION vl22j1( gp, aaa, tmp, z1, z, zp, rxy )
    !-------------------------------------------------------------------------
    !
    USE kinds,  ONLY : DP
    !
    IMPLICIT NONE
    !
    REAL(DP) :: vl22j1
    REAL(DP), INTENT(IN) :: gp, aaa, tmp, z1, z, zp, rxy
    !
    ! ... local variables
    !
    REAL(DP) :: alpha, beta, xi, chi, lambda, t1, t2, t3, arg000, &
                arg004, arg005, arg008, arg009, arg010, arg011,   &
                arg012, arg104, arg105, arg107, arg108, arg109,   &
                arg110, arg111, arg112, arg113, arg114
    !
    alpha = aaa+gp+SQRT(aaa**2+gp**2)
    beta  = aaa+gp-SQRT(aaa**2+gp**2)
    xi    = aaa   +SQRT(aaa**2+gp**2)
    chi   = aaa   -SQRT(aaa**2+gp**2)
    lambda=        SQRT(aaa**2+gp**2)
    !
    arg000 = -xi*(z+zp-2.d0*z1)
    arg004 =  gp*(z-z1)+xi*(z1-zp)
    arg005 = -gp*(z1-zp)-xi*(z-z1)
    arg008 =  aaa/2.d0/tmp**2* xi+ xi*(z1-zp)-chi*(z-z1)
    arg009 =  aaa/2.d0/tmp**2* xi+ xi*(z1-zp)- xi*(z-z1)
    arg010 =  aaa/2.d0/tmp**2* xi+chi*(z1-zp)- xi*(z-z1)
    arg011 =  aaa/2.d0/tmp**2*chi+chi*(z1-zp)- xi*(z-z1)
    arg012 =  aaa/2.d0/tmp**2*chi+ xi*(z1-zp)-chi*(z-z1)
    arg104 =   gp/2.d0/tmp-tmp*(z1-zp)
    arg105 =   gp/2.d0/tmp+tmp*( z-z1)
    arg107 =   xi/2.d0/tmp+tmp*( z-zp)
    arg108 =   xi/2.d0/tmp-tmp*( z-zp)
    arg109 =   xi/2.d0/tmp+tmp*(z1-zp)
    arg110 =   xi/2.d0/tmp-tmp*( z-z1)
    arg111 =  chi/2.d0/tmp+tmp*( z-zp)
    arg112 =  chi/2.d0/tmp-tmp*( z-zp)
    arg113 =  chi/2.d0/tmp+tmp*(z1-zp)
    arg114 =  chi/2.d0/tmp-tmp*( z-z1)
    !
    t1 = -EXP(arg000)*beta/alpha
    t2 =  exp_erfc(arg011,arg113)-exp_erfc(arg011,arg111) &
         +exp_erfc(arg008,arg107)-beta/alpha*exp_erfc(arg009,arg109) &
         +exp_erfc(arg005,arg104)*2.d0*lambda/alpha
    t3 =  exp_erfc(arg012,arg114)-exp_erfc(arg012,arg112) &
         +exp_erfc(arg010,arg108)-beta/alpha*exp_erfc(arg009,arg110) &
         +exp_erfc(arg004,arg105)*2.d0*lambda/alpha
    !
    vl22j1 = gp**2*(t1-(t2+t3)/4.d0)*dbesj1(gp*rxy)/lambda
    !
    RETURN
    !
  END FUNCTION vl22j1
  !
  !
  !-----------------------------------------------------------------------
  FUNCTION vl11( gp, aaa, tmp, z1, z, zp, rxy )
    !---------------------------------------------------------------------
    !
    USE kinds,  ONLY : DP
    !
    IMPLICIT NONE
    !
    REAL(DP) :: vl11
    REAL(DP), INTENT(IN) :: gp, aaa, tmp, z1, z, zp, rxy
    !
    ! ... local variables
    !
    REAL(DP) :: alpha, kappa, xi, t1, t2, t3, arg001, arg002,   &
                arg003, arg006, arg101, arg102, arg103, arg104, &
                arg109
    !
    alpha = aaa+gp+SQRT(aaa**2+gp**2)
    kappa = aaa-gp+SQRT(aaa**2+gp**2)
    xi    = aaa   +SQRT(aaa**2+gp**2)
    !
    arg001 =  gp*(z-zp)
    arg002 = -gp*(z-zp)
    arg003 =  gp*(z+zp-2.d0*z1)
    arg006 =  aaa/2.d0/tmp**2*xi+gp*(z-z1)+xi*(z1-zp)
    arg101 =   gp/2.d0/tmp+tmp*( z-zp)
    arg102 =   gp/2.d0/tmp-tmp*( z-zp)
    arg103 =   gp/2.d0/tmp+tmp*(z1-zp)
    arg104 =   gp/2.d0/tmp-tmp*(z1-zp)
    arg109 =   xi/2.d0/tmp+tmp*(z1-zp)
    !
    t1 = -EXP(arg003)*kappa/alpha
    t2 =  exp_erfc(arg001,arg101)-exp_erfc(arg001,arg103) &
         +exp_erfc(arg002,arg102)-kappa/alpha*exp_erfc(arg003,arg104)
    t3 =  exp_erfc(arg006,arg109)/alpha
    !
    vl11 = t1/2.d0-(t2/4.d0+gp*t3/2.d0)
    !
    RETURN
    !
  END FUNCTION vl11
  !
  !
  !---------------------------------------------------------------------
  FUNCTION vl22( gp, aaa, tmp, z1, z, zp, rxy )
    !------------------------------------------------------------------
    !
    USE kinds,  ONLY : DP
    !
    IMPLICIT NONE
    !
    REAL(DP) :: vl22
    REAL(DP), INTENT(IN) :: gp, aaa, tmp, z1, z, zp, rxy
    !
    ! ... local variables
    !
    REAL(DP) :: alpha, beta, xi, chi, lambda, t1, t2, t3, arg000, &
                arg005, arg008, arg009, arg011, arg104, arg107, &
                arg109, arg111, arg113
    !
    alpha = aaa+gp+SQRT(aaa**2+gp**2)
    beta  = aaa+gp-SQRT(aaa**2+gp**2)
    xi    = aaa   +SQRT(aaa**2+gp**2)
    chi   = aaa   -SQRT(aaa**2+gp**2)
    lambda=        SQRT(aaa**2+gp**2)
    !
    arg000 =-xi*(z+zp-2.d0*z1)
    arg005 =-gp*(z1-zp)-xi*(z-z1)
    arg008 = aaa/2.d0/tmp**2*xi + xi*(z1-zp)-chi*(z-z1)
    arg009 = aaa/2.d0/tmp**2*xi + xi*(z1-zp)- xi*(z-z1)
    arg011 = aaa/2.d0/tmp**2*chi+chi*(z1-zp)- xi*(z-z1)
    arg104 =  gp/2.d0/tmp-tmp*(z1-zp)
    arg107 =  xi/2.d0/tmp+tmp*( z-zp)
    arg109 =  xi/2.d0/tmp+tmp*(z1-zp)
    arg111 = chi/2.d0/tmp+tmp*( z-zp)
    arg113 = chi/2.d0/tmp+tmp*(z1-zp)
    !
    t1 = -EXP(arg000)*beta/alpha
    t2 =  exp_erfc(arg011,arg113)-exp_erfc(arg011,arg111) &
         +exp_erfc(arg008,arg107)-beta/alpha*exp_erfc(arg009,arg109)
    t3 =  exp_erfc(arg005,arg104)/alpha
    !
    vl22 = gp*t1/2.d0/lambda-gp*(t2/4.d0/lambda+t3/2.d0)
    !
    RETURN
    !
  END FUNCTION vl22
  !
  !
  !--------------------------------------------------------------------------
  FUNCTION dvl11( gp, aaa, tmp, z1, z, zp, rxy )
    !------------------------------------------------------------------------
    !
    USE kinds,  ONLY : DP
    !
    IMPLICIT NONE
    !
    REAL(DP) :: dvl11
    REAL(DP), INTENT(IN) :: gp, aaa, tmp, z1, z, zp, rxy
    !
    ! ... local variables
    !
    REAL(DP) :: alpha, kappa, xi, t1, t2, t3, arg001, arg002,   &
                arg003, arg006, arg007, arg101, arg102, arg103, &
                arg104, arg105, arg106, arg109, arg110
    !
    alpha = aaa+gp+SQRT(aaa**2+gp**2)
    kappa = aaa-gp+SQRT(aaa**2+gp**2)
    xi    = aaa   +SQRT(aaa**2+gp**2)
    !
    arg001 =  gp*(z-zp)
    arg002 = -gp*(z-zp)
    arg003 =  gp*(z+zp-2.d0*z1)
    arg006 =  aaa/2.d0/tmp**2*xi+gp*( z-z1)+xi*(z1-zp)
    arg007 =  aaa/2.d0/tmp**2*xi-gp*(z1-zp)-xi*( z-z1)
    arg101 =   gp/2.d0/tmp+tmp*( z-zp)
    arg102 =   gp/2.d0/tmp-tmp*( z-zp)
    arg103 =   gp/2.d0/tmp+tmp*(z1-zp)
    arg104 =   gp/2.d0/tmp-tmp*(z1-zp)
    arg105 =   gp/2.d0/tmp+tmp*( z-z1)
    arg106 =   gp/2.d0/tmp-tmp*( z-z1)
    arg109 =   xi/2.d0/tmp+tmp*(z1-zp)
    arg110 =   xi/2.d0/tmp-tmp*( z-z1)
    !
    t1 = -EXP(arg003)*kappa/alpha
    t2 =  exp_erfc(arg001,arg103)-exp_erfc(arg001,arg101) &
         +exp_erfc(arg002,arg102)-exp_erfc(arg003,arg104)*kappa/alpha &
         -exp_erfc(arg006,arg109)*xi/alpha*2.d0
    t3 =  exp_erfc(arg002,arg102)-exp_erfc(arg002,arg106) &
         -exp_erfc(arg001,arg101)-exp_erfc(arg003,arg105)*kappa/alpha &
         +exp_erfc(arg007,arg110)*gp/alpha*2.d0
    !
    dvl11 = gp*(t1-(t2+t3)/4.d0)
    !
    RETURN
    !
  END FUNCTION dvl11
  !
  !
  !--------------------------------------------------------------------------
  FUNCTION dvl22( gp, aaa, tmp, z1, z, zp, rxy )
    !------------------------------------------------------------------------
    !
    USE kinds, ONLY : DP
    !
    IMPLICIT NONE
    !
    REAL(DP) :: dvl22
    REAL(DP), INTENT(IN) :: gp, aaa, tmp, z1, z, zp, rxy
    !
    ! ... local variables
    !
    REAL(DP) :: alpha, beta, kappa, xi, chi, lambda, arg000,    &
                arg004, arg005, arg008, arg009, arg010, arg011, &
                arg012, arg104, arg105, arg107, arg108, arg109, &
                arg110, arg111, arg112, arg113, arg114, t1, t2, t3
    !
    alpha = aaa+gp+SQRT(aaa**2+gp**2)
    beta  = aaa+gp-SQRT(aaa**2+gp**2)
    kappa = aaa-gp+SQRT(aaa**2+gp**2)
    xi    = aaa   +SQRT(aaa**2+gp**2)
    chi   = aaa   -SQRT(aaa**2+gp**2)
    lambda=        SQRT(aaa**2+gp**2)
    !
    arg000 = -xi*(z+zp-2.d0*z1)
    arg004 =  gp*( z-z1)+xi*(z1-zp)
    arg005 = -gp*(z1-zp)-xi*( z-z1)
    arg008 =  aaa/2.d0/tmp**2*xi +xi* (z1-zp)-chi*(z-z1)
    arg009 =  aaa/2.d0/tmp**2*xi +xi* (z1-zp)- xi*(z-z1)
    arg010 =  aaa/2.d0/tmp**2*xi +chi*(z1-zp)- xi*(z-z1)
    arg011 =  aaa/2.d0/tmp**2*chi+chi*(z1-zp)- xi*(z-z1)
    arg012 =  aaa/2.d0/tmp**2*chi+ xi*(z1-zp)-chi*(z-z1)
    arg104 =   gp/2.d0/tmp-tmp*(z1-zp)
    arg105 =   gp/2.d0/tmp+tmp*( z-z1)
    arg107 =   xi/2.d0/tmp+tmp*( z-zp)
    arg108 =   xi/2.d0/tmp-tmp*( z-zp)
    arg109 =   xi/2.d0/tmp+tmp*(z1-zp)
    arg110 =   xi/2.d0/tmp-tmp*( z-z1)
    arg111 =  chi/2.d0/tmp+tmp*( z-zp)
    arg112 =  chi/2.d0/tmp-tmp*( z-zp)
    arg113 =  chi/2.d0/tmp+tmp*(z1-zp)
    arg114 =  chi/2.d0/tmp-tmp*( z-z1)
    !
    t1 = EXP(arg000)*beta*xi/alpha/lambda
    t2 = (exp_erfc(arg011,arg111)-exp_erfc(arg011,arg113))*chi/lambda &
         -exp_erfc(arg008,arg107)*xi/lambda  &
         +exp_erfc(arg009,arg109)*xi*beta/alpha/lambda &
         +exp_erfc(arg005,arg104)*gp/alpha*2.d0
    t3 = (exp_erfc(arg012,arg112)-exp_erfc(arg012,arg114))*xi/lambda  &
         -exp_erfc(arg010,arg108)*chi/lambda &
         +exp_erfc(arg009,arg110)*xi*beta/alpha/lambda &
         -exp_erfc(arg004,arg105)*xi/alpha*2.d0
    !
    dvl22 = gp*(t1-(t2+t3)/4.d0)
    !
    RETURN
    !
  END FUNCTION dvl22
  !
  !
  !---------------------------------------------------------------------
  FUNCTION dvl11j0( gp, aaa, tmp, z1, z, zp, rxy )
    !--------------------------------------------------------------------
    !
    USE kinds,  ONLY : DP
    !
    IMPLICIT NONE
    !
    REAL(DP) :: dvl11j0
    REAL(DP), INTENT(IN) :: gp, aaa, tmp, z1, z, zp, rxy
    !
    ! ... local variables
    !
    REAL(DP) :: alpha, kappa, xi, t1, t2, t3, arg001, arg002,   &
                arg003, arg006, arg007, arg101, arg102, arg103, &
                arg104, arg105, arg106, arg109, arg110
    !
    alpha = aaa+gp+SQRT(aaa**2+gp**2)
    kappa = aaa-gp+SQRT(aaa**2+gp**2)
    xi    = aaa   +SQRT(aaa**2+gp**2)
    !
    arg001 =  gp*(z-zp)
    arg002 = -gp*(z-zp)
    arg003 =  gp*(z+zp-2.d0*z1)
    arg006 =  aaa/2.d0/tmp**2*xi+gp*(z-z1)+xi*(z1-zp)
    arg007 =  aaa/2.d0/tmp**2*xi-gp*(z1-zp)-xi*(z-z1)
    arg101 =   gp/2.d0/tmp+tmp*(z-zp)
    arg102 =   gp/2.d0/tmp-tmp*(z-zp)
    arg103 =   gp/2.d0/tmp+tmp*(z1-zp)
    arg104 =   gp/2.d0/tmp-tmp*(z1-zp)
    arg105 =   gp/2.d0/tmp+tmp*(z-z1)
    arg106 =   gp/2.d0/tmp-tmp*(z-z1)
    arg109 =   xi/2.d0/tmp+tmp*(z1-zp)
    arg110 =   xi/2.d0/tmp-tmp*(z-z1)
    !
    t1 = -EXP(arg003)*kappa/alpha
    t2 =  exp_erfc(arg001,arg103)-exp_erfc(arg001,arg101) &
         +exp_erfc(arg002,arg102)-exp_erfc(arg003,arg104)*kappa/alpha &
         -exp_erfc(arg006,arg109)*xi/alpha*2.d0
    t3 =  exp_erfc(arg002,arg102)-exp_erfc(arg002,arg106) &
         -exp_erfc(arg001,arg101)-exp_erfc(arg003,arg105)*kappa/alpha &
         +exp_erfc(arg007,arg110)*gp/alpha*2.d0
    !
    dvl11j0 = gp*(t1-(t2+t3)/4.d0)*dbesj0(gp*rxy)
    !
    RETURN
    !
  END FUNCTION dvl11j0
  !
  !
  !------------------------------------------------------------------------
  FUNCTION dvl12j0( gp, aaa, tmp, z1, z, zp, rxy )
    !----------------------------------------------------------------------
    !
    USE kinds,  ONLY : DP
    !
    IMPLICIT NONE
    !
    REAL(DP)            :: dvl12j0
    REAL(DP), INTENT(IN) :: gp, aaa, tmp, z1, z, zp, rxy
    !
    ! ... local variables
    !
    REAL(DP) :: alpha, beta, kappa, xi, chi, lambda, arg001,    &
                arg002, arg003, arg004, arg006, arg009, arg010, &
                arg012, arg101, arg102, arg103, arg104, arg105, &
                arg108, arg109, arg110, arg112, arg114, t1, t2, t3
    !
    alpha = aaa+gp+SQRT(aaa**2+gp**2)
    beta  = aaa+gp-SQRT(aaa**2+gp**2)
    kappa = aaa-gp+SQRT(aaa**2+gp**2)
    xi    = aaa   +SQRT(aaa**2+gp**2)
    chi   = aaa   -SQRT(aaa**2+gp**2)
    lambda=        SQRT(aaa**2+gp**2)
    !
    arg001 =  gp*(z-zp)
    arg002 = -gp*(z-zp)
    arg003 =  gp*(z+zp-2.d0*z1)
    arg004 =  gp*(z-z1)+xi*(z1-zp)
    arg006 =  aaa/2.d0/tmp**2*xi + gp*(z-z1) +xi*(z1-zp)
    arg009 =  aaa/2.d0/tmp**2*xi + xi*(z1-zp)-xi*( z-z1)
    arg010 =  aaa/2.d0/tmp**2*xi +chi*(z1-zp)-xi*( z-z1)
    arg012 =  aaa/2.d0/tmp**2*chi+ xi*(z1-zp)-chi*(z-z1)
    arg101 =   gp/2.d0/tmp+tmp*( z-zp)
    arg102 =   gp/2.d0/tmp-tmp*( z-zp)
    arg103 =   gp/2.d0/tmp+tmp*(z1-zp)
    arg104 =   gp/2.d0/tmp-tmp*(z1-zp)
    arg105 =   gp/2.d0/tmp+tmp*( z-z1)
    arg108 =   xi/2.d0/tmp-tmp*( z-zp)
    arg109 =   xi/2.d0/tmp+tmp*(z1-zp)
    arg110 =   xi/2.d0/tmp-tmp*( z-z1)
    arg112 =  chi/2.d0/tmp-tmp*( z-zp)
    arg114 =  chi/2.d0/tmp-tmp*( z-z1)
    !
    t1 = -EXP(arg004)*xi/alpha
    t2 =  exp_erfc(arg001,arg103)-exp_erfc(arg001,arg101) &
         +exp_erfc(arg002,arg102)-exp_erfc(arg003,arg104)*kappa/alpha &
         -exp_erfc(arg006,arg109)*xi/alpha*2.d0
    t3 =  (exp_erfc(arg012,arg112)-exp_erfc(arg012,arg114))*xi/lambda &
         -exp_erfc(arg010,arg108)*chi/lambda &
         +exp_erfc(arg009,arg110)*xi*beta/alpha/lambda &
         -exp_erfc(arg004,arg105)*xi/alpha*2.d0
    !
    dvl12j0 = gp*(2.d0*t1-(t2+t3)/4.d0)*dbesj0(gp*rxy)
    !
    RETURN
    !
  END FUNCTION dvl12j0
  !
  !
  !--------------------------------------------------------------------------
  FUNCTION dvl21j0( gp, aaa, tmp, z1, z, zp, rxy )
    !------------------------------------------------------------------------
    !
    USE kinds,  ONLY : DP
    !
    IMPLICIT NONE
    !
    REAL(DP)            :: dvl21j0
    REAL(DP), INTENT(IN) :: gp, aaa, tmp, z1, z, zp, rxy
    !
    ! ... local variables
    !
    REAL(DP) :: alpha, beta, kappa, xi, chi, lambda, arg001,    &
                arg002, arg003, arg005, arg007, arg008, arg009, &
                arg011, arg101, arg102, arg104, arg105, arg106, &
                arg107, arg109, arg110, arg111, arg113, t1, t2, t3
    !
    alpha = aaa+gp+SQRT(aaa**2+gp**2)
    beta  = aaa+gp-SQRT(aaa**2+gp**2)
    kappa = aaa-gp+SQRT(aaa**2+gp**2)
    xi    = aaa   +SQRT(aaa**2+gp**2)
    chi   = aaa   -SQRT(aaa**2+gp**2)
    lambda=        SQRT(aaa**2+gp**2)
    !
    arg001 =  gp*(z-zp)
    arg002 = -gp*(z-zp)
    arg003 =  gp*(z+zp-2.d0*z1)
    arg005 = -gp*(z1-zp)-xi*(z-z1)
    arg007 =  aaa/2.d0/tmp**2*xi - gp*(z1-zp)- xi*(z-z1)
    arg008 =  aaa/2.d0/tmp**2*xi + xi*(z1-zp)-chi*(z-z1)
    arg009 =  aaa/2.d0/tmp**2*xi + xi*(z1-zp)- xi*(z-z1)
    arg011 =  aaa/2.d0/tmp**2*chi+chi*(z1-zp)- xi*(z-z1)
    arg101 =   gp/2.d0/tmp+tmp*( z-zp)
    arg102 =   gp/2.d0/tmp-tmp*( z-zp)
    arg104 =   gp/2.d0/tmp-tmp*(z1-zp)
    arg105 =   gp/2.d0/tmp+tmp*( z-z1)
    arg106 =   gp/2.d0/tmp-tmp*( z-z1)
    arg107 =   xi/2.d0/tmp+tmp*( z-zp)
    arg109 =   xi/2.d0/tmp+tmp*(z1-zp)
    arg110 =   xi/2.d0/tmp-tmp*( z-z1)
    arg111 =  chi/2.d0/tmp+tmp*( z-zp)
    arg113 =  chi/2.d0/tmp+tmp*(z1-zp)
    !
    t1 =  EXP(arg005)*gp/alpha
    t2 = (exp_erfc(arg011,arg111)-exp_erfc(arg011,arg113))*chi/lambda &
         -exp_erfc(arg008,arg107)*xi/lambda &
         +exp_erfc(arg009,arg109)*xi*beta/alpha/lambda    &
         +exp_erfc(arg005,arg104)*gp/alpha*2.d0
    t3 =  exp_erfc(arg002,arg102)-exp_erfc(arg002,arg106) &
         -exp_erfc(arg001,arg101)-exp_erfc(arg003,arg105)*kappa/alpha &
         +exp_erfc(arg007,arg110)*gp/alpha*2.d0
    !
    dvl21j0 = gp*(2.d0*t1-(t2+t3)/4.d0)*dbesj0(gp*rxy)
    !
    RETURN
    !
  END FUNCTION dvl21j0
  !
  !
  !-------------------------------------------------------------------------
  FUNCTION dvl22j0( gp, aaa, tmp, z1, z, zp, rxy )
    !------------------------------------------------------------------------
    !
    USE kinds,  ONLY : DP
    !
    IMPLICIT NONE
    !
    REAL(DP) :: dvl22j0
    REAL(DP), INTENT(IN) :: gp, aaa, tmp, z1, z, zp, rxy
    !
    ! ... local variables
    !
    REAL(DP) :: alpha, beta, kappa, xi, chi, lambda, arg000,    &
                arg004, arg005, arg008, arg009, arg010, arg011, &
                arg012, arg104, arg105, arg107, arg108, arg109, &
                arg110, arg111, arg112, arg113, arg114, t1, t2, t3
    !
    alpha = aaa+gp+SQRT(aaa**2+gp**2)
    beta  = aaa+gp-SQRT(aaa**2+gp**2)
    kappa = aaa-gp+SQRT(aaa**2+gp**2)
    xi    = aaa   +SQRT(aaa**2+gp**2)
    chi   = aaa   -SQRT(aaa**2+gp**2)
    lambda=        SQRT(aaa**2+gp**2)
    !
    arg000 = -xi*(z+zp-2.d0*z1)
    arg004 =  gp*( z-z1)+xi*(z1-zp)
    arg005 = -gp*(z1-zp)-xi*( z-z1)
    arg008 =  aaa/2.d0/tmp**2*xi + xi*(z1-zp)-chi*(z-z1)
    arg009 =  aaa/2.d0/tmp**2*xi + xi*(z1-zp)- xi*(z-z1)
    arg010 =  aaa/2.d0/tmp**2*xi +chi*(z1-zp)- xi*(z-z1)
    arg011 =  aaa/2.d0/tmp**2*chi+chi*(z1-zp)- xi*(z-z1)
    arg012 =  aaa/2.d0/tmp**2*chi+ xi*(z1-zp)-chi*(z-z1)
    arg104 =   gp/2.d0/tmp-tmp*(z1-zp)
    arg105 =   gp/2.d0/tmp+tmp*( z-z1)
    arg107 =   xi/2.d0/tmp+tmp*( z-zp)
    arg108 =   xi/2.d0/tmp-tmp*( z-zp)
    arg109 =   xi/2.d0/tmp+tmp*(z1-zp)
    arg110 =   xi/2.d0/tmp-tmp*( z-z1)
    arg111 =  chi/2.d0/tmp+tmp*( z-zp)
    arg112 =  chi/2.d0/tmp-tmp*( z-zp)
    arg113 =  chi/2.d0/tmp+tmp*(z1-zp)
    arg114 =  chi/2.d0/tmp-tmp*( z-z1)
    !
    t1 =  EXP(arg000)*beta*xi/alpha/lambda
    t2 = (exp_erfc(arg011,arg111)-exp_erfc(arg011,arg113))*chi/lambda &
         -exp_erfc(arg008,arg107)*xi/lambda  &
         +exp_erfc(arg009,arg109)*xi*beta/alpha/lambda &
         +exp_erfc(arg005,arg104)*gp/alpha*2.d0
    t3 = (exp_erfc(arg012,arg112)-exp_erfc(arg012,arg114))*xi/lambda  &
         -exp_erfc(arg010,arg108)*chi/lambda &
         +exp_erfc(arg009,arg110)*xi*beta/alpha/lambda &
         -exp_erfc(arg004,arg105)*xi/alpha*2.d0
    !
    dvl22j0 = gp*(t1-(t2+t3)/4.d0)*dbesj0(gp*rxy)
    !
    RETURN
    !
  END FUNCTION dvl22j0
  !
  !
  !-----------------------------------------------------------------------
  SUBROUTINE qromb( func, aaa, tmp, z1, z, zp, rxy, b, ss )
    !---------------------------------------------------------------------
    !
    USE kinds,  ONLY : DP
    !
    REAL(DP), INTENT(IN) :: aaa, tmp, z1, z, zp, rxy, b
    REAL(DP), INTENT(OUT) :: ss
    ! ss=int_a^b func(gp,aaa,tmp,z1,z,zp,rxy) dgp
    !
    ! ... local variables
    !
    INTEGER, PARAMETER :: jmax=20, jmaxp=jmax+1, k=5, km=k-1
    INTEGER :: j
    REAL(DP), PARAMETER :: a=0.0_DP, eps=1.e-12
    REAL(DP) :: dss, h(jmaxp), s(jmaxp)
    REAL(DP), EXTERNAL :: func
    !
    h(1) = 1.0_DP
    !
    DO j = 1, jmax
      !
      CALL trapzd( func, aaa, tmp, z1, z, zp, rxy, a, b, s(j), j )
      !
      IF (j >= k) THEN
        CALL polint(h(j-km),s(j-km),k,0.0_DP,ss,dss)
        IF (ABS(ss) <=1.e-8 ) RETURN
        IF (ABS(dss)<=eps*ABS(ss)) RETURN
      ENDIF
      !
      s(j+1) = s(j)
      h(j+1) = 0.25*h(j)
      !
    ENDDO
    !
    STOP 'too many steps in qromb'
    !
  END SUBROUTINE qromb
  !
  !
  !---------------------------------------------------------------------
  SUBROUTINE trapzd( func, aaa, tmp, z1, z, zp, rxy, a, b, s, n )
    !-------------------------------------------------------------------
    !
    USE kinds,  ONLY : DP
    !
    INTEGER, INTENT(IN) :: n
    REAL(DP), INTENT(IN) :: aaa, tmp, z1, z, zp, rxy, a, b
    REAL(DP), INTENT(INOUT) :: s
    REAL(DP), EXTERNAL :: func
    !
    ! ... local variables
    !
    INTEGER :: it,j
    REAL(DP) :: del, sum, tnm, x
    !
    IF (n == 1) THEN
      s = 0.5*(b-a)*(func(a,aaa,tmp,z1,z,zp,rxy)+func(b,aaa,tmp,z1,z,zp,rxy))
    ELSE
      it = 2**(n-2)
      tnm = it
      del = (b-a)/tnm
      x = a+0.5*del
      sum = 0.
      !
      DO j = 1, it
        sum = sum+func(x,aaa,tmp,z1,z,zp,rxy)
        x = x+del
      ENDDO
      !
      s = 0.5*(s+(b-a)*sum/tnm)
    ENDIF
    !
    RETURN
    !
  END SUBROUTINE trapzd
  !
  !
  !--------------------------------------------------------------------
  SUBROUTINE polint( xa, ya, n, x, y, dy )
    !-------------------------------------------------------------------
    !
    USE kinds,  ONLY : DP
    !
    INTEGER, INTENT(IN) :: n
    REAL(DP), INTENT(IN) :: x, xa(n), ya(n)
    REAL(DP), INTENT(OUT) :: y, dy
    !
    ! ... local variables
    !
    INTEGER :: i, m, ns
    INTEGER, PARAMETER :: nmax=10
    REAL(DP) :: den, dif, dift, ho, hp, w, c(nmax), d(nmax)
    !
    ns = 1
    dif = ABS(x-xa(1))
    !
    DO i = 1, n
       dift = ABS(x-xa(i))
       IF (dift < dif) THEN
         ns = i
         dif = dift
       ENDIF
       c(i) = ya(i)
       d(i) = ya(i)
    ENDDO
    !
    y = ya(ns)
    ns = ns-1
    !
    DO m = 1, n-1
       !
       DO i = 1, n-m
          ho = xa(i)-x
          hp = xa(i+m)-x
          w = c(i+1)-d(i)
          den = ho-hp
          IF (den == 0.) STOP 'failure in polint'
          den = w/den
          d(i) = hp*den
          c(i) = ho*den
       ENDDO
       !
       IF (2*ns < n-m) THEN
         dy = c(ns+1)
       ELSE
         dy = d(ns)
         ns = ns-1
       ENDIF
       !
       y = y+dy
       !
    ENDDO
    !
    RETURN
    !
  END SUBROUTINE polint
  !
  !
  !-------------------------------------------------------------------------
  REAL(8) FUNCTION dbesj0( x )
   !------------------------------------------------------------------------
   !! Bessel J_0(x) function in double precision.
   !
   IMPLICIT NONE
   !
   REAL(8), INTENT(IN) :: x
   REAL(8), PARAMETER  :: pi4 = 0.78539816339744830962d0
   REAL(8), PARAMETER  :: a(0:7) = (/ &
      -0.0000000000023655394d0, 0.0000000004708898680d0, &
      -0.0000000678167892231d0, 0.0000067816840038636d0, &
      -0.0004340277777716935d0, 0.0156249999999992397d0, &
      -0.2499999999999999638d0, 0.9999999999999999997d0 /)
   REAL(8), PARAMETER  :: b(0:12, 0:4) = RESHAPE( (/ &
       0.0000000000626681117d0, -0.0000000022270614428d0, & 
       0.0000000662981656302d0, -0.0000016268486502196d0, & 
       0.0000321978384111685d0, -0.0005005237733315830d0, & 
       0.0059060313537449816d0, -0.0505265323740109701d0, & 
       0.2936432097610503985d0, -1.0482565081091638637d0, & 
       1.9181123286040428113d0, -1.1319199475221700100d0, & 
      -0.1965480952704682000d0, &
       0.0000000000457457332d0, -0.0000000015814772025d0, & 
       0.0000000455487446311d0, -0.0000010735201286233d0, & 
       0.0000202015179970014d0, -0.0002942392368203808d0, & 
       0.0031801987726150648d0, -0.0239875209742846362d0, & 
       0.1141447698973777641d0, -0.2766726722823530233d0, & 
       0.1088620480970941648d0,  0.5136514645381999197d0, & 
      -0.2100594022073706033d0, &
       0.0000000000331366618d0, -0.0000000011119090229d0, & 
       0.0000000308823040363d0, -0.0000006956602653104d0, & 
       0.0000123499947481762d0, -0.0001662951945396180d0, & 
       0.0016048663165678412d0, -0.0100785479932760966d0, & 
       0.0328996815223415274d0, -0.0056168761733860688d0, & 
      -0.2341096400274429386d0,  0.2551729256776404262d0, & 
       0.2288438186148935667d0, & 
       0.0000000000238007203d0, -0.0000000007731046439d0, & 
       0.0000000206237001152d0, -0.0000004412291442285d0, & 
       0.0000073107766249655d0, -0.0000891749801028666d0, & 
       0.0007341654513841350d0, -0.0033303085445352071d0, & 
       0.0015425853045205717d0,  0.0521100583113136379d0, & 
      -0.1334447768979217815d0, -0.1401330292364750968d0, & 
       0.2685616168804818919d0, &
       0.0000000000169355950d0, -0.0000000005308092192d0, & 
       0.0000000135323005576d0, -0.0000002726650587978d0, & 
       0.0000041513240141760d0, -0.0000443353052220157d0, & 
       0.0002815740758993879d0, -0.0004393235121629007d0, & 
      -0.0067573531105799347d0,  0.0369141914660130814d0, & 
       0.0081673361942996237d0, -0.2573381285898881860d0, & 
       0.0459580257102978932d0 /), (/13, 5/) )
   REAL(8), PARAMETER  :: c(0:13, 0:4) = RESHAPE( (/ &
      -0.00000000003009451757d0, -0.00000000014958003844d0, & 
       0.00000000506854544776d0,  0.00000001863564222012d0, & 
      -0.00000060304249068078d0, -0.00000147686259937403d0, & 
       0.00004714331342682714d0,  0.00006286305481740818d0, & 
      -0.00214137170594124344d0, -0.00089157336676889788d0, & 
       0.04508258728666024989d0, -0.00490362805828762224d0, & 
      -0.27312196367405374426d0,  0.04193925184293450356d0,  &
      -0.00000000000712453560d0, -0.00000000041170814825d0, & 
       0.00000000138012624364d0,  0.00000005704447670683d0, & 
      -0.00000019026363528842d0, -0.00000533925032409729d0, & 
       0.00001736064885538091d0,  0.00030692619152608375d0, & 
      -0.00092598938200644367d0, -0.00917934265960017663d0, & 
       0.02287952522866389076d0,  0.10545197546252853195d0, & 
      -0.16126443075752985095d0, -0.19392874768742235538d0,  &
       0.00000000002128344556d0, -0.00000000031053910272d0, & 
      -0.00000000334979293158d0,  0.00000004507232895050d0, & 
       0.00000036437959146427d0, -0.00000446421436266678d0, & 
      -0.00002523429344576552d0,  0.00027519882931758163d0, & 
       0.00097185076358599358d0, -0.00898326746345390692d0, & 
      -0.01665959196063987584d0,  0.11456933464891967814d0, & 
       0.07885001422733148815d0, -0.23664819446234712621d0,  & 
       0.00000000003035295055d0,  0.00000000005486066835d0, & 
      -0.00000000501026824811d0, -0.00000000501246847860d0, & 
       0.00000058012340163034d0,  0.00000016788922416169d0, & 
      -0.00004373270270147275d0,  0.00001183898532719802d0, & 
       0.00189863342862291449d0, -0.00113759249561636130d0, & 
      -0.03846797195329871681d0,  0.02389746880951420335d0, & 
       0.22837862066532347461d0, -0.06765394811166522844d0,  &
       0.00000000001279875977d0,  0.00000000035925958103d0, & 
      -0.00000000228037105967d0, -0.00000004852770517176d0, & 
       0.00000028696428000189d0,  0.00000440131125178642d0, & 
      -0.00002366617753349105d0, -0.00024412456252884129d0, & 
       0.00113028178539430542d0,  0.00708470513919789080d0, & 
      -0.02526914792327618386d0, -0.08006137953480093426d0, & 
       0.16548380461475971846d0,  0.14688405470042110229d0/), (/14, 5/) )
   REAL(8), PARAMETER  :: d(0:12, 0:3) = RESHAPE( (/ &
       1.059601355592185731d-14, -2.71150591218550377d-13, & 
       8.6514809056201638d-12,   -4.6264028554286627d-10, & 
       5.0815403835647104d-8,    -1.76722552048141208d-5, & 
       0.16286750396763997378d0,  2.949651820598278873d-13, & 
      -8.818215611676125741d-12,  3.571119876162253451d-10, & 
      -2.631924120993717060d-8,   4.709502795656698909d-6, & 
      -5.208333333333283282d-3, & 
       7.18344107717531977d-15,  -2.51623725588410308d-13, & 
       8.6017784918920604d-12,   -4.6256876614290359d-10, & 
       5.0815343220437937d-8,    -1.76722551764941970d-5, & 
       0.16286750396763433767d0,  2.2327570859680094777d-13, & 
      -8.464594853517051292d-12,  3.563766464349055183d-10, & 
      -2.631843986737892965d-8,   4.709502342288659410d-6, & 
      -5.2083333332278466225d-3, & 
       5.15413392842889366d-15,  -2.27740238380640162d-13, & 
       8.4827767197609014d-12,   -4.6224753682737618d-10, & 
       5.0814848128929134d-8,    -1.76722547638767480d-5, & 
       0.16286750396748926663d0,  1.7316195320192170887d-13, & 
      -7.971122772293919646d-12,  3.544039469911895749d-10, & 
      -2.631443902081701081d-8,   4.709498228695400603d-6, & 
      -5.2083333315143653610d-3, & 
       3.84653681453798517d-15,  -2.04464520778789011d-13, & 
       8.3089298605177838d-12,   -4.6155016158412096d-10, & 
       5.0813263696466650d-8,    -1.76722528311426167d-5, & 
       0.16286750396650065930d0,  1.3797879972460878797d-13, & 
      -7.448089381011684812d-12,  3.512733797106959780d-10, & 
      -2.630500895563592722d-8,   4.709483934775839193d-6, & 
      -5.2083333227940760113d-3 /), (/13, 4/) )
   REAL(8) :: w, t, y, v, theta
   INTEGER :: k, i
   !
   w = ABS(x)
   IF ( w < 1.0d0 ) THEN
      t = w * w
      y = a(0)
      DO i = 1, 7
         y = y * t + a(i)
      ENDDO
   ELSEIF( w < 8.5d0 ) THEN
      t = w * w * 0.0625d0
      k = INT(t)
      t = t - (k + 0.5d0)
      y = b(0,k)
      DO i = 1, 12
         y = y * t + b(i,k)
      ENDDO
   ELSEIF( w < 12.5d0 ) THEN
      k = INT(w)
      t = w - ( k + 0.5d0 )
      k = k - 8
      y = c(0,k)
      DO i = 1, 13
         y = y * t + c(i,k)
      ENDDO
   ELSE
      v = 24.0d0 / w
      t = v * v
      k = INT(t)
      y = d(0,k)
      DO i = 1, 6
         y = y * t + d(i,k)
      ENDDO
      y = y * SQRT(v)
      theta = d(7,k)
      DO i = 8, 12
         theta = theta * t + d(i,k)
      ENDDO
      theta = theta * v - pi4
      y = y * COS( w + theta )
   ENDIF
   !
   dbesj0 = y
   !
  END FUNCTION dbesj0
  !
  !
  !--------------------------------------------------------------------------
  REAL(8) FUNCTION dbesj1( x )
   !-------------------------------------------------------------------------
   !! Bessel J_1(x) function in double precision.
   !
   IMPLICIT NONE
   !
   REAL(8), INTENT(IN) :: x
   REAL(8), PARAMETER :: pi4 = 0.78539816339744830962d0
   REAL(8), PARAMETER :: a(0:7) = (/ &
      -0.00000000000014810349d0,  0.00000000003363594618d0, &
      -0.00000000565140051697d0,  0.00000067816840144764d0, &
      -0.00005425347222188379d0,  0.00260416666666662438d0, &
      -0.06249999999999999799d0,  0.49999999999999999998d0 /)
   REAL(8), PARAMETER :: b(0:12, 0:4) = RESHAPE( (/ &
       0.00000000000243721316d0, -0.00000000009400554763d0, &
       0.00000000306053389980d0, -0.00000008287270492518d0, &
       0.00000183020515991344d0, -0.00003219783841164382d0, &
       0.00043795830161515318d0, -0.00442952351530868999d0, &
       0.03157908273375945955d0, -0.14682160488052520107d0, &
       0.39309619054093640008d0, -0.47952808215101070280d0, &
       0.14148999344027125140d0, &
       0.00000000000182119257d0, -0.00000000006862117678d0, &
       0.00000000217327908360d0, -0.00000005693592917820d0, &
       0.00000120771046483277d0, -0.00002020151799736374d0, &
       0.00025745933218048448d0, -0.00238514907946126334d0, &
       0.01499220060892984289d0, -0.05707238494868888345d0, &
       0.10375225210588234727d0, -0.02721551202427354117d0, &
      -0.06420643306727498985d0, &
       0.000000000001352611196d0, -0.000000000049706947875d0, &
       0.000000001527944986332d0, -0.000000038602878823401d0, &
       0.000000782618036237845d0, -0.000012349994748451100d0, &
       0.000145508295194426686d0, -0.001203649737425854162d0, &
       0.006299092495799005109d0, -0.016449840761170764763d0, &
       0.002106328565019748701d0,  0.058527410006860734650d0, &
      -0.031896615709705053191d0, &
       0.000000000000997982124d0, -0.000000000035702556073d0, &
       0.000000001062332772617d0, -0.000000025779624221725d0, &
       0.000000496382962683556d0, -0.000007310776625173004d0, &
       0.000078028107569541842d0, -0.000550624088538081113d0, &
       0.002081442840335570371d0, -0.000771292652260286633d0, &
      -0.019541271866742634199d0,  0.033361194224480445382d0, &
       0.017516628654559387164d0, &
       0.000000000000731050661d0, -0.000000000025404499912d0, &
       0.000000000729360079088d0, -0.000000016915375004937d0, &
       0.000000306748319652546d0, -0.000004151324014331739d0, &
       0.000038793392054271497d0, -0.000211180556924525773d0, &
       0.000274577195102593786d0,  0.003378676555289966782d0, &
      -0.013842821799754920148d0, -0.002041834048574905921d0, &
       0.032167266073736023299d0 /), (/13, 5/) )
   REAL(8), PARAMETER :: c(0:13, 0:4) = RESHAPE( (/ &
       -0.00000000001185964494d0,  0.00000000039110295657d0, &
        0.00000000180385519493d0, -0.00000005575391345723d0, &
       -0.00000018635897017174d0,  0.00000542738239401869d0, &
        0.00001181490114244279d0, -0.00033000319398521070d0, &
       -0.00037717832892725053d0,  0.01070685852970608288d0, &
        0.00356629346707622489d0, -0.13524776185998074716d0, &
        0.00980725611657523952d0,  0.27312196367405374425d0,  &
       -0.00000000003029591097d0,  0.00000000009259293559d0, &
        0.00000000496321971223d0, -0.00000001518137078639d0, &
       -0.00000057045127595547d0,  0.00000171237271302072d0, &
        0.00004271400348035384d0, -0.00012152454198713258d0, &
       -0.00184155714921474963d0,  0.00462994691003219055d0, &
        0.03671737063840232452d0, -0.06863857568599167175d0, &
       -0.21090395092505707655d0,  0.16126443075752985095d0,  & 
       -0.00000000002197602080d0, -0.00000000027659100729d0, &
        0.00000000374295124827d0,  0.00000003684765777023d0, &
       -0.00000045072801091574d0, -0.00000327941630669276d0, &
        0.00003571371554516300d0,  0.00017664005411843533d0, &
       -0.00165119297594774104d0, -0.00485925381792986774d0, &
        0.03593306985381680131d0,  0.04997877588191962563d0, &
       -0.22913866929783936544d0, -0.07885001422733148814d0,  &
        0.00000000000516292316d0, -0.00000000039445956763d0, &
       -0.00000000066220021263d0,  0.00000005511286218639d0, &
        0.00000005012579400780d0, -0.00000522111059203425d0, &
       -0.00000134311394455105d0,  0.00030612891890766805d0, &
       -0.00007103391195326182d0, -0.00949316714311443491d0, &
        0.00455036998246516948d0,  0.11540391585989614784d0, &
       -0.04779493761902840455d0, -0.22837862066532347460d0,  &
        0.00000000002697817493d0, -0.00000000016633326949d0, &
       -0.00000000433134860350d0,  0.00000002508404686362d0, &
        0.00000048528284780984d0, -0.00000258267851112118d0, &
       -0.00003521049080466759d0,  0.00016566324273339952d0, &
        0.00146474737522491617d0, -0.00565140892697147306d0, &
       -0.02833882055679300400d0,  0.07580744376982855057d0, &
        0.16012275906960187978d0, -0.16548380461475971845d0 /), (/14, 5/) )
   REAL(8), PARAMETER :: d(0:12, 0:3) = RESHAPE( (/ &
       -1.272346002224188092d-14,  3.370464692346669075d-13, &
       -1.144940314335484869d-11,  6.863141561083429745d-10, &
       -9.491933932960924159d-8,   5.301676561445687562d-5,  &
        0.1628675039676399740d0,  -3.652982212914147794d-13, &
        1.151126750560028914d-11, -5.165585095674343486d-10, &
        4.657991250060549892d-8,  -1.186794704692706504d-5,  &
        1.562499999999994026d-2, &
       -8.713069680903981555d-15,  3.140780373478474935d-13, &
       -1.139089186076256597d-11,  6.862299023338785566d-10, &
       -9.491926788274594674d-8,   5.301676558106268323d-5,  &
        0.1628675039676466220d0,  -2.792555727162752006d-13, &
        1.108650207651756807d-11, -5.156745588549830981d-10, &
        4.657894859077370979d-8,  -1.186794650130550256d-5,  &
        1.562499999987299901d-2, & 
       -6.304859171204770696d-15,  2.857249044208791652d-13, &
       -1.124956921556753188d-11,  6.858482894906716661d-10, &
       -9.491867953516898460d-8,   5.301676509057781574d-5,  &
        0.1628675039678191167d0,  -2.185193490132496053d-13, &
        1.048820673697426074d-11, -5.132819367467680132d-10, &
        4.657409437372994220d-8,  -1.186794150862988921d-5,  &
        1.562499999779270706d-2, &
       -4.740417209792009850d-15,  2.578715253644144182d-13, &
       -1.104148898414138857d-11,  6.850134201626289183d-10, &
       -9.491678234174919640d-8,   5.301676277588728159d-5,  &
        0.1628675039690033136d0,  -1.755122057493842290d-13, &
        9.848723331445182397d-12, -5.094535425482245697d-10, &
        4.656255982268609304d-8,  -1.186792402114394891d-5,  &
        1.562499998712198636d-2 /), (/13, 4/) )
   REAL(8) :: w, t, y, v, theta
   INTEGER :: k, i
   !
   w = ABS(x)
   IF ( w < 1.0d0 ) THEN
      t = w * w
      y = a(0)
      DO i = 1, 7
         y = y * t + a(i)
      ENDDO
      y = y * w
   ELSEIF( w < 8.5d0 ) THEN
      t = w * w * 0.0625d0
      k = INT(t)
      t = t - (k + 0.5d0)
      y = b(0,k)
      DO i = 1, 12
         y = y * t + b(i,k)
      ENDDO
      y = y * w
   ELSEIF( w < 12.5d0 ) THEN
      k = INT(w)
      t = w - (k + 0.5d0)
      k = k - 8
      y = c(0,k)
      DO i = 1, 13
         y = y * t + c(i,k)
      ENDDO
   ELSE
      v = 24.0d0 / w
      t = v * v
      k = INT(t)
      y = d(0,k)
      DO i = 1, 6
         y = y * t + d(i,k)
      ENDDO
      y = y * SQRT(v)
      theta = d(7,k)
      DO i = 8, 12
         theta = theta * t + d(i,k)
      ENDDO
      theta = theta * v - pi4
      y = y * SIN(w + theta)
   ENDIF
   !
   IF ( x < 0.0d0 ) y = -y
   !
   dbesj1 = y
   !
  END FUNCTION dbesj1
  !
  !
  !---------------------------------------------------------------------
  REAL(8) FUNCTION exp_erfc( x, y )
    !--------------------------------------------------------------------
    !! EXP(x) * erfc(y): this function is to avoid INFINITY * ZERO for 
    !! large positive x and y.
    !
    IMPLICIT NONE
    !
    REAL(8), INTENT(IN) :: x, y
    REAL(8) :: ym, ym2, nume, deno
    REAL(8), PARAMETER :: rtpim = 0.564189583547756279d0 ! 1/SQRT(PI)
    REAL(8), PARAMETER :: r(0:4) = (/ &
       -2.99610707703542174d-3, -4.94730910623250734d-2, &
       -2.26956593539686930d-1, -2.78661308609647788d-1, &
       -2.23192459734184686d-2 /)
    REAL(8), PARAMETER :: s(0:4) = (/ &
       1.06209230528467918d-2, 1.91308926107829841d-1, &
       1.05167510706793207d0,  1.98733201817135256d0,  &
       1.00000000000000000d0 /)
    !
    IF ( x < 709.0d0 .OR. y < 4.0d0 ) THEN
       exp_erfc = EXP(x) * qe_erfc(y)
    ELSE
       ym  = 1d0 / y
       ym2 = ym**2
       nume = ( ( ( r(4) * ym2 + r(3) ) * ym2 + r(2) ) * ym2 + r(1) ) * ym2 + r(0)
       deno = ( ( ( s(4) * ym2 + s(3) ) * ym2 + s(2) ) * ym2 + s(1) ) * ym2 + s(0)
       exp_erfc = EXP( - y**2 + x) * ym * ( rtpim + ym2 * nume / deno )
    ENDIF
    !   
    RETURN
    !
  END FUNCTION exp_erfc
  !
  !
END MODULE esm

