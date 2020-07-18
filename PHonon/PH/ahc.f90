!
! Copyright (C) 2019-2019 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE ahc
!----------------------------------------------------------------------------
!! In this subroutine, the matrix elements required for AHC temperature-
!! dependent electronic structure calculation are calculated and written to
!! file. Three different quantities are computed.
!!
!! 1. ahc_gkk(ib, jb, imode)
!!   = <\psi_ib(k+q)|dV/du_{q,imode}|\psi_jb(k)>
!!   (Eq.(2) of PHonon/Doc/dfpt_self_energy.pdf)
!!
!!   1 <= ib <= nbnd, ib_ahc_min <= jb <= ib_ahc_max
!!   Needed to calculate the static or dynamic Fan term.
!!
!! 2. ahc_upfan(ib,jb,imode,jmode)
!!   = <P_{c,k+q}^+ d\psi_ib(k+q)/du_{q,imode} | dV/du_{q,jmode} | \psi_jb(k)>
!!   (Eq.(4) of PHonon/Doc/dfpt_self_energy.pdf)
!!
!!   ib_ahc_min <= ib, jb <= ib_ahc_max
!!   Here, P_{c,k+q} is the orthogonalization to the nbnd lowest-lying bands
!!   in the k+q subspace. (nbnd may differ from the number of bands used in
!!   the SCF and phonon calculations.)
!!
!!   Needed to compute the "upper Fan" term which approximates the contribution
!!   of high-energy (unoccpied) bands in the Fan term.
!!   Ref: X. Gonze, P. Boulanger, and M. Cote, Ann. Phys. 523, 168 (2011)
!!
!! 3. ahc_dw(ib, jb, imode, jdir)
!!   = i * <\psi_ib(k)|[dV_{SCF}/du_{Gamma,imode}, p_jdir]|\psi_jb(k)>
!!   (Eq.(3) of PHonon/Doc/dfpt_self_energy.pdf)
!!
!!   ib_ahc_min <= ib, jb <= ib_ahc_min+ahc_nbnd-1
!!   Here, p_jdir = -i * d/dr_jdir is the momentum operator.
!!   Computed only for q = Gamma.
!!
!!   Needed to calculate the Debye-Waller term.
!!   We use the generalized acoustic sum rule for electron-phonon matrix,
!!   which gives both the diagonal the off-diagonal matrix elements of the
!!   Debye-Waller term (in the electron eigenbasis).
!!   Ref: J.-M. Lihm, and C.-H. Park, Phys. Rev. B, 101, 121102(R) (2020)
!!
!! In all cases, the imode index is in the Cartesian basis, so that only one
!! atom with index iatm is displaced along Cartesian direction idir. In this
!! case, the mode index is imode = 3 * (iatm - 1) + idir.
!!
!! ib_ahc_min = ahc_nbndskip + 1
!! ib_ahc_max = ahc_nbndskip + ahc_nbnd
!!
!! Eigenvalues of the ahc_nbnd bands should be well-separated with the nbnd-th
!! band to avoid problems in solving the Sternheimer equation.
!!
!! Not implemented (or not tested) for the following cases:
!!   - USPP
!!   - PAW
!!   - DFPT + U
!!   - magnetism (both collinear and noncollinear)
!!   - 2d Coulomb cutoff
!!
!----------------------------------------------------------------------------
  USE kinds, ONLY :  DP
  !
  IMPLICIT NONE
  !
  ! Input parameters
  !
  CHARACTER(LEN=256) :: ahc_dir
  !! Directory where the output binary files for AHC e-ph coupling are written
  INTEGER :: ahc_nbnd
  !! Number of bands for which the electron self-energy is to be computed.
  INTEGER :: ahc_nbndskip
  !! Number of bands to exclude when computing the self-energy. The
  !! self-energy is computed for ibnd from ahc_nbndskip + 1
  !! to ahc_nbndskip + ahc_nbnd.
  LOGICAL :: skip_upperfan = .FALSE.
  !! If .true., skip the calculation of upper Fan self-energy,
  !! which involves solving the Sternheimer equation.
  !
  ! Public variables
  !
  LOGICAL :: elph_ahc
  !! If .true., calculate ahc e-ph variables.
  INTEGER :: ahc_nbnd_gauge
  !! Number of bands to compute self-energy or bands degenerate with those.
  INTEGER :: ib_ahc_gauge_min
  !! Minimum band index to compute dvpsi. Needed to deal with degeneracy.
  INTEGER :: ib_ahc_gauge_max
  !! Maximum band index to compute dvpsi. Needed to deal with degeneracy.
  !
  ! Local variables
  !
  INTEGER :: ib_ahc_min
  !! Mininum band index to compute electron self-energy
  INTEGER :: ib_ahc_max
  !! Maximum band index to compute electron self-energy
  INTEGER :: iungkk
  !! Unit for ahc_gkk (dV_q matrix element) output
  INTEGER :: iunupfan
  !! Unit for ahc_upfan (upper Fan by Sternheimer) output
  INTEGER :: iundw
  !! Unit for ahc_dw (Debye-Waller) output
  INTEGER :: nbase_ik
  !! The position in the list of the first point that belong to this npool - 1
  REAL(DP) :: e_degen_thr = 1.d-4
  !! threshold for degeneracy
  COMPLEX(DP), ALLOCATABLE :: ahc_gkk(:,:,:)
  !! (nbnd, ahc_nbnd, nmodes) dV_q matrix element
  COMPLEX(DP), ALLOCATABLE :: ahc_upfan(:,:,:,:)
  !! (ahc_nbnd, ahc_nbnd, nmodes, nmodes) upper Fan self-energy by Sternheimer
  COMPLEX(DP), ALLOCATABLE :: ahc_dw(:,:,:,:)
  !! (ahc_nbnd, ahc_nbnd, nmodes, 3)
  !! [dV,p] matrix element for Debye-Waller
  COMPLEX(DP), ALLOCATABLE :: dvpsi_cart(:,:,:)
  !! (npwx*npol, ahc_nbnd, nmodes) dV/du_{q, imode} * psi_nk in
  !! the Cartesian atomic displacement.
  COMPLEX(DP), ALLOCATABLE :: psi_gauge(:,:)
  !! (ahc_nbnd_gauge, ahc_nbnd) <psi_nk(at q)|psi_mk(at q=gamma)>
  !! Used to fix gauge of psi_nk computed at different q points
  !
CONTAINS
!------------------------------------------------------------------------------
SUBROUTINE ahc_do_upperfan(ik)
!------------------------------------------------------------------------------
!!
!! Compute ahc_upperfan(ib,jb,imode,jmode)
!! = <P_{c,k+q}^+ d\psi_ib(k+q)/du_{q,imode} | dV/du_{q,jmode} | \psi_jb(k)>
!!
!! Implements Eq.(4) of PHonon/Doc/dfpt_self_energy.pdf
!!
!! To do so, compute d\psi_ib(k+q)/du{q,imode} by solving Sternheimer equation
!! (H - e_ib) dpsi_ib = - P_{c,k+q}^+ dV/du_{q,imode} psi_ib.
!!
!! Adapted from solve_linter.f90 with modifications by Jae-Mo Lihm
!!
!------------------------------------------------------------------------------
  USE kinds,            ONLY : DP
  USE io_global,        ONLY : stdout, ionode
  USE mp,               ONLY : mp_sum
  USE mp_pools,         ONLY : intra_pool_comm, me_pool, root_pool
  USE wvfct,            ONLY : npwx, nbnd, et
  USE klist,            ONLY : ngk, igk_k, xk
  USE wavefunctions,    ONLY : evc
  USE noncollin_module, ONLY : npol
  USE uspp,             ONLY : vkb
  USE buffers,          ONLY : get_buffer
  USE qpoint,           ONLY : ikks, ikqs
  USE modes,            ONLY : nmodes
  USE eqv,              ONLY : dvpsi, dpsi, evq
  USE units_lr,         ONLY : lrwfc, iuwfc
  USE control_ph,       ONLY : tr2_ph
  USE control_lr,       ONLY : lgamma
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ik
  !! k point index where dvpsi is calculated
  !
  INTEGER :: ikk, ikq, npw, npwq, nrec, ibnd, imode, jmode
  LOGICAL :: conv_root
  !! true if linear system is converged
  INTEGER :: lter
  !! Counter on iterations of Sternheimer equation
  REAL(DP) :: thresh
  !! convergence threshold for solving Sternheimer equation
  REAL(DP) :: anorm
  !! the norm of the error of Sternheimer equation
  REAL(DP) :: error_sum
  !! error for checking the solution of Sternheimer equation
  REAL(DP) , ALLOCATABLE :: h_diag (:,:)
  !! diagonal part of the Hamiltonian
  COMPLEX(DP), ALLOCATABLE :: dpsi_cart(:,:,:)
  !! Solution of the Sternheimer equation
  EXTERNAL h_psi
  !
  EXTERNAL ch_psi_all, cg_psi
  !
  WRITE(stdout, '(5x,a,I8)') 'Computing ahc_upperfan for ik = ', ik
  !
  CALL start_clock('ahc_upfan')
  !
  ! Set Sternheimer threshold, following that of solve_linter
  thresh = MIN(1.d-1 * SQRT(tr2_ph / npol), 1.d-2)
  !
  ALLOCATE(h_diag(npwx*npol, nbnd))
  ALLOCATE(dpsi_cart(npwx*npol, ahc_nbnd, nmodes))
  dpsi_cart = (0.d0, 0.d0)
  !
  ikk = ikks(ik)
  ikq = ikqs(ik)
  npw = ngk(ikk)
  npwq = ngk(ikq)
  !
  ! Read unperturbed wavefunctions psi(k) and psi(k+q) from file
  !
  IF (lgamma) THEN
    CALL get_buffer(evc, lrwfc, iuwfc, ikk)
  ELSE
    CALL get_buffer(evc, lrwfc, iuwfc, ikk)
    CALL get_buffer(evq, lrwfc, iuwfc, ikq)
  ENDIF
  !
  ! Setup for Sternheimer solver
  !
  ! compute beta functions and kinetic energy for k-point ikq
  ! needed by h_psi, called by ch_psi_all, called by cgsolve_all
  !
  CALL init_us_2(npwq, igk_k(1, ikq), xk(1, ikq), vkb)
  CALL g2_kin(ikq)
  !
  ! compute preconditioning matrix h_diag used by cgsolve_all
  CALL h_prec(ik, evq, h_diag)
  !
  DO imode = 1, nmodes
    dvpsi = (0.d0, 0.d0)
    dvpsi(:,ib_ahc_min:ib_ahc_max) = dvpsi_cart(:,:,imode)
    !
    ! Apply -P_c^+ : orthogonalize dvpsi to all nbnd states
    CALL orthogonalize(dvpsi, evq, ikk, ikq, dpsi, npwq, .FALSE.)
    !
    conv_root = .true.
    dpsi(:,:) = (0.d0, 0.d0)
    !
    ! Iteratively solve the Sternheimer equation
    CALL cgsolve_all(ch_psi_all, cg_psi, et(ib_ahc_min, ikk), &
        dvpsi(1, ib_ahc_min), dpsi, h_diag, npwx, npwq, thresh, ik, &
        lter, conv_root, anorm, ahc_nbnd, npol)
    !
    IF (.NOT. conv_root) THEN
      WRITE( stdout, '(5x,"kpoint",i4," mode",i4,  &
            & " ahc_upperfan: cgsolve_all root not converged ", &
            & es10.3)') ik , imode, anorm
    ENDIF
    !
    dpsi_cart(:,:,imode) = dpsi(:,1:ahc_nbnd)
    !
  ENDDO ! imode
  !
  ! Calculate ahc_upfan(:,:,ik,:,:)
  !
  ahc_upfan = (0.d0, 0.d0)
  DO jmode = 1, nmodes
     DO imode = 1, nmodes
        CALL ZGEMM('C', 'N', ahc_nbnd, ahc_nbnd, npwx * npol, &
           (1.d0, 0.d0), dpsi_cart(1, 1, imode), npwx * npol, &
                        dvpsi_cart(1, 1, jmode), npwx * npol, &
           (0.d0, 0.d0), ahc_upfan(1, 1, imode, jmode), ahc_nbnd)
     ENDDO
  ENDDO
  !
  CALL mp_sum(ahc_upfan, intra_pool_comm)
  !
  ! Write ahc_upfan to file
  !
  IF (me_pool == root_pool) THEN
     nrec = ik + nbase_ik
     WRITE(iunupfan, REC=nrec) ahc_upfan
  ENDIF
  !
  ! Check: (H - e(ib)) dpsi_cart(:,ib,imode) = - P_c dvpsi_cart(:,ib,imode)
  imode = 2
  dvpsi = (0.d0, 0.d0)
  dvpsi(:, ib_ahc_min:ib_ahc_max) = dvpsi_cart(:, :, imode)
  CALL orthogonalize(dvpsi, evq, ikk, ikq, dpsi, npwq, .FALSE.)
  dpsi = (0.d0, 0.d0)
  !
  CALL h_psi(npwx, npwq, ahc_nbnd, dpsi_cart(1,1,imode), dpsi)
  DO ibnd = 1, ahc_nbnd
    dpsi(:,ibnd) = dpsi(:,ibnd) - et(ibnd + ahc_nbndskip, ikk) * dpsi_cart(:,ibnd,imode)
  ENDDO
  !
  ! Check total error is small
  !
  error_sum = SUM(ABS(dpsi(:,1:ahc_nbnd) - dvpsi(:,ib_ahc_min:ib_ahc_max)))
  CALL mp_sum(error_sum, intra_pool_comm)
  IF (ionode) THEN
    IF (error_sum > 1.d-4 * REAL(ahc_nbnd, DP)) THEN
      WRITE(stdout, *) 'WARNING: Sternheimer equation error is large.'
      WRITE(stdout, *) 'Consider increasing nbnd of the NSCF calculation.'
      WRITE(stdout, *) 'sum of absolute error = ', error_sum
      WRITE(stdout, *) 'sum of absolute value = ', &
          SUM(ABS(dvpsi(:,ib_ahc_min:ib_ahc_max)))
    ENDIF
  ENDIF
  !
  ! To print the error of each bands, uncomment the following block.
  !
  ! DO ibnd = 1, nbnd
  !   WRITE(stdout, *) 'Sternheimer error ', ibnd, &
  !       SUM(ABS(dpsi(:,ibnd) - dvpsi(:,ibnd)))
  ! ENDDO
  !
  DEALLOCATE(h_diag)
  DEALLOCATE(dpsi_cart)
  !
  CALL stop_clock('ahc_upfan')
  !
END SUBROUTINE ahc_do_upperfan
!----------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
SUBROUTINE ahc_do_dw(ik)
!------------------------------------------------------------------------------
!!
!! Compute and write to file the matrix elements of [dV,p]. This quantity is
!! needed to compute Debye-Waller self-energy within rigid-ion approximation.
!!
!! ahc_dw(ib, jb, imode, jdir)
!!   = i * <\psi_ib(k)|[dV_{SCF}/du_{Gamma, imode}, p_jdir]|\psi_jb(k)>
!!
!! Implements Eq.(3) of PHonon/Doc/dfpt_self_energy.pdf
!!
!! Here, the "operator-generalized acoustic sum rule" is used to represent
!! Debye-Waller self-energy as a simple matrix element.
!! See Eq.(13) of the following reference:
!! Jae-Mo Lihm and Cheol-Hwan Park, Phys. Rev. B, 101, 121102(R) (2020).
!!
!------------------------------------------------------------------------------
  USE kinds,            ONLY : DP
  USE io_global,        ONLY : stdout
  USE mp,               ONLY : mp_sum
  USE mp_pools,         ONLY : intra_pool_comm, me_pool, root_pool
  USE wavefunctions,    ONLY : evc
  USE wvfct,            ONLY : npwx
  USE klist,            ONLY : ngk, igk_k, xk
  USE gvect,            ONLY : g
  USE noncollin_module, ONLY : noncolin
  USE cell_base,        ONLY : tpiba
  USE noncollin_module, ONLY : npol
  USE buffers,          ONLY : get_buffer
  USE qpoint,           ONLY : ikks, ikqs
  USE modes,            ONLY : nmodes
  USE units_lr,         ONLY : lrwfc, iuwfc
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ik
  !! k point index where dvpsi is calculated
  !
  INTEGER :: ikk
  !! Counter on k-point
  INTEGER :: ikq
  !! Counter on k+q-point
  INTEGER :: imode
  !! Counter on modes
  INTEGER :: ibnd
  !! Counter on bands
  INTEGER :: idir
  !! Counter on Cartesian directions
  INTEGER :: ig
  !! Counter on plane waves
  INTEGER :: npw
  !! number of plane waves at k
  INTEGER :: npwq
  !! number of plane waves at k+q
  INTEGER :: nrec
  !! Record to write file
  COMPLEX(DP), ALLOCATABLE :: p_psi(:, :, :)
  !! wavefunction multiplied by the momentum operator p
  COMPLEX(DP), ALLOCATABLE :: p_psi_gauged(:, :, :)
  !! wavefunction multiplied by the momentum operator p after gauge fixing
  !
  CALL start_clock('ahc_dw')
  !
  WRITE(stdout, '(5x,a,I8)') 'Computing ahc_dw for ik = ', ik
  !
  ALLOCATE(p_psi(npwx * npol, ahc_nbnd_gauge, 3))
  ALLOCATE(p_psi_gauged(npwx * npol, ahc_nbnd, 3))
  !
  ikk = ikks(ik)
  ikq = ikqs(ik)
  npw = ngk(ikk)
  npwq = ngk(ikq)
  !
  ! Read unperturbed wavefunctions psi(k) from file
  ! No need of k+q because q == 0.
  !
  CALL get_buffer(evc, lrwfc, iuwfc, ikk)
  !
  ! Compute p_psi(:,ibnd,idir) = p_idir |evc(:,ibnd)>
  ! p_idir = -i d/dr_idir
  !
  p_psi = (0.d0, 0.d0)
  DO idir = 1, 3
    DO ibnd = 1, ahc_nbnd_gauge
      DO ig = 1, npw
        p_psi(ig, ibnd, idir) = evc(ig, ibnd + ib_ahc_gauge_min - 1) &
          * ( xk(idir, ik) + g(idir, igk_k(ig, ik)) ) * tpiba
      ENDDO
      IF (noncolin) THEN
        DO ig = 1, npw
          p_psi(ig + npwx, ibnd, idir) &
            = evc(ig + npwx, ibnd + ib_ahc_gauge_min - 1) &
            * ( xk(idir, ik) + g(idir, igk_k(ig, ik)) ) * tpiba
        ENDDO
      ENDIF
    ENDDO
  ENDDO
  !
  ! Fix gauge of p_psi: dvpsi_gauged = p_psi * psi_gauge
  !
  DO idir = 1, 3
    CALL ZGEMM('N', 'N', npwx*npol, ahc_nbnd, ahc_nbnd_gauge, &
               (1.d0,0.d0), p_psi(1, 1, idir), npwx * npol, &
                            psi_gauge, ahc_nbnd_gauge, &
               (0.d0,0.d0), p_psi_gauged(1, 1, idir), npwx * npol)
  ENDDO
  !
  ! Calculate ahc_dw(ib,jb,imode,jdir)
  ! = i * <\psi_ib(k)|[dV_{SCF}/du^Gamma_{imode}, p_jdir]|\psi_jb(k)>
  !
  ahc_dw = (0.d0, 0.d0)
  DO idir = 1, 3
    DO imode = 1, nmodes
      CALL ZGEMM('C', 'N', ahc_nbnd, ahc_nbnd, npwx * npol, &
        (0.d0, 1.d0), dvpsi_cart(1,1,imode), npwx * npol, &
                      p_psi_gauged(1,1,idir), npwx * npol, &
        (1.d0, 0.d0), ahc_dw(1,1,imode,idir), ahc_nbnd)
      !
      CALL ZGEMM('C', 'N', ahc_nbnd, ahc_nbnd, npwx * npol, &
        (0.d0,-1.d0), p_psi_gauged(1,1,idir), npwx * npol, &
                      dvpsi_cart(1,1,imode), npwx * npol, &
        (1.d0, 0.d0), ahc_dw(1,1,imode,idir), ahc_nbnd)
    ENDDO
  ENDDO
  !
  CALL mp_sum(ahc_dw, intra_pool_comm)
  !
  ! Write ahc_dw to file
  IF (me_pool == root_pool) THEN
    nrec = ik + nbase_ik
    WRITE(iundw, REC=nrec) ahc_dw
  ENDIF
  !
  DEALLOCATE(p_psi)
  DEALLOCATE(p_psi_gauged)
  !
  CALL stop_clock('ahc_dw')
  !
!------------------------------------------------------------------------------
END SUBROUTINE ahc_do_dw
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
SUBROUTINE ahc_do_gkk(ik)
!------------------------------------------------------------------------------
!!
!! Calculate and write to file ahc_gkk.
!! ahc_gkk(ib, jb, imode) = <\psi_ib(k+q)|dV/du_{q, imode}|\psi_jb(k)>
!! 1 <= ib <= nbnd, ib_ahc_min <= jb <= ib_ahc_max
!!
!! Implements Eq.(2) of PHonon/Doc/dfpt_self_energy.pdf
!!
!------------------------------------------------------------------------------
   USE kinds,        ONLY : DP
   USE io_global,    ONLY : stdout
   USE mp,           ONLY : mp_sum
   USE mp_pools,     ONLY : intra_pool_comm, me_pool, root_pool
   USE modes,        ONLY : nmodes
   USE wvfct,        ONLY : npwx, nbnd
   USE qpoint,       ONLY : ikqs
   USE buffers,      ONLY : get_buffer
   USE noncollin_module, ONLY : npol
   USE units_lr,     ONLY : lrwfc, iuwfc
   USE eqv,          ONLY : evq
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN) :: ik
   !! k point index where dvpsi is calculated
   !
   INTEGER :: imode, nrec
   !
   CALL start_clock('ahc_gkk')
   !
   WRITE(stdout, '(5x,a,I8)') 'Computing ahc_gkk for ik = ', ik
   !
   ahc_gkk = (0.d0, 0.d0)
   CALL get_buffer(evq, lrwfc, iuwfc, ikqs(ik))
   !
   DO imode = 1, nmodes
      CALL ZGEMM('C', 'N', nbnd, ahc_nbnd, npwx*npol, &
         (1.d0,0.d0), evq, npwx*npol, dvpsi_cart(1, 1, imode), npwx*npol, &
         (0.d0,0.d0), ahc_gkk(1, 1, imode), nbnd)
   ENDDO
   !
   CALL mp_sum(ahc_gkk, intra_pool_comm)
   !
   ! Write ahc_gkk to file
   !
   IF (me_pool == root_pool) THEN
      nrec = ik + nbase_ik
      WRITE(iungkk, REC=nrec) ahc_gkk
   ENDIF
   !
   CALL stop_clock('ahc_gkk')
   !
END SUBROUTINE ahc_do_gkk
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
SUBROUTINE compute_psi_gauge(ik)
!------------------------------------------------------------------------------
!! Compute psi_gauge(n, m) = <psi_nk(q) | psi_ref_mk>.
!! n = ib_ahc_gauge_min, ..., ib_ahc_gauge_max
!! m = ib_ahc_min, ..., ib_ahc_max
!!
!! The range of m is wider than the range of n because of possible degeneracy
!! at ib_ahc_min or ib_ahc_max.
!!
!! The overlap is computed indirectly by comparing psi_nk at a few G vectors.
!! G vectors with miller indices (Gx, Gy, Gz) with -3 <= Gx, Gy, Gz <= 3
!! are used.
!!
!! P: projection operator to the selected G vectors.
!! P|psi_ref_m> = P|psi_n> * psi_gauge(n, m)
!! S_nm = <psi_ref_n|P|psi_ref_m>
!! A_nm = <psi_n|P|psi_ref_m>
!! psi_gauge(n, m) = (A * inv(S))_nm
!!
!! Implemented only for NCPP case.
!------------------------------------------------------------------------------
  USE kinds,            ONLY : DP
  USE io_files,         ONLY : seqopn
  USE mp,               ONLY : mp_bcast, mp_sum
  USE mp_pools,         ONLY : me_pool, root_pool, intra_pool_comm
  USE buffers,          ONLY : get_buffer, save_buffer
  USE matrix_inversion, ONLY : invmat
  USE wavefunctions,    ONLY : evc
  USE gvect,            ONLY : mill, ngm
  USE wvfct,            ONLY : npwx, et
  USE klist,            ONLY : igk_k, ngk
  USE noncollin_module, ONLY : noncolin, npol
  USE qpoint,           ONLY : ikks
  USE disp,             ONLY : nqs
  USE control_ph,       ONLY : current_iq
  USE units_lr,         ONLY : iuwfc, lrwfc
  USE units_ph,         ONLY : iugauge
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ik
  !! Current k point index
  !
  INTEGER, PARAMETER :: gmax_search = 3
  !! Max |G_x,y,z| to use for gauge fixing.
  LOGICAL :: is_last_group
  !! true if at the last degenerate group
  INTEGER :: info
  !! Info output of zheev
  INTEGER :: found_g
  !! set to 1 if given G vector is found in mill
  INTEGER :: ikk
  !! Counter on k point
  INTEGER :: igx
  !! Counter on plane waves
  INTEGER :: igy
  !! Counter on plane waves
  INTEGER :: igz
  !! Counter on plane waves
  INTEGER :: ng_gauge
  !! Number of G vectors to use in computing overlap
  INTEGER :: ig
  !! Counter on plane waves
  INTEGER :: jg
  !! Counter on plane waves
  INTEGER :: igroup
  !! Counter on degenerate groups
  INTEGER :: ndegen
  !! Degree of degeneracy
  INTEGER :: ndegen_ref
  !! Degree of degeneracy
  INTEGER :: npw
  !! Number of plane waves
  INTEGER :: ibnd
  !! Counter on bands
  INTEGER :: ibnd_ref
  !! Counter on bands
  INTEGER :: jbnd
  !! Counter on bands
  INTEGER :: lwork
  !! dimension of workspace for diagonalization
  INTEGER, ALLOCATABLE :: gauge_mill_tmp(:, :)
  !! G-points to use for computing gauge
  INTEGER, ALLOCATABLE :: gauge_mill(:, :)
  !! G-points to use for computing gauge
  INTEGER, ALLOCATABLE :: ndegen_list(:)
  !! Degree of degeneracy for each degenerate groups.
  REAL(DP), ALLOCATABLE :: rwork(:)
  !! workspace for diagonalization
  REAL(DP), ALLOCATABLE :: eigvals(:)
  !! eigenvalues of smat
  COMPLEX(DP), ALLOCATABLE :: evc_select_ref(:, :, :)
  !! Reference wavefunction at selected G vectors
  COMPLEX(DP), ALLOCATABLE :: evc_select(:, :, :)
  !! Wavefunction at selected G vectors
  COMPLEX(DP), ALLOCATABLE :: smat(:, :)
  !! overlap matrix of psi_ref
  COMPLEX(DP), ALLOCATABLE :: smat_inv(:, :)
  !! inverse of smat
  COMPLEX(DP), ALLOCATABLE :: amat(:, :)
  !! overlap matrix of psi_ref and psi
  COMPLEX(DP), ALLOCATABLE :: work(:)
  !! workspace for diagonalization
  !
  CALL start_clock('ahc_gauge')
  !
  ! If there is only a single q point, no need of gauge fixing.
  ! Set psi_gauge to identity and return.
  !
  IF (nqs == 1) THEN
    psi_gauge = (0.d0, 0.d0)
    !
    DO ibnd = 1, ahc_nbnd
       psi_gauge(ibnd - ib_ahc_gauge_min + ib_ahc_min, ibnd) = (1.d0, 0.d0)
    ENDDO
    !
    CALL stop_clock('ahc_gauge')
    RETURN
    !
  ENDIF
  !
  psi_gauge = (0.d0, 0.d0)
  !
  ikk = ikks(ik)
  npw = ngk(ikk)
  !
  CALL get_buffer(evc, lrwfc, iuwfc, ikk)
  !
  ! If iq == 1,
  ! 1) Set gauge_mill with the G vectors to use to compute overlap.
  ! 2) Group degenerate bands and set ndegen_list.
  ! 3) Find selected G vectors from each pools and set evc_select_ref
  ! 4) Write gauge_mill, ndegen_list, and evc_select_ref to file.
  ! 5) For each degenerate group, compute inv(S) matrix and write to file.
  ! 6) Set psi_gauge to identity and return.
  !
  IF (current_iq == 1) THEN
    !
    ! 1) Set gauge_mill.
    !
    ALLOCATE(gauge_mill_tmp(3, (2 * gmax_search + 1)**3))
    gauge_mill_tmp = 0
    ng_gauge = 0
    !
    DO igx = -gmax_search, gmax_search
      DO igy = -gmax_search, gmax_search
        DO igz = -gmax_search, gmax_search
          !
          ! Check if (igx, igy, igz) is in mill_g
          !
          found_g = 0
          DO ig = 1, ngm
            IF (mill(1, ig) == igx .AND. &
                mill(2, ig) == igy .AND. &
                mill(3, ig) == igz) THEN
              found_g = 1
              EXIT
            ENDIF
          ENDDO
          !
          CALL mp_sum(found_g, intra_pool_comm)
          !
          ! found_g must be 0 (not found) or 1 (found)
          !
          IF (found_g /= 0 .AND. found_g /= 1) CALL errore(&
            'compute_psi_gauge', 'problem setting gauge_mill', 1)
          !
          IF (found_g == 1) THEN
            ng_gauge = ng_gauge + 1
            gauge_mill_tmp(:, ng_gauge) = (/ igx, igy, igz /)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    !
    ! Make a smaller array containing found G vectors only.
    !
    ALLOCATE(gauge_mill(3, ng_gauge))
    gauge_mill(:, 1:ng_gauge) = gauge_mill_tmp(:, 1:ng_gauge)
    DEALLOCATE(gauge_mill_tmp)
    !
    ! 2) Group degenerate bands and set ndegen_list
    !
    ALLOCATE(ndegen_list(ahc_nbnd))
    ndegen_list = 0
    igroup = 1
    !
    DO ibnd = ib_ahc_min, ib_ahc_max
      !
      ! If nondegenerate, increase igroup
      !
      IF (ibnd > ib_ahc_min) THEN
        IF (et(ibnd, ikk) - et(ibnd - 1, ikk) > e_degen_thr) THEN
          igroup = igroup + 1
        ENDIF
      ENDIF
      !
      ndegen_list(igroup) = ndegen_list(igroup) + 1
      !
    ENDDO ! ibnd
    !
    ! 3) Find selected G vectors from each pools and set evc_select_ref.
    !
    ! In case of PW parallelization, evc_select_ref is set by the processor that
    ! contains the G vector and is later summed by mp_sum.
    !
    ALLOCATE(evc_select_ref(ng_gauge, npol, ahc_nbnd))
    evc_select_ref = (0.d0, 0.d0)
    !
    DO ig = 1, ng_gauge
      DO jg = 1, npw
        IF (mill(1, igk_k(jg, ikk)) == gauge_mill(1, ig) .AND. &
            mill(2, igk_k(jg, ikk)) == gauge_mill(2, ig) .AND. &
            mill(3, igk_k(jg, ikk)) == gauge_mill(3, ig)) THEN
          !
          DO ibnd = 1, ahc_nbnd
            jbnd = ibnd - 1 + ib_ahc_min
            evc_select_ref(ig, 1, ibnd) = evc(jg, jbnd)
            IF (noncolin) THEN
              evc_select_ref(ig, 2, ibnd) = evc(jg + npwx, jbnd)
            ENDIF
          ENDDO
          !
          EXIT
          !
        ENDIF
      ENDDO ! jg
    ENDDO ! ig
    !
    CALL mp_sum(evc_select_ref, intra_pool_comm)
    !
    ! 4) Write gauge_mill, ndegen_list, and evc_select_ref to file.
    !
    IF (me_pool == root_pool) THEN
      WRITE(iugauge) ng_gauge
      WRITE(iugauge) gauge_mill
      WRITE(iugauge) ndegen_list
      WRITE(iugauge) evc_select_ref
    ENDIF
    !
    ! 5) For each degenerate group, compute inv(S) matrix and write to file.
    !
    IF (me_pool == root_pool) THEN
      !
      ibnd = 1
      DO igroup = 1, ahc_nbnd
        !
        IF (ndegen_list(igroup) == 0) EXIT
        ndegen = ndegen_list(igroup)
        !
        ! S_nm = sum_G [evc_select_ref(ig, n)]* evc_select_ref(ig, m)
        !
        ALLOCATE(smat(ndegen, ndegen))
        ALLOCATE(smat_inv(ndegen, ndegen))
        !
        CALL ZGEMM('C', 'N', ndegen, ndegen, npol * ng_gauge, &
          (1.d0, 0.d0), evc_select_ref(1, 1, ibnd), npol * ng_gauge, &
                        evc_select_ref(1, 1, ibnd), npol * ng_gauge, &
          (0.d0, 0.d0), smat, ndegen)
        !
        CALL invmat(ndegen, smat, smat_inv)
        !
        ! Check whether smat is singular by computing the eigenvalues
        !
        lwork = 2 * ndegen - 1
        ALLOCATE(work(lwork))
        ALLOCATE(eigvals(ndegen))
        ALLOCATE(rwork(3 * ndegen - 2))
        CALL ZHEEV('N', 'U', ndegen, smat, ndegen, eigvals, work, lwork, rwork, info)
        IF (info /= 0) THEN
          CALL errore('compute_psi_gauge', 'problem diagonalizing smat', info)
        ENDIF
        !
        IF (eigvals(1) < 1.d-2) THEN
          CALL errore('compute_psi_gauge', 'smat close to singluar. &
            &Try increasing gmax_search.', 1)
        ENDIF
        !
        WRITE(iugauge) smat_inv
        !
        DEALLOCATE(work)
        DEALLOCATE(eigvals)
        DEALLOCATE(rwork)
        DEALLOCATE(smat)
        DEALLOCATE(smat_inv)
        !
        ibnd = ibnd + ndegen
        !
      ENDDO
    ENDIF
    !
    IF (ibnd /= ahc_nbnd + 1) CALL errore('compute_psi_gauge', &
    'ibnd /= ahc_nbnd + 1 after loop over degenreate groups at first iq', 1)
    !
    ! 6) Set psi_gauge to identity and return.
    !
    psi_gauge = (0.d0, 0.d0)
    !
    DO ibnd = 1, ahc_nbnd
       psi_gauge(ibnd - ib_ahc_gauge_min + ib_ahc_min, ibnd) = (1.d0, 0.d0)
    ENDDO
    !
    DEALLOCATE(gauge_mill)
    DEALLOCATE(ndegen_list)
    DEALLOCATE(evc_select_ref)
    !
    CALL stop_clock('ahc_gauge')
    !
    RETURN
    !
  ENDIF
  !
  ! If iq > 1,
  ! 1) Read wfcgauge file from root_pool and bcast.
  ! 2) Find selected G vectors from each pools and set evc_select.
  ! 3) For each degenerate group:
  !   3-1) Read inv(S) matrix from file
  !   3-2) Compute A matrix
  !   3-3) Compute psi_gauge.
  !
  ! 1) Read wfcgauge file from root_pool and bcast.
  !
  IF (me_pool == root_pool) READ(iugauge) ng_gauge
  CALL mp_bcast(ng_gauge, root_pool, intra_pool_comm)
  !
  ALLOCATE(ndegen_list(ahc_nbnd))
  ALLOCATE(gauge_mill(3, ng_gauge))
  ALLOCATE(evc_select_ref(ng_gauge, npol, ahc_nbnd))
  ALLOCATE(evc_select(ng_gauge, npol, ahc_nbnd_gauge))
  !
  IF (me_pool == root_pool) THEN
    READ(iugauge) gauge_mill
    READ(iugauge) ndegen_list
    READ(iugauge) evc_select_ref
  ENDIF
  !
  CALL mp_bcast(gauge_mill, root_pool, intra_pool_comm)
  CALL mp_bcast(ndegen_list, root_pool, intra_pool_comm)
  CALL mp_bcast(evc_select_ref, root_pool, intra_pool_comm)
  !
  ! 2) Find selected G vectors from each pools and set evc_select.
  !
  evc_select = (0.d0, 0.d0)
  !
  DO ig = 1, ng_gauge
    DO jg = 1, npw
      IF (mill(1, igk_k(jg, ikk)) == gauge_mill(1, ig) .AND. &
          mill(2, igk_k(jg, ikk)) == gauge_mill(2, ig) .AND. &
          mill(3, igk_k(jg, ikk)) == gauge_mill(3, ig)) THEN
        !
        DO ibnd = 1, ahc_nbnd_gauge
          jbnd = ibnd - 1 + ib_ahc_gauge_min
          evc_select(ig, 1, ibnd) = evc(jg, jbnd)
          IF (noncolin) THEN
            evc_select(ig, 2, ibnd) = evc(jg + npwx, jbnd)
          ENDIF
        ENDDO
        !
        EXIT
        !
      ENDIF
    ENDDO ! jg
  ENDDO ! ig
  !
  CALL mp_sum(evc_select, intra_pool_comm)
  !
  ! 3) For each degenerate group:
  !
  ibnd_ref = 1
  ibnd = 1
  !
  DO igroup = 1, ahc_nbnd
    !
    IF (ndegen_list(igroup) == 0) EXIT
    ndegen_ref = ndegen_list(igroup)
    !
    is_last_group = .FALSE.
    IF (igroup == ahc_nbnd) THEN
      is_last_group = .TRUE.
    ELSEIF (ndegen_list(igroup + 1) == 0) THEN
      is_last_group = .TRUE.
    ENDIF
    !
    ! Enlarge ndegen when at the first or last degenerate group.
    !
    IF (igroup == 1) THEN
      ! First group.
      ndegen = ndegen_ref + ib_ahc_min - ib_ahc_gauge_min
    ELSEIF (is_last_group) THEN
      ! Last group.
      ndegen = ndegen_ref - ib_ahc_max + ib_ahc_gauge_max
    ELSE
      ndegen = ndegen_ref
    ENDIF
    !
    ALLOCATE(smat_inv(ndegen_ref, ndegen_ref))
    ALLOCATE(amat(ndegen, ndegen_ref))
    !
    ! 3-1) Read inv(S) matrix from file
    ! S_nm = sum_G [evc_select_ref(ig, n)]* evc_select_ref(ig, m)
    !
    IF (me_pool == root_pool) READ(iugauge) smat_inv
    CALL mp_bcast(smat_inv, root_pool, intra_pool_comm)
    !
    ! 3-2) Compute A matrix
    ! A_nm = sum_G [evc_select(ig, n)]* evc_select_ref(ig, m)
    !
    CALL ZGEMM('C', 'N', ndegen, ndegen_ref, npol * ng_gauge, &
        (1.d0, 0.d0), evc_select(1, 1, ibnd), npol * ng_gauge, &
                      evc_select_ref(1, 1, ibnd_ref), npol * ng_gauge, &
        (0.d0, 0.d0), amat, ndegen)
    !
    ! 3-3) Compute psi_gauge
    ! psi_gauge(n, m) = (A * inv(S))_nm
    !
    CALL ZGEMM('N', 'N', ndegen, ndegen_ref, ndegen_ref, &
        (1.d0, 0.d0), amat, ndegen, smat_inv, ndegen_ref, &
        (0.d0, 0.d0), psi_gauge(ibnd, ibnd_ref), ahc_nbnd_gauge)
    !
    DEALLOCATE(smat_inv)
    DEALLOCATE(amat)
    ibnd_ref = ibnd_ref + ndegen_ref
    ibnd = ibnd + ndegen
    !
  ENDDO ! igroup
  !
  IF (ibnd_ref /= ahc_nbnd + 1) CALL errore('compute_psi_gauge', &
    'ibnd_ref /= ahc_nbnd + 1 after loop over degenreate groups', 1)
  IF (ibnd /= ahc_nbnd_gauge + 1) CALL errore('compute_psi_gauge', &
    'ibnd /= ahc_nbnd_gauge + 1 after loop over degenreate groups', 1)
  !
  DEALLOCATE(gauge_mill)
  DEALLOCATE(ndegen_list)
  DEALLOCATE(evc_select)
  DEALLOCATE(evc_select_ref)
  !
  CALL stop_clock('ahc_gauge')
  !
END SUBROUTINE compute_psi_gauge
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
SUBROUTINE elph_ahc_setup()
!------------------------------------------------------------------------------
!! Initialize variables for ahc calculation.
!! Look for degeneracy around ib_ahc_min or ib_ahc_max and set
!! ahc_nbnd_gauge, the number of bands to compute dvpsi
!------------------------------------------------------------------------------
  USE kinds,        ONLY : DP
  USE io_global,    ONLY : ionode, ionode_id
  USE io_files,     ONLY : create_directory
  USE mp,           ONLY : mp_min, mp_max, mp_bcast
  USE mp_images,    ONLY : intra_image_comm
  USE mp_pools,     ONLY : intra_pool_comm, root_pool, me_pool
  USE wvfct,        ONLY : nbnd, et
  USE qpoint,       ONLY : nksq, ikks, nksqtot
  USE control_lr,   ONLY : lgamma
  USE control_ph,   ONLY : current_iq
  !
  IMPLICIT NONE
  !
  LOGICAL :: exst
  !! True if folder exists
  CHARACTER(LEN=256) :: filoutetk
  !! Filename for e_n(k) energy eigenvalue output
  CHARACTER(LEN=256) :: filoutetq
  !! Filename for e_n(k+q) energy eigenvalue output
  INTEGER :: ik
  !! Counter for k points
  INTEGER :: ikk
  !! Counter for k points
  INTEGER :: ibnd
  !! Counter for bands
  INTEGER :: ib_ahc_gauge_min_ik
  !! Counter for bands
  INTEGER :: ib_ahc_gauge_max_ik
  !! Counter for bands
  INTEGER :: recl
  !! length of the files to be written
  INTEGER :: iunetk
  !! Unit for e_n(k) energy eigenvalue output
  INTEGER :: iunetq
  !! Unit for e_n(k+q) energy eigenvalue output
  REAL(DP), ALLOCATABLE :: et_collect(:, :, :)
  !! energy eigenvalues at k and k+q for all k points
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  INTEGER, EXTERNAL :: find_free_unit
  !
  IF (ahc_nbndskip + ahc_nbnd > nbnd) CALL errore('elph_ahc_setup', &
     'ahc_nbndskip + ahc_nbnd cannot be greater than nbnd', 1)
  !
  ! Create output directory
  !
  IF (ionode) INQUIRE(FILE=TRIM(ahc_dir), EXIST=exst)
  CALL mp_bcast(exst, ionode_id, intra_image_comm)
  IF (.NOT. exst) CALL create_directory(ahc_dir)
  !
  ib_ahc_min = ahc_nbndskip + 1
  ib_ahc_max = ahc_nbndskip + ahc_nbnd
  !
  ! Search degeneracy around ib_ahc_min or ib_ahc_max to set
  ! ib_ahc_gauge_min, ib_ahc_gauge_max, and ahc_nbnd_gauge.
  !
  ! Compute at each root_pool and bcast to other nodes.
  ! ib_ahc_gauge_min/max are different for each pool.
  !
  IF (me_pool == root_pool) THEN
    ib_ahc_gauge_min = ib_ahc_min
    ib_ahc_gauge_max = ib_ahc_max
    !
    DO ik = 1, nksq
      ikk = ikks(ik)
      !
      ib_ahc_gauge_min_ik = ib_ahc_min
      IF (ib_ahc_min > 1) THEN
        DO ibnd = ib_ahc_min - 1, 1, -1
          IF (ABS(et(ibnd, ikk) - et(ib_ahc_min, ikk)) < e_degen_thr) THEN
            ib_ahc_gauge_min_ik = ibnd
          ELSE
            EXIT
          ENDIF
        ENDDO
      ENDIF
      !
      ib_ahc_gauge_max_ik = ib_ahc_max
      IF (ib_ahc_max < nbnd) THEN
        DO ibnd = ib_ahc_max + 1, nbnd
          IF (ABS(et(ibnd, ikk) - et(ib_ahc_max, ikk)) < e_degen_thr) THEN
            ib_ahc_gauge_max_ik = ibnd
          ELSE
            EXIT
          ENDIF
        ENDDO
      ENDIF
      !
      ib_ahc_gauge_min = MIN(ib_ahc_gauge_min_ik, ib_ahc_gauge_min)
      ib_ahc_gauge_max = MAX(ib_ahc_gauge_max_ik, ib_ahc_gauge_max)
    ENDDO
  ENDIF ! root_pool
  !
  CALL mp_bcast(ib_ahc_gauge_min, root_pool, intra_pool_comm)
  CALL mp_bcast(ib_ahc_gauge_max, root_pool, intra_pool_comm)
  !
  ahc_nbnd_gauge = ib_ahc_gauge_max - ib_ahc_gauge_min + 1
  !
  ! Write energy eigenvalues to file
  !
  IF (ionode) THEN
    !
    ALLOCATE(et_collect(nbnd, nksqtot, 2))
    !
    filoutetk = TRIM(ahc_dir) // 'ahc_etk_iq' // TRIM(int_to_char(current_iq)) // '.bin'
    filoutetq = TRIM(ahc_dir) // 'ahc_etq_iq' // TRIM(int_to_char(current_iq)) // '.bin'
    INQUIRE(IOLENGTH=recl) et_collect(:, :, 1)
    !
    iunetk = find_free_unit()
    OPEN(UNIT=iunetk, FILE=TRIM(filoutetk), FORM='unformatted', &
        ACCESS='direct', RECL=recl)
    !
    iunetq = find_free_unit()
    OPEN(UNIT=iunetq, FILE=TRIM(filoutetq), FORM='unformatted', &
        ACCESS='direct', RECL=recl)
    !
    DO ik = 1, nksqtot
      IF (lgamma) THEN
        et_collect(:, ik, 1) = et(:, ik)
        et_collect(:, ik, 2) = et(:, ik)
      ELSE
        et_collect(:, ik, 1) = et(:, 2*ik-1)
        et_collect(:, ik, 2) = et(:, 2*ik)
      ENDIF
    ENDDO
    !
    WRITE(iunetk, REC=1) et_collect(:,:,1)
    WRITE(iunetq, REC=1) et_collect(:,:,2)
    !
    CLOSE(iunetk, STATUS='KEEP')
    CLOSE(iunetq, STATUS='KEEP')
    !
    DEALLOCATE(et_collect)
    !
  ENDIF
  !
END SUBROUTINE elph_ahc_setup
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
SUBROUTINE elph_do_ahc()
!------------------------------------------------------------------------------
!! The main driver of the AHC calculation.
!------------------------------------------------------------------------------
  USE kinds,        ONLY : DP
  USE io_global,    ONLY : stdout
  USE io_files,     ONLY : diropn
  USE mp_pools,     ONLY : npool, my_pool_id, me_pool, root_pool
  USE wvfct,        ONLY : npwx, nbnd
  USE klist,        ONLY : lgauss
  USE noncollin_module, ONLY : npol
  USE buffers,      ONLY : get_buffer
  USE qpoint,       ONLY : nksq, nksqtot
  USE modes,        ONLY : u, nmodes
  USE control_lr,   ONLY : lgamma, nbnd_occ
  USE units_ph,     ONLY : lrdvpsi, iudvpsi
  USE eqv,          ONLY : dvpsi
  USE control_ph,   ONLY : current_iq
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=6) :: iq_name
  !! int_to_char(current_iq)
  CHARACTER(LEN=256) :: filoutgkk
  !! Filename for ahc_gkk (dV_q matrix element) output
  CHARACTER(LEN=256) :: filoutdw
  !! Filename for ahc_dw (Debye-Waller) output
  CHARACTER(LEN=256) :: filoutupfan
  !! Filename for ahc_upfan (upper Fan by Sternheimer) output
  INTEGER :: ik
  !! Counter for k points
  INTEGER :: imode
  !! Counter for atomic displacements
  INTEGER :: jmode
  !! Counter for atomic displacements
  INTEGER :: recl
  !! length of the files to be written
  INTEGER :: nrec
  !! record number
  INTEGER :: rest
  !! integer for calculating nbase_ik
  INTEGER :: nks1
  !! integer for calculating nbase_ik
  ! CHARACTER(LEN=256) :: filename_wf
  COMPLEX(DP), ALLOCATABLE :: dvpsi_gauged(:, :)
  !! Temporary storage of dvpsi for rotating to Cartesian
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  INTEGER, EXTERNAL :: find_free_unit
  !
  WRITE(stdout, *) ""
  WRITE(stdout, '(5x,a)') "Begin electron-phonon calculation for AHC theory"
  !
  CALL start_clock('ahc_elph')
  !
  ! For k point parallelization, compute the global ik index
  ! ik_global = ik_local + nbase_ik
  ! Adapted from el_ph_collect.f90
  !
  IF (npool == 1) THEN
    nbase_ik = 0
  ELSE
    nks1 = nksqtot / npool
    !
    rest = nksqtot - nks1 * npool
    !
    IF ((my_pool_id + 1) <= rest) nks1 = nks1 + 1
    !
    IF (nks1 /= nksq) CALL errore('elph_do_ahc', 'problem with nks1', 1)
    !
    nbase_ik = nksq * my_pool_id
    !
    IF ((my_pool_id + 1) > rest) nbase_ik = nbase_ik + rest
  ENDIF ! npool
  !
  ALLOCATE(ahc_gkk(nbnd, ahc_nbnd, nmodes))
  ALLOCATE(ahc_upfan(ahc_nbnd, ahc_nbnd, nmodes, nmodes))
  ALLOCATE(ahc_dw(ahc_nbnd, ahc_nbnd, nmodes, 3))
  ALLOCATE(dvpsi_cart(npwx*npol, ahc_nbnd, nmodes))
  ALLOCATE(psi_gauge(ahc_nbnd_gauge, ahc_nbnd))
  !
  ! Open units for binary output. Each root_pool writes different records.
  !
  IF (me_pool == root_pool) THEN
    iq_name = int_to_char(current_iq)
    filoutgkk = TRIM(ahc_dir) // 'ahc_gkk_iq' // TRIM(iq_name) // '.bin'
    filoutupfan = TRIM(ahc_dir) // 'ahc_upfan_iq' // TRIM(iq_name) // '.bin'
    filoutdw = TRIM(ahc_dir) // 'ahc_dw.bin'
    !
    iungkk = find_free_unit()
    INQUIRE(IOLENGTH=recl) ahc_gkk
    OPEN(UNIT=iungkk, FILE=TRIM(filoutgkk), FORM='unformatted', &
         ACCESS='direct', RECL=recl)
    !
    IF (.NOT. skip_upperfan) THEN
      iunupfan = find_free_unit()
      INQUIRE(IOLENGTH=recl) ahc_upfan
      OPEN(UNIT=iunupfan, FILE=TRIM(filoutupfan), FORM='unformatted', &
           ACCESS='direct', RECL=recl)
    ENDIF
    !
    IF (lgamma) THEN
      iundw = find_free_unit()
      INQUIRE(IOLENGTH=recl) ahc_dw
      OPEN(UNIT=iundw, FILE=TRIM(filoutdw), FORM='unformatted', &
          ACCESS='direct', RECL=recl)
    ENDIF
  ENDIF ! root_pool
  !
  ! We need to change nbnd_occ because we solve the Sternheimer
  ! equation for all bands, not only for occupied bands.
  nbnd_occ(:) = nbnd
  !
  ! Metallic case: we do not want lgauss anymore.
  lgauss = .FALSE.
  !
  ! alpha_pv for Sternheimer equation changes because nbnd_occ is changed.
  CALL setup_alpha_pv()
  !
  DO ik = 1, nksq
    !
    ! Compute phase of wavefunctions psi_nk at q, using psi_nk at q=gamma
    !
    CALL compute_psi_gauge(ik)
    !
    dvpsi_cart = (0.d0, 0.d0)
    !
    ! Read dV_{SCF} psi from file and rotate from pattern to Cartesian basis
    ! Note that dvpsi is computed in elphon only for the relevant bands:
    ! ibnd = ib_ahc_gauge_min,..., ib_ahc_gauge_max (in total, ahc_nbnd_gauge)
    !
    ALLOCATE(dvpsi_gauged(npwx*npol, ahc_nbnd))
    !
    DO imode = 1, nmodes
      !
      nrec = (ik - 1) * nmodes + imode
      CALL get_buffer(dvpsi, lrdvpsi, iudvpsi, nrec)
      !
      ! Fix gauge of dvpsi: dvpsi_gauged = dvpsi * psi_gauge
      !
      CALL ZGEMM('N', 'N', npwx*npol, ahc_nbnd, ahc_nbnd_gauge, &
                 (1.d0,0.d0), dvpsi(1, 1), npwx*npol, &
                              psi_gauge, ahc_nbnd_gauge, &
                 (0.d0,0.d0), dvpsi_gauged, npwx*npol)
      !
      ! Rotate dvpsi from pattern basis (u) to Cartesian basis
      !
      DO jmode = 1, nmodes
         CALL ZAXPY(npwx*npol*ahc_nbnd, CONJG(u(jmode, imode)), &
                    dvpsi_gauged, 1, dvpsi_cart(:, :, jmode), 1)
      ENDDO
      !
    ENDDO ! imode
    !
    DEALLOCATE(dvpsi_gauged)
    !
    CALL ahc_do_gkk(ik)
    !
    IF (.NOT. skip_upperfan) THEN
      CALL ahc_do_upperfan(ik)
    ENDIF
    !
    ! If q = Gamma, compute Debye-Waller matrix elements
    !
    IF (lgamma) CALL ahc_do_dw(ik)
    !
  ENDDO ! ik
  !
  WRITE(stdout, '(5x,a)') "AHC e-ph calculation done"
  !
  DEALLOCATE(ahc_gkk)
  DEALLOCATE(ahc_upfan)
  DEALLOCATE(ahc_dw)
  DEALLOCATE(dvpsi_cart)
  DEALLOCATE(psi_gauge)
  !
  ! Close output file units
  !
  IF (me_pool == root_pool) THEN
    CLOSE(iungkk, STATUS='KEEP')
    IF (.NOT. skip_upperfan) CLOSE(iunupfan, STATUS='KEEP')
    IF (lgamma) CLOSE(iundw, STATUS='KEEP')
  ENDIF ! root_pool
  !
  CALL stop_clock('ahc_elph')
  !
END SUBROUTINE elph_do_ahc
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
END MODULE ahc
!------------------------------------------------------------------------------
