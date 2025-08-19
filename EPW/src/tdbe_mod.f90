  !
  ! Copyright (C) 2016-2023 EPW-Collaboration
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !-------------------------------------------------------------------------
  MODULE tdbe_mod
  !-------------------------------------------------------------------------
  !!
  !! 
  !! This module contains all the subroutines for real-time simulation with
  !! time-dependent Boltzmann Equation. 
  !! Authored by Yiming Pan (pan@physik.uni-kiel.de) 
  !! and Fabio Caruso (caruso@physik.uni-kiel.de) (Apr. 2025)
  !! wanninterp is adapted from use_wannier, to perform wan2bloch of electron 
  !! and phonon band structure. 
  !! read_g2_tdbe is adpated from read_eigenvalues, read_ephmat and read_kqmap
  !! (see io/io_supercond.f90), to read g2 into a grid defined by nkf1d, 
  !! nkf2d and nkf3d.
  !!
  !
  IMPLICIT NONE
  !
  CONTAINS
    !
    !-----------------------------------------------------------------------
    SUBROUTINE double_grid_map()
    !-----------------------------------------------------------------------
    !! This subroutine map k-point defined on a fine mesh with nkf1d, nkf2d, 
    !! and nkf3d to a finer grid defined with nkf1, nkf2, nkf3
    !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE input,         ONLY : nkf1, nkf2, nkf3, nqf1, nqf2, nqf3, &
                              nkf1d, nkf2d, nkf3d, nqf1d, nqf2d,  &
                              nqf3d
    USE tdbe_common,   ONLY : indx_mapk, indx_mapq  
    !
    IMPLICIT NONE
    !
    INTEGER :: ierr
    !! Error info
    INTEGER :: i
    !! index of k point 
    INTEGER :: j
    !! index of k point 
    INTEGER :: k
    !! index of k point 
    INTEGER :: rs
    !! index of k point 
    INTEGER :: rss
    !! variable saving the result of quotient
    INTEGER :: pk(3)
    !! nkf1 / nkf1d
    INTEGER :: pq(3)
    !! nqf1 / nqf1d
    INTEGER :: ii
    !! Index of direction b1
    INTEGER :: jj 
    !! Index of direction b2
    INTEGER :: kk
    !! Index of direction b3
    INTEGER :: nktotff
    !! nkf1 * nkf2 * nkf3
    INTEGER :: nqtotff
    !! nqf1 * nqf2 * nqf3
    INTEGER :: ik
    !! k index
    INTEGER :: iq
    !! q index
    !
    nktotff = nkf1 * nkf2 * nkf3
    !
    nqtotff = nqf1 * nqf2 * nqf3
    !
    pk(1) = INT(nkf1 / nkf1d)
    IF (pk(1) * nkf1d /= nkf1) CALL errore('double_grid_map','nkf1 should be a multiple of nkf1d', 1)
    pk(2) = INT(nkf2 / nkf2d)
    IF (pk(2) * nkf2d /= nkf2) CALL errore('double_grid_map','nkf2 should be a multiple of nkf2d', 1)
    pk(3) = INT(nkf3 / nkf3d)
    IF (pk(3) * nkf3d /= nkf3) CALL errore('double_grid_map','nkf3 should be a multiple of nkf3d', 1)
    pq(1) = INT(nqf1 / nqf1d)
    IF (pq(1) * nqf1d /= nqf1) CALL errore('double_grid_map', 'nqf1 should be a multiple of nqf1d', 1)
    pq(2) = INT(nqf2 / nqf2d)
    IF (pq(2) * nqf2d /= nqf2) CALL errore('double_grid_map', 'nqf2 should be a multiple of nqf2d', 1)
    pq(3) = INT(nqf3 / nqf3d)
    IF (pq(3) * nqf3d /= nqf3) CALL errore('double_grid_map', 'nqf3 should be a multiple of nqf3d', 1)
    !
    ALLOCATE(indx_mapk(nktotff), STAT = ierr)
    IF (ierr /= 0) CALL errore('double_grid_map', 'Error allocating indx_mapk', 1)
    ALLOCATE(indx_mapq(nqtotff), STAT = ierr)
    IF (ierr /= 0) CALL errore('double_grid_map', 'Error allocating indx_mapq', 1)
    indx_mapk(:) = 0
    indx_mapq(:) = 0
    DO ik = 1, nktotff
      k = 1 + MODULO(ik - 1, nkf3)
      rs = INT((ik - 1) / nkf3)
      j = 1 + MODULO(rs, nkf2)
      rss = INT(rs / nkf2)
      i = 1 + MODULO(rss, nkf1)
      ii = indxmap(i, pk(1), nkf1d)
      jj = indxmap(j, pk(2), nkf2d)
      kk = indxmap(k, pk(3), nkf3d)
      indx_mapk(ik) = (ii - 1) * nkf2d * nkf3d + (jj - 1) * nkf3d + kk
    ENDDO
    !
    DO iq = 1, nqtotff
      k = 1 + MODULO(iq - 1, nqf3)
      rs = INT((iq - 1) / nqf3)
      j = 1 + MODULO(rs, nqf2)
      rss = INT(rs / nqf2)
      i = 1 + MODULO(rss, nqf1)
      ii = indxmap(i, pq(1), nqf1d)
      jj = indxmap(j, pq(2), nqf2d)
      kk = indxmap(k, pq(3), nqf3d)
      indx_mapq(iq) = (ii - 1) * nqf2d * nqf3d + (jj - 1) * nqf3d + kk
    ENDDO
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE double_grid_map 
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    FUNCTION indxmap(inxf,p,nkfd)
    !-----------------------------------------------------------------------
    !!
    !-----------------------------------------------------------------------
    IMPLICIT NONE
    !
    INTEGER :: indxmap
    !! index on 
    INTEGER, INTENT(IN) :: inxf
    !! 
    INTEGER, INTENT(IN) :: p
    !!
    INTEGER, INTENT(IN) :: nkfd
    !!
    !
    indxmap = INT(DBLE((inxf - 1)) / DBLE(p) + 0.5) + 1
    IF (indxmap > nkfd) THEN
        indxmap = indxmap - nkfd
    ENDIF
    RETURN
    !
    !-----------------------------------------------------------------------
    END FUNCTION indxmap
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE map2nkfs() 
    !-----------------------------------------------------------------------
    !!
    !-----------------------------------------------------------------------
    USE kinds,             ONLY : DP
    USE input,             ONLY : nkf1, nkf2, nkf3, mp_mesh_k
    USE supercond_common,  ONLY : ixkff, xkfs, nkfs
    USE tdbe_common,       ONLY : ikfsdx, nkfsden, indx_mapk
    USE io_global,         ONLY : stdout, ionode_id, meta_ionode
    USE mp_global,         ONLY : inter_pool_comm,npool, my_pool_id
    USE mp,                ONLY : mp_bcast, mp_barrier, mp_sum, mp_max, mp_min
    USE parallelism,       ONLY : fkbounds
    !
    IMPLICIT NONE
    !
    INTEGER :: ierr
    !! Error info
    INTEGER :: nktotff
    !! nkf1 * nkf2 * nkf3
    INTEGER :: ik
    !! k point index
    INTEGER :: ipool
    !! index of pool
    INTEGER :: lower_bound
    !! lower bound of k point parallelization
    INTEGER :: upper_bound
    !! upper bound of k point parallelization
    INTEGER :: fermicount_dense
    !! number of k points on dense grid in fermi window
    INTEGER :: nkfsden_all(npool)
    !! k point counted in each pool
    !
    CALL fkbounds(nkfs, lower_bound, upper_bound)
    IF (mp_mesh_k) THEN
      WRITE(stdout, '(/5x,a,i8,a/)') 'scatter nkfs (', nkfs ,') irr- k-points within energy window into cpus'
    ELSE
      WRITE(stdout, '(/5x,a,i8,a/)') 'scatter nkfs (', nkfs ,') k points within energy window into cpus'
    ENDIF
    nkfsden = 0
    nkfsden_all(:) = 0
    nktotff = nkf1 * nkf2 * nkf3

    DO ik = 1, nktotff
      IF (ixkff(indx_mapk(ik)) >= lower_bound .AND. ixkff(indx_mapk(ik)) <= upper_bound) THEN
         nkfsden = nkfsden + 1
      ENDIF
    ENDDO
    nkfsden_all(my_pool_id + 1) = nkfsden
    CALL mp_sum(nkfsden_all, inter_pool_comm)
    WRITE(stdout, '(5x, a)') '     Number of kpoints of the fine mesh within energy window determined by nkfs :'
    WRITE(stdout, *) nkfsden_all(:)
    WRITE(stdout,'(/5x, a/)') REPEAT('=',67)
    ALLOCATE(ikfsdx(nkfsden), STAT = ierr)
    IF (ierr /= 0) CALL errore('map2nkfs', 'Error allocating ikfsdx', 1)
    ikfsdx(:) = 0
    fermicount_dense = 0
    DO ik = 1, nktotff
      IF (ixkff(indx_mapk(ik)) >= lower_bound .AND. ixkff(indx_mapk(ik)) <= upper_bound) THEN
         fermicount_dense = fermicount_dense + 1
         ikfsdx(fermicount_dense) = ik
      ENDIF
    ENDDO
    RETURN
    !-----------------------------------------------------------------------
    END SUBROUTINE map2nkfs 
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE qwindow_tdbe(exst)
    !-----------------------------------------------------------------------
    !!
    !! This routine pre-computes the q-points that falls within the fstichk.
    !! If at least 1 k-point is such that at least one k+q eigenenergy falls
    !! within the user-defined fstichk, then the q-point is taken.
    !!
    !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE global_var,         ONLY : nqf, xqf, xkf, nkf, nqtotf, nktotf, &
                              nbndfst
    USE tdbe_common,       ONLY : enk_all, totq, selecq  
    USE io_global,     ONLY : ionode_id, stdout
    USE io_var,        ONLY : iunselecq
    USE mp_global,     ONLY : npool, world_comm, my_pool_id
    USE mp,            ONLY : mp_sum, mp_bcast
    USE bzgrid,        ONLY : kpmq_map
!    USE ep_constants, ONLY : twopi, ci, zero, eps6, ryd2ev, czero
    USE input,         ONLY : fsthick, mp_mesh_k, restart_step
    USE pwcom,         ONLY : ef
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(in) :: exst
    !! Local variable
    INTEGER :: iq
    !! Counter on coarse q-point grid
    INTEGER :: ik, ikk
    !! Counter on coarse k-point grid
    INTEGER :: found(npool)
    !! Indicate if a q-point was found within the window
    INTEGER :: ind1
    !! Index of the k point from the full grid.
    INTEGER :: ind2
    !! Index of the k+q point from the full grid.
    INTEGER :: ierr, ios
    !! Error status
    INTEGER :: nqtot
    !! total nb of q read from fmt file
    REAL(KIND = DP) :: xxq(3)
    !! Current q-point
    REAL(KIND = DP) :: xkk(3)
    !! Current k-point on the fine grid
    REAL(KIND = DP) :: xkq(3)
    !! Current k-point on the fine grid
    REAL(KIND = DP) :: gammapoint(3) = (/0.0, 0.0, 0.0/)
    !
    IF(exst) THEN
      WRITE(stdout, '(5x, a)') 'Read seleced q points from selecq.fmt'
      IF (my_pool_id == ionode_id) THEN
        OPEN(UNIT = iunselecq, FILE = 'selecq_tdbe.fmt', STATUS = 'old', IOSTAT = ios)
        READ(iunselecq,*) totq
        ALLOCATE(selecq(totq), STAT = ierr)
        IF (ierr /= 0) CALL errore('qwindow_tdbe', 'Error allocating selecq', 1)
        selecq(:) = 0
        READ(iunselecq,*) nqtot
        READ(iunselecq,*) selecq(:)
        CLOSE(iunselecq)
      ENDIF
      CALL mp_bcast(totq, ionode_id, world_comm)
      IF (my_pool_id /= ionode_id) ALLOCATE(selecq(totq))
      CALL mp_bcast(nqtot , ionode_id, world_comm)
      CALL mp_bcast(selecq, ionode_id, world_comm)
      IF (nqtot /= nqtotf) THEN
        CALL errore('qwindow_tdbe', 'Cannot read from selecq_tdbe.fmt', 1 )
      ENDIF
      !      
    ELSE
      ALLOCATE(selecq(nqf), STAT = ierr)
      IF (ierr /= 0) CALL errore('qwindow_tdbe', 'Error allocating selecq', 1)
      selecq(:) = 0
      !
      DO iq = 1, nqf
        xxq = xqf(:, iq)
        !
        found(:) = 0
        DO ik = 1, nkf
          ikk = 2 * ik - 1
          xkk = xkf(:, ikk)
          xkq = xkk + xxq
          !
          !CALL kpmq_map(xkk, (/0d0,0d0,0d0/), 1, ind1)
          CALL kpmq_map(xkk, gammapoint, 1, ind1)
          CALL kpmq_map(xkk, xxq, 1, ind2)
          !
          IF ((MINVAL(ABS(enk_all(:, ind1) - ef)) < fsthick) .AND. &
              (MINVAL(ABS(enk_all(:, ind2) - ef)) < fsthick)) THEN
              found(my_pool_id + 1) = 1
              EXIT ! exit the loop
          ENDIF

        ENDDO ! k-loop
        ! If found on any k-point from the pools
        CALL mp_sum(found, world_comm)
        !
        IF (SUM(found) > 0) THEN
          totq = totq + 1
          selecq(totq) = iq
          !
          IF (MODULO(totq, restart_step) == 0) THEN
            WRITE(stdout, '(5x,a,i15,i15)')'Number selected, total', totq, iq
          ENDIF
        ENDIF
      ENDDO ! iq
      IF (my_pool_id == ionode_id) THEN
        OPEN(UNIT = iunselecq, FILE = 'selecq_tdbe.fmt', ACTION = 'write')
        WRITE(iunselecq, *) totq    ! Selected number of q-points
        WRITE(iunselecq, *) nqtotf  ! Total number of q-points
        WRITE(iunselecq, *) selecq(1:totq)
        CLOSE(iunselecq)
      ENDIF
    ENDIF  
    !-----------------------------------------------------------------------
    END SUBROUTINE qwindow_tdbe
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE read_g2_tdbe()
    !-----------------------------------------------------------------------
    !
    !! read the frequencies, energy, e-ph matrix element
    !! obtained from a previous epw run
    !
    USE kinds,                ONLY : DP
    USE io_global,            ONLY : stdout, ionode_id
    USE io_var,               ONLY : iufilfreq, iufilegnv, iuselecqfd,      &
                                     iufilikmap, iufileph, iurpa
    USE io_files,             ONLY : prefix, tmp_dir
    USE input,                ONLY : nkf1, nkf2, nkf3, nqf1, nqf2, nqf3,    &
                                     nqf1d, nqf2d, nqf3d, fsthick, nkf1d,   &
                                     nkf2d, nkf3d, degaussw, eps_acoustic                     
    USE ep_constants,         ONLY : zero, ryd2ev, eps8, one, two, czero
    USE mp_global,            ONLY : inter_pool_comm, npool, my_pool_id
    USE mp,                   ONLY : mp_bcast, mp_barrier, mp_sum
    USE supercond_common,     ONLY : dosef,efm0 => ef0, ekfs, wkfs,ixkf, g2, &
                                     ixkff, xkff, xkfs, nkfs, ixkqf, ixqfs,  &
                                     nbndfs, nqfs, memlt_pool, ekfs_all,     &
                                     wkfs_all, xkfs_all, nkfs_all, ixkf_inv, &
                                     nbndfs_all, ibnd_kfs_to_kfs_all,        &
                                     ibnd_kfs_all_to_kfs
    USE tdbe_common,          ONLY : iq2nqfs
    USE parallelism,          ONLY : fkbounds
    USE low_lvl,              ONLY : set_ndnmbr, mem_size_eliashberg
    USE input,                ONLY : mp_mesh_k, dg_tdbe, ephmat_dir,     &
                                     lscreen_tdbe, scr_typ 
    USE global_var,           ONLY : epsi, wf
    USE screening,            ONLY : rpa_epsilon, tf_epsilon
    USE modes,                ONLY : nmodes
    USE control_flags,        ONLY : iverbosity
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = 256) :: filfreq, filegnv, filikmap, filephmat
    !! file name
    CHARACTER(LEN = 256) :: dirname
    !! Name of the directory where ikmap egnv freq ephmat files are saved
    CHARACTER(LEN = 4) :: filelab
    !! File name
    INTEGER :: iq, iqd
    !! Counter on q points
    INTEGER :: imode
    !! Counter on modes
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: ibnd, jbnd
    !! Counter on bands
    INTEGER :: ikfs
    !! Counter on k-points for ekfs
    INTEGER :: ibndfs
    !! Counter on bands for ekfs
    INTEGER :: iqq
    !! Q-point index from selecq.fmt window
    INTEGER :: totqfd
    !! Total number of q-points inside fsthick
    INTEGER :: ipool 
    !! Counter on pools
    INTEGER :: npool_g2
    !! Number of pools used to calculate g2
    INTEGER :: tmp_pool_id
    !! Pool index read from file
    INTEGER :: nqtotfd
    !! Number of q points
    INTEGER :: nkftotfd
    !! Number of k points
    !INTEGER :: nmode
    !! number of phonon branches
    INTEGER :: i, j, k, nk, n, rs, rss
    !! Counter on k points
    INTEGER :: ii, jj, kk
    !! Counter on q points
    INTEGER :: nkq
    !! Index of k+sign*q on the fine k-mesh
    INTEGER :: lower_bnd, upper_bnd
    !! Lower/upper bound index after k parallelization
    INTEGER :: nks
    !! Number of non-equivalent k points
    !! Counter on k points within the Fermi shell
    INTEGER :: ns
    !! Counter on rotation operations
    INTEGER(8) :: imelt
    !! Memory allocated
    INTEGER, ALLOCATABLE :: index_(:, :)
    !! Index of q-point on the full q-mesh for which k+sign*q is within the Fermi shell
    INTEGER :: nqf1d_, nqf2d_, nqf3d_
    !! Temporary variable for number of q-points along each direction
    INTEGER :: nkf1d_, nkf2d_, nkf3d_
    !! Temporary variable for number of k-points along each direction
    INTEGER, ALLOCATABLE :: selecqfd(:)    
    !! List of selected q-points
    INTEGER :: m
    !! Band indexes
    INTEGER :: ios
    !! IO error message
    INTEGER :: ierr
    !! Error status
    REAL(KIND = DP) :: ebnd
    !! Local variable for energy
    REAL(KIND = DP) :: ebndmax1
    !! temporarily store the maximum value of ABS(ekfs-ekfs_all)
    REAL(KIND = DP) :: ebndmax2
    !! temporarily store the maximum value of ABS(ekfs-ekfs_all)
    REAL(KIND = DP), ALLOCATABLE :: wfd(:, :) 
    !! phonon frequency and q-points on the double mesh
    REAL(KIND = DP), ALLOCATABLE :: wqfd(:)
    !! weight of q-points read from filfreq
    REAL(KIND = DP), ALLOCATABLE :: xqfd(:, :)
    !! coordinates of q-points read from filfreq
    REAL(KIND = DP) :: degaussw0
    !! smearing parameter read from filegnv
    REAL(KIND = DP) ::  fsthick0
    !! Fermi thickness read from filegnv
    REAL(KIND = DP) ::  eftmp
    !! Fermi level read from filegnv
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! The derivative of wgauss:  an approximation to the delta function
    REAL(KIND = DP) :: xk(3)
    !! coordinates of k points
    REAL(KIND = DP) :: xq(3), qcart(3)
    !! coordinates of q points
    INTEGER :: nmin, nmax
    !! Lower/upper bound index for .ephmat file read in current pool
!    INTEGER :: nmin_all(npool), nmax_all(npool)
!    !! Lower/upper bound index for .ephmat file read in current pool
    INTEGER :: nnk
    !! Number of k-points within the Fermi shell
    INTEGER, ALLOCATABLE ::nnk_all(:)
    !! Number of k-points within the Fermi shell
    INTEGER, ALLOCATABLE :: nnq(:)
    !! Number of k+q points within the Fermi shell for a given k-point
    INTEGER, ALLOCATABLE :: nkpool(:)
    !! nkpool(ipool) - sum of nr. of k points from pool 1 to pool ipool
    REAL(KIND = DP) :: gmat
    !! Electron-phonon matrix element square
    INTEGER :: ikdum
    !! Dummy
    INTEGER :: nkf_mesh
    !! Total number of k points
    COMPLEX(KIND = DP), ALLOCATABLE :: eps_rpa(:)
    !
    IF (my_pool_id == ionode_id) THEN
      OPEN(UNIT = iuselecqfd, FILE = 'selecq.fmt', STATUS = 'old', IOSTAT = ios)
      ! ENDIF
      READ(iuselecqfd, *) totqfd
      ALLOCATE(selecqfd(totqfd), STAT = ierr)
      IF (ierr /= 0) CALL errore('read_g2_tdbe', 'Error allocating selecqfd', 1)
      selecqfd(:) = 0
      READ(iuselecqfd, *) nqtotfd
      IF (nqtotfd /= nqf1d * nqf2d * nqf3d) &
        CALL errore('read_g2_tdbe', 'selecq.fmt is not calculated on the nqf1d, nqf2d, nqf3d mesh', 1)
      READ(iuselecqfd, *) selecqfd(:)
      CLOSE(iuselecqfd)
      !
      ! read frequencies from file
      dirname = TRIM(ephmat_dir ) // TRIM(prefix) // '.ephmat'
      npool_g2=0
      CALL SYSTEM('ls '//TRIM(dirname)//'| grep ephmat | wc -l > npool_g2')
      OPEN(123,file='npool_g2',IOSTAT = ios)
      IF (ios /=0) CALL errore('read_g2_tdbe','error executing a shell command',1)
      read(123,*) npool_g2
      IF (npool_g2 == 0) CALL errore('read_g2_tdbe','npool_g2 can not be zero',1)
      CLOSE(123)
      WRITE(stdout,*) 'npool_g2 = ', npool_g2
      call SYSTEM('rm -rf npool_g2')
    ENDIF
    CALL mp_bcast(totqfd, ionode_id, inter_pool_comm)
    IF (my_pool_id /= ionode_id) ALLOCATE(selecqfd(totqfd))
    CALL mp_bcast(selecqfd, ionode_id, inter_pool_comm)
    nqtotfd = nqf1d * nqf2d  * nqf3d
    CALL mp_bcast(npool_g2, ionode_id, inter_pool_comm)
    ALLOCATE(nnk_all(npool_g2), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_g2_tdbe', 'Error allocating nnk_all', 1)
    ALLOCATE(nkpool(npool_g2), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_g2_tdbe', 'Error allocating nkpool', 1)
    CALL mp_barrier(inter_pool_comm)
    !
    ALLOCATE(wfd(nmodes, nqtotfd), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_g2_tdbe', 'Error allocating wfd', 1)
    ALLOCATE(wqfd(nqtotfd), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_g2_tdbe', 'Error allocating wqfd', 1)
    ALLOCATE(xqfd(3, nqtotfd), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_g2_tdbe', 'Error allocating xqfd', 1)
    wfd(:, :) = zero
    wqfd(:) = 1.d0 / DBLE(nqtotfd)
    xqfd(:, :) = zero
    !
    !IF (my_pool_id == ionode_id) THEN
    DO iqq = 1, totqfd ! loop over q-points in fsthick
      iq = selecqfd(iqq)
      !
      k = 1 + MODULO(iq - 1, nqf3d)
      rs = INT((iq - 1) / nqf3d) 
      j = 1 + MODULO(rs, nqf2d)
      rss = INT(rs / nqf2d)
      i = 1 + MODULO(rss, nqf1d)
      xqfd(1, iq) = DBLE(i - 1) / DBLE(nqf1d)
      xqfd(2, iq) = DBLE(j - 1) / DBLE(nqf2d)
      xqfd(3, iq) = DBLE(k - 1) / DBLE(nqf3d)
      ii = INT(xqfd(1, iq) * nqf1) + 1
      jj = INT(xqfd(2, iq) * nqf2) + 1
      kk = INT(xqfd(3, iq) * nqf3) + 1
      iqd = (ii - 1) * nkf2 * nkf3 + (jj - 1) * nkf3 + kk
      DO imode = 1, nmodes
        wfd(imode, iq) = wf(imode, iqd)
      ENDDO
    ENDDO
    CALL mp_barrier(inter_pool_comm)
    WRITE(stdout,'(/5x,a/)') 'Finished reading freq file'
    !! Now start reading eigenvalues of states in the fermi window
    IF (my_pool_id == ionode_id) THEN
      !
      ! SP: Needs to be initialized
      nbndfs = 0
      nkfs = 0
      !
      ! read eigenvalues on the irreducible fine k-mesh
      !
      dirname = TRIM(ephmat_dir ) // TRIM(prefix) // '.ephmat'
      filegnv = TRIM(dirname) // '/' // 'egnv'
      OPEN(UNIT = iufilegnv, FILE = filegnv, STATUS = 'unknown', &
           FORM = 'unformatted', ACCESS = 'stream', IOSTAT = ios)
      IF (ios /= 0) CALL errore('read_g2_tdbe', 'error opening file '//filegnv, iufilegnv)
      READ(iufilegnv) nkftotfd, nkf1d_, nkf2d_, nkf3d_, nkfs
      READ(iufilegnv) nbndfs, eftmp, efm0, dosef, degaussw0,fsthick0
      IF (nkf1d /= nkf1d_ .OR. nkf2d /= nkf2d_ .OR. nkf3d /= nkf3d) &
        CALL errore('read_g2_tdbe', 'e-ph mat elements were not calculated on the nkf1d, nkf2d, nkf3d mesh', 1)
      !
      !!! Since in real time carrier dynamic simulation, the user usually manually
      !!! defined a fermi level, given by efm0, one needs to read the exact fermi energy
      !!! calculated from nscf calculation as ef
      !
      WRITE(stdout, '(5x, a32, ES20.10)') 'Fermi level defined by user(eV) = ', efm0 * ryd2ev
      WRITE(stdout, '(5x, a, ES20.10)') 'Fermi level ef  written in the previous calculation (eV) = ', eftmp * ryd2ev
      WRITE(stdout, '(5x, a32, ES20.10)') 'DOS(states/spin/eV/Unit Cell) = ', dosef / ryd2ev
      WRITE(stdout, '(5x, a32, ES20.10)') 'Electron smearing (eV) = ', degaussw * ryd2ev
      WRITE(stdout,'(5x,a32,ES20.10)') 'Fermi window (eV) = ', fsthick0 * ryd2ev
      WRITE(stdout, '(5x, i7, a/)') nbndfs, ' bands within the Fermi window'
      IF (mp_mesh_k) THEN
        WRITE(stdout, '(5x, a, i9, a, i9/)') 'Nr irreducible k-points within the Fermi shell = ', nkfs, ' out of ', nkftotfd
      ELSE
        WRITE(stdout, '(5x, a, i9, a, i9/)') 'Nr k-points within the Fermi shell = ', nkfs, ' out of ', nkftotfd
      ENDIF
    ENDIF
    !
    ! first node broadcasts everything to all nodes
    CALL mp_bcast(nkf1d, ionode_id, inter_pool_comm)
    CALL mp_bcast(nkf2d, ionode_id, inter_pool_comm)
    CALL mp_bcast(nkf3d, ionode_id, inter_pool_comm)
    CALL mp_bcast(nkfs, ionode_id, inter_pool_comm)
    CALL mp_bcast(nbndfs, ionode_id, inter_pool_comm)
    CALL mp_bcast(efm0, ionode_id, inter_pool_comm)
    CALL mp_bcast(dosef, ionode_id, inter_pool_comm)
    CALL mp_bcast(fsthick0, ionode_id, inter_pool_comm)
    !
    !
    ALLOCATE(wkfs(nkfs), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_g2_tdbe', 'Error allocating wkfs', 1)
    ALLOCATE(xkfs(3, nkfs), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_g2_tdbe', 'Error allocating xkfs', 1)
    wkfs(:) = zero
    xkfs(:, :) = zero
    !
    IF (my_pool_id == ionode_id) THEN
      ! HM: nbndfs_all should be same with nbndfst.
      nbndfs_all = nbndfs
      nkfs_all = nkftotfd
    ENDIF
    CALL mp_bcast(nbndfs_all, ionode_id, inter_pool_comm)
    CALL mp_bcast(nkfs_all, ionode_id, inter_pool_comm)
    ALLOCATE(ekfs_all(nbndfs_all, nkfs_all), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_g2_tdbe', 'Error allocating ekfs_all', 1)
    ALLOCATE(wkfs_all(nkfs_all), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_g2_tdbe', 'Error allocating wkfs_all', 1)
    ALLOCATE(xkfs_all(3, nkfs_all), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_g2_tdbe', 'Error allocating xkfs_all', 1)
    xkfs_all(:, :) = zero
    ! sanity choice
    wkfs_all(:) = zero
    ekfs_all(:, :) = zero
    ALLOCATE(ixkf(nkfs_all), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_g2_tdbe', 'Error allocating ixkf', 1)
    ALLOCATE(ixkf_inv(nkfs), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_g2_tdbe', 'Error allocating ixkf_inv', 1)
    ixkf(:) = 0
    ixkf_inv(:) = 0
    !
    IF (my_pool_id == ionode_id) THEN
      !
      ! nbndfs - nr of bands within the Fermi shell
      !
      nbndfs = 0
      DO ik = 1, nkfs_all ! loop over irreducible k-points
        READ(iufilegnv) wkfs_all(ik), xkfs_all(:, ik), ikdum, ixkf(ik)
        IF (ik /= ikdum) CALL errore('read_eigenvalues', 'error reading file '//filegnv, ik)
        IF (ixkf(ik) /= 0) THEN
          ixkf_inv(ixkf(ik)) = ik
        ENDIF
        n = 0
        DO ibnd = 1, nbndfs_all
          READ(iufilegnv) ekfs_all(ibnd, ik)
          ! ! go from Ryd to eV
          ! ekfs_all(ibnd, ik) = ekfs_all(ibnd, ik) * ryd2ev
          IF (ABS(ekfs_all(ibnd, ik) - efm0) < fsthick) THEN
            n = n + 1
            IF (nbndfs < n) nbndfs = n
          ENDIF
        ENDDO
      ENDDO
      !
      WRITE(stdout, '(5x, i7, a/)') nbndfs, ' bands within the Fermi window'
      CLOSE(iufilegnv)
      !
    ENDIF
    ! first node broadcasts everything to all nodes
    CALL mp_bcast(nbndfs, ionode_id, inter_pool_comm)
    CALL mp_bcast(ekfs_all, ionode_id, inter_pool_comm)
    CALL mp_bcast(wkfs_all, ionode_id, inter_pool_comm)
    CALL mp_bcast(xkfs_all, ionode_id, inter_pool_comm)
    CALL mp_bcast(ixkf_inv, ionode_id, inter_pool_comm)
    CALL mp_bcast(ixkf, ionode_id, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    !  
    ALLOCATE(ekfs(nbndfs, nkfs), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_g2_tdbe', 'Error allocating ekfs', 1)
    ekfs(:, :) = efm0 - 10.d0 * fsthick
    !
    ALLOCATE(ibnd_kfs_all_to_kfs(nbndfs_all, nkfs), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_g2_tdbe', 'Error allocating ibnd_kfs_all_to_kfs', 1)
    ALLOCATE(ibnd_kfs_to_kfs_all(nbndfs, nkfs), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_g2_tdbe', 'Error allocating ibnd_kfs_to_kfs_all', 1)
    ibnd_kfs_all_to_kfs(:, :) = 0
    ibnd_kfs_to_kfs_all(:, :) = 0
    !
    !
    IF (my_pool_id == ionode_id) THEN
      DO ik = 1, nkfs_all
        IF (ixkf(ik) == 0) CYCLE
        ikfs = ixkf(ik)
        wkfs(ikfs) = wkfs_all(ik)
        xkfs(:, ikfs) = xkfs_all(:, ik)
        n = 0
        DO ibnd = 1, nbndfs_all
          IF (ABS(ekfs_all(ibnd, ik) - efm0) < fsthick) THEN
            n = n + 1
            ekfs(n, ikfs) = ekfs_all(ibnd, ik)
            ibnd_kfs_all_to_kfs(ibnd, ikfs) = n
            ibnd_kfs_to_kfs_all(n, ikfs) = ibnd
          ENDIF
        ENDDO
      ENDDO
    ENDIF
    !
    ! first node broadcasts everything to all nodes
    CALL mp_bcast(wkfs, ionode_id, inter_pool_comm)
    CALL mp_bcast(xkfs, ionode_id, inter_pool_comm)
    CALL mp_bcast(ekfs, ionode_id, inter_pool_comm)
    CALL mp_bcast(ibnd_kfs_all_to_kfs, ionode_id, inter_pool_comm)
    CALL mp_bcast(ibnd_kfs_to_kfs_all, ionode_id, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    !
    WRITE(stdout,'(/5x,a/)') 'Finished reading egnv file '
    !
    ebndmax1 = 0.0d0
    ebndmax2 = 0.0d0
    DO ik = 1, nkfs_all
      IF (ixkf(ik) == 0) CYCLE
      ikfs = ixkf(ik)
      DO ibnd = 1, nbndfs_all
        IF (ABS(ekfs_all(ibnd, ik) - efm0) < fsthick) THEN
          ibndfs = ibnd_kfs_all_to_kfs(ibnd, ikfs)
          ebnd = ekfs(ibndfs, ikfs) - ekfs_all(ibnd, ik)
          ebnd = ABS(ebnd)
          ebndmax1 = MAX(ebnd, ebndmax1)
          IF (iverbosity == 5) THEN
            WRITE(stdout, '(5x, f15.7, 1x, f15.7)') &
                  ekfs(ibndfs, ikfs), ekfs_all(ibnd, ik)
          ENDIF
        ENDIF
      ENDDO
    ENDDO
    !
    DO ikfs = 1, nkfs
      ik = ixkf_inv(ikfs)
      DO ibndfs = 1, nbndfs
        IF (ABS(ekfs(ibndfs, ikfs) - efm0) < fsthick) THEN
          ibnd = ibnd_kfs_to_kfs_all(ibndfs, ikfs)
          ebnd = ekfs(ibndfs, ikfs) - ekfs_all(ibnd, ik)
          ebnd = ABS(ebnd)
          ebndmax2 = MAX(ebnd, ebndmax2)
          !!!!!!!! FOR DEBUG !!!!!!!!
          !WRITE(stdout, '(5x, f15.7, 1x, f15.7)') &
          !ekfs(ibndfs, ikfs), ekfs_all(ibnd, ik)
          !!!!!!!! FOR DEBUG !!!!!!!!
        ENDIF
      ENDDO
    ENDDO
    !
    IF ((ebndmax1 > eps8) .OR. (ebndmax2 > eps8)) THEN
      CALL errore('read_g2_tdbe', 'ekfs_all is not equal to ekfs.', 1)
    ENDIF
    !
    WRITE(stdout,'(/5x,a/)') 'Start reading e-ph matrix elements '
    !
    IF (lscreen_tdbe .AND. scr_typ <= 1 ) THEN
      ALLOCATE(eps_rpa(nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('read_g2_tdbe', 'Error allocating eps_rpa', 1)
      eps_rpa(:)           = czero
    ENDIF
    !
    !! Write dielectric function
    IF (lscreen_tdbe .AND. scr_typ <= 1 .AND. my_pool_id == ionode_id) THEN
      OPEN(UNIT = iurpa, FILE = 'eps_rpa.dat')
        DO iq = 1, nqtotfd ! loop over q-points 
          xq = xqfd(:, iq)
          CALL shift_q(xq, qcart)  ! shift q point to the 1st BZ
          IF (scr_typ == 0) CALL rpa_epsilon(xq, wfd(:, iq), nmodes, epsi, eps_rpa)
          IF (scr_typ == 1) CALL tf_epsilon(xq, nmodes, epsi, eps_rpa)
          WRITE(iurpa, *) qcart(:)
          DO imode = 1, nmodes
            WRITE(iurpa, '(3E22.14)') wfd(imode, iq), eps_rpa(imode)
          ENDDO
        ENDDO
      CLOSE(iurpa)
    ENDIF
    ! 
    ALLOCATE(memlt_pool(npool), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_g2_tdbe', 'Error allocating memlt_pool', 1)
    memlt_pool(:) = zero
    !
    ! get the size of arrays for frequency and eigenvalue variables allocated in
    ! read_frequencies and read_eigenvalues
    imelt = (nmodes + 4) * nqtotfd  + (4 + 2 * nbndfs) * nkfs
    CALL mem_size_eliashberg(2, imelt)
    !
    nkftotfd = nkf1d * nkf2d * nkf3d
    !
    ALLOCATE(ixkff(nkftotfd), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_g2_tdbe', 'Error allocating ixkff', 1)
    ixkff(:) = 0
    !
    IF (my_pool_id == ionode_id) THEN
      !
      dirname = TRIM(ephmat_dir) // TRIM(prefix) // '.ephmat'
      filikmap = TRIM(dirname) // '/' // 'ikmap'
      !OPEN(UNIT = iufilikmap, FILE = filikmap, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
      OPEN(UNIT = iufilikmap, FILE = filikmap, STATUS = 'unknown', &
           FORM = 'unformatted', ACCESS = 'stream', IOSTAT = ios)
      IF (ios /= 0) CALL errore('read_g2_tdbe', 'error opening file ' // filikmap, iufilikmap)
      !
      !READ(iufilikmap, *) ixkff(1:nkftot)
      READ(iufilikmap) ixkff(1: nkftotfd)
      !
      CLOSE(iufilikmap)
    ENDIF
    CALL mp_bcast(ixkff, ionode_id, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    !
    CALL fkbounds(nkfs, lower_bnd, upper_bnd)
    ! get the size of required memory for ixkqf, nqfs, index_
    imelt = (2*nqtotfd + 1) * nkfs + (upper_bnd - lower_bnd + 1) * nqtotfd
    CALL mem_size_eliashberg(1, imelt)
    !
    ALLOCATE(ixkqf(nkfs, nqtotfd), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_g2_tdbe', 'Error allocating ixkqf', 1)
    ALLOCATE(nqfs(nkfs), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_g2_tdbe', 'Error allocating nqfs', 1)
    ALLOCATE(index_(lower_bnd:upper_bnd, nqtotfd), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_g2_tdbe', 'Error allocating index_', 1)
    allocate(iq2nqfs(nkfs,nqtotfd), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_g2_tdbe', 'Error allocating iq2nqfs', 1)
    ixkqf(:, :) = 0
    nqfs(:) = 0
    index_(:, :) = 0
    iq2nqfs(:, :) = 0
    ! find the index of k+sign*q on the fine k-mesh
    ! nkfs - total nr. of k-points within the Fermi shell (fine mesh)
    !      - these are irreducible k-points if mp_mesh_k=.TRUE.
    ! nqtotf - total nr of q-points on the fine mesh
    !
    DO ik = lower_bnd, upper_bnd
      DO iqq = 1, totqfd
        iq = selecqfd(iqq)
        xk(:) = xkfs(:, ik)
        xq(:) = xqfd(:, iq)
        !
        !  nkq - index of k+sign*q on the full fine k-mesh.
        !
        CALL kpmq_map_dg(xk, xq, +1, nkq)
        !
        !  ixkqf(ik,iq) - index of k+sign*q on the fine k-mesh within the Fermi shell
        !
        ixkqf(ik, iq) = ixkff(nkq)
        !
        ! nqfs(ik) - nr of q-points at each k-point for which k+sign*q is within the Fermi shell
        ! index_   - index q-point on the full q-mesh for which k+sign*q is within the Fermi shell
        !
        IF (ixkqf(ik, iq) > 0) THEN
          nqfs(ik) = nqfs(ik) + 1
          index_(ik, nqfs(ik)) = iq
          iq2nqfs(ik,iq) = nqfs(ik)
        ENDIF
      ENDDO
    ENDDO
    !
    ! collect contributions from all pools (sum over k-points)
    CALL mp_sum(ixkqf, inter_pool_comm)
    CALL mp_sum(nqfs,  inter_pool_comm)
    CALL mp_sum(iq2nqfs, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    !
    ! get the size of required memory for ixqfs
    imelt = nkfs * MAXVAL(nqfs(:))
    CALL mem_size_eliashberg(1, imelt)
    !
    ALLOCATE(ixqfs(nkfs, MAXVAL(nqfs(:))), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_g2_tdbe', 'Error allocating ixqfs', 1)
    ixqfs(:, :) = 0
    !
    DO ik = lower_bnd, upper_bnd
      DO iq = 1, nqfs(ik)
        !
        ! ixqfs - index q-point on the full q-mesh for which k+sign*q is within the Fermi shell
        !
        ixqfs(ik,iq) = index_(ik,iq)
      ENDDO
    ENDDO
    !
    ! collect contributions from all pools (sum over k-points)
    CALL mp_sum(ixqfs, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    !
    DEALLOCATE(index_, STAT = ierr)
    IF (ierr /= 0) CALL errore('read_g2_tdbe', 'Error deallocating index_', 1)
    !
    ! remove memory allocated for index_
    imelt = nqtotfd * (upper_bnd - lower_bnd + 1)
    CALL mem_size_eliashberg(1, -imelt)
    !
    DEALLOCATE(selecqfd, STAT = ierr)
    IF (ierr /= 0) CALL errore('read_g2_tdbe', 'Error allocating selecq', 1)
    !
    WRITE(stdout, '(/5x, a, i9/)') 'Max nr of q-points = ', MAXVAL(nqfs(:))
    WRITE(stdout, '(/5x, a/)') 'Finished reading ikmap files'
    !
    ! get the size of the e-ph matrices that need to be stored in each pool
    imelt = (upper_bnd - lower_bnd + 1) * MAXVAL(nqfs(:)) * nbndfs**2 * nmodes
    CALL mem_size_eliashberg(2, imelt)
    !
    ALLOCATE(nnq(nkfs), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_g2_tdbe', 'Error allocating nnq', 1)
    ALLOCATE(g2(lower_bnd:upper_bnd , MAXVAL(nqfs(:)), nbndfs, nbndfs, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_g2_tdbe', 'Error allocating g2', 1)
    g2(:, :, :, :, :) = zero   
    ! eps_acoustic is given in units of cm-1 in the input file and converted to Ryd in epw_readin
    WRITE(stdout, '(/5x, a/)') 'Start reading .ephmat files'
    !
    dirname = TRIM(ephmat_dir) // TRIM(prefix) // '.ephmat'
    !
    DO ipool = 1, npool_g2 ! nr of pools
      CALL set_ndnmbr(0, ipool, 1, npool_g2, filelab)
#if defined(__MPI)
      filephmat = TRIM(dirname) // '/' // 'ephmat' // filelab
#else
      filephmat = TRIM(dirname) // '/' // 'ephmat'
#endif
      OPEN(UNIT = iufileph, FILE = filephmat, STATUS = 'unknown', &
           FORM = 'unformatted', ACCESS = 'stream', IOSTAT = ios)
      IF (ios /= 0) CALL errore('read_g2_tdbe', 'error opening file ' // filephmat, iufileph)

      READ(iufileph) tmp_pool_id, nkpool(ipool)
      IF (ipool /= tmp_pool_id)  & 
          CALL errore('read_g2_tdbe', &
                 'npool_g2 should be equal to the number of .ephmat files', 1)
      IF (ipool > 1) &
        nkpool(ipool) = nkpool(ipool) + nkpool(ipool - 1)

      CLOSE(iufileph)
    ENDDO
    CALL mp_barrier(inter_pool_comm)
    !
    ! since the nkfs k-points within the Fermi shell are not evenly distrubed
    ! among the .ephmat files, we re-distribute them here among the npool-pools
    nmin = npool_g2
    nmax = npool_g2
    DO ipool = npool_g2, 1, -1
      IF (lower_bnd <= nkpool(ipool)) THEN
        nmin = ipool
      ENDIF
      IF (upper_bnd <= nkpool(ipool)) THEN
        nmax = ipool
      ENDIF
    ENDDO
    !
    nnk = 0
    nnk_all(:)=0
    nnq(:) = 0
    DO ipool = 1, npool_g2 ! nr of pools
      CALL set_ndnmbr(0, ipool, 1, npool_g2, filelab)
#if defined(__MPI)
      filephmat = TRIM(dirname) // '/' // 'ephmat' // filelab
#else
      filephmat = TRIM(dirname) // '/' // 'ephmat'
#endif
      OPEN(UNIT = iufileph, FILE = filephmat, STATUS = 'unknown', &
           FORM ='unformatted', ACCESS = 'stream', IOSTAT = ios)
      IF (ios /= 0) CALL errore('read_g2_tdbe', 'error opening file ' // filephmat, iufileph)
      READ(iufileph) tmp_pool_id, nks
      IF (ipool >= nmin .AND. ipool <= nmax) THEN
        DO iq = 1, nqtotfd ! loop over q-points
          IF (lscreen_tdbe .AND. scr_typ <= 1 ) THEN
            xq = xqfd(:, iq)
            CALL shift_q(xq, qcart)  ! shift q point to the 1st BZ
            IF (scr_typ == 0) CALL rpa_epsilon(xq, wfd(:, iq), nmodes, epsi, eps_rpa)
            IF (scr_typ == 1) CALL tf_epsilon(xq, nmodes, epsi, eps_rpa)
          ENDIF
          DO ik = 1, nks ! loop over k-points in the pool
            IF (ixkqf(ik+nnk, iq) > 0) THEN
              nnq(ik + nnk) = nnq(ik + nnk) + 1
              DO imode = 1, nmodes ! loop over phonon modes
                DO ibnd = 1, nbndfs ! loop over iband's
                  IF (ABS(ekfs(ibnd, ik+nnk) - efm0) < fsthick) THEN
                     DO jbnd = 1, nbndfs ! loop over jband's
                      IF (ABS(ekfs(jbnd, ixkqf(ik + nnk, iq)) - efm0) < fsthick) THEN
                         !READ(iufileph, '(ES20.10)') gmat
                         READ(iufileph) gmat
                         IF (ik+nnk >= lower_bnd .AND. ik+nnk <= upper_bnd) THEN
                           IF (wfd(imode, iq) > eps_acoustic) THEN
                             !g2(ik + nnk, nnq(ik + nnk), ibnd, jbnd, imode) = gmat * ryd2ev * ryd2ev
                             g2(ik + nnk, nnq(ik + nnk), ibnd, jbnd, imode) = gmat
                             IF (lscreen_tdbe .AND. scr_typ <= 1) g2(ik + nnk, nnq(ik + nnk), ibnd, jbnd, imode)  &
                                        = g2(ik + nnk, nnq(ik + nnk), ibnd, jbnd, imode)/ &
                                          ABS(eps_rpa(imode))**two 
                           ELSE
                            g2(ik + nnk, nnq(ik + nnk), ibnd, jbnd, imode) = zero
                           ENDIF
                         ENDIF
                       ENDIF ! ekq
                     ENDDO ! jbnd
                   ENDIF ! ekk
                ENDDO ! ibnd
              ENDDO ! imode
            ENDIF ! ekk and ekq
          ENDDO ! ik
        ENDDO ! iq
        CLOSE(iufileph)
      ENDIF ! ipool
      nnk = nnk + nks
      IF (ipool == npool_g2 .AND. nnk /= nkfs)  CALL errore('read_ephmat_dg', &
          'nnk should be equal to nkfs', 1)
    ENDDO ! ipool
    !
    CALL mp_barrier(inter_pool_comm)
    !
    WRITE(stdout, '(/5x, a/)') 'Finished reading .ephmat files'
    !
    IF (lscreen_tdbe .AND. scr_typ <= 1 ) THEN
      DEALLOCATE(eps_rpa, STAT = ierr)
      IF (ierr /= 0) CALL errore('read_g2_tdbe', 'Error deallocating eps_rpa', 1)
    ENDIF
    !
    DEALLOCATE(nnk_all, STAT = ierr)
    IF (ierr /= 0) CALL errore('read_g2_tdbe', 'Error deallocating nnk_all', 1)
    DEALLOCATE(nkpool, STAT = ierr)
    IF (ierr /= 0) CALL errore('read_g2_tdbe', 'Error deallocating nkpool', 1)    
    DEALLOCATE(wfd, STAT = ierr)
    IF (ierr /= 0) CALL errore('read_g2_tdbe', 'Error deallocating wfd', 1)
    DEALLOCATE(wqfd, STAT = ierr)
    IF (ierr /= 0) CALL errore('read_g2_tdbe', 'Error deallocating wqfd', 1)
    DEALLOCATE(xqfd, STAT = ierr)
    IF (ierr /= 0) CALL errore('read_g2_tdbe', 'Error deallocating xqfd', 1)
    DEALLOCATE(wkfs, STAT = ierr)
    IF (ierr /= 0) CALL errore('read_g2_tdbe', 'Error deallocating wkfs', 1)
    DEALLOCATE(nnq, STAT = ierr)
    IF (ierr /= 0) CALL errore('read_g2_tdbe', 'Error deallocating nnq', 1)
    DEALLOCATE(ekfs_all, STAT = ierr)
    IF (ierr /= 0) CALL errore('read_g2_tdbe', 'Error deallocating ekfs_all', 1)
    DEALLOCATE(wkfs_all, STAT = ierr)
    IF (ierr /= 0) CALL errore('read_g2_tdbe', 'Error deallocating wkfs_all', 1)
    DEALLOCATE(xkfs_all, STAT = ierr)
    IF (ierr /= 0) CALL errore('read_g2_tdbe', 'Error deallocating wkfs_all', 1)
    DEALLOCATE(ixkf, STAT = ierr)
    IF (ierr /= 0) CALL errore('read_g2_tdbe', 'Error deallocating ixkf', 1)
    DEALLOCATE(ixkf_inv, STAT = ierr)
    IF (ierr /= 0) CALL errore('read_g2_tdbe', 'Error deallocating ixkf_inv', 1)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE read_g2_tdbe
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE shift_q(q, qcart)
    !-----------------------------------------------------------------------
    !! shift the q point to the wigner-seitz cell
    !! 
    USE kinds,            ONLY : DP 
    USE cell_base,        ONLY : at, bg, alat, omega
    !   
    REAL(KIND = DP) , INTENT(inout) :: q(3)
    !! q-coordinate
    REAL(KIND = DP) , INTENT(out) :: qcart(3) 
    INTEGER :: m1, m2, m3
    !! Counter for G vectors
    INTEGER :: m1f, m2f, m3f
    !! index of G that minimizes q+G
    REAL :: g1, g2, g3
    !! Cartesian coordinates of q+G
    REAL :: val
    !! test value for finding G that minimizes q+G
    REAL :: qmG
    !! |q+G|
    CALL cryst_to_cart(1, q, bg, 1)
    ! Look for G that minimizes q+G --> enforce periodicity
    !
    val = 1.0d24
    m1f = 0
    m2f = 0
    m3f = 0
    !
    !
    DO m1 = -2, 2
      DO m2 = -2, 2
        DO m3 = -2, 2
          g1 = q(1) + (m1 * bg(1, 1) + m2 * bg(1, 2) + m3 * bg(1, 3))
          g2 = q(2) + (m1 * bg(2, 1) + m2 * bg(2, 2) + m3 * bg(2, 3))
          g3 = q(3) + (m1 * bg(3, 1) + m2 * bg(3, 2) + m3 * bg(3, 3))
          qmG = SQRT(g1**2+g2**2+g3**2)
          IF (qmG < val) THEN
            val = qmG
            m1f = m1
            m2f = m2
            m3f = m3
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    !
    q(1) = q(1) + m1f * bg(1, 1) + m2f * bg(1, 2) + m3f * bg(1, 3)
    q(2) = q(2) + m1f * bg(2, 1) + m2f * bg(2, 2) + m3f * bg(2, 3)
    q(3) = q(3) + m1f * bg(3, 1) + m2f * bg(3, 2) + m3f * bg(3, 3)
    qcart(:) = q(:)
    CALL cryst_to_cart(1, q, at, -1)
    ! 
    !-----------------------------------------------------------------------
    END SUBROUTINE shift_q
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE kpmq_map_dg(xk, xq, sign1, nkq)
    !-----------------------------------------------------------------------
    !!
    !! this routine finds the index of k+q or k-q point on the fine k-mesh
    !!
    USE kinds,          ONLY : DP
    USE input,         ONLY : nkf1d, nkf2d, nkf3d
    USE ep_constants,  ONLY : eps5
    USE mp,             ONLY : mp_bcast, mp_barrier
    USE kfold,          ONLY : backtoBZ
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: sign1
    !! +1 for searching k+q, -1 for k-q
    INTEGER, INTENT(out) :: nkq
    !! the index of k+sign*q
    !
    REAL(KIND = DP), INTENT(in) :: xk(3)
    !! coordinates of k points
    REAL(KIND = DP), INTENT(in) :: xq(3)
    !! coordinates of q points
    !
    ! Local variables
    LOGICAL :: in_the_list
    !! Check if k point is in the list
    !
    REAL(KIND = DP) :: xx, yy, zz
    !! Temporary variables
    REAL(KIND = DP) :: xxk(3)
    !! k + (sign1) * q
    !
    xxk(:) = xk(:) + DBLE(sign1) * xq(:)
    xx = xxk(1) * nkf1d
    yy = xxk(2) * nkf2d
    zz = xxk(3) * nkf3d
    in_the_list = ABS(xx - NINT(xx)) <= eps5 .AND. &
                  ABS(yy - NINT(yy)) <= eps5 .AND. &
                  ABS(zz - NINT(zz)) <= eps5
    IF (.NOT. in_the_list) CALL errore('kpmq_map', 'k+q does not fall on k-grid', 1)
    !
    !  find the index of this k+q or k-q in the k-grid
    !  make sure xx, yy, zz are in the 1st BZ
    !
    CALL backtoBZ(xx, yy, zz, nkf1d, nkf2d, nkf3d)
    !
    ! since k- and q- meshes are commensurate, nkq can be easily found
    !
    nkq = NINT(xx) * nkf2d * nkf3d + NINT(yy) * nkf3d + NINT(zz) + 1
    !
    !  Now nkq represents the index of k+sign*q on the fine k-grid.
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE kpmq_map_dg
    !-----------------------------------------------------------------------
    !
    !----------------------------------------------------------------------
    SUBROUTINE deallocate_tdbe()
    !-----------------------------------------------------------------------
    !
      !!  deallocates the variables allocated by read_frequencies, 
      !!  read_eigenvalues, read_kqmap, read_ephmat, and tdbte

      USE input,        ONLY : mp_mesh_k, phph_tdbe, vme, &
                                lscreen_tdbe, scr_typ, assume_metal
      USE global_var,         ONLY : wf, wqf, xqf,adapt_smearing, zstar, &
                                epsi 
      USE supercond_common, ONLY : ekfs, xkfs, wkfs, g2, ixkff,        &
                                ixkqf, ixqfs, nqfs,  memlt_pool
      USE tdbe_common,       ONLY : enk_all, wkf_all, xkf_all, vnk_all, & 
                                vnuq_all, uf_all, selecq,           &
                                bztoibz_dg, s_bztoibz_dg, ie_ph,    &
                                iph_e, iph_ph, ie_e, felec, nphon,  &
                                nphon_pre, iq2nqfs, indx_mapk,      &
                                indx_mapq, ikfsdx, unk_all, Amnk_all               
      USE ions_base,     ONLY : ityp, tau
      !
      IMPLICIT NONE
      !
      INTEGER :: ierr
      !! Error status
      !! 
      IF(ALLOCATED(epsi))     DEALLOCATE(epsi)
      IF(ALLOCATED(zstar))    DEALLOCATE(zstar)
      IF(ALLOCATED(ityp))     DEALLOCATE(ityp)
      IF(ALLOCATED(tau))      DEALLOCATE(tau)
      ! read_ephmat
      DEALLOCATE(g2, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_tdbe', 'Error deallocating g2', 1)
      ! read_frequencies
      DEALLOCATE(wf, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_tdbe', 'Error deallocating wf', 1)
      DEALLOCATE(wqf, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_tdbe', 'Error deallocating wqf', 1)
      DEALLOCATE(xqf, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_tdbe', 'Error deallocating xqf', 1)
      ! read_eigenvalues
      DEALLOCATE(ekfs, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_tdbe', 'Error deallocating ekfs', 1)
      DEALLOCATE(xkfs, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_tdbe', 'Error deallocating xkfs', 1)
      ! read_kqmap
      DEALLOCATE(ixkff, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_tdbe', 'Error deallocating ixkff', 1)
      DEALLOCATE(ixkqf, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_tdbe', 'Error deallocating ixkqf', 1)
      DEALLOCATE(ixqfs, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_tdbe', 'Error deallocating ixqfs', 1)
      DEALLOCATE(nqfs, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_tdbe', 'Error deallocating nqfs', 1)
      DEALLOCATE(memlt_pool, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_tdbe', 'Error deallocating memlt_pool', 1)
      DEALLOCATE(iq2nqfs, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_tdbe', 'Error deallocating iq2nqfs', 1)
      DEALLOCATE(indx_mapk, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_tdbe', 'Error deallocating indx_mapk', 1)
      DEALLOCATE(indx_mapq, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_tdbe', 'Error deallocating indx_mapq', 1)
      DEALLOCATE(ikfsdx, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_tdbe', 'Error deallocating indx_mapq', 1)
      !!! Deallocate variables defined in module tdbe_mod.f90
      DEALLOCATE(enk_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_tdbe','Error deallocating enk_all',1)
      DEALLOCATE(xkf_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_tdbe','Error deallocating xkf_all',1)
      DEALLOCATE(wkf_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_tdbe','Error deallocating wkf_all',1)
      IF (adapt_smearing) THEN
        DEALLOCATE(vnk_all, STAT = ierr)
        IF (ierr /= 0) CALL errore('deallocate_tdbe','Error deallocating vnuq_all',1)
        DEALLOCATE(vnuq_all, STAT = ierr)
        IF (ierr /= 0) CALL errore('deallocate_tdbe','Error deallocating vnuq_all',1)
      ENDIF
      IF (lscreen_tdbe .AND. scr_typ > 1 .AND. (.NOT. assume_metal)) THEN
        DEALLOCATE(unk_all, STAT = ierr) 
        IF (ierr /= 0) CALL errore('wanninterp', 'Error deallocating cufkk_all', 1)
        IF (vme == 'wannier') THEN 
          DEALLOCATE(Amnk_all, STAT = ierr)
          IF (ierr /= 0) CALL errore('wanninterp', 'Error allocating Amnk_all', 1)
        ENDIF
      ENDIF
      DEALLOCATE(uf_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_tdbe','Error deallocating uf',1)
      DEALLOCATE(selecq, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_tdbe','Error deallocating selecq',1)
      IF (mp_mesh_k) THEN
      DEALLOCATE(bztoibz_dg, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_tdbe','Error deallocating bztoibz_dg',1)
      DEALLOCATE(s_bztoibz_dg , STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_tdbe','Error deallocating s_bztoibz_dg',1)
      ENDIF
      DEALLOCATE(ie_ph, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_tdbe','Error deallocating ie_ph',1)
      DEALLOCATE(iph_e, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_tdbe','Error deallocating iph_e',1)
      IF (phph_tdbe) THEN
        DEALLOCATE(iph_ph, STAT = ierr)
        IF (ierr /= 0) CALL errore('deallocate_tdbe','Error deallocating iph_ph',1)
        DEALLOCATE(nphon_pre, STAT = ierr)
        IF (ierr /= 0) CALL errore('deallocate_tdbe','Error deallocating nphon_pre',1)
      ENDIF
      !
      DEALLOCATE(felec, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_tdbe','Error deallocating felec',1)
      DEALLOCATE(nphon, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_tdbe','Error deallocating nphon',1)
      ! evaluate_a2f_lambda
      !
      RETURN
    !-----------------------------------------------------------------------
    END SUBROUTINE deallocate_tdbe
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE rotate_kpoint(ns, ik, iks, nkc1, nkc2, nkc3)
    !-----------------------------------------------------------------------
    !! find out the index of the k point index (iks) which is reached by
    !! rotating the k point ik by symmetry operation ns 
    !-----------------------------------------------------------------------
    USE kinds,            ONLY : DP
    USE symm_base,        ONLY : s, t_rev, invs
    USE ep_constants,    ONLY : eps6
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: ns
    !! index of symmetry operation
    INTEGER, INTENT(in) :: ik 
    !! index of k point
    INTEGER :: nkc1, nkc2, nkc3
    !! size of the k grid
    INTEGER, INTENT(out) :: iks
    !! index of k point after rotation
    INTEGER :: nsinv
    !! inversion operation contained (or not) in ns
    REAL(KIND = DP) :: xk(3), xks(3)
    !! coordinates of k point before and after rotation
    REAL(KIND = DP) :: xx, yy, zz
    !! temporary variable for coordinates
    INTEGER :: i, j , k , rs, rss
    !! !! Counters
    LOGICAL :: in_the_list
    !! After rotation the k point is on the k mesh defined by nkc 


    ! i = 1 + INT(ik /(nkc2*nkc3))
    ! j = 1 + INT((ik-(i-1)*nkc2*nkc3)/nkc3)
    ! k = ik - (i - 1) * nkc2*nkc3 - (j -1)*nkc3
    k = 1 + MODULO(ik - 1, nkc3)
    rs = INT((ik-1)/nkc3)
    j = 1 + MODULO(rs, nkc2)
    rss= INT(rs / nkc2)
    i = 1 +MODULO(rss,nkc1)
    xk(1) = DBLE(i - 1) / DBLE(nkc1)
    xk(2) = DBLE(j - 1) / DBLE(nkc2)
    xk(3) = DBLE(k - 1) / DBLE(nkc3)
    nsinv = invs(ns)
    DO i = 1, 3
      xks(i) = s(i, 1, nsinv) * xk(1) &
             + s(i, 2, nsinv) * xk(2) &
             + s(i, 3, nsinv) * xk(3) 
      xks(i) = xks(i) - NINT(xks(i)) 
    ENDDO
    IF (t_rev(nsinv) == 1) xks = -xks
    xx = xks(1) * nkc1
    yy = xks(2) * nkc2
    zz = xks(3) * nkc3

    in_the_list = ABS(xx - NINT(xx)) <= eps6 .AND. &
                  ABS(yy - NINT(yy)) <= eps6 .AND. &
                  ABS(zz - NINT(zz)) <= eps6 

    IF (in_the_list) THEN
      i = MOD(NINT(xks(1) * nkc1 + 2 * nkc1), nkc1) + 1  
      j = MOD(NINT(xks(2) * nkc2 + 2 * nkc2), nkc2) + 1
      k = MOD(NINT(xks(3) * nkc3 + 2 * nkc3), nkc3) + 1
      iks = (k - 1) + (j - 1) * nkc3 + (i - 1) * nkc2 * nkc3 + 1
    ELSE
      iks = 0
    ENDIF
    END SUBROUTINE rotate_kpoint
    !------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------
    SUBROUTINE allocate_tdbe
    !------------------------------------------------------------------------------ 
    !! Allocate space to time-dependent variables 
    !-----------------------------------------------------------------------
    !
    USE kinds,            ONLY : DP 
    USE input,            ONLY : solver_tdbe, phph_tdbe
    USE tdbe_common,      ONLY : b_tab, c_tab, a_tab, nstg, felec,  &
                                 nphon, ie_ph, iph_e, iph_ph, ie_e, &
                                 nphon_pre
    USE global_var,       ONLY : nktotf, nqtotf, nbndfst
    USE modes,            ONLY : nmodes
    !
    IMPLICIT NONE
    !
    INTEGER :: ierr
    !! Error status
    !
    b_tab = 0.0E0_DP
    c_tab = 0.0E0_DP
    a_tab = 0.0E0_DP
    !
    SELECT CASE(solver_tdbe)
    CASE ('euler')
      nstg = 1
      b_tab(1) = 1.0E0_DP
    CASE ('heun')
      nstg = 2
      b_tab(1) = 0.5E0_DP
      b_tab(2) = 0.5E0_DP
      c_tab(2) = 1.0E0_DP
      a_tab(2, 1) = 1.0E0_DP
    CASE ('rk_4')
      nstg = 4
      b_tab(1) = 1.0E0_DP/6.0E0_DP
      b_tab(2) = 1.0E0_DP/3.0E0_DP
      b_tab(3) = 1.0E0_DP/3.0E0_DP
      b_tab(4) = 1.0E0_DP/6.0E0_DP
      c_tab(2) = 0.5E0_DP
      c_tab(3) = 0.5E0_DP
      c_tab(4) = 1.0E0_DP
      a_tab(2, 1) = 0.5E0_DP
      a_tab(3, 1) = 0.0E0_DP
      a_tab(3, 2) = 0.5E0_DP
      a_tab(4, 1) = 0.0E0_DP
      a_tab(4, 2) = 0.0E0_DP
      a_tab(4, 3) = 1.0E0_DP
    CASE DEFAULT  
      CALL errore('allocate_tdbe', 'solver_tdbe is either euler, heun or rk_4', 1)
    END SELECT

    !! Allocate space to electron and phonon distributions
    ALLOCATE(nphon(nmodes, nqtotf), STAT = ierr)
    IF (ierr /= 0) CALL errore('allocate_tdbe', 'Error allocating nphon', 1)
    ALLOCATE(felec(nbndfst, nktotf), STAT = ierr)
    IF (ierr /= 0) CALL errore('allocate_tdbe', 'Error allocating felec', 1)
    !! Allocate space to collision integrals
    ALLOCATE(ie_ph(nbndfst, nktotf, nstg), STAT = ierr)
    IF (ierr /= 0) CALL errore('allocate_tdbe', 'Error allocating ie_ph', 1)
    ALLOCATE(iph_e(nmodes, nqtotf, nstg), STAT = ierr)
    IF (ierr /= 0) CALL errore('allocate_tdbe', 'Error allocating iph_e', 1)
    IF (phph_tdbe) THEN
      allocate(iph_ph(nmodes, nqtotf, nstg), STAT = ierr)
      IF (ierr /= 0) CALL errore('allocate_tdbe', 'Error allocating iph_ph', 1)
      allocate(nphon_pre(nmodes, nqtotf), STAT = ierr)
      IF (ierr /= 0) CALL errore('allocate_tdbe', 'Error allocating nphon_pre', 1)
    ENDIF
    ! Initialize
    felec(:, :) = 0.0E0_DP
    nphon(:, :) = 0.0E0_DP
    ie_ph(:, :, :) = 0.0E0_DP
    iph_e(:, :, :) = 0.0E0_DP
    IF (phph_tdbe) THEN
      nphon_pre(:, :) = 0.0E0_DP
      iph_ph(:, :, :) = 0.0E0_DP
    ENDIF
    !
    RETURN
    !
    !------------------------------------------------------------------------------
    END SUBROUTINE allocate_tdbe
    !------------------------------------------------------------------------------
    !    
    !------------------------------------------------------------------------------
    SUBROUTINE crystal_fmt_read(w_centers)
    !-------------------------------------------------------------------------------
    !!
    !! Routine to read the real space quantities for fine grid interpolation
    !------------------------------------------------------------------------------
    !!
    USE kinds,     ONLY : DP
    USE input,    ONLY : system_2d, nbndsub
    USE pwcom,     ONLY : nelec
    USE ions_base, ONLY : nat, tau, ityp, amass
    USE modes,     ONLY : nmodes
    USE cell_base,     ONLY : at, bg, alat, omega, tpiba
    USE noncollin_module, ONLY : noncolin
    USE io_global, ONLY : stdout
    USE io_var,    ONLY : crystal
    USE io_global, ONLY : ionode_id
    USE mp,        ONLY : mp_barrier, mp_bcast
    USE mp_world,  ONLY : world_comm, mpime
    USE mp_global, ONLY : my_pool_id, inter_pool_comm
    USE global_var,     ONLY : L, do_cutoff_2D_epw, area, nbndskip
    !
    IMPLICIT NONE
    !
    !
    REAL(KIND = DP), INTENT(OUT) :: w_centers(3, nbndsub)
    !! wannier centers
    INTEGER :: ios
    !! Status of files
    INTEGER :: ierr
    !! Error status
    !
    IF (mpime == ionode_id) THEN
      !
      OPEN(UNIT = crystal, FILE = 'crystal.fmt', STATUS = 'old', IOSTAT = ios)
      IF (ios /= 0) CALL errore('crystal_fmt_read', 'error opening crystal.fmt', crystal)
      READ(crystal, *) nat
      READ(crystal, *) nmodes
      READ(crystal, *) nelec, nbndskip
      READ(crystal, *) at
      READ(crystal, *) bg
      READ(crystal, *) omega
      READ(crystal, *) alat
      ALLOCATE(tau(3, nat), STAT = ierr)
      IF (ierr /= 0) CALL errore('crystal_fmt_read', 'Error allocating tau', 1)
      READ(crystal, *) tau
      READ(crystal, *) amass
      ALLOCATE(ityp(nat), STAT = ierr)
      IF (ierr /= 0) CALL errore('crystal_fmt_read', 'Error allocating ityp', 1)
      READ(crystal, *) ityp
      READ(crystal, *) noncolin
      READ(crystal, *) do_cutoff_2D_epw
      READ(crystal, *) w_centers
      READ(crystal, *) L      
    ENDIF
    CALL mp_bcast(nat      , ionode_id, world_comm)
    IF (mpime /= ionode_id) ALLOCATE(ityp(nat))
    CALL mp_bcast(nmodes   , ionode_id, world_comm)
    CALL mp_bcast(nelec    , ionode_id, world_comm)
    CALL mp_bcast(nbndskip , ionode_id, world_comm)
    CALL mp_bcast(at       , ionode_id, world_comm)
    CALL mp_bcast(bg       , ionode_id, world_comm)
    CALL mp_bcast(omega    , ionode_id, world_comm)
    CALL mp_bcast(alat     , ionode_id, world_comm)
    IF (mpime /= ionode_id) ALLOCATE(tau(3, nat))
    CALL mp_bcast(tau      , ionode_id, world_comm)
    CALL mp_bcast(amass    , ionode_id, world_comm)
    CALL mp_bcast(ityp     , ionode_id, world_comm)
    CALL mp_bcast(noncolin , ionode_id, world_comm)
    CALL mp_bcast(w_centers, ionode_id, world_comm)
    CALL mp_bcast(do_cutoff_2D_epw, ionode_id, world_comm)
    CALL mp_bcast(L        , ionode_id, world_comm)
    IF (mpime == ionode_id) THEN
      CLOSE(crystal)
    ENDIF
    IF (system_2d /= 'no') THEN
      area = omega * bg(3, 3) / alat
      WRITE(stdout,'(5x,a,f10.6)') 'Area is [Bohr^2] ', area
    ENDIF
    RETURN
    !
    !------------------------------------------------------------------------------
    END SUBROUTINE crystal_fmt_read
    !------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------
    SUBROUTINE epwdata_fmt_read(nrr_k,nrr_q)
    !--------------------------------------------------------------------------------
    !!
    !! Routine to read the real space quantities for fine grid interpolation
    !------------------------------------------------------------------------------
    USE kinds,          ONLY : DP
    USE input,          ONLY : nbndsub, vme, lifc, nqc1, nqc2, nqc3
    USE ep_constants,   ONLY : zero, czero, cone
    USE pwcom,          ONLY : ef
    USE global_var,     ONLY : chw, rdw,zstar, crrw, epsi,          &
                               adapt_smearing, cvmew, ifc
    USE ions_base,      ONLY : nat
    USE modes,          ONLY : nmodes
    USE io_global,      ONLY : stdout
    USE io_var,         ONLY : epwdata, crystal, iunvmedata
#if defined(__NAG)
    USE f90_unix_io,    ONLY : flush
#endif
    USE io_global,      ONLY : ionode_id
    USE mp,             ONLY : mp_barrier, mp_bcast
    USE mp_global,      ONLY : inter_pool_comm, my_pool_id
    USE mp_world,       ONLY : world_comm, mpime
    USE io,             ONLY : read_ifc_epw
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(OUT) :: nrr_k
    !! Number of WS vectors for the electrons
    INTEGER, INTENT(OUT) :: nrr_q
    !! Number of WS vectors for the phonons
    !
    !
    INTEGER :: ibnd, jbnd, ipol
    !! Band index
    INTEGER :: jmode, imode
    !! Mode index
    INTEGER :: irk, irq
    !! WS vector looping index on electron, phonons and el-ph
    INTEGER :: ios
    !! Status of files
    INTEGER :: ierr
    !! Error status
    !
    WRITE(stdout,'(/5x,"Reading Hamiltonian, Dynamical matrix in Wann rep from epwdata.fmt and crystal.fmt"/)')
    FLUSH(stdout)
    !
    ! This is important in restart mode as zstar etc has not been allocated
    ALLOCATE(zstar(3, 3, nat), STAT = ierr)
    IF (ierr /= 0) CALL errore('epwdata_fmt_read', 'Error allocating zstar', 1)
    ALLOCATE(epsi(3, 3), STAT = ierr)
    IF (ierr /= 0) CALL errore('epwdata_fmt_read', 'Error allocating epsi', 1)
    !
    IF (mpime == ionode_id) THEN
      !
      OPEN(UNIT = epwdata, FILE = 'epwdata.fmt', STATUS = 'old', IOSTAT = ios)
      IF (ios /= 0) CALL errore ('epwdata_fmt_read', 'error opening epwdata.fmt', epwdata)
      READ(epwdata,*) ef
      READ(epwdata,*) nbndsub, nrr_k, nmodes, nrr_q
      READ(epwdata,*) zstar, epsi
      !
    ENDIF
    CALL mp_bcast(ef,      ionode_id, world_comm)
    CALL mp_bcast(nbndsub, ionode_id, world_comm)
    CALL mp_bcast(nrr_k,   ionode_id, world_comm)
    CALL mp_bcast(nmodes,  ionode_id, world_comm)
    CALL mp_bcast(nrr_q,   ionode_id, world_comm)
    CALL mp_bcast(zstar,   ionode_id, world_comm)
    CALL mp_bcast(epsi,    ionode_id, world_comm)
    !
    ALLOCATE(chw(nbndsub, nbndsub, nrr_k), STAT = ierr)
    IF (ierr /= 0) CALL errore('epwdata_fmt_read', 'Error allocating chw', 1)
    ALLOCATE(rdw(nmodes, nmodes,  nrr_q), STAT = ierr)
    IF (ierr /= 0) CALL errore('epwdata_fmt_read', 'Error allocating rdw', 1)
    !
    IF (mpime == ionode_id) THEN
      !
      DO ibnd = 1, nbndsub
        DO jbnd = 1, nbndsub
          DO irk = 1, nrr_k
            READ(epwdata,*) chw(ibnd, jbnd, irk)
          ENDDO
        ENDDO
      ENDDO
      !
      IF (.NOT. lifc) THEN
        DO imode = 1, nmodes
          DO jmode = 1, nmodes
            DO irq = 1, nrr_q
              READ(epwdata,*) rdw(imode, jmode, irq)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      !
     !
    ENDIF
    IF (lifc) THEN
      ALLOCATE(ifc(nqc1, nqc2, nqc3, 3, 3, nat, nat), STAT = ierr)
      IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error allocating ifc', 1)
      ifc(:, :, :, :, :, :, :) = zero
      CALL read_ifc_epw
    ENDIF
    !
    CALL mp_bcast(chw, ionode_id, world_comm)
    CALL mp_bcast(rdw, ionode_id, world_comm)
    !
    !
    !CALL mp_barrier(inter_pool_comm)
    IF (mpime == ionode_id) THEN
      CLOSE(epwdata)
    ENDIF
    !
    WRITE(stdout, '(/5x,"Finished reading epwdata.fmt and crystal.fmt"/)')
    !
    !! Read velocity operator in wannier basis if 
    !! adaptive smearing is turned on
    ! IF (adapt_smearing) THEN
    IF (vme /= 'wannier') CALL errore('epwdata_fmt_read', &
         'velocity operator in wannier basis is required for adaptive smearing',1)
    IF (vme == 'wannier') THEN
      IF (mpime == ionode_id) THEN
        OPEN(UNIT = iunvmedata, FILE = 'vmedata.fmt', STATUS = 'old', IOSTAT = ios)
        IF (ierr /= 0) CALL errore('epwdata_fmt_read', 'Error allocating cvmew', 1)
      ENDIF
      ALLOCATE(cvmew(3, nbndsub, nbndsub, nrr_k), STAT = ierr)
      IF (ierr /= 0) CALL errore('epwdata_fmt_read', 'Error allocating cvmew', 1)
      ALLOCATE(crrw(3, nbndsub, nbndsub, nrr_k), STAT = ierr)
      IF (ierr /= 0) CALL errore('epwdata_fmt_read', 'Error allocating crrw', 1)
      cvmew(:, :, :, :) = czero
      crrw(:, :, :, :)  = czero
      IF (mpime == ionode_id) THEN
        DO ibnd = 1, nbndsub
          DO jbnd = 1, nbndsub
            DO irk = 1, nrr_k
              DO ipol = 1,3
              READ(iunvmedata, *) cvmew(ipol, ibnd, jbnd, irk)
              READ(iunvmedata, *) crrw(ipol, ibnd, jbnd, irk)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        CLOSE(iunvmedata)
      ENDIF
      WRITE(stdout, '(/5x,"Finished reading vmedata.fmt"/)')
      CALL mp_bcast(cvmew, ionode_id, world_comm) 
      CALL mp_bcast(crrw, ionode_id, world_comm) 
    ENDIF
    !
    !------------------------------------------------------------------------------
    END SUBROUTINE epwdata_fmt_read
    !------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------
    SUBROUTINE quadrupole_fmt_read()
    !--------------------------------------------------------------------------------
    !!
    !! Routine to read the real space quantities for fine grid interpolation
    !------------------------------------------------------------------------------
    USE kinds,            ONLY : DP
    USE ep_constants,     ONLY : zero
    USE ions_base,        ONLY : nat
    USE io_global,        ONLY : stdout
    USE io_var,           ONLY : iuquad
    USE io_global,        ONLY : ionode_id
    USE mp,               ONLY : mp_bcast
    USE mp_global,        ONLY : world_comm, inter_pool_comm, my_pool_id
    USE global_var,       ONLY : Qmat, qrpl
    USE input,            ONLY : lpolar
    USE mp_world,         ONLY : mpime
    !
    IMPLICIT NONE
    !
    INTEGER :: na
    !! Counter on atom
    INTEGER :: idir
    !! Cartesian direction
    INTEGER :: i
    !! Counter of modes
    REAL(KIND = DP) :: Qxx, Qyy, Qzz, Qyz, Qxz, Qxy
    !! Specific quadrupole value read from file.
    CHARACTER(LEN = 256) :: dummy
    !! Dummy variable
    INTEGER :: ios
    !! Status of files
    INTEGER :: ierr
    !! Error status
    LOGICAL :: exst
    !! Find if a file exists.
    IF (mpime == ionode_id) THEN
      INQUIRE(FILE = 'quadrupole.fmt', EXIST = exst)
    ENDIF
    CALL mp_bcast(exst, ionode_id, world_comm)
    qrpl = .FALSE.
    ALLOCATE(Qmat(nat, 3, 3, 3), STAT = ierr)
    IF (ierr /= 0) CALL errore('quadrupole_fmt_read', 'Error allocating Qmat', 1)
    Qmat(:, :, :, :) = zero
    IF (exst) THEN
      qrpl = .TRUE.
      IF (mpime == ionode_id) THEN
        OPEN(UNIT = iuquad, FILE = 'quadrupole.fmt', STATUS = 'old', IOSTAT = ios)
        READ(iuquad, *) dummy
        DO i = 1, 3 * nat
          READ(iuquad, *) na, idir, Qxx, Qyy, Qzz, Qyz, Qxz, Qxy
          Qmat(na, idir, 1, 1) = Qxx
          Qmat(na, idir, 2, 2) = Qyy
          Qmat(na, idir, 3, 3) = Qzz
          Qmat(na, idir, 2, 3) = Qyz
          Qmat(na, idir, 3, 2) = Qyz
          Qmat(na, idir, 1, 3) = Qxz
          Qmat(na, idir, 3, 1) = Qxz
          Qmat(na, idir, 1, 2) = Qxy
          Qmat(na, idir, 2, 1) = Qxy
        ENDDO
        CLOSE(iuquad)
      ENDIF ! my_pool_id == ionode_id
      CALL mp_bcast(Qmat, ionode_id, world_comm)
      WRITE(stdout, '(a)') '     '
      WRITE(stdout, '(a)') '     ------------------------------------ '
      WRITE(stdout, '(a)') '     Quadrupole tensor is correctly read: '
      WRITE(stdout, '(a)') '     ------------------------------------ '
      WRITE(stdout, '(a)') '     atom   dir        Qxx       Qyy      Qzz        Qyz       Qxz       Qxy'
      DO na = 1, nat
        WRITE(stdout, '(i8, a,6f10.5)' ) na, '        x    ', Qmat(na, 1, 1, 1), Qmat(na, 1, 2, 2), Qmat(na, 1, 3, 3), &
                                                              Qmat(na, 1, 2, 3), Qmat(na, 1, 1, 3), Qmat(na, 1, 1, 2)
        WRITE(stdout, '(i8, a,6f10.5)' ) na, '        y    ', Qmat(na, 2, 1, 1), Qmat(na, 2, 2, 2), Qmat(na, 2, 3, 3), &
                                                              Qmat(na, 2, 2, 3), Qmat(na, 2, 1, 3), Qmat(na, 2, 1, 2)
        WRITE(stdout, '(i8, a,6f10.5)' ) na, '        z    ', Qmat(na, 3, 1, 1), Qmat(na, 3, 2, 2), Qmat(na, 3, 3, 3), &
                                                              Qmat(na, 3, 2, 3), Qmat(na, 3, 1, 3), Qmat(na, 3, 1, 2)
      ENDDO
      WRITE(stdout, '(a)') '     '
    ENDIF ! exst
    !
    IF (lpolar) THEN
      WRITE(stdout, '(/,5x,a)' ) 'Computes the analytic long-range interaction for polar materials [lpolar]'
      WRITE(stdout, '(5x,a)' )   ' '
    ENDIF
    IF (.NOT. lpolar .AND. qrpl) THEN
      WRITE(stdout, '(/,5x,a)' ) 'Computes the analytic quadrupole long-range interaction for non-polar materials [Q1,Q2]'
      WRITE(stdout, '(5x,a)' )   ' '
    ENDIF
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE quadrupole_fmt_read
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE check_dg()
    !-----------------------------------------------------------------------
    !! Check the k points in case of mp_mesh_k = .true.
    !------------------------------------------------------------------------------
    USE kinds,             ONLY : DP
    USE input,             ONLY : nkf1d, nkf2d, nkf3d, nqf1d, nqf2d, nqf3d
    USE symm_base,         ONLY : nrot,s, t_rev, time_reversal, nsym
    USE symmetry,          ONLY : s_k, nsym_k, t_rev_k
    USE supercond_common,  ONLY : ixkff, xkfs, nkfs
    USE tdbe_common,       ONLY : bztoibz_dg, s_bztoibz_dg, nkirr_dg,       &
                                  nkfsden, ikfsdx, indx_mapk
    USE control_flags,     ONLY : iverbosity
    USE mp_global,         ONLY : my_pool_id
    USE io_global,         ONLY : ionode_id
    USE global_var,        ONLY : bztoibz, s_bztoibz
    USE bzgrid,            ONLY : kpoint_grid_epw
    !
    IMPLICIT NONE
    !
    INTEGER  :: ierr
    !! Error info
    INTEGER :: ikfs
    !! Counter of k points within Fermi window
    INTEGER ::  ik
    !! Index of k point on the full mesh
    INTEGER :: ix, iy, iz
    !! Index of k point in x, y,z direction
    INTEGER :: ikd
    !! Index of k point on the grid defined by nkf1d
    INTEGER :: ikibz, ikibz2
    !! Index of irrducible k point 
    REAL(KIND = DP) :: xxk(3)    
    !! coordinates
    !
    ALLOCATE(bztoibz_dg(nkf1d * nkf2d * nkf3d), STAT = ierr)
    IF (ierr /= 0) CALL errore('check_dg', 'error allocating bztoibz_dg',1)
    ALLOCATE(s_bztoibz_dg(nkf1d * nkf2d * nkf3d), STAT = ierr)
    IF (ierr /= 0) CALL errore('check_dg', 'error allocating s_bztoibz_dg',1)
    bztoibz_dg(:) = 0
    s_bztoibz_dg(:) = 0
    IF (ALLOCATED(bztoibz)) THEN
      DEALLOCATE(bztoibz, STAT = ierr)
      IF (ierr /= 0) CALL errore('check_dg', 'error deallocating bztoibz', 1) 
    ENDIF
    IF (ALLOCATED(s_bztoibz)) THEN
      DEALLOCATE(s_bztoibz, STAT = ierr)
      IF (ierr /= 0) CALL errore('check_dg', 'error deallocating s_bztoibz', 1) 
    ENDIF
    CALL kpoint_grid_epw(nsym_k, .FALSE., s_k, t_rev_k, nkf1d, nkf2d, &
                          nkf3d, nkirr_dg, .FALSE.)
    bztoibz_dg = bztoibz
    s_bztoibz_dg = s_bztoibz
    DEALLOCATE(bztoibz, STAT = ierr)
    IF (ierr /= 0) CALL errore('check_dg', 'error deallocating bztoibz', 1) 
    DEALLOCATE(s_bztoibz, STAT = ierr)
    IF (ierr /= 0) CALL errore('check_dg', 'error deallocating s_bztoibz', 1)
    !!! Now I want to check xkfs(ixkff(ikd,:)) and xkfd(s_bztoibz(ikd),:) are the same point
    IF (iverbosity == 5) THEN
      IF (my_pool_id == ionode_id) THEN
        OPEN(unit = 565, file = 'BZtoIBZ_dg')
        OPEN(unit = 566, file = 's_BZtoIBZ_dg')
        DO ik = 1 , nkf1d*nkf2d*nkf3d
          WRITE(565,*) bztoibz_dg(ik)
          WRITE(566,*) s_bztoibz_dg(ik)
        ENDDO
        CLOSE(565)
        CLOSE(566)
      ENDIF
    ENDIF
    !
    DO ikfs = 1, nkfsden
      ik = ikfsdx(ikfs)
      ikd = indx_mapk(ik)
      IF (ixkff(ikd) /= 0) THEN
        xxk(:) = xkfs(:,ixkff(ikd))
        ikibz = bztoibz_dg(ikd)
        ix = MOD(NINT(xxk(1) * nkf1d + 2 * nkf1d), nkf1d) + 1
        iy = MOD(NINT(xxk(2) * nkf2d + 2 * nkf2d), nkf2d) + 1
        iz = MOD(NINT(xxk(3) * nkf3d + 2 * nkf3d), nkf3d) + 1
        ikibz2 = (iz - 1) + (iy - 1) * nkf3d + (ix - 1) * nkf2d * nkf3d + 1
        IF (ikibz /= ikibz2) CALL errore('check_dg', &
      'k points in IBZ are not consistent, turn off mp_mesh_k and recalculate',1) 
      ENDIF
    ENDDO
  RETURN
  !
  !-----------------------------------------------------------------------
  END SUBROUTINE check_dg
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------  
  SUBROUTINE restart_tdbe_read()
  !-----------------------------------------------------------------------
  !! Restart time-dependent simulation by reading the electron and phonon
  !! distribution 
  !-----------------------------------------------------------------------
  USE kinds,            ONLY : DP
  USE input,            ONLY : phph_tdbe, lscreen_tdbe, scr_typ, &
                               assume_metal
  USE modes,            ONLY : nmodes
  USE io_global,        ONLY : stdout, ionode_id, meta_ionode
  USE io_var,           ONLY : iufelecrestart, iunphonrestart
  USE global_var,       ONLY : nqtotf, nktotf, nbndfst
  USE tdbe_common,      ONLY : tstart, felec, nphon, &
                               nphon_pre, epsil_ipa
  USE mp,               ONLY : mp_bcast
  USE mp_global,        ONLY : inter_pool_comm, my_pool_id
  !
  IMPLICIT NONE
  !
  INTEGER :: ierr
  !! Error status
  INTEGER :: ios
  !! INTEGER variable for I/O control
  INTEGER :: ik
  !! Counter on k points
  INTEGER :: nu 
  !! Counter on phonon mode index
  INTEGER :: iq
  !! Counter on q points 
  INTEGER :: ibnd 
  !! Counter on band index
  REAL(KIND = DP) :: trestart
  !! Restart time
  REAL(KIND = DP) :: qcart(3), q2
  !!
  tstart = 0.d0
    WRITE(stdout,'(5x,a)') 'TDBE simulation restarting from an interupted calculation or a user defined inital condition.'
    WRITE(stdout,'(5x,a)') 'electron and phonon distribution read from felec_restart and nphon_restart'
    IF (my_pool_id == ionode_id) THEN
      OPEN(UNIT = iufelecrestart, FILE = 'felec_restart',STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
      IF (ios /= 0) CALL errore('restart_tdbe_read', 'error opening file felec_restart', iufelecrestart)
      READ (iufelecrestart, *) trestart
      DO ik = 1, nktotf
        DO ibnd = 1, nbndfst
          READ(iufelecrestart, *) felec(ibnd, ik)
        ENDDO
      ENDDO      
      OPEN(UNIT = iunphonrestart, FILE = 'nphon_restart',STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
      IF (ios /= 0) CALL errore('restart_tdbe_read', 'error opening file nphon_restart', iunphonrestart)
      READ (iunphonrestart, *) trestart
      DO iq = 1, nqtotf
        DO nu = 1 , nmodes
          READ(iunphonrestart, *) nphon(nu, iq)
        ENDDO
      ENDDO
      IF (phph_tdbe) THEN
        nphon_pre(:, :) =  nphon(:, :)
      ENDIF
      WRITE(stdout,'(5x,a)') 'Fininsh reading electron and phonon distribution, remove files felec_restart and nphon_restart'
      CLOSE(iufelecrestart)
      CLOSE(iunphonrestart)
!       CLOSE(unit=iufelecrestart,STATUS = 'delete')
!       CLOSE(unit=iunphonrestart,STATUS = 'delete')
    ENDIF
    CALL mp_bcast(trestart, ionode_id, inter_pool_comm)
    tstart = trestart
    CALL mp_bcast(felec, ionode_id, inter_pool_comm)
    CALL mp_bcast(nphon, ionode_id, inter_pool_comm)
    !
    IF (phph_tdbe)   CALL mp_bcast(nphon_pre, ionode_id, inter_pool_comm) 
    !
    IF (lscreen_tdbe .AND. scr_typ > 1 .AND. (.NOT. assume_metal)) THEN
      !ALLOCATE(epsil_ipa(nqtotf), STAT = ierr)
      !IF (ierr /= 0) CALL errore('epsil_rpa', 'Error allocating epsil_ipa', 1)
      !
      epsil_ipa = 0.d0
      !
      IF (my_pool_id == ionode_id) THEN
      OPEN(UNIT = 480, FILE = 'eps_ipa_eh.dat')
        DO iq = 1, nqtotf ! loop over q-points
        READ(480, *) qcart(:), q2
        DO nu = 1, nmodes
          READ(480, *) epsil_ipa(nu, iq)
        ENDDO
        ENDDO
      CLOSE(480)
      ENDIF
      CALL mp_bcast(epsil_ipa, ionode_id, inter_pool_comm)
      !
    ENDIF
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE restart_tdbe_read   
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE set_initial_dist()
    !-----------------------------------------------------------------------
    !! Set up initial distribution of electron and phonon
    !-----------------------------------------------------------------------     
    USE kinds,            ONLY : DP
    USE input,            ONLY : init_type_tdbe, init_sigma_tdbe, carr_dyn_tdbe, &
                                 ef_c_tdbe, ef_v_tdbe, temp_el_tdbe, temp_ph_tdbe, eps_acoustic,    &
                                 phph_tdbe, lscreen_tdbe, scr_typ, assume_metal
    USE ep_constants,     ONLY : kelvin2Ry, one, two, zero, ryd2ev
    USE modes,            ONLY : nmodes
    USE io_global,        ONLY : stdout, ionode_id, meta_ionode
    USE mp_global,        ONLY : inter_pool_comm, world_comm,npool, my_pool_id
    USE global_var,       ONLY : wf, nqtotf, nktotf, nbndfst, upper_bnd, lower_bnd
    USE tdbe_common,      ONLY : tstart, enk_all, ef0, felec, nphon, nphon_pre
    USE mp,               ONLY : mp_bcast, mp_barrier, mp_sum
    !
    IMPLICIT NONE
    !
    INTEGER :: ierr
    !! Error status
    INTEGER :: ios
    !! INTEGER variable for I/O control
    REAL(KIND = DP) :: inv_temp_el, inv_temp_ph
    !! Inverse temperature for carriers and phonon
    INTEGER :: ik, ibnd, iq, nu
    !! Counter on k, q grid, band index and phonon
    !! branch index
    REAL(KIND = DP), PARAMETER :: sqrtpi = 1.77245385090551602729_DP
    REAL(KIND = DP), EXTERNAL :: wgauss
    !! Fermi-Dirac distribution function (when -99)
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! The derivative of wgauss:  an approximation to the delta function
    !! invert temperature of electron and phonon, used for generating initial distribution


    tstart = zero
    WRITE(stdout, '(5x, a, f8.3, a/)' ) 'Initial distribution of phonon is simulated with BE distribution at ', & 
                                     temp_ph_tdbe, ' K temperarue'
    IF (init_type_tdbe == 'FD') THEN
      inv_temp_el = one / (kelvin2Ry * temp_el_tdbe)
      WRITE(stdout, '(5x, a, f8.3, a)' ) 'Initial distribution of electrons is simulated with FD distribution at ', &
      temp_el_tdbe, ' K (temperarue)'
      IF (carr_dyn_tdbe == 1) THEN   
        WRITE(stdout,'(5x,a)') 'Real time dynamics for electron, initial distribution determined by FD distribution'
        WRITE(stdout,'(5x,a,f8.6,a)') 'chemical potential of electron is', ef_c_tdbe,' eV'
        DO ik = lower_bnd, upper_bnd
          DO ibnd = 1, nbndfst
            IF(enk_all(ibnd, ik) < ef0) CALL errore('set_initial_dist', 'Please keep only &
                                                    conduction bands in the energy window', 1)
            felec(ibnd,ik) = wgauss(-(enk_all(ibnd, ik) - (ef_c_tdbe / ryd2ev)) * inv_temp_el, -99)
          ENDDO
        ENDDO
      ELSEIF (carr_dyn_tdbe == 2) THEN
        WRITE(stdout,'(5x,a)') 'Real time dynamics for holes, initial distribution determined by FD distribution'
        WRITE(stdout,'(5x,a,f8.6,a)') 'chemical potential of hole is', ef_v_tdbe,' eV'
        DO ik = lower_bnd, upper_bnd
          DO ibnd = 1, nbndfst
            IF(enk_all(ibnd, ik) > ef0) CALL errore('set_initial_dist', 'Please keep only &
                                                  valence bands in the energy window', 1)
            felec(ibnd,ik) = wgauss(-(enk_all(ibnd, ik) - (ef_v_tdbe / ryd2ev)) * inv_temp_el,-99)
          ENDDO
        ENDDO
      ELSEIF (carr_dyn_tdbe  == 3) THEN
        WRITE(stdout,'(5x,a)') 'Real time dynamics for electron and holes, initial distribution determined &
                                by FD distribution'
        WRITE(stdout,'(5x,a,f10.6,a)') 'chemical potential of electron is', ef_c_tdbe,' eV'
        WRITE(stdout,'(5x,a,f10.6,a)') 'chemical potential of hole is', ef_v_tdbe,' eV'
        DO ik = lower_bnd, upper_bnd
          DO ibnd = 1, nbndfst
            IF(enk_all(ibnd, ik) > ef0) THEN
              felec(ibnd, ik) = wgauss(-(enk_all(ibnd, ik) - (ef_c_tdbe / ryd2ev)) * inv_temp_el, -99)
            ELSE
              felec(ibnd, ik) = wgauss(-(enk_all(ibnd, ik) - (ef_v_tdbe / ryd2ev)) * inv_temp_el, -99)
            ENDIF
          ENDDO
        ENDDO
      ELSE
        CALL errore('set_initial_dist', 'carr_dyn_tdbe could only be 1, 2 or 3', 1)
      ENDIF
    ELSEIF (init_type_tdbe == 'gaussian') THEN   
      IF (carr_dyn_tdbe == 1) THEN
        WRITE(stdout,'(5x,a)') 'Real time dynamics for electron, initial distribtuion given by gaussian function.'
        WRITE(stdout,'(5x,a,f10.6,a)') 'initial energy of hot electron is', ef_c_tdbe,' eV'
        WRITE(stdout,'(5x,a,f10.6,a)') 'gaussian spread for initial distribution is', init_sigma_tdbe, ' eV'
        DO ik = lower_bnd, upper_bnd
          DO ibnd = 1, nbndfst
            IF(enk_all(ibnd, ik) < ef0) CALL errore('set_initial_dist', &
              'Please keep only conduction bands in the energy window', 1)
            felec(ibnd,ik) = sqrtpi * w0gauss((enk_all(ibnd, ik) - (ef_c_tdbe / ryd2ev)) / &
                             (init_sigma_tdbe / ryd2ev),0)
          ENDDO
        ENDDO
      ELSEIF (carr_dyn_tdbe == 2) THEN
        WRITE(stdout,'(5x,a)') 'Real time dynamics for holes, initial distribtuion given by gaussian function.'
        WRITE(stdout,'(5x,a,f10.6,a)') 'initial energy of hot hole is', ef_v_tdbe,' eV'
        WRITE(stdout,'(5x,a,f10.6,a)') 'gaussian spread for initial distribution is', init_sigma_tdbe,' eV'
        DO ik = lower_bnd, upper_bnd
          DO ibnd = 1, nbndfst
            IF(enk_all(ibnd, ik) > ef0) CALL errore('set_initial_dist', &
              'Please keep only valence bands in the energy window', 1)
             felec(ibnd, ik) = 1.d0 - sqrtpi * w0gauss((enk_all(ibnd, ik) -  &
                              (ef_v_tdbe / ryd2ev)) / (init_sigma_tdbe / ryd2ev), 0)
          ENDDO
        ENDDO 
      ELSEIF (carr_dyn_tdbe  == 3) THEN
        WRITE(stdout,'(5x,a)') 'Real time dynamics for electron and hole, initial distribtuion &
                                given by gaussian function.'
        WRITE(stdout,'(5x,a,f10.6,a)') 'initial energy of hot electron is', ef_c_tdbe,' eV'
        WRITE(stdout,'(5x,a,f10.6,a)') 'initial energy of hot hole is', ef_v_tdbe,' eV'
        WRITE(stdout,'(5x,a,f10.6,a)') 'gaussian spread for initial distribution is', init_sigma_tdbe,' eV'
        DO ik = lower_bnd, upper_bnd
          DO ibnd = 1, nbndfst
            IF(enk_all(ibnd, ik) > ef0) THEN
              felec(ibnd, ik) = sqrtpi * w0gauss((enk_all(ibnd, ik) - (ef_c_tdbe / ryd2ev)) / &
                               (init_sigma_tdbe / ryd2ev), 0)
            ELSE
              felec(ibnd, ik) = 1.d0 - sqrtpi * w0gauss((enk_all(ibnd, ik) - (ef_v_tdbe / ryd2ev)) &
                               / (init_sigma_tdbe / ryd2ev), 0)
            ENDIF
          ENDDO
        ENDDO
      ELSE
        CALL errore('set_initial_dist', 'carr_dyn_tdbe could only be 1, 2 or 3', 1)
      ENDIF
    ELSE
      CALL errore('set_initial_dist', 'init_type_tdbe is either FD or gaussian', 1)
    ENDIF
    CALL mp_sum(felec, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    inv_temp_ph = one / (kelvin2Ry * temp_ph_tdbe)
    DO iq = 1, nqtotf
      DO nu = 1 , nmodes
        !IF (wf(nu,iq) > zero) THEN
        IF (wf(nu, iq) > eps_acoustic) THEN
          nphon(nu, iq) = wgauss(-wf(nu, iq) * inv_temp_ph, -99)
          nphon(nu, iq) = nphon(nu, iq) / (one - two * nphon(nu, iq))
        ENDIF
      ENDDO
    ENDDO
    IF (phph_tdbe) THEN
      nphon_pre(:, :) = nphon(:, :)
    ENDIF
    !
    !! Calculate the screening induced by the electron-hole pairs
    IF (lscreen_tdbe .AND. scr_typ > 1 .AND. (.NOT. assume_metal)) CALL eps_ipa_eh(ef0)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE set_initial_dist
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE write_dist(initial, t_in_fs)
    !-----------------------------------------------------------------------
    !! Write the distribution of electrons and phonons 
    !-----------------------------------------------------------------------
    USE kinds,            ONLY : DP
    USE io_var,           ONLY : iufilfelec, iufilnphon
    USE io_global,        ONLY : stdout, ionode_id, meta_ionode
    USE global_var,       ONLY : nqtotf, nktotf, nbndfst, upper_bnd, lower_bnd, &
                                 wf
    USE tdbe_common,      ONLY : enk_all, ef0, felec, nphon, nphon_pre
    USE ep_constants,     ONLY : ryd2ev, ryd2mev
    USE modes,            ONLY : nmodes
    USE mp_global,        ONLY : my_pool_id, inter_pool_comm
    !
    IMPLICIT NONE
    !
    INTEGER :: ierr
    !! Error status
    INTEGER :: ios 
    !! INTEGER variable for I/O control
    LOGICAL,INTENT(IN) :: initial
    REAL(KIND = DP), INTENT(IN) :: t_in_fs
    INTEGER :: ik , iq
    !! Counter for k and q points
    INTEGER :: ibnd, nu
    !! Counter for band and phonon branch indexes

    IF (my_pool_id == ionode_id) THEN
      IF (initial) THEN
        OPEN(UNIT = iufilfelec, FILE = 'felec',STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
        IF (ios /= 0) CALL errore('tdbe', 'error opening file ' // 'felec', iufilfelec)
        WRITE(iufilfelec, '(2x, a)')'Initial Distribution of electrons'
        DO ik = 1,  nktotf
          DO ibnd = 1, nbndfst
            WRITE(iufilfelec,*) felec(ibnd, ik) , enk_all(ibnd, ik) * ryd2ev
          ENDDO
        ENDDO
        CLOSE(iufilfelec) 
        OPEN(UNIT = iufilnphon, FILE = 'nphon',STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)  
        IF (ios /= 0) CALL errore('tdbe', 'error opening file' // 'nphon', iufilnphon)
        WRITE(iufilnphon, '(2x, a)')'Initial Distribution of phonons'
        DO iq = 1, nqtotf
        DO nu = 1 , nmodes
          WRITE(iufilnphon,*) nphon(nu, iq),  wf(nu, iq) * ryd2ev
        ENDDO
        ENDDO  
        CLOSE(iufilnphon)     
      ELSE
        OPEN(UNIT = iufilfelec, FILE = 'felec',STATUS = 'unknown', POSITION = 'append', FORM = 'formatted', IOSTAT = ios)
        IF (ios /= 0) CALL errore('tdbe', 'error opening file ' // 'felec', iufilfelec)
        WRITE(iufilfelec, '(5x, a, f15.6 , a)' ) 'Distribution of electron at ', t_in_fs, ' fs'
        DO ik = 1,  nktotf
          DO ibnd = 1, nbndfst
            WRITE(iufilfelec, *) felec(ibnd, ik)
          ENDDO
        ENDDO
        CLOSE(iufilfelec) 
        OPEN(UNIT = iufilnphon, FILE = 'nphon',STATUS = 'unknown', POSITION = 'append', FORM = 'formatted', IOSTAT = ios)
        IF (ios /= 0) CALL errore('tdbe', 'error opening file ' // 'nphon', iufilnphon)
        WRITE(iufilnphon, '(5x, a, f15.6 , a)' ) 'Distribution of phonon at ', t_in_fs, ' fs'
        DO iq = 1, nqtotf
        DO nu = 1 , nmodes
          WRITE(iufilnphon, *) nphon(nu, iq)
        ENDDO
        ENDDO
        CLOSE(iufilnphon) 
      ENDIF  
      CLOSE(iufilfelec) 
      !
    ENDIF
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE write_dist
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE  write_points()
    !-----------------------------------------------------------------------
    !! write the k and q points as well as eigenvalues 
    !-----------------------------------------------------------------------
    USE kinds,            ONLY : DP
    USE io_var,           ONLY : iukmeshf, iuqmeshf 
    USE io_global,        ONLY : stdout, ionode_id, meta_ionode
    USE global_var,       ONLY : nqtotf, nktotf, nbndfst,wf, xqf
    USE tdbe_common,      ONLY : enk_all, ef0, felec, nphon, nphon_pre, &
                                 xkf_all
    USE ep_constants,     ONLY : ryd2ev
    USE modes,            ONLY : nmodes
    USE mp_global,        ONLY : my_pool_id
    !
    IMPLICIT NONE
    !
    INTEGER :: ierr 
    !! Error status
    INTEGER :: ios 
    !! INTEGER variable for I/O control
    INTEGER :: ik, iq
    !! Counter on  k and q point grid
    INTEGER :: ibnd, nu
    !! Counter on band and mode
    IF (my_pool_id == ionode_id) THEN
      OPEN(UNIT = iukmeshf, FILE = 'kpoints',STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
      IF (ios /= 0) CALL errore('tdbe', 'error opening file kpoints ', iukmeshf)
      DO ik = 1, nktotf
       WRITE(iukmeshf, *) xkf_all(:, ik)
        DO ibnd = 1, nbndfst
          WRITE(iukmeshf,*) enk_all(ibnd, ik) * ryd2ev
        ENDDO
      ENDDO
      CLOSE(iukmeshf)
      !
      OPEN(UNIT = iuqmeshf, FILE = 'qpoints',STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
      IF (ios /= 0) CALL errore('tdbe', 'error opening file qpoints', iuqmeshf)
      !WRITE(iuqmeshf,'(/2x,a/)') 'Phonon frequencies (in eV) on the fine q-mesh'
      DO iq = 1, nqtotf
        WRITE(iuqmeshf, *) xqf(:, iq)
        DO nu = 1, nmodes
          WRITE(iuqmeshf, *) wf(nu, iq) * ryd2ev
        ENDDO
      ENDDO
      CLOSE(iuqmeshf)
    ENDIF
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE write_points
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE wanninterp()
    !-----------------------------------------------------------------------
    !! Wannier interpolate electron and phonon eigenvalues
    !-----------------------------------------------------------------------
    USE kinds,            ONLY : DP                               
    USE global_var,       ONLY : wf, wqf, xqf, xkf, nqtotf, nqf, nkf, nkqf, &
                                 nktotf, nkqtotf, etf, chw, ibndmax, area,  &
                                 qrpl, ibndmin, nbndfst,     lower_bnd,     &
                                 upper_bnd, nbndskip, wkf, rdw, wscache,    &
                                 adapt_smearing, zstar, crrw, cvmew,        &
                                 bztoibz, s_bztoibz
    USE input,            ONLY : nkf1, nkf2, nkf3, nqf1, nqf2, nqf3,        &
                                 eps_acoustic, nbndsub,  fsthick, dt_tdbe,  & 
                                 use_ws, nqc1, nqc2, nqc3, nkc1, nkc2, nkc3,&
                                 efermi_read, fermi_energy, scissor, vme,   &
                                 restart_tdbe, lifc, phph_tdbe,     &
                                 mp_mesh_k, fixsym                               
    USE tdbe_common,      ONLY : enk_all, wnuq_all, wkf_all, xkf_all,       &
                                 xqf_all, wqf_all, vnk_all, vnuq_all,       &
                                 totq, selecq, uf_all, ef0, unk_all,        &
                                 Amnk_all
    USE ep_constants,     ONLY : zero, czero, eps6, eps4, ryd2ev, ci, one,  &
                                 two, ryd2mev, twopi, ryd2ghz , bohr2ang     
    USE noncollin_module, ONLY : noncolin                        
    USE wannier2bloch,    ONLY : hamwan2bloch, dynwan2bloch,vmewan2bloch,   &
                                 vmewan2blochp, dynifc2blochf, rrwan2bloch 
    USE io_global,        ONLY : stdout, ionode_id, meta_ionode
    USE mp_global,        ONLY : inter_pool_comm, world_comm, npool,        &
                                 my_pool_id
    USE parallelism,      ONLY : fkbounds
    USE utilities,        ONLY : fermiwindow 
    USE mp,               ONLY : mp_bcast, mp_barrier, mp_sum, mp_max, mp_min
    USE cell_base,        ONLY : at, bg, alat, omega
    USE ions_base,        ONLY : nat, amass, ityp, tau
    USE symm_base,        ONLY : nrot,s, t_rev, time_reversal, set_sym_bl,  &
                                 find_sym, nsym
    USE symmetry,         ONLY : kpoints_time_reversal_init, s_k, nsym_k, t_rev_k
    USE modes,            ONLY : nmodes
    USE wigner,           ONLY : wigner_seitz_wrap, wigner_divide_ndegen
    USE bzgrid,           ONLY : loadqmesh_serial, loadkmesh_para, kpoint_grid_epw
    USE pwcom,            ONLY : ef,nbnd,nelec
    USE low_lvl,          ONLY : fix_sym   
    USE input,            ONLY : lscreen_tdbe, scr_typ, assume_metal                      
    !
    IMPLICIT NONE
    !
    INTEGER :: ierr
    !! Error status
    INTEGER :: ios
    !! INTEGER variable for I/O control
    INTEGER :: ik, ikk
    !! Counter on k points
    INTEGER :: ibnd, jbnd
    !! Counter on band index
    INTEGER :: iq 
    !! Counter on q points
    INTEGER :: mu
    !! counter on mode
    INTEGER :: nu
    !! counter on mode
    INTEGER :: ir
    !! Counter for WS loop
    INTEGER :: iw
    !! Counter on bands when use_ws == .TRUE.
    INTEGER :: iw2
    !! Counter on bands when use_ws == .TRUE.
    INTEGER :: nrr_k
    !! Number of WS points for electrons
    INTEGER :: nrr_q
    !! Number of WS points for phonons
    INTEGER :: dims
    !! Dims is either nbndsub if use_ws or 1 if not
    INTEGER :: dims2
    !! Dims is either nat if use_ws or 1 if not
    REAL(KIND = DP) :: dummy(3)
    !! Dummy variable
    REAL(KIND = DP) :: xxk(3)
    !! Current k-point
    REAL(KIND = DP) :: xxq(3), qcart(3), zbg(3)
    !! Current q-point
    LOGICAL :: mp_mesh_k_tmp
    !! mask the flag mp_mesh_k 
    REAL(KIND = DP) :: ediff_max
    !! 
    REAL(DP) ::  mdum(3,nat)
    !! 
    INTEGER :: i, j ,k
    INTEGER :: bztoibz_local(nkf1 * nkf2 * nkf3)
    !! BZ to IBZ mapping
    INTEGER :: s_bztoibz_local(nkf1 * nkf2 * nkf3)
    !! symmetry matrix for each k-point from the full BZ
    INTEGER :: ikibz, iqibz
    !! Index of irrducible k point 
    INTEGER :: ns
    !! Index on symmetry operations
    INTEGER :: nkirr
    !! local electron and phonon velocity
    REAL(KIND = DP) :: vmek_av(3), vmeq_av(3), vdiff(3)
    REAL(KIND = DP) :: rws(0:3, 200)
    !! Real-space wigner-Seitz vectors
    REAL(KIND = DP) :: atws(3, 3)
    !! Maximum vector: at*nq 
    INTEGER :: n_av
    !! Count the degeracy
    INTEGER :: nrws
    !! Number of real-space Wigner-Seitz
    LOGICAL :: already_skipped
    !! Skipping band during the Wannierization
    INTEGER, PARAMETER :: nrwsx = 200
    !! Maximum number of real-space Wigner-Seitz
    !! Number of real-space Wigner-Seitz
    INTEGER, ALLOCATABLE :: irvec_k(:, :)
    !! INTEGER components of the ir-th Wigner-Seitz grid point in the basis
    !! of the lattice vectors for electrons
    INTEGER, ALLOCATABLE :: irvec_q(:, :)
    !! INTEGER components of the ir-th Wigner-Seitz grid point for phonons
    INTEGER, ALLOCATABLE :: irvec_g(:, :)
    !! INTEGER components of the ir-th Wigner-Seitz grid point for electron-phonon
    INTEGER, ALLOCATABLE :: ndegen_k(:, :, :)
    !! Wigner-Seitz number of degenerescence (weights) for the electrons grid
    INTEGER, ALLOCATABLE :: ndegen_q(:, :, :)
    !! Wigner-Seitz weights for the phonon grid that depend on atomic positions $R + \tau(nb) - \tau(na)$
    INTEGER, ALLOCATABLE :: ndegen_g(:, :, :)
    !! Wigner-Seitz weights for the electron-phonon grid that depend on
    !! atomic positions $R - \tau(na)$  
    REAL(KIND = DP), ALLOCATABLE :: w2(:)
    !! Square of phonon frequency
    REAL(KIND = DP), ALLOCATABLE :: wf_tmp(:,:)
    !! Save phonon frequencies for symmetrisation
    REAL(KIND = DP), ALLOCATABLE :: irvec_r(:, :)
    !! Wigner-Size supercell vectors, store in real instead of integer
    REAL(KIND = DP), ALLOCATABLE :: rdotk(:)
    !! $r\cdot k$
    REAL(KIND = DP), ALLOCATABLE :: wslen_k(:)
    !! real-space length for electrons, in units of alat
    REAL(KIND = DP), ALLOCATABLE :: wslen_q(:)
    !! real-space length for phonons, in units of alat
    REAL(KIND = DP), ALLOCATABLE :: wslen_g(:)
    !! real-space length for electron-phonons, in units of alat
    COMPLEX(KIND = DP), ALLOCATABLE :: cufkk(:, :)
    !! Rotation matrix, fine mesh, points k
    COMPLEX(KIND = DP), ALLOCATABLE :: cufkq(:, :)
    !! Rotation matrix, fine mesh, points k
    COMPLEX(KIND = DP), ALLOCATABLE :: cufkk_all(:, :, :)
    !! Save eigenvectors for all k points
    COMPLEX(KIND = DP), ALLOCATABLE :: rrf_all(:, :, :, :)
    !! Save Berry connection for all kpoints in the pool
    REAL(KIND = DP), ALLOCATABLE :: vnk(:,:,:)
    COMPLEX(KIND = DP), ALLOCATABLE :: vmef(:,:,:)
    !! velocity operator in wannier basis
    COMPLEX(KIND = DP), ALLOCATABLE :: vmefp(:, :, :)
    !! velocity operator in wannier basis
    COMPLEX(KIND = DP), ALLOCATABLE :: uf0(:, :)
    !! Rotation matrix for phonons
    COMPLEX(KIND = DP), ALLOCATABLE :: cfac(:)
    !! Used to store $e^{2\pi r \cdot k}$ exponential
    REAL(KIND = DP), ALLOCATABLE :: enk_all_nosym(:,:)
    !! electronic eigenvalues 
    REAL(KIND = DP), ALLOCATABLE :: w_centers(:, :)
    !! wannier centers
    !
    ALLOCATE(w_centers(3, nbndsub), STAT = ierr)
    IF (ierr /= 0) CALL errore('crystal_fmt_read', 'Error allocating w_centers', 1)
    CALL crystal_fmt_read(w_centers)
    CALL quadrupole_fmt_read()
  ! Center the WS at Gamma for electonic part, the phonon part and el-ph part
    IF (use_ws) THEN
      dims = nbndsub
      dims2 = nat 
      CALL wigner_seitz_wrap(nkc1, nkc2, nkc3, nqc1, nqc2, nqc3, irvec_k, irvec_q,irvec_g, &
                             ndegen_k, ndegen_q, ndegen_g, wslen_k, wslen_q, wslen_g,      &
                             w_centers, dims, tau, dims2 )
    ELSE
      dims  = 1
      dims2 = 1
      dummy(:) = (/0.0, 0.0, 0.0/)
      CALL wigner_seitz_wrap(nkc1, nkc2, nkc3, nqc1, nqc2, nqc3, irvec_k, irvec_q,irvec_g,&
                             ndegen_k, ndegen_q, ndegen_g, wslen_k, wslen_q,wslen_g,      &
                             dummy, dims, dummy,dims2 )
    ENDIF
    !! On exiting, ndegen_k contains the degeneracies of each pairs 
    !! of wannier centers while irvec_k contains the minimal communal sets of WS vectors.
    !!            ndegen_q contains the degeneracies of each pairs
    !! of atoms while irvec_q contains the minimal communal sets of WS vectors.
    !! deallocate since not used
    DEALLOCATE(irvec_g, STAT = ierr)
    IF (ierr /= 0) CALL errore('wanninterp', 'Error deallocating irvec_g', 1)
    DEALLOCATE(ndegen_g, STAT = ierr)
    IF (ierr /= 0) CALL errore('wanninterp', 'Error deallocating ndegen_g', 1)  
    DEALLOCATE(wslen_g, STAT = ierr)
    IF (ierr /= 0) CALL errore('wanninterp', 'Error deallocating wslen_g', 1)    
    !
    nrr_k = SIZE(irvec_k(1, :))
    nrr_q = SIZE(irvec_q(1, :))
    IF (use_ws)THEN
      WRITE(stdout, '(5x,a)' )    'Construct the Wigner-Seitz cell using Wannier centers and atomic positions '
      WRITE(stdout, '(5x,a,i8)' ) 'Number of WS vectors for electrons ', nrr_k
      WRITE(stdout, '(5x,a,i8)' ) 'Number of WS vectors for phonons ', nrr_q
    ELSE
      WRITE(stdout, '(5x,a)' )    'Use zone-centred Wigner-Seitz cells '
      WRITE(stdout, '(5x,a,i8)' ) 'Number of WS vectors for electrons ', nrr_k
      WRITE(stdout, '(5x,a,i8)' ) 'Number of WS vectors for phonons ', nrr_q
      WRITE(stdout, '(5x,a)' )    'Results may be improved by using use_ws == .TRUE. '
    ENDIF ! use_ws

    !! Read Hamiltonian in wannier basis
    CALL epwdata_fmt_read(nrr_k, nrr_q)
    CALL wigner_divide_ndegen(chw, 1, nbndsub, nrr_k, 1, ndegen_k, dims)
    CALL wigner_divide_ndegen(cvmew, 3, nbndsub, nrr_k, 1, ndegen_k, dims)
    CALL wigner_divide_ndegen(crrw, 3, nbndsub, nrr_k, 1, ndegen_k, dims)
    ef0 = ef
    WRITE(stdout,'(/5x,a,f10.6,a)') 'Fermi energy coarse grid = ', ef0 * ryd2ev, ' eV'
    !! mask the symmetry flag in order to obtain the full k grid
    mp_mesh_k_tmp = .FALSE.
    IF (mp_mesh_k) THEN
      mp_mesh_k_tmp = .TRUE.
      mp_mesh_k = .FALSE.
    ENDIF
    !! Generate the homogeneous k and q grid
    CALL loadqmesh_serial
    CALL loadkmesh_para
    !! Defines the total number of k-points
    nktotf = nkqtotf / 2
    IF (nqtotf /= nqf) CALL errore('wanninterp', 'uniform q grid must be loaded in sequential', 1)
    ALLOCATE(wf(nmodes, nqf), STAT = ierr) ! defined in global_var, global
    IF (ierr /= 0) CALL errore('wanninterp', 'Error allocating wnuq_all', 1)
    ALLOCATE(etf(nbndsub, nkqf), STAT = ierr) ! defined in global_var, global
    IF (ierr /= 0) CALL errore('wanninterp', 'Error allocating etf', 1)
    ALLOCATE(cufkk(nbndsub, nbndsub), STAT = ierr) ! defined in this subroutine, local
    IF (ierr /= 0) CALL errore('wanninterp', 'Error allocating cufkk', 1)
    ALLOCATE(cufkq(nbndsub, nbndsub), STAT = ierr) ! defined in this subroutine, local
    !
    IF (lscreen_tdbe .AND. scr_typ > 1 .AND. (.NOT. assume_metal)) THEN
      ALLOCATE(cufkk_all(nbndsub, nbndsub, nkf), STAT = ierr) 
      IF (ierr /= 0) CALL errore('wanninterp', 'Error allocating cufkk_all', 1)
      cufkk_all = czero
      IF (vme == 'wannier') THEN 
        ALLOCATE(rrf_all(3, nbndsub, nbndsub, nkf), STAT = ierr)
        IF (ierr /= 0) CALL errore('wanninterp', 'Error allocating rrf', 1)
        rrf_all(:, :, :, :) = czero
      ENDIF
    ENDIF
    !
    IF (adapt_smearing) THEN
      ALLOCATE(vnk(3, nbndsub, nkqf), STAT = ierr)
      IF (ierr /= 0) CALL errore('wanninterp', 'Error allocating vnk', 1)
      ALLOCATE(vmef(3, nbndsub, nbndsub), STAT = ierr)
      IF (ierr /= 0) CALL errore('wanninterp', 'Error allocating vmef', 1)
      ALLOCATE(vmefp(3, nmodes, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('wanninterp', 'Error allocating vmefp', 1)
      vmefp(:, :, :) = czero
      vnk(:, :, :) = zero
      vmef(:, :, :) = czero
    ENDIF
    !! To include ph-ph interaction, one needs to 
    !! save eigenvects for all q points
    IF (phph_tdbe) THEN
      ALLOCATE(uf_all(nqtotf, nmodes, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('wanninterp', 'Error allocating uf', 1)
    ELSE
      ALLOCATE(uf_all(1, 1, 1), STAT = ierr)
    ENDIF
    uf_all(:, :, :) = czero
    ALLOCATE(uf0(nmodes, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('wanninterp', 'Error allocating uf0', 1)
    uf0(:, :)            = czero
    wf(:, :)             = zero
    etf(:, :)            = zero
    cufkk(:, :)          = czero
    cufkq(:, :)          = czero  
    !! restore the original value of mp_mesh_k
    IF (mp_mesh_k_tmp) THEN
      mp_mesh_k = .TRUE.
    ENDIF

    !! Find out the symmetry operations, the map from BZ to IBZ and 
    !! the symmetry operation that links a BZ point to its IBZ friend
    bztoibz_local(:)   = 0
    s_bztoibz_local(:) = 0
    s(:, :, :) = 0 !! Symmetry in crystal axis with dim: 3,3,48
    CALL set_sym_bl()
    !
    ! Setup Bravais lattice symmetry
    WRITE(stdout,'(5x,a,i3)') "Symmetries of Bravais lattice: ", nrot
    !
    !! Setup crystal symmetry
    CALL find_sym(nat, tau, ityp, .FALSE., mdum)
    IF (fixsym) CALL fix_sym(.FALSE.)
    WRITE(stdout, '(5x, a, i3)') "Symmetries of crystal:         ", nsym
    !
    CALL kpoints_time_reversal_init()
    ! CALL kpoint_grid_epw(nsym, time_reversal, s, t_rev, nkf1, nkf2, &
    !                      nkf3, nkirr, bztoibz, s_bztoibz)
    IF (ALLOCATED(bztoibz)) THEN
      DEALLOCATE(bztoibz, STAT = ierr)
      IF (ierr /= 0) CALL errore('wanninterp', 'error deallocating bztoibz', 1)
    ENDIF
    IF (ALLOCATED(s_bztoibz)) THEN 
      DEALLOCATE(s_bztoibz, STAT = ierr)
      IF (ierr /= 0) CALL errore('wanninterp', 'error deallocating s_bztoibz', 1)
    ENDIF
    !
    CALL kpoint_grid_epw(nsym_k, .FALSE., s_k, t_rev_k, nkf1, nkf2, &
                         nkf3, nkirr, .FALSE.)
    bztoibz_local = bztoibz
    s_bztoibz_local = s_bztoibz
    DEALLOCATE(bztoibz, STAT = ierr)
    IF (ierr /= 0) CALL errore('wanninterp', 'error deallocating bztoibz', 1)
    DEALLOCATE(s_bztoibz, STAT = ierr)
    IF (ierr /= 0) CALL errore('wanninterp', 'error deallocating s_bztoibz', 1)
    !
    write(stdout,'(5x,a,i5)') "Number of irreducible k points is: ", nkirr
    !
    !! exponential factor and 
    ALLOCATE(cfac(nrr_k), STAT = ierr)
    IF (ierr /= 0) CALL errore('wanninterp', 'Error allocating cfac', 1)
    ALLOCATE(rdotk(nrr_k), STAT = ierr)
    IF (ierr /= 0) CALL errore('wanninterp', 'Error allocating rdotk', 1)
    ! This is simply because dgemv take only real number (not integer)
    ALLOCATE(irvec_r(3, nrr_k), STAT = ierr)
    IF (ierr /= 0) CALL errore('wanninterp', 'Error allocating irvec_r', 1)
    irvec_r = REAL(irvec_k, KIND = DP)
    !
    !! Zeroing everything 
    cfac(:)  = czero
    rdotk(:)       = zero
    !! Hamiltonian : Wannier -> Bloch 
    ! ------------------------------------------------------
    !
    xxq = zero
    !
    !! nkqf is the number of kpoints in the pool

    DO ik = 1, nkqf
      !
      xxk = xkf(:, ik)
      !
      IF (2 * (ik / 2) == ik) THEN
        !
        !  this is a k+q point : redefine as xkf (:, ik-1) + xxq
        !
        CALL cryst_to_cart(1, xxq, at, -1)
        xxk = xkf(:, ik - 1) + xxq
        CALL cryst_to_cart(1, xxq, bg, 1)
        !
      ENDIF
      !
      ! SP: Compute the cfac only once here since the same are use in both hamwan2bloch and dmewan2bloch
      ! + optimize the 2\pi r\cdot k with Blas
      CALL DGEMV('t', 3, nrr_k, twopi, irvec_r, 3, xxk, 1, 0.0_DP, rdotk, 1)
      !
      cfac(:) = EXP(ci * rdotk(:))
      !
      CALL hamwan2bloch(nbndsub, nrr_k, cufkk, etf(:, ik), chw, cfac)
      IF (lscreen_tdbe .AND. scr_typ > 1 .AND. (2 * (ik / 2) /= ik) .AND. (.NOT. assume_metal)) THEN
        ikk = ik / 2 + 1
        cufkk_all(:, :, ikk) = cufkk(:, :)
        IF (vme == 'wannier') THEN          
          CALL rrwan2bloch(nbndsub, nrr_k, cfac, rrf_all(:, :, :, ikk))
        ENDIF
      ENDIF
      !! If we use adaptive broadening parameter in delta function, we need 
      !! additionally the band velocity
      IF (adapt_smearing) THEN
        CALL vmewan2bloch(nbndsub, nrr_k, irvec_k, cufkk, vmef, etf(:, ik), etf(:, ik), chw, cfac)
        DO ibnd = 1, nbndsub
          vnk(:, ibnd, ik) = REAL(vmef(:, ibnd, ibnd))
        ENDDO
      ENDIF
    ENDDO
    !
    IF (ABS(scissor) > eps6) THEN
      WRITE(stdout, '(5x,"Applying a scissor shift of ",f9.5," eV to the CB ")' ) scissor * ryd2ev
      DO  ik = 1, nkf
        DO ibnd = 1, nbndsub
          IF (etf(ibnd, 2 * ik - 1) > ef0) THEN
            etf(ibnd, 2 * ik - 1) = etf(ibnd, 2 * ik - 1) + scissor
          ENDIF
        ENDDO
      ENDDO
    ENDIF
    !
    IF (efermi_read) THEN
      !
      ef = fermi_energy
      WRITE(stdout,'(/5x, a)') REPEAT('=',67)
      WRITE(stdout, '(/5x, a, f10.6, a)') 'Fermi energy is read from the input file: Ef = ', ef * ryd2ev, ' eV'
      !! SP: even when reading from input the number of electron needs to be correct
      already_skipped = .FALSE.
      IF (nbndskip > 0) THEN
        IF (.NOT. already_skipped) THEN
          IF (noncolin) THEN
            nelec = nelec - one * nbndskip
          ELSE
            nelec = nelec - two * nbndskip
          ENDIF
          already_skipped = .TRUE.
          WRITE(stdout, '(/5x,"Skipping the first ", i4, " bands:")') nbndskip
          WRITE(stdout, '(/5x,"The Fermi level will be determined with ", f9.5, " electrons")') nelec
        ENDIF
      ENDIF    
      WRITE(stdout,'(/5x,a)') REPEAT('=',67)
      WRITE(stdout, '(/5x,"There are ", f9.5, " electrons per unit cell")') nelec
      WRITE(stdout, '(/5x,"Total nb of band ", i4)') nbndsub
      !
    ENDIF
    !
    CALL fermiwindow()
    nbndfst = ibndmax - ibndmin + 1
    !
    !! Define it only once for the full run.
    CALL fkbounds(nktotf, lower_bnd, upper_bnd)
    !! electronic energies within the Fermi window
    ALLOCATE(enk_all_nosym(nbndfst, nktotf), STAT = ierr)
    IF (ierr /= 0) CALL errore('wanninterp', 'Error allocating enk_all', 1)
    ALLOCATE(enk_all(nbndfst, nktotf), STAT = ierr)
    IF (ierr /= 0) CALL errore('wanninterp', 'Error allocating enk_all', 1)
    ALLOCATE(wkf_all(nktotf), STAT = ierr)
    IF (ierr /= 0) CALL errore('wanninterp', 'Error allocating wkf_all', 1)
    ALLOCATE(xkf_all(3, nktotf), STAT = ierr)
    IF (ierr /= 0) CALL errore('wanninterp', 'Error allocating xkf_all', 1)
    xkf_all(:, :) = zero
    !
    IF (lscreen_tdbe .AND. scr_typ > 1 .AND. (.NOT. assume_metal)) THEN
      ALLOCATE(unk_all(nbndsub, nbndsub, nktotf), STAT = ierr) 
      IF (ierr /= 0) CALL errore('wanninterp', 'Error allocating unk_all', 1)
      IF (vme == 'wannier') THEN
        ALLOCATE(Amnk_all(3, nbndsub, nbndsub, nktotf), STAT = ierr)
        IF (ierr /= 0) CALL errore('wanninterp','Error deallocating Amnk_all',1)
        Amnk_all = czero
      ENDIF
      unk_all = czero
      DO ik = 1, nkf
        unk_all(:, :, ik + lower_bnd - 1) = cufkk_all(:, :, ik)
        IF (vme == 'wannier') THEN
          Amnk_all(:, :, :, ik + lower_bnd - 1) = rrf_all(:, :, :, ik)
        ENDIF
      ENDDO
      !
      CALL mp_sum(unk_all, inter_pool_comm)
      IF (vme == 'wannier') CALL mp_sum(Amnk_all, inter_pool_comm)
      CALL mp_barrier(inter_pool_comm)
      DEALLOCATE(cufkk_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('wanninterp','Error deallocating cufkk_all',1)
      IF (vme == 'wannier') THEN
        DEALLOCATE(rrf_all, STAT = ierr)
        IF (ierr /= 0) CALL errore('wanninterp','Error deallocating rrf_all',1)
      ENDIF
      !
      ! CALL eps_ipa( ef0)
    ENDIF 
    !
    DO ik = 1, nkf
      xkf_all(:,ik + lower_bnd - 1) = xkf(:,2 * ik - 1)
    ENDDO
    !
    CALL mp_sum(xkf_all, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    !
    ! ALLOCATE(xqf_all(3, nktotf), STAT = ierr)
    ! IF (ierr /= 0) CALL errore('wanninterp', 'Error allocating xqf_all', 1)
    ! xqf_all(:,:) = zero
    ! !
    ! DO iq = 1, nqtotf
    !   xqf_all(:,iq) = xqf(iq)
    ! ENDDO
    ! 
    enk_all_nosym(:, :) = zero
    enk_all(:, :) = zero
    wkf_all(:) = zero
    !
    DO ik = 1, nkf
      wkf_all(ik + lower_bnd - 1) = one / DBLE(nktotf)
      DO ibnd = 1, nbndfst
        enk_all_nosym(ibnd, ik + lower_bnd - 1) = etf(ibnd + ibndmin - 1, 2 * ik - 1)
      ENDDO
    ENDDO
    !
    IF (noncolin) wkf_all(:) = wkf_all(:) / two
    CALL mp_sum(enk_all_nosym, inter_pool_comm)
    CALL mp_sum(wkf_all, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    ! IF (my_pool_id == ionode_id) then 
    !   open(unit=566,file='enk_all_nosym')
    !   do ik = 1, nktotf
    !     do ibnd = 1, nbndfst
    !       write(566,*) enk_all_nosym(ibnd,ik)
    !     enddo
    !   enddo
    !   close(566)
    ! ENDIF
    ediff_max = zero! Ry
    !! The electron energies are symmetrized
    DO ik = 1, nkf
      ikibz = bztoibz_local(ik + lower_bnd - 1)
      DO ibnd = 1, nbndfst
        ediff_max = ABS(enk_all_nosym(ibnd, ik + lower_bnd - 1)-enk_all_nosym(ibnd, ikibz))
        IF (ediff_max > one / ryd2ev) THEN
          WRITE(stdout,'(5x, a)') 'waring: the difference between enk and enk_irr &
          are greater than 1.0 eV.'
        ENDIF
        enk_all(ibnd, ik + lower_bnd-1) = enk_all_nosym(ibnd, ikibz)
      ENDDO
    ENDDO
    !
    CALL mp_sum(enk_all, inter_pool_comm)
    DEALLOCATE(enk_all_nosym, STAT = ierr)
    IF (ierr /= 0) CALL errore('wanninterp','Error deallocating enk_all_nosym',1)
    !
    IF (adapt_smearing) THEN
      ALLOCATE(vnk_all(3, nbndfst, nktotf), STAT = ierr)
      IF (ierr /= 0) CALL errore('wanninterp','Error allocating vnk_all',1)
      vnk_all(:,:,:) = zero
      DO ik = 1, nkf
        DO ibnd = 1, nbndfst
          vmek_av(:) =  zero
          n_av = 0
          DO jbnd = 1, nbndfst
            IF (ABS(etf(ibnd + ibndmin - 1,2 * ik - 1)-           &
                etf(jbnd + ibndmin - 1,2 * ik - 1)) < eps4) THEN
              n_av = n_av + 1
              vmek_av(:) = vmek_av(:) + vnk(:, jbnd + ibndmin - 1,2 * ik -1)
            ENDIF
          ENDDO
          vnk_all(:,ibnd,ik + lower_bnd - 1) = vmek_av(:) / FLOAT(n_av)
          !vnk_all(:,ibnd,ik+lower_bnd-1)=vnk(:,ibnd+ibndmin-1,2*ik-1)
        ENDDO
      ENDDO
      CALL mp_sum(vnk_all, inter_pool_comm)
      DEALLOCATE(vnk, STAT = ierr)
      IF (ierr /= 0) CALL errore('wanninterp','Error deallocating vnk',1)
      DEALLOCATE(vmef, STAT = ierr)
      IF (ierr /= 0) CALL errore('wanninterp','Error deallocating vmef',1)
      ALLOCATE(vnuq_all(3, nmodes, nqtotf), STAT = ierr)
      IF (ierr /= 0) CALL errore('wanninterp','Error allocating vnuq_all',1)
      vnuq_all(:,:,:)=zero
      ! ALLOCATE(vmeq(3,nmodes), STAT = ierr)
      ! IF (ierr /= 0) CALL errore('wanninterp','Error deallocating vmeq',1)
      ! vmeq(:,:) = zero
    ELSE
      ALLOCATE(vnk_all(1, 1, 1), STAT = ierr)
      IF (ierr /= 0) CALL errore('wanninterp','Error allocating vnk_all',1)    
      ALLOCATE(vnuq_all(1, 1, 1), STAT = ierr)
      IF (ierr /= 0) CALL errore('wanninterp','Error allocating vnuq_all',1)
      vnk_all = zero 
      vnuq_all = zero
    ENDIF
    !
    totq = 0
    IF (restart_tdbe) THEN
      CALL qwindow_tdbe(.true.)
    ELSE
      CALL qwindow_tdbe(.false.)
    ENDIF
    WRITE(stdout, '(5x, a, i8, a)') 'We need to compute ', totq, ' q-points'
    WRITE(stdout, '(5x, a)')' '
    !
    ALLOCATE(w2(3 * nat), STAT = ierr)
    IF (ierr /= 0) CALL errore('wanninterp','Error allocating w2', 1)
    w2(:) = zero
    !
    ! In the case of crystal ASR
    IF (lifc) THEN
      !
      ! build the WS cell corresponding to the force constant grid
      !
      atws(:, 1) = at(:, 1) * DBLE(nqc1)
      atws(:, 2) = at(:, 2) * DBLE(nqc2)
      atws(:, 3) = at(:, 3) * DBLE(nqc3)
      rws(:, :)  = zero
      nrws       = 0
      ! initialize WS r-vectors
      CALL wsinit(rws, nrwsx, nrws, atws)
      ALLOCATE(wscache(-2 * nqc3:2 * nqc3, -2 * nqc2:2 * nqc2, -2 * nqc1:2 * nqc1, nat, nat), STAT = ierr)
      IF (ierr /= 0) CALL errore('wanninterp', 'Error allocating wscache', 1)
      wscache(:, :, :, :, :) = zero
    ELSE
      atws(:, :) = zero
      rws(:, :)  = zero
      nrws       = 0
    ENDIF
    !
    CALL mp_barrier(inter_pool_comm) 
!  DO iqq = 1, totq
    DO iq = 1 , nqtotf
      xxq = xqf(:, iq)
      IF (.NOT. lifc) THEN
        CALL dynwan2bloch(nmodes, nrr_q, irvec_q, ndegen_q, xxq, uf0, w2)
      ELSE 
        CALL dynifc2blochf(nmodes, rws, nrws, xxq, uf0, w2)
      ENDIF
      IF (phph_tdbe) uf_all(iq, :, :) = uf0
      DO nu = 1, nmodes
        IF (w2(nu) > 0.d0) THEN
          wf(nu, iq) =  SQRT(ABS(w2(nu)))
        ELSE
          wf(nu, iq) = -SQRT(ABS(w2(nu)))
        ENDIF
      ENDDO
      IF (adapt_smearing) THEN
        CALL vmewan2blochp(xxq, nmodes, nrr_q, irvec_q, ndegen_q, uf0, vmefp, wf(:, iq), rws, nrws)
        DO nu = 1, nmodes
          vmeq_av(:) = zero
          n_av = 0
          DO mu = 1, nmodes
            IF (ABS(wf(mu, iq) - wf(nu, iq)) < eps6) THEN
              n_av = n_av + 1
              vmeq_av(:) = vmeq_av(:) + REAL(vmefp(:, mu, mu))
            ENDIF
          ENDDO
          vnuq_all(:, nu, iq) = vmeq_av(:) / FLOAT(n_av)
          vnuq_all(:, nu, iq) = vnuq_all(:, nu, iq) / (two * wf(nu, iq))
        ENDDO
      ENDIF
    ENDDO
    !! The phonon frequencies need to be symmetrized
    ALLOCATE(wf_tmp(nmodes, nqf), STAT = ierr)
    wf_tmp =zero
    IF (ierr /= 0) CALL errore('wanninterp', 'Error allocating wf', 1)
    DO iq = 1, nqtotf
      iqibz = bztoibz_local(iq)
      DO nu = 1, nmodes
        ediff_max = ABS(wf(nu, iq)-wf(nu, iqibz))
        IF (ediff_max >  1E-2 / ryd2ev) THEN
          WRITE(stdout,'(5x, a)') 'waring: the difference between wnuq and wnuq_irr are greater than 10 meV, &
          check the dynamic matrix'
        ENDIF
        wf_tmp(nu, iq) = wf(nu, iqibz)
      ENDDO
    ENDDO
    wf(:,:) = wf_tmp(:,:)
    DEALLOCATE(wf_tmp, STAT = ierr)
    IF (ierr /= 0 ) CALL errore('wanninterp', 'Error allocating wf_tmp', 1)
    !
    WRITE(stdout,'(5x,a,f10.6)') 'threshold of phonon freqency(mev):',eps_acoustic * ryd2mev
    !
    WRITE(stdout, '(5x, a, f10.6, a)' ) 'Fermi Surface thickness = ', fsthick * ryd2ev, ' eV'
    !! Now deallocating
    IF (lifc) THEN
      DEALLOCATE(wscache, STAT = ierr)
      IF (ierr /= 0) CALL errore('wanninterp', 'Error deallocating wscache', 1)
    ENDIF    
    !
    IF (phph_tdbe) THEN
      IF (my_pool_id == ionode_id) THEN 
      OPEN(UNIT = 560, FILE = 'ph_freq', STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
      IF (ios /= 0) CALL errore('wanninterp', 'error opening file ph_freq',1)
      do iq = 1, nqtotf 
        do nu = 1, nmodes
          write(560,'(E22.14)') wf(nu, iq)
          do mu = 1, nmodes
            write(560,'(2E22.14)') uf_all(iq ,mu, nu)
          enddo
          !read(560,*)
        enddo
      enddo 
      close(560)
      ENDIF
    ENDIF
    !
    !! Now calculate the dielectric function
    IF (lscreen_tdbe .AND. scr_typ > 1 .AND. (.NOT. assume_metal)) CALL eps_ipa( ef0)
    !!
    !
    DEALLOCATE(chw, STAT = ierr)
    IF (ierr /= 0) CALL errore('wanninterp', 'Error deallocating chw', 1)
    DEALLOCATE(rdw, STAT = ierr)
    IF (ierr /= 0) CALL errore('wanninterp', 'Error deallocating rdw', 1)
    DEALLOCATE(cufkk, STAT = ierr)
    IF (ierr /= 0) CALL errore('wanninterp', 'Error deallocating cufkk', 1)
    DEALLOCATE(cufkq, STAT = ierr)
    IF (ierr /= 0) CALL errore('wanninterp', 'Error deallocating cufkq', 1)
    DEALLOCATE(uf0, STAT = ierr)
    IF (ierr /= 0) CALL errore('wanninterp', 'Error deallocating uf0', 1)
    DEALLOCATE(w2, STAT = ierr)
    IF (ierr /= 0) CALL errore('wanninterp', 'Error deallocating w2', 1)
    DEALLOCATE(cfac, STAT = ierr)
    IF (ierr /= 0) CALL errore('wanninterp', 'Error deallocating cfac', 1)
     DEALLOCATE(rdotk, STAT = ierr)
    IF (ierr /= 0) CALL errore('wanninterp', 'Error deallocating rdotk', 1)
    DEALLOCATE(irvec_r, STAT = ierr)
    IF (ierr /= 0) CALL errore('wanninterp', 'Error deallocating irvec_r', 1)
    DEALLOCATE(irvec_k, STAT = ierr)
    IF (ierr /= 0) CALL errore('wanninterp', 'Error deallocating irvec_k', 1)
    DEALLOCATE(irvec_q, STAT = ierr)
    IF (ierr /= 0) CALL errore('wanninterp', 'Error deallocating irvec_q', 1)
    DEALLOCATE(ndegen_k, STAT = ierr)
    IF (ierr /= 0) CALL errore('wanninterp', 'Error deallocating ndegen_k', 1)
    DEALLOCATE(ndegen_q, STAT = ierr)
    IF (ierr /= 0) CALL errore('wanninterp', 'Error deallocating ndegen_q', 1)
    DEALLOCATE(wslen_k, STAT = ierr)
    IF (ierr /= 0) CALL errore('wanninterp', 'Error deallocating wslen_k', 1)
    DEALLOCATE(wslen_q, STAT = ierr)
    IF (ierr /= 0) CALL errore('wanninterp', 'Error deallocating wslen_q', 1)
    DEALLOCATE(w_centers, STAT = ierr)
    IF (ierr /= 0) CALL errore('wanninterp','Error deallocating w_centers',1)
    IF (adapt_smearing) THEN
        DEALLOCATE(vmefp, STAT = ierr)
      IF (ierr /= 0) CALL errore('wanninterp','Error deallocating vmefp',1)
    ENDIF
    CALL mp_barrier(inter_pool_comm)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE wanninterp
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE eps_ipa( ef0)
    !-----------------------------------------------------------------------
    !! calculate the dielectric function of the updoped material
    !-----------------------------------------------------------------------
    USE kinds,            ONLY : DP
    USE global_var,       ONLY : wf, nkf, etf, wkf, nqtotf, lower_bnd, upper_bnd, &
                                 nbndskip, xqf, xkf, nktotf
    USE input,            ONLY : nbndsub, vme
    USE cell_base,        ONLY : alat, omega
    USE tdbe_common,      ONLY : unk_all, Amnk_all, epsil_ipa
    USE modes,            ONLY : nmodes
    USE io_global,        ONLY : stdout, ionode_id 
    USE mp,               ONLY : mp_barrier, mp_sum, mp_max
    USE mp_global,        ONLY : inter_pool_comm, npool, my_pool_id
    USE ep_constants,     ONLY : czero, eps5, zero, fpi, twopi, ci, cone
    USE longrange,        ONLY : compute_umn_f
    USE noncollin_module, ONLY : noncolin
    USE bzgrid,           ONLY : kpmq_map
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(IN) :: ef0 
    !! The Fermi energy of the undoped material
    COMPLEX(KIND = DP) :: cufkk(nbndsub, nbndsub)
    !! Rotation matrix, fine mesh, points k
    COMPLEX(KIND = DP) :: cufkq(nbndsub, nbndsub)
    !! Rotation matrix, fine mesh, points k+q
    COMPLEX(KIND = DP) :: cufkk_d(nbndsub, nbndsub)
    !! complex conjugate of cufkk
    COMPLEX(KIND = DP) :: bmatf(nbndsub, nbndsub)
    !! Electron wavefunction in wannier basis
    COMPLEX(KIND = DP) :: Amat(nbndsub, nbndsub)
    !! Temporary matrix saving A^W U_k, where A^W is the Berry 
    !! connection in wannier basis
    COMPLEX(KIND = DP) :: qdotA(nbndsub, nbndsub)
    !! Save i {\bf q} \cdot A
    REAL(KIND = DP) :: etf_all(nbndsub, nktotf)
    !! overlap U_k+q U_k^\dagger in smooth Bloch basis, fine mesh
    INTEGER :: ierr 
    !! error info
    INTEGER :: ios 
    !! io error info
    INTEGER :: ibnd, jbnd
    !! band index
    INTEGER :: imode
    !! phonon branch index
    INTEGER :: ivbm
    !! index of VBM
    INTEGER :: ik, iq 
    !! k and q index
    INTEGER :: ikq
    !! index of k+q
    REAL(KIND = DP) :: freq 
    !! phonon energy
    REAL(KIND = DP) :: xxk(3)
    !! Crystal coordinates of k-point
    REAL(KIND = DP) :: xxq(3), xxq_s(3)
    !! Crystal coordinates of q-point
    REAL(KIND = DP) :: qcart(3)
    !! Cartesian coordinates of q-point
    REAL(KIND = DP) :: q2
    !! norm of q point
    REAL(KIND = DP) :: fac
    !
    ivbm = 1
    !
    fac = 2.0 / omega  !! 2.0 accounts for the spin degeneracy
    IF (noncolin) fac = fac / 2.0 !! In case of noncolin, spin degeneracy is 1.0
    !
    fac = fac * 2.0 !! 2.0 = e^2, where e is the electron charge in Rydberg atomic unit
    !
    etf_all = zero
    DO ik = 1, nkf 
      etf_all(:, ik + lower_bnd - 1) = etf(:, 2 * ik - 1)
    ENDDO
    CALL mp_sum(etf_all, inter_pool_comm)
    !
    DO ik = 1, nkf 
      DO ibnd = 1, nbndsub
        IF (etf(ibnd, 2 * ik - 1) < ef0) ivbm = ibnd
      ENDDO
    ENDDO    
    !
    CALL mp_max(ivbm, inter_pool_comm)
    !
    WRITE(stdout, '(5x,"Index of maximum VB ",i6)' ) ivbm
    !
    ALLOCATE(epsil_ipa(nmodes, nqtotf), STAT = ierr)
    IF (ierr /= 0) CALL errore('epsil_rpa', 'Error allocating epsil_ipa', 1)
    !
    epsil_ipa = zero
    bmatf = czero
    !
    DO iq = 1, nqtotf
      xxq_s(:) = xqf(:, iq)
      xxq(:)   = xqf(:, iq)
      CALL shift_q(xxq_s, qcart)  ! shift q point to the 1st BZ
      q2 = (qcart(1)**2 + qcart(2)**2 + qcart(3)**2) * (twopi / alat) ** 2
      IF (ABS(q2) < eps5) CYCLE
      DO ik = 1, nkf
        xxk(:) = xkf(:, 2 * ik - 1)
        CALL kpmq_map(xxk, xxq, 1, ikq)
        cufkk(:, :) = unk_all(:, :, ik + lower_bnd - 1)
        cufkq(:, :) = unk_all(:, :, ikq)
        CALL compute_umn_f(nbndsub, cufkk, cufkq, bmatf)
        !! Now include the Berry connection term
        qdotA = czero
        Amat  = czero
        IF (vme == 'wannier') THEN
          qdotA(:, :) = qcart(1) * Amnk_all(1, :, :, ik + lower_bnd - 1) + &
                        qcart(2) * Amnk_all(2, :, :, ik + lower_bnd - 1) + &
                        qcart(3) * Amnk_all(3, :, :, ik + lower_bnd - 1)
          qdotA = qdotA * (twopi / alat)
          qdotA = ci * qdotA
          ! cufkk_d = conjg(transpose(cufkk))
          ! Amat  = matmul(qdotA,  cufkk_d)
          ! qdotA = matmul(cufkq, Amat)
          CALL ZGEMM('n', 'n', nbndsub, nbndsub, nbndsub, cone, cufkq, nbndsub, qdotA, nbndsub, czero, Amat, nbndsub)
          CALL ZGEMM('n', 'c', nbndsub, nbndsub, nbndsub, cone, Amat, nbndsub, cufkk, nbndsub, czero, qdotA, nbndsub)
        ENDIF
        !
        DO ibnd = 1, ivbm
          DO jbnd = ivbm + 1, nbndsub
            DO imode = 1, nmodes
              IF (wf(imode, iq) < 1e-5) THEN
                freq = 0.d0
              ELSE 
                freq = wf(imode, iq)
              ENDIF
              epsil_ipa(imode, iq) = epsil_ipa(imode, iq) - fac * fpi/ q2 / nktotf /   &
              (etf_all(jbnd, ikq) - etf_all(ibnd, ik + lower_bnd - 1) - freq) * &
              ABS(bmatf(jbnd, ibnd) + qdotA(jbnd, ibnd)) ** 2
            ENDDO
          ENDDO
        ENDDO
        DO ibnd = ivbm + 1, nbndsub
          DO jbnd = 1, ivbm
            DO imode = 1, nmodes
              IF (wf(imode, iq) < 1e-5) THEN
                freq = 0.d0
              ELSE 
                freq = wf(imode, iq)
              ENDIF
              epsil_ipa(imode, iq) = epsil_ipa(imode, iq) + fac * fpi/ q2 / nktotf /   &
              (etf_all(jbnd, ikq) - etf_all(ibnd, ik + lower_bnd - 1) - freq) * &
              ABS(bmatf(jbnd, ibnd) + qdotA(jbnd, ibnd)) ** 2
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !
    CALL mp_sum(epsil_ipa, inter_pool_comm)
    !
    epsil_ipa(:, :) = 1.d0 - epsil_ipa(:, :) 
    !! Start writing the dielectric constants
    IF (my_pool_id == ionode_id) THEN
      OPEN(UNIT = 480, FILE = 'eps_ipa.dat')
        DO iq = 1, nqtotf ! loop over q-points 
          xxq_s(:) = xqf(:, iq)
          CALL shift_q(xxq_s, qcart)  ! shift q point to the 1st BZ
          q2 = (qcart(1)**2 + qcart(2)**2 + qcart(3)**2) * (twopi / alat) ** 2
          WRITE(480, *) qcart(:), q2
          DO imode = 1, nmodes
            WRITE(480, *) epsil_ipa(imode, iq)
          ENDDO
        ENDDO
      CLOSE(480)
    ENDIF
    CALL mp_barrier(inter_pool_comm)
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE eps_ipa
    !-----------------------------------------------------------------------  
    !
    !-----------------------------------------------------------------------
    SUBROUTINE eps_ipa_eh( ef0)
    !-----------------------------------------------------------------------
    !! calculate the dielectric function of the updoped material
    !-----------------------------------------------------------------------
    USE kinds,            ONLY : DP
    USE global_var,       ONLY : wf, nkf, etf, wkf, nqtotf, lower_bnd, upper_bnd, &
                                 nbndskip, xqf, xkf, nktotf, ibndmin, nbndfst
    USE input,            ONLY : nbndsub, vme
    USE cell_base,        ONLY : alat, omega
    USE tdbe_common,      ONLY : unk_all, Amnk_all, epsil_ipa, enk_all, felec, uf_all
    USE modes,            ONLY : nmodes
    USE io_global,        ONLY : stdout, ionode_id 
    USE mp,               ONLY : mp_barrier, mp_sum, mp_max
    USE mp_global,        ONLY : inter_pool_comm, npool, my_pool_id
    USE ep_constants,     ONLY : czero, eps5, zero, fpi, twopi, ci, cone
    USE longrange,        ONLY : compute_umn_f
    USE noncollin_module, ONLY : noncolin
    USE bzgrid,           ONLY : kpmq_map
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(IN) :: ef0 
    !! The Fermi energy of the undoped material
    COMPLEX(KIND = DP) :: cufkk(nbndsub, nbndsub)
    !! Rotation matrix, fine mesh, points k
    COMPLEX(KIND = DP) :: cufkq(nbndsub, nbndsub)
    !! Rotation matrix, fine mesh, points k+q
    COMPLEX(KIND = DP) :: cufkk_d(nbndsub, nbndsub)
    !! complex conjugate of cufkk
    COMPLEX(KIND = DP) :: bmatf(nbndsub, nbndsub)
    !! Electron wavefunction in wannier basis
    COMPLEX(KIND = DP) :: Amat(nbndsub, nbndsub)
    !! Temporary matrix saving A^W U_k, where A^W is the Berry 
    !! connection in wannier basis
    COMPLEX(KIND = DP) :: qdotA(nbndsub, nbndsub)
    !! Save i {\bf q} \cdot A
    REAL(KIND = DP) :: etf_all(nbndsub, nktotf)
    !! overlap U_k+q U_k^\dagger in smooth Bloch basis, fine mesh
    REAL(KIND = DP) :: dchi(nmodes, nqtotf)
    !! \delta\chi -> independent particle polarizability contributed by free carriers 
    REAL(KIND = DP) :: fik, fjkq 
    !! electronic occupation in undoped semicoonductors
    REAL(KIND = DP) :: etrans
    !! etf_all(jbnd + ibndmin - 1, ikq) - etf_all(ibnd + ibndmin - 1, ik + lower_bnd - 1)
    INTEGER :: ierr 
    !! error info
    INTEGER :: ios 
    !! io error info
    INTEGER :: ibnd, jbnd
    !! band index
    INTEGER :: imode 
    !! phonon branch index
    INTEGER :: ivbm
    !! index of VBM
    INTEGER :: ik, iq 
    !! k and q index
    INTEGER :: ikq
    !! index of k+q
    REAL(KIND = DP) :: xxk(3)
    !! Crystal coordinates of k-point
    REAL(KIND = DP) :: xxq(3), xxq_s(3)
    !! Crystal coordinates of q-point
    REAL(KIND = DP) :: qcart(3)
    !! Cartesian coordinates of q-point
    REAL(KIND = DP) :: q2
    !! norm of q point
    REAL(KIND = DP) :: fac
    !!
    REAL(KIND = DP) :: freq 
    !! phonon energy
    !
    ivbm = 1
    !
    !
    fac = 2.0 * fpi / nktotf / omega !! 2.0 = e^2, where e is the electron charge in Rydberg atomic unit
    ! fpi = 1/epsilon_0
    !
    IF (.NOT. noncolin) fac = fac * 2.0 !! In case of noncolin, spin degeneracy is 1.0, otherwise is 2.0
    !
    etf_all = zero
    !
    DO ik = 1, nkf 
      etf_all(:, ik + lower_bnd - 1) = etf(:, 2 * ik - 1)
    ENDDO
    CALL mp_sum(etf_all, inter_pool_comm)
    !
!     DO ik = 1, nkf 
!       DO ibnd = 1, nbndsub
!         IF (etf(ibnd, 2 * ik - 1) < ef0) ivbm = ibnd
!       ENDDO
!     ENDDO    
!     !
!     CALL mp_max(ivbm, inter_pool_comm)
!     !
!     WRITE(stdout, '(5x,"Index of maximum VB ",i6)' ) ivbm
!     !
!     ALLOCATE(epsil_ipa(nqtotf), STAT = ierr)
!     IF (ierr /= 0) CALL errore('epsil_rpa', 'Error allocating epsil_ipa', 1)
    !
!     epsil_ipa = zero
    dchi  = zero
    bmatf = czero
    !
    DO iq = 1, nqtotf
      xxq_s(:) = xqf(:, iq)
      xxq(:)   = xqf(:, iq)
      CALL shift_q(xxq_s, qcart)  ! shift q point to the 1st BZ
      q2 = (qcart(1)**2 + qcart(2)**2 + qcart(3)**2) * (twopi / alat) ** 2
      IF (ABS(q2) < eps5) CYCLE
      DO ik = 1, nkf
        xxk(:) = xkf(:, 2 * ik - 1)
        CALL kpmq_map(xxk, xxq, 1, ikq)
        cufkk(:, :) = unk_all(:, :, ik + lower_bnd - 1)
        cufkq(:, :) = unk_all(:, :, ikq)
        CALL compute_umn_f(nbndsub, cufkk, cufkq, bmatf)
        !! Now include the Berry connection term
        qdotA = czero
        Amat  = czero
        IF (vme == 'wannier') THEN
          qdotA(:, :) = qcart(1) * Amnk_all(1, :, :, ik + lower_bnd - 1) + &
                        qcart(2) * Amnk_all(2, :, :, ik + lower_bnd - 1) + &
                        qcart(3) * Amnk_all(3, :, :, ik + lower_bnd - 1)
          qdotA = qdotA * (twopi / alat)
          qdotA = ci * qdotA
          ! cufkk_d = conjg(transpose(cufkk))
          ! Amat  = matmul(qdotA,  cufkk_d)
          ! qdotA = matmul(cufkq, Amat)
          CALL ZGEMM('n', 'n', nbndsub, nbndsub, nbndsub, cone, cufkq, nbndsub, qdotA, nbndsub, czero, Amat, nbndsub)
          CALL ZGEMM('n', 'c', nbndsub, nbndsub, nbndsub, cone, Amat, nbndsub, cufkk, nbndsub, czero, qdotA, nbndsub)
        ENDIF
        !
        DO ibnd = 1, nbndfst
          IF (etf_all(ibnd + ibndmin - 1, ik + lower_bnd - 1) < ef0) THEN 
            fik = 1.d0
          ELSE
            fik = 0.d0
          ENDIF
          DO jbnd = 1, nbndfst
            IF (etf_all(jbnd + ibndmin - 1, ikq) < ef0) THEN 
              fjkq = 1.d0
            ELSE 
              fjkq = 0.d0 
            ENDIF
            DO imode = 1, nmodes 
              IF (wf(imode, iq) < 1.0e-5) THEN 
                freq = 0.0
              ELSE 
                freq = wf(imode, iq) 
              ENDIF
              etrans = etf_all(jbnd + ibndmin - 1, ikq) - etf_all(ibnd + ibndmin - 1, ik + lower_bnd - 1) &
               - freq
              IF (ABS(etrans) < eps5) CYCLE
              dchi(imode, iq) = dchi(imode, iq) + fac / q2  /  etrans * &
              (felec(jbnd, ikq) - fjkq - felec(ibnd, ik + lower_bnd - 1) + fik) * &
              ABS(bmatf(jbnd + ibndmin - 1, ibnd + ibndmin - 1) +     &
              qdotA(jbnd + ibndmin - 1, ibnd + ibndmin - 1)) ** 2
            ENDDO
          ENDDO
        ENDDO
      ENDDO ! ik
    ENDDO ! iq
    !
    CALL mp_sum(dchi, inter_pool_comm)
    !
    epsil_ipa(:, :) = 1.d0 - dchi(:, :) / epsil_ipa(:, :) 
    !! Start writing the dielectric constants
    IF (my_pool_id == ionode_id) THEN
      OPEN(UNIT = 480, FILE = 'eps_ipa_eh.dat')
        DO iq = 1, nqtotf ! loop over q-points 
        xxq_s(:) = xqf(:, iq)
        CALL shift_q(xxq_s, qcart)  ! shift q point to the 1st BZ
        q2 = (qcart(1)**2 + qcart(2)**2 + qcart(3)**2) * (twopi / alat) ** 2
        WRITE(480, *) qcart(:), q2
        DO imode = 1, nmodes
          WRITE(480, *) epsil_ipa(imode, iq)
        ENDDO
        ENDDO
      CLOSE(480)
    ENDIF
    CALL mp_barrier(inter_pool_comm)
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE eps_ipa_eh
    !-----------------------------------------------------------------------  
    !      
    !-----------------------------------------------------------------------
    SUBROUTINE time_propagation()
    !-----------------------------------------------------------------------
    !! subroutine controling the time propagation of electron and phonon
    !! distribution
    !! Yiming Pan
    !-----------------------------------------------------------------------
    USE kinds,            ONLY : DP
    USE ep_constants,     ONLY : zero, one, two, twopi, hbar, ryd2ev,      &
                                 eps40, ryd2mev, czero
    USE input,            ONLY : degaussw, dt_tdbe, twrite_tdbe,        &
                                 nt_tdbe, mp_mesh_k, fsthick, nkf1d, nkf2d, &
                                 nkf3d, eps_acoustic, nkf1, nkf2, nkf3,    &
                                 nqf1, nqf2, nqf3 , solver_tdbe,  phph_tdbe, &
                                 dph_tdbe, lscreen_tdbe,     &
                                 scr_typ, assume_metal
    USE global_var,       ONLY : adapt_smearing, wf, xqf, wqf, nktotf,     &
                                 nqtotf, nbndfst, zstar
    USE tdbe_common,      ONLY : nstg, a_tab, b_tab, c_tab, nkfsden,       &
                                 indx_mapk, indx_mapq, xkf_all, iq2nqfs,   &
                                 enk_all, wkf_all, ie_ph, iph_e, iph_ph,   &
                                 ie_e, totq, selecq, ef0, tstart, fs2sec,  &
                                 felec, nphon, nphon_pre, vnk_all,         &
                                 vnuq_all, s_bztoibz_dg, bztoibz_dg,       &
                                 dt_in_ps, ikfsdx, epsil_ipa, uf_all
    USE supercond_common, ONLY : ixkff, g2, ibnd_kfs_to_kfs_all,           &
                                 ibnd_kfs_all_to_kfs
    USE io_global,        ONLY : stdout, ionode_id, meta_ionode
    USE cell_base,        ONLY : at, bg, alat, omega
    USE mp,               ONLY : mp_bcast, mp_barrier, mp_sum
    USE mp_global,        ONLY : inter_pool_comm, my_pool_id
    USE io_var,           ONLY : iufilfelec, iufilnphon
    USE modes,            ONLY : nmodes
    USE pwcom,            ONLY : ef
    USE bzgrid,           ONLY : kpmq_map
    USE phph,             ONLY : iphph_calc, iphph_rta_calc
    USE ions_base,        ONLY : amass, ityp, nat 
    !
    IMPLICIT NONE
    !
    INTEGER :: ierr
    !! Error status
    INTEGER :: ios
    !! INTEGER variable for I/O control
    INTEGER :: itt      
    !! Counter on time step 
    INTEGER :: iq, iqd
    !! Counter on coarse q-point grid
    INTEGER :: iqq
    !! Counter on coarse q-point grid
    INTEGER :: ik, ikfs, ikd, jkqd 
    !! Counter on coarse k-point grid
    INTEGER :: ikq
    !! index of  k+q on the dense k-point grid
    INTEGER :: ibnd, ibndfs
    !! Counter on band
    INTEGER :: jbnd, jbndfs
    !! Counter on band
    INTEGER :: ir
    !! Counter for WS loop
    INTEGER :: nu
    !! counter on mode
    INTEGER :: istg, pp
    !! counter
    INTEGER :: ns
    !! Index on symmetry operations
    REAL(KIND = DP):: fdos1, fdos2, dt
    !! time-dependent joint density of state
    REAL(KIND = DP) :: fik, fjkq, nqv
    !! electron and phonon population
    REAL(KIND = DP) :: fact
    !! two*pi/hbar
    REAL(KIND = DP) :: inv_degaussw
    !! 1.0/degaussw. Defined for efficiency reasons.
    !INTEGER :: nt_tdbe
    REAL(KIND = DP) :: delta1 , delta2
    !! delta function
    REAL(KIND = DP) :: etmp1
    !! Temporary variable to store etmp1 = ek - ekq + wq
    REAL(KIND = DP) :: etmp2
    !! Index of lowest conduction band for every kpts in the fermi shell
    REAL(KIND = DP) :: xxk(3), xxq(3), qcart(3), zbg(3), qabs
    !! Current k/q-point
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! Delta function
    REAL(KIND = DP) :: vdiff(3)
    !! vnuq- vjk+q
    LOGICAL  :: screen_g2
    !
    INTEGER  ::  iat, ipol
    !
    IF (adapt_smearing) THEN
      inv_degaussw = zero
    ELSE
      inv_degaussw = one / degaussw
    ENDIF
    fact = twopi / (hbar / ryd2ev)
    dt = dt_tdbe * fs2sec
    WRITE(stdout,'(5x,a,a10,a)') 'Solving TDBE with ', solver_tdbe ,' method.'
    ! WRITE(stdout,'(5x,a)') 'The Butcher Tableaux for this solver are: '
    ! WRITE(stdout,'(5x,a,i)') 'nstg: ', nstg
    ! WRITE(stdout,'(5x,a,6f5.3)') 'b_tab : ', b_tab(1:6)
    ! WRITE(stdout,'(5x,a,6f5.3)') 'c_tab : ', c_tab(1:6)
    ! WRITE(stdout,'(5x,a)') 'a_tab : '
    ! WRITE(stdout,*) '      ', a_tab(:,:)
    screen_g2 = .FALSE.
    IF (lscreen_tdbe .AND. scr_typ > 1 .AND. (.NOT. assume_metal) ) THEN 
      WRITE(stdout, '(5x, a)') 'The e-ph matrix elements are screened by light-induced carriers'
      WRITE(stdout, '(5x, a)') 'The dielectric function is calculated under IPA and carriers &
      are fixed with its initial distribution.'
      WRITE(stdout, '(5x, a)') 'Time-dependent screening is still under development'
      screen_g2 = .TRUE.
      !! inverse dielectric function
      epsil_ipa(:, :) = one / epsil_ipa(:, :)
      epsil_ipa(:, :) = epsil_ipa(:, :) ** two
      WRITE(stdout, '(5x, a , f16.10)') 'The maximum value of 1/epsilon(nu, iq)**2 is : ', MAXVAL(epsil_ipa)
      WRITE(stdout, '(5x, a , f16.10)') 'The minimum value of 1/epsilon(nu, iq)**2 is : ', MINVAL(epsil_ipa) 
    ENDIF
    DO itt = 1, nt_tdbe
      CALL start_clock('TDBE: dt')
      ie_ph(:,:,:) = 0.d0
      iph_e(:,:,:) = 0.d0
      DO istg = 1, nstg
        DO ikfs = 1, nkfsden
          ik = ikfsdx(ikfs)
          ikd = indx_mapk(ik)
          xxk = xkf_all(:, ik)
          IF (ixkff(ikd) == 0) WRITE(stdout, '(5x, a)') 'warning : ixkff(ikd) == 0'
          DO iqq = 1, totq
            iq = selecq(iqq)
            IF (iq == 0) CALL errore('tdbe','iq should not be 0',1)
            IF (mp_mesh_k) THEN
              ns = s_bztoibz_dg(ikd)
              CALL rotate_kpoint(ns, indx_mapq(iq), iqd, nkf1d, nkf2d, nkf3d)
              if (iqd == 0 .and. itt == 1) THEN
                WRITE(stdout,'(5x,a,i5,a,i5,a)') 'map of q point ', indx_mapq(iq),' with operation', ns, 'fails'
                CYCLE
              ENDIF
            ELSE
              iqd = indx_mapq(iq)
            ENDIF
            xxq = xqf(:, iq)
            IF (ixkff(ikd) /= 0 .AND. iq2nqfs(ixkff(ikd), iqd) > 0) THEN
              CALL kpmq_map(xxk, xxq, +1, ikq)
              jkqd = indx_mapk(ikq)
              IF (ixkff(jkqd) == 0) WRITE(stdout, '(5x, a)') 'warning : ixkff(jkqd) == 0'
              IF (ikq == 0) CALL errore('tdbe','ikq should not be 0',1)
              DO ibnd = 1, nbndfst
                IF (ABS(enk_all(ibnd, ik) - ef) < fsthick) THEN
                  DO jbnd = 1, nbndfst
                    IF (ABS(enk_all(jbnd, ikq) - ef) < fsthick) THEN
                      ibndfs = ibnd_kfs_all_to_kfs(ibnd, ixkff(ikd))
                      jbndfs = ibnd_kfs_all_to_kfs(jbnd, ixkff(jkqd))
                      IF (ibndfs == 0 .OR. jbndfs == 0) CYCLE
                      DO nu = 1, nmodes
                        IF (wf(nu, iq) > eps_acoustic) THEN
                          fik = felec(ibnd, ik)
                          fjkq = felec(jbnd, ikq)
                          nqv = nphon(nu, iq)
                          IF (istg > 1) THEN 
                            DO pp = 1, istg-1
                              fik = fik + a_tab(istg, pp) * ie_ph(ibnd, ik, pp) * dt
                              fjkq = fjkq + a_tab(istg, pp) * ie_ph(jbnd, ikq, pp) * dt
                              nqv = nqv + a_tab(istg, pp) * iph_e(nu, iq, pp) * dt
                            ENDDO
                          ENDIF
                          ! delta[E_k - E_k+q + w_q] and delta[E_k - E_k+q - w_q]
                          !! absorption of a phonon
                          etmp1  = enk_all(ibnd, ik) - enk_all(jbnd, ikq) + wf(nu, iq)
                          !! emission of a phonon
                          etmp2  = enk_all(ibnd, ik) - enk_all(jbnd, ikq) - wf(nu, iq)
                          IF (adapt_smearing) THEN
                            vdiff(:) = vnuq_all(:, nu, iq) - vnk_all(:, jbnd, ikq)
                            inv_degaussw = one / eta(vdiff)
                            delta1 = w0gauss(etmp1 * inv_degaussw,0) * inv_degaussw
                            delta2 = w0gauss(etmp2 * inv_degaussw,0) * inv_degaussw
                          ELSE
                            delta1 = w0gauss(etmp1 * inv_degaussw,0) * inv_degaussw
                            delta2 = w0gauss(etmp2 * inv_degaussw,0) * inv_degaussw
                          ENDIF
                          fdos1 = (one + nqv) * (one - fik) *     &
                                 & fjkq - nqv * fik *  &
                                 & (one-fjkq)                      
                          fdos2 = -fik * (one - fjkq) *        &
                                 & (one + nqv) + (one - fik) *    &
                                 & fjkq * nqv
                          IF (screen_g2) THEN
                            ie_ph(ibnd, ik, istg)=ie_ph(ibnd, ik, istg) + fact * wqf(iq) *     &
                                         & (delta1 * fdos1 + delta2 * fdos2) *       &
                                         & g2(ixkff(ikd), iq2nqfs(ixkff(ikd), iqd), ibndfs, jbndfs, nu) &
                                         & * epsil_ipa(nu, iq) 
                            iph_e(nu, iq, istg)=iph_e(nu, iq, istg) + two * fact * wkf_all(ik)*   &
                                         & g2(ixkff(ikd), iq2nqfs(ixkff(ikd), iqd), ibndfs, jbndfs, nu) * delta1 * fdos1 &
                                         & * epsil_ipa(nu, iq) 
                          ELSE
                            ie_ph(ibnd, ik, istg)=ie_ph(ibnd, ik, istg) + fact * wqf(iq) *     &
                                         & (delta1 * fdos1 + delta2 * fdos2) *       &
                                         & g2(ixkff(ikd), iq2nqfs(ixkff(ikd), iqd), ibndfs, jbndfs, nu)
                            iph_e(nu, iq, istg)=iph_e(nu, iq, istg) + two * fact * wkf_all(ik)*   &
                                         & g2(ixkff(ikd), iq2nqfs(ixkff(ikd), iqd), ibndfs, jbndfs, nu) * delta1 * fdos1
                          ENDIF ! screen_g2                                                          
                         ENDIF
                      ENDDO ! nu
                    ENDIF ! ekq
                  ENDDO ! jbnd
                ENDIF ! ekk
              ENDDO ! ibnd
            ENDIF ! iq2nqfs
          !iq2nqfs(nkfs,nqtotfd)
          ENDDO ! iqq
        ENDDO ! ik    
        CALL mp_barrier(inter_pool_comm)
        CALL mp_sum(ie_ph(:, :, istg), inter_pool_comm)
        CALL mp_sum(iph_e(:, :, istg), inter_pool_comm)
        CALL mp_barrier(inter_pool_comm)
      ENDDO !  istg
      !
      !! If include ph-ph interaction, calculate the corresponding collsion integral
      IF (phph_tdbe .AND. MODULO(itt, dph_tdbe) == 0) THEN
        IF (.FALSE.) THEN
          CALL iphph_rta_calc()
        ELSE
          CALL iphph_calc()
        ENDIF
      ENDIF
      !! Now update the distribution 
      !! e-ph interaction
      DO istg = 1, nstg
        felec(:, :) = felec(:, :) + dt * b_tab(istg) * ie_ph(:, :, istg)
        !!! electron-electron interaction within RTA
        nphon(:, :) = nphon(:, :) + dt * b_tab(istg) * iph_e(:, :, istg)
      ENDDO
      !
      !! ph-ph interaction
      IF (phph_tdbe .AND. MODULO(itt, dph_tdbe) == 0) THEN
        DO istg = 1, nstg
          nphon(:, :) = nphon(:, :) + dt_in_ps * b_tab(istg) * iph_ph(:, :, istg)
        ENDDO
        nphon_pre(:, :) = nphon(:, :)
      ENDIF      

      !! output 
      IF (itt == 1) THEN
        IF (my_pool_id == ionode_id) THEN
          OPEN (unit = 108, file = 'ie_ph.dat')
          DO ik = 1, nktotf
            DO ibnd = 1, nbndfst
              WRITE(108,*) enk_all(ibnd,ik), ie_ph(ibnd,ik,1)
            ENDDO
          ENDDO
          CLOSE(108)
          OPEN (unit = 109, file = 'iph_e.dat')
          DO iq = 1, nqtotf 
            DO nu = 1, nmodes
              WRITE(109,*) wf(nu,iq), iph_e(nu,iq,1)
            ENDDO
          ENDDO
          CLOSE(109)
        ENDIF
      ENDIF
      IF (MODULO(itt, twrite_tdbe) == 0) THEN
        WRITE(stdout, '(/5x, a, f15.6 , a)' ) 'Electron and phonon distributions at ',itt*dt_tdbe+tstart, ' fs are written to files'
        CALL write_dist(.false.,itt*dt_tdbe+tstart)
      ENDIF   
      CALL stop_clock('TDBE: dt')
    ENDDO

    CONTAINS
      !
      FUNCTION eta(v)
      IMPLICIT NONE
      !
      REAL(KIND = DP), INTENT(in) :: v(3)
      REAL(KIND = DP) :: eta_tmp(3)
      REAL(KIND = DP) :: eta
      !
      IF (SQRT(DOT_PRODUCT(v, v)) < eps40) THEN
        eta= 1.0d0 / ryd2mev
      ELSE
        eta_tmp(1) = (twopi / alat) * ABS(DOT_PRODUCT(v(:), bg(:, 1)) / DBLE(nqf1))
        eta_tmp(2) = (twopi / alat) * ABS(DOT_PRODUCT(v(:), bg(:, 2)) / DBLE(nqf2))
        eta_tmp(3) = (twopi / alat) * ABS(DOT_PRODUCT(v(:), bg(:, 3)) / DBLE(nqf3))
        eta = 0.5d0 * DSQRT(eta_tmp(1)**two + eta_tmp(2)**two + &
                             eta_tmp(3) ** two) / SQRT(12.0d0)
      ENDIF
      ! If the smearing is too small, set 1 meV. Too small value are numerically unstable.
      IF (eta * ryd2mev < 1.0d0) eta = 1.0d0 / ryd2mev
      !
      RETURN
        !
      END FUNCTION eta
    !-----------------------------------------------------------------------
    END SUBROUTINE time_propagation
    !-----------------------------------------------------------------------
  !
  !-------------------------------------------------------------------------
  END MODULE tdbe_mod
  !-------------------------------------------------------------------------
