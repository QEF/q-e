  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
  MODULE kfold
  !-----------------------------------------------------------------------
  !!
  !! This module contains routine to fold k+q grids back on the k-grid
  !!
  !  
  USE kinds, ONLY : DP 
  ! 
  IMPLICIT NONE
  !
  INTEGER :: g0vec_all(3, 125)
  !! G-vectors needed to fold the k+q grid into the k grid 
  INTEGER :: ng0vec           
  !! number of inequivalent such translations (125)
  INTEGER, ALLOCATABLE :: shift(:)
  !! for every k+q, index of the G_0-vector needed to fold k+q into k+q+G0 
  INTEGER, ALLOCATABLE :: gmap(:, :)        
  !! the map G --> G-G_0 in the large (density) set, for every G_0 (125 at most)
  REAL(KIND = DP) :: g0vec_all_r(3, 125)
  !! G-vectors needed to fold the k+q grid into the k grid, cartesian coord. 
  ! 
  CONTAINS
    ! 
    !-----------------------------------------------------------------------
    SUBROUTINE createkmap(xq)
    !-----------------------------------------------------------------------
    !!  
    !! This SUBROUTINE is called from elphon_shuffle_wrap for each
    !! nqc1*nqc2*nqc3 phonon on the coarse mesh.    
    !!
    !! It folds the k+q mesh into the k mesh using 5^3 G_0 translations 
    !!
    !! SP - 2016 - iverbosity cannot be tested here. Generates Tb of data ... 
    !
    USE kinds,         ONLY : DP
    USE cell_base,     ONLY : at, bg
    USE klist,         ONLY : nkstot, xk
    USE epwcom,        ONLY : nkc1, nkc2, nkc3
    USE io_files,      ONLY : prefix
    USE io_var,        ONLY : iukmap
    USE klist_epw,     ONLY : kmap
    USE io_global,     ONLY : meta_ionode
    USE mp,            ONLY : mp_barrier
    USE elph2,         ONLY : xkq
    USE constants_epw, ONLY : eps5, zero
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: xq(3)
    !! Coords. of q-point 
    !
    ! Local variables
    LOGICAL :: in_the_list
    !! Is the point in the list
    LOGICAL :: found
    !! Is the point found
    INTEGER :: ik
    !! k-point index
    INTEGER :: jk
    !! k-point index
    INTEGER :: i
    !! Index of the k+q on the k-grid
    INTEGER :: j
    !! Index of the k+q on the k-grid
    INTEGER :: k
    !! Index of the k+q on the k-grid
    INTEGER :: ig1, ig2, ig3
    !! Counter on G_0-translations
    INTEGER :: g0vec(3)
    !! G_0 in INTEGER coords.
    INTEGER :: ig0
    !! Index of G_0 such that k+q+G_0 belongs to the 1st BZ
    INTEGER :: n
    !! Mapping index of k+q on k
    REAL(KIND = DP) :: xk_q(3)
    !! Coords. of k+q-point
    REAL(KIND = DP) :: xx_c(nkstot), yy_c(nkstot), zz_c(nkstot)
    !! k-points in crystal coords. in multiple of nkc1, nkc2, nkc3
    REAL(KIND = DP) :: xx, yy, zz
    !! k+q in crystal coords. in multiple of nkc1, nkc2, nkc3
    REAL(KIND = DP) :: xx_n, yy_n, zz_n
    !! k+q in crystal coords. in multiple of nkc1, nkc2, nkc3 in 1st BZ
    !
    IF (meta_ionode) THEN
      !
      !  Now fold k+q back into the k-grid for wannier interpolation.
      !  Since this is done before divide and impera, every pool has all the kpoints.
      !
      !  bring q in crystal coordinates and check commensuration
      !  loosy tolerance: not important since k+q is defined through NINT() 
      !  bring q-point from cartesian to crystal coords.  
      !
      CALL cryst_to_cart(1, xq, at, -1)
      !
      xx = xq(1) * nkc1 
      yy = xq(2) * nkc2 
      zz = xq(3) * nkc3 
      in_the_list = ABS(xx - NINT(xx)) <= eps5 .AND. &
                    ABS(yy - NINT(yy)) <= eps5 .AND. &
                    ABS(zz - NINT(zz)) <= eps5
      IF (.NOT. in_the_list) CALL errore('createkmap', 'q-vec not commensurate', 1)
      !
      ng0vec = 0
      DO ig1 = -2, 2
        DO ig2 = -2, 2
          DO ig3 = -2, 2
            ng0vec = ng0vec + 1
            g0vec_all(1, ng0vec) = ig1
            g0vec_all(2, ng0vec) = ig2
            g0vec_all(3, ng0vec) = ig3
          ENDDO
        ENDDO
      ENDDO
      !
      !  bring all the k points from cartesian to crystal coordinates
      !
      CALL cryst_to_cart(nkstot, xk, at, -1)
      !
      DO ik = 1, nkstot
        !
        !  check that the k's are actually on a uniform mesh centered at gamma
        !
        xx_c(ik) = xk(1, ik) * nkc1
        yy_c(ik) = xk(2, ik) * nkc2
        zz_c(ik) = xk(3, ik) * nkc3
        in_the_list = ABS(xx_c(ik) - NINT(xx_c(ik))) <= eps5 .AND. &
                      ABS(yy_c(ik) - NINT(yy_c(ik))) <= eps5 .AND. &
                      ABS(zz_c(ik) - NINT(zz_c(ik))) <= eps5
        IF (.NOT. in_the_list) CALL errore('createkmap', 'is this a uniform k-mesh?', 1)
        !
        IF ((xx_c(ik) < -eps5) .OR. (yy_c(ik) < -eps5) .OR. (zz_c(ik) < -eps5)) THEN
          CALL errore('createkmap', 'coarse k-mesh needs to be strictly positive in 1st BZ', 1)
        ENDIF 
      ENDDO
      !
      DO ik = 1, nkstot
        !
        ! Now add the phonon wavevector and check that k+q falls again on the k grid
        !
        xk_q(:) = xk(:, ik) + xq(:)
        !
        xx = xk_q(1) * nkc1
        yy = xk_q(2) * nkc2
        zz = xk_q(3) * nkc3
        in_the_list = ABS(xx - NINT(xx)) <= eps5 .AND. &
                      ABS(yy - NINT(yy)) <= eps5 .AND. &
                      ABS(zz - NINT(zz)) <= eps5
        IF (.NOT. in_the_list) CALL errore('createkmap', 'k+q does not fall on k-grid', 1)
        !
        ! Find the index of this k+q in the k-grid
        !
        i = MOD(NINT(xx + 2 * nkc1), nkc1) 
        j = MOD(NINT(yy + 2 * nkc2), nkc2) 
        k = MOD(NINT(zz + 2 * nkc3), nkc3) 
        !
        xx_n = xx
        yy_n = yy
        zz_n = zz
        !
        !  make sure xx_n, yy_n and zz_n are in 1st BZ
        !
        CALL backtoBZ(xx_n, yy_n, zz_n, nkc1, nkc2, nkc3)
        !
        n = 0
        found = .FALSE.
        DO jk = 1, nkstot
          found = NINT(xx_c(jk)) == NINT(xx_n) .AND. &
                  NINT(yy_c(jk)) == NINT(yy_n) .AND. &
                  NINT(zz_c(jk)) == NINT(zz_n)
          IF (found) THEN
            n = jk
            EXIT
          ENDIF
        ENDDO
        !
        !  26/06/2012 RM
        !  since coarse k- and q- meshes are commensurate, one can easily find n
        !  n = NINT(xx_n) * nk2 * nk3 + NINT(yy_n) * nk3 + NINT(zz_n) + 1
        !  n represents the index of k+q on the coarse k-grid.
        !
        IF (n == 0) CALL errore('createkmap','problem indexing k+q', 1)
        !
        kmap(ik) = n
        !
        ! Determine the G_0 such that k+q+G_0 belongs to the first BZ
        !
        g0vec(1) = (i - NINT(xx)) / nkc1
        g0vec(2) = (j - NINT(yy)) / nkc2
        g0vec(3) = (k - NINT(zz)) / nkc3
        !
        !  now store the shift for this k+q point
        !
        in_the_list = .FALSE.
        ig0 = 0
        DO WHILE((ig0 <= ng0vec) .AND. (.NOT. in_the_list))
          ig0 = ig0 + 1
          in_the_list = ((ABS(g0vec(1) - g0vec_all(1,ig0)) <= eps5) .AND. &
                         (ABS(g0vec(2) - g0vec_all(2,ig0)) <= eps5) .AND. &
                         (ABS(g0vec(3) - g0vec_all(3,ig0)) <= eps5))
        ENDDO
        shift(ik) = ig0
        !
        IF (.NOT. in_the_list) CALL errore('createkmap', 'cannot find the folding vector in the list', 1)
        !
        ! very important: now redefine k+q through the corresponding kpoint on the k mesh
        ! Note that this will require using the periodic gauge in the calculation of the
        ! electron-phonon matrix elements (factor e^iG_0*r if G_0 is the vector used for refolding)
        xkq(:, ik) = xk(:, n) 
        !
      ENDDO
      !
      ! bring k-points, q-point, and G_0-vectors back to cartesian coordinates 
      !
      CALL cryst_to_cart(1, xq, bg, 1) 
      CALL cryst_to_cart(nkstot, xk, bg, 1)
      !
      g0vec_all_r = DBLE(g0vec_all)
      CALL cryst_to_cart(ng0vec, g0vec_all_r, bg, 1)
      !
      !  the unit with kmap(ik) and shift(ik)
      ! 
      OPEN(iukmap, FILE = TRIM(prefix) // '.kmap', FORM = 'formatted')
      DO ik = 1, nkstot
        WRITE(iukmap,'(3i6)') ik, kmap(ik), shift(ik) 
      ENDDO
      CLOSE(iukmap)
      !
    ENDIF
    !CALL mp_barrier(world_comm)
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE createkmap
    !-----------------------------------------------------------------------
    ! 
    !-----------------------------------------------------------------------
    SUBROUTINE createkmap2(xxq)
    !-----------------------------------------------------------------------
    !!
    !! generate the map k+q --> k for folding the rotation matrix U(k+q) 
    !! 
    !! in parallel case, this SUBROUTINE must be called only by first proc 
    !! (which has all the kpoints)
    !!
    !
    USE kinds,         ONLY : DP
    USE cell_base,     ONLY : at, bg
    USE klist,         ONLY : nkstot, xk
    USE klist_epw,     ONLY : kmap  
    USE epwcom,        ONLY : nkc1, nkc2, nkc3
    USE elph2,         ONLY : xkq
    USE constants_epw, ONLY : eps5, zero
    ! 
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: xxq(3)
    !! The current q-point 
    ! 
    ! Local variables
    LOGICAL :: in_the_list
    !! Is the file in the list
    LOGICAL :: found
    !! Is the file found
    INTEGER :: ik
    !! K-point index
    INTEGER :: jk
    !! Another k-point index
    INTEGER :: n
    !! Mapping index of k+q on k
    REAL(KIND = DP) :: xx_c(nkstot), yy_c(nkstot), zz_c(nkstot)
    !! k-points in crystal coords. in multiple of nkc1, nkc2, nkc3
    REAL(KIND = DP) :: xx, yy, zz
    !! k+q in crystal coords. in multiple of nkc1, nkc2, nkc3
    !
    ! The first proc keeps a copy of all kpoints !
    ! bring q from cartesian to crystal coordinates and check commensuration
    ! 
    CALL cryst_to_cart(1, xxq, at, -1)
    !
    xx = xxq(1) * nkc1 
    yy = xxq(2) * nkc2 
    zz = xxq(3) * nkc3 
    in_the_list = ABS(xx - NINT(xx)) <= eps5 .AND. &
                  ABS(yy - NINT(yy)) <= eps5 .AND. &
                  ABS(zz - NINT(zz)) <= eps5
    IF (.NOT. in_the_list) CALL errore('createkmap2', 'q-vec not commensurate', 1)
    !
    !  bring all the k-points from cartesian to crystal coordinates 
    !
    CALL cryst_to_cart(nkstot, xk, at, -1)
    !
    DO ik = 1, nkstot
      !
      !  check that the k's are actually on a uniform mesh centered at gamma
      !
      xx_c(ik) = xk(1, ik) * nkc1
      yy_c(ik) = xk(2, ik) * nkc2
      zz_c(ik) = xk(3, ik) * nkc3
      in_the_list = ABS(xx_c(ik) - NINT(xx_c(ik))) <= eps5 .AND. &
                    ABS(yy_c(ik) - NINT(yy_c(ik))) <= eps5 .AND. &
                    ABS(zz_c(ik) - NINT(zz_c(ik))) <= eps5
      IF (.NOT. in_the_list) CALL errore('createkmap2', 'is this a uniform k-mesh?', 1)
      !
      IF ((xx_c(ik) < -eps5) .OR. (yy_c(ik) < -eps5) .OR. (zz_c(ik) < -eps5)) THEN
        CALL errore('createkmap2', 'coarse k-mesh needs to be strictly positive in 1st BZ', 1)
      ENDIF 
    ENDDO
    !
    DO ik = 1, nkstot
      !
      ! now add the phonon wavevector and check that k+q falls again on the k grid
      !
      xkq(:, ik) = xk(:, ik) + xxq(:)
      !
      xx = xkq(1, ik) * nkc1
      yy = xkq(2, ik) * nkc2
      zz = xkq(3, ik) * nkc3
      in_the_list = ABS(xx - NINT(xx)) <= eps5 .AND. &
                    ABS(yy - NINT(yy)) <= eps5 .AND. &
                    ABS(zz - NINT(zz)) <= eps5
      IF (.NOT. in_the_list) CALL errore('createkmap2', 'k+q does not fall on k-grid', 1)
      !
      ! Find the index of this k+q in the k-grid
      !
      ! make sure xx, yy and zz are in 1st BZ
      !
      CALL backtoBZ(xx, yy, zz, nkc1, nkc2, nkc3)
      !
      n = 0
      found = .FALSE.
      DO jk = 1, nkstot
        !
        found = NINT(xx_c(jk)) == NINT(xx) .AND. &
                NINT(yy_c(jk)) == NINT(yy) .AND. &
                NINT(zz_c(jk)) == NINT(zz)
        IF (found) THEN
          n = jk
          EXIT
        ENDIF
      ENDDO
      !
      IF (n == 0) CALL errore('createkmap2', 'problem indexing k+q', 1)
      !
      kmap(ik) = n
      !
    ENDDO
    !
    ! bring everybody back to cartesian coordinates 
    !
    CALL cryst_to_cart(1, xxq, bg, 1) 
    CALL cryst_to_cart(nkstot, xk, bg, 1)
    !
    RETURN
    !-------------------------------------------------------------------------
    END SUBROUTINE createkmap2
    !-------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------
    SUBROUTINE createkmap_pw2()
    !-------------------------------------------------------------------------
    !!
    !! Creates the first instance of [prefix].kgmap. 
    !!
    USE kinds,         ONLY : DP
    USE cell_base,     ONLY : at, bg
    USE epwcom,        ONLY : nkc1, nkc2, nkc3
    USE pwcom,         ONLY : nkstot
    USE klist_epw,     ONLY : xk_cryst
    USE io_global,     ONLY : stdout, meta_ionode
    USE io_files,      ONLY : prefix
    USE io_var,        ONLY : iukgmap
    USE gvect,         ONLY : ngm, ngm_g, gcutm
    USE fft_base,      ONLY : dfftp
    USE fft_types,     ONLY : fft_stick_index
    USE fft_ggen,      ONLY : fft_set_nl
    USE constants,     ONLY : eps8
    USE constants_epw, ONLY : eps5
#if defined(__NAG) 
    USE f90_unix_io,   ONLY : flush
#endif
    USE mp_global,     ONLY : inter_pool_comm, inter_image_comm
    USE mp,            ONLY : mp_barrier
    !
    IMPLICIT NONE
    !
    ! Local variables
    LOGICAL :: is_local
    !! Local variable
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: ig1, ig2, ig3
    !! Counter on G_0-translations
    INTEGER :: ig0
    !! Index of G_0-vector at the origin
    INTEGER :: ngm_save, ngm_max, ngm_local
    !! Local variables for the nr. of G-vectors
    INTEGER :: ni, nj, nk
    !! Max Miller indices
    INTEGER ::  i, j, k
    !! Counter on Miller indices
    INTEGER :: istart, jstart, kstart
    !! Starting counter on Miller indices
    INTEGER :: ng
    !! Counter on G-vectors
    INTEGER :: ierr
    !! Error status
    INTEGER, ALLOCATABLE :: ig_l2g(:)
    !! Converts a local G-vector index into the global index 
    INTEGER, ALLOCATABLE :: g2l(:)
    !! Local index of G-vector 
    INTEGER, ALLOCATABLE :: mill_unsorted(:, :)
    !! Array of unsorted Miller indices of G-vectors
    INTEGER, ALLOCATABLE :: mill(:, :)
    !! Array of sorted Miller indices of G-vectors in increasing order of G^2 
    INTEGER, ALLOCATABLE :: igsrt(:)
    !! Array of G-vector indices in the initial (unsorted) list 
    !! (index of i-th G-vector in the unsorted list of G-vectors) 
    INTEGER, ALLOCATABLE :: jtoi(:)
    !! For the i-th G-vector in the sorted list, jtoi(i) 
    !! returns its index in the unsorted list
    INTEGER, ALLOCATABLE :: itoj(:)
    !! itoj(i) returns the index of the G-vector in the sorted list 
    !! that was at i-th position in the unsorted list
    REAL(KIND = DP) :: xx, yy, zz
    !! k-point in crystal coords. in multiple of nkc1, nkc2, nkc3
    REAL(KIND = DP) :: tx(3), ty(3), t(3)
    !! Reciprocal lattice vectors
    REAL(KIND = DP), ALLOCATABLE :: gg(:)
    !! G^2 in increasing order
    REAL(KIND = DP), ALLOCATABLE :: g(:, :)
    !! G-vectors cartesian components in increasing order of G^2
    REAL(KIND = DP), ALLOCATABLE :: g2sort_g(:)
    !! G-vectors for the current processor
    REAL(KIND = DP), ALLOCATABLE :: tt(:)
    !! Temporal array
    ! 
    IF (meta_ionode) THEN
      !
      WRITE(stdout, '(/5x,a)') 'Calculating kgmap'
      FLUSH(stdout)
      !
      OPEN(iukgmap, FILE = TRIM(prefix) // '.kgmap', FORM = 'formatted')
      ! 
      ! the 5^3 possible G_0 translations
      ng0vec = 0
      DO ig1 = -2, 2
        DO ig2 = -2, 2
          DO ig3 = -2, 2
            ng0vec = ng0vec + 1
            g0vec_all(1, ng0vec) = ig1
            g0vec_all(2, ng0vec) = ig2
            g0vec_all(3, ng0vec) = ig3
          ENDDO
        ENDDO
      ENDDO
      ig0 = NINT(DBLE(ng0vec) / 2)
      !
      ALLOCATE(shift(nkstot), STAT = ierr)
      IF (ierr /= 0) CALL errore('createkmap_pw2', 'Error allocating shift', 1)
      ! 
      DO ik = 1, nkstot       
        !
        xx = xk_cryst(1, ik) * nkc1
        yy = xk_cryst(2, ik) * nkc2
        zz = xk_cryst(3, ik) * nkc3
        ! check that the k-mesh was defined in the positive region of 1st BZ
        !
        IF ((xx < -eps5) .OR. (yy < -eps5) .OR. (zz < -eps5)) THEN 
          CALL errore('createkmap_pw2','coarse k-mesh needs to be strictly positive in 1st BZ', 1)
        ENDIF 
        ! 
        shift(ik) = ig0
        WRITE(iukgmap,'(3i6)') ik, shift(ik)
        !
      ENDDO
      DEALLOCATE(shift, STAT = ierr)
      IF (ierr /= 0) CALL errore('createkmap_pw2', 'Error deallocating shift', 1)
      !
      g0vec_all_r = DBLE(g0vec_all)
      ! bring G_0 vectors from crystal to cartesian coordinates
      CALL cryst_to_cart(ng0vec, g0vec_all_r, bg, 1)
      !
      WRITE(iukgmap,'(i5)') ng0vec
      DO ig0 = 1, ng0vec
        WRITE(iukgmap, '(3f20.15)') g0vec_all_r(:, ig0)
      ENDDO
      !
    ENDIF
    !   
    CALL mp_barrier(inter_pool_comm)
    CALL mp_barrier(inter_image_comm)
    !
    ! RM: The following is adapted from ggen SUBROUTINE in Modules/recvec_subs.f90
    !
    ngm_max = ngm_g
    !
    ! save current value of ngm
    !
    ngm_save = ngm
    !
    ngm = 0
    ngm_local = 0
    !
    ! and computes all the g vectors inside a sphere
    !
    ALLOCATE(mill_unsorted(3, ngm_save), STAT = ierr)
    IF (ierr /= 0) CALL errore('createkmap_pw2', 'Error allocating mill_unsorted', 1)
    ALLOCATE(igsrt(ngm_max), STAT = ierr)
    IF (ierr /= 0) CALL errore('createkmap_pw2', 'Error allocating igsrt', 1)
    ALLOCATE(g2l(ngm_max), STAT = ierr)
    IF (ierr /= 0) CALL errore('createkmap_pw2', 'Error allocating g2l', 1)
    ALLOCATE(g2sort_g(ngm_max), STAT = ierr)
    IF (ierr /= 0) CALL errore('createkmap_pw2', 'Error allocating g2sort_g', 1)
    ALLOCATE(ig_l2g(ngm_max), STAT = ierr)
    IF (ierr /= 0) CALL errore('createkmap_pw2', 'Error allocating ig_l2g', 1)
    ALLOCATE(mill(3, ngm_max), STAT = ierr)
    IF (ierr /= 0) CALL errore('createkmap_pw2', 'Error allocating mill', 1)
    ALLOCATE(jtoi(ngm_max), STAT = ierr)
    IF (ierr /= 0) CALL errore('createkmap_pw2', 'Error allocating jtoi', 1)
    ALLOCATE(itoj(ngm_max), STAT = ierr)
    IF (ierr /= 0) CALL errore('createkmap_pw2', 'Error allocating itoj', 1)
    ALLOCATE(g(3, ngm_max), STAT = ierr)
    IF (ierr /= 0) CALL errore('createkmap_pw2', 'Error allocating g', 1)
    ALLOCATE(gg(ngm_max), STAT = ierr)
    IF (ierr /= 0) CALL errore('createkmap_pw2', 'Error allocating gg', 1)
    !
    ! Set the total number of FFT mesh points and and initial value of gg.
    ! The choice of gcutm is due to the fact that we have to order the
    ! vectors after computing them
    !
    gg(:) = gcutm + 1.d0
    !
    g2sort_g(:) = 1.0d20
    !
    ! Allocate temporal array
    !
    ALLOCATE(tt(dfftp%nr3), STAT = ierr)
    IF (ierr /= 0) CALL errore('createkmap_pw2', 'Error allocating tt', 1)
    !
    ! max miller indices (same convention as in module stick_set)
    !
    ni = (dfftp%nr1 - 1) / 2
    nj = (dfftp%nr2 - 1) / 2
    nk = (dfftp%nr3 - 1) / 2
    !
    istart = -ni
    DO i = istart, ni
      !
      tx(1:3) = i * bg(1:3, 1)
      !
      jstart = -nj
      DO j = jstart, nj
        !
        IF (dfftp%lpara .AND. fft_stick_index(dfftp, i, j) == 0) THEN
          is_local = .FALSE.
        ELSE
          is_local = .TRUE.
        ENDIF
        !
        ty(1:3) = tx(1:3) + j * bg(1:3, 2)
        !
        ! Compute all the norm square
        !
        kstart = -nk
        DO k = kstart, nk
          !
          t(1) = ty(1) + k * bg(1, 3)
          t(2) = ty(2) + k * bg(2, 3)
          t(3) = ty(3) + k * bg(3, 3)
          tt(k - kstart + 1) = t(1)**2 + t(2)**2 + t(3)**2
        ENDDO
        !
        ! Save all the norm square within cutoff
        !
        DO k = kstart, nk
          IF (tt(k - kstart + 1) <= gcutm) THEN
            ngm = ngm + 1
            IF (ngm > ngm_max) CALL errore('createkmap_pw2 1', 'too many g-vectors', ngm)
            IF (tt(k - kstart + 1) > eps8) THEN
              g2sort_g(ngm) = tt(k - kstart + 1)
            ELSE
              g2sort_g(ngm) = 0.d0
            ENDIF
            IF (is_local) THEN
              ngm_local = ngm_local + 1
              mill_unsorted(:, ngm_local) = (/ i, j, k /)
              g2l(ngm) = ngm_local
            ELSE
              g2l(ngm) = 0
            ENDIF
          ENDIF
        ENDDO
      ENDDO !jloop
    ENDDO !iloop
    IF (ngm /= ngm_max)  CALL errore('createkmap_pw2', 'g-vectors missing !', ABS(ngm - ngm_max))
    !
    igsrt(1) = 0
    CALL hpsort_eps(ngm_g, g2sort_g, igsrt, eps8)
    DEALLOCATE(g2sort_g, STAT = ierr)
    IF (ierr /= 0) CALL errore('createkmap_pw2', 'Error deallocating g2sort_g', 1)
    DEALLOCATE(tt, STAT = ierr)
    IF (ierr /= 0) CALL errore('createkmap_pw2', 'Error deallocating tt', 1)
    !  
    ngm = 0
    !
    DO ng = 1, ngm_max
      !
      IF (g2l(igsrt(ng)) > 0) THEN
        ! fetch the indices
        i = mill_unsorted(1, g2l(igsrt(ng)))
        j = mill_unsorted(2, g2l(igsrt(ng)))
        k = mill_unsorted(3, g2l(igsrt(ng)))
        !
        ngm = ngm + 1
        !
        ig_l2g(ngm) = ng
        ! 
        g(1:3, ngm) = i * bg(:, 1) + j * bg(:, 2) + k * bg(:, 3)
        gg(ngm) = SUM(g(1:3, ngm)**2)
      ENDIF
    ENDDO !ngloop
    DEALLOCATE(g2l, STAT = ierr)
    IF (ierr /= 0) CALL errore('createkmap_pw2', 'Error deallocating g2l', 1)
    !
    IF (ngm /= ngm_save) CALL errore('createkmap_pw2', 'g-vectors (ngm) missing !', ABS(ngm - ngm_save))
    !
    CALL fft_set_nl(dfftp, at, g, mill)
    !
    DO i = 1, ngm_g
      jtoi(i) = igsrt(i)
    ENDDO !
    !
    DO i = 1, ngm_g
      itoj(jtoi(i)) = i
    ENDDO
    !
    CALL refold(ngm_g, mill, itoj, jtoi)
    !
    !CALL mp_barrier(inter_pool_comm)
    !CALL mp_barrier(inter_image_comm)
    !
    DEALLOCATE(ig_l2g, STAT = ierr)
    IF (ierr /= 0) CALL errore('createkmap_pw2', 'Error deallocating ig_l2g', 1)
    DEALLOCATE(mill, STAT = ierr)
    IF (ierr /= 0) CALL errore('createkmap_pw2', 'Error deallocating mill', 1)
    DEALLOCATE(mill_unsorted, STAT = ierr)
    IF (ierr /= 0) CALL errore('createkmap_pw2', 'Error deallocating mill_unsorted', 1)
    DEALLOCATE(igsrt, STAT = ierr)
    IF (ierr /= 0) CALL errore('createkmap_pw2', 'Error deallocating igsrt', 1)
    DEALLOCATE(jtoi, STAT = ierr)
    IF (ierr /= 0) CALL errore('createkmap_pw2', 'Error deallocating jtoi', 1)
    DEALLOCATE(itoj, STAT = ierr)
    IF (ierr /= 0) CALL errore('createkmap_pw2', 'Error deallocating itoj', 1)
    DEALLOCATE(g, STAT = ierr)
    IF (ierr /= 0) CALL errore('createkmap_pw2', 'Error deallocating g', 1)
    DEALLOCATE(gg, STAT = ierr)
    IF (ierr /= 0) CALL errore('createkmap_pw2', 'Error deallocating gg', 1)
    !
    RETURN
    !
    !-------------------------------------------------------------------------
    END SUBROUTINE createkmap_pw2
    !-------------------------------------------------------------------------
    ! 
    !-----------------------------------------------------------------------
    SUBROUTINE refold(ngm_g, mill_g, itoj, jtoi)
    !----------------------------------------------------------------------
    !!
    !! Map the indices of G+G_0 into those of G 
    !! this is used to calculate electron-phonon matrix elements by
    !! refolding the k+q points into the first BZ (original k grid)
    !!
    !! No parallelization on G-vecs at the moment  
    !! (actually this is done on the global array, but in elphel2.f90
    !! every processor has just a chunk of the array, I may need some
    !! communication)
    !!
    !! I use the rule : if not found then gmap = 0 
    !! Note that the map will be used only up to npwx (small sphere), 
    !! while the G-vectors lost in the process are on the surface of 
    !! the large sphere (density set).
    !!
    !-----------------------------------------------------------------
    USE io_global, ONLY : stdout, meta_ionode
    USE io_var,    ONLY : iukgmap
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: ngm_g
    !! Counter on G-vectors
    INTEGER, INTENT(in) :: mill_g(3, ngm_g)
    !!  Array of Miller indices of G-vectors in increasing order of G^2
    INTEGER, INTENT(in) :: jtoi(ngm_g)
    !! For the i-th G-vector in the sorted list, jtoi(i)
    !! returns its index in the unsorted list
    INTEGER, INTENT(in) :: itoj(ngm_g)
    !! itoj(i) returns the index of the G-vector in the sorted list
    !! that was at i-th position in the unsorted list
    ! 
    ! Local variables
    LOGICAL :: tfound
    !! Found
    INTEGER :: ig0 
    !! Counter on G_0 vectors
    INTEGER :: ig1, ig2
    !! Counter on G vectors
    INTEGER :: i, j, k
    !! Miller indices for G+G_0 vector
    INTEGER :: ig1_use
    !! Temporary G-vectors indices
    INTEGER :: ig2_use
    !!  
    INTEGER :: ig2_guess
    !!  
    INTEGER :: notfound
    !!  
    INTEGER :: guess_skip
    !!  
    INTEGER :: indnew
    !!
    INTEGER :: indold
    !! Counter on G_0 vectors indices for writing to file
    INTEGER :: ierr
    !! Error status
    !
    ALLOCATE(gmap(ngm_g, ng0vec), STAT = ierr)
    IF (ierr /= 0) CALL errore('refold', 'Error allocating gmap', 1)
    gmap(:, :) = 0
    guess_skip = 0
    !
    !  Loop on the inequivalent G_0 vectors
    !
    DO ig0 = 1, ng0vec
      !
      IF (ig0 == 1) THEN
        WRITE(stdout, '(/5x,"Progress kgmap: ")', ADVANCE = 'no')
        indold = 0
      ENDIF
      indnew = NINT(DBLE(ig0) / DBLE(ng0vec) * 40)
      IF (indnew /= indold) WRITE(stdout, '(a)', ADVANCE = 'no') '#'
      indold = indnew
      !
      notfound = 0
      DO ig1 = 1, ngm_g
        !
        ig1_use = itoj(ig1)
        !
        ! the initial G vector
        i = mill_g(1, ig1_use)
        j = mill_g(2, ig1_use)
        k = mill_g(3, ig1_use)
        !
        ! the final G+G_0 vector
        i = i + g0vec_all(1, ig0)
        j = j + g0vec_all(2, ig0)
        k = k + g0vec_all(3, ig0)
        !
        ig2 = 0
        tfound = .FALSE.
        !
        ! try to guess next index
        !
        ig2_guess = jtoi(ig1_use) + guess_skip
        !
        IF ((ig2_guess > 0) .AND. (ig2_guess < ngm_g + 1)) THEN
          !
          ig2_guess = itoj(ig2_guess)
          !
          IF ((i == mill_g(1, ig2_guess)) .AND. (j == mill_g(2, ig2_guess)) .AND. (k == mill_g(3, ig2_guess))) THEN
            !
            ig2_use = ig2_guess
            tfound = .TRUE.
            !
          ENDIF
          !
        ENDIF
        !
        DO WHILE ((.NOT.  tfound) .AND. (ig2 < ngm_g))
          !
          ig2 = ig2 + 1
          ig2_use = itoj(ig2)
          tfound = (i == mill_g(1, ig2_use)) .AND. (j == mill_g(2, ig2_use)) .AND. (k == mill_g(3, ig2_use))
          !
        ENDDO
        !
        IF (tfound) THEN
          gmap(ig1_use, ig0) = ig2_use
          guess_skip = jtoi(ig2_use) - jtoi(ig1_use)
        ELSE
          gmap(ig1_use, ig0) = 0
          notfound = notfound + 1
        ENDIF
      ENDDO
    ENDDO ! ng0vec
    ! 
    !  output on file for electron-phonon matrix elements
    !
    IF (.NOT. meta_ionode) iukgmap = stdout
    !
    DO ig1 = 1, ngm_g
      WRITE(iukgmap, '(9i10)') (gmap(ig1, ig0), ig0 = 1, ng0vec)
    ENDDO
    !
    IF (iukgmap /= stdout) CLOSE(iukgmap)
    WRITE(stdout, *)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE refold
    !-----------------------------------------------------------------------
    ! 
    !---------------------------------
    SUBROUTINE backtoBZ(xx, yy, zz, n1, n2, n3)
    !---------------------------------
    !!
    !!  Brings xx, yy, and zz  into first BZ 
    !!
    !---------------------------------
    !
    USE kinds,  ONLY : DP
    !
    IMPLICIT NONE
    ! 
    INTEGER, INTENT(in) :: n1, n2, n3
    !! cell size
    REAL(KIND = DP), INTENT(inout) :: xx, yy, zz
    !! kgrid
    ! 
    ! Local variables
    INTEGER :: ib
    !! Size of replicas
    !
    ! More translations are needed to go back to the first BZ when the unit cell
    ! is far from cubic
    !
    DO ib = -2, 0
      IF (NINT(xx) < ib * n1) xx = xx + (-ib + 1) * n1
      IF (NINT(yy) < ib * n2) yy = yy + (-ib + 1) * n2
      IF (NINT(zz) < ib * n3) zz = zz + (-ib + 1) * n3
    ENDDO
    DO ib = 2, 1, -1
      IF (NINT(xx) >= ib * n1) xx = xx - ib * n1
      IF (NINT(yy) >= ib * n2) yy = yy - ib * n2
      IF (NINT(zz) >= ib * n3) zz = zz - ib * n3
    ENDDO
    !
    !-------------------------------------------
    END SUBROUTINE backtoBZ
    !-------------------------------------------
    ! 
    !--------------------------------------------------------
    SUBROUTINE ktokpmq(xk, xq, sign, ipool, nkq, nkq_abs)
    !--------------------------------------------------------
    !!
    !! For a given k point in cart coord, find the index 
    !! of the corresponding (k + sign*q) point
    !!
    !! In the parallel case, determine also the pool number
    !! nkq is the in-pool index, nkq_abs is the absolute
    !! index
    !!
    !--------------------------------------------------------
    !
    USE kinds,          only : DP
    use pwcom,          ONLY : nkstot
    USE cell_base,      ONLY : at
    USE epwcom,         ONLY : nkc1, nkc2, nkc3
    use klist_epw,      ONLY : xk_cryst
    USE mp_global,      ONLY : nproc_pool, npool
    USE mp_images,      ONLY : nproc_image
    USE mp,             ONLY : mp_barrier, mp_bcast
    USE constants_epw,  ONLY : eps5
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: sign
    !! +1 for searching k+q, -1 for k-q
    INTEGER, INTENT(out) :: nkq
    !! The pool hosting the k+-q point    
    INTEGER, INTENT(out) :: nkq_abs
    !! the index of k+sign*q
    INTEGER, INTENT(out) :: ipool
    !! The pool hosting the k+sign*q point
    REAL(KIND = DP), INTENT(in) :: xk(3)
    !! coordinates of k points
    REAL(KIND = DP), INTENT(in) :: xq(3)
    !! Coordinates of q point
    !
    ! Local variables
    LOGICAL :: in_the_list
    !! Is it in the list
    LOGICAL :: found
    !! Is it found
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: n
    !! Mapping index of k+q on k
    INTEGER :: iks
    !! 
    INTEGER :: nkl
    !! 
    INTEGER :: nkr
    !! Nb of kpt per pool
    INTEGER :: jpool
    !! 
    INTEGER :: kunit
    !! 
    REAL(KIND = DP) :: xxk(3)
    !! Coords. of k-point
    REAL(KIND = DP) :: xxq(3)
    !! Coords. of q-point
    REAL(KIND = DP) :: xx, yy, zz
    !! current k and k+q points in crystal coords. in multiple of nkc1, nkc2, nkc3
    REAL(KIND = DP) :: xx_c, yy_c, zz_c
    !! k-points in crystal coords. in multiple of nkc1, nkc2, nkc3
    !
    IF (ABS(sign) /= 1) CALL errore('ktokpmq', 'sign must be +1 or -1', 1)
    !
    ! bring k and q in crystal coordinates
    !
    xxk = xk
    xxq = xq
    !
    CALL cryst_to_cart(1, xxk, at, -1)
    CALL cryst_to_cart(1, xxq, at, -1)
    !
    !  check that k is actually on a uniform mesh centered at gamma
    !
    xx = xxk(1) * nkc1
    yy = xxk(2) * nkc2
    zz = xxk(3) * nkc3
    in_the_list = ABS(xx - NINT(xx)) <= eps5 .AND. &
                  ABS(yy - NINT(yy)) <= eps5 .AND. &
                  ABS(zz - NINT(zz)) <= eps5
    IF (.NOT. in_the_list) CALL errore('ktokpmq', 'is this a uniform k-mesh?', 1)
    !
    IF (xx < -eps5 .OR. yy < -eps5 .OR. zz < -eps5) THEN
      CALL errore('ktokpmq', 'coarse k-mesh needs to be strictly positive in 1st BZ', 1)
    ENDIF
    !
    ! now add the phonon wavevector and check that k+q falls again on the k grid
    !
    xxk = xxk + DBLE(sign) * xxq
    !
    xx = xxk(1) * nkc1
    yy = xxk(2) * nkc2
    zz = xxk(3) * nkc3
    in_the_list = ABS(xx - NINT(xx)) <= eps5 .AND. &
                  ABS(yy - NINT(yy)) <= eps5 .AND. &
                  ABS(zz - NINT(zz)) <= eps5
    IF (.NOT. in_the_list) CALL errore('ktokpmq', 'k+q does not fall on k-grid', 1)
    !
    ! Find the index of this k+q in the k-grid
    ! make sure xx, yy and zz are in 1st BZ
    !
    CALL backtoBZ(xx, yy, zz, nkc1, nkc2, nkc3)
    !
    n = 0
    found = .FALSE.
    DO ik = 1, nkstot
      xx_c = xk_cryst(1, ik) * nkc1
      yy_c = xk_cryst(2, ik) * nkc2
      zz_c = xk_cryst(3, ik) * nkc3
      !
      ! Check that the k-mesh was defined in the positive region of 1st BZ
      !
      IF (xx_c < -eps5 .OR. yy_c < -eps5 .OR. zz_c < -eps5) THEN
        CALL errore('ktokpmq', 'coarse k-mesh needs to be strictly positive in 1st BZ', 1)
      ENDIF
      !
      found = NINT(xx_c) == NINT(xx) .AND. &
              NINT(yy_c) == NINT(yy) .AND. &
              NINT(zz_c) == NINT(zz)
      IF (found) THEN  
        n = ik
        EXIT
      ENDIF
    ENDDO
    !
    ! 26/06/2012 RM
    ! since coarse k- and q- meshes are commensurate, one can easily find n
    ! n = NINT(xx) * nk2 * nk3 + NINT(yy) * nk3 + NINT(zz) + 1
    !
    IF (n == 0) CALL errore('ktokpmq', 'problem indexing k+q', 1)
    !
    ! Now n represents the index of k+sign*q in the original k grid.
    ! In the parallel case we have to find the corresponding pool 
    ! and index in the pool
    !
#if defined(__MPI)
    !
    npool = nproc_image / nproc_pool
    kunit = 1
    !
    DO jpool = 0, npool - 1
      !
      nkl = kunit * (nkstot / npool)
      nkr = (nkstot - nkl * npool) / kunit
      !
      ! The reminder goes to the first nkr pools (0...nkr-1)
      !
      IF (jpool < nkr ) nkl = nkl + kunit
      !
      !  the index of the first k point in this pool
      !
      iks = nkl * jpool + 1
      IF (jpool >= nkr ) iks = iks + nkr * kunit
      !
      IF (n >= iks) THEN
        ipool = jpool + 1
        nkq = n - iks + 1
      ENDIF
      !
    ENDDO
    !
#else
    ipool = 1
    nkq = n
#endif
    nkq_abs = n
    !---------------------------------------------------------------------
    END SUBROUTINE ktokpmq
    !---------------------------------------------------------------------
    ! 
  !-----------------------------------------------------------------------
  END MODULE kfold
  !-----------------------------------------------------------------------
