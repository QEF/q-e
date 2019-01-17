  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                        
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
  SUBROUTINE createkmap( xq )
  !-----------------------------------------------------------------------
  !!  
  !!  This subroutine is called from elphon_shuffle_wrap for each
  !!  nq1*nq2*nq3 phonon on the coarse mesh.    
  !!
  !!  It folds the k+q mesh into the k mesh using 5^3 G_0 translations 
  !!
  !!  SP - 2016 - iverbosity cannot be tested here. Generates Tb of data ... 
  !
  USE kinds,         ONLY : DP
  USE cell_base,     ONLY : at, bg
  USE klist,         ONLY : nkstot, xk
  USE start_k,       ONLY : nk1, nk2, nk3
  USE io_files,      ONLY : prefix
  USE io_epw,        ONLY : iukmap
  USE klist_epw,     ONLY : kmap
  USE kfold,         ONLY : g0vec_all, ng0vec, shift, g0vec_all_r
  USE io_global,     ONLY : meta_ionode
  USE mp,            ONLY : mp_barrier
  USE mp_world,      ONLY : world_comm
  USE elph2,         ONLY : xkq
  USE constants_epw, ONLY : eps5, zero
  !
  IMPLICIT NONE
  !
  REAL(kind=DP), INTENT(in) :: xq(3)
  !! Coords. of q-point 
  !
  ! Local variables
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
  !! G_0 in integer coords.
  INTEGER :: ig0
  !! Index of G_0 such that k+q+G_0 belongs to the 1st BZ
  INTEGER :: n
  !! Mapping index of k+q on k
  !
  REAL(kind=DP) :: xk_q(3)
  !! Coords. of k+q-point
  REAL(kind=DP) :: xx_c(nkstot), yy_c(nkstot), zz_c(nkstot)
  !! k-points in crystal coords. in multiple of nk1, nk2, nk3
  REAL(kind=DP) :: xx, yy, zz
  !! k+q in crystal coords. in multiple of nk1, nk2, nk3
  REAL(kind=DP) :: xx_n, yy_n, zz_n
  !! k+q in crystal coords. in multiple of nk1, nk2, nk3 in 1st BZ
  !
  LOGICAL :: in_the_list, found
  !
  IF (.not. ALLOCATED(xkq) ) ALLOCATE( xkq(3,nkstot) )
  xkq(:,:) = zero
  !
  IF (meta_ionode) THEN
    !
    !  the first proc keeps a copy of all kpoints !
    !
    IF ( .not. ALLOCATED(shift) ) ALLOCATE( shift(nkstot) )
    shift(:) = 0
    !
    !  Now fold k+q back into the k-grid for wannier interpolation.
    !  Since this is done before divide and impera, every pool has all the kpoints.
    !
    !  bring q in crystal coordinates and check commensuration
    !  loosy tolerance: not important since k+q is defined through nint() 
    !  bring q-point from cartesian to crystal coords.  
    !
    CALL cryst_to_cart(1, xq, at, -1)
    !
    xx = xq(1) * nk1 
    yy = xq(2) * nk2 
    zz = xq(3) * nk3 
    in_the_list = abs(xx-nint(xx)) .le. eps5 .AND. &
                  abs(yy-nint(yy)) .le. eps5 .AND. &
                  abs(zz-nint(zz)) .le. eps5
    IF (.not.in_the_list) CALL errore('createkmap','q-vec not commensurate',1)
    !
    ng0vec = 0
    DO ig1 = -2, 2
      DO ig2 = -2, 2
        DO ig3 = -2, 2
          ng0vec = ng0vec + 1
          g0vec_all(1,ng0vec) = ig1
          g0vec_all(2,ng0vec) = ig2
          g0vec_all(3,ng0vec) = ig3
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
      xx_c(ik) = xk(1,ik) * nk1
      yy_c(ik) = xk(2,ik) * nk2
      zz_c(ik) = xk(3,ik) * nk3
      in_the_list = abs(xx_c(ik)-nint(xx_c(ik))) .le. eps5 .AND. &
                    abs(yy_c(ik)-nint(yy_c(ik))) .le. eps5 .AND. &
                    abs(zz_c(ik)-nint(zz_c(ik))) .le. eps5
      IF (.not.in_the_list) CALL errore('createkmap','is this a uniform k-mesh?',1)
      !
      IF ( (xx_c(ik) .lt. -eps5) .OR. (yy_c(ik) .lt. -eps5) .OR. (zz_c(ik) .lt. -eps5) ) &
        CALL errore('createkmap','coarse k-mesh needs to be strictly positive in 1st BZ',1)
    ENDDO
    !
    DO ik = 1, nkstot
      !
      !  now add the phonon wavevector and check that k+q falls again on the k grid
      !
      xk_q(:) = xk(:,ik) + xq(:)
      !
      xx = xk_q(1) * nk1
      yy = xk_q(2) * nk2
      zz = xk_q(3) * nk3
      in_the_list = abs(xx-nint(xx)) .le. eps5 .AND. &
                    abs(yy-nint(yy)) .le. eps5 .AND. &
                    abs(zz-nint(zz)) .le. eps5
      IF (.not.in_the_list) CALL errore('createkmap','k+q does not fall on k-grid',1)
      !
      !  find the index of this k+q in the k-grid
      !
      i = mod( nint(xx + 2*nk1), nk1 ) 
      j = mod( nint(yy + 2*nk2), nk2 ) 
      k = mod( nint(zz + 2*nk3), nk3 ) 
      !
      xx_n = xx
      yy_n = yy
      zz_n = zz
      !
      !  make sure xx_n, yy_n and zz_n are in 1st BZ
      !
      CALL backtoBZ( xx_n, yy_n, zz_n, nk1, nk2, nk3 )
      !
      n = 0
      found = .false.
      DO jk = 1, nkstot
         !
         found = nint(xx_c(jk)) .eq. nint(xx_n) .AND. &
                 nint(yy_c(jk)) .eq. nint(yy_n) .AND. &
                 nint(zz_c(jk)) .eq. nint(zz_n)
         IF (found) THEN
            n = jk
            EXIT
         ENDIF
      ENDDO
      !
      !  26/06/2012 RM
      !  since coarse k- and q- meshes are commensurate, one can easily find n
      !  n = nint(xx_n) * nk2 * nk3 + nint(yy_n) * nk3 + nint(zz_n) + 1
      !  n represents the index of k+q on the coarse k-grid.
      !
      IF (n .eq. 0) CALL errore('createkmap','problem indexing k+q',1)
      !
      kmap(ik) = n
      !
      !  determine the G_0 such that k+q+G_0 belongs to the first BZ
      !
      g0vec(1) = ( i - nint(xx) ) / nk1
      g0vec(2) = ( j - nint(yy) ) / nk2
      g0vec(3) = ( k - nint(zz) ) / nk3
      !
      !  now store the shift for this k+q point
      !
      in_the_list = .false.
      ig0 = 0
      DO WHILE ( (ig0.le.ng0vec) .AND. (.not.in_the_list) )
        ig0 = ig0 + 1
        in_the_list = ( (abs(g0vec(1) - g0vec_all(1,ig0)) .le. eps5) .AND. &
                        (abs(g0vec(2) - g0vec_all(2,ig0)) .le. eps5) .AND. &
                        (abs(g0vec(3) - g0vec_all(3,ig0)) .le. eps5))
      ENDDO
      shift(ik) = ig0
      !
      IF (.not.in_the_list) CALL errore &
         ('createkmap','cannot find the folding vector in the list',1)
      !
      !  obsolete:
      !
      !  very important: now redefine k+q through the corresponding kpoint on the k mesh
      !  Note that this will require using the periodic gauge in the calculation of the
      !  electron-phonon matrix elements (factor e^iG_0*r if G_0 is the vector used for 
      !  refolding)
      !
      xkq(:,ik) = xk(:,n) 
      !
    ENDDO
    !
    ! bring k-points, q-point, and G_0-vectors back to cartesian coordinates 
    !
    CALL cryst_to_cart(1, xq, bg, 1) 
    CALL cryst_to_cart(nkstot, xk, bg, 1)
    !
    g0vec_all_r = dble(g0vec_all)
    CALL cryst_to_cart(ng0vec, g0vec_all_r, bg, 1)
    !
    !  the unit with kmap(ik) and shift(ik)
    ! 
    OPEN(iukmap, file = TRIM(prefix)//'.kmap', form = 'formatted')
    DO ik = 1, nkstot
      WRITE(iukmap,'(3i6)') ik, kmap(ik), shift(ik) 
    ENDDO
    CLOSE(iukmap)
    !
  ENDIF
  CALL mp_barrier(world_comm)
  !
  RETURN
  !
  END SUBROUTINE createkmap

  !-----------------------------------------------------------------------
  SUBROUTINE createkmap2( xxq )
  !-----------------------------------------------------------------------
  !!
  !!  generate the map k+q --> k for folding the rotation matrix U(k+q) 
  !!  
  !!  in parallel case, this subroutine must be called only by first proc 
  !!  (which has all the kpoints)
  !!
  !-----------------------------------------------------------------------
  !
  USE kinds,         ONLY : DP
  USE cell_base,     ONLY : at, bg
  USE klist,         ONLY : nkstot, xk
  USE klist_epw,     ONLY : kmap  
  USE start_k,       ONLY : nk1, nk2, nk3
  USE elph2,         ONLY : xkq
  USE constants_epw, ONLY : eps5, zero

  IMPLICIT NONE
  !
  REAL(kind=DP), INTENT(in) :: xxq(3)
  !! The current q-point 
  ! 
  ! Local variables
  INTEGER :: ik
  !! K-point index
  INTEGER :: jk
  !! Another k-point index
  INTEGER :: n
  !! Mapping index of k+q on k
  !
  REAL(kind=DP) :: xx_c(nkstot), yy_c(nkstot), zz_c(nkstot)
  !! k-points in crystal coords. in multiple of nk1, nk2, nk3
  REAL(kind=DP) :: xx, yy, zz
  !! k+q in crystal coords. in multiple of nk1, nk2, nk3
  !
  LOGICAL :: in_the_list, found
  !
  ! the first proc keeps a copy of all kpoints !
  !
  ! bring q from cartesian to crystal coordinates and check commensuration
  ! 
  CALL cryst_to_cart(1, xxq, at, -1)
  !
  xx = xxq(1) * nk1 
  yy = xxq(2) * nk2 
  zz = xxq(3) * nk3 
  in_the_list = abs(xx-nint(xx)) .le. eps5 .AND. &
                abs(yy-nint(yy)) .le. eps5 .AND. &
                abs(zz-nint(zz)) .le. eps5
  IF (.not.in_the_list) CALL errore('createkmap2','q-vec not commensurate',1)
  IF (.not. ALLOCATED(xkq) ) ALLOCATE( xkq(3,nkstot) )
  xkq(:,:) = zero
  !
  !  bring all the k-points from cartesian to crystal coordinates 
  !
  CALL cryst_to_cart(nkstot, xk, at, -1)
  !
  DO ik = 1, nkstot
    !
    !  check that the k's are actually on a uniform mesh centered at gamma
    !
    xx_c(ik) = xk(1,ik) * nk1
    yy_c(ik) = xk(2,ik) * nk2
    zz_c(ik) = xk(3,ik) * nk3
    in_the_list = abs(xx_c(ik)-nint(xx_c(ik))) .le. eps5 .AND. &
                  abs(yy_c(ik)-nint(yy_c(ik))) .le. eps5 .AND. &
                  abs(zz_c(ik)-nint(zz_c(ik))) .le. eps5
    IF (.not.in_the_list) CALL errore('createkmap2','is this a uniform k-mesh?',1)
    !
    IF ( (xx_c(ik) .lt. -eps5) .OR. (yy_c(ik) .lt. -eps5) .OR. (zz_c(ik) .lt. -eps5) ) &
      CALL errore('createkmap2','coarse k-mesh needs to be strictly positive in 1st BZ',1)
  ENDDO
  !
  DO ik = 1, nkstot
    !
    !  now add the phonon wavevector and check that k+q falls again on the k grid
    !
    xkq(:,ik) = xk(:,ik) + xxq(:)
    !
    xx = xkq(1,ik) * nk1
    yy = xkq(2,ik) * nk2
    zz = xkq(3,ik) * nk3
    in_the_list = abs(xx-nint(xx)) .le. eps5 .AND. &
                  abs(yy-nint(yy)) .le. eps5 .AND. &
                  abs(zz-nint(zz)) .le. eps5
    IF (.not.in_the_list) CALL errore('createkmap2','k+q does not fall on k-grid',1)
    !
    !  find the index of this k+q in the k-grid
    !
    ! make sure xx, yy and zz are in 1st BZ
    !
    CALL backtoBZ( xx, yy, zz, nk1, nk2, nk3 )
    !
    n = 0
    found = .false.
    DO jk = 1, nkstot
       !
       found = nint(xx_c(jk)) .eq. nint(xx) .and. &
               nint(yy_c(jk)) .eq. nint(yy) .and. &
               nint(zz_c(jk)) .eq. nint(zz)
       IF (found) THEN
          n = jk
          EXIT
       ENDIF
    ENDDO
    !
    IF (n .eq. 0) CALL errore('createkmap2','problem indexing k+q',1)
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
  END SUBROUTINE createkmap2
  !
  !-------------------------------------------------------------------------
  SUBROUTINE createkmap_pw2
  !-------------------------------------------------------------------------
  !!
  !! Creates the first instance of [prefix].kgmap. 
  !!
  !-------------------------------------------------------------------------
  USE kinds,         ONLY : DP
  USE cell_base,     ONLY : at, bg
  USE start_k,       ONLY : nk1, nk2, nk3
  USE pwcom,         ONLY : nkstot
  USE epwcom,        ONLY : xk_cryst
  USE io_global,     ONLY : stdout, meta_ionode
  USE io_files,      ONLY : prefix
  USE io_epw,        ONLY : iukgmap
  USE gvect,         ONLY : ngm, ngm_g, gcutm
  USE fft_base,      ONLY : dfftp
  USE fft_types,     ONLY : fft_stick_index
  USE fft_ggen,      ONLY : fft_set_nl
  USE constants,     ONLY : eps8
  USE kfold
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
  !
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
  !
  INTEGER, ALLOCATABLE :: ig_l2g(:)
  !! Converts a local G-vector index into the global index 
  INTEGER, ALLOCATABLE :: g2l(:)
  !! Local index of G-vector 
  INTEGER, ALLOCATABLE :: mill_unsorted(:,:)
  !! Array of unsorted Miller indices of G-vectors
  INTEGER, ALLOCATABLE :: mill(:,:)
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
  !
  REAL(kind=DP) :: xx, yy, zz
  !! k-point in crystal coords. in multiple of nk1, nk2, nk3
  REAL(kind=DP), ALLOCATABLE :: gg(:)
  !! G^2 in increasing order
  REAL(kind=DP), ALLOCATABLE :: g(:,:)
  !! G-vectors cartesian components in increasing order of G^2
  REAL(kind=DP), ALLOCATABLE :: g2sort_g(:)
  !! G-vectors for the current processor
  REAL(kind=DP) :: tx(3), ty(3), t(3)
  REAL(kind=DP), ALLOCATABLE :: tt(:)
  !
  LOGICAL :: in_the_list
  !! Variable in the list
  LOGICAL :: is_local
  !! Local variable
  ! 
  IF (meta_ionode) THEN
    !
    WRITE(stdout, '(/5x,a)') 'Calculating kgmap'
    CALL flush(stdout)
    !
    OPEN(iukgmap,file = TRIM(prefix)//'.kgmap',form='formatted')
    ! 
    ! the 5^3 possible G_0 translations
    ng0vec = 0
    DO ig1 = -2, 2
      DO ig2 = -2, 2
        DO ig3 = -2, 2
          ng0vec = ng0vec + 1
          g0vec_all(1,ng0vec) = ig1
          g0vec_all(2,ng0vec) = ig2
          g0vec_all(3,ng0vec) = ig3
        ENDDO
      ENDDO
    ENDDO
    ig0 = nint( dble(ng0vec) / 2 )
    !
    IF (.not. ALLOCATED(shift)) ALLOCATE( shift(nkstot) )
    ! 
    DO ik = 1, nkstot       
      !
      xx = xk_cryst(1,ik) * nk1
      yy = xk_cryst(2,ik) * nk2
      zz = xk_cryst(3,ik) * nk3
      ! check that the k-mesh was defined in the positive region of 1st BZ
      !
      IF ( (xx .lt. -eps5) .OR. (yy .lt. -eps5) .OR. (zz .lt. -eps5) ) &
         CALL errore('createkmap_pw2','coarse k-mesh needs to be strictly positive in 1st BZ',1)
      ! 
      shift(ik) = ig0
      WRITE(iukgmap,'(3i6)') ik, shift(ik)
      !
    ENDDO
    IF (ALLOCATED(shift)) DEALLOCATE(shift)
    !
    g0vec_all_r = dble(g0vec_all)
    ! bring G_0 vectors from crystal to cartesian coordinates
    CALL cryst_to_cart(ng0vec, g0vec_all_r, bg, 1)
    !
    WRITE(iukgmap,'(i5)') ng0vec
    DO ig0 = 1, ng0vec
      WRITE(iukgmap,'(3f20.15)') g0vec_all_r(:,ig0)
    ENDDO
    !
  ENDIF
  !   
  CALL mp_barrier(inter_pool_comm)
  CALL mp_barrier(inter_image_comm)
  !
  ! RM: The following is adapted from ggen subroutine in Modules/recvec_subs.f90
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
  !    and computes all the g vectors inside a sphere
  !
  ALLOCATE( mill_unsorted(3,ngm_save) )
  ALLOCATE( igsrt(ngm_max) )
  ALLOCATE( g2l(ngm_max) )
  ALLOCATE( g2sort_g(ngm_max) )
  ALLOCATE( ig_l2g(ngm_max) )
  ALLOCATE( mill(3,ngm_max) )
  ALLOCATE( jtoi(ngm_max) )
  ALLOCATE( itoj(ngm_max) )
  ALLOCATE( g(3,ngm_max) )
  ALLOCATE( gg(ngm_max) )
  !
  !    Set the total number of FFT mesh points and and initial value of gg.
  !    The choice of gcutm is due to the fact that we have to order the
  !    vectors after computing them
  !
  gg(:) = gcutm + 1.d0
  !
  g2sort_g(:) = 1.0d20
  !
  ! allocate temporal array
  !
  ALLOCATE( tt(dfftp%nr3) )
  !
  ! max miller indices (same convention as in module stick_set)
  !
  ni = (dfftp%nr1-1)/2
  nj = (dfftp%nr2-1)/2
  nk = (dfftp%nr3-1)/2
  !
  istart = -ni
  DO i = istart, ni
    !
    tx(1:3) = i * bg(1:3,1)
    !
    jstart = -nj
    DO j = jstart, nj
      !
      IF ( dfftp%lpara .AND. fft_stick_index(dfftp,i,j) == 0) THEN
        is_local = .FALSE.
      ELSE
        is_local = .TRUE.
      ENDIF
      !
      ty(1:3) = tx(1:3) + j * bg(1:3,2)
      !
      !  compute all the norm square
      !
      kstart = -nk
      DO k = kstart, nk
        !
        t(1) = ty(1) + k * bg(1,3)
        t(2) = ty(2) + k * bg(2,3)
        t(3) = ty(3) + k * bg(3,3)
        tt(k-kstart+1) = t(1)**2 + t(2)**2 + t(3)**2
      ENDDO
      !
      !  save all the norm square within cutoff
      !
      DO k = kstart, nk
        IF (tt(k-kstart+1) <= gcutm) THEN
          ngm = ngm + 1
          IF (ngm > ngm_max) CALL errore('createkmap_pw2 1', 'too many g-vectors', ngm)
          IF ( tt(k-kstart+1) > eps8 ) THEN
            g2sort_g(ngm) = tt(k-kstart+1)
          ELSE
            g2sort_g(ngm) = 0.d0
          ENDIF
          IF (is_local) THEN
            ngm_local = ngm_local + 1
            mill_unsorted(:,ngm_local) = (/ i,j,k /)
            g2l(ngm) = ngm_local
          ELSE
            g2l(ngm) = 0
          ENDIF
        ENDIF
      ENDDO
    ENDDO !jloop
  ENDDO !iloop
  IF (ngm  /= ngm_max) &
      CALL errore('createkmap_pw2', 'g-vectors missing !', abs(ngm - ngm_max))
  !
  igsrt(1) = 0
  CALL hpsort_eps( ngm_g, g2sort_g, igsrt, eps8 )
  DEALLOCATE( g2sort_g, tt )
  !  
  ngm = 0
  !
  DO ng = 1, ngm_max
     !
     IF ( g2l(igsrt(ng))>0 ) THEN
        ! fetch the indices
        i = mill_unsorted(1,g2l(igsrt(ng)))
        j = mill_unsorted(2,g2l(igsrt(ng)))
        k = mill_unsorted(3,g2l(igsrt(ng)))
        !
        ngm = ngm + 1
        !
        ig_l2g(ngm) = ng
        ! 
        g(1:3,ngm) = i * bg(:,1) + j * bg(:,2) + k * bg(:,3)
        gg(ngm) = sum(g(1:3,ngm)**2)
     ENDIF
  ENDDO !ngloop
  DEALLOCATE( g2l )
  !
  IF (ngm /= ngm_save) &
     CALL errore('createkmap_pw2', 'g-vectors (ngm) missing !', abs(ngm - ngm_save))
  !
  CALL fft_set_nl( dfftp, at, g, mill )
  !
  DO i = 1, ngm_g
     jtoi(i) = igsrt(i)
  ENDDO !
  !
  DO i = 1, ngm_g
     itoj(jtoi(i)) = i
  ENDDO
  !
  CALL refold( ngm_g, mill, itoj, jtoi )
  !
  CALL mp_barrier(inter_pool_comm)
  CALL mp_barrier(inter_image_comm)
  !
  DEALLOCATE( ig_l2g, mill, mill_unsorted, igsrt, jtoi, itoj, g, gg )
  !
  RETURN
  !
  END SUBROUTINE createkmap_pw2
  !-------------------------------------------------------------------------
  ! 
  !-----------------------------------------------------------------------
  SUBROUTINE refold( ngm_g, mill_g, itoj, jtoi )
  !----------------------------------------------------------------------
  !
  !   Map the indices of G+G_0 into those of G 
  !   this is used to calculate electron-phonon matrix elements by
  !   refolding the k+q points into the first BZ (original k grid)
  !
  !   No parallelization on G-vecs at the moment  
  !   (actually this is done on the global array, but in elphel2.f90
  !   every processor has just a chunk of the array, I may need some
  !   communication)
  !
  !   No ultrasoft now
  !
  !   I use the rule : if not found then gmap = 0 
  !   Note that the map will be used only up to npwx (small sphere), 
  !   while the G-vectors lost in the process are on the surface of 
  !   the large sphere (density set).
  !
  !-----------------------------------------------------------------
  USE io_global,     ONLY : stdout, meta_ionode
  USE io_epw,        ONLY : iukgmap
! SP: Sucidal. Produce too much data. Only use for debugging. 
!  USE control_flags, ONLY : iverbosity
  USE kfold
  !
  IMPLICIT NONE
  !
  INTEGER :: ngm_g
  !! Counter on G-vectors
  INTEGER :: mill_g(3,ngm_g)
  !!  Array of Miller indices of G-vectors in increasing order of G^2
  INTEGER :: jtoi(ngm_g)
  !! For the i-th G-vector in the sorted list, jtoi(i)
  !! returns its index in the unsorted list
  INTEGER :: itoj(ngm_g)
  !! itoj(i) returns the index of the G-vector in the sorted list
  !! that was at i-th position in the unsorted list
  INTEGER :: ig0 
  !! Counter on G_0 vectors
  INTEGER :: ig1, ig2
  !! Counter on G vectors
  INTEGER :: i, j, k
  !! Miller indices for G+G_0 vector
  INTEGER :: ig1_use, ig2_use, ig2_guess, notfound, guess_skip
  !! Temporary G-vectors indices
  INTEGER :: indold, indnew
  !! Counter on G_0 vectors indices for writing to file
  !
  LOGICAL :: tfound
  !
  ALLOCATE( gmap(ngm_g,ng0vec) )
  gmap(:,:) = 0
  guess_skip = 0
  !
  !  Loop on the inequivalent G_0 vectors
  !
  DO ig0 = 1, ng0vec
    !
    IF (ig0 .eq. 1) THEN
      WRITE(stdout,'(/5x,"Progress kgmap: ")',advance='no')
      indold = 0
    ENDIF
    indnew = nint( dble(ig0) / dble(ng0vec) * 40 )
    IF (indnew.ne.indold) WRITE(stdout,'(a)',advance='no') '#'
    indold = indnew
    !
    notfound = 0
    DO ig1 = 1, ngm_g
      !
      ig1_use = itoj(ig1)
      !
      !  the initial G vector
      !
      i = mill_g(1,ig1_use)
      j = mill_g(2,ig1_use)
      k = mill_g(3,ig1_use)
      !
      !  the final G+G_0 vector
      !
      i = i + g0vec_all(1,ig0)
      j = j + g0vec_all(2,ig0)
      k = k + g0vec_all(3,ig0)
      !
      ig2 = 0
      tfound = .false.
      !
      ! try to guess next index
      !
      ig2_guess = jtoi(ig1_use) + guess_skip
      !
      IF ((ig2_guess .gt. 0) .AND. (ig2_guess .lt. ngm_g+1)) THEN
        !
        ig2_guess = itoj(ig2_guess)
        !
        IF ((i .eq. mill_g(1,ig2_guess)) .AND. (j .eq. mill_g(2,ig2_guess)) .AND. (k .eq. mill_g(3,ig2_guess))) THEN
          !
          ig2_use = ig2_guess
          tfound = .true.
          !
        ENDIF
        !
      ENDIF
      !
      DO WHILE ((.not. tfound) .AND. (ig2 .lt. ngm_g))
        !
        ig2 = ig2 + 1
        ig2_use = itoj(ig2)
        tfound = (i .eq. mill_g(1,ig2_use)) .AND. & 
                 (j .eq. mill_g(2,ig2_use)) .AND. & 
                 (k .eq. mill_g(3,ig2_use))
        !
      ENDDO
      !
      IF (tfound) THEN
        gmap(ig1_use,ig0) = ig2_use
        guess_skip = jtoi(ig2_use) - jtoi(ig1_use)
      ELSE
        gmap(ig1_use,ig0) = 0
        notfound = notfound + 1
      ENDIF
      !
    ENDDO
    !
  ENDDO
  ! 
  !  output on file for electron-phonon matrix elements
  !
  IF (.NOT. meta_ionode) iukgmap = stdout
  !
  DO ig1 = 1, ngm_g
    WRITE(iukgmap,'(9i10)') (gmap(ig1,ig0), ig0 = 1, ng0vec)
  ENDDO
  !
  IF (iukgmap .ne. stdout) CLOSE(iukgmap)
  WRITE(stdout,*)
  !
  RETURN
  !
  END SUBROUTINE refold
