  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
SUBROUTINE loadkmesh_para
  !-----------------------------------------------------------------------
  !!
  !!  load fine k mesh and distribute among pools
  !!
  !-----------------------------------------------------------------------
  USE io_global, ONLY : ionode_id
  USE mp_global, ONLY : inter_pool_comm, my_pool_id, npool
  USE mp,        ONLY : mp_bcast, mp_sum
  USE mp_world,  ONLY : mpime 
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout
  USE epwcom,    ONLY : filkf, nkf1, nkf2, nkf3, &
                        rand_k, rand_nk, mp_mesh_k, system_2d
  USE elph2,     ONLY : nkqtotf, nkqf, xkf, wkf, nkf
  USE cell_base, ONLY : at, bg
  USE symm_base, ONLY : s, t_rev, time_reversal, set_sym_bl, nrot
  USE io_epw,    ONLY : iunkf

  USE noncollin_module, ONLY : noncolin
  !
  implicit none
  !
  real(kind=DP), ALLOCATABLE :: xkf_(:,:), wkf_(:), xkf_tmp(:,:), &
       wkf_tmp(:)
  integer :: ik, ikk, ikq, lower_bnd, upper_bnd, i, j, k, ios
  !
  !
  integer :: rest
  !
  IF (mpime .eq. ionode_id) THEN
    IF (filkf .ne. '') THEN ! load from file (crystal coordinates)
       !
       WRITE(stdout, *) '     Using k-mesh file: ', trim(filkf)
       OPEN( unit = iunkf, file = filkf, status = 'old', form = 'formatted',err=100, iostat=ios)
100    CALL errore('loadkmesh_para','opening file '//filkf,abs(ios))
       READ(iunkf, *) nkqtotf 
       !
       ALLOCATE (xkf_ (3, 2*nkqtotf), wkf_(2*nkqtotf))
       !
       DO ik = 1, nkqtotf
          !
          ikk = 2 * ik - 1
          ikq = ikk + 1
          !
          READ (iunkf, *) xkf_ (:, ikk ), wkf_ (ikk)
          !
          ! SP: This is so we can input a weight of 1 to random file 
          !     This way you can feed the same file for the k and q grid  
          wkf_ (ikk) = wkf_ (ikk)*2.d0  
          !
          !  bring the k point to crystal coordinates
          ! CALL cryst_to_cart ( 1, xkf_ (:,ikk), at, -1)
          !
          xkf_ (:, ikq) = xkf_ (:, ikk) 
          wkf_ ( ikq ) = 0.d0
          !
       ENDDO
       CLOSE(iunkf)
       !
       ! redefine nkqtotf to include the k+q points
       !
       nkqtotf = 2 * nkqtotf
       !
    ELSEIF ( (nkf1.ne.0) .and. (nkf2.ne.0) .and. (nkf3.ne.0) ) THEN ! generate grid
       IF (mp_mesh_k) THEN
           ! get size of the mp_mesh in the irr wedge 
          WRITE (stdout, '(a,3i4)') '     Using uniform MP k-mesh: ', nkf1, nkf2, nkf3
          call set_sym_bl ( )
          !
          ALLOCATE ( xkf_ (3, 2*nkf1*nkf2*nkf3), wkf_(2*nkf1*nkf2*nkf3) )
          ! the result of this call is just nkqtotf
          CALL kpoint_grid ( nrot, time_reversal, .false., s, t_rev, bg, nkf1*nkf2*nkf3, &
               0,0,0, nkf1,nkf2,nkf3, nkqtotf, xkf_, wkf_)
          DEALLOCATE (xkf_, wkf_)
          ALLOCATE ( xkf_ (3, 2*nkqtotf), wkf_(2*nkqtotf)) 
          ALLOCATE (xkf_tmp (3,nkqtotf), wkf_tmp(nkqtotf))
          CALL kpoint_grid ( nrot, time_reversal, .false., s, t_rev, bg, nkf1*nkf2*nkf3, &
               0,0,0, nkf1,nkf2,nkf3, nkqtotf, xkf_tmp, wkf_tmp)
          !  
          ! assign to k and k+q for xkf and wkf 
          ! 
          DO ik = 1, nkqtotf
             ikk = 2 * ik - 1
             ikq = ikk + 1
             xkf_(:,ikk) = xkf_tmp(:,ik)
             xkf_(:,ikq) = xkf_tmp(:,ik)
             wkf_(ikk)   = 2.d0 * wkf_tmp(ik)
             wkf_(ikq)   = 0.d0
          ENDDO
          DEALLOCATE (xkf_tmp, wkf_tmp)
          !       
          ! bring the k point to crystal coordinates       
          CALL cryst_to_cart (2*nkqtotf, xkf_, at, -1)
          !
          ! redefine nkqtotf to include the k+q points
          !
          nkqtotf = 2 * nkqtotf
          !
       ELSE
          !
          WRITE (stdout, '(a,3i4)') '     Using uniform k-mesh: ', nkf1, nkf2, nkf3
          !
          nkqtotf = 2 * nkf1 * nkf2 * nkf3
          ALLOCATE ( xkf_ (3, nkqtotf), wkf_(nkqtotf) )
          wkf_(:) = 0.d0
          DO ik = 1, nkf1 * nkf2 * nkf3
             wkf_(2*ik-1) = 2.d0/(dble(nkqtotf/2))
          ENDDO
          DO i = 1, nkf1
             DO j = 1, nkf2
                DO k = 1, nkf3
                   ik = (i-1)*nkf2*nkf3 + (j-1)*nkf3 + k
                   ikk = 2 * ik - 1
                   ikq = ikk + 1
                   xkf_(1, ikk) = dble(i-1)/dble(nkf1)
                   xkf_(2, ikk) = dble(j-1)/dble(nkf2)
                   xkf_(3, ikk) = dble(k-1)/dble(nkf3)
                   xkf_(1, ikq) = xkf_(1, ikk)
                   xkf_(2, ikq) = xkf_(2, ikk)
                   xkf_(3, ikq) = xkf_(3, ikk) 
                ENDDO
             ENDDO
          ENDDO
          !
       ENDIF
       !
    ELSEIF (rand_k) THEN  ! random points
       !
       WRITE (stdout, *) '    Using random k-mesh: ', rand_nk
       !
       nkqtotf = rand_nk
       ALLOCATE (xkf_ (3, 2*nkqtotf), wkf_(2*nkqtotf))
       !
       CALL init_random_seed()
       !
       DO ik = 1, nkqtotf
          !
          ikk = 2 * ik - 1
          ikq = ikk + 1
          !
          wkf_(ikk) = 2.d0/ dble(nkqtotf)
          wkf_(ikq) = 0.d0
          !
          IF ( system_2d ) THEN
             CALL random_number(xkf_(1:2,ikk))
             xkf_(3,ikk) = 0.d0
          ELSE
             CALL random_number(xkf_(:,ikk))
          ENDIF
          !
          xkf_(:,ikq) = xkf_(:,ikk)
          !  bring the k point to cartesian coordinates
          !CALL cryst_to_cart ( 1, xkf_ (:,ikk), bg, 1)
          !WRITE(stdout,'(4f12.7)') xkf_(:,ikk), wkf_(ikk)
          !  bring the k point to crystal coordinates
          !CALL cryst_to_cart ( 1, xkf_ (:,ikk), at, -1)
       ENDDO
       !WRITE(stdout,'(a)') '  '
       !
       ! redefine nkqtotf to include the k+q points
       !
       nkqtotf = 2 * nkqtotf
       !
    ELSE ! don't know how to get grid
       CALL errore('loadkmesh_para', "Cannot load fine k points", 1)
    ENDIF
 ENDIF
 !
#if defined(__MPI)
 CALL mp_bcast (nkqtotf, ionode_id, inter_pool_comm)
 !
 !  scatter the k points of the fine mesh across the pools
 !
 nkqf = 2 * ( nkqtotf / 2 / npool )
 rest = ( nkqtotf - nkqf * npool ) / 2
 IF (my_pool_id < rest ) THEN
    nkqf = nkqf + 2
    lower_bnd = my_pool_id*nkqf + 1
    upper_bnd = lower_bnd + nkqf - 1
 ELSE
    lower_bnd = rest*(nkqf+2)+(my_pool_id-rest)*nkqf + 1
    upper_bnd = lower_bnd + nkqf - 1
 ENDIF
 !
 nkf = nkqf / 2 
 IF (.not.ALLOCATED(xkf_)) ALLOCATE (xkf_(3,nkqtotf))
 IF (.not.ALLOCATED(wkf_)) ALLOCATE (wkf_(  nkqtotf))
 CALL mp_bcast(xkf_, ionode_id, inter_pool_comm)
 CALL mp_bcast(wkf_, ionode_id, inter_pool_comm)
 !
#else
 !
 ! In serial the definitions are much easier 
 !
 nkqf = nkqtotf
 nkf = nkqf / 2 
 lower_bnd = 1
 upper_bnd = nkqf
 !
#endif
 !
 !  Assign the weights and vectors to the correct bounds
 !
 ALLOCATE(xkf(3,nkqf))
 ALLOCATE(wkf(  nkqf))
 xkf(:,:) = xkf_ (:, lower_bnd:upper_bnd)
 IF (noncolin) THEN 
    wkf(  :) = wkf_ ( lower_bnd:upper_bnd)/2.d0
 ELSE
    wkf(  :) = wkf_ ( lower_bnd:upper_bnd)
 ENDIF  
 !
 IF (abs(sum (wkf_ (:)) - 2.d0) .gt. 1.d-4 ) &
    WRITE(stdout,'(5x,"WARNING: k-point weigths do not add up to 1 [loadkmesh_para]")')
 !
 WRITE( stdout, '(5x,"Size of k point mesh for interpolation: ",i10)' ) nkqtotf 
 WRITE( stdout, '(5x,"Max number of k points per pool:",7x,i10)' ) nkqf 
 !
 IF (ALLOCATED(xkf_)) DEALLOCATE(xkf_)
 IF (ALLOCATED(wkf_)) DEALLOCATE(wkf_)
 !
END SUBROUTINE loadkmesh_para
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE loadkmesh_serial
!-----------------------------------------------------------------------
!!
!!  Load fine k mesh
!!
!-----------------------------------------------------------------------
  USE io_global, ONLY : ionode_id
  USE mp_global, ONLY : inter_pool_comm
  USE mp,        ONLY : mp_bcast
  USE mp_world,  ONLY : mpime
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout
  USE epwcom,    ONLY : filkf, nkf1, nkf2, nkf3, &
                        rand_k, rand_nk, mp_mesh_k, system_2d
  USE elph2,     ONLY : xkf, wkf, nkqtotf, nkf, nkqf
  USE cell_base, ONLY : at, bg
  USE symm_base, ONLY : s, t_rev, time_reversal, set_sym_bl, nrot
  USE io_epw,    ONLY : iunkf
  implicit none
  !
  integer :: ik, i, j, k, ios, ikk, ikq
  real(kind=DP), ALLOCATABLE ::  xkf_tmp(:,:), wkf_tmp(:)
  !
  IF (mpime .eq. ionode_id) THEN
    IF (filkf .ne. '') THEN ! load from file (crystal coordinates)
       !
       ! Each pool gets its own copy from the action=read statement
       !
       WRITE (stdout, *) '     Using k-mesh file: ', trim(filkf)
       OPEN ( unit = iunkf, file = filkf, status = 'old', form = 'formatted', err=100, iostat=ios)
100    CALL errore('loadkmesh_serial','opening file '//filkf,abs(ios))
       READ(iunkf, *) nkqtotf
       ALLOCATE (xkf(3, 2*nkqtotf), wkf(2*nkqtotf))
       !DO ik = 1, nkqtotf
       !   READ (iunkf, *) xkf (:, ik), wkf(ik)
       !ENDDO
       DO ik = 1, nkqtotf
         !
         ikk = 2 * ik - 1
         ikq = ikk + 1
         !
         READ (iunkf, *) xkf (:, ikk ), wkf (ikk)
         !
         ! SP: This is so we can input a weight of 1 to random file 
         !     This way you can feed the same file for the k and q grid  
         wkf (ikk) = wkf (ikk)*2.d0
         !
         !  bring the k point to crystal coordinates
         ! CALL cryst_to_cart ( 1, xkf_ (:,ikk), at, -1)
         !
         xkf (:, ikq) = xkf (:, ikk)
         wkf ( ikq ) = 0.d0
         !
       ENDDO
       CLOSE(iunkf)
       !
       ! redefine nkqtotf to include the k+q points
       !
       nkqtotf = 2 * nkqtotf
       !
       !
       ! bring xkf in crystal coordinates
       ! CALL cryst_to_cart (nkqtotf, xkf, at, -1)
       !
    ELSEIF ( (nkf1.ne.0) .and. (nkf2.ne.0) .and. (nkf3.ne.0) ) THEN ! generate grid
       IF (mp_mesh_k) THEN
          ! get size of the mp_mesh in the irr wedge 
          WRITE (stdout, '(a,3i4)') '     Using uniform k-mesh: ', nkf1, nkf2, nkf3
          call set_sym_bl ( )
          !                                         
          ALLOCATE ( xkf (3, 2*nkf1*nkf2*nkf3), wkf(2*nkf1*nkf2*nkf3) )
          ! the result of this call is just nkqtotf
          CALL kpoint_grid ( nrot, time_reversal, s, t_rev, bg, nkf1*nkf2*nkf3, &
               0,0,0, nkf1,nkf2,nkf3, nkqtotf, xkf, wkf)
          DEALLOCATE ( xkf, wkf) 
          ALLOCATE ( xkf(3, 2*nkqtotf), wkf(2*nkqtotf))
          ALLOCATE (xkf_tmp (3,nkqtotf), wkf_tmp(nkqtotf))
          CALL kpoint_grid ( nrot, time_reversal, s, t_rev, bg, nkf1*nkf2*nkf3, &
               0,0,0, nkf1,nkf2,nkf3, nkqtotf, xkf_tmp, wkf_tmp)
          !  
          ! assign to k and k+q for xkf and wkf 
          ! 
          DO ik = 1, nkqtotf
             ikk = 2 * ik - 1
             ikq = ikk + 1
             xkf(:,ikk) = xkf_tmp(:,ik)
             xkf(:,ikq) = xkf_tmp(:,ik)
             wkf(ikk)   = 2.d0 * wkf_tmp(ik)
             wkf(ikq)   = 0.d0
          ENDDO
          DEALLOCATE (xkf_tmp, wkf_tmp)
          !       
          ! bring the k point to crystal coordinates       
          CALL cryst_to_cart (2*nkqtotf, xkf, at, -1)
          !
          ! redefine nkqtotf to include the k+q points
          !
          nkqtotf = 2 * nkqtotf
          !
       ELSE
          WRITE (stdout, '(a,3i4)') '     Using uniform k-mesh: ', nkf1, nkf2, nkf3
          !
          nkqtotf = 2 * nkf1 * nkf2 * nkf3
          ALLOCATE ( xkf(3, nkqtotf), wkf(nkqtotf) )
          wkf(:) = 0.d0
          DO ik = 1, nkf1 * nkf2 * nkf3
             wkf(2*ik-1) = 2.d0/(dble(nkqtotf/2))
          ENDDO
          DO i = 1, nkf1
             DO j = 1, nkf2
                DO k = 1, nkf3
                   ik = (i-1)*nkf2*nkf3 + (j-1)*nkf3 + k
                   ikk = 2 * ik - 1
                   ikq = ikk + 1
                   xkf(1, ikk) = dble(i-1)/dble(nkf1)
                   xkf(2, ikk) = dble(j-1)/dble(nkf2)
                   xkf(3, ikk) = dble(k-1)/dble(nkf3)
                   xkf(1, ikq) = xkf(1, ikk)
                   xkf(2, ikq) = xkf(2, ikk)
                   xkf(3, ikq) = xkf(3, ikk)
                ENDDO
             ENDDO
          ENDDO
          !
       ENDIF
    ELSEIF (rand_k) THEN  ! random points
      WRITE (stdout, *) '    Using random k-mesh: ', rand_nk
      !
      nkqtotf = rand_nk
      ALLOCATE (xkf(3, 2*nkqtotf), wkf(2*nkqtotf))
      !
      CALL init_random_seed()
      !
      DO ik = 1, nkqtotf
         !
         ikk = 2 * ik - 1
         ikq = ikk + 1
         !
         wkf(ikk) = 2.d0/ dble(nkqtotf)
         wkf(ikq) = 0.d0
         !
         IF ( system_2d ) THEN
            CALL random_number(xkf(1:2,ikk))
            xkf(3,ikk) = 0.d0
         ELSE
            CALL random_number(xkf(:,ikk))
         ENDIF
         !
         xkf(:,ikq) = xkf(:,ikk)
         !
      ENDDO
      !
      ! redefine nkqtotf to include the k+q points
      !
      nkqtotf = 2 * nkqtotf
      ! 
    ELSE ! don't know how to get grid
       CALL errore('loadkmesh_serial', "Cannot load fine k points", 1)
    ENDIF
    !
    ! Serial
    nkf = nkqtotf/2
    nkqf = nkqtotf
    !
  ENDIF
  !
  CALL mp_bcast (nkf, ionode_id, inter_pool_comm)
  CALL mp_bcast (nkqf, ionode_id, inter_pool_comm)
  CALL mp_bcast (nkqtotf, ionode_id, inter_pool_comm)
  IF (.not.ALLOCATED(xkf)) ALLOCATE (xkf(3,nkqtotf))
  IF (.not.ALLOCATED(wkf)) ALLOCATE (wkf(  nkqtotf))
  CALL mp_bcast(xkf, ionode_id, inter_pool_comm)
  CALL mp_bcast(wkf, ionode_id, inter_pool_comm)
  !
  IF (abs(sum (wkf) - 2.d0) .gt. 1.d-4 ) &
    WRITE(stdout,'(5x,"WARNING: k-point weigths do not add up to 1 [loadkmesh_serial]")') 
  !
  WRITE( stdout, '(5x,"Size of k point mesh for interpolation: ",i10)' ) nkqtotf
  !
END SUBROUTINE loadkmesh_serial
!
!-----------------------------------------------------------------------
!
SUBROUTINE init_random_seed()
  !
  INTEGER :: i, n, clock
  INTEGER, DIMENSION(:), ALLOCATABLE :: seed
  !     
  CALL RANDOM_SEED(size = n)
  ALLOCATE(seed(n))
  !      
  CALL SYSTEM_CLOCK(COUNT=clock)
  !        
  seed = clock + 37 * (/ (i - 1, i = 1, n) /)
  CALL RANDOM_SEED(PUT = seed)
  !        
  DEALLOCATE(seed)
  !
END SUBROUTINE init_random_seed
!-----------------------------------------------------------------------
