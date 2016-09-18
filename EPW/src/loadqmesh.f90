  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
SUBROUTINE loadqmesh_para
  !-----------------------------------------------------------------------
  !!
  !!  Load fine q mesh and distribute among pools
  !!
  !-----------------------------------------------------------------------
  USE io_global, ONLY : ionode_id
  USE mp_global, ONLY : inter_pool_comm, my_pool_id, npool
  USE mp,        ONLY : mp_bcast, mp_sum
  USE mp_world,  ONLY : mpime 
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout
  USE epwcom,    ONLY : filqf, nqf1, nqf2, nqf3, &
                        rand_q, rand_nq, mp_mesh_q, system_2d
  USE elph2,     ONLY : xqf, wqf, nqf, nqtotf
  USE cell_base, ONLY : at, bg
  USE symm_base, ONLY : s, t_rev, time_reversal, set_sym_bl, nrot
  USE io_epw,    ONLY : iunqf

  USE noncollin_module, ONLY : noncolin
  !
  implicit none
  !
  real(kind=DP), ALLOCATABLE :: xqf_(:,:), wqf_(:)
  integer :: iq, lower_bnd, upper_bnd, i, j, k, ios
  !
  !
  integer :: rest
  !
  IF (mpime .eq. ionode_id) THEN
    IF (filqf .ne. '') THEN ! load from file (crystal coordinates)
       !
       WRITE(stdout, *) '     Using q-mesh file: ', trim(filqf)
       OPEN( unit = iunqf, file = filqf, status = 'old', form = 'formatted',err=100, iostat=ios)
100    CALL errore('loadkmesh_para','opening file '//filqf,abs(ios))
       READ(iunqf, *) nqtotf
       !
       ALLOCATE (xqf_ (3, nqtotf), wqf_(nqtotf))
       !
       DO iq = 1, nqtotf
          !
          READ (iunqf, *) xqf_ (:, iq ), wqf_ (iq)
          !
       ENDDO
       CLOSE(iunqf)
       !
    ELSEIF ( (nqf1.ne.0) .and. (nqf2.ne.0) .and. (nqf3.ne.0) ) THEN ! generate grid
       IF (mp_mesh_q) THEN
           ! get size of the mp_mesh in the irr wedge 
          WRITE (stdout, '(a,3i4)') '     Using uniform MP q-mesh: ', nqf1, nqf2, nqf3
          call set_sym_bl ( )
          !
          ALLOCATE ( xqf_ (3, nqf1*nqf2*nqf3), wqf_(nqf1*nqf2*nqf3) )
          ! the result of this call is just nkqtotf
          CALL kpoint_grid ( nrot, time_reversal, .false., s, t_rev, bg, nqf1*nqf2*nqf3, &
               0,0,0, nqf1,nqf2,nqf3, nqtotf, xqf_, wqf_)
          DEALLOCATE (xqf_, wqf_)
          ALLOCATE ( xqf_ (3, nqtotf), wqf_(nqtotf)) 
          CALL kpoint_grid ( nrot, time_reversal, .false., s, t_rev, bg, nqf1*nqf2*nqf3, &
               0,0,0, nqf1,nqf2,nqf3, nqtotf, xqf_, wqf_)
          !  
          ! bring the k point to crystal coordinates       
          CALL cryst_to_cart (nqtotf, xqf_, at, -1)
          !
       ELSE
          !
          WRITE (stdout, '(a,3i4)') '     Using uniform q-mesh: ', nqf1, nqf2, nqf3
          !
          nqtotf =  nqf1 * nqf2 * nqf3
          ALLOCATE ( xqf_ (3, nqtotf), wqf_(nqtotf) )
          wqf_(:) = 0.d0
          DO iq = 1, nqf1 * nqf2 * nqf3
             wqf_(iq) = 1.d0/(dble(nqtotf))
          ENDDO
          DO i = 1, nqf1
             DO j = 1, nqf2
                DO k = 1, nqf3
                   iq = (i-1)*nqf2*nqf3 + (j-1)*nqf3 + k
                   xqf_(1, iq) = dble(i-1)/dble(nqf1)
                   xqf_(2, iq) = dble(j-1)/dble(nqf2)
                   xqf_(3, iq) = dble(k-1)/dble(nqf3)
                ENDDO
             ENDDO
          ENDDO
          !
       ENDIF
       !
    ELSEIF (rand_q) THEN  ! random points
       !
       WRITE (stdout, *) '    Using random q-mesh: ', rand_nq
       !
       nqtotf = rand_nq
       ALLOCATE (xqf_ (3, nqtotf), wqf_(nqtotf))
       !
       CALL init_random_seed()
       !
       DO iq = 1, nqtotf
          !
          !
          wqf_(iq) = 1.d0/ dble(nqtotf)
          !
          IF ( system_2d ) THEN
             CALL random_number(xqf_(1:2,iq))
             xqf_(3,iq) = 0.d0
          ELSE
             CALL random_number(xqf_(:,iq))
          ENDIF
          !
       ENDDO
       !
    ELSE ! don't know how to get grid
       CALL errore('loadqmesh_para', "Cannot load fine q points", 1)
    ENDIF
 ENDIF
 !
#if defined(__MPI)
 CALL mp_bcast (nqtotf, ionode_id, inter_pool_comm)
 !
 !  scatter the q points of the fine mesh across the pools
 !
 nqf = ( nqtotf / npool )
 rest = ( nqtotf - nqf * npool ) / 2
 IF (my_pool_id < rest ) THEN
    nqf = nqf + 2
    lower_bnd = my_pool_id*nqf + 1
    upper_bnd = lower_bnd + nqf - 1
 ELSE
    lower_bnd = rest*(nqf+2)+(my_pool_id-rest)*nqf + 1
    upper_bnd = lower_bnd + nqf - 1
 ENDIF
 !
 IF (.not.ALLOCATED(xqf_)) ALLOCATE (xqf_(3,nqtotf))
 IF (.not.ALLOCATED(wqf_)) ALLOCATE (wqf_(  nqtotf))
 CALL mp_bcast(xqf_, ionode_id, inter_pool_comm)
 CALL mp_bcast(wqf_, ionode_id, inter_pool_comm)
 !
 !CALL FLUSH()
 ! 
#else
 !
 ! In serial the definitions are much easier 
 !
 nqf = nqtotf  
 lower_bnd = 1
 upper_bnd = nqf
 !
#endif
 !
 !  Assign the weights and vectors to the correct bounds
 !
 ALLOCATE(xqf(3,nqf))
 ALLOCATE(wqf(  nqf))
 xqf(:,:) = xqf_ (:, lower_bnd:upper_bnd)
 IF (noncolin) THEN 
    wqf(  :) = wqf_ ( lower_bnd:upper_bnd)/2.d0
 ELSE
    wqf(  :) = wqf_ ( lower_bnd:upper_bnd)
 ENDIF  
 !
 IF (abs(sum (wqf_ (:)) - 1.d0) .gt. 1.d-4 ) &
    WRITE(stdout,'(5x,"WARNING: q-point weigths do not add up to 1 [loadqmesh_para]")')
 !
 WRITE( stdout, '(5x,"Size of q point mesh for interpolation: ",i10)' ) nqtotf 
 WRITE( stdout, '(5x,"Max number of q points per pool:",7x,i10)' ) nqf 
 !
 IF (ALLOCATED(xqf_)) DEALLOCATE(xqf_)
 IF (ALLOCATED(wqf_)) DEALLOCATE(wqf_)
 !
END SUBROUTINE loadqmesh_para
!-----------------------------------------------------------------------
!

!-----------------------------------------------------------------------
SUBROUTINE loadqmesh_serial
!-----------------------------------------------------------------------
!!
!!  Load fine q mesh
!!
!-----------------------------------------------------------------------
  USE io_global, ONLY : ionode_id
  USE mp_global, ONLY : inter_pool_comm
  USE mp,        ONLY : mp_bcast
  USE mp_world,  ONLY : mpime
  USE io_global, ONLY : stdout
  USE epwcom,    ONLY : filqf, nqf1, nqf2, nqf3, &
                        rand_q, rand_nq, mp_mesh_q, system_2d
  USE elph2,     ONLY : xqf, wqf, nqtotf, nqf
  USE pwcom,     ONLY : at, bg
  USE symm_base, ONLY : s, t_rev, time_reversal, set_sym_bl, nrot
  USE io_epw,    ONLY : iunqf
  implicit none
  !
  integer :: iq, i, j, k, ios
  !
  IF (mpime .eq. ionode_id) THEN
    IF (filqf .ne. '') THEN ! load from file (crystal coordinates)
       !
       ! Each pool gets its own copy from the action=read statement
       !
       WRITE (stdout, *) '     Using q-mesh file: ', trim(filqf)
       OPEN ( unit = iunqf, file = filqf, status = 'old', form = 'formatted', err=100, iostat=ios)
100    CALL errore('loadqmesh_serial','opening file '//filqf,abs(ios))
       READ(iunqf, *) nqtotf
       ALLOCATE (xqf(3, nqtotf), wqf(nqtotf))
       DO iq = 1, nqtotf
          READ (iunqf, *) xqf (:, iq), wqf(iq)
       ENDDO
       CLOSE(iunqf)
       !
       ! bring xqf in crystal coordinates
       ! CALL cryst_to_cart (nqtotf, xqf, at, -1)
       !
    ELSEIF ( (nqf1.ne.0) .and. (nqf2.ne.0) .and. (nqf3.ne.0) ) THEN ! generate grid
       IF (mp_mesh_q) THEN
          ! get size of the mp_mesh in the irr wedge 
          WRITE (stdout, '(a,3i4)') '     Using uniform q-mesh: ', nqf1, nqf2, nqf3
          call set_sym_bl ( )
          !                                         
          ALLOCATE ( xqf (3, nqf1*nqf2*nqf3), wqf(nqf1*nqf2*nqf3) )
          ! the result of this call is just nkqtotf
          CALL kpoint_grid ( nrot, time_reversal, s, t_rev, bg, nqf1*nqf2*nqf3, &
               0,0,0, nqf1,nqf2,nqf3, nqtotf, xqf, wqf)
          DEALLOCATE ( xqf, wqf) 
          ALLOCATE ( xqf(3, nqtotf), wqf(nqtotf)) 
          CALL kpoint_grid ( nrot, time_reversal, s, t_rev, bg, nqf1*nqf2*nqf3, &
               0,0,0, nqf1,nqf2,nqf3, nqtotf, xqf, wqf)
          !
          ! bring xqf in crystal coordinates       
          CALL cryst_to_cart (nqtotf, xqf, at, -1)
          !
       ELSE
          ! currently no offset.  
          ! q's are in crystal coordinates in xqf
          WRITE (stdout, '(a,3i4)') '     Using uniform q-mesh: ', nqf1, nqf2, nqf3
          !
          nqtotf = nqf1 * nqf2 * nqf3
          ALLOCATE ( xqf (3, nqtotf), wqf(nqtotf) )
          wqf(:) = 1.d0 / (dble(nqtotf))
          DO i = 1, nqf1
             DO j = 1, nqf2
                DO k = 1, nqf3
                   iq = (i-1)*nqf2*nqf3 + (j-1)*nqf3 + k
                   xqf(1, iq) = dble(i-1)/dble(nqf1)
                   xqf(2, iq) = dble(j-1)/dble(nqf2)
                   xqf(3, iq) = dble(k-1)/dble(nqf3)
                ENDDO
             ENDDO
          ENDDO
          !
          !WRITE(stdout,'(a)') '  '
          !DO iq = 1, nqtotf
             !WRITE(stdout,'(4f12.7)') xqf(:,iq), wqf(iq)
          !ENDDO
          !WRITE(stdout,'(a)') '  '
          !
       ENDIF
    ELSEIF (rand_q) THEN  ! random points
       ! random grid
       WRITE (stdout, *) '    Using random q-mesh: ', rand_nq
       !
       nqtotf = rand_nq
       ALLOCATE (xqf(3, nqtotf), wqf(nqtotf))
       !WRITE(stdout,'(a)') '  '
       wqf(:) = 1.d0/(dble(nqtotf))
       !
       CALL init_random_seed()
       !
       DO iq = 1, nqtotf
          !
          IF ( system_2d ) THEN
             CALL random_number(xqf(1:2,iq))
             xqf(3,iq) = 0.d0
          ELSE
             CALL random_number(xqf(:,iq))
          ENDIF
          !
          !WRITE(stdout,'(4f12.7)') xqf(:,iq), wqf(iq)
          !
       ENDDO
       !WRITE(stdout,'(a)') '  '
       !
    ELSE ! don't know how to get grid
       CALL errore('loadqmesh_serial', "Cannot load fine q points", 1)
    ENDIF
    !
    ! Since serial 
    nqf = nqtotf 
  ENDIF
  !
  CALL mp_bcast (nqf, ionode_id, inter_pool_comm)
  CALL mp_bcast (nqtotf, ionode_id, inter_pool_comm)
  IF (.not.ALLOCATED(xqf)) ALLOCATE (xqf(3,nqtotf))
  IF (.not.ALLOCATED(wqf)) ALLOCATE (wqf(  nqtotf))
  CALL mp_bcast(xqf, ionode_id, inter_pool_comm)
  CALL mp_bcast(wqf, ionode_id, inter_pool_comm)
  !
  IF (abs(sum (wqf) - 1.d0) .gt. 1.d-4 ) &
    WRITE(stdout,'(5x,"WARNING: q-point weigths do not add up to 1 [loadqmesh_serial]")') 
  !
  WRITE( stdout, '(5x,"Size of q point mesh for interpolation: ",i10)' ) nqtotf
  !
END SUBROUTINE loadqmesh_serial
!
!-----------------------------------------------------------------------
