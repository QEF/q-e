subroutine read_export (pp_file,kunit,uspp_spsi, ascii, single_file, raw)
  !-----------------------------------------------------------------------
  !
  use iotk_module


  use kinds,          ONLY : DP 
  use pwcom  
  use cell_base,      ONLY : tpiba2, bg
  use gvect,          ONLY : ngm, ngm_g, ig_l2g, mill, g 
  use control_flags,  ONLY : gamma_only  
  use becmod,         ONLY : bec_type, becp, calbec, &
                             allocate_bec_type, deallocate_bec_type
!  use symme,          ONLY : nsym, s, invsym, sname, irt, ftau
!  use symme,          ONLY : nsym, s, invsym, irt, ftau
!  use char,           ONLY : sname
! occhio sname is in symme which is now outside pwcom
  use  uspp,          ONLY : nkb, vkb
  use wavefunctions_module,  ONLY : evc
  use io_files,       ONLY : nd_nmbr, prefix, iunwfc, nwordwfc, iunsat, nwordatwfc
  use io_files,       ONLY : pseudo_dir, psfile
  use io_global,      ONLY : ionode, stdout
  USE ions_base,      ONLY : atm, nat, ityp, tau, nsp
  use mp_world,      ONLY : nproc,  mpime
  use mp_pools,      ONLY : my_pool_id, intra_pool_comm, inter_pool_comm, nproc_pool
  USE mp_world,             ONLY : world_comm
  use mp,             ONLY : mp_sum, mp_max
!  use ldaU,           ONLY : swfcatom, lda_plus_u
  use ldaU,           ONLY :  lda_plus_u
  USE gvecw,              ONLY :  ecutwfc
  USE klist, ONLY : igk_k

  implicit none

  CHARACTER(5), PARAMETER :: fmt_name="QEXPT"
  CHARACTER(5), PARAMETER :: fmt_version="1.1.0"

  integer, intent(in) :: kunit
  character(80), intent(in) :: pp_file
  logical, intent(in) :: uspp_spsi, ascii, single_file, raw

  integer :: i, j, k, ig, ik, ibnd, na, ngg,ig_, ierr
  integer, allocatable :: kisort(:)
  real(DP) :: xyz(3), tmp(3)
  integer :: ike, iks, npw_g, npwx_g, ispin, local_pw
  integer, allocatable :: ngk_g( : )
  integer, allocatable :: itmp_g( :, : )
  real(DP),allocatable :: rtmp_g( :, : )
  real(DP),allocatable :: rtmp_gg( : )
  integer, allocatable :: itmp1( : )
  integer, allocatable :: igwk( :, : )
  integer, allocatable :: l2g_new( : )
  integer, allocatable :: igk_l2g( :, : )
  integer, external :: global_kpoint_index

  real(DP) :: wfc_scal 
  logical :: twf0, twfm
  character(iotk_attlenx) :: attr
  complex(DP), allocatable :: sevc (:,:)

  call start_clock('read_export')
  write(stdout,*) "nkstot=", nkstot

  IF( nkstot > 0 ) THEN

     IF( ( kunit < 1 ) .OR. ( MOD( nkstot, kunit ) /= 0 ) ) &
       CALL errore( ' write_export ',' wrong kunit ', 1 )

     IF( ( nproc_pool > nproc ) .OR. ( MOD( nproc, nproc_pool ) /= 0 ) ) &
       CALL errore( ' write_export ',' nproc_pool ', 1 )

     iks = global_kpoint_index (nkstot, 1)
     ike = iks + nks - 1

  END IF

  write(stdout,*) "after first init"

  ! find out the global number of G vectors: ngm_g  
  ngm_g = ngm
  call mp_sum( ngm_g , world_comm )

  ! collect all G vectors across processors within the pools
  ! and compute their modules
  !
  allocate( itmp_g( 3, ngm_g ) )
  allocate( rtmp_g( 3, ngm_g ) )
  allocate( rtmp_gg( ngm_g ) )

  itmp_g = 0
  do  ig = 1, ngm
    itmp_g( 1, ig_l2g( ig ) ) = mill(1, ig )
    itmp_g( 2, ig_l2g( ig ) ) = mill(2, ig )
    itmp_g( 3, ig_l2g( ig ) ) = mill(3, ig )
  end do
  call mp_sum( itmp_g , world_comm )
  !
  ! here we are in crystal units
  rtmp_g(1:3,1:ngm_g) = REAL( itmp_g(1:3,1:ngm_g) )
  !
  ! go to cartesian units (tpiba)
  call cryst_to_cart( ngm_g, rtmp_g, bg , 1 )
  !
  ! compute squared moduli
  do  ig = 1, ngm_g 
     rtmp_gg(ig) = rtmp_g(1,ig)**2 + rtmp_g(2,ig)**2 + rtmp_g(3,ig)**2 
  enddo
  deallocate( rtmp_g )

  ! build the G+k array indexes
  allocate ( igk_l2g ( npwx, nks ) )
  allocate ( kisort( npwx ) )
  do ik = 1, nks
     kisort = 0
     npw = npwx
     call gk_sort (xk (1, ik+iks-1), ngm, g, ecutwfc / tpiba2, npw, kisort(1), g2kin)
     !
     ! mapping between local and global G vector index, for this kpoint
     !
     DO ig = 1, npw
        !
        igk_l2g(ig,ik) = ig_l2g( kisort(ig) )
        !
     END DO
     !
     igk_l2g( npw+1 : npwx, ik ) = 0
     !
     ngk (ik) = npw
  end do
  deallocate (kisort)

  ! compute the global number of G+k vectors for each k point
  allocate( ngk_g( nkstot ) )
  ngk_g = 0
  ngk_g( iks:ike ) = ngk( 1:nks )
  CALL mp_sum( ngk_g, world_comm )

  ! compute the Maximum G vector index among all G+k and processors
  npw_g = MAXVAL( igk_l2g(:,:) )
  CALL mp_max( npw_g, world_comm )

  ! compute the Maximum number of G vector among all k points
  npwx_g = MAXVAL( ngk_g( 1:nkstot ) )

  deallocate(rtmp_gg)

  allocate( igwk( npwx_g,nkstot ) )

  write(stdout,*) "after g stuff"

! wfc grids

  DO ik = 1, nkstot
    igwk(:,ik) = 0
    !
    ALLOCATE( itmp1( npw_g ), STAT= ierr )
    IF ( ierr/=0 ) CALL errore('pw_export','allocating itmp1', ABS(ierr) )
    itmp1 = 0
    ! 
    IF( ik >= iks .AND. ik <= ike ) THEN 
      DO  ig = 1, ngk( ik-iks+1 )
        itmp1( igk_l2g( ig, ik-iks+1 ) ) = igk_l2g( ig, ik-iks+1 ) 
      END DO
    END IF
    !
    CALL mp_sum( itmp1, world_comm )
    !
    ngg = 0
    DO  ig = 1, npw_g
      IF( itmp1( ig ) == ig ) THEN
        ngg = ngg + 1
        igwk( ngg , ik) = ig
      END IF
    END DO
    IF( ngg /= ngk_g( ik ) ) THEN
      WRITE( stdout,*) ' ik, ngg, ngk_g = ', ik, ngg, ngk_g( ik )
    END IF
    !
    DEALLOCATE( itmp1 )
    !
  ENDDO
  !
  deallocate( itmp_g )
  
  write(stdout,*)"after wfc waves"

  call poolrecover (et, nbnd, nkstot, nks)
 
  wfc_scal = 1.0d0
  twf0 = .true.
  twfm = .false.

  do ik = 1, nkstot
     local_pw = 0
     IF( (ik >= iks) .AND. (ik <= ike) ) THEN

       call davcio (evc, 2*nwordwfc, iunwfc, (ik-iks+1), - 1)
!       IF ( lda_plus_u ) CALL davcio( swfcatom, nwordatwfc, iunsat, (ik-iks+1), -1 )
       local_pw = ngk(ik-iks+1)

     ENDIF


     allocate(l2g_new(local_pw))

     l2g_new = 0
     do ig = 1, local_pw
       ngg = igk_l2g(ig,ik-iks+1)
       do ig_ = 1, ngk_g(ik)
         if(ngg == igwk(ig_,ik)) then
           l2g_new(ig) = ig_
           exit
         end if
       end do
     end do


     ispin = isk( ik )
     !  WRITE(0,*) ' ### ', ik,nkstot,iks,ike,kunit,nproc,nproc_pool 
     deallocate(l2g_new)
  end do
  !  

  write(stdout,*) "after davcio"

  ! If specified and if USPP are used the wfcs S_psi are written  
  ! | spsi_nk > = \hat S | psi_nk >  
  ! where S is the overlap operator of US PP 
  !  
  IF ( uspp_spsi .AND. nkb > 0 ) THEN

       ALLOCATE( sevc(npwx,nbnd), STAT=ierr )
       IF (ierr/=0) CALL errore( ' read_export ',' Unable to allocate SEVC ', ABS(ierr) )

       CALL init_us_1
       CALL init_at_1

       CALL allocate_bec_type (nkb,nbnd,becp)

       do ik = 1, nkstot
 
           local_pw = 0
           IF( (ik >= iks) .AND. (ik <= ike) ) THEN
               
               CALL gk_sort (xk (1, ik+iks-1), ngm, g, ecutwfc / tpiba2, npw, igk_k(1,1), g2kin)
               CALL davcio (evc, 2*nwordwfc, iunwfc, (ik-iks+1), - 1)

               CALL init_us_2(npw, igk_k(1,1), xk(1, ik), vkb)
               local_pw = ngk(ik-iks+1)
                            
               IF ( gamma_only ) THEN
                  if(nkb>0) CALL calbec ( ngk_g(ik), vkb, evc, becp )
               ELSE
                  CALL calbec ( npw, vkb, evc, becp )
               ENDIF
               CALL s_psi(npwx, npw, nbnd, evc, sevc)
           ENDIF

           ALLOCATE(l2g_new(local_pw))

           l2g_new = 0
           DO ig = 1, local_pw
             ngg = igk_l2g(ig,ik-iks+1)
             DO ig_ = 1, ngk_g(ik)
               IF(ngg == igwk(ig_,ik)) THEN
                 l2g_new(ig) = ig_
                 EXIT
               ENDIF
             ENDDO
           ENDDO

           ispin = isk( ik )
           DEALLOCATE(l2g_new)
       ENDDO
      
       DEALLOCATE( sevc, STAT=ierr )
       IF ( ierr/= 0 ) CALL errore('read_export','Unable to deallocate SEVC',ABS(ierr))
       CALL deallocate_bec_type ( becp )
  ENDIF

  DEALLOCATE( igk_l2g )
  DEALLOCATE( igwk )
  DEALLOCATE ( ngk_g )
  call stop_clock('read_export')

end subroutine read_export
