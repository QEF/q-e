! Copyright (C) 2005-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------
MODULE exx_band
  !--------------------------------------
  !
  ! Variables and subroutines for band parallelization over pairs of bands:
  ! see T. Barnes, T. Kurth, P. Carrier, N. Wichmann, D. Prendergast,
  ! P.R.C. Kent, J. Deslippe, Computer Physics Communications 2017,
  ! dx.doi.org/10.1016/j.cpc.2017.01.008
  !
  USE kinds,                ONLY : DP
  USE noncollin_module,     ONLY : npol
  USE io_global,            ONLY : ionode, stdout
  !
  USE control_flags,        ONLY : gamma_only
  USE fft_types,            ONLY : fft_type_descriptor
  USE stick_base,           ONLY : sticks_map
  !
  IMPLICIT NONE
  SAVE
  !!
  COMPLEX(DP), ALLOCATABLE :: evc_exx(:,:)
  COMPLEX(DP), ALLOCATABLE :: psi_exx(:,:), hpsi_exx(:,:)
  INTEGER :: lda_original, n_original
  INTEGER :: nwordwfc_exx
  INTEGER, ALLOCATABLE :: igk_exx(:,:)
  !
  ! mapping for the data structure conversion
  !
  TYPE comm_packet
     INTEGER :: size
     INTEGER, ALLOCATABLE :: indices(:)
     COMPLEX(DP), ALLOCATABLE :: msg(:,:,:)
  END TYPE comm_packet
  TYPE(comm_packet), ALLOCATABLE :: comm_recv(:,:), comm_send(:,:)
  TYPE(comm_packet), ALLOCATABLE :: comm_recv_reverse(:,:)
  TYPE(comm_packet), ALLOCATABLE :: comm_send_reverse(:,:,:)
  INTEGER, ALLOCATABLE :: lda_local(:,:)
  INTEGER, ALLOCATABLE :: lda_exx(:,:)
  INTEGER, ALLOCATABLE :: ngk_local(:), ngk_exx(:)
  INTEGER :: npwx_local = 0
  INTEGER :: npwx_exx = 0
  INTEGER :: n_local = 0
  LOGICAL :: first_data_structure_change = .TRUE.

  INTEGER :: ngm_loc, ngm_g_loc, gstart_loc
  INTEGER, ALLOCATABLE :: ig_l2g_loc(:)
  REAL(DP), ALLOCATABLE :: g_loc(:,:), gg_loc(:)
  INTEGER, ALLOCATABLE :: mill_loc(:,:), nl_loc(:)
  INTEGER :: ngms_loc, ngms_g_loc
  INTEGER, ALLOCATABLE :: nls_loc(:)
  INTEGER, ALLOCATABLE :: nlm_loc(:)
  INTEGER, ALLOCATABLE :: nlsm_loc(:)

  INTEGER :: ngm_exx, ngm_g_exx, gstart_exx
  INTEGER, ALLOCATABLE :: ig_l2g_exx(:)
  REAL(DP), ALLOCATABLE :: g_exx(:,:), gg_exx(:)
  INTEGER, ALLOCATABLE :: mill_exx(:,:), nl_exx(:)
  INTEGER :: ngms_exx, ngms_g_exx
  INTEGER, ALLOCATABLE :: nls_exx(:)
  INTEGER, ALLOCATABLE :: nlm_exx(:)
  INTEGER, ALLOCATABLE :: nlsm_exx(:)

  TYPE(fft_type_descriptor) :: dfftp_loc, dffts_loc
  TYPE(fft_type_descriptor) :: dfftp_exx, dffts_exx
  TYPE (sticks_map) :: smap_exx ! Stick map descriptor

 CONTAINS
  !
  !-----------------------------------------------------------------------
  SUBROUTINE transform_evc_to_exx(type)
  !-----------------------------------------------------------------------
    !
    USE wvfct,                ONLY : npwx, nbnd
    USE io_files,             ONLY : nwordwfc, iunwfc, iunwfc_exx
    USE klist,                ONLY : nks, ngk, igk_k
    USE uspp,                 ONLY : nkb, vkb
    USE wavefunctions, ONLY : evc
    USE control_flags,        ONLY : io_level
    USE buffers,              ONLY : open_buffer, get_buffer, save_buffer
    USE mp_exx,               ONLY : max_ibands, negrp
    !
    !
    IMPLICIT NONE
    !
    INTEGER, intent(in) :: type
    INTEGER :: lda, n, ik
    LOGICAL :: exst_mem, exst_file
    !
    IF (negrp.eq.1) THEN
       !
       ! no change in data structure is necessary
       ! just copy all of the required data
       !
       ! get evc_exx
       !
       IF(.not.allocated(evc_exx))THEN
          ALLOCATE(evc_exx(npwx*npol,nbnd))
       END IF
       evc_exx = evc
       !
       ! get igk_exx
       !
       IF(.not.allocated(igk_exx)) THEN
          ALLOCATE( igk_exx( npwx, nks ) )
          igk_exx = igk_k
       END IF
       !
       ! get the wfc buffer is used
       !
       iunwfc_exx = iunwfc
       nwordwfc_exx = nwordwfc
       !
       RETURN
       !
    END IF
    !
    ! change the data structure of evc and igk
    !
    lda = npwx
    n = npwx 
    npwx_local = npwx
    IF( .not.allocated(ngk_local) ) allocate(ngk_local(nks))
    ngk_local = ngk
    !
    IF ( .not.allocated(comm_recv) ) THEN
       !
       ! initialize all of the conversion maps and change the data structure
       !
       CALL initialize_local_to_exact_map(lda, nbnd)
    ELSE
       !
       ! change the data structure
       !
       CALL change_data_structure(.TRUE.)
    END IF
    !
    lda = npwx
    n = npwx 
    npwx_exx = npwx
    IF( .not.allocated(ngk_exx) ) allocate(ngk_exx(nks))
    ngk_exx = ngk
    !
    ! get evc_exx
    !
    IF(.not.allocated(evc_exx))THEN
       ALLOCATE(evc_exx(lda*npol,max_ibands+2))
       !
       ! ... open files/buffer for wavefunctions (nwordwfc set in openfil)
       ! ... io_level > 1 : open file, otherwise: open buffer
       !
       nwordwfc_exx  = size(evc_exx)
       CALL open_buffer( iunwfc_exx, 'wfc_exx', nwordwfc_exx, io_level, &
            exst_mem, exst_file )
    END IF
    !
    DO ik=1, nks
       !
       ! read evc for the local data structure
       !
       IF ( nks > 1 ) CALL get_buffer(evc, nwordwfc, iunwfc, ik)
       !
       ! transform evc to the EXX data structure
       !
       CALL transform_to_exx(lda, n, nbnd, nbnd, ik, evc, evc_exx, type)
       !
       ! save evc to a buffer
       !
       IF ( nks > 1 ) CALL save_buffer ( evc_exx, nwordwfc_exx, iunwfc_exx, ik )
    END DO
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE transform_evc_to_exx
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE transform_psi_to_exx(lda, n, m, psi)
  !-----------------------------------------------------------------------
    USE wvfct,        ONLY : current_k, npwx, nbnd
    USE mp_exx,       ONLY : negrp, nibands, my_egrp_id, max_ibands
    !
    !
    IMPLICIT NONE
    !
    Integer, INTENT(in) :: lda
    INTEGER, INTENT(in) :: m
    INTEGER, INTENT(inout) :: n
    COMPLEX(DP), INTENT(in) :: psi(lda*npol,m) 
    !
    ! change to the EXX data strucutre
    !
    npwx_local = npwx
    n_local = n
    IF ( .not.allocated(comm_recv) ) THEN
       !
       ! initialize all of the conversion maps and change the data structure
       !
       CALL initialize_local_to_exact_map(lda, m)
    ELSE
       !
       ! change the data structure
       !
       CALL change_data_structure(.TRUE.)
    END IF
    npwx_exx = npwx
    n = ngk_exx(current_k)
    !
    ! get igk
    !
    CALL update_igk(.TRUE.)
    !
    ! get psi_exx
    !
    CALL transform_to_exx(lda, n, m, max_ibands, &
         current_k, psi, psi_exx, 0)
    !
    ! zero hpsi_exx
    !
    hpsi_exx = 0.d0
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE transform_psi_to_exx
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE transform_hpsi_to_local(lda, n, m, hpsi)
  !-----------------------------------------------------------------------
    USE mp_exx,         ONLY : iexx_istart, iexx_iend, my_egrp_id
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: lda
    INTEGER, INTENT(in) :: m
    INTEGER, INTENT(inout) :: n
    COMPLEX(DP), INTENT(out) :: hpsi(lda_original*npol,m)
    INTEGER :: m_exx
    !
    ! change to the local data structure
    !
    CALL change_data_structure(.FALSE.)
    !
    ! get igk
    !
    CALL update_igk(.FALSE.)
    n = n_local
    !
    ! transform hpsi_exx to the local data structure
    !
    m_exx = iexx_iend(my_egrp_id+1) - iexx_istart(my_egrp_id+1) + 1
    CALL transform_to_local(m,m_exx,hpsi_exx,hpsi)
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE transform_hpsi_to_local
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE initialize_local_to_exact_map(lda, m)
  !-----------------------------------------------------------------------
    !
    ! determine the mapping between the local and EXX data structures
    !
    USE wvfct,          ONLY : npwx, nbnd
    USE klist,          ONLY : nks, igk_k
    USE mp_exx,         ONLY : nproc_egrp, negrp, my_egrp_id, me_egrp, &
                               intra_egrp_comm, inter_egrp_comm, &
                               ibands, nibands, init_index_over_band, &
                               iexx_istart, iexx_iend, max_ibands
    USE mp_pools,       ONLY : nproc_pool, me_pool, intra_pool_comm
    USE mp,             ONLY : mp_sum
    USE gvect,          ONLY : ig_l2g
    USE uspp,           ONLY : nkb
    !
    !
    IMPLICIT NONE
    !
    Integer :: lda
    INTEGER :: n, m
    INTEGER, ALLOCATABLE :: local_map(:,:), exx_map(:,:)
    INTEGER, ALLOCATABLE :: l2e_map(:,:), e2l_map(:,:,:)
    INTEGER, ALLOCATABLE :: psi_source(:), psi_source_exx(:,:)
    INTEGER :: current_index
    INTEGER :: i, j, k, ik, im, ig, count, iproc, prev, iegrp
    INTEGER :: total_lda(nks), prev_lda(nks)
    INTEGER :: total_lda_exx(nks), prev_lda_exx(nks)
    INTEGER :: n_local
    INTEGER :: lda_max_local, lda_max_exx
    INTEGER :: max_lda_egrp
    INTEGER :: request_send(nproc_egrp), request_recv(nproc_egrp)
    INTEGER :: ierr
    INTEGER :: egrp_base, total_lda_egrp(nks), prev_lda_egrp(nks)
    INTEGER :: igk_loc(npwx)
    !
    ! initialize the pair assignments
    !
    CALL init_index_over_band(inter_egrp_comm,nbnd,m)
    !
    ! allocate bookeeping arrays
    !
    IF ( .not.allocated(comm_recv) ) THEN
       ALLOCATE(comm_recv(nproc_egrp,nks),comm_send(nproc_egrp,nks))
    END IF
    IF ( .not.allocated(lda_local) ) THEN
       ALLOCATE(lda_local(nproc_pool,nks))
       ALLOCATE(lda_exx(nproc_egrp,nks))
    END IF
    !
    ! store the original values of lda and n
    !
    lda_original = lda
    n_original = n
    !
    ! construct the local map
    !
    lda_local = 0
    DO ik = 1, nks
       igk_loc = igk_k(:,ik)
       n = 0
       DO i = 1, size(igk_loc)
          IF(igk_loc(i).gt.0)n = n + 1
       END DO
       lda_local(me_pool+1,ik) = n
       CALL mp_sum(lda_local(:,ik),intra_pool_comm)
       total_lda(ik) = sum(lda_local(:,ik))
       prev_lda(ik) = sum(lda_local(1:me_pool,ik))
    END DO
    ALLOCATE(local_map(maxval(total_lda),nks))
    local_map = 0
    DO ik = 1, nks
       local_map(prev_lda(ik)+1:prev_lda(ik)+lda_local(me_pool+1,ik),ik) = &
            ig_l2g(igk_k(1:lda_local(me_pool+1,ik),ik))
    END DO
    CALL mp_sum(local_map,intra_pool_comm)
    !
    !-----------------------------------------!
    ! Switch to the exx data structure        !
    !-----------------------------------------!
    !
    CALL change_data_structure(.TRUE.)
    !
    ! construct the exx map
    !
    lda_exx = 0
    DO ik = 1, nks
       n = 0
       DO i = 1, size(igk_exx(:,ik))
          IF(igk_exx(i,ik).gt.0)n = n + 1
       END DO
       lda_exx(me_egrp+1,ik) = n
       CALL mp_sum(lda_exx(:,ik),intra_egrp_comm)
       total_lda_exx(ik) = sum(lda_exx(:,ik))
       prev_lda_exx(ik) = sum(lda_exx(1:me_egrp,ik))
    END DO
    ALLOCATE(exx_map(maxval(total_lda_exx),nks))
    exx_map = 0
    DO ik = 1, nks
       exx_map(prev_lda_exx(ik)+1:prev_lda_exx(ik)+lda_exx(me_egrp+1,ik),ik) = &
            ig_l2g(igk_exx(1:lda_exx(me_egrp+1,ik),ik))    
    END DO
    CALL mp_sum(exx_map,intra_egrp_comm)
    !
    ! construct the l2e_map
    !
    allocate( l2e_map(maxval(total_lda_exx),nks) )
    l2e_map = 0
    DO ik = 1, nks
       DO ig = 1, lda_exx(me_egrp+1,ik)
          DO j=1, total_lda(ik)
             IF( local_map(j,ik).EQ.exx_map(ig+prev_lda_exx(ik),ik) ) exit
          END DO
          l2e_map(ig+prev_lda_exx(ik),ik) = j
       END DO
    END DO
    CALL mp_sum(l2e_map,intra_egrp_comm)
    !
    ! plan communication for the data structure change
    !
    lda_max_local = maxval(lda_local)
    lda_max_exx = maxval(lda_exx)
    allocate(psi_source(maxval(total_lda_exx)))
    !
    DO ik = 1, nks
       !
       ! determine which task each value will come from
       !
       psi_source = 0
       DO ig = 1, lda_exx(me_egrp+1,ik)
          j = 1
          DO i = 1, nproc_pool
             j = j + lda_local(i,ik)
             IF( j.gt.l2e_map(ig+prev_lda_exx(ik),ik) ) exit
          END DO
          psi_source(ig+prev_lda_exx(ik)) = i-1
       END DO
       CALL mp_sum(psi_source,intra_egrp_comm)
       !
       ! allocate communication packets to receive psi and hpsi
       !
       DO iproc=0, nproc_egrp-1
          !
          ! determine how many values need to come from iproc
          !
          count = 0
          DO ig=1, lda_exx(me_egrp+1,ik)
             IF ( MODULO(psi_source(ig+prev_lda_exx(ik)),nproc_egrp)&
                  .eq.iproc ) THEN
                count = count + 1
             END IF
          END DO
          !
          ! allocate the communication packet
          !
          comm_recv(iproc+1,ik)%size = count
          IF (count.gt.0) THEN
             IF (.not.ALLOCATED(comm_recv(iproc+1,ik)%msg)) THEN
                ALLOCATE(comm_recv(iproc+1,ik)%indices(count))
                ALLOCATE(comm_recv(iproc+1,ik)%msg(count,npol,max_ibands+2))
             END IF
          END IF
          !
          ! determine which values need to come from iproc
          !
          count = 0
          DO ig=1, lda_exx(me_egrp+1,ik)
             IF ( MODULO(psi_source(ig+prev_lda_exx(ik)),nproc_egrp)&
                  .eq.iproc ) THEN
                count = count + 1
                comm_recv(iproc+1,ik)%indices(count) = ig
             END IF
          END DO
          !
       END DO
       !
       ! allocate communication packets to send psi and hpsi
       !
       prev = 0
       DO iproc=0, nproc_egrp-1
          !
          ! determine how many values need to be sent to iproc
          !
          count = 0
          DO ig=1, lda_exx(iproc+1,ik)
             IF ( MODULO(psi_source(ig+prev),nproc_egrp).eq.me_egrp ) THEN
                count = count + 1
             END IF
          END DO
          !
          ! allocate the communication packet
          !
          comm_send(iproc+1,ik)%size = count
          IF (count.gt.0) THEN
             IF (.not.ALLOCATED(comm_send(iproc+1,ik)%msg)) THEN
                ALLOCATE(comm_send(iproc+1,ik)%indices(count))
                ALLOCATE(comm_send(iproc+1,ik)%msg(count,npol,max_ibands+2))
             END IF
          END IF
          !
          ! determine which values need to be sent to iproc
          !
          count = 0
          DO ig=1, lda_exx(iproc+1,ik)
             IF ( MODULO(psi_source(ig+prev),nproc_egrp).eq.me_egrp ) THEN
                count = count + 1
                comm_send(iproc+1,ik)%indices(count) = l2e_map(ig+prev,ik)
             END IF
          END DO
          !
          prev = prev + lda_exx(iproc+1,ik)
          !
       END DO
       !
    END DO
    !
    ! allocate psi_exx and hpsi_exx
    !
    IF(allocated(psi_exx))DEALLOCATE(psi_exx)
    ALLOCATE(psi_exx(npwx*npol, max_ibands ))
    IF(allocated(hpsi_exx))DEALLOCATE(hpsi_exx)
    ALLOCATE(hpsi_exx(npwx*npol, max_ibands ))
    !
    ! allocate communication arrays for the exx to local transformation
    !
    IF ( .not.allocated(comm_recv_reverse) ) THEN
       ALLOCATE(comm_recv_reverse(nproc_egrp,nks))
       ALLOCATE(comm_send_reverse(nproc_egrp,negrp,nks))
    END IF
    !
    egrp_base = my_egrp_id*nproc_egrp
    DO ik = 1, nks
       total_lda_egrp(ik) = &
            sum( lda_local(egrp_base+1:(egrp_base+nproc_egrp),ik) )
       prev_lda_egrp(ik) = &
            sum( lda_local(egrp_base+1:(egrp_base+me_egrp),ik) )
    END DO
    !
    max_lda_egrp = 0
    DO j = 1, negrp
       DO ik = 1, nks
          max_lda_egrp = max(max_lda_egrp, &
               sum( lda_local((j-1)*nproc_egrp+1:j*nproc_egrp,ik) ) )
       END DO
    END DO
    !
    ! determine the e2l_map
    !
    allocate( e2l_map(max_lda_egrp,nks,negrp) )
    e2l_map = 0
    DO ik = 1, nks
       DO ig = 1, lda_local(me_pool+1,ik)
          DO j=1, total_lda_exx(ik)
             IF( local_map(ig+prev_lda(ik),ik).EQ.exx_map(j,ik) ) exit
          END DO
          e2l_map(ig+prev_lda_egrp(ik),ik,my_egrp_id+1) = j
       END DO
    END DO
    CALL mp_sum(e2l_map(:,:,my_egrp_id+1),intra_egrp_comm)
    CALL mp_sum(e2l_map,inter_egrp_comm)
    !
    ! plan communication for the local to EXX data structure transformation
    !
    allocate(psi_source_exx( max_lda_egrp, negrp ))
    DO ik = 1, nks
       !
       ! determine where each value is coming from
       !
       psi_source_exx = 0
       DO ig = 1, lda_local(me_pool+1,ik)
          j = 1
          DO i = 1, nproc_egrp
             j = j + lda_exx(i,ik)
             IF( j.gt.e2l_map(ig+prev_lda_egrp(ik),ik,my_egrp_id+1) ) exit
          END DO
          psi_source_exx(ig+prev_lda_egrp(ik),my_egrp_id+1) = i-1
       END DO
       CALL mp_sum(psi_source_exx(:,my_egrp_id+1),intra_egrp_comm)
       CALL mp_sum(psi_source_exx,inter_egrp_comm)
       !
       ! allocate communication packets to receive psi and hpsi (reverse)
       !
       DO iegrp=my_egrp_id+1, my_egrp_id+1
          DO iproc=0, nproc_egrp-1
             !
             ! determine how many values need to come from iproc
             !
             count = 0
             DO ig=1, lda_local(me_pool+1,ik)
                IF ( psi_source_exx(ig+prev_lda_egrp(ik),iegrp).eq.iproc ) THEN
                   count = count + 1
                END IF
             END DO
             !
             ! allocate the communication packet
             !
             comm_recv_reverse(iproc+1,ik)%size = count
             IF (count.gt.0) THEN
                IF (.not.ALLOCATED(comm_recv_reverse(iproc+1,ik)%msg)) THEN
                   ALLOCATE(comm_recv_reverse(iproc+1,ik)%indices(count))
                   ALLOCATE(comm_recv_reverse(iproc+1,ik)%msg(count,npol,m+2))
                END IF
             END IF
             !
             ! determine which values need to come from iproc
             !
             count = 0
             DO ig=1, lda_local(me_pool+1,ik)
                IF ( psi_source_exx(ig+prev_lda_egrp(ik),iegrp).eq.iproc ) THEN
                   count = count + 1
                   comm_recv_reverse(iproc+1,ik)%indices(count) = ig
                END IF
             END DO
             
          END DO
       END DO
       !
       ! allocate communication packets to send psi and hpsi
       !
       DO iegrp=1, negrp
          prev = 0
          DO iproc=0, nproc_egrp-1
             !
             ! determine how many values need to be sent to iproc
             !
             count = 0
             DO ig=1, lda_local(iproc+(iegrp-1)*nproc_egrp+1,ik)
                IF ( psi_source_exx(ig+prev,iegrp).eq.me_egrp ) THEN
                   count = count + 1
                END IF
             END DO
             !
             ! allocate the communication packet
             !
             comm_send_reverse(iproc+1,iegrp,ik)%size = count
             IF (count.gt.0) THEN
                IF (.not.ALLOCATED(comm_send_reverse(iproc+1,iegrp,ik)%msg))THEN
                   ALLOCATE(comm_send_reverse(iproc+1,iegrp,ik)%indices(count))
                   ALLOCATE(comm_send_reverse(iproc+1,iegrp,ik)%msg(count,npol,&
                        iexx_iend(my_egrp_id+1)-iexx_istart(my_egrp_id+1)+3))
                END IF
             END IF
             !
             ! determine which values need to be sent to iproc
             !
             count = 0
             DO ig=1, lda_local(iproc+(iegrp-1)*nproc_egrp+1,ik)
                IF ( psi_source_exx(ig+prev,iegrp).eq.me_egrp ) THEN
                   count = count + 1
                   comm_send_reverse(iproc+1,iegrp,ik)%indices(count) = &
                        e2l_map(ig+prev,ik,iegrp)
                END IF
             END DO
             !
             prev = prev + lda_local( (iegrp-1)*nproc_egrp+iproc+1, ik )
             !
          END DO
       END DO
       !
    END DO
    !
    ! deallocate arrays
    !
    DEALLOCATE( local_map, exx_map )
    DEALLOCATE( l2e_map, e2l_map )
    DEALLOCATE( psi_source, psi_source_exx )
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE initialize_local_to_exact_map
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE transform_to_exx(lda, n, m, m_out, ik, psi, psi_out, type)
  !-----------------------------------------------------------------------
    !
    ! transform psi into the EXX data structure, and place the result in psi_out
    !
    USE wvfct,        ONLY : nbnd
    USE mp,           ONLY : mp_sum
    USE mp_pools,     ONLY : nproc_pool, me_pool
    USE mp_exx,       ONLY : intra_egrp_comm, inter_egrp_comm, &
         nproc_egrp, me_egrp, negrp, my_egrp_id, nibands, ibands, &
         max_ibands, all_start, all_end
    USE parallel_include
    !
    !
    IMPLICIT NONE
    !
    Integer :: lda
    INTEGER :: n, m, m_out
    COMPLEX(DP) :: psi(npwx_local*npol,m) 
    COMPLEX(DP) :: psi_out(npwx_exx*npol,m_out)
    INTEGER, INTENT(in) :: type

    COMPLEX(DP), ALLOCATABLE :: psi_work(:,:,:,:), psi_gather(:,:)
    INTEGER :: i, j, im, iproc, ig, ik, iegrp
    INTEGER :: prev, lda_max_local

    INTEGER :: request_send(nproc_egrp), request_recv(nproc_egrp)
    INTEGER :: ierr
    INTEGER :: current_ik
    !INTEGER, EXTERNAL :: find_current_k
    INTEGER :: recvcount(negrp)
    INTEGER :: displs(negrp)
#if defined(__MPI)
    INTEGER :: istatus(MPI_STATUS_SIZE)
#endif
    INTEGER :: requests(max_ibands+2,negrp)
    !
    INTEGER :: ipol, my_lda, lda_offset, count
    lda_max_local = maxval(lda_local)
    current_ik = ik
    !
    !-------------------------------------------------------!
    !Communication Part 1
    !-------------------------------------------------------!
    !
    allocate(psi_work(lda_max_local,npol,m,negrp))
    allocate(psi_gather(lda_max_local*npol,m))
    DO im=1, m
       my_lda = lda_local(me_pool+1,current_ik)
       DO ipol=1, npol
          lda_offset = lda_max_local*(ipol-1)
          DO ig=1, my_lda
             psi_gather(lda_offset+ig,im) = psi(npwx_local*(ipol-1)+ig,im)
          END DO
       END DO
    END DO
    recvcount = lda_max_local*npol
    count = lda_max_local*npol
#if defined (__MPI_NONBLOCKING)
    IF ( type.eq.0 ) THEN
       DO iegrp=1, negrp
          displs(iegrp) = (iegrp-1)*(count*m)
       END DO
       DO iegrp=1, negrp
          DO im=1, nibands(iegrp)
             IF ( my_egrp_id.eq.(iegrp-1) ) THEN
                DO j=1, negrp
                   displs(j) = (j-1)*(count*m) + count*(ibands(im,iegrp)-1)
                END DO
             END IF
#if defined(__MPI)
             CALL MPI_IGATHERV( psi_gather(:, ibands(im,iegrp) ), &
                  count, MPI_DOUBLE_COMPLEX, &
                  psi_work, &
                  recvcount, displs, MPI_DOUBLE_COMPLEX, &
                  iegrp-1, &
                  inter_egrp_comm, requests(im,iegrp), ierr )
#endif
          END DO
       END DO

    ELSE IF(type.eq.1) THEN
#elif defined(__MPI)
       CALL MPI_ALLGATHER( psi_gather, &
            count*m, MPI_DOUBLE_COMPLEX, &
            psi_work, &
            count*m, MPI_DOUBLE_COMPLEX, &
            inter_egrp_comm, ierr )
#endif

#if defined (__MPI_NONBLOCKING)
    ELSE IF(type.eq.2) THEN !evc2

       DO iegrp=1, negrp
          displs(iegrp) = (iegrp-1)*(count*m)
       END DO
       DO iegrp=1, negrp
          DO im=1, all_end(iegrp) - all_start(iegrp) + 1
             !
             IF(all_start(iegrp).eq.0) CYCLE
             !
             IF ( my_egrp_id.eq.(iegrp-1) ) THEN
                DO j=1, negrp
                   displs(j) = (j-1)*(count*m) + count*(im+all_start(iegrp)-2)
                END DO
             END IF
#if defined(__MPI)
             CALL MPI_IGATHERV( psi_gather(:, im+all_start(iegrp)-1 ), &
                  count, MPI_DOUBLE_COMPLEX, &
                  psi_work, &
                  recvcount, displs, MPI_DOUBLE_COMPLEX, &
                  iegrp-1, &
                  inter_egrp_comm, requests(im,iegrp), ierr )
#endif
          END DO
       END DO
          
    END IF
    !
    IF(type.eq.0)THEN
       DO iegrp=1, negrp
          DO im=1, nibands(iegrp)
             CALL MPI_WAIT(requests(im,iegrp), istatus, ierr)
          END DO
       END DO
    ELSEIF(type.eq.2)THEN
       DO iegrp=1, negrp
          DO im=1, all_end(iegrp) - all_start(iegrp) + 1
             IF(all_start(iegrp).eq.0) CYCLE
             CALL MPI_WAIT(requests(im,iegrp), istatus, ierr)
          END DO
       END DO
    END IF
#endif
    !
    !-------------------------------------------------------!
    !Communication Part 2
    !-------------------------------------------------------!
    !
    !send communication packets
    !
    DO iproc=0, nproc_egrp-1
       IF ( comm_send(iproc+1,current_ik)%size.gt.0) THEN
          DO i=1, comm_send(iproc+1,current_ik)%size
             ig = comm_send(iproc+1,current_ik)%indices(i)
             !
             !determine which egrp this corresponds to
             !
             prev = 0
             DO j=1, nproc_pool
                IF ((prev+lda_local(j,current_ik)).ge.ig) THEN 
                   ig = ig - prev
                   exit
                END IF
                prev = prev + lda_local(j,current_ik)
             END DO
             !
             ! prepare the message
             !
             DO ipol=1, npol
             IF ( type.eq.0 ) THEN !psi or hpsi
                DO im=1, nibands(my_egrp_id+1)
                   comm_send(iproc+1,current_ik)%msg(i,ipol,im) = &
                        psi_work(ig,ipol,ibands(im,my_egrp_id+1),1+(j-1)/nproc_egrp)
                END DO
             ELSE IF (type.eq.1) THEN !evc
                DO im=1, m
                   comm_send(iproc+1,current_ik)%msg(i,ipol,im) = &
                        psi_work(ig,ipol,im,1+(j-1)/nproc_egrp)
                END DO
             ELSE IF ( type.eq.2 ) THEN !evc2
                IF(all_start(my_egrp_id+1).gt.0) THEN
                   DO im=1, all_end(my_egrp_id+1) - all_start(my_egrp_id+1) + 1
                      comm_send(iproc+1,current_ik)%msg(i,ipol,im) = &
                           psi_work(ig,ipol,im+all_start(my_egrp_id+1)-1,1+(j-1)/nproc_egrp)
                   END DO
                END IF
             END IF
             END DO
             !
          END DO
          !
          ! send the message
          !
#if defined(__MPI)
          IF ( type.eq.0 ) THEN !psi or hpsi
             CALL MPI_ISEND( comm_send(iproc+1,current_ik)%msg, &
                  comm_send(iproc+1,current_ik)%size*npol*nibands(my_egrp_id+1), &
                  MPI_DOUBLE_COMPLEX, &
                  iproc, 100+iproc*nproc_egrp+me_egrp, &
                  intra_egrp_comm, request_send(iproc+1), ierr )
          ELSE IF (type.eq.1) THEN !evc
             CALL MPI_ISEND( comm_send(iproc+1,current_ik)%msg, &
                  comm_send(iproc+1,current_ik)%size*npol*m, MPI_DOUBLE_COMPLEX, &
                  iproc, 100+iproc*nproc_egrp+me_egrp, &
                  intra_egrp_comm, request_send(iproc+1), ierr )
          ELSE IF (type.eq.2) THEN !evc2
             CALL MPI_ISEND( comm_send(iproc+1,current_ik)%msg, &
                  comm_send(iproc+1,current_ik)%size*npol*(all_end(my_egrp_id+1)-all_start(my_egrp_id+1)+1), &
                  MPI_DOUBLE_COMPLEX, &
                  iproc, 100+iproc*nproc_egrp+me_egrp, &
                  intra_egrp_comm, request_send(iproc+1), ierr )
          END IF
#endif
          !
       END IF
    END DO
    !
    ! begin receiving the messages
    !
    DO iproc=0, nproc_egrp-1
       IF ( comm_recv(iproc+1,current_ik)%size.gt.0) THEN
          !
          ! receive the message
          !
#if defined(__MPI)
          IF (type.eq.0) THEN !psi or hpsi
             CALL MPI_IRECV( comm_recv(iproc+1,current_ik)%msg, &
                  comm_recv(iproc+1,current_ik)%size*npol*nibands(my_egrp_id+1), &
                  MPI_DOUBLE_COMPLEX, &
                  iproc, 100+me_egrp*nproc_egrp+iproc, &
                  intra_egrp_comm, request_recv(iproc+1), ierr )
          ELSE IF (type.eq.1) THEN !evc
             CALL MPI_IRECV( comm_recv(iproc+1,current_ik)%msg, &
                  comm_recv(iproc+1,current_ik)%size*npol*m, MPI_DOUBLE_COMPLEX, &
                  iproc, 100+me_egrp*nproc_egrp+iproc, &
                  intra_egrp_comm, request_recv(iproc+1), ierr )
          ELSE IF (type.eq.2) THEN !evc2
             CALL MPI_IRECV( comm_recv(iproc+1,current_ik)%msg, &
                  comm_recv(iproc+1,current_ik)%size*npol*(all_end(my_egrp_id+1)-all_start(my_egrp_id+1)+1), &
                  MPI_DOUBLE_COMPLEX, &
                  iproc, 100+me_egrp*nproc_egrp+iproc, &
                  intra_egrp_comm, request_recv(iproc+1), ierr )
          END IF
#endif
          !
       END IF
    END DO
    !
    ! assign psi_out
    !
    DO iproc=0, nproc_egrp-1
       IF ( comm_recv(iproc+1,current_ik)%size.gt.0 ) THEN
          !
          ! wait for the message to be received
          !
#if defined(__MPI)
          CALL MPI_WAIT(request_recv(iproc+1), istatus, ierr)
#endif
          !
          DO i=1, comm_recv(iproc+1,current_ik)%size
             ig = comm_recv(iproc+1,current_ik)%indices(i)
             !
             ! place the message into the correct elements of psi_out
             !
             IF (type.eq.0) THEN !psi or hpsi
                DO im=1, nibands(my_egrp_id+1)
                   DO ipol=1, npol
                      psi_out(ig+npwx_exx*(ipol-1),im) = &
                           comm_recv(iproc+1,current_ik)%msg(i,ipol,im)
                   END DO
                END DO
             ELSE IF (type.eq.1) THEN !evc
                DO im=1, m
                   DO ipol=1, npol
                      psi_out(ig+npwx_exx*(ipol-1),im) = &
                           comm_recv(iproc+1,current_ik)%msg(i,ipol,im)
                   END DO
                END DO
             ELSE IF (type.eq.2) THEN !evc2
                DO im=1, all_end(my_egrp_id+1) - all_start(my_egrp_id+1) + 1
                   DO ipol=1, npol
                      psi_out(ig+npwx_exx*(ipol-1),im) = comm_recv(iproc+1,current_ik)%msg(i,ipol,im)
                   END DO
                END DO
             END IF
             !
          END DO
          !
       END IF
    END DO
    !
    ! wait for everything to finish sending
    !
#if defined(__MPI)
    DO iproc=0, nproc_egrp-1
       IF ( comm_send(iproc+1,current_ik)%size.gt.0 ) THEN
          CALL MPI_WAIT(request_send(iproc+1), istatus, ierr)
       END IF
    END DO
#endif
    !
    ! deallocate arrays
    !
    DEALLOCATE( psi_work, psi_gather )
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE transform_to_exx
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE change_data_structure(is_exx)
  !-----------------------------------------------------------------------
    !
    ! change between the local and EXX data structures
    ! is_exx = .TRUE. - change to the EXX data structure
    ! is_exx = .FALSE. - change to the local data strucutre
    !
    USE cell_base,      ONLY : at, bg, tpiba2
    USE cellmd,         ONLY : lmovecell
    USE wvfct,          ONLY : npwx
    USE gvect,          ONLY : gcutm, ig_l2g, g, gg, ngm, ngm_g, mill, &
                               gstart, gvect_init, deallocate_gvect_exx
    USE gvecs,          ONLY : gcutms, ngms, ngms_g, gvecs_init
    USE gvecw,          ONLY : gkcut, ecutwfc, gcutw
    USE klist,          ONLY : xk, nks, ngk
    USE mp_bands,       ONLY : intra_bgrp_comm, ntask_groups, nyfft
    USE mp_exx,         ONLY : intra_egrp_comm, me_egrp, exx_mode, nproc_egrp, &
                               negrp, root_egrp
    USE io_global,      ONLY : stdout
    USE fft_base,       ONLY : dfftp, dffts, smap, fft_base_info
    USE fft_types,      ONLY : fft_type_init
    USE recvec_subs,    ONLY : ggen, ggens
    !
    IMPLICIT NONE
    !
    LOGICAL, intent(in) :: is_exx
    INTEGER, EXTERNAL  :: n_plane_waves
    COMPLEX(DP), ALLOCATABLE :: work_space(:)
    INTEGER :: ik, i
    INTEGER :: ngm_, ngs_, ngw_
    LOGICAL exst
#if defined (__MPI)
  LOGICAL :: lpara = .true.
#else
  LOGICAL :: lpara = .false.
#endif
    !
    IF (negrp.eq.1) RETURN
    !
    IF (first_data_structure_change) THEN
       allocate( ig_l2g_loc(ngm), g_loc(3,ngm), gg_loc(ngm) )
       allocate( mill_loc(3,ngm), nl_loc(ngm) )
       allocate( nls_loc(size(dffts%nl)) )
       IF ( gamma_only) THEN
          allocate( nlm_loc(size(dfftp%nlm)) )
          allocate( nlsm_loc(size(dffts%nlm)) )
       END IF
       ig_l2g_loc = ig_l2g
       g_loc = g
       gg_loc = gg
       mill_loc = mill
       nl_loc = dfftp%nl
       nls_loc = dffts%nl
       IF ( gamma_only) THEN
          nlm_loc = dfftp%nlm
          nlsm_loc = dffts%nlm
       END IF
       ngm_loc = ngm
       ngm_g_loc = ngm_g
       gstart_loc = gstart
       ngms_loc = ngms
       ngms_g_loc = ngms_g
    END IF
    !
    ! generate the gvectors for the new data structure
    !
    IF (is_exx) THEN
       exx_mode = 1
       IF(first_data_structure_change)THEN
          dfftp_loc = dfftp
          dffts_loc = dffts

          CALL fft_type_init( dffts_exx, smap_exx, "wave", gamma_only, &
               lpara, intra_egrp_comm, at, bg, gkcut, gcutms/gkcut, &
               nyfft=ntask_groups )
          CALL fft_type_init( dfftp_exx, smap_exx, "rho", gamma_only, &
               lpara, intra_egrp_comm, at, bg,  gcutm, nyfft=nyfft )
          CALL fft_base_info( ionode, stdout )
          ngs_ = dffts_exx%ngl( dffts_exx%mype + 1 )
          ngm_ = dfftp_exx%ngl( dfftp_exx%mype + 1 )
          IF( gamma_only ) THEN
             ngs_ = (ngs_ + 1)/2
             ngm_ = (ngm_ + 1)/2
          END IF
          dfftp = dfftp_exx
          dffts = dffts_exx
          ngm = ngm_
          ngms = ngs_
       ELSE
          dfftp = dfftp_exx
          dffts = dffts_exx
          ngm = ngm_exx
          ngms = ngms_exx
       END IF
       call deallocate_gvect_exx()
       call gvect_init( ngm , intra_egrp_comm )
       call gvecs_init( ngms , intra_egrp_comm )
    ELSE
       exx_mode = 2
       dfftp = dfftp_loc
       dffts = dffts_loc
       ngm = ngm_loc
       ngms = ngms_loc
       call deallocate_gvect_exx()
       call gvect_init( ngm , intra_bgrp_comm )
       call gvecs_init( ngms , intra_bgrp_comm )
       exx_mode = 0
    END IF
    !
    IF (first_data_structure_change) THEN
       CALL ggen ( dfftp, gamma_only, at, bg, gcutm, ngm_g, ngm, &
            g, gg, mill, ig_l2g, gstart )
       CALL ggens( dffts, gamma_only, at, g, gg, mill, gcutms, ngms )
       allocate( ig_l2g_exx(ngm), g_exx(3,ngm), gg_exx(ngm) )
       allocate( mill_exx(3,ngm), nl_exx(ngm) )
       allocate( nls_exx(size(dffts%nl)) )
       allocate( nlm_exx(size(dfftp%nlm) ) )
       allocate( nlsm_exx(size(dffts%nlm) ) )
       ig_l2g_exx = ig_l2g
       g_exx = g
       gg_exx = gg
       mill_exx = mill
       nl_exx = dfftp%nl
       nls_exx = dffts%nl
       IF( gamma_only ) THEN
          nlm_exx = dfftp%nlm
          nlsm_exx = dffts%nlm
       END IF
       ngm_exx = ngm
       ngm_g_exx = ngm_g
       gstart_exx = gstart
       ngms_exx = ngms
       ngms_g_exx = ngms_g
    ELSE IF ( is_exx ) THEN
       ig_l2g = ig_l2g_exx
       g = g_exx
       gg = gg_exx
       mill = mill_exx
       ! workaround: here dfft?%nl* are unallocated
       ! some compilers go on and allocate, some others crash
       IF ( .NOT. ALLOCATED(dfftp%nl) ) ALLOCATE (dfftp%nl(size(nl_exx)))
       IF ( .NOT. ALLOCATED(dffts%nl) ) ALLOCATE (dffts%nl(size(nls_exx)))
       IF ( gamma_only .AND. .NOT.ALLOCATED(dfftp%nlm) ) ALLOCATE (dfftp%nlm(size(nlm_exx)))
       IF ( gamma_only .AND. .NOT.ALLOCATED(dffts%nlm) ) ALLOCATE (dffts%nlm(size(nlsm_exx)))
       ! end workaround. FIXME: this part of code must disappear ASAP
       dfftp%nl = nl_exx
       dffts%nl = nls_exx
       IF( gamma_only ) THEN
          dfftp%nlm = nlm_exx
          dffts%nlm = nlsm_exx
       ENDIF
       ngm = ngm_exx
       ngm_g = ngm_g_exx
       gstart = gstart_exx
       ngms = ngms_exx
       ngms_g = ngms_g_exx
    ELSE ! not is_exx
       ig_l2g = ig_l2g_loc
       g = g_loc
       gg = gg_loc
       mill = mill_loc
       dfftp%nl = nl_loc
       dffts%nl = nls_loc
       IF( gamma_only ) THEN
          dfftp%nlm = nlm_loc
          dffts%nlm = nlsm_loc
       END IF
       ngm = ngm_loc
       ngm_g = ngm_g_loc
       gstart = gstart_loc
       ngms = ngms_loc
       ngms_g = ngms_g_loc
    END IF
    !
    ! get npwx and ngk
    !
    IF ( is_exx.and.npwx_exx.gt.0 ) THEN
       npwx = npwx_exx
       ngk = ngk_exx
    ELSE IF ( .not.is_exx.and.npwx_local.gt.0 ) THEN
       npwx = npwx_local
       ngk = ngk_local
    ELSE
       npwx = n_plane_waves (gcutw, nks, xk, g, ngm)
    END IF
    !
    ! get igk
    !
    IF( first_data_structure_change ) THEN
       allocate(igk_exx(npwx,nks),work_space(npwx))
       first_data_structure_change = .FALSE.
       IF ( nks.eq.1 ) THEN
          CALL gk_sort( xk, ngm, g, ecutwfc / tpiba2, ngk, igk_exx, work_space )
       END IF
       IF ( nks > 1 ) THEN
          !
          DO ik = 1, nks
             CALL gk_sort( xk(1,ik), ngm, g, ecutwfc / tpiba2, ngk(ik), &
                  igk_exx(1,ik), work_space )
          END DO
          !
       END IF
       DEALLOCATE( work_space )
    END IF
    !
    ! generate ngl and igtongl
    !
    CALL gshells( lmovecell )
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE change_data_structure
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE update_igk(is_exx)
  !-----------------------------------------------------------------------
    USE cell_base,      ONLY : tpiba2
    USE gvect,          ONLY : ngm, g
    USE gvecw,          ONLY : ecutwfc
    USE wvfct,          ONLY : npwx, current_k
    USE klist,          ONLY : xk, igk_k
    USE mp_exx,         ONLY : negrp
    !
    !
    IMPLICIT NONE
    !
    LOGICAL, intent(in) :: is_exx
    COMPLEX(DP), ALLOCATABLE :: work_space(:)
    INTEGER :: comm
    INTEGER :: npw, ik, i
    LOGICAL exst
    !
    IF (negrp.eq.1) RETURN
    !
    ! get igk
    !
    allocate(work_space(npwx))
    ik = current_k
    IF(is_exx) THEN
       CALL gk_sort( xk(1,ik), ngm, g, ecutwfc / tpiba2, npw, igk_exx(1,ik), &
            work_space )
    ELSE
       CALL gk_sort( xk(1,ik), ngm, g, ecutwfc / tpiba2, npw, igk_k(1,ik), &
            work_space )
    END IF
    !
    DEALLOCATE( work_space )
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE update_igk
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE result_sum (n, m, data)
  !-----------------------------------------------------------------------
    USE parallel_include
    USE mp_exx,       ONLY : iexx_start, iexx_end, inter_egrp_comm, &
                               intra_egrp_comm, my_egrp_id, negrp, &
                               max_pairs, egrp_pairs, max_contributors, &
                               contributed_bands, all_end, &
                               iexx_istart, iexx_iend, band_roots
    USE mp,           ONLY : mp_sum, mp_bcast

    INTEGER, INTENT(in) :: n, m
    COMPLEX(DP), INTENT(inout) :: data(n,m)
#if defined(__MPI)
    INTEGER :: istatus(MPI_STATUS_SIZE)
    COMPLEX(DP), ALLOCATABLE :: recvbuf(:,:)
    COMPLEX(DP) :: data_sum(n,m), test(negrp)
    INTEGER :: im, iegrp, ibuf, i, j, nsending(m)
    INTEGER :: ncontributing(m)
    INTEGER :: contrib_this(negrp,m), displs(negrp,m)

    INTEGER sendcount, sendtype, ierr, root, request(m)
    INTEGER sendc(negrp), sendd(negrp)
#endif
    !
    IF (negrp.eq.1) RETURN
    !
#if defined(__MPI)
    ! gather data onto the correct nodes
    !
    ALLOCATE( recvbuf( n*max_contributors, max(1,iexx_end-iexx_start+1) ) )
    displs = 0
    ibuf = 0
    nsending = 0
    contrib_this = 0
    !
    DO im=1, m
       !
       IF(contributed_bands(im,my_egrp_id+1)) THEN
          sendcount = n
       ELSE
          sendcount = 0
       END IF
       !
       root = band_roots(im)
       !
       IF(my_egrp_id.eq.root) THEN
          !
          ! determine the number of sending processors
          !
          ibuf = ibuf + 1
          ncontributing(im) = 0
          DO iegrp=1, negrp
             IF(contributed_bands(im,iegrp)) THEN
                ncontributing(im) = ncontributing(im) + 1
                contrib_this(iegrp,im) = n
                IF(iegrp.lt.negrp) displs(iegrp+1,im) = displs(iegrp,im) + n
                nsending(im) = nsending(im) + 1
             ELSE
                contrib_this(iegrp,im) = 0
                IF(iegrp.lt.negrp) displs(iegrp+1,im) = displs(iegrp,im)
             END IF
          END DO
       END IF
       !
#if defined(__MPI_NONBLOCKING)
       CALL MPI_IGATHERV(data(:,im), sendcount, MPI_DOUBLE_COMPLEX, &
            recvbuf(:,max(1,ibuf)), contrib_this(:,im), &
            displs(:,im), MPI_DOUBLE_COMPLEX, &
            root, inter_egrp_comm, request(im), ierr)
#else
       CALL MPI_GATHERV(data(:,im), sendcount, MPI_DOUBLE_COMPLEX, &
            recvbuf(:,max(1,ibuf)), contrib_this(:,im), &
            displs(:,im), MPI_DOUBLE_COMPLEX, &
            root, inter_egrp_comm, ierr)
#endif
       !
    END DO
    !
#if defined(__MPI_NONBLOCKING)
    DO im=1, m
       CALL MPI_WAIT(request(im), istatus, ierr)
    END DO
#endif
    !
    ! perform the sum
    !
    DO im=iexx_istart(my_egrp_id+1), iexx_iend(my_egrp_id+1)
       IF(im.eq.0)exit
       data(:,im) = 0._dp
       ibuf = im - iexx_istart(my_egrp_id+1) + 1
       DO j=1, nsending(im)
          DO i=1, n
             data(i,im) = data(i,im) + recvbuf(i+n*(j-1),ibuf)
          END DO
       END DO
    END DO
#endif
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE result_sum
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE transform_to_local(m, m_exx, psi, psi_out)
  !-----------------------------------------------------------------------
    USE mp,           ONLY : mp_sum
    USE mp_pools,     ONLY : nproc_pool, me_pool, intra_pool_comm
    USE mp_exx,       ONLY : intra_egrp_comm, inter_egrp_comm, &
         nproc_egrp, me_egrp, negrp, my_egrp_id, iexx_istart, iexx_iend
    USE parallel_include
    USE wvfct,        ONLY : current_k
    !
    !
    IMPLICIT NONE
    !
    INTEGER :: m, m_exx
    COMPLEX(DP) :: psi(npwx_exx*npol,m_exx)
    COMPLEX(DP) :: psi_out(npwx_local*npol,m)
    !
    INTEGER :: i, j, im, iproc, ig, ik, current_ik, iegrp
    INTEGER :: prev, lda_max_local, prev_lda_exx
    INTEGER :: my_bands, recv_bands, tag
    !
    INTEGER :: request_send(nproc_egrp,negrp), request_recv(nproc_egrp,negrp)
    INTEGER :: ierr
    INTEGER :: ipol
#if defined(__MPI)
    INTEGER :: istatus(MPI_STATUS_SIZE)
#endif
    !INTEGER, EXTERNAL :: find_current_k
    !
    current_ik = current_k
    prev_lda_exx = sum( lda_exx(1:me_egrp,current_ik) )
    !
    my_bands = iexx_iend(my_egrp_id+1) - iexx_istart(my_egrp_id+1) + 1
    !
    ! send communication packets
    !
    IF ( iexx_istart(my_egrp_id+1).gt.0 ) THEN
       DO iegrp=1, negrp
          DO iproc=0, nproc_egrp-1
             IF ( comm_send_reverse(iproc+1,iegrp,current_ik)%size.gt.0) THEN
                DO i=1, comm_send_reverse(iproc+1,iegrp,current_ik)%size
                   ig = comm_send_reverse(iproc+1,iegrp,current_ik)%indices(i)
                   ig = ig - prev_lda_exx
                   !
                   DO im=1, my_bands
                      DO ipol=1, npol
                         comm_send_reverse(iproc+1,iegrp,current_ik)%msg(i,ipol,im) = &
                              psi(ig+npwx_exx*(ipol-1),im)
                      END DO
                   END DO
                   !
                END DO
                !
                ! send the message
                !
                tag = 0
#if defined(__MPI)
                CALL MPI_ISEND( comm_send_reverse(iproc+1,iegrp,current_ik)%msg, &
                     comm_send_reverse(iproc+1,iegrp,current_ik)%size*npol*my_bands, &
                     MPI_DOUBLE_COMPLEX, &
                     iproc+(iegrp-1)*nproc_egrp, &
                     tag, &
                     intra_pool_comm, request_send(iproc+1,iegrp), ierr )
#endif
                !
             END IF
          END DO
       END DO
    END IF
    !
    ! begin receiving the communication packets
    !
    DO iegrp=1, negrp
       !
       IF ( iexx_istart(iegrp).le.0 ) CYCLE
       !
       recv_bands = iexx_iend(iegrp) - iexx_istart(iegrp) + 1
       !
       DO iproc=0, nproc_egrp-1
          IF ( comm_recv_reverse(iproc+1,current_ik)%size.gt.0) THEN
             !
             !receive the message
             !
             tag = 0
#if defined(__MPI)
             CALL MPI_IRECV( comm_recv_reverse(iproc+1,current_ik)%msg(:,:,iexx_istart(iegrp)), &
                  comm_recv_reverse(iproc+1,current_ik)%size*npol*recv_bands, &
                  MPI_DOUBLE_COMPLEX, &
                  iproc+(iegrp-1)*nproc_egrp, &
                  tag, &
                  intra_pool_comm, request_recv(iproc+1,iegrp), ierr )
#endif
             !
          END IF
       END DO
    END DO
    !
    ! assign psi
    !
    DO iproc=0, nproc_egrp-1
       IF ( comm_recv_reverse(iproc+1,current_ik)%size.gt.0 ) THEN
          !
#if defined(__MPI)
          DO iegrp=1, negrp
             IF ( iexx_istart(iegrp).le.0 ) CYCLE
             CALL MPI_WAIT(request_recv(iproc+1,iegrp), istatus, ierr)
          END DO
#endif
          !
          DO i=1, comm_recv_reverse(iproc+1,current_ik)%size
             ig = comm_recv_reverse(iproc+1,current_ik)%indices(i)
             !
             ! set psi_out
             !
             DO im=1, m
                DO ipol=1, npol
                   psi_out(ig+npwx_local*(ipol-1),im) = psi_out(ig+npwx_local*(ipol-1),im) + &
                        comm_recv_reverse(iproc+1,current_ik)%msg(i,ipol,im)
                END DO
             END DO
             !
          END DO
          !
       END IF
    END DO
    !
    ! wait for everything to finish sending
    !
#if defined(__MPI)
    IF ( iexx_istart(my_egrp_id+1).gt.0 ) THEN
       DO iproc=0, nproc_egrp-1
          DO iegrp=1, negrp
             IF ( comm_send_reverse(iproc+1,iegrp,current_ik)%size.gt.0 ) THEN
                CALL MPI_WAIT(request_send(iproc+1,iegrp), istatus, ierr)
             END IF
          END DO
       END DO
    END IF
#endif
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE transform_to_local
  !-----------------------------------------------------------------------

END MODULE exx_band
!-----------------------------------------------------------------------
