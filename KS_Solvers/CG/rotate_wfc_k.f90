!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE rotate_wfc_k( h_psi, s_psi, overlap, &
                         npwx, npw, nstart, nbnd, npol, psi, evc, e )
  !----------------------------------------------------------------------------
  !
  ! ... Serial version of rotate_wfc for colinear, k-point calculations
  !
  USE cg_param,      ONLY : DP
  USE mp_bands_util, ONLY : intra_bgrp_comm, inter_bgrp_comm, root_bgrp_id,&
          nbgrp, my_bgrp_id, me_bgrp, root_bgrp
  USE mp,            ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  INCLUDE 'laxlib.fh'
  !
  ! ... I/O variables
  !
  INTEGER, INTENT(IN) :: npw, npwx, nstart, nbnd, npol
    ! dimension of the matrix to be diagonalized
    ! leading dimension of matrix psi, as declared in the calling pgm unit
    ! input number of states
    ! output number of states
    ! number of spin polarizations
  LOGICAL :: overlap
    ! if .FALSE. : S|psi> not needed
  COMPLEX(DP) :: psi(npwx*npol,nstart), evc(npwx*npol,nbnd)
    ! input and output eigenvectors (may overlap)
  REAL(DP) :: e(nbnd)
    ! eigenvalues
  !
  ! ... local variables
  !
  INTEGER :: kdim, kdmx
  COMPLEX(DP), ALLOCATABLE :: aux(:,:)
  COMPLEX(DP), ALLOCATABLE :: hc(:,:), sc(:,:), vc(:,:)
  REAL(DP),    ALLOCATABLE :: en(:)
  INTEGER :: n_start, n_end, my_n
  !
  EXTERNAL  h_psi,    s_psi
    ! h_psi(npwx,npw,nvec,psi,hpsi)
    !     calculates H|psi>
    ! s_psi(npwx,npw,nvec,spsi)
    !     calculates S|psi> (if needed)
    !     Vectors psi,hpsi,spsi are dimensioned (npwx,npol,nvec)

  IF ( npol == 1 ) THEN
     !
     kdim = npw
     kdmx = npwx
     !
  ELSE
     !
     kdim = npwx*npol
     kdmx = npwx*npol
     !
  END IF
  !
  ALLOCATE( aux(kdmx, nstart ) )    
  ALLOCATE( hc( nstart, nstart) )    
  ALLOCATE( sc( nstart, nstart) )    
  ALLOCATE( vc( nstart, nstart) )    
  ALLOCATE( en( nstart ) )
  call start_clock('rotwfck'); !write(*,*) 'start rotwfck';FLUSH(6)
  !
  ! ... Set up the Hamiltonian and Overlap matrix on the subspace :
  !
  ! ...      H_ij = <psi_i| H |psi_j>     S_ij = <psi_i| S |psi_j>
  !
  call start_clock('rotwfck:hpsi'); !write(*,*) 'start rotwfck:hpsi';FLUSH(6)
  CALL h_psi( npwx, npw, nstart, psi, aux )
  call stop_clock('rotwfck:hpsi') ; !write(*,*) 'stop rotwfck:hpsi';FLUSH(6)
  !
  call start_clock('rotwfck:hc'); !write(*,*) 'start rotwfck:hc';FLUSH(6)
  hc=(0.D0,0.D0)
  CALL divide(inter_bgrp_comm,nstart,n_start,n_end)
  my_n = n_end - n_start + 1; !write (*,*) nstart,n_start,n_end
  if (n_start .le. n_end) &
  call ZGEMM( 'C','N', nstart, my_n, kdim, (1.D0,0.D0), psi, kdmx, aux(1,n_start), kdmx, (0.D0,0.D0), hc(1,n_start), nstart )
  CALL mp_sum( hc, inter_bgrp_comm )
  !            
  CALL mp_sum( hc, intra_bgrp_comm )
  !
  sc=(0.D0,0.D0)
  IF ( overlap ) THEN
     !
     CALL s_psi( npwx, npw, nstart, psi, aux )
     if (n_start .le. n_end) &
     CALL ZGEMM( 'C','N', nstart, my_n, kdim, (1.D0,0.D0), psi, kdmx, aux(1,n_start), kdmx, (0.D0,0.D0), sc(1,n_start), nstart )
     !
  ELSE
     !
     if (n_start .le. n_end) &
     CALL ZGEMM( 'C','N', nstart, my_n, kdim, (1.D0,0.D0), psi, kdmx, psi(1,n_start), kdmx, (0.D0,0.D0), sc(1,n_start), nstart )
     !  
  END IF
  CALL mp_sum( sc, inter_bgrp_comm )
  !
  CALL mp_sum( sc, intra_bgrp_comm )
  call stop_clock('rotwfck:hc'); !write(*,*) 'stop rotwfck:hc';FLUSH(6)
  !
  ! ... Diagonalize
  !
  call start_clock('rotwfck:diag');  !write(*,*) 'start rotwfck:diag';FLUSH(6)
  CALL diaghg( nstart, nbnd, hc, sc, nstart, en, vc, me_bgrp, root_bgrp, intra_bgrp_comm )
  call stop_clock('rotwfck:diag');  !write(*,*) 'stop rotwfck:diag';FLUSH(6)
  call start_clock('rotwfck:evc'); !write(*,*) 'start rotwfck:evc';FLUSH(6)
  !
  e(:) = en(1:nbnd)
  !
  ! ...  update the basis set
  !  
  aux=(0.D0,0.D0)
  if (n_start .le. n_end) &
  CALL ZGEMM( 'N','N', kdim, nbnd, my_n, (1.D0,0.D0), psi(1,n_start), kdmx, vc(n_start,1), nstart, (0.D0,0.D0), aux, kdmx )
  CALL mp_sum( aux, inter_bgrp_comm )
  !     
  evc(:,:) = aux(:,1:nbnd)
  call stop_clock('rotwfck:evc') ; !write(*,*) 'start rotwfck;evc';FLUSH(6)
  !
  DEALLOCATE( en )
  DEALLOCATE( vc )
  DEALLOCATE( sc )
  DEALLOCATE( hc )
  DEALLOCATE( aux )
  call stop_clock('rotwfck'); !write(*,*) 'stop rotwfck';FLUSH(6)
  !call print_clock('rotwfck')
  !call print_clock('rotwfck:hpsi')
  !call print_clock('rotwfck:hc')
  !call print_clock('rotwfck:diag')
  !call print_clock('rotwfck:evc')
  !
  RETURN
  !
END SUBROUTINE rotate_wfc_k
!
!
!----------------------------------------------------------------------------
SUBROUTINE protate_wfc_k( h_psi, s_psi, overlap, &
                          npwx, npw, nstart, nbnd, npol, psi, evc, e )
  !----------------------------------------------------------------------------
  !
  ! ... Parallel version of rotate_wfc for colinear, k-point calculations
  ! ... Subroutine with distributed matrices, written by Carlo Cavazzoni
  !
  USE cg_param,         ONLY : DP
  USE mp_bands_util,    ONLY : intra_bgrp_comm, inter_bgrp_comm, root_bgrp_id,&
          nbgrp, my_bgrp_id
  USE mp,               ONLY : mp_bcast, mp_root_sum, mp_sum, mp_barrier
  !
  IMPLICIT NONE
  !
  include 'laxlib.fh'
  !
  ! ... I/O variables
  !
  INTEGER :: npw, npwx, nstart, nbnd, npol
    ! dimension of the matrix to be diagonalized
    ! leading dimension of matrix psi, as declared in the calling pgm unit
    ! input number of states
    ! output number of states
    ! number of spin polarizations
  LOGICAL :: overlap
    ! if .FALSE. : S|psi> not needed
  COMPLEX(DP) :: psi(npwx*npol,nstart), evc(npwx*npol,nbnd)
    ! input and output eigenvectors (may overlap)
  REAL(DP) :: e(nbnd)
    ! eigenvalues
  !
  ! ... local variables
  !
  INTEGER :: kdim, kdmx
  COMPLEX(DP), ALLOCATABLE :: aux(:,:)
  COMPLEX(DP), ALLOCATABLE :: hc(:,:), sc(:,:), vc(:,:)
  REAL(DP),    ALLOCATABLE :: en(:)
  !
  INTEGER :: idesc(LAX_DESC_SIZE)
    ! matrix distribution descriptors
  INTEGER :: nx
    ! maximum local block dimension
  LOGICAL :: la_proc
    ! flag to distinguish procs involved in linear algebra
  INTEGER, ALLOCATABLE :: idesc_ip( :, :, : )
  INTEGER, ALLOCATABLE :: rank_ip( :, : )
  !
  INTEGER :: ortho_comm, np_ortho(2), me_ortho(2), ortho_comm_id, leg_ortho, &
             ortho_parent_comm, ortho_cntx
  LOGICAL :: do_distr_diag_inside_bgrp
  !
  EXTERNAL  h_psi,    s_psi
    ! h_psi(npwx,npw,nvec,psi,hpsi)
    !     calculates H|psi>
    ! s_psi(npwx,npw,nvec,spsi)
    !     calculates S|psi> (if needed)
    !     Vectors psi,hpsi,spsi are dimensioned (npwx,npol,nvec)

  call start_clock('protwfck')
  !
  !
  CALL desc_init( nstart, nx, la_proc, idesc, idesc_ip )
  !
  IF ( npol == 1 ) THEN
     !
     kdim = npw
     kdmx = npwx
     !
  ELSE
     !
     kdim = npwx*npol
     kdmx = npwx*npol
     !
  END IF
  !
  ALLOCATE( aux(kdmx, nstart ) )    
  ALLOCATE( hc( nx, nx) )    
  ALLOCATE( sc( nx, nx) )    
  ALLOCATE( vc( nx, nx) )    
  ALLOCATE( en( nstart ) )

  aux=(0.0_DP,0.0_DP)
  !
  ! ... Set up the Hamiltonian and Overlap matrix on the subspace :
  !
  ! ...      H_ij = <psi_i| H |psi_j>     S_ij = <psi_i| S |psi_j>
  !
  call start_clock('protwfck:hpsi')
  CALL h_psi( npwx, npw, nstart, psi, aux )
  call stop_clock('protwfck:hpsi')
  !
  call start_clock('protwfck:hc')
  CALL compute_distmat( hc, psi, aux ) 
  !            
  IF ( overlap ) THEN
     !
     CALL s_psi( npwx, npw, nstart, psi, aux )
     CALL compute_distmat( sc, psi, aux )
     !
  ELSE
     !
     CALL compute_distmat( sc, psi, psi )
     !  
  END IF
  call stop_clock('protwfck:hc')
  !
  ! ... Diagonalize
  !
  call start_clock('protwfck:diag')
  IF ( do_distr_diag_inside_bgrp ) THEN ! NB on output of pdiaghg en and vc are the same across ortho_parent_comm
     ! only the first bgrp performs the diagonalization
     IF( my_bgrp_id == root_bgrp_id ) CALL pdiaghg( nstart, hc, sc, nx, en, vc, idesc )
     IF( nbgrp > 1 ) THEN ! results must be brodcast to the other band groups
       CALL mp_bcast( vc, root_bgrp_id, inter_bgrp_comm )
       CALL mp_bcast( en, root_bgrp_id, inter_bgrp_comm )
     ENDIF
  ELSE
     CALL pdiaghg( nstart, hc, sc, nx, en, vc, idesc )
  END IF
  call stop_clock('protwfck:diag')
  !
  e(:) = en(1:nbnd)
  !
  ! ...  update the basis set
  !  
  call start_clock('protwfck:evc')
  CALL refresh_evc()
  !     
  evc(:,:) = aux(:,1:nbnd)
  call stop_clock('protwfck:evc')
  !
  DEALLOCATE( en )
  DEALLOCATE( vc )
  DEALLOCATE( sc )
  DEALLOCATE( hc )
  DEALLOCATE( aux )
  !
  DEALLOCATE( idesc_ip )
  DEALLOCATE( rank_ip )
  call stop_clock('protwfck')
  !call print_clock('protwfck')
  !call print_clock('protwfck:hpsi')
  !call print_clock('protwfck:hc')
  !call print_clock('protwfck:diag')
  !call print_clock('protwfck:evc')
  !
  RETURN
  !
  !
CONTAINS
  !
  SUBROUTINE desc_init( nsiz, nx, la_proc, idesc, idesc_ip )
     !
     INTEGER, INTENT(IN)  :: nsiz
     INTEGER, INTENT(OUT) :: nx
     LOGICAL, INTENT(OUT) :: la_proc
     INTEGER, INTENT(OUT) :: idesc(:)
     INTEGER, INTENT(OUT), ALLOCATABLE :: idesc_ip(:,:,:)
     INTEGER :: i, j, rank
     INTEGER :: coor_ip( 2 )
     !
     CALL laxlib_getval( np_ortho = np_ortho, me_ortho = me_ortho, &
          ortho_comm = ortho_comm, leg_ortho = leg_ortho, &
          ortho_comm_id = ortho_comm_id, ortho_parent_comm = ortho_parent_comm,&
          ortho_cntx = ortho_cntx, do_distr_diag_inside_bgrp = do_distr_diag_inside_bgrp )
     !
     ALLOCATE( idesc_ip( LAX_DESC_SIZE, np_ortho(1), np_ortho(2) ) )
     ALLOCATE( rank_ip( np_ortho(1), np_ortho(2) ) )
     CALL laxlib_init_desc( idesc, nsiz, nsiz, np_ortho, me_ortho, ortho_comm,&
          ortho_cntx, ortho_comm_id )
     !
     nx = idesc(LAX_DESC_NRCX)
     !
     DO j = 0, idesc(LAX_DESC_NPC) - 1
        DO i = 0, idesc(LAX_DESC_NPR) - 1
           coor_ip( 1 ) = i
           coor_ip( 2 ) = j
           CALL laxlib_init_desc( idesc_ip(:,i+1,j+1), idesc(LAX_DESC_N), &
                idesc(LAX_DESC_NX), np_ortho, coor_ip, ortho_comm, ortho_cntx, 1 )
           CALL GRID2D_RANK( 'R', idesc(LAX_DESC_NPR), idesc(LAX_DESC_NPC), i, j, rank )
           rank_ip( i+1, j+1 ) = rank * leg_ortho
        END DO
     END DO
     !
     la_proc = .FALSE.
     IF( idesc(LAX_DESC_ACTIVE_NODE) > 0 ) la_proc = .TRUE.
     !
     RETURN
  END SUBROUTINE desc_init
  !
  !
  SUBROUTINE compute_distmat( dm, v, w )
     !
     !  This subroutine compute <vi|wj> and store the
     !  result in distributed matrix dm 
     !
     INTEGER :: ipc, ipr
     INTEGER :: nr, nc, ir, ic, root
     COMPLEX(DP), INTENT(OUT) :: dm( :, : )
     COMPLEX(DP) :: v(:,:), w(:,:)
     COMPLEX(DP), ALLOCATABLE :: work( :, : )
     !
     ALLOCATE( work( nx, nx ) )
     !
     work = ( 0.0_DP, 0.0_DP )
     !
     DO ipc = 1, idesc(LAX_DESC_NPC) !  loop on column procs 
        !
        nc = idesc_ip( LAX_DESC_NC, 1, ipc )
        ic = idesc_ip( LAX_DESC_IC, 1, ipc )
        !
        DO ipr = 1, ipc ! desc%npr ! ipc ! use symmetry for the loop on row procs
           !
           nr = idesc_ip( LAX_DESC_NR, ipr, ipc )
           ir = idesc_ip( LAX_DESC_IR, ipr, ipc )
           !
           !  rank of the processor for which this block (ipr,ipc) is destinated
           !
           root = rank_ip( ipr, ipc )

           ! use blas subs. on the matrix block

           CALL ZGEMM( 'C', 'N', nr, nc, kdim, ( 1.D0, 0.D0 ),  v(1,ir), kdmx, w(1,ic), kdmx, ( 0.D0, 0.D0 ), work, nx )

           ! accumulate result on dm of root proc.
           CALL mp_root_sum( work, dm, root, ortho_parent_comm )

        END DO
        !
     END DO
     if (ortho_parent_comm.ne.intra_bgrp_comm .and. nbgrp > 1) dm = dm/nbgrp
     !
     CALL laxlib_zsqmher( nstart, dm, nx, idesc )
     !
     DEALLOCATE( work )
     !
     RETURN
  END SUBROUTINE compute_distmat


  SUBROUTINE refresh_evc( )
     !
     INTEGER :: ipc, ipr
     INTEGER :: nr, nc, ir, ic, root
     COMPLEX(DP), ALLOCATABLE :: vtmp( :, : )
     COMPLEX(DP) :: beta

     ALLOCATE( vtmp( nx, nx ) )
     !
     DO ipc = 1, idesc(LAX_DESC_NPC)
        !
        nc = idesc_ip( LAX_DESC_NC, 1, ipc )
        ic = idesc_ip( LAX_DESC_IC, 1, ipc )
        !
        IF( ic <= nbnd ) THEN
           !
           nc = min( nc, nbnd - ic + 1 )
           !
           beta = ( 0.D0, 0.D0 )

           DO ipr = 1, idesc(LAX_DESC_NPR)
              !
              nr = idesc_ip( LAX_DESC_NR, ipr, ipc )
              ir = idesc_ip( LAX_DESC_IR, ipr, ipc )
              !
              root = rank_ip( ipr, ipc )

              IF( ipr-1 == idesc(LAX_DESC_MYR) .AND. ipc-1 == idesc(LAX_DESC_MYC) .AND. la_proc ) THEN
                 !
                 !  this proc sends his block
                 ! 
                 CALL mp_bcast( vc(:,1:nc), root, ortho_parent_comm )
                 CALL ZGEMM( 'N', 'N', kdim, nc, nr, ( 1.D0, 0.D0 ), psi(1,ir), kdmx, vc, nx, beta, aux(1,ic), kdmx )
              ELSE
                 !
                 !  all other procs receive
                 ! 
                 CALL mp_bcast( vtmp(:,1:nc), root, ortho_parent_comm )
                 CALL ZGEMM( 'N', 'N', kdim, nc, nr, ( 1.D0, 0.D0 ), psi(1,ir), kdmx, vtmp, nx, beta, aux(1,ic), kdmx )
              END IF
              ! 

              beta = ( 1.D0, 0.D0 )

           END DO
           !
        END IF
        !
     END DO
     !
     DEALLOCATE( vtmp )

     RETURN
  END SUBROUTINE refresh_evc

END SUBROUTINE protate_wfc_k
