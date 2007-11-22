!
! Copyright (C) 2001-2007 Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine rotate_wfc_nc &
     (npwx, npw, nstart, nbnd, psi, npol, overlap, evc, e)
  !----------------------------------------------------------------------
  !
  !   Serial version of rotate_wfc for noncollinear caluculations
  !
#include "f_defs.h"
  USE kinds
  implicit none
  !
  integer :: npw, npwx, nstart, nbnd, npol
  ! dimension of the matrix to be diagonalized
  ! leading dimension of matrix psi, as declared in the calling pgm unit
  ! input number of states
  ! output number of states
  ! number of spin polarizations
  logical :: overlap
  ! if .false. : S|psi> not needed

  complex(DP) :: psi (npwx, npol, nstart), evc (npwx, npol, nbnd)

  ! input and output eigenvectors (may overlap)

  real(DP) :: e (nbnd)
  ! eigenvalues

  ! auxiliary variables:
  integer :: i, j, idx
  ! counters
  complex(DP), allocatable :: hpsi (:,:,:), spsi (:,:,:), &
       hc (:,:), sc (:,:), vc (:,:) 
  real(DP), allocatable :: en (:)
  COMPLEX(DP), EXTERNAL :: ZDOTC
  !
  allocate (hpsi(  npwx ,npol, nstart))    
  if (overlap) allocate (spsi(  npwx , npol, nstart))    
  allocate (hc( nstart , nstart))    
  allocate (sc( nstart , nstart))    
  allocate (vc( nstart , nstart))    
  allocate (en( nstart ))
  !!!!!!!!!!!!!!!!! 
  hpsi=(0.d0,0.d0)
  if (overlap) spsi=(0.d0,0.d0)
  hc=(0.d0,0.d0)
  sc=(0.d0,0.d0)
  vc=(0.d0,0.d0)
  en=0.d0
  !!!!!!!!!!!!!!!!
  !
  ! Set up the Hamiltonian and Overlap matrix
  !
  call h_psi_nc (npwx, npw, nstart, psi(1,1,1), hpsi(1,1,1))
  if (overlap) call s_psi_nc (npwx, npw, nstart, psi(1,1,1), spsi(1,1,1))

  if (npol.eq.1) then
     call ZGEMM ('c','n',nstart,nstart,npw,(1.d0,0.d0),psi(1,1,1),npwx, &
       hpsi(1,1,1), npwx, (0.d0, 0.d0) , hc, nstart)
  else
     ! FIXME: this should be replaced by a call to ZGEMM
     do j=1,nstart
        do i=1,nstart
           hc(i,j) = ZDOTC(npw,psi(1,1,i),1,hpsi(1,1,j),1) + &
                     ZDOTC(npw,psi(1,2,i),1,hpsi(1,2,j),1)
        enddo
     enddo
  endif
#ifdef __PARA
  call reduce (2 * nstart * nstart, hc)
#endif
  if (overlap) then
     if (npol.eq.1) then
        call ZGEMM ('c','n',nstart,nstart,npw,(1.d0,0.d0),psi(1,1,1),npwx, &
             spsi(1,1,1),npwx,(0.d0,0.d0),sc,nstart)
     else
        ! FIXME: this should be replaced by a call to ZGEMM
        do j=1,nstart
           do i=1,nstart
              sc(i,j) = ZDOTC(npw,psi(:,1,i),1,spsi(:,1,j),1) + &
                        ZDOTC(npw,psi(:,2,i),1,spsi(:,2,j),1)
              ! WRITE( stdout,*) ' sc ', i, j, sc(i,j)
           enddo
        enddo
     end if
  else
     if (npol.eq.1) then
        call ZGEMM ('c','n',nstart,nstart,npw,(1.d0,0.d0),psi(1,1,1),npwx, &
             psi(1,1,1),npwx,(0.d0,0.d0),sc,nstart)
     else
        ! FIXME: this should be replaced by a call to ZGEMM
        do j=1,nstart
           do i=1,nstart
              sc(i,j) = ZDOTC(npw,psi(:,1,i),1,psi(:,1,j),1) + &
                   ZDOTC(npw,psi(:,2,i),1,psi(:,2,j),1)
              ! WRITE( stdout,*) ' sc ', i, j, sc(i,j)
           enddo
        enddo
     endif
  end if
  
#ifdef __PARA
  call reduce (2 * nstart * nstart, sc)
#endif
  !
  ! Diagonalize
  !
  call cdiaghg (nstart, nbnd, hc, sc, nstart, en, vc)
  !
  e (:) = en(1:nbnd)
  !
  !   update the basis set
  !  
  if (npol.eq.1) then
     call ZGEMM ('n', 'n', npw, nbnd, nstart, (1.d0, 0.d0) , psi(1,1,1), npwx, &
                 vc, nstart, (0.d0, 0.d0) , hpsi(1,1,1), npwx)
  else
     do j=1,nbnd
        do i=1,npw
           hpsi(i,1,j) = (0.d0,0.d0)
           hpsi(i,2,j) = (0.d0,0.d0)
           do idx=1,nstart
              hpsi(i,1,j) = hpsi(i,1,j) + psi(i,1,idx)*vc(idx,j)
              hpsi(i,2,j) = hpsi(i,2,j) + psi(i,2,idx)*vc(idx,j)
           enddo
        enddo
     enddo
  endif
  evc(:, :,  :) = hpsi(:,:, 1:nbnd)

  deallocate (en)
  deallocate (vc)
  deallocate (sc)
  deallocate (hc)
  if (overlap) deallocate (spsi)
  deallocate (hpsi)

  return
end subroutine rotate_wfc_nc

!
!----------------------------------------------------------------------
subroutine protate_wfc_nc &
     (npwx, npw, nstart, nbnd, psi, npol, overlap, evc, e)
  !----------------------------------------------------------------------
  !
  ! ... Parallel version of rotate_wfc for noncollinear calculations
  ! ... Subroutine with distributed matrices, written by Carlo Cavazzoni
  !
#include "f_defs.h"
  USE kinds
  USE mp_global,        ONLY : npool, nproc_pool, me_pool, root_pool, &
                               intra_pool_comm, init_ortho_group, me_image, &
                               ortho_comm, np_ortho, me_ortho, ortho_comm_id,&
                               leg_ortho
  USE descriptors,      ONLY : descla_siz_ , descla_init , lambda_node_ , &
                               la_nx_ , la_nrl_ , la_n_ , &
                               ilac_ , ilar_ , nlar_ , nlac_ , la_npc_ , &
                               la_npr_ , la_me_ , la_comm_ , &
                               la_myr_ , la_myc_ , nlax_
  USE parallel_toolkit, ONLY : zsqmred, zsqmher, zsqmdst
  USE mp,               ONLY : mp_bcast, mp_root_sum, mp_sum, mp_barrier, mp_end

  implicit none
  !
  integer :: npw, npwx, nstart, nbnd, npol
  ! dimension of the matrix to be diagonalized
  ! leading dimension of matrix psi, as declared in the calling pgm unit
  ! input number of states
  ! output number of states
  ! number of spin polarizations
  logical :: overlap
  ! if .false. : S|psi> not needed

  complex(DP) :: psi (npwx, npol, nstart), evc (npwx, npol, nbnd)
  ! input and output eigenvectors (may overlap)

  real(DP) :: e (nbnd)
  ! eigenvalues

  ! auxiliary variables:
  integer :: i, j, idx
  ! counters
  complex(DP), allocatable :: hpsi (:,:,:), spsi (:,:,:), &
       hc (:,:), sc (:,:), vc (:,:) 
  real(DP), allocatable :: en (:)
  !
  INTEGER :: desc( descla_siz_ )
    ! matrix distribution descriptors
  INTEGER :: nx
    ! maximum local block dimension
  LOGICAL :: la_proc
    ! flag to distinguish procs involved in linear algebra
  INTEGER, ALLOCATABLE :: desc_ip( :, :, : )
  INTEGER, ALLOCATABLE :: rank_ip( :, : )

  COMPLEX(DP), external :: ZDOTC
  !
  ALLOCATE( desc_ip( descla_siz_ , np_ortho(1), np_ortho(2) ) )
  ALLOCATE( rank_ip( np_ortho(1), np_ortho(2) ) )
  !
  CALL desc_init( nstart, desc, desc_ip )
  !
  allocate (hpsi(  npwx ,npol, nstart))    
  if (overlap) allocate (spsi(  npwx , npol, nstart))    
  allocate (hc( nx , nx))    
  allocate (sc( nx , nx))    
  allocate (vc( nx , nx))    
  allocate (en( nstart ))
  !!!!!!!!!!!!!!!!! 
  hpsi=(0.d0,0.d0)
  if (overlap) spsi=(0.d0,0.d0)
  hc=(0.d0,0.d0)
  sc=(0.d0,0.d0)
  vc=(0.d0,0.d0)
  en=0.d0
  !!!!!!!!!!!!!!!!
  !
  ! Set up the Hamiltonian and Overlap matrix
  !
  call h_psi_nc (npwx, npw, nstart, psi(1,1,1), hpsi(1,1,1))
  if (overlap) call s_psi_nc (npwx, npw, nstart, psi(1,1,1), spsi(1,1,1))

  call compute_distmat( hc, psi, hpsi )

  if (overlap) then

     call compute_distmat( sc, psi, spsi )

  else

     call compute_distmat( sc, psi, psi )

  end if
  !
  ! Diagonalize
  !
  CALL pcdiaghg( nstart, hc, sc, nx, en, vc, desc )
  !
  e (:) = en(1:nbnd)
  !
  !   update the basis set
  !  
  call refresh_evc()

  evc(:, :,  :) = hpsi(:,:, 1:nbnd)

  deallocate( en )
  deallocate( vc )
  deallocate( sc )
  deallocate( hc )
  if (overlap) deallocate (spsi)
  deallocate (hpsi)

  deallocate( desc_ip )
  deallocate( rank_ip )

  call stop_clock ('wfcrot')
  return

contains 

  SUBROUTINE desc_init( nsiz, desc, desc_ip )
     !
     INTEGER, INTENT(IN)  :: nsiz
     INTEGER, INTENT(OUT) :: desc(:)
     INTEGER, INTENT(OUT) :: desc_ip(:,:,:)
     INTEGER :: i, j, rank
     INTEGER :: coor_ip( 2 )
     !
     CALL descla_init( desc, nsiz, nsiz, np_ortho, me_ortho, ortho_comm, ortho_comm_id )
     !
     nx = desc( nlax_ )
     !
     DO j = 0, desc( la_npc_ ) - 1
        DO i = 0, desc( la_npr_ ) - 1
           coor_ip( 1 ) = i
           coor_ip( 2 ) = j
           CALL descla_init( desc_ip(:,i+1,j+1), desc( la_n_ ), desc( la_nx_ ), &
                             np_ortho, coor_ip, ortho_comm, 1 ) 
           CALL GRID2D_RANK( 'R', desc( la_npr_ ), desc( la_npc_ ), i, j, rank )
           rank_ip( i+1, j+1 ) = rank * leg_ortho
        END DO
     END DO
     !
     la_proc = .FALSE.
     IF( desc( lambda_node_ ) > 0 ) la_proc = .TRUE.
     !
     RETURN
  END SUBROUTINE desc_init

  !
  SUBROUTINE compute_distmat( dm, v, w )
     !
     !  This subroutine compute <vi|wj> and store the
     !  result in distributed matrix dm 
     !
     INTEGER :: ipc, ipr
     INTEGER :: nr, nc, ir, ic, root
     COMPLEX(DP), INTENT(OUT) :: dm( :, : )
     COMPLEX(DP) :: v(:,:,:), w(:,:,:)
     COMPLEX(DP), ALLOCATABLE :: work( :, : )
     !
     ALLOCATE( work( nx, nx ) )
     !
     work = ( 0.0_DP, 0.0_DP )
     !
     DO ipc = 1, desc( la_npc_ ) !  loop on column procs 
        !
        nc = desc_ip( nlac_ , 1, ipc )
        ic = desc_ip( ilac_ , 1, ipc )
        !
        DO ipr = 1, ipc ! desc( la_npr_ ) ! ipc ! use symmetry for the loop on row procs
           !
           nr = desc_ip( nlar_ , ipr, ipc )
           ir = desc_ip( ilar_ , ipr, ipc )
           !
           !  rank of the processor for which this block (ipr,ipc) is destinated
           !
           root = rank_ip( ipr, ipc )

           ! use blas subs. on the matrix block

           if (npol.eq.1) then
              call ZGEMM ('c','n',nr,nc,npw,(1.d0,0.d0),v(1,1,ir),npwx, &
                          w(1,1,ic), npwx, (0.d0, 0.d0) , work, nx )
           else
              ! FIXME: this should be replaced by a call to ZGEMM
              do j=1,nc
                 do i=1,nr
                    work(i,j) = ZDOTC(npw,v(1,1,i+ir-1),1,w(1,1,j+ic-1),1) + &
                                ZDOTC(npw,v(1,2,i+ir-1),1,w(1,2,j+ic-1),1)
                 enddo
              enddo
           endif

           ! accumulate result on dm of root proc.

           CALL mp_root_sum( work, dm, root, intra_pool_comm )

        END DO
        !
     END DO
     !
     CALL zsqmher( nstart, dm, nx, desc )
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
     DO ipc = 1, desc( la_npc_ )
        !
        nc = desc_ip( nlac_ , 1, ipc )
        ic = desc_ip( ilac_ , 1, ipc )
        !
        IF( ic <= nbnd ) THEN
           !
           nc = min( nc, nbnd - ic + 1 )
           !
           beta = ( 0.D0, 0.D0 )

           DO ipr = 1, desc( la_npr_ )
              !
              nr = desc_ip( nlar_ , ipr, ipc )
              ir = desc_ip( ilar_ , ipr, ipc )
              !
              root = rank_ip( ipr, ipc )

              IF( ipr-1 == desc( la_myr_ ) .AND. ipc-1 == desc( la_myc_ ) .AND. la_proc ) THEN
                 !
                 !  this proc sends his block
                 ! 
                 CALL mp_bcast( vc(:,1:nc), root, intra_pool_comm )
                 !
                 if (npol.eq.1) then
                    call ZGEMM ('n', 'n', npw, nc, nr, (1.d0, 0.d0) , psi(1,1,ir), npwx, &
                                vc, nx, beta, hpsi(1,1,ic), npwx)
                 else
                    do j=1,nc
                       do i=1,npw
                          IF( beta == ( 0.D0, 0.D0 ) ) THEN
                             hpsi(i,1,j+ic-1) = (0.d0,0.d0)
                             hpsi(i,2,j+ic-1) = (0.d0,0.d0)
                          ENDIF
                          do idx=1,nr
                             hpsi(i,1,j+ic-1) = hpsi(i,1,j+ic-1) + psi(i,1,idx+ir-1)*vc(idx,j)
                             hpsi(i,2,j+ic-1) = hpsi(i,2,j+ic-1) + psi(i,2,idx+ir-1)*vc(idx,j)
                          enddo
                       enddo
                    enddo
                 endif

              ELSE
                 !
                 !  all other procs receive
                 ! 
                 CALL mp_bcast( vtmp(:,1:nc), root, intra_pool_comm )
                 !
                 if (npol.eq.1) then
                    call ZGEMM ('n', 'n', npw, nc, nr, (1.d0, 0.d0) , psi(1,1,ir), npwx, &
                                vtmp, nx, beta, hpsi(1,1,ic), npwx)
                 else
                    do j=1,nc
                       do i=1,npw
                          IF( beta == ( 0.D0, 0.D0 ) ) THEN
                             hpsi(i,1,j+ic-1) = (0.d0,0.d0)
                             hpsi(i,2,j+ic-1) = (0.d0,0.d0)
                          ENDIF
                          do idx=1,nr
                             hpsi(i,1,j+ic-1) = hpsi(i,1,j+ic-1) + psi(i,1,idx+ir-1)*vtmp(idx,j)
                             hpsi(i,2,j+ic-1) = hpsi(i,2,j+ic-1) + psi(i,2,idx+ir-1)*vtmp(idx,j)
                          enddo
                       enddo
                    enddo
                 endif
                 !
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


end subroutine protate_wfc_nc

