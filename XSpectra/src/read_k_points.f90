    !<MCB>


    subroutine read_k_points
      USE start_k,            ONLY : nk1, nk2, nk3, k1, k2, k3
      USE io_global,          ONLY : ionode_id 
      USE noncollin_module,     ONLY : noncolin
      USE klist,              ONLY : npk, xk, wk, nks
      USE lsda_mod,           ONLY : isk, lsda, nspin
      USE constants,          ONLY : degspin
      USE parser,             ONLY : read_line
      USE cell_base,       ONLY : bg, at, celldm
      USE mp,                 ONLY : mp_bcast
      USE mp_world,           ONLY : world_comm ! not sure about this
      
      implicit none
      INTEGER               :: npool, nkl, nkr, nkbl, iks, ike
      CHARACTER(LEN=256)         :: input_line
      INTEGER i,j,k,n
      !   Define new k-point mesh
      
      CALL read_line( input_line )
      READ(input_line, *) nk1, nk2, nk3, k1, k2 ,k3
      IF ( k1 < 0 .OR. k1 > 1 .OR. &
           k2 < 0 .OR. k2 > 1 .OR. &
           k3 < 0 .OR. k3 > 1 ) CALL errore &
           ('card_kpoints', 'invalid offsets: must be 0 or 1', 1)
      IF ( nk1 <= 0 .OR. nk2 <= 0 .OR. nk3 <= 0 ) CALL errore &
           ('card_kpoints', 'invalid values for nk1, nk2, nk3', 1)
      
      CALL mp_bcast( k1,  ionode_id, world_comm )
      CALL mp_bcast( k2,  ionode_id, world_comm )
      CALL mp_bcast( k3,  ionode_id, world_comm )
      CALL mp_bcast( nk1,  ionode_id,world_comm )
      CALL mp_bcast( nk2,  ionode_id,world_comm )
      CALL mp_bcast( nk3,  ionode_id,world_comm )
      
      !     if(nsym.eq.1) then
      nks=nk1*nk2*nk3
      do i=1,nk1
         do j=1,nk2
            do k=1,nk3
               !  this is nothing but consecutive ordering
               n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
               !  xkg are the components of the complete grid in crystal axis
               xk(1,n) = DBLE(i-1)/nk1 + DBLE(k1)/2/nk1
               xk(2,n) = DBLE(j-1)/nk2 + DBLE(k2)/2/nk2
               xk(3,n) = DBLE(k-1)/nk3 + DBLE(k3)/2/nk3
            end do
         end do
      end do
      wk(1:nks)=1.0/DBLE(nks)
      call cryst_to_cart(nks,xk,bg,1)
      !     else
      !       CALL kpoint_grid( nsym, s, bg, npk, k1, k2, k3, &
      !          nk1, nk2, nk3, nks, xk, wk )
      !     endif
      
      if(lsda) then
         CALL set_kup_and_kdw( xk, wk, isk, nks, npk )
      ELSE IF ( noncolin ) THEN
         call errore('define_and_distribute_k_points', &
              'noncolinear not implemented', 1 )
      else 
         isk(1:nks)=1
         wk(1:nks)=2.d0/dfloat(nks)
      endif
      
      
    end subroutine read_k_points

    !</MCB>


