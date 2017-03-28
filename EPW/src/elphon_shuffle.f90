  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino  
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                        
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
  SUBROUTINE elphon_shuffle ( iq_irr, nqc_irr, iq, gmapsym, eigv, isym, xq0, timerev )
  !-----------------------------------------------------------------------
  !!
  !! Electron-phonon calculation from data saved in fildvscf
  !! Shuffle2 mode (shuffle on electrons + load all phonon q's)
  !!
  !! no ultrasoft yet
  !!
  !! RM - Nov/Dec 2014
  !! Imported the noncolinear case implemented by xlzhang
  !!
  !
  !-----------------------------------------------------------------------
  !
  USE mp,        ONLY : mp_barrier, mp_sum
  USE mp_global, ONLY : my_pool_id, nproc_pool,npool,kunit,&
                          inter_pool_comm
  USE mp_images, ONLY : nproc_image
  USE pwcom,     ONLY : nbnd, nks, nkstot
  USE gvect,     ONLY : ngm
  USE gvecs,     ONLY : doublegrid
  USE kinds,     ONLY : DP
  USE modes,     ONLY : nmodes, nirr, npert, u
  USE elph2,     ONLY : epmatq, el_ph_mat
  USE constants_epw, ONLY : czero, cone
  USE fft_base,  ONLY : dfftp, dffts
  USE noncollin_module,     ONLY : nspin_mag
!  USE noncollin_module,     ONLY : noncolin
  !
  implicit none
  !
  integer :: irr, imode0, ipert, is, iq, iq_irr, nqc_irr
  ! counter on the representations
  ! counter on the modes
  ! the change of Vscf due to perturbations
  ! the current qpoint in the uniform grid
  ! the current ireducible qpoint
  ! the total number of irreducible qpoints in the list
  complex(kind=DP), POINTER :: dvscfin(:,:,:), dvscfins (:,:,:)
  logical :: timerev
  !  true if we are using time reversal
  !
  integer :: tmp_pool_id, ik0, ik, ibnd, jbnd
  integer :: iks, nkl, nkr
  integer :: gmapsym ( ngm, 48 ), isym
  ! the correspondence G-->S(G)
  ! the symmetry which generates the current q in the star
  complex(kind=DP) :: eigv (ngm, 48)
  ! e^{ iGv} for 1...nsym ( v the fractional translation)
  real(kind=DP) :: xq0(3)
  ! the first q in the star (cartesian)
  !
  CALL start_clock ('elphon')
  !
  ik0 = 0
  tmp_pool_id = 0
  !
  npool =  nproc_image / nproc_pool
  IF (npool.gt.1) THEN
    !
    ! number of kpoint blocks, kpoints per pool and reminder
    kunit = 1 
    nkl   = kunit * ( nkstot / npool )
    nkr   = ( nkstot - nkl * npool ) / kunit
    ! the reminder goes to the first nkr pools
    IF ( my_pool_id < nkr ) nkl = nkl + kunit
    !
    iks = nkl * my_pool_id + 1
    IF ( my_pool_id >= nkr ) iks = iks + nkr * kunit
    !
    !  the index of the first k point block in this pool - 1
    !  (I will need the index of ik, not ikk)
    !
    ik0 = ( iks - 1 ) / kunit
    !
  ENDIF
  !
  ! read Delta Vscf and calculate electron-phonon coefficients
  !
  imode0 = 0
  DO irr = 1, nirr
     ALLOCATE (dvscfin ( dfftp%nnr , nspin_mag , npert(irr)) )
!DBSP
!     if (noncolin) then
!       ALLOCATE (dvscfin ( dfftp%nnr , nspin_mag , npert(irr)) )
!     endif
!END
     !
     !   read the <prefix>.dvscf_q[iq] files
     !
     dvscfin = (0.d0,0.d0)
     IF (my_pool_id.eq.0) THEN
        DO ipert = 1, npert (irr)
           CALL readdvscf ( dvscfin(1,1,ipert), imode0 + ipert, iq_irr, nqc_irr )
        ENDDO
     ENDIF
     CALL mp_sum(dvscfin,inter_pool_comm)
     !
     !
     IF (doublegrid) THEN
        ALLOCATE (dvscfins ( dffts%nnr , nspin_mag , npert(irr)) )
        DO is = 1, nspin_mag
           DO ipert = 1, npert(irr)
              CALL cinterpolate (dvscfin(1,is,ipert),dvscfins(1,is,ipert),-1)
           ENDDO 
        ENDDO
     ELSE
        dvscfins => dvscfin
     ENDIF
     !
     CALL elphel2_shuffle (npert(irr), imode0, dvscfins, gmapsym, eigv, isym, xq0, timerev)
     !
     imode0 = imode0 + npert (irr)
     IF (doublegrid) DEALLOCATE (dvscfins)
     DEALLOCATE (dvscfin)
  ENDDO
  !
  CALL mp_barrier(inter_pool_comm)
  !
  !  the output e-p matrix in the pattern representation
  !  must be transformed in the cartesian basis
  !  epmat_{CART} = conjg ( U ) * epmat_{PATTERN}
  !
  !  note it is not U^\dagger ! Have a look to symdyn_munu.f90 
  !  for comparison
  !
  DO ibnd = 1, nbnd
    DO jbnd = 1, nbnd
      DO ik = 1, nks
        ! 
        ! Here is where we calculate epmatq, it appears to be
        ! epmatq = cone * conjug(u) * el_ph_mat + czero  
        IF ( timerev ) THEN
          CALL zgemv ('n', nmodes, nmodes, cone, u , nmodes, &
            el_ph_mat (ibnd,jbnd,ik,:), 1, czero, epmatq (ibnd,jbnd,ik,:,iq), 1)
        ELSE
        CALL zgemv ('n', nmodes, nmodes, cone, CONJG ( u ), nmodes, &
          el_ph_mat (ibnd,jbnd,ik,:), 1, czero, epmatq (ibnd,jbnd,ik,:,iq), 1 )
        ENDIF 
        !
      ENDDO
    ENDDO
  ENDDO
  !DBSP
  !write(*,*)'epmatq(:,:,215,:,iq)**2',SUM((REAL(REAL(epmatq(:,:,215,:,iq))))**2)+&
  !        SUM((REAL(AIMAG(epmatq(:,:,215,:,iq))))**2)
  !END
  !
  CALL stop_clock ('elphon')
  !
  END SUBROUTINE elphon_shuffle
