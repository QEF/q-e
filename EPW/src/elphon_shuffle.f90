  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino  
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                        
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
  SUBROUTINE elphon_shuffle( iq_irr, nqc_irr, iq, gmapsym, eigv, isym, xq0, timerev )
  !-----------------------------------------------------------------------
  !!
  !! Electron-phonon calculation from data saved in fildvscf
  !! Shuffle2 mode (shuffle on electrons + load all phonon q's)
  !!
  !! RM - Nov/Dec 2014
  !! Imported the noncolinear case implemented by xlzhang
  !!
  !! Roxana Margine - Jan 2019: Updated based on QE 6.3 for US
  !!
  !-----------------------------------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE mp,        ONLY : mp_barrier, mp_sum
  USE mp_pools,  ONLY : my_pool_id, npool, inter_pool_comm
  USE ions_base, ONLY : nat
  USE pwcom,     ONLY : nbnd, nks, nkstot
  USE gvect,     ONLY : ngm
  USE gvecs,     ONLY : doublegrid
  USE modes,     ONLY : nmodes, nirr, npert, u
  USE elph2,     ONLY : epmatq, el_ph_mat
  USE lrus,      ONLY : int3, int3_nc
  USE uspp,      ONLY : okvan
  USE lsda_mod,  ONLY : nspin
  USE fft_base,  ONLY : dfftp, dffts
  USE uspp_param,       ONLY : nhm
  USE constants_epw,    ONLY : czero, cone
  USE fft_interfaces,   ONLY : fft_interpolate
  USE noncollin_module, ONLY : nspin_mag, noncolin
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: iq_irr
  !! Current ireducible q-point
  INTEGER, INTENT(in) :: nqc_irr
  !! Total number of irreducible q-points in the list
  INTEGER, INTENT(in) :: iq
  !! Current q-point in the star of iq_irr q-point
  INTEGER, INTENT(in) :: gmapsym(ngm,48)
  !! Correspondence G-->S(G)
  INTEGER, INTENT(in) :: isym
  !! The symmetry which generates the current q in the star
  REAL(DP), INTENT(in) :: xq0(3)
  !! The first q-point in the star (cartesian coords.)
  COMPLEX(DP), INTENT(in) :: eigv(ngm,48)
  !! e^{iGv} for 1...nsym (v the fractional translation)
  LOGICAL, INTENT(in) :: timerev
  !!  true if we are using time reversal
  !
  ! Local variables
  !
  INTEGER :: irr 
  !! Counter on representations
  INTEGER :: imode0
  !! Counter on modes
  INTEGER :: ipert
  !! Change of Vscf due to perturbations
  INTEGER :: npe
  !! Number of perturbations for irr representation
  INTEGER :: is
  !! Counter on spin
  INTEGER :: ik
  !! Counter on k-points in the pool
  INTEGER :: ibnd
  !! Counter on bands
  INTEGER :: jbnd
  !! Counter on bands
  !
  COMPLEX(kind=DP), POINTER :: dvscfin(:,:,:)
  !! Change of the scf potential 
  COMPLEX(kind=DP), POINTER :: dvscfins(:,:,:)
  !! Change of the scf potential (smooth)
  !
  CALL start_clock('elphon_shuffle')
  !
  ! read Delta Vscf and calculate electron-phonon coefficients
  !
  imode0 = 0
  DO irr = 1, nirr
     npe = npert(irr)
     ALLOCATE( dvscfin(dfftp%nnr, nspin_mag, npe) )
     IF (okvan) THEN
        ALLOCATE( int3(nhm, nhm, nat, nspin_mag, npe) )
        IF (noncolin) ALLOCATE( int3_nc(nhm, nhm, nat, nspin, npe) )
     ENDIF
     !
     !   read the <prefix>.dvscf_q[iq] files
     !
     dvscfin = czero
     IF ( my_pool_id.eq.0 ) THEN
        DO ipert = 1, npe
           CALL readdvscf( dvscfin(1,1,ipert), imode0 + ipert, iq_irr, nqc_irr )
        ENDDO
     ENDIF
     CALL mp_sum(dvscfin,inter_pool_comm)
     !
     IF (doublegrid) THEN
       ALLOCATE( dvscfins(dffts%nnr, nspin_mag, npe) )
       DO is = 1, nspin_mag
         DO ipert = 1, npe
           CALL fft_interpolate(dfftp, dvscfin(:,is,ipert), dffts, dvscfins(:,is,ipert))
         ENDDO 
       ENDDO
     ELSE
       dvscfins => dvscfin
     ENDIF
     !
     CALL newdq2( dvscfin, npe, xq0, timerev )
     CALL elphel2_shuffle( npe, imode0, dvscfins, gmapsym, eigv, isym, xq0, timerev )
     !
     imode0 = imode0 + npe
     IF (doublegrid) DEALLOCATE(dvscfins)
     DEALLOCATE(dvscfin)
     IF (okvan) THEN
        DEALLOCATE(int3)
        IF (noncolin) DEALLOCATE(int3_nc)
     ENDIF
  ENDDO
  !
  CALL mp_barrier(inter_pool_comm)
  !
  !  the output e-p matrix in the pattern representation
  !  must be transformed in the cartesian basis
  !  epmat_{CART} = conjg ( U ) * epmat_{PATTERN}
  !
  !  note it is not U^\dagger but u_pattern! 
  !  Have a look to symdyn_munu.f90 for comparison
  !
  DO ibnd = 1, nbnd
    DO jbnd = 1, nbnd
      DO ik = 1, nks
        ! 
        ! Here is where we calculate epmatq, it appears to be
        ! epmatq = cone * conjug(u) * el_ph_mat + czero  
        IF ( timerev ) THEN
          CALL zgemv( 'n', nmodes, nmodes, cone, u, nmodes, &
            el_ph_mat(ibnd,jbnd,ik,:), 1, czero, epmatq(ibnd,jbnd,ik,:,iq), 1 )
        ELSE
          CALL zgemv( 'n', nmodes, nmodes, cone, conjg(u), nmodes, &
            el_ph_mat(ibnd,jbnd,ik,:), 1, czero, epmatq(ibnd,jbnd,ik,:,iq), 1 )
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
  CALL stop_clock('elphon_shuffle')
  !
  END SUBROUTINE elphon_shuffle
