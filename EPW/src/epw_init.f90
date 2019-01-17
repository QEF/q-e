  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  ! Adapted from the code PH/phq_init - Quantum-ESPRESSO group                 
  !--------------------------------------------------------------------
  SUBROUTINE epw_init(first_run)
  !----------------------------------------------------------------------------
  !
  !!     This initialization is done nqc_irr times from elphon_shuffle_wrap
  !!     not all of the following code is necessary.  More adaptation from
  !!     phq_init is needed   
  !!
  !!     Roxana Margine - Dec 2018: Updated based on QE 6.3
  !!
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp, tau
  USE becmod,               ONLY : calbec
  USE phus,                 ONLY : alphap
  USE lrus,                 ONLY : becp1
  USE uspp,                 ONLY : vkb
  USE pwcom,                ONLY : npwx, nbnd, nks, lsda, current_spin, &
                                   isk, xk
  USE constants,            ONLY : tpi
  USE constants_epw,        ONLY : zero, czero, cone
  USE cell_base,            ONLY : tpiba2, tpiba, bg, omega
  USE klist,                ONLY : ngk, igk_k, nkstot
  USE gvect,                ONLY : g, ngm
  USE atom,                 ONLY : msh, rgrid
  USE wavefunctions,        ONLY : evc
  USE noncollin_module,     ONLY : noncolin, npol
  USE uspp_param,           ONLY : upf
  USE m_gth,                ONLY : setlocq_gth
  USE units_lr,             ONLY : lrwfc, iuwfc
  USE phcom,                ONLY : vlocq
  USE qpoint,               ONLY : xq, eigqts
  USE nlcc_ph,              ONLY : drc                           
  USE uspp,                 ONLY : nlcc_any
  USE elph2,                ONLY : igk_k_all, ngk_all
  USE mp,                   ONLY : mp_barrier
  USE mp_global,            ONLY : inter_pool_comm
  !
  IMPLICIT NONE
  !
  LOGICAL :: first_run
  !
  ! ... Local variables
  !
  INTEGER :: nt
  !! counter on atom types
  INTEGER :: ik
  !! counter on k points
  INTEGER :: ipol
  !! counter on polarizations
  INTEGER :: ibnd
  !! counter on bands
  INTEGER :: na
  !! counter on atoms
  INTEGER :: ig
  !! counter on G vectors
  !
  REAL(DP) :: arg
  !! the argument of the phase
  !
  COMPLEX(DP), ALLOCATABLE :: aux1(:,:)
  !! used to compute alphap
  !
  !
  CALL start_clock( 'epw_init' )
  !
  DO na = 1, nat
    !
    ! xq here is the first q of the star
    arg = ( xq(1) * tau(1,na) + &
            xq(2) * tau(2,na) + &
            xq(3) * tau(3,na) ) * tpi
    !        
    eigqts(na) = CMPLX( COS( arg ), - SIN( arg ), kind=DP )
    !
  END DO
  !
  ! ... a0) compute rhocore for each atomic-type if needed for nlcc
  !
  IF ( nlcc_any ) CALL set_drhoc( xq, drc )
  !
  ! ... b) the fourier components of the local potential at q+G
  !
  vlocq(:,:) = zero
  !
  DO nt = 1, ntyp
    !
    IF (upf(nt)%is_gth) then
      CALL setlocq_gth( nt, xq, upf(nt)%zp, tpiba2, ngm, g, omega, vlocq(1,nt) )
    ELSE
      CALL setlocq( xq, rgrid(nt)%mesh, msh(nt), rgrid(nt)%rab, rgrid(nt)%r, &
                    upf(nt)%vloc(1), upf(nt)%zp, tpiba2, ngm, g, omega, &
                    vlocq(1,nt) )
    ENDIF
    !
  END DO
  !
  ALLOCATE( aux1( npwx*npol, nbnd ) )
  !
  DO ik = 1, nks
    !
    !
    IF ( lsda ) current_spin = isk( ik )
    !
    ! ... d) The functions vkb(k+G)
    !
    CALL init_us_2( ngk(ik), igk_k(1,ik), xk(1,ik), vkb )
    !
    ! ... read the wavefunctions at k
    !
    CALL davcio( evc, lrwfc, iuwfc, ik, -1 )
    !
    ! ... e) we compute the becp terms which are used in the rest of
    ! ...    the code
    !
    CALL calbec( ngk(ik), vkb, evc, becp1(ik) )
    !
    ! we compute the derivative of the becp term with respect to an
    !   atomic displacement
    !
    DO ipol = 1, 3
      aux1 = czero
      DO ibnd = 1, nbnd
        DO ig = 1, ngk(ik)
          aux1(ig,ibnd) = evc(ig,ibnd) * tpiba * cone * & 
                          ( xk(ipol,ik) + g(ipol,igk_k(ig,ik)) )
        ENDDO
        IF (noncolin) THEN
          DO ig = 1, ngk(ik)
            aux1(ig+npwx,ibnd) = evc(ig+npwx,ibnd) * tpiba *cone *& 
                      ( xk(ipol,ik) + g(ipol,igk_k(ig,ik)) )
          ENDDO
        ENDIF
      ENDDO
      CALL calbec( ngk(ik), vkb, aux1, alphap(ipol,ik) )
    ENDDO
    !
    !
  ENDDO
  !
  DEALLOCATE( aux1 )
  !
  IF(.not. ALLOCATED(igk_k_all)) ALLOCATE(igk_k_all(npwx,nkstot))
  IF(.not. ALLOCATED(ngk_all))   ALLOCATE(ngk_all(nkstot))
  !
#if defined(__MPI)
  !
  CALL poolgather_int(npwx, nkstot, nks, igk_k(:,1:nks), igk_k_all) 
  CALL poolgather_int1(nkstot, nks, ngk(1:nks), ngk_all) 
  CALL mp_barrier(inter_pool_comm)
  !
#else
  !
  igk_k_all = igk_k
  ngk_all = ngk
  !
#endif
  !
  IF (.not.first_run) CALL dvanqq2()
  !
  CALL stop_clock( 'epw_init' )
  !
  END SUBROUTINE epw_init
