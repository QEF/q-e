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
  !     This initialization is done nqc_irr times from elphon_shuffle_wrap
  !     not all of the following code is necessary.  More adaptation from
  !     phq_init is needed   
  !
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp, tau
  USE becmod,               ONLY : calbec
  USE phus,                 ONLY : alphap
  USE lrus,                 ONLY : becp1
  USE uspp,                 ONLY : vkb
  USE pwcom,                ONLY : npwx, nbnd, nks, lsda, current_spin,&
                                   g2kin, isk, xk, strf
  USE constants,            ONLY : tpi
  USE cell_base,            ONLY : tpiba2, tpiba, bg, omega
  USE klist,                ONLY : ngk, igk_k, nkstot
  USE constants_epw,        ONLY : zero
  USE gvecw,                ONLY : ecutwfc
  USE gvect,                ONLY : eigts1, eigts2, eigts3, g, ngm
  USE atom,                 ONLY : msh, rgrid
  USE wavefunctions_module, ONLY : evc
  USE kinds,                ONLY : DP
  USE noncollin_module,     ONLY : noncolin, npol
  USE uspp_param,           ONLY : upf
  USE phcom,                ONLY : lrwfc, iuwfc, vlocq
  USE qpoint,               ONLY : xq, eigqts, npwq
  USE nlcc_ph,              ONLY : drc                           
  USE uspp,                 ONLY : nlcc_any
  USE fft_base,             ONLY : dfftp
  USE elph2,                ONLY : igk_k_all, ngk_all
  USE mp,                   ONLY : mp_barrier
  USE mp_global,            ONLY : inter_pool_comm
  !
  IMPLICIT NONE
  !
  ! ... local variables
  !
  INTEGER :: nt, ik, ipol, ibnd, na, ig
    ! counter on atom types
    ! counter on k points
    ! counter on k+q points
    ! counter on polarizations
    ! counter on bands
    ! index for wavefunctions at k
    ! counter on atoms
    ! counter on G vectors
  REAL(DP) :: arg
    ! the argument of the phase
  COMPLEX(DP), ALLOCATABLE :: aux1(:,:)
    ! used to compute alphap
  logical :: first_run
  !
  !
  !
  CALL start_clock( 'epw_init' )
  !
  ALLOCATE( aux1( npwx*npol, nbnd ) )    
  !
  ! ... initialize structure factor array
  !
  CALL struc_fact( nat, tau, ntyp, ityp, ngm, g, bg, dfftp%nr1,dfftp%nr2, dfftp%nr3, &
                   strf, eigts1, eigts2, eigts3 )
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
  ! compute rhocore for each atomic-type if needed for nlcc
  !
  IF ( nlcc_any ) CALL set_drhoc( xq, drc )
  !
  ! the fourier components of the local potential for each |G|
  !
  CALL init_vloc()
  !
  ! the fourier components of the local potential at q+G
  !
  vlocq(:,:) = 0.D0
  !
  DO nt = 1, ntyp
     !
     CALL setlocq( xq, rgrid(nt)%mesh, msh(nt), rgrid(nt)%rab, rgrid(nt)%r,&
                   upf(nt)%vloc(1), upf(nt)%zp, tpiba2, ngm, g, omega, &
                   vlocq(1,nt) )
     !
  END DO
  !
  ! the parameters defining the pseudopotential
  !
  ! we compute the denominators of the KB types, or the
  ! parameters which define the non-local pseudopotential and
  ! which are independent of the k point for the US case
  !
  CALL init_us_1()
  !
  DO ik = 1, nks
     !
     !
     IF ( lsda ) current_spin = isk( ik )
     !
     ! g2kin is used here as work space
     !
     CALL gk_sort( xk(1,ik), ngm, g, ( ecutwfc / tpiba2 ), ngk(ik), igk_k(1,ik), g2kin )
     ! 
     ! if there is only one k-point evc, evq, npw, igk stay in memory
     !
     npwq = ngk(ik)
     !
     ! The functions vkb(k+G)
     !
     CALL init_us_2( ngk(ik), igk_k(1,ik), xk(1,ik), vkb )
     !
     ! ... read the wavefunctions at k
     !
     CALL davcio( evc, lrwfc, iuwfc, ik, -1 )
     !
     ! we compute the becp terms which are used in the rest of
     ! the code
     !
     IF (noncolin) THEN
        CALL calbec (ngk(ik), vkb, evc, becp1(ik)%nc(:,:,:) )
     ELSE
        CALL calbec (ngk(ik), vkb, evc, becp1(ik)%k(:,:))
     ENDIF
     !
     ! we compute the derivative of the becp term with respect to an
     !   atomic displacement
     !
     DO ipol = 1, 3
        aux1=(0.d0,0.d0)
        DO ibnd = 1, nbnd
           DO ig = 1, ngk(ik)
              aux1(ig,ibnd) = evc(ig,ibnd) * tpiba * ( 0.D0, 1.D0 ) * & 
                              ( xk(ipol,ik) + g(ipol,igk_k(ig,ik)) )
           ENDDO
           IF (noncolin) THEN
              DO ig = 1, ngk(ik)
                 aux1(ig+npwx,ibnd) = evc(ig+npwx,ibnd)*tpiba*(0.D0,1.D0)*& 
                           ( xk(ipol,ik) + g(ipol,igk_k(ig,ik)) )
              ENDDO
           ENDIF
        ENDDO
        IF (noncolin) THEN
           CALL calbec (ngk(ik), vkb, aux1, alphap(ipol,ik)%nc(:,:,:) )
        ELSE
           CALL calbec (ngk(ik), vkb, aux1, alphap(ipol,ik)%k(:,:) )
        ENDIF
     ENDDO
     !
     !
  ENDDO
  !
  IF(.not. ALLOCATED(igk_k_all) ) ALLOCATE(igk_k_all( npwx, nkstot))
  IF(.not. ALLOCATED(ngk_all) ) ALLOCATE(ngk_all(nkstot))
  !
#if defined(__MPI)
  !
  CALL poolgather_int (npwx, nkstot, nks, igk_k(:,1:nks), igk_k_all ) 
  CALL poolgather_int1 (nkstot, nks, ngk(1:nks), ngk_all ) 
  CALL mp_barrier(inter_pool_comm)
  !
#else
  !
  igk_k_all = igk_k
  ngk_all = ngk
  !
#endif
  !
  DEALLOCATE( aux1 )
  !
  IF (.not.first_run) CALL dvanqq2()
  !
  !
  CALL stop_clock( 'epw_init' )
  !
  END SUBROUTINE epw_init
