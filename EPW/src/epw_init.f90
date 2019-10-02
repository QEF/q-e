  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  ! Adapted from the code PH/phq_init - Quantum-ESPRESSO group                 
  !----------------------------------------------------------------------------
  SUBROUTINE epw_init(first_run)
  !----------------------------------------------------------------------------
  !
  !! This initialization is done nqc_irr times from elphon_shuffle_wrap
  !! not all of the following code is necessary.  More adaptation from
  !! phq_init is needed   
  !!
  !! Roxana Margine - Dec 2018: Updated based on QE 6.3
  !!
  !
  USE kinds,            ONLY : DP
  USE ions_base,        ONLY : nat, ntyp => nsp, tau
  USE becmod,           ONLY : calbec, becp, allocate_bec_type
  USE lrus,             ONLY : becp1
  USE uspp,             ONLY : vkb, nlcc_any, okvan, nkb
  USE pwcom,            ONLY : npwx, nbnd, nks
  USE klist_epw,        ONLY : xk_loc, isk_loc
  USE constants,        ONLY : tpi
  USE constants_epw,    ONLY : zero, czero, cone
  USE cell_base,        ONLY : tpiba2, tpiba, omega
  USE klist,            ONLY : ngk, igk_k, nkstot
  USE gvect,            ONLY : g, ngm
  USE atom,             ONLY : msh, rgrid
  USE wavefunctions,    ONLY : evc
  USE noncollin_module, ONLY : noncolin, npol, nspin_mag
  USE uspp_param,       ONLY : upf, nhm
  USE m_gth,            ONLY : setlocq_gth
  USE units_lr,         ONLY : lrwfc, iuwfc
  USE phcom,            ONLY : vlocq
  USE qpoint,           ONLY : xq, eigqts
  USE nlcc_ph,          ONLY : drc                           
  USE elph2,            ONLY : igk_k_all, ngk_all
  USE mp,               ONLY : mp_barrier
  USE mp_global,        ONLY : inter_pool_comm, my_pool_id
  USE spin_orb,         ONLY : lspinorb
  USE lsda_mod,         ONLY : nspin, lsda, current_spin
  USE phus,             ONLY : int1, int1_nc, int2, int2_so,        &
                               int4, int4_nc, int5, int5_so, alphap
  USE poolgathering,    ONLY : poolgather_int, poolgather_int1
  USE io_epw,           ONLY : readwfc
  USE dvqpsi,           ONLY : dvanqq2
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(in) :: first_run
  !
  ! Local variables
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
  INTEGER :: ierr
  !! Error status
  REAL(KIND = DP) :: arg
  !! the argument of the phase
  COMPLEX(KIND = DP), ALLOCATABLE :: aux1(:, :)
  !! used to compute alphap
  !
  CALL start_clock('epw_init')
  ! 
  IF (first_run) THEN
    ALLOCATE(vlocq(ngm, ntyp), STAT = ierr)
    IF (ierr /= 0) CALL errore('epw_init', 'Error allocating vlocq', 1)
    ALLOCATE(eigqts(nat), STAT = ierr)
    IF (ierr /= 0) CALL errore('epw_init', 'Error allocating eigqts', 1)
    IF (okvan) THEN
      ALLOCATE(int1(nhm, nhm, 3, nat, nspin_mag), STAT = ierr)
      IF (ierr /= 0) CALL errore('epw_init', 'Error allocating int1', 1)
      ALLOCATE(int2(nhm, nhm, 3, nat, nat), STAT = ierr)
      IF (ierr /= 0) CALL errore('epw_init', 'Error allocating int2', 1)
      ALLOCATE(int4(nhm * (nhm + 1) / 2, 3, 3, nat, nspin_mag), STAT = ierr)
      IF (ierr /= 0) CALL errore('epw_init', 'Error allocating int4(nhm * ', 1)
      ALLOCATE(int5(nhm * (nhm + 1) / 2, 3, 3, nat , nat), STAT = ierr)
      IF (ierr /= 0) CALL errore('epw_init', 'Error allocating int5(nhm * ', 1)
      IF (noncolin) THEN
        ALLOCATE(int1_nc(nhm, nhm, 3, nat, nspin), STAT = ierr)
        IF (ierr /= 0) CALL errore('epw_init', 'Error allocating int1_nc', 1)
        ALLOCATE(int4_nc(nhm, nhm, 3, 3, nat, nspin), STAT = ierr)
        IF (ierr /= 0) CALL errore('epw_init', 'Error allocating int4_nc', 1)
        IF (lspinorb) THEN
          ALLOCATE(int2_so(nhm, nhm, 3, nat, nat, nspin), STAT = ierr)
          IF (ierr /= 0) CALL errore('epw_init', 'Error allocating int2_so', 1)
          ALLOCATE(int5_so(nhm, nhm, 3, 3, nat, nat, nspin), STAT = ierr)
          IF (ierr /= 0) CALL errore('epw_init', 'Error allocating int5_so', 1)
        ENDIF
      ENDIF ! noncolin
    ENDIF ! okvan
    !  
    ALLOCATE(becp1(nks), STAT = ierr)
    IF (ierr /= 0) CALL errore('epw_init', 'Error allocating becp1', 1)
    ALLOCATE(alphap(3, nks), STAT = ierr)
    IF (ierr /= 0) CALL errore('epw_init', 'Error allocating alphap', 1)
    ! 
    DO ik = 1, nks
      CALL allocate_bec_type(nkb, nbnd, becp1(ik))
      DO ipol = 1, 3
        CALL allocate_bec_type(nkb, nbnd, alphap(ipol,ik))
      ENDDO
    ENDDO
    CALL allocate_bec_type(nkb, nbnd, becp)
  ENDIF
  ! 
  DO na = 1, nat
    !
    ! xq here is the first q of the star
    arg = (xq(1) * tau(1, na) + &
           xq(2) * tau(2, na) + &
           xq(3) * tau(3, na)) * tpi
    !        
    eigqts(na) = CMPLX(COS(arg), - SIN(arg), KIND = DP)
    !
  END DO
  !
  ! ... a0) compute rhocore for each atomic-type if needed for nlcc
  !
  IF (nlcc_any) CALL set_drhoc(xq, drc)
  !
  ! ... b) the fourier components of the local potential at q+G
  !
  vlocq(:, :) = zero
  !
  DO nt = 1, ntyp
    !
    IF (upf(nt)%is_gth) THEN
      CALL setlocq_gth(nt, xq, upf(nt)%zp, tpiba2, ngm, g, omega, vlocq(1, nt))
    ELSE
      CALL setlocq(xq, rgrid(nt)%mesh, msh(nt), rgrid(nt)%rab, rgrid(nt)%r, &
                   upf(nt)%vloc(1), upf(nt)%zp, tpiba2, ngm, g, omega, vlocq(1, nt))
    ENDIF
    !
  END DO
  !
  ALLOCATE(aux1(npwx * npol, nbnd), STAT = ierr)
  IF (ierr /= 0) CALL errore('epw_init', 'Error allocating aux1', 1)
  ! 
  DO ik = 1, nks
    !
    !
    IF (lsda) current_spin = isk_loc(ik)
    !
    ! ... d) The functions vkb(k+G)
    !
    CALL init_us_2(ngk(ik), igk_k(1, ik), xk_loc(1, ik), vkb)
    !
    ! ... read the wavefunctions at k
    !
    CALL readwfc(my_pool_id + 1, ik, evc)
    !
    ! ... e) we compute the becp terms which are used in the rest of
    ! ...    the code
    !
    CALL calbec(ngk(ik), vkb, evc, becp1(ik))
    !
    ! we compute the derivative of the becp term with respect to an
    !   atomic displacement
    !
    DO ipol = 1, 3
      aux1 = czero
      DO ibnd = 1, nbnd
        DO ig = 1, ngk(ik)
          aux1(ig, ibnd) = evc(ig, ibnd) * tpiba * cone * & 
                          (xk_loc(ipol, ik) + g(ipol, igk_k(ig, ik)))
        ENDDO
        IF (noncolin) THEN
          DO ig = 1, ngk(ik)
            aux1(ig + npwx, ibnd) = evc(ig + npwx, ibnd) * tpiba *cone *& 
                      (xk_loc(ipol, ik) + g(ipol, igk_k(ig, ik)) )
          ENDDO
        ENDIF
      ENDDO
      CALL calbec(ngk(ik), vkb, aux1, alphap(ipol, ik)) 
    ENDDO
    !
  ENDDO
  !
  DEALLOCATE(aux1, STAT = ierr)
  IF (ierr /= 0) CALL errore('epw_init', 'Error deallocating aux1', 1)
  !
  IF (first_run) THEN
    ALLOCATE(igk_k_all(npwx, nkstot), STAT = ierr)
    IF (ierr /= 0) CALL errore('epw_init', 'Error allocating igk_k_all', 1)
    ALLOCATE(ngk_all(nkstot), STAT = ierr)
    IF (ierr /= 0) CALL errore('epw_init', 'Error allocating ngk_all', 1)
  ENDIF
  !
#if defined(__MPI)
  !
  CALL poolgather_int(npwx, nkstot, nks, igk_k(:, 1:nks), igk_k_all) 
  CALL poolgather_int1(nkstot, nks, ngk(1:nks), ngk_all) 
  !CALL mp_barrier(inter_pool_comm)
  !
#else
  !
  igk_k_all = igk_k
  ngk_all = ngk
  !
#endif
  !
  IF (.NOT. first_run) THEN
    CALL dvanqq2()
  ENDIF
  !
  CALL stop_clock('epw_init')
  !
  !----------------------------------------------------------------------------
  END SUBROUTINE epw_init
  !----------------------------------------------------------------------------
