!
! Copyright (C) 2001-2023 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE hp_dnsq (lmetq0, iter, conv_root, dnsq) 
  !-----------------------------------------------------------------------
  ! 
  ! This routine calculates, at each SCF iteration, 
  ! the variation of the occupation matrix dnsq for a fixed q.
  ! See Eq. (43) in Ref. [1].
  ! [1] Phys. Rev. B 98, 085127 (2018)
  ! 
  ! dnsq(m1,m2,ispin,I) = 
  !  = \sum_{n,k} [ <psi(n,k,ispin)| S(k)*phi_(k,I,m1)> * 
  !                   * < S(k+q)*phi(k+q,I,m2)| dpsi(n,k+q,ispin)> + 
  !               + <psi(n,k,ispin)| S(k)*phi_(k,I,m2)> * 
  !                   * < S(k+q)*phi(k+q,I,m1)| dpsi(n,k+q,ispin)> ]
  !      + the term due to shift of the Fermi energy (only if lmetq0=.true.)
  !
  USE kinds,                ONLY : DP
  USE io_files,             ONLY : nwordwfcU
  USE ions_base,            ONLY : nat, ityp
  USE klist,                ONLY : xk, wk, degauss, ngauss, ngk
  USE wvfct,                ONLY : npwx, wg, nbnd, et 
  USE uspp_param,           ONLY : nh 
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE wavefunctions,        ONLY : evc
  USE ener,                 ONLY : ef
  USE uspp,                 ONLY : okvan 
  USE buffers,              ONLY : get_buffer
  USE mp_pools,             ONLY : intra_pool_comm, inter_pool_comm       
  USE mp,                   ONLY : mp_sum
  USE io_global,            ONLY : stdout
  USE constants,            ONLY : rytoev
  USE control_flags,        ONLY : iverbosity
  USE qpoint,               ONLY : nksq, ikks, ikqs
  USE control_lr,           ONLY : lgamma, nbnd_occ
  USE units_lr,             ONLY : iuwfc, lrwfc, iuatswfc, iudwf, lrdwf
  USE lr_symm_base,         ONLY : nsymq
  USE ldaU,                 ONLY : Hubbard_lmax, Hubbard_l, is_hubbard, offsetU, nwfcU
  USE ldaU_hp,              ONLY : conv_thr_chi, trace_dns_tot_old, &
                                   conv_thr_chi_best, iter_best
  USE ldaU_lr,              ONLY : swfcatomk, swfcatomkpq
  USE efermi_shift,         ONLY : def
  USE noncollin_module,     ONLY : noncolin, npol, nspin_mag, domag
  USE qpoint_aux,           ONLY : ikmks
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN) :: lmetq0 
  ! Change of the Fermi energy (only for q=0 and metal) 
  INTEGER, INTENT(IN) :: iter
  ! Number of the current iteration
  LOGICAL, INTENT(OUT) :: conv_root
  ! .True. if convergence has been reached
  COMPLEX(DP), INTENT(OUT) :: dnsq(2*Hubbard_lmax+1, 2*Hubbard_lmax+1, nspin, nat)
  !
  !  the local variables
  !
  INTEGER ::  ik, ikk, ikq, i, j, k, ios, icar, isp1, isp2, &
              counter, nt, na, l, ih, jkb2, n, &
              ihubst, ihubst1, ihubst2, m, m1, m2, ibnd, is, ldim, &
              ikmk, isolv, nsolv, nrec, is1, is2
  INTEGER :: npw, npwq
  ! number of plane waves at k and k+q
  COMPLEX(DP), ALLOCATABLE :: dpsi(:, :), proj1(:,:), proj2(:,:)
  REAL(DP) :: weight, wdelta, w1
  REAL(DP),    EXTERNAL :: w0gauss ! an approximation to the delta function
  COMPLEX(DP) :: trace_dns(2)
  COMPLEX(DP), ALLOCATABLE :: trace_dns_tot(:)
  LOGICAL, ALLOCATABLE :: conv(:) ! If true the root is converged
  REAL(DP) :: diff_max
  REAL(DP), ALLOCATABLE :: diff(:)
  ! 
  CALL start_clock( 'hp_dnsq' )
  !
  ios = 0 
  ! 
  ALLOCATE (dpsi(npwx*npol,nbnd))
  ALLOCATE (proj1(nbnd,nwfcU))
  ALLOCATE (proj2(nbnd,nwfcU))
  ALLOCATE (trace_dns_tot(nat))
  ALLOCATE (conv(nat))
  ALLOCATE (diff(nat))
  !
  dpsi = (0.d0, 0.d0) 
  !
  ldim = 2 * Hubbard_lmax + 1
  ! 
  ! At each iteration dnsq is set to zero, because it is recomputed
  !
  dnsq(:,:,:,:) = (0.d0, 0.d0) 
  !
  DO na = 1, nat
     nt = ityp(na)
     IF (is_hubbard(nt)) THEN
        conv(na) = .false.
     ELSE
        conv(na) = .true.
     ENDIF
  ENDDO
  ! 
  nsolv = 1 
  IF (noncolin.AND.domag) nsolv = 2
  !
  DO ik = 1, nksq
     !
     ikk  = ikks(ik)
     ikq  = ikqs(ik)
     npw  = ngk(ikk)
     npwq = ngk(ikq)
     !
     DO isolv = 1, nsolv
        !
        IF (isolv == 2) THEN
           ikmk = ikmks(ik)
        ELSE
           ikmk = ikk
        ENDIF
        !
        IF (lsda) current_spin = isk (ikk)
        !
        ! Read unperturbed KS wfc's psi(k)
        IF (nksq.GT.1) CALL get_buffer (evc, lrwfc, iuwfc, ikmk)
        !
        ! Read atomic wfc's : S(k)*phi(k) and S(k+q)*phi(k+q) 
        CALL get_buffer (swfcatomk, nwordwfcU, iuatswfc, ikk)
        IF (.NOT.lgamma) CALL get_buffer (swfcatomkpq, nwordwfcU, iuatswfc, ikq)
        !
        ! At each SCF iteration for each ik read dpsi from file
        nrec = ik + (isolv - 1) * nksq 
        CALL get_buffer (dpsi, lrdwf, iudwf, nrec)
        ! 
        ! Loop on Hubbard atoms
        proj1 = (0.d0, 0.d0) 
        proj2 = (0.d0, 0.d0) 
        DO na = 1, nat
           nt = ityp(na)
           IF (is_hubbard(nt)) THEN
              ldim = (2*Hubbard_l(nt) + 1) * npol
              DO m = 1, ldim
                 ihubst = offsetU(na) + m   ! I m index
                 !
                 !  proj1(ibnd, ihubst) = < S(k)*phi(k) | psi(k) >
                 !  proj2(ibnd, ihubst) = < S(k+q)*phi(k+q) | dpsi(k+q) > 
                 ! FIXME: use ZGEMM instead of dot_product
                 DO ibnd = 1, nbnd_occ(ikk)
                    proj1(ibnd, ihubst) = DOT_PRODUCT( swfcatomk(1:npwx*npol,ihubst), &
                                                       evc(1:npwx*npol,ibnd) )
                    proj2(ibnd, ihubst) = DOT_PRODUCT( swfcatomkpq(1:npwx*npol,ihubst), &
                                                       dpsi(1:npwx*npol,ibnd) )
                 ENDDO
              ENDDO 
           ENDIF
           !
        ENDDO
        !
        CALL mp_sum(proj1, intra_pool_comm)  
        CALL mp_sum(proj2, intra_pool_comm)
        !
        IF (noncolin) THEN
            DO na = 1, nat
               nt = ityp(na)
               IF ( is_hubbard(nt) ) THEN
                  ldim = 2*Hubbard_l(nt)+1
                  DO is1 = 1, npol
                     DO is2 = 1, npol
                        is = npol*(is1-1) + is2 
                        isp1 = is1 + (is2 - is1)*(isolv - 1)
                        isp2 = is2 + (is1 - is2)*(isolv - 1)
                        DO m1 = 1, ldim
                           DO m2 = 1, ldim
                              ihubst2 = offsetU(na) + m2 + (m1 - m2)*(isolv - 1)
                              ihubst1 = offsetU(na) + m1 + (m2 - m1)*(isolv - 1)
                              DO ibnd = 1, nbnd_occ(ikk) 
                                 dnsq(m1,m2,is,na) = dnsq(m1,m2,is,na) + &
                                 wk(ikk) * CONJG(proj1(ibnd,ihubst1+ldim*(isp1-1)) ) * &
                                                 proj2(ibnd,ihubst2+ldim*(isp2-1))
                                 IF (.not.domag) THEN 
                                      dnsq(m1,m2,is,na) = dnsq(m1,m2,is,na) + &
                                      wk(ikk) * CONJG(proj1(ibnd,ihubst2+ldim*(isp2-1)) ) * &
                                                      proj2(ibnd,ihubst1+ldim*(isp1-1))
                                 ENDIF
               
                                 ! 
                                 ! correction for metals at q=0
                                 IF (lmetq0.and.isolv==1) then
                                    !
                                    ! wdelta: smeared delta function at the Fermi level
                                    wdelta = w0gauss ( (ef-et(ibnd,ikk)) / degauss, ngauss) / degauss
                                    weight = wk(ikk)  
                                    w1 = weight * wdelta
                                    !
                                    dnsq(m1, m2, is, na) = dnsq(m1, m2, is, na) +  &
                                                w1 * def(1) * CONJG(proj1(ibnd,ihubst1+ldim*(is1-1))) &
                                                                  * proj1(ibnd,ihubst2+ldim*(is2-1))
                                 ENDIF
                              ENDDO
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDIF
            ENDDO
        ELSE
            DO na = 1, nat
               nt = ityp(na)
               IF (is_hubbard(nt)) THEN
                  DO m1 = 1, 2 * Hubbard_l(nt) + 1
                     ihubst1 = offsetU(na) + m1
                     DO m2 = m1, 2 * Hubbard_l(nt) + 1
                        ihubst2 = offsetU(na) + m2
                        DO ibnd = 1, nbnd_occ(ikk)
                           dnsq(m1, m2, current_spin, na) = dnsq(m1, m2, current_spin, na) +    &
                                    wk(ikk) * ( CONJG(proj1(ibnd,ihubst1)) * proj2(ibnd,ihubst2) + &
                                                CONJG(proj1(ibnd,ihubst2)) * proj2(ibnd,ihubst1) )
                           !
                           ! Correction for metals at q=0                                                
                           IF (lmetq0) THEN
                              !
                              ! wdelta: smeared delta function at the Fermi level
                              wdelta = w0gauss ( (ef-et(ibnd,ikk)) / degauss, ngauss) / degauss
                              weight = wk(ikk)  
                              w1 = weight * wdelta
                              !
                              dnsq(m1, m2, current_spin, na) = dnsq(m1, m2, current_spin, na) +  &
                              w1 * def(1) * CONJG(proj1(ibnd,ihubst1)) * proj1(ibnd,ihubst2)
                              !
                           ENDIF
                        ENDDO 
                     ENDDO 
                  ENDDO
               ENDIF 
            ENDDO
         ENDIF
      ENDDO ! isolv 
     ! 
  ENDDO ! ik 
  !
#if defined (__MPI)
  CALL mp_sum(dnsq, inter_pool_comm) 
#endif
  !
  ! dnsq is symmetric: filling the m1 m2 lower triangular part 
  !
  IF (.NOT.noncolin) THEN
     ldim = 2 * Hubbard_lmax + 1
     DO m1 = 2, ldim
        DO m2 = 1, m1-1
           dnsq(m1,m2,:,:) = dnsq(m2,m1,:,:)
        ENDDO
     ENDDO
  ENDIF
  !
  ! If nspin=1, k-point's weight wk is normalized to 2 el/band in the whole BZ,
  ! but we are interested in dnsq of one spin component
  !
  IF (nspin.EQ.1) dnsq = 0.5d0 * dnsq
  ! 
  IF ( iverbosity > 3 ) THEN
     WRITE(stdout,'(6x,"RESPONSE OCCUPATION MATRICES:")')
     DO na = 1, nat
        nt = ityp(na)
        IF (is_hubbard(nt)) THEN
         ldim = 2*Hubbard_l(nt)+1
           DO is = 1, nspin
              WRITE(stdout,'(5x,a,1x,i2,2x,a,1x,i2)') ' Hubbard atom', na, 'spin', is
              DO m1 = 1, ldim 
                 WRITE(stdout,'(3x,10(f15.10,1x))') ((dnsq(m1,m2,is,na)), &
                                                         m2=1,ldim)
              ENDDO
           ENDDO
        ENDIF
     ENDDO
  ENDIF
  !
  ! Symmetrization of dnsq
  ! 
  CALL hp_symdnsq (dnsq)
  !
  IF ( nsymq > 1 .AND. iverbosity > 3 ) THEN
     WRITE(stdout,'(6x,"RESPONSE OCCUPATION MATRICES (SYMMETRIZED):")')
     DO na = 1, nat
        nt = ityp(na)
        IF (is_hubbard(nt)) THEN
         ldim = 2*Hubbard_l(nt)+1
           DO is = 1, nspin
              WRITE(stdout,'(5x,a,1x,i2,2x,a,1x,i2)') ' Hubbard atom', na, 'spin', is
              DO m1 = 1, ldim
                 WRITE(stdout,'(3x,10(f15.10,1x))') (DBLE(dnsq(m1,m2,is,na)), &
                                                         m2=1,ldim)
              ENDDO
           ENDDO
        ENDIF
     ENDDO
  ENDIF
  !
  ! Write the trace of dnsq for a given iteration
  !
  trace_dns_tot(:) = 0.d0
  !
  conv_root = .true.
  !
  diff_max = 0.0d0
  diff(:)  = 0.0d0
  !
  DO na = 1, nat
     !
     nt = ityp(na)
     ! 
     IF (is_hubbard(nt)) THEN 
        ldim = 2*Hubbard_l(nt)+1
        !
        ! Divide by rytoev -> conversion of units from Ry to eV
        !
        trace_dns(:) = 0.d0
        DO is = 1, nspin/npol
           DO m = 1, ldim
              trace_dns(is) = trace_dns(is) + dnsq(m,m,is**npol,na)/rytoev
           ENDDO
           trace_dns_tot(na) = trace_dns_tot(na) + trace_dns(is)
        ENDDO
        !
        ! If nspin=1, multiply back by a factor of 2 due to spin degeneracy
        !
        IF (nspin.EQ.1) trace_dns_tot(na) = 2.0d0 * trace_dns_tot(na) 
        !
        ! Check for convergence
        !
        diff(na) = DABS( DBLE(trace_dns_tot(na)) - DBLE(trace_dns_tot_old(na)) )
        !
        IF ( diff(na) < conv_thr_chi ) conv(na) = .true.
        !
        conv_root = conv_root .AND. conv(na) 
        !
        ! Write the trace of dnsq to file
        !
        IF ( iter == 1 ) THEN
           WRITE(stdout,'(6x,"chi:",1x,i3,1x,f14.10,3x)') &
                                & na, DBLE(trace_dns_tot(na))
        ELSE
           WRITE(stdout,'(6x,"chi:",1x,i3,1x,f14.10,3x,"residue:"2x,f14.10)') &
                                & na, DBLE(trace_dns_tot(na)), diff(na)
        ENDIF
        !
        ! Save for comparison in the next iteration
        !
        trace_dns_tot_old(na) = trace_dns_tot(na)
        !
     ENDIF
     !
  ENDDO
  !
  ! Find the largest error among diff(:) for the current iteration
  ! 
  diff_max = MAXVAL(diff(:))
  !
  IF ( iter == 1 ) THEN
     !
     conv_root = .false.
     conv_thr_chi_best = 1000.0d0
     iter_best = 1
     GO TO 101
     !
  ENDIF
  !
  IF ( diff_max < conv_thr_chi_best ) THEN
     !
     conv_thr_chi_best = diff_max
     iter_best = iter
     !
  ENDIF
  !
101 CONTINUE
  !
  DEALLOCATE (proj1)
  DEALLOCATE (proj2)
  DEALLOCATE (dpsi)
  DEALLOCATE (trace_dns_tot)
  DEALLOCATE (conv)
  DEALLOCATE (diff)
  !
  CALL stop_clock( 'hp_dnsq' )
  ! 
  RETURN
  !
END SUBROUTINE hp_dnsq
