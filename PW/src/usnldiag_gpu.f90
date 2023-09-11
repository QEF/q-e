!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE usnldiag_gpu( npw, h_diag_d, s_diag_d )
  !-----------------------------------------------------------------------
  !! Add nonlocal pseudopotential term to diagonal part of Hamiltonian.  
  !! Compute the diagonal part of the S matrix.
  !
  !    Routine splitted to improve performance
  !
  USE kinds,            ONLY: DP
  USE ions_base,        ONLY: nat, ityp, ntyp => nsp
  USE wvfct,            ONLY: npwx
  USE uspp,             ONLY: ofsbeta, deeq, qq_at, qq_so, &
                              deeq_nc
  USE uspp_param,       ONLY: upf, nh
  USE noncollin_module, ONLY: noncolin, npol, lspinorb
  USE device_memcpy_m,  ONLY: dev_memset
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: npw
  ! number of plane waves
  REAL(dp), INTENT(inout) :: h_diag_d (npwx,npol)
  ! the diagonal part of the hamiltonian
  REAL(dp), INTENT(out)   :: s_diag_d (npwx,npol)
  ! the diagonal part of the S matrix
#if defined(__CUDA)
  attributes(DEVICE) :: h_diag_d, s_diag_d
#endif
  !
  INTEGER :: ig, ipol
  !
  ! initialise s_diag
  !
  CALL dev_memset( s_diag_d, 1.d0 )
  !
  IF (lspinorb) THEN
     CALL usnldiag_spinorb()
  ELSEIF (noncolin) THEN
     CALL usnldiag_noncollinear()
  ELSE
     CALL usnldiag_collinear()
  END IF
  
CONTAINS
  
  SUBROUTINE usnldiag_collinear()
     USE lsda_mod, ONLY: current_spin
     USE uspp,     ONLY: deeq, qq_at, vkb
     
     IMPLICIT NONE
     !
     INTEGER :: ijkb_start, ikb, jkb, ih, jh, na, nt, ig, ipol
     COMPLEX(DP) :: ps1_1, ps1_2, ps2_1, ps2_2, ar, cv
     !
     INTEGER :: nh_
     REAL(DP) :: sum_h, sum_s
     !
     !    multiply on projectors
     !
     DO nt = 1, ntyp
        IF ( upf(nt)%tvanp .or. upf(nt)%is_multiproj ) THEN
           DO na = 1, nat
              IF (ityp (na) == nt) THEN
                   ijkb_start = ofsbeta(na)
                   nh_ = nh(nt)
                   !$acc data present(vkb(:,:),deeq) deviceptr(h_diag_d(:,:), s_diag_d(:,:))
                   !$acc parallel vector_length(32)
                   !$acc loop gang reduction(+:sum_h,sum_s)
                   DO ig = 1, npw 
                      sum_h = 0.d0
                      sum_s = 0.d0
                      !$acc loop vector collapse(2) private(ikb,cv,jkb,ar) reduction(+:sum_h,sum_s)
                      DO ih = 1, nh_
                         DO jh = 1, nh_
                            ikb = ijkb_start + ih
                            cv = vkb(ig,ikb)
                            jkb = ijkb_start + jh
                            ar = cv*conjg(vkb(ig,jkb))
                            sum_h = sum_h + dble(deeq(ih,jh,na,current_spin) * ar)
                            sum_s = sum_s + dble(qq_at(ih,jh,na) * ar)
                         END DO
                      END DO
                      !$acc atomic update 
                      h_diag_d (ig,1) = h_diag_d (ig,1) + sum_h
                      !$acc end atomic 
                      !$acc atomic update 
                      s_diag_d (ig,1) = s_diag_d (ig,1) + sum_s
                      !$acc end atomic 
                   ENDDO
                   !$acc end parallel
                   !$acc end data
              END IF
           END DO
        ELSE
           DO na = 1, nat
              IF (ityp (na) == nt) THEN
                   ijkb_start = ofsbeta(na)
                   nh_ = nh(nt)
                   !$acc data present(vkb(:,:),deeq) deviceptr(h_diag_d(:,:), s_diag_d(:,:))
                   !$acc parallel vector_length(32)
                   !$acc loop gang reduction(+:sum_h,sum_s)
                   DO ig = 1, npw 
                      sum_h = 0.d0
                      sum_s = 0.d0
                      !$acc loop vector private(ikb,ar) reduction(+:sum_h,sum_s)
                      DO ih = 1, nh_
                         ikb = ijkb_start + ih
                         ar = vkb(ig,ikb)*conjg(vkb(ig,ikb))
                         sum_h = sum_h + dble(deeq(ih,ih,na,current_spin) * ar)
                         sum_s = sum_s + dble(qq_at(ih,ih,na) * ar)
                      END DO
                      !$acc atomic update
                      h_diag_d (ig,1) = h_diag_d (ig,1) + sum_h
                      !$acc end atomic 
                      !$acc atomic update
                      s_diag_d (ig,1) = s_diag_d (ig,1) + sum_s
                      !$acc end atomic 
                   ENDDO
                   !$acc end parallel
                   !$acc end data
              END IF
           END DO
        END IF
     END DO
  END SUBROUTINE usnldiag_collinear
  !
  SUBROUTINE usnldiag_noncollinear()
     USE lsda_mod,  ONLY: current_spin
     USE uspp,      ONLY: vkb, qq_at, qq_so, deeq_nc
     
     IMPLICIT NONE
     !
     INTEGER :: ijkb_start, ikb, jkb, ih, jh, na, nt, ig, ipol
     COMPLEX(DP) :: ps1_1, ps1_2, ps2_1, ps2_2, ar, cv
     !
     INTEGER :: nh_
     REAL(DP) :: sum_h1, sum_h4, sum_s
     !
     !    multiply on projectors
     !
     DO nt = 1, ntyp
        IF ( upf(nt)%tvanp .or. upf(nt)%is_multiproj ) THEN
           DO na = 1, nat
              IF (ityp (na) == nt) THEN
                   ijkb_start = ofsbeta(na)
                   nh_ = nh(nt)
                   !$acc data present(vkb(:,:),deeq_nc) deviceptr(h_diag_d(:,:), s_diag_d(:,:))
                   !$acc parallel vector_length(32) 
                   !$acc loop gang reduction(+:sum_h1,sum_h4,sum_s)
                   DO ig = 1, npw   ! change this to 2*npw ?
                      sum_h1 = 0.d0
                      sum_h4 = 0.d0
                      sum_s = 0.d0
                      !$acc loop vector collapse(2) private(ikb,cv,jkb,ar) reduction(+:sum_h1,sum_h4,sum_s)
                      DO ih = 1, nh_
                         DO jh = 1, nh_
                            ikb = ijkb_start + ih
                            cv = vkb (ig, ikb)
                            jkb = ijkb_start + jh
                            ar = cv*conjg(vkb(ig,jkb))
                            sum_h1 = sum_h1 + dble(deeq_nc(ih,jh,na,1) * ar)
                            sum_h4 = sum_h4 + dble(deeq_nc(ih,jh,na,4) * ar)
                            sum_s  = sum_s  + dble(qq_at(ih,jh,na) * ar)
                         END DO
                      END DO
                      !
                      ! OPTIMIZE HERE : this scattered assign is bad!
                      !
                      !$acc atomic update
                      h_diag_d (ig,1) = h_diag_d (ig,1) + sum_h1
                      !$acc end atomic 
                      !$acc atomic update
                      h_diag_d (ig,2) = h_diag_d (ig,2) + sum_h4
                      !$acc end atomic 
                      !$acc atomic update
                      s_diag_d (ig,1) = s_diag_d (ig,1) + sum_s
                      !$acc end atomic 
                      !$acc atomic update
                      s_diag_d (ig,2) = s_diag_d (ig,2) + sum_s
                      !$acc end atomic 
                   ENDDO
                   !$acc end parallel
                   !$acc end data 
              END IF
           END DO
        ELSE
           DO na = 1, nat
              IF (ityp (na) == nt) THEN
                   ijkb_start = ofsbeta(na)
                   nh_ = nh(nt)
                   !$acc data present(vkb(:,:),deeq_nc) deviceptr(h_diag_d(:,:), s_diag_d(:,:))
                   !$acc parallel vector_length(32) 
                   !$acc loop gang reduction(+:sum_h1,sum_h4,sum_s)
                   DO ig = 1, npw 
                      sum_h1 = 0.d0
                      sum_h4 = 0.d0
                      sum_s = 0.d0
                      !$acc loop vector private(ikb,ar) reduction(+:sum_h1,sum_h4,sum_s)
                      DO ih = 1, nh_
                         ikb = ijkb_start + ih
                         ar = vkb (ig, ikb)*conjg(vkb(ig,ikb))
                         sum_h1 = sum_h1 + dble(deeq_nc(ih,ih,na,1) * ar)
                         sum_h4 = sum_h4 + dble(deeq_nc(ih,ih,na,4) * ar)
                         sum_s = sum_s + dble(qq_at(ih,ih,na) * ar)
                      END DO
                      !
                      ! OPTIMIZE HERE : this scattered assign is bad!
                      !
                      !$acc atomic update
                      h_diag_d (ig,1) = h_diag_d (ig,1) + sum_h1
                      !$acc end atomic 
                      !$acc atomic update
                      h_diag_d (ig,2) = h_diag_d (ig,2) + sum_h4
                      !$acc end atomic 
                      !$acc atomic update
                      s_diag_d (ig,1) = s_diag_d (ig,1) + sum_s
                      !$acc end atomic 
                      !$acc atomic update
                      s_diag_d (ig,2) = s_diag_d (ig,2) + sum_s
                      !$acc end atomic 
                   ENDDO
                   !$acc end parallel
                   !$acc end data 
              END IF
           END DO
        END IF
     END DO
  END SUBROUTINE usnldiag_noncollinear
  !
  SUBROUTINE usnldiag_spinorb()
     USE lsda_mod, ONLY: current_spin
     USE uspp,     ONLY: vkb, qq_at, qq_so, deeq_nc

     IMPLICIT NONE
     !
     INTEGER :: ijkb_start, ikb, jkb, ih, jh, na, nt, ig, ipol
     COMPLEX(DP) :: ps1_1, ps1_2, ps2_1, ps2_2, ar, cv
     !
     INTEGER :: nh_
     REAL(DP) :: sum_h1, sum_h4, sum_s1, sum_s4
     !
     !    multiply on projectors
     !
     DO nt = 1, ntyp
        IF ( upf(nt)%tvanp .or. upf(nt)%is_multiproj ) THEN
           DO na = 1, nat
              IF (ityp (na) == nt) THEN
                   ijkb_start = ofsbeta(na)
                   nh_ = nh(nt)
                   !$acc data present(vkb(:,:),deeq_nc) deviceptr(h_diag_d(:,:), s_diag_d(:,:))
                   !$acc parallel vector_length(32)
                   !$acc loop gang reduction(+:sum_h1,sum_h4,sum_s1,sum_s4)
                   DO ig = 1, npw   ! change this to 2*npw ?
                      sum_h1 = 0.d0
                      sum_h4 = 0.d0
                      sum_s1 = 0.d0
                      sum_s4 = 0.d0
                      !$acc loop vector collapse(2) private(ikb,cv,jkb,ar) reduction(+:sum_h1,sum_h4,sum_s1,sum_s4) 
                      DO ih = 1, nh_
                         DO jh = 1, nh_
                            ikb = ijkb_start + ih
                            cv = vkb(ig,ikb)
                            jkb = ijkb_start + jh
                            ar = cv*conjg(vkb(ig,jkb))
                            sum_h1 = sum_h1 + dble(deeq_nc(ih,jh,na,1) * ar)
                            sum_h4 = sum_h4 + dble(deeq_nc(ih,jh,na,4) * ar)
                            sum_s1 = sum_s1 + dble(qq_so(ih,jh,1,nt) * ar)
                            sum_s4 = sum_s4 + dble(qq_so(ih,jh,4,nt) * ar)
                         END DO
                      END DO
                      !
                      ! OPTIMIZE HERE : this scattered assign is bad!
                      !
                      !$acc atomic update
                      h_diag_d (ig,1) = h_diag_d (ig,1) + sum_h1
                      !$acc end atomic 
                      !$acc atomic update
                      h_diag_d (ig,2) = h_diag_d (ig,2) + sum_h4
                      !$acc end atomic 
                      !$acc atomic update
                      s_diag_d (ig,1) = s_diag_d (ig,1) + sum_s1
                      !$acc end atomic 
                      !$acc atomic update
                      s_diag_d (ig,2) = s_diag_d (ig,2) + sum_s4
                      !$acc end atomic 
                   ENDDO
                   !$acc end parallel
                   !$acc end data 
              END IF
           END DO
        ELSE
           DO na = 1, nat
              IF (ityp (na) == nt) THEN
                   ijkb_start = ofsbeta(na)
                   nh_ = nh(nt)
                   !$acc data present(vkb(:,:),deeq_nc) deviceptr(h_diag_d(:,:), s_diag_d(:,:))
                   !$acc parallel vector_length(32)
                   !$acc loop gang reduction(+:sum_h1,sum_h4,sum_s1,sum_s4)
                   DO ig = 1, npw 
                      sum_h1 = 0.d0
                      sum_h4 = 0.d0
                      sum_s1 = 0.d0
                      sum_s4 = 0.d0
                      !$acc loop vector private(ikb,ar) reduction(+:sum_h1,sum_h4,sum_s1,sum_s4) 
                      DO ih = 1, nh_
                         ikb = ijkb_start + ih
                         ar = vkb (ig, ikb)*conjg(vkb (ig, ikb))
                         sum_h1 = sum_h1 + dble(deeq_nc(ih,ih,na,1) * ar)
                         sum_h4 = sum_h4 + dble(deeq_nc(ih,ih,na,4) * ar)
                         sum_s1 = sum_s1 + dble(qq_so(ih,ih,1,nt) * ar)
                         sum_s4 = sum_s4 + dble(qq_so(ih,ih,4,nt) * ar)
                      END DO
                      !
                      ! OPTIMIZE HERE : this scattered assign is bad!
                      !
                      !$acc atomic update
                      h_diag_d (ig,1) = h_diag_d (ig,1) + sum_h1
                      !$acc end atomic 
                      !$acc atomic update
                      h_diag_d (ig,2) = h_diag_d (ig,2) + sum_h4
                      !$acc end atomic 
                      !$acc atomic update
                      s_diag_d (ig,1) = s_diag_d (ig,1) + sum_s1
                      !$acc end atomic 
                      !$acc atomic update
                      s_diag_d (ig,2) = s_diag_d (ig,2) + sum_s4
                      !$acc end atomic 
                   ENDDO
                   !$acc end parallel
                   !$acc end data 
              END IF
           END DO
        END IF
     END DO
  END SUBROUTINE usnldiag_spinorb
END SUBROUTINE usnldiag_gpu
