!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE usnldiag_gpu (npw, h_diag_d, s_diag_d)
  !-----------------------------------------------------------------------
  !
  !    add nonlocal pseudopotential term to diagonal part of Hamiltonian
  !    compute the diagonal part of the S matrix.
  !
  !    Routine splitted for improving performance
  !
  USE kinds, ONLY: DP
  USE ions_base,  ONLY : nat, ityp, ntyp => nsp
  USE wvfct, ONLY: npwx
  USE uspp,  ONLY: indv_ijkb0
  USE uspp_gpum,  ONLY: deeq_d, vkb_d, qq_at_d, qq_so_d, deeq_nc_d
  USE uspp_param, ONLY: upf, nh, newpseudo
  USE spin_orb, ONLY: lspinorb
  USE noncollin_module, ONLY: noncolin, npol
  !
  USE uspp_gpum, ONLY : using_vkb_d, using_indv_ijkb0, using_deeq_d, using_deeq_nc_d, &
                        using_qq_at_d, using_qq_so_d
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
  CALL using_vkb_d(0)
  CALL using_indv_ijkb0(0)
  CALL using_deeq_d(0)
  IF (lspinorb .or. noncolin) CALL using_deeq_nc_d(0)
  IF (.not. lspinorb)         CALL using_qq_at_d(0)
  IF (lspinorb)               CALL using_qq_so_d(0)
  !
  ! initialise s_diag
  !
  s_diag_d = 1.d0
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
     USE uspp_gpum,  ONLY: deeq_d, vkb_d, qq_at_d
     
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
        IF ( upf(nt)%tvanp .or.newpseudo (nt) ) THEN
           DO na = 1, nat
              IF (ityp (na) == nt) THEN
                   ijkb_start = indv_ijkb0(na)
                   nh_ = nh(nt)
!$cuf kernel do(1) <<<*,*>>>
                   DO ig = 1, npw 
                      sum_h = 0.d0
                      sum_s = 0.d0
                      
                      DO ih = 1, nh_
                         ikb = ijkb_start + ih
                         cv = vkb_d (ig, ikb)
                         DO jh = 1, nh_
                            jkb = ijkb_start + jh
                            
                            ar = cv*conjg(vkb_d (ig, jkb))
                            
                            sum_h = sum_h + dble(deeq_d (ih, jh, na, current_spin) * ar)
                            sum_s = sum_s + dble(qq_at_d (ih, jh, na) * ar)
                         END DO
                      END DO
                       
                      h_diag_d (ig,1) = h_diag_d (ig,1) + sum_h
                      s_diag_d (ig,1) = s_diag_d (ig,1) + sum_s
                   ENDDO
              END IF
           END DO
        ELSE
           DO na = 1, nat
              IF (ityp (na) == nt) THEN
                   ijkb_start = indv_ijkb0(na)
                   nh_ = nh(nt)
!$cuf kernel do(1) <<<*,*>>>
                   DO ig = 1, npw 
                      sum_h = 0.d0
                      sum_s = 0.d0
                      
                      DO ih = 1, nh_
                         ikb = ijkb_start + ih
                         
                         ar = vkb_d (ig, ikb)*conjg(vkb_d (ig, ikb))
                         
                         sum_h = sum_h + dble(deeq_d (ih, ih, na, current_spin) * ar)
                         sum_s = sum_s + dble(qq_at_d (ih, ih, na) * ar)
                      END DO
                       
                      h_diag_d (ig,1) = h_diag_d (ig,1) + sum_h
                      s_diag_d (ig,1) = s_diag_d (ig,1) + sum_s
                   ENDDO
              END IF
           END DO
        END IF
     END DO
  END SUBROUTINE usnldiag_collinear
  !
  SUBROUTINE usnldiag_noncollinear()
     USE lsda_mod, ONLY: current_spin
     USE uspp_gpum,  ONLY: vkb_d, qq_at_d, qq_so_d, deeq_nc_d
     
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
        IF ( upf(nt)%tvanp .or.newpseudo (nt) ) THEN
           DO na = 1, nat
              IF (ityp (na) == nt) THEN
                   ijkb_start = indv_ijkb0(na)
                   nh_ = nh(nt)
!$cuf kernel do(1) <<<*,*>>>
                   DO ig = 1, npw   ! change this to 2*npw ?
                      sum_h1 = 0.d0
                      sum_h4 = 0.d0
                      sum_s = 0.d0
                      
                      DO ih = 1, nh_
                         ikb = ijkb_start + ih
                         cv = vkb_d (ig, ikb)
                         DO jh = 1, nh_
                            jkb = ijkb_start + jh
                            
                            ar = cv*conjg(vkb_d (ig, jkb))
                            
                            sum_h1 = sum_h1 + dble(deeq_nc_d (ih, jh, na, 1) * ar)
                            sum_h4 = sum_h4 + dble(deeq_nc_d (ih, jh, na, 4) * ar)
                            sum_s  = sum_s  + dble(qq_at_d (ih, jh, na) * ar)
                         END DO
                      END DO
                      !
                      ! OPTIMIZE HERE : this scattered assign is bad!
                      !
                      h_diag_d (ig,1) = h_diag_d (ig,1) + sum_h1
                      h_diag_d (ig,2) = h_diag_d (ig,2) + sum_h4
                      s_diag_d (ig,1) = s_diag_d (ig,1) + sum_s
                      s_diag_d (ig,2) = s_diag_d (ig,2) + sum_s
                   ENDDO
              END IF
           END DO
        ELSE
           DO na = 1, nat
              IF (ityp (na) == nt) THEN
                   ijkb_start = indv_ijkb0(na)
                   nh_ = nh(nt)
!$cuf kernel do(1) <<<*,*>>>
                   DO ig = 1, npw 
                      sum_h1 = 0.d0
                      sum_h4 = 0.d0
                      sum_s = 0.d0
                      
                      DO ih = 1, nh_
                         ikb = ijkb_start + ih
                         
                         ar = vkb_d (ig, ikb)*conjg(vkb_d (ig, ikb))
                         
                         sum_h1 = sum_h1 + dble(deeq_nc_d (ih, ih, na, 1) * ar)
                         sum_h4 = sum_h4 + dble(deeq_nc_d (ih, ih, na, 4) * ar)
                         sum_s = sum_s + dble(qq_at_d (ih, ih, na) * ar)
                      END DO
                      !
                      ! OPTIMIZE HERE : this scattered assign is bad!
                      !
                      h_diag_d (ig,1) = h_diag_d (ig,1) + sum_h1
                      h_diag_d (ig,2) = h_diag_d (ig,2) + sum_h4
                      s_diag_d (ig,1) = s_diag_d (ig,1) + sum_s
                      s_diag_d (ig,2) = s_diag_d (ig,2) + sum_s
                   ENDDO
              END IF
           END DO
        END IF
     END DO
  END SUBROUTINE usnldiag_noncollinear
  !
  SUBROUTINE usnldiag_spinorb()
     USE lsda_mod, ONLY: current_spin
     USE uspp_gpum,  ONLY: vkb_d, qq_at_d, qq_so_d, deeq_nc_d

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
        IF ( upf(nt)%tvanp .or.newpseudo (nt) ) THEN
           DO na = 1, nat
              IF (ityp (na) == nt) THEN
                   ijkb_start = indv_ijkb0(na)
                   nh_ = nh(nt)
!$cuf kernel do(1) <<<*,*>>>
                   DO ig = 1, npw   ! change this to 2*npw ?
                      sum_h1 = 0.d0
                      sum_h4 = 0.d0
                      sum_s1 = 0.d0
                      sum_s4 = 0.d0
                      
                      DO ih = 1, nh_
                         ikb = ijkb_start + ih
                         cv = vkb_d (ig, ikb)
                         DO jh = 1, nh_
                            jkb = ijkb_start + jh
                            
                            ar = cv*conjg(vkb_d (ig, jkb))
                            
                            sum_h1 = sum_h1 + dble(deeq_nc_d (ih, jh, na, 1) * ar)
                            sum_h4 = sum_h4 + dble(deeq_nc_d (ih, jh, na, 4) * ar)
                            sum_s1 = sum_s1  + dble(qq_so_d (ih, jh, 1, nt) * ar)
                            sum_s4 = sum_s4  + dble(qq_so_d (ih, jh, 4, nt) * ar)
                         END DO
                      END DO
                      !
                      ! OPTIMIZE HERE : this scattered assign is bad!
                      !
                      h_diag_d (ig,1) = h_diag_d (ig,1) + sum_h1
                      h_diag_d (ig,2) = h_diag_d (ig,2) + sum_h4
                      s_diag_d (ig,1) = s_diag_d (ig,1) + sum_s1
                      s_diag_d (ig,2) = s_diag_d (ig,2) + sum_s4
                   ENDDO
              END IF
           END DO
        ELSE
           DO na = 1, nat
              IF (ityp (na) == nt) THEN
                   ijkb_start = indv_ijkb0(na)
                   nh_ = nh(nt)
!$cuf kernel do(1) <<<*,*>>>
                   DO ig = 1, npw 
                      sum_h1 = 0.d0
                      sum_h4 = 0.d0
                      sum_s1 = 0.d0
                      sum_s4 = 0.d0
                      
                      DO ih = 1, nh_
                         ikb = ijkb_start + ih
                         
                         ar = vkb_d (ig, ikb)*conjg(vkb_d (ig, ikb))
                         
                         sum_h1 = sum_h1 + dble(deeq_nc_d (ih, ih, na, 1) * ar)
                         sum_h4 = sum_h4 + dble(deeq_nc_d (ih, ih, na, 4) * ar)
                         sum_s1 = sum_s1 + dble(qq_so_d (ih, ih, 1, nt) * ar)
                         sum_s4 = sum_s4 + dble(qq_so_d (ih, ih, 4, nt) * ar)
                      END DO
                      !
                      ! OPTIMIZE HERE : this scattered assign is bad!
                      !
                      h_diag_d (ig,1) = h_diag_d (ig,1) + sum_h1
                      h_diag_d (ig,2) = h_diag_d (ig,2) + sum_h4
                      s_diag_d (ig,1) = s_diag_d (ig,1) + sum_s1
                      s_diag_d (ig,2) = s_diag_d (ig,2) + sum_s4
                   ENDDO
              END IF
           END DO
        END IF
     END DO
  END SUBROUTINE usnldiag_spinorb
END SUBROUTINE usnldiag_gpu
