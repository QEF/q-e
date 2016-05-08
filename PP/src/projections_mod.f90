!
! Copyright (C) 2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE projections
  USE kinds, ONLY : DP
  
  TYPE wfc_label
     INTEGER na, n, l, m, ind
     REAL (DP) jj
  END TYPE wfc_label
  TYPE(wfc_label), ALLOCATABLE :: nlmchi(:)
  
  REAL (DP),    ALLOCATABLE :: proj (:,:,:)
  COMPLEX (DP), ALLOCATABLE :: proj_aux (:,:,:)
  COMPLEX (DP), ALLOCATABLE :: ovps_aux (:,:,:)
  
  CONTAINS
    !
    SUBROUTINE fill_nlmchi ( natomwfc, nwfc, lmax_wfc )
      !
      USE ions_base, ONLY : ityp, nat
      USE uspp_param, ONLY: upf
      USE spin_orb, ONLY: lspinorb
      USE noncollin_module, ONLY: noncolin
      !
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: natomwfc
      INTEGER, INTENT (OUT) :: nwfc, lmax_wfc 
      !
      INTEGER :: na, nt, n, n1, n2, l, m, ind
      REAL(dp) :: jj, fact(2)
      REAL(dp), EXTERNAL :: spinor
      !
      ALLOCATE (nlmchi(natomwfc))
      nwfc=0
      lmax_wfc = 0
      DO na = 1, nat
         nt = ityp (na)
         n2 = 0
         DO n = 1, upf(nt)%nwfc
            IF (upf(nt)%oc (n) >= 0.d0) THEN
               l = upf(nt)%lchi (n)
               lmax_wfc = max (lmax_wfc, l )
               IF (lspinorb) THEN
                  IF (upf(nt)%has_so) THEN
                     jj = upf(nt)%jchi (n)
                     ind = 0
                     DO m = -l-1, l
                        fact(1) = spinor(l,jj,m,1)
                        fact(2) = spinor(l,jj,m,2)
                        IF (abs(fact(1)) > 1.d-8 .or. abs(fact(2)) > 1.d-8) THEN
                           nwfc = nwfc + 1
                           ind = ind + 1
                           nlmchi(nwfc)%na = na
                           nlmchi(nwfc)%n  =  n
                           nlmchi(nwfc)%l  =  l
                           nlmchi(nwfc)%m  =  m
                           nlmchi(nwfc)%ind  =  ind
                           nlmchi(nwfc)%jj  =  jj
                        ENDIF
                     ENDDO
                  ELSE
                     DO n1 = l, l+1
                        jj= dble(n1) - 0.5d0
                        ind = 0
                        IF (jj>0.d0)  THEN
                           n2 = n2 + 1
                           DO m = -l-1, l
                              fact(1) = spinor(l,jj,m,1)
                              fact(2) = spinor(l,jj,m,2)
                              IF (abs(fact(1)) > 1.d-8 .or. abs(fact(2)) > 1.d-8) THEN
                                 nwfc = nwfc + 1
                                 ind = ind + 1
                                 nlmchi(nwfc)%na = na
                                 nlmchi(nwfc)%n  =  n2
                                 nlmchi(nwfc)%l  =  l
                                 nlmchi(nwfc)%m  =  m
                                 nlmchi(nwfc)%ind  =  ind
                                 nlmchi(nwfc)%jj  =  jj
                              ENDIF
                           ENDDO
                        ENDIF
                     ENDDO
                  ENDIF
               ELSE
                  DO m = 1, 2 * l + 1
                     nwfc=nwfc+1
                     nlmchi(nwfc)%na = na
                     nlmchi(nwfc)%n  =  n
                     nlmchi(nwfc)%l  =  l
                     nlmchi(nwfc)%m  =  m
                     nlmchi(nwfc)%ind=  m
                     nlmchi(nwfc)%jj =  0.d0
                  ENDDO
                  IF ( noncolin) THEN
                     DO m = 1, 2 * l + 1
                        nwfc=nwfc+1
                        nlmchi(nwfc)%na = na
                        nlmchi(nwfc)%n  =  n
                        nlmchi(nwfc)%l  =  l
                        nlmchi(nwfc)%m  =  m
                        nlmchi(nwfc)%ind=  m+2*l+1
                        nlmchi(nwfc)%jj =  0.d0
                     END DO
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
      ENDDO
      !
      IF (lmax_wfc > 3) CALL errore ('fill_nlmchi', 'l > 3 not yet implemented',1)
      IF (nwfc /= natomwfc) CALL errore ('fill_nlmchi','wrong # of atomic wfcs',1)
      
    END SUBROUTINE fill_nlmchi
    !
END MODULE projections
