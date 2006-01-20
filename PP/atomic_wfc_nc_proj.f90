!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!-----------------------------------------------------------------------
SUBROUTINE atomic_wfc_nc_proj (ik, wfcatom)
  !-----------------------------------------------------------------------
  !
  ! This routine computes the  superposition of atomic wavefunctions for a
  ! given k-point in the non-collinear case.
  ! If lspinorb=.TRUE. it generates eigenstates of the total angular 
  ! momentum j, else eigenstates of l, m and s_z (the z-component of the spin) 
  ! are generated.  
  !
  USE kinds,      ONLY : DP
  USE parameters, ONLY : nchix
  USE atom,       ONLY : nchi, lchi, jchi, chi, oc, r, rab, msh
  USE constants,  ONLY : tpi, fpi, pi
  USE cell_base,  ONLY : omega, tpiba
  USE ions_base,  ONLY : nat, ntyp => nsp, ityp, tau
  USE basis,      ONLY : natomwfc
  USE gvect,      ONLY : ig1, ig2, ig3, eigts1, eigts2, eigts3, g
  USE klist,      ONLY : xk
  USE wvfct,      ONLY : npwx, npw, nbndx, nbnd, igk
  USE us,         ONLY : tab_at, dq
  USE noncollin_module, ONLY : noncolin, npol
  USE spin_orb,   ONLY : lspinorb, so, rot_ylm, fcoef, lmaxx
  !
  IMPLICIT NONE
  !
  INTEGER :: ik
  ! input: k-point
  COMPLEX(DP) :: wfcatom (npwx, npol, natomwfc) ! output: atomic wavefunctions
  !
  INTEGER :: n_starting_wfc, lmax_wfc, nt, l, nb, na, m, lm, ig, iig, &
             i0, i1, i2, i3
  !
  REAL(DP), ALLOCATABLE :: qg(:), ylm (:,:), chiq (:,:,:), aux (:), &
       gk (:,:), vchi (:)
  COMPLEX(DP), ALLOCATABLE :: sk (:), aux_so(:)
  COMPLEX(DP) :: kphase, lphase, scalar
  REAL(DP) :: vqint, arg, px, ux, vx, wx, alpha, gamman, &
                   fact(2), j, spinor
  INTEGER ::  ind, ind1, n1, n2, is, sph_ind

  CALL start_clock ('atomic_wfc')
  IF (.NOT.noncolin) CALL errore('atomic_wfc_nc_proj', &
                                 'called in the wrong case',1)

  ALLOCATE ( qg(npw), chiq(npw,nchix,ntyp), gk(3,npw), sk(npw))
  IF (lspinorb) ALLOCATE(aux_so(npw))

  ! calculate max angular momentum required in wavefunctions
  lmax_wfc = 0
  DO nt = 1, ntyp
     DO nb = 1, nchi (nt)
        lmax_wfc = max (lmax_wfc, lchi (nb, nt) )
     ENDDO
  ENDDO
  !
  ALLOCATE(ylm (npw,(lmax_wfc+1)**2) )
  !
  DO ig = 1, npw
     gk (1,ig) = xk(1, ik) + g(1, igk(ig) )
     gk (2,ig) = xk(2, ik) + g(2, igk(ig) )
     gk (3,ig) = xk(3, ik) + g(3, igk(ig) )
     qg(ig) = gk(1, ig)**2 +  gk(2, ig)**2 + gk(3, ig)**2
  ENDDO
  !
  !  ylm = spherical harmonics
  !
  CALL ylmr2 ((lmax_wfc+1)**2, npw, gk, qg, ylm)
  !
  ! set now q=|k+G| in atomic units
  !
  DO ig = 1, npw
     qg(ig) = sqrt(qg(ig))*tpiba
  ENDDO
  !
  n_starting_wfc = 0
  !
  ! chiq = radial fourier transform of atomic orbitals chi
  !
  DO nt = 1, ntyp
     DO nb = 1, nchi (nt)
        IF ( oc (nb, nt) >= 0.d0) THEN
           DO ig = 1, npw
              px = qg (ig) / dq - int (qg (ig) / dq)
              ux = 1.d0 - px
              vx = 2.d0 - px
              wx = 3.d0 - px
              i0 = qg (ig) / dq + 1
              i1 = i0 + 1
              i2 = i0 + 2
              i3 = i0 + 3
              chiq (ig, nb, nt) = &
                     tab_at (i0, nb, nt) * ux * vx * wx / 6.d0 + &
                     tab_at (i1, nb, nt) * px * vx * wx / 2.d0 - &
                     tab_at (i2, nb, nt) * px * ux * wx / 2.d0 + &
                     tab_at (i3, nb, nt) * px * ux * vx / 6.d0
           ENDDO
        ENDIF
     ENDDO
  ENDDO
  alpha = 0.d0
  gamman = 0.5d0*pi
  DO na = 1, nat
     arg = (xk(1,ik)*tau(1,na) + xk(2,ik)*tau(2,na) + xk(3,ik)*tau(3,na)) * tpi
     kphase = CMPLX (cos (arg), - sin (arg) )
     !
     !     sk is the structure factor
     !
     DO ig = 1, npw
        iig = igk (ig)
        sk (ig) = kphase * eigts1 (ig1 (iig), na) * eigts2 (ig2 (iig), na) * &
                           eigts3 (ig3 (iig), na)
     ENDDO
     !
     nt = ityp (na)
     DO nb = 1, nchi (nt)
        IF (oc (nb, nt) >= 0.d0) THEN
           l = lchi (nb, nt)
           lphase = (0.d0,1.d0)**l
           !  the factor i^l MUST BE PRESENT in order to produce
           !  wavefunctions for k=0 that are real in real space
           IF (lspinorb) THEN
             IF (so(nt)) THEN
               j = jchi (nb, nt)
               DO m = -l-1, l
                 fact(1) = spinor(l,j,m,1)
                 fact(2) = spinor(l,j,m,2)
                 IF (abs(fact(1)).gt.1.d-8.or.abs(fact(2)).gt.1.d-8) THEN
                    n_starting_wfc = n_starting_wfc + 1
                    IF (n_starting_wfc.gt.natomwfc) &
                       CALL errore ('atomic_wfc_nc_proj', 'too many wfcs', 1)
                    DO is=1,2
                       IF (abs(fact(is)).gt.1.d-8) THEN
                          ind=lmaxx+1+sph_ind(l,j,m,is)
                          aux_so=(0.d0,0.d0)
                          DO n1=1,2*l+1
                             ind1=l**2+n1
                             IF (abs(rot_ylm(ind,n1)).gt.1.d-8) &
                             aux_so(:)=aux_so(:)+rot_ylm(ind,n1)*ylm(:,ind1)
                          ENDDO
                          DO ig=1,npw
                             wfcatom (ig,is,n_starting_wfc) = lphase*fact(is)*&
                             sk(ig)*aux_so(ig)*chiq (ig, nb, nt)
                          ENDDO
                       ENDIF
                    ENDDO
                 ENDIF
               ENDDO
             ELSE
               DO n2 = l, l + 1
                 j= DBLE(n2) - 0.5d0
                 IF (j.gt.0.d0)  THEN                          
                   DO m = -l-1, l
                     fact(1) = spinor(l,j,m,1)
                     fact(2) = spinor(l,j,m,2)
                     IF (abs(fact(1)).gt.1.d-8.or.abs(fact(2)).gt.1.d-8) THEN
                       n_starting_wfc = n_starting_wfc + 1
                       IF (n_starting_wfc.gt.natomwfc) &
                         CALL errore ('atomic_wfc_nc_proj', 'too many wfcs', 1)
                       DO is=1,2
                         IF (abs(fact(is)).gt.1.d-8) THEN
                           ind=lmaxx+1+sph_ind(l,j,m,is)
                           aux_so=(0.d0,0.d0)
                           DO n1=1,2*l+1
                             ind1=l**2+n1
                             IF (abs(rot_ylm(ind,n1)).gt.1.d-8) &
                             aux_so(:)=aux_so(:)+rot_ylm(ind,n1)*ylm(:,ind1)
                           ENDDO
                           DO ig=1,npw
                             wfcatom (ig,is,n_starting_wfc) = lphase*fact(is)*&
                             sk(ig)*aux_so(ig)*chiq (ig, nb, nt)
                           ENDDO
                         ENDIF
                       ENDDO
                     ENDIF
                   ENDDO 
                 ENDIF
               ENDDO 
             ENDIF 
           ELSE
             DO m = 1, 2 * l + 1
                lm = l**2 + m
                n_starting_wfc = n_starting_wfc + 1
                IF (n_starting_wfc.GT.natomwfc) &
                   CALL errore ('atomic_wfc_nc_proj', 'too many wfcs', 1)
                IF (n_starting_wfc+2*l+1 .GT. nbndx) &
                   CALL errore('atomic_wfc_nc_proj','too many wfcs',1)
                DO ig=1,npw
                   scalar = sk(ig)*ylm(ig,lm)*chiq(ig,nb,nt)
                   wfcatom(ig,1,n_starting_wfc) = scalar
                   wfcatom(ig,2,n_starting_wfc) = (0.d0,0.d0)

                   wfcatom(ig,1,n_starting_wfc+2*l+1) =(0.d0,0.d0)
                   wfcatom(ig,2,n_starting_wfc+2*l+1) = scalar
                ENDDO
             ENDDO
             n_starting_wfc = n_starting_wfc + 2*l+1
           ENDIF
        ENDIF
     ENDDO
  ENDDO

  IF (n_starting_wfc.ne.natomwfc) CALL errore ('atomic_wfc_nc_proj', &
       'something wrong', 1)

  DEALLOCATE(qg, chiq ,gk, sk ,ylm)
  IF (lspinorb) DEALLOCATE(aux_so)

  CALL stop_clock ('atomic_wfc')
  RETURN
END SUBROUTINE atomic_wfc_nc_proj
