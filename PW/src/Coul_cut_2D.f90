!
! Copyright (C) 2007-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This module contains the variables and routines necessary to the implementation
! of the two-dimensional Coulomb cutoff. Details of the implementation can be found in:
!
! Sohier, T., Calandra, M., & Mauri, F. (2017), 
! "Density functional perturbation theory for gated two-dimensional heterostructures: 
! Theoretical developments and application to flexural phonons in graphene." 
! Physical Review B, 96(7), 75448. https://doi.org/10.1103/PhysRevB.96.075448
!
!----------------------------------------------------------------------------
MODULE Coul_cut_2D
  !----------------------------------------------------------------------------
  !! This module contains the variables and subroutines needed for the 
  !! 2D Coulomb cutoff.
  !
  USE kinds,       ONLY : DP
  USE constants,   ONLY : tpi, pi
  !
  SAVE
  !
  LOGICAL :: do_cutoff_2D = .FALSE.
  !! flag for 2D cutoff. If true, the cutoff is active.
  REAL(DP) :: lz
  !! The distance in the out-plne direction after which potential are cut off.
  !
  REAL(DP), ALLOCATABLE :: cutoff_2D(:)
  !! The factor appended to the Coulomb Kernel to cut off potentials
  REAL(DP), ALLOCATABLE :: lr_Vloc(:,:)
  !! The long-range part of the local part of the ionic potential
  !
CONTAINS
!
!----------------------------------------------------------------------
SUBROUTINE cutoff_fact()
  !----------------------------------------------------------------------
  !! This routine calculates the cutoff factor in G-space and stores it in 
  !! a vector called \(\text{cutoff_2D}(:)\), to be re-used in various routines.  
  !! See Eq.(24) of PRB 96, 075448
  !
  USE kinds
  USE io_global,    ONLY : stdout
  USE gvect,        ONLY : g, ngm, ngmx
  USE cell_base,    ONLY : alat, celldm, at
  !
  IMPLICIT NONE
  !
  ! ... local variables
  !
  INTEGER :: ng, i
  ! counter over G vectors, cartesian coord.
  REAL(DP) :: Gzlz, Gplz
  !
  ALLOCATE( cutoff_2D(ngmx) ) 
  !
  ! Message to indicate that the cutoff is active. 
  WRITE(stdout, *) "----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D"
  WRITE(stdout, *) " The code is running with the 2D cutoff"
  WRITE(stdout, *) " Please refer to:"
  WRITE(stdout, *) " Sohier, T., Calandra, M., & Mauri, F. (2017), "
  WRITE(stdout, *) " Density functional perturbation theory for gated two-dimensional heterostructures:" 
  WRITE(stdout, *) " Theoretical developments and application to flexural phonons in graphene." 
  WRITE(stdout, *) " Physical Review B, 96(7), 75448. https://doi.org/10.1103/PhysRevB.96.075448"
  WRITE(stdout, *) "----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D"
  !  at(:,i) are the lattice vectors of the simulation cell, a_i,
  !  in alat units: a_i(:) = at(:,i)/alat
  !  Check that material is in the x-y plane
  DO i = 1, 2
     IF (ABS(at(3,i))>1d-8) WRITE(stdout, *) "2D CODE WILL NOT WORK, 2D MATERIAL NOT IN X-Y PLANE!!"
  ENDDO
  ! define cutoff distnce and compute cutoff factor
  lz = 0.5d0*at(3,3)*alat
  DO ng = 1, ngm
     Gplz = SQRT( g(1,ng)**2 + g(2,ng)**2 )*tpi*lz/alat
     Gzlz = g(3,ng)*tpi*lz/alat
     cutoff_2D(ng) = 1.0d0 - EXP(-Gplz)*COS(Gzlz)
  ENDDO
  !
  RETURN
  !
END SUBROUTINE cutoff_fact
!
!----------------------------------------------------------------------
SUBROUTINE cutoff_lr_Vloc( )
  !----------------------------------------------------------------------
  !! This routine calculates the long-range part of \(\text{vloc}(g)\) for
  !! 2D calculations.  
  !! See Eq. (32) of PRB 96, 075448.
  !
  USE kinds
  USE constants,    ONLY : fpi, e2, eps8
  USE fft_base,     ONLY : dfftp
  USE gvect,        ONLY : ngm, gg, g, ngmx
  ! gg is G^2 in increasing order (in units of tpiba2=(2pi/a)^2)
  USE ions_base,    ONLY : zv, nsp
  USE uspp_param,   ONLY : upf
  USE cell_base,    ONLY : omega, tpiba2
  !
  IMPLICIT NONE
  !
  ! ... local variables
  !
  INTEGER :: ng, nt, ng0 
  REAL(DP) ::fac
  !
  IF (.NOT. ALLOCATED(lr_Vloc)) ALLOCATE( lr_Vloc(ngmx,nsp) )
  !
  lr_Vloc(:,:) = 0.0d0
  ! set G=0 value to zero
  IF (gg(1)<eps8) THEN
     lr_Vloc(1,:) = 0.0d0
     ng0 = 2
  ELSE
     ng0 = 1
  ENDIF
  ! Set g.neq.0 values
  DO nt = 1, nsp
     fac = upf(nt)%zp * e2 / tpiba2
     DO ng = ng0, ngm
        lr_Vloc(ng,nt) = - fpi / omega* fac * cutoff_2D(ng)* &
                       & EXP( -gg(ng) * tpiba2 * 0.25d0) / gg(ng)
     ENDDO
  ENDDO
  !
  RETURN
  !
END SUBROUTINE cutoff_lr_Vloc
!
!----------------------------------------------------------------------
SUBROUTINE cutoff_local( aux )
  !----------------------------------------------------------------------
  !! This subroutine is called to re-add the long-range part of the local
  !! part of the ionic potential, using \(\text{lr_Vloc}\) computed in
  !! routine \(\texttt{cutoff_lr_Vloc}\).  
  !! See Eq. (33) of PRB 96, 075448
  !
  USE kinds
  USE fft_base,   ONLY : dfftp
  USE gvect,      ONLY : ngm
  USE vlocal,     ONLY : strf
  USE ions_base,  ONLY : nsp
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(INOUT):: aux(dfftp%nnr)
  !! input: local part of ionic potential 
  !
  ! ... local variables
  !
  INTEGER :: ng, nt 
  !
  DO nt = 1, nsp
     DO ng = 1, ngm
        aux(dfftp%nl(ng)) = aux(dfftp%nl(ng)) + lr_Vloc(ng,nt) * strf(ng,nt)
     ENDDO
  ENDDO
  !
  RETURN
  !
END SUBROUTINE cutoff_local
!
!
!----------------------------------------------------------------------
SUBROUTINE cutoff_hartree( rhog, aux1, ehart )
  !----------------------------------------------------------------------
  !! This subroutine cuts off the Hartree potential and defines Hartree
  !! energy accordingly in G-space.  
  !! See Eq. (34) and (41) of PRB 96, 075448
  !
  USE kinds
  USE gvect,       ONLY : ngm, gg , gstart
  USE io_global,   ONLY : stdout
  !
  IMPLICIT NONE
  ! 
  COMPLEX(DP), INTENT(IN) :: rhog(ngm)
  !! local potential
  REAL(DP), INTENT(INOUT) :: aux1(2,ngm)
  !! Hartree potential
  REAL(DP), INTENT(INOUT) :: ehart
  !! Hartree energy
  !
  ! ... local variables
  !
  INTEGER :: ig
  REAL(DP) :: fac
  REAL(DP) :: rgtot_re, rgtot_im
  !  
  DO ig = gstart, ngm
     !
     fac = 1.D0 / gg(ig) * cutoff_2D(ig)
     !
     rgtot_re = REAL( rhog(ig) )
     rgtot_im = AIMAG( rhog(ig) )
     !
     ehart = ehart + ( rgtot_re**2 + rgtot_im**2 ) * fac
     !
     aux1(1,ig) = rgtot_re * fac
     aux1(2,ig) = rgtot_im * fac
     !
  ENDDO
  !
  RETURN
  !
END SUBROUTINE cutoff_hartree
!
!
!----------------------------------------------------------------------
SUBROUTINE cutoff_ewald( alpha, ewaldg, omega )
  !----------------------------------------------------------------------
  !! This subroutine defines computes the cutoff version of the 
  !! Ewald sum in G space.
  !! See Eq. (46) of PRB 96, 075448
  !
  USE kinds
  USE gvect,      ONLY : ngm, gg, gstart
  USE ions_base,  ONLY : zv, nsp, nat, ityp
  USE cell_base,  ONLY : tpiba2, alat
  USE vlocal,     ONLY : strf
  USE io_global,  ONLY : stdout
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: alpha
  !! tuning parameter for Ewald LR/SR separation
  REAL(DP), INTENT(INOUT) :: ewaldg
  !! Ewald sum
  REAL(DP), INTENT(IN) :: omega
  !! unit-cell volume
  !
  ! ... local variables
  !
  INTEGER :: ng, nt, na, nr, ir, iz, nz, rmax
  COMPLEX(DP) :: rhon
  REAL(DP) :: rp, z
  !
  ! The G=0 component of the long-ranged local part of the 
  ! pseudopotential minus the Hartree potential is set to 0.
  ! This is equivalent to substracting the finite non-singular
  ! part of the ionic potential at G=0. See Appendix D.2 of PRB 96, 075448.
  ! In practice, with respect to the 3D ewald sum, we must subtract the 
  ! G=0 energy term - 4 pi/omega * e**2/2 * charge**2 / alpha / 4.0d0.
  ! That is, we simply set setting ewaldg(G=0)=0.
  !
  ewaldg = 0.0d0
  ! now the G.neq.0 terms
  DO ng = gstart, ngm
     rhon = (0.d0, 0.d0)
     DO nt = 1, nsp
        rhon = rhon + zv(nt) * CONJG(strf(ng,nt) )
     ENDDO
     ewaldg = ewaldg +  ABS(rhon)**2 * EXP( - gg(ng) * tpiba2 /&
              alpha / 4.d0) / gg(ng)*cutoff_2D(ng) / tpiba2
  ENDDO
  ewaldg = 2.d0 * tpi / omega * ewaldg
  !
  ! ... here add the other constant term (Phi_self)
  !
  IF (gstart == 2) THEN
     DO na = 1, nat
        ewaldg = ewaldg - zv(ityp(na))**2 * SQRT(8.d0 / tpi * &
                 alpha)
     ENDDO
  ENDIF
  !  
  RETURN
  !
END SUBROUTINE cutoff_ewald
!
!
!----------------------------------------------------------------------
SUBROUTINE cutoff_force_ew( aux, alpha )
  !----------------------------------------------------------------------
  !! This subroutine cuts off the Ewald contribution to the forces. More
  !! precisely, it cuts off the LR ion-ion potential that is then used to 
  !! compute the Ewald forces.  
  !! See Eq. (55) of PRB 96, 075448 (note that Eq. (56), derived from Eq. (55),
  !! looks somewhat different from what is implemented in the code, but it is
  !! equivalent).
  !
  USE kinds
  USE gvect,        ONLY : ngm, gg , gstart
  USE cell_base,    ONLY : tpiba2, alat
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(INOUT) :: aux(ngm)
  !! long-range part of the ionic potential
  REAL(DP), INTENT(IN) :: alpha
  !! tuning parameter for the LR/SR separation
  !
  ! ... local variables
  !
  INTEGER :: ig
  !
  DO ig = gstart, ngm
     aux(ig) = aux(ig) * EXP( - gg(ig) * tpiba2 / alpha / 4.d0) &
               / (gg(ig) * tpiba2) * cutoff_2D(ig)
  ENDDO  
  !
  RETURN
  !
END SUBROUTINE cutoff_force_ew
!
!
!----------------------------------------------------------------------
SUBROUTINE cutoff_force_lc( aux, forcelc )
  !----------------------------------------------------------------------
  !! This subroutine re-adds the cutoff contribution from the long-range 
  !! local part of the ionic potential to the forces. In the 2D code, this 
  !! contribution is missing from the \(\text{Vloc}\).  
  !! See Eq. (54) of PRB 96, 075448.
  !
  USE kinds
  USE gvect,         ONLY : ngm, gg, g , gstart
  USE constants,     ONLY : fpi, e2, eps8, tpi
  USE uspp_param,    ONLY : upf 
  USE cell_base,     ONLY : tpiba2, alat, omega
  USE ions_base,     ONLY : nat, zv, tau, ityp
  USE io_global,     ONLY : stdout
  USE fft_base,      ONLY : dfftp
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(IN) :: aux(dfftp%nnr)
  !! local ionic potential
  REAL(DP), INTENT(INOUT) :: forcelc(3,nat)
  !! corresponding force contribution
  !
  ! ... local variables
  !
  REAL(DP) :: arg
  INTEGER :: ig, na, ipol
  !
  DO na = 1, nat
     DO ig = gstart, ngm 
        arg = (g(1,ig) * tau(1,na) + g(2,ig) * tau(2,na) + &
               g(3,ig) * tau(3,na) ) * tpi
        DO ipol = 1, 3
           forcelc(ipol,na) = forcelc (ipol,na) + tpi / alat * &
                 g(ipol,ig) * lr_Vloc(ig, ityp(na)) * omega  * &
                ( SIN(arg)*DBLE(aux(dfftp%nl(ig))) + COS(arg)*AIMAG(aux(dfftp%nl(ig))) )
        ENDDO
     ENDDO 
  ENDDO
  !
  RETURN
  !
END SUBROUTINE cutoff_force_lc
!
!
!----------------------------------------------------------------------
SUBROUTINE cutoff_stres_evloc( psic_G, strf, evloc )
  !----------------------------------------------------------------------
  !! This subroutine adds the contribution from the cutoff long-range part
  !! of the local part of the ionic potential to \(\text{evloc}\).  
  !! evloc corresponds to the delta term in Eq. (63) of PRB 96, 075448.
  !! It is the energy of the electrons in the local ionic potential.  
  !! Note that it is not calculated as such (by itself) in the standard code.
  !! Indeed, it is "hidden" in the sum of KS eigenvalues. That is why we need 
  !! to re-compute it here for the stress.
  !
  USE kinds
  USE ions_base,  ONLY : ntyp => nsp
  USE gvect,      ONLY : ngm , gstart
  USE io_global,  ONLY : stdout
  USE fft_base,   ONLY : dfftp
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(IN) :: psic_G(dfftp%nnr)
  !! charge density in G space
  COMPLEX(DP), INTENT(IN) :: strf(ngm,ntyp)
  !! the structure factor
  REAL(DP), INTENT(INOUT) :: evloc
  !! the energy of the electrons in the local ionic potential
  !
  ! ... local variables
  !
  INTEGER :: ng, nt
  !
  ! If gstart=2, it means g(1) is G=0, but we have nothing to add for G=0
  ! So we start at gstart.
  DO nt = 1, ntyp
     DO ng = gstart, ngm
        evloc = evloc + DBLE( CONJG(psic_G(dfftp%nl(ng))) * strf(ng,nt) ) &
                        * lr_Vloc(ng,nt) 
     ENDDO
  ENDDO
  !
  RETURN
  !
END SUBROUTINE cutoff_stres_evloc
!
!
!----------------------------------------------------------------------
SUBROUTINE cutoff_stres_sigmaloc( psic_G, strf, sigmaloc )
  !----------------------------------------------------------------------
  !! This subroutine adds the contribution from the cutoff long-range part 
  !! of the local part of the ionic potential to the rest of the 
  !! \(\text{sigmaloc}\). That is, the rest of Eq. (63) of PRB 96, 075448.
  !
  USE kinds
  USE ions_base,   ONLY : ntyp => nsp
  USE constants,   ONLY : eps8
  USE gvect,       ONLY : ngm, g, gg, gstart
  USE cell_base,   ONLY : tpiba, tpiba2, alat, omega
  USE io_global,   ONLY : stdout
  USE fft_base,    ONLY : dfftp
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(IN) :: psic_G(dfftp%nnr)
  !! charge density in G space
  COMPLEX(DP), INTENT(IN) :: strf(ngm,ntyp)
  !! the structure factor
  REAL(DP), INTENT(INOUT) :: sigmaloc(3,3)
  !! stress contribution for the local ionic potential
  !
  ! ... local variables
  !
  INTEGER :: ng, nt, l, m
  REAL(DP) :: Gp, G2lzo2Gp, beta, dlr_Vloc
  !
  ! no G=0 contribution
  DO nt = 1, ntyp
     DO ng = gstart, ngm
        !
        Gp = SQRT( g(1,ng)**2 + g(2,ng)**2 )*tpiba
        ! below is a somewhat cumbersome way to define beta of Eq. (61) of PRB 96, 075448
        IF (Gp < eps8) THEN
           ! G^2*lz/2|Gp|
           G2lzo2Gp = 0.0d0
           beta = 0.0d0
        ELSE
           G2lzo2Gp = gg(ng)*tpiba2*lz/2.0d0/Gp
           beta = G2lzo2Gp*(1.0d0-cutoff_2D(ng))/cutoff_2D(ng)
        ENDIF
        ! dlr_vloc corresponds to the derivative of the long-range local ionic potential
        ! with respect to G
        DO l = 1, 3
           IF (l == 3) THEN
              dlr_Vloc = - 1.0d0/ (gg(ng)*tpiba2) * lr_Vloc(ng,nt)  &
                               * (1.0d0+ gg(ng)*tpiba2/4.0d0)
           ELSE
              dlr_Vloc = - 1.0d0/ (gg(ng)*tpiba2) * lr_Vloc(ng,nt)  &
                               * (1.0d0- beta + gg(ng)*tpiba2/4.0d0)
           ENDIF
           !
           DO m = 1, l
              sigmaloc(l,m) = sigmaloc(l,m) +  DBLE( CONJG( psic_G(dfftp%nl(ng) ) ) &
                              * strf(ng,nt) ) * 2.0d0 * dlr_Vloc  &
                              * tpiba2 * g(l,ng) * g(m,ng) 
           ENDDO
        ENDDO
        !
     ENDDO
  ENDDO
  !
  RETURN
  !
END SUBROUTINE cutoff_stres_sigmaloc
!
!
!----------------------------------------------------------------------
SUBROUTINE cutoff_stres_sigmahar( psic_G, sigmahar )
  !----------------------------------------------------------------------
  !! This subroutine cuts off the Hartree part of the stress.  
  !! See Eq. (62) of PRB 96, 075448.
  !
  USE kinds
  USE gvect,      ONLY : ngm, g, gg, gstart
  USE constants,  ONLY : eps8
  USE cell_base,  ONLY : tpiba2, alat, tpiba
  USE io_global,  ONLY : stdout
  USE fft_base,   ONLY : dfftp
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(IN) :: psic_G(dfftp%nnr)
  !! charge density in G-space
  REAL(DP), INTENT(INOUT) :: sigmahar(3,3)
  !! hartree contribution to stress
  !
  ! ... local variables
  !
  INTEGER :: ng, nt, l, m
  REAL(DP) :: Gp, G2lzo2Gp, beta, shart, g2, fact
  !
  DO ng = gstart, ngm
     Gp = SQRT( g(1,ng)**2 + g(2,ng)**2 )*tpiba
     IF (Gp < eps8) THEN
        G2lzo2Gp = 0.0d0
        beta = 0.0d0
     ELSE
        G2lzo2Gp = gg(ng)*tpiba2*lz/2.0d0/Gp
        beta = G2lzo2Gp*(1.0d0-cutoff_2D(ng))/cutoff_2D(ng)
     ENDIF
     g2 = gg (ng) * tpiba2
     shart = psic_G(dfftp%nl(ng)) * CONJG(psic_G(dfftp%nl(ng))) / g2 * cutoff_2D(ng)
     DO l = 1, 3
        IF (l == 3) THEN
           fact = 1.0d0
        ELSE
           fact = 1.0d0 - beta
        ENDIF
        DO m = 1, l
           sigmahar(l,m) = sigmahar(l,m) + shart * tpiba2 * 2 * &
                           g(l,ng) * g(m,ng) / g2  * fact
        ENDDO
     ENDDO
  ENDDO
  !sigma is multiplied by 0.5*fpi*e2 after
  !
  RETURN
  !
END SUBROUTINE cutoff_stres_sigmahar
!
!
!----------------------------------------------------------------------
SUBROUTINE cutoff_stres_sigmaewa( alpha, sdewald, sigmaewa )
  !----------------------------------------------------------------------
  !! This subroutine cuts off the Ewald part of the stress.  
  !! See Eq. (64) in PRB 96 075448
  !
  USE kinds
  USE ions_base,   ONLY : nat, zv, tau, ityp
  USE constants,   ONLY : e2, eps8
  USE gvect,       ONLY : ngm, g, gg, gstart
  USE cell_base,   ONLY : tpiba2, alat, omega, tpiba
  USE io_global,   ONLY : stdout
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: alpha
  !! tuning param for LR/SR separation
  REAL(DP), INTENT(INOUT) :: sigmaewa(3,3)
  !! ewald contribution to stress
  REAL(DP), INTENT(INOUT) :: sdewald
  !! constant and diagonal terms
  !
  ! ... local variables
  !
  INTEGER :: ng, na, l, m
  REAL(DP) :: Gp, G2lzo2Gp, beta, sewald, g2, g2a, arg, fact
  COMPLEX(DP) :: rhostar
  !
  ! g(1) is a problem if it's G=0, because we divide by G^2. 
  ! So start at gstart.
  ! fact=1.0d0, gamma_only not implemented
  ! G=0 componenent of the long-range part of the local part of the 
  ! pseudopotminus the Hartree potential is set to 0.
  ! in other words, sdewald=0.  
  ! sdewald is the last term in equation B1 of PRB 32 3792.
  ! See also similar comment for ewaldg in cutoff_ewald routine
  !
  sdewald = 0.0D0
  DO ng = gstart, ngm
     Gp = SQRT( g(1,ng)**2 + g(2,ng)**2 )*tpiba
     IF (Gp < eps8) THEN
        G2lzo2Gp = 0.0d0
        beta = 0.0d0
     ELSE
        G2lzo2Gp = gg(ng)*tpiba2*lz/2.0d0/Gp
        beta = G2lzo2Gp*(1.0d0-cutoff_2D(ng))/cutoff_2D(ng)
     ENDIF
     g2 = gg(ng) * tpiba2
     g2a = g2 / 4.d0 / alpha
     rhostar = (0.d0, 0.d0)
     DO na = 1, nat
        arg = (g(1,ng) * tau(1,na) + g(2,ng) * tau(2,na) + &
               g(3,ng) * tau(3,na) ) * tpi
        rhostar = rhostar + zv (ityp(na) ) * CMPLX(COS(arg), SIN(arg), KIND=DP)
     ENDDO
     rhostar = rhostar / omega
     sewald = tpi * e2 * EXP(-g2a) / g2* cutoff_2D(ng) * ABS(rhostar)**2
     ! ... sewald is an other diagonal term that is similar to the diagonal terms 
     ! in the other stress contributions. It basically gives a term prop to 
     ! the ewald energy
     sdewald = sdewald-sewald
     DO l = 1, 3
        IF (l == 3) THEN
           fact = (g2a + 1.0d0)
        ELSE
           fact = (1.0d0+g2a-beta)
        ENDIF
        !
        DO m = 1, l
           sigmaewa(l,m) = sigmaewa(l,m) + sewald * tpiba2 * 2.d0 * &
                 g(l,ng) * g(m,ng) / g2 * fact
        ENDDO
     ENDDO
     !
  ENDDO
  !
  RETURN
  !
END SUBROUTINE cutoff_stres_sigmaewa
!
!
END MODULE Coul_cut_2D
