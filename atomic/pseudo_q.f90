! Copyright (C) 2007 QUANTUM-espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This routine derives from the routine us2paw in atomic_paw.
!
SUBROUTINE pseudo_q (qfunc, qfuncl)       
USE kinds, ONLY : DP
USE constants, ONLY : eps8
USE ld1_parameters, ONLY : nwfsx
USE ld1inc, ONLY : which_augfun, rmatch_augfun, lls, &
                   grid, ndmx, lmx2, nbeta, ikk
USE atomic_paw, ONLY : find_bes_qi
IMPLICIT NONE
!
REAL(DP), INTENT(IN) :: qfunc(ndmx,nwfsx,nwfsx)
REAL(DP), INTENT(OUT) :: qfuncl(ndmx,nwfsx,nwfsx,0:lmx2) 
REAL(DP),  EXTERNAL :: int_0_inf_dr
!
! variables for aug. functions generation
! 
INTEGER  :: irc, ns, ns1, l, l1, l2, l3, lll, mesh, ircm, ir, nc, iok 
REAL(DP) :: aux(ndmx), raux
REAL(DP) :: qc(2), xc(2), b1(2), b2(2), energy(5,3), twosigma2, rm, gaussian(ndmx)  
REAL(DP), ALLOCATABLE :: j1(:,:), augmom(:,:,:)
!
ALLOCATE(augmom(nbeta,nbeta,0:lmx2))

mesh = grid%mesh
irc = maxval(ikk(1:nbeta))+8
!
qfuncl = 0.0_dp
augmom(:,:,:) = 0.0_dp
!
!   Compute the multipoles of the function
!    
DO ns=1,nbeta
   l1=lls(ns)
   DO ns1=ns,nbeta
      l2=lls(ns1)
      DO l3 = abs(l1-l2), l1+l2
         aux(1:irc) = qfunc(1:irc,ns,ns1) * grid%r(1:irc)**l3
         lll = l1 + l2 + 2 + l3
         augmom(ns,ns1,l3)=int_0_inf_dr(aux(1:irc),grid,irc,lll)
         augmom(ns1,ns,l3)=augmom(ns,ns1,l3)
      END DO
      WRITE (*,'(a,2i3,a,2i3,a,i3,a,i3,10f8.4)') " MULTIPOLE",ns,l1,":",ns1,l2,&
                   " l3=",ABS(l1-l2)," - ",l1+l2, &
                         (augmom(ns,ns1,l3), l3=abs(l1-l2),l1+l2)
   END DO
END DO
!
CALL infomsg ('pseudo_q','You have specified: '//TRIM(which_augfun))  

DO ns=1,nbeta
   l1 = lls(ns)
   DO ns1=ns,nbeta
      l2 = lls(ns1)
        !
      SELECT CASE (which_augfun)
      CASE ('QVAN')
!
!   This is the original case. The qvan are generated from norm conserving
!   wavefunctions
!
         qfuncl(1:mesh,ns,ns1,0) = qfunc(1:mesh,ns,ns1)
      CASE ('GAUSS')
         ! define a "gaussian" (with a volume element)
         rm = rmatch_augfun
         twosigma2 = 2.0_DP*(rm/3.0_dp)**2
         DO ir=1,mesh
            IF (grid%r(ir) <= rm) THEN
               gaussian(ir) = ( EXP(-grid%r(ir)**2/twosigma2) + &
                                EXP(-(grid%r(ir)-2.0_DP*rm)**2/twosigma2 ) - &
                                2.0_DP*EXP(-rm**2/twosigma2) ) * grid%r2(ir)
            ELSE
               gaussian(ir) = 0.0_dp
            END IF
         END DO
         DO l3 = max (l1-l2,l2-l1), l1+l2
            ! 
            aux(1:irc) = gaussian(1:irc) * grid%r(1:irc)**l3
            ! calculate the corresponding multipole
            raux=int_0_inf_dr(aux,grid,irc,l3+2)
            IF (ABS(raux) < eps8) CALL errore &
               ('ld1_to_paw','norm of augm. func. is too small',ns*100+ns1)
            raux=augmom(ns,ns1,l3)/raux
            qfuncl(1:mesh,ns,ns1,l3) = raux*gaussian(1:mesh)
            qfuncl(1:mesh,ns1,ns,l3) = raux*gaussian(1:mesh)
            ! 
         END DO
         !
      CASE ('BESSEL')
         ALLOCATE (j1 (grid%mesh,2))
         ! ... or linear combination of Bessel functions?
         DO ir=1,irc
            IF (grid%r(ir)<rmatch_augfun) ircm=ir
         END DO
         DO l3 = max(l1-l2,l2-l1), l1+l2 
            ! 
            CALL find_bes_qi(qc,grid%r(ircm),l3,2,iok)
            IF (iok.ne.0) CALL errore('pseudo_q','problems with find_aug_qi',1)
            DO nc = 1, 2
               !
               CALL sph_bes(irc,grid%r(1),qc(nc),l3,j1(1,nc))
               aux(1:irc) = j1(1:irc,nc) * grid%r(1:irc)**(l3+2)
               b1(nc) = j1(ircm,nc)
               b2(nc) = int_0_inf_dr(aux,grid,ircm,l3+2)
               !
            ENDDO
            xc(1) = b1(2) / (b1(2) * b2(1) - b1(1) * b2(2))
            xc(2) = - b1(1) * xc(1) / b1(2)
            qfuncl(1:ircm,ns,ns1,l3) = &
                   augmom(ns,ns1,l3) * grid%r2(1:ircm) * &
                   (xc(1) * j1(1:ircm,1) + xc(2) * j1(1:ircm,2)) 
            qfuncl(1:mesh,ns1,ns,l3)=qfuncl(1:mesh,ns,ns1,l3)
            !
         END DO 
         DEALLOCATE (j1)
         !
      CASE DEFAULT
         !
         CALL errore ('pseudo_q',&
                 'Specified augmentation functions not allowed or coded',1)
         !
      END SELECT
   END DO
END DO
DEALLOCATE(augmom)
!
RETURN
END SUBROUTINE pseudo_q
