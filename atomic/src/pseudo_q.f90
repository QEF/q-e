! Copyright (C) 2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
SUBROUTINE pseudo_q (qfunc, qfuncl)       
USE kinds, ONLY : DP
USE io_global, ONLY : stdout
USE ld1_parameters, ONLY : nwfsx
USE ld1inc, ONLY : rcut, lls,  grid, ndmx, lmx2, nbeta, ikk, ecutrho, &
            rmatch_augfun, rmatch_augfun_nc
IMPLICIT NONE
!
REAL(DP), INTENT(IN) :: qfunc(ndmx,nwfsx,nwfsx)
REAL(DP), INTENT(OUT) :: qfuncl(ndmx,nwfsx,nwfsx,0:lmx2) 
REAL(DP),  EXTERNAL :: int_0_inf_dr
!
! variables for aug. functions generation
! 
INTEGER  :: irc, ns, ns1, l1, l2, l3, lll, mesh, n, ik
INTEGER  :: l1_e, l2_e
REAL(DP) :: aux(ndmx)
REAL(DP) :: augmom, ecutrhoq, rmatch

ecutrho=0.0_DP
mesh = grid%mesh
qfuncl=0.0_DP
do ns=1,nbeta
   l1 = lls(ns)
   do ns1=ns,nbeta
      l2 = lls(ns1)
      !
      !  Find the matching point
      !
      ik=0
      IF (rmatch_augfun_nc) THEN
         rmatch=min(rcut(ns),rcut(ns1))
      ELSE
         rmatch=rmatch_augfun
      ENDIF
      do n=1,mesh
         if (grid%r(n)>rmatch) then
            ik=n
            exit
         endif
      enddo
      if (ik==0.or.ik>mesh-20) call errore('pseudo_q','wrong rmatch_augfun',1)
      !
      ! Do the pseudization
      !
      do l3 = abs(l1-l2), l1+l2, 2
         CALL compute_q_3bess(l3,l1+l2,ik,qfunc(1,ns,ns1), &
                                          qfuncl(1,ns,ns1,l3),ecutrhoq)
         IF (ecutrhoq>ecutrho) then
            ecutrho=ecutrhoq
            l1_e=l1
            l2_e=l2
         ENDIF
         qfuncl(1:mesh,ns1,ns,l3)=qfuncl(1:mesh,ns,ns1,l3)
      end do
   end do
end do
!
!  Check that multipoles have not changed
!
irc = maxval(ikk(1:nbeta))+8
augmom=0.0_DP
DO ns=1,nbeta
   l1=lls(ns)
   DO ns1=ns,nbeta
      l2=lls(ns1)
      DO l3 = abs(l1-l2), l1+l2, 2
         aux(1:irc) = (qfuncl(1:irc,ns,ns1,l3)-qfunc(1:irc,ns,ns1)) &
                                               * grid%r(1:irc)**l3
         lll = l1 + l2 + 2 + l3
         augmom=int_0_inf_dr(aux(1:irc),grid,irc,lll)
      
         IF (abs(augmom)>1.d-5) WRITE (stdout,'(5x,a,2i3,a,2i3,a,i3,f15.7)') &
              " Problem with multipole",ns,l1,":",ns1,l2, " l3=",l3, augmom
      END DO
   END DO
END DO
WRITE(stdout,'(/,5x, "Q pseudized with Bessel functions")')
WRITE(stdout,'(5x,"Expected ecutrho= ",f12.4," due to l1=",i3,"   l2=",i3)') &
                            ecutrho, l1_e, l2_e
RETURN
END SUBROUTINE pseudo_q
