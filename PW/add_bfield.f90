SUBROUTINE add_bfield (v)
!--------------------------------------------------------------------
!
! If noncolinear is set, one can calculate constrains either on 
! the local magnetization, calculated in get_locals or on the 
! total magnetization.
! 
! To this end, a "penalty term" of the form
! E_p = lambda * ( m_loc - m_loc_constr)^2
! is added to the energy. Here we calculate the resulting
! "constraining B-field" and add it to v(ir,2..4)
! Moreover there is also the possibility to add a fixed
! magnetic field. 
! 
! NB: So far, the contribution of the orbital currents 
!     to the magnetization is not included.
!
!
USE kinds,            ONLY : dp
USE ions_base,        ONLY : nat, ntyp => nsp, ityp
USE cell_base,        ONLY : omega
USE gvect,            ONLY : nr1, nr2, nr3, nrxx 
USE lsda_mod,         ONLY : nspin
USE noncollin_module, ONLY : magtot_nc, bfield, lambda, i_cons, mcons, &
                             pointlist, pointnum, factlist, m_loc, r_loc 
IMPLICIT NONE
  !
REAL(KIND=dp) :: v(nrxx, nspin)
REAL(KIND=dp) :: ma, xx, xxx, fact, m1(3)

INTEGER :: ir, ipol, nt, na
! counters

IF (i_cons==0) RETURN
!
! get the actual values of the local integrated quantities
IF (i_cons.LT.3) THEN
   CALL get_locals(r_loc, m_loc)

   DO na = 1,ntyp
      nt = ityp(na)
      IF (i_cons==1) THEN
         ! i_cons = 1 means that the 3 components of the magnetization
         ! are constrained, they are given in the input-file
         DO ipol = 1,3
            m1(ipol) = m_loc(ipol,nt) - mcons(ipol,nt)
         END DO
      ELSE IF (i_cons==2) THEN
         ! i_cons = 2 means that the angle theta between the local
         ! magn. moment and the z-axis is constrained
         ! mcons (1,nt) is the cos of the constraining angle theta
         ! the penalty functional in this case is
         ! lambda*(acos(m_loc(z)/|m_loc|) - theta )^2
         ma = dsqrt(m_loc(1,nt)**2+m_loc(2,nt)**2+m_loc(3,nt)**2)
         xx = m_loc(3,nt)/ma
         IF (ABS(xx - 1.D0).GT.1.D-10) THEN
            xxx = - (ACOS(xx) - mcons(1,nt))/SQRT(1.D0 - xx*xx)
         ELSE
            xxx = -1.D0
         END IF
         m1(1) = - xxx * m_loc(1,nt)*m_loc(3,nt) / (ma*ma*ma)
         m1(2) = - xxx * m_loc(2,nt)*m_loc(3,nt) / (ma*ma*ma)
         m1(3) =   xxx * (ma*ma-m_loc(3,nt)*m_loc(3,nt)) / (ma*ma*ma)
      END IF

      DO ir = 1,pointnum(na)
         fact = 2.D0*lambda*factlist(ir,na)*omega/(nr1*nr2*nr3)
         DO ipol = 1,3
            v(pointlist(ir,na),ipol+1) = v(pointlist(ir,na),ipol+1) &
                                       + fact*m1(ipol)
         END DO       ! ipol
      END DO      ! points
   END DO      ! na
ELSE IF (i_cons==3) THEN
   fact = 2.D0*lambda
   DO ipol=1,3
      bfield(ipol)=-fact*(magtot_nc(ipol)-mcons(ipol,1))
   END DO
   write(6,'(5x," External magnetic field: ", 3f13.5)') &
                                          (bfield(ipol),ipol=1,3)
   DO ipol = 2,4
      fact=bfield(ipol-1)
      DO ir =1,nrxx
         v(ir,ipol) = v(ir,ipol)-fact
      END DO            
   END DO              
ELSE IF (i_cons==4) THEN
   bfield(:)=mcons(:,1)
   write(6,'(5x," External magnetic field fixed: ", 3f13.5)') &
                                         (bfield(ipol),ipol=1,3)
   DO ipol = 2,4
      fact=bfield(ipol-1)
      DO ir =1,nrxx
         v(ir,ipol) = v(ir,ipol)-fact
      END DO            
   END DO
ELSE
   CALL errore('add_bfield','i_cons not programmed',1)
END IF

RETURN
END SUBROUTINE add_bfield
