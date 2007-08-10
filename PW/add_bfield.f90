!
! Copyright (C) 2006 Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE add_bfield (v,rho)
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
  USE constants,        ONLY : pi
  USE io_global,        ONLY : stdout
  USE ions_base,        ONLY : nat, ntyp => nsp, ityp
  USE cell_base,        ONLY : omega
  USE gvect,            ONLY : nr1, nr2, nr3, nrxx 
  USE lsda_mod,         ONLY : nspin
  USE noncollin_module, ONLY : magtot_nc, bfield, lambda, i_cons, mcons, &
       pointlist, pointnum, factlist 
  IMPLICIT NONE
  !
  REAL(dp) :: v(nrxx, nspin), rho(nrxx,nspin)
  REAL(dp) :: ma, xx, fact, m1(3), m_loc(3,nat), r_loc(nat)

  INTEGER :: ir, ipol, nt, na
  REAL(DP) :: etcon


  IF (i_cons==0) RETURN
  !
  ! get the actual values of the local integrated quantities
  etcon=0.d0
  IF (i_cons.LT.3) THEN
     CALL get_locals(r_loc, m_loc, rho)

     DO na = 1,nat
        nt = ityp(na)
        IF (i_cons==1) THEN
           ! i_cons = 1 means that the 3 components of the magnetization
           ! are constrained, they are given in the input-file
           DO ipol = 1,3
              m1(ipol) = m_loc(ipol,na) - mcons(ipol,nt)
           END DO
           etcon = etcon + &
                lambda * (m1(1)*m1(1) + m1(2)*m1(2) + m1(3)*m1(3) )
        ELSE IF (i_cons==2) THEN
           ! i_cons = 2 means that the angle theta between the local
           ! magn. moment and the z-axis is constrained
           ! mcons (3,nt) is the cos of the constraining angle theta
           ! the penalty functional in this case is
           ! lambda*(m_loc(z)/|m_loc| - cos(theta) )^2
           ma = dsqrt(m_loc(1,na)**2+m_loc(2,na)**2+m_loc(3,na)**2)
           if (ma.lt.1.d-30) call errore('add_bfield', &
                'local magnetization is zero',1)
           xx=(m_loc(3,na)/ma - mcons(3,nt))
           m1(1) = - xx*m_loc(1,na)*m_loc(3,na) / (ma*ma*ma)
           m1(2) = - xx*m_loc(2,na)*m_loc(3,na) / (ma*ma*ma)
           m1(3) =   xx*(-m_loc(3,na)*m_loc(3,na) / (ma*ma*ma) + 1.d0/ma)
           etcon = etcon + &
                lambda * (m_loc(3,na)/ma - mcons(3,nt))**2

        END IF

        DO ir = 1,pointnum(na)
           fact = 2.D0*lambda*factlist(ir,na)*omega/(nr1*nr2*nr3)
           DO ipol = 1,3
              v(pointlist(ir,na),ipol+1) = v(pointlist(ir,na),ipol+1) &
                   + fact*m1(ipol)
           END DO       ! ipol
        END DO      ! points
     END DO      ! na
     write (stdout,'(4x,a,F15.8)' ) " constraint energy (Ryd) = ", etcon
  ELSE IF (i_cons==3.or.i_cons==6) THEN
     m1 = 0.d0
     DO ir = 1,nrxx
        DO ipol = 1, 3
           m1(ipol) = m1(ipol) + rho(ir,ipol+1)
        END DO
     END DO
     CALL reduce( 3, m1 )
     DO ipol = 1, 3
        m1(ipol) = m1(ipol) * omega / ( nr1 * nr2 * nr3 )
     END DO

     IF (i_cons==3) THEN
       fact = 2.D0*lambda
       DO ipol=1,3
          bfield(ipol)=-fact*(m1(ipol)-mcons(ipol,1))
       END DO
       write(6,'(5x," External magnetic field: ", 3f13.5)') &
            (bfield(ipol),ipol=1,3)
       DO ipol = 2,4
          fact=bfield(ipol-1)
          DO ir =1,nrxx
             v(ir,ipol) = v(ir,ipol)-fact
          END DO
       END DO
     END IF

     IF (i_cons==6) THEN
       fact = TAN(mcons(3,1)/180.d0*pi)
       write(6,'(5x,"-angle- ", f13.5)') atan(m1(1)/m1(3))*180.d0/pi
       write(6,'(5x," E_constrain: ", f13.9)') (m1(1)-m1(3)*fact)**2*lambda
       write(6,'(5x," lambda: ", f10.3)') lambda
       bfield(1) = 2.d0*lambda*(m1(1)-m1(3)*fact)
       bfield(2) = 0.d0
       bfield(3) = -2.d0*lambda*fact*(m1(1)-m1(3)*fact)
       write(6,'(5x," External magnetic field fixed: ", 3f13.5)') &
          (bfield(ipol),ipol=1,3)
       DO ipol = 2,4
          fact=bfield(ipol-1)
          DO ir =1,nrxx
             v(ir,ipol) = v(ir,ipol)+fact
          END DO
       END DO

     END IF

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
