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
  USE mp_global,        ONLY : intra_pool_comm
  USE mp,               ONLY : mp_sum
  USE noncollin_module, ONLY : magtot_nc, bfield, lambda, i_cons, mcons, &
       pointlist, factlist, noncolin
  IMPLICIT NONE
  !
  REAL(dp) :: v(nrxx, nspin), rho(nrxx,nspin)
  REAL(dp) :: ma, xx, fact, m1(3)
  REAL(dp), allocatable :: m2(:,:), m_loc(:,:), r_loc(:)

  INTEGER :: ir, ipol, nt, na, npol
  REAL(DP) :: etcon

  IF (nspin ==1 .or. i_cons==0 ) RETURN

  npol = nspin - 1  ! number of relevant magnetic components
                    ! 3 for non-collinear case; 1 for collinear case
  !
  ! get the actual values of the local integrated quantities
  etcon=0.d0
  IF (i_cons.LT.3) THEN
     allocate ( m2(npol,nat), m_loc(npol,nat), r_loc(nat) )

     CALL get_locals(r_loc, m_loc, rho)

     DO na = 1,nat
        nt = ityp(na)
        IF (i_cons==1) THEN
           ! i_cons = 1 means that the npol components of the magnetization
           ! are constrained, they are given in the input-file
           m2(1:npol,na) = m_loc(1:npol,na) - mcons(1:npol,nt)
           do ipol=1,npol
              etcon = etcon + lambda * m2(ipol,na)*m2(ipol,na) 
           end do
        ELSE IF (i_cons==2) THEN
           ! i_cons = 2 means that the angle theta between the local
           ! magn. moment and the z-axis is constrained
           ! mcons (3,nt) is the cos of the constraining angle theta
           ! the penalty functional in this case is
           ! lambda*(m_loc(z)/|m_loc| - cos(theta) )^2
           if (.not. noncolin) &
              call errore('add_bfield','this magnetic constraint only applies to non collinear calculations',2)
           ma = dsqrt(m_loc(1,na)**2+m_loc(2,na)**2+m_loc(3,na)**2)
           if (ma.lt.1.d-30) call errore('add_bfield', &
                'local magnetization is zero',1)
           xx=(m_loc(3,na)/ma - mcons(3,nt))
           m2(1,na) = - xx*m_loc(1,na)*m_loc(3,na) / (ma*ma*ma)
           m2(2,na) = - xx*m_loc(2,na)*m_loc(3,na) / (ma*ma*ma)
           m2(3,na) =   xx*(-m_loc(3,na)*m_loc(3,na) / (ma*ma*ma) + 1.d0/ma)
           etcon = etcon + &
                lambda * (m_loc(3,na)/ma - mcons(3,nt))**2

        END IF
     END DO ! na

     if (noncolin) then
        DO ir = 1, nrxx
           if (pointlist(ir) .eq. 0 ) cycle
           fact = 2.D0*lambda*factlist(ir)*omega/(nr1*nr2*nr3)
           DO ipol = 1,3
              v(ir,ipol+1) = v(ir,ipol+1) + fact*m2(ipol,pointlist(ir))
           END DO       ! ipol
        END DO      ! points
     else
        DO ir = 1, nrxx
           if (pointlist(ir) .eq. 0 ) cycle
           fact = 2.D0*lambda*factlist(ir)*omega/(nr1*nr2*nr3)
           v(ir,1) = v(ir,1) + fact*m2(1,pointlist(ir))
           v(ir,2) = v(ir,2) - fact*m2(1,pointlist(ir))
        END DO      ! points
     end if
     deallocate (m2, m_loc, r_loc)

     write (stdout,'(4x,a,F15.8)' ) " constraint energy (Ryd) = ", etcon
  ELSE IF (i_cons==3.or.i_cons==6) THEN
     m1 = 0.d0
     IF (npol==1) THEN
        DO ir = 1,nrxx
           m1(1) = m1(1) + rho(ir,1) - rho(ir,2)
        END DO
        m1(1) = m1(1) * omega / ( nr1 * nr2 * nr3 )
     ELSE
        DO ipol = 1, 3
           DO ir = 1,nrxx
              m1(ipol) = m1(ipol) + rho(ir,ipol+1)
           END DO
           m1(ipol) = m1(ipol) * omega / ( nr1 * nr2 * nr3 )
        END DO
     END IF
     CALL mp_sum( m1, intra_pool_comm )

     IF (i_cons==3) THEN
       IF (npol==1) THEN
          fact = 2.D0*lambda
          bfield(1)=-fact*(m1(1)-mcons(1,1))
          DO ir =1,nrxx
             v(ir,1) = v(ir,1)-bfield(1)
             v(ir,2) = v(ir,2)+bfield(1)
          END DO
       ELSE
          fact = 2.D0*lambda
          DO ipol=1,3
             bfield(ipol)=-fact*(m1(ipol)-mcons(ipol,1))
             DO ir =1,nrxx
                v(ir,ipol+1) = v(ir,ipol+1)-bfield(ipol)
             END DO
          END DO
       END IF
       write(6,'(5x," External magnetic field: ", 3f13.5)') &
            (bfield(ipol),ipol=1,npol)
     END IF

     IF (i_cons==6) THEN
        if (.not. noncolin) &
           call errore('add_bfield','this magnetic constraint only applies to non collinear calculations',6)
       fact = TAN(mcons(3,1)/180.d0*pi)
       write(6,'(5x,"-angle- ", f13.5)') atan(m1(1)/m1(3))*180.d0/pi
       write(6,'(5x," E_constrain: ", f13.9)') (m1(1)-m1(3)*fact)**2*lambda
       write(6,'(5x," lambda: ", f10.3)') lambda
       bfield(1) = 2.d0*lambda*(m1(1)-m1(3)*fact)
       bfield(2) = 0.d0
       bfield(3) =-2.d0*lambda*fact*(m1(1)-m1(3)*fact)
       write(6,'(5x," External magnetic field fixed: ", 3f13.5)') &
          (bfield(ipol),ipol=1,3)
       DO ipol = 1,3
          DO ir =1,nrxx
             v(ir,ipol+1) = v(ir,ipol+1)+bfield(ipol)
          END DO
       END DO
     END IF
  ELSE IF (i_cons==4) THEN
     bfield(1:npol)=mcons(1:npol,1)
     write(6,'(5x," External magnetic field fixed: ", 3f13.5)') &
             (bfield(ipol),ipol=1,npol)
     IF (npol==1) THEN
        DO ir =1,nrxx
           v(ir,1) = v(ir,1)-bfield(ipol)
           v(ir,2) = v(ir,2)+bfield(ipol)
        END DO
     ELSE
        DO ipol = 1,3
           DO ir =1,nrxx
              v(ir,ipol+1) = v(ir,ipol+1)-bfield(ipol)
           END DO
        END DO
     END IF
  ELSE
     CALL errore('add_bfield','i_cons not programmed',1)
  END IF

  RETURN
END SUBROUTINE add_bfield
