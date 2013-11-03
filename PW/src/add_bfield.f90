!
! Copyright (C) 2006 Quantum ESPRESSO group
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
  ! magnetic field (apparently disabled at the moment).
  ! 
  ! NB: So far, the contribution of the orbital currents 
  !     to the magnetization is not included.
  !
  !
  USE kinds,            ONLY : DP
  USE constants,        ONLY : pi
  USE io_global,        ONLY : stdout
  USE ions_base,        ONLY : nat, ntyp => nsp, ityp
  USE cell_base,        ONLY : omega
  USE fft_base,         ONLY : dfftp
  USE lsda_mod,         ONLY : nspin
  USE mp_bands,         ONLY : intra_bgrp_comm
  USE mp,               ONLY : mp_sum
  USE noncollin_module, ONLY : bfield, lambda, i_cons, mcons, &
                               pointlist, factlist, noncolin
  IMPLICIT NONE
  ! input/outpt variables
  REAL(DP), INTENT(IN) :: rho(dfftp%nnr,nspin)
  REAL(DP), INTENT(INOUT) :: v(dfftp%nnr, nspin)
  ! local variables
  REAL(DP) :: ma, mperp, xx, fact, m1(3), etcon, fact1(3)
  REAL(DP), allocatable :: m2(:,:), m_loc(:,:), r_loc(:)

  INTEGER :: ir, ipol, nt, na, npol

  etcon=0.D0

  IF (nspin ==1 .or. i_cons==0)  RETURN
  ! i_cons==0, no constraint
  npol = nspin - 1  ! number of relevant magnetic components
                    ! 3 for non-collinear case; 1 for collinear case
  !
  ! get the actual values of the local integrated quantities
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
           IF (.NOT. noncolin) CALL errore('add_bfield', &
              'this magnetic constraint only applies to non collinear calculations',2)
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
        DO ir = 1, dfftp%nnr
           if (pointlist(ir) .eq. 0 ) cycle
           fact = 2.D0*lambda*factlist(ir)*omega/(dfftp%nr1*dfftp%nr2*dfftp%nr3)
           DO ipol = 1,3
              v(ir,ipol+1) = v(ir,ipol+1) + fact*m2(ipol,pointlist(ir))
           END DO       ! ipol
        END DO      ! points
     else
        DO ir = 1, dfftp%nnr
           if (pointlist(ir) .eq. 0 ) cycle
           fact = 2.D0*lambda*factlist(ir)*omega/(dfftp%nr1*dfftp%nr2*dfftp%nr3)
           v(ir,1) = v(ir,1) + fact*m2(1,pointlist(ir))
           v(ir,2) = v(ir,2) - fact*m2(1,pointlist(ir))
        END DO      ! points
     end if
     deallocate (m2, m_loc, r_loc)

     write (stdout,'(4x,a,F15.8)' ) " constraint energy (Ryd) = ", etcon
  ELSE IF (i_cons==3.or.i_cons==6) THEN
     m1 = 0.d0
     IF (npol==1) THEN
        DO ir = 1,dfftp%nnr
           m1(1) = m1(1) + rho(ir,1) - rho(ir,2)
        END DO
        m1(1) = m1(1) * omega / ( dfftp%nr1 * dfftp%nr2 * dfftp%nr3 )
     ELSE
        DO ipol = 1, 3
           DO ir = 1,dfftp%nnr
              m1(ipol) = m1(ipol) + rho(ir,ipol+1)
           END DO
           m1(ipol) = m1(ipol) * omega / ( dfftp%nr1 * dfftp%nr2 * dfftp%nr3 )
        END DO
     END IF
     CALL mp_sum( m1, intra_bgrp_comm )

     IF (i_cons==3) THEN
       IF (npol==1) THEN
          fact = 2.D0*lambda
          bfield(1)=-fact*(m1(1)-mcons(1,1))
          DO ir =1,dfftp%nnr
             v(ir,1) = v(ir,1)-bfield(1)
             v(ir,2) = v(ir,2)+bfield(1)
          END DO
       ELSE
          fact = 2.D0*lambda
          DO ipol=1,3
             bfield(ipol)=-fact*(m1(ipol)-mcons(ipol,1))
             DO ir =1,dfftp%nnr
                v(ir,ipol+1) = v(ir,ipol+1)-bfield(ipol)
             END DO
          END DO
       END IF
       write(stdout,'(5x," External magnetic field: ", 3f13.5)') &
            (bfield(ipol),ipol=1,npol)
     END IF

     IF (i_cons==6) THEN
       !
       IF (.NOT. noncolin)  CALL errore('add_bfield', &
         'this magnetic constraint only applies to non collinear calculations',6)
       !
       ! penalty functional: E = lambda*(arccos(m_z/|m|) - theta)^2
       !
       ! modulus and azimuthal component of the magnetization:
       ma = SQRT(m1(1)**2 + m1(2)**2 + m1(3)**2)
       mperp = SQRT(m1(1)**2 + m1(2)**2)
       IF (ma < 1.D-12)  CALL errore('add_bfield', &
           'magnetization too small, cannot constrain polar angle', 1)
       fact = ACOS(m1(3)/ma)
       xx = fact - mcons(3,1)/180.D0*pi
       IF (mperp < 1.D-14) THEN
          fact1(1:2) = 0.D0
          ! when m is along z, in order to allow the magnetization to rotate
          ! add a tiny B_ext along x (when required, because of theta-target > 0)
          IF (mcons(3,1) > 0.D0) fact1(1) = 1.D-14
       ELSE
          fact1(1:2) = m1(1:2)/mperp * m1(3)/ma/ma
       ENDIF
       fact1(3) = - SQRT(1.D0 - (m1(3)/ma)**2)/ma
       etcon = lambda * xx**2
       bfield(:) = 2.D0 * lambda * xx * fact1(:)
       DO ipol = 1,3
          DO ir =1,dfftp%nnr
             v(ir,ipol+1) = v(ir,ipol+1)+bfield(ipol)
          END DO
       END DO
       !
       write(stdout,'(/,5x,"Constraint on the polar angle of the magnetization")')
       ! N.B.: since the magnetization is here computed starting from the mixed
       ! rho (i.e. the input rho for the next scf iteration), as all the other
       ! contributions to the potential for the next iteration, it will differ
       ! from the magnetization written on the output, since that is calculated
       ! with the output rho of the current iteration. At convergence the two
       ! magnetizations will coincide (and so will do the polar angles).
       write(stdout,'(5x,"theta (target): ",F10.5,"     (",F10.5,")")') &
             ACOS(m1(3)/ma)*180.d0/pi, mcons(3,1)
       write(stdout,'(5x,"E_constraint:  ",F15.9," (lambda:",F15.9,")")') etcon, lambda
       write(stdout,'(5x,"External magnetic field: ", 3F12.6)') bfield(1:npol)
       !write(stdout,'(5x,"Magnetization          : ", 3F12.6)') m1(1:npol)
       !
     END IF
  ELSE IF (i_cons==4) THEN
     write(stdout,'(5x," External magnetic field: ", 3f13.5)') &
             (bfield(ipol),ipol=1,npol)
     IF (npol==1) THEN
        DO ir =1,dfftp%nnr
           v(ir,1) = v(ir,1)-bfield(ipol)
           v(ir,2) = v(ir,2)+bfield(ipol)
        END DO
     ELSE
        DO ipol = 1,3
           DO ir =1,dfftp%nnr
              v(ir,ipol+1) = v(ir,ipol+1)-bfield(ipol)
           END DO
        END DO
     END IF
  ELSE
     CALL errore('add_bfield','i_cons not programmed',1)
  END IF

  RETURN
END SUBROUTINE add_bfield
