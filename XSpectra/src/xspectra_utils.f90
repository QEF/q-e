!
! Copyright (C) 2009-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
SUBROUTINE stop_xspectra
  !----------------------------------------------------------------------------
  !
  ! Synchronize processes before stopping. This is a copy of stop_pp.
  !
  USE io_files, ONLY: iunwfc
  USE mp_global, ONLY: mp_global_end
  !
#if defined(__MPI)

  INTEGER :: info
  LOGICAL :: op

  INQUIRE ( iunwfc, opened = op )

  IF ( op ) CLOSE (unit = iunwfc, status = 'delete')

  CALL mp_global_end()

#endif

  STOP
END SUBROUTINE stop_xspectra

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
FUNCTION lastterm(a,b,g)
  !----------------------------------------------------------------------------
  USE kinds, ONLY: dp
  IMPLICIT NONE
  REAL(dp)    :: a, b, g, y1, y2, z1, z2, r
  COMPLEX(dp) :: lastterm

  y1 = a*a - g*g - 4*b
  y2 = -2*a*g
  r  = 0.5*SQRT(y1*y1 + y2*y2)

  IF (g<0) THEN
     z1 =  a/2 + 0.5*SIGN(SQRT(y1/2 + r),y2)
     z2 = -g/2 + 0.5*SQRT(-y1/2 + r)
  ELSE
     z1 =  a/2 - 0.5*SIGN(SQRT(y1/2 + r),y2)
     z2 = -g/2 - 0.5*SQRT(-y1/2+r)
  ENDIF
 
 lastterm = CMPLX(z1,z2,kind=DP)

END FUNCTION lastterm

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
FUNCTION paste_fermi(e,e0,a,b,gamma,m,term,first)
  !----------------------------------------------------------------------------
  USE kinds,     ONLY: dp
  USE constants, ONLY: tpi
  USE xspectra,  ONLY: xnitermax, xcheck_conv
  USE cut_valence_green, ONLY: cut_ierror, cut_stepu, cut_stepl, &
                               cut_startt, cut_tinf, cut_tsup,   &
                               cut_nmemu, cut_nmeml, memu, meml
  IMPLICIT NONE

  REAL(dp) :: paste_fermi
  ! Arguments
  INTEGER  :: m
  REAL(dp) :: a(xnitermax)
  REAL(dp) :: b(xnitermax)
  REAL(dp) :: gamma
  REAL(dp) :: e, e0
  LOGICAL  :: term, first
  ! Local
  COMPLEX(dp) :: green, y, dy, c1, c2, e1, e2
  REAL(dp)    :: t, dt, t1, ta, tb
  INTEGER, save :: n1
  INTEGER, save :: n2
  INTEGER :: nn1, nn2

  IF (first) THEN
     memu(:,:) = (0.d0,0.d0)
     meml(:,:) = (0.d0,0.d0) 
     n1 = 0
     n2 = 0
     first = .false.
  ENDIF

  dy = cut_ierror + 1.0
  y = 0.d0

  nn1 = 1
  nn2 = 1

  t1 = 0.5773502692

  t = cut_startt

  DO WHILE ((abs(dy)>cut_ierror).OR.(t<cut_tsup))
     dt = cut_stepu*t
     ta = t + dt*(1-t1)/2
     tb = t + dt*(1+t1)/2
     e1 = CMPLX(e0,ta,kind=DP)
     e2 = CMPLX(e0,tb,kind=DP)

     IF (nn1>n1) THEN
        c1 = green(a,b,e1,m,term)
        c2 = green(a,b,e2,m,term)
        IF (nn1<cut_nmemu) THEN
           memu(nn1,1) = c1
           memu(nn1,2) = c2
           n1 = nn1
        ENDIF
     ELSE
        c1 = memu(nn1,1)
        c2 = memu(nn1,2)
     ENDIF

     dy = (dt/2) * &
          ( c1/CMPLX(e0-e,ta-gamma,kind=DP)          &
           + CONJG(c1)/CMPLX(e0-e,-ta-gamma,kind=DP) &
           + c2/CMPLX(e0-e,tb-gamma,kind=DP)         &
           + CONJG(c2)/CMPLX(e0-e,-tb-gamma,kind=DP) )
     y = y + dy
     t = t + dt
     nn1 = nn1 + 1
  ENDDO

  t = cut_startt
  dy = cut_ierror + 1

  DO WHILE((abs(dy)>cut_ierror).OR.(t>cut_tinf))
     dt = cut_stepl * t
     ta = t - dt*(1-t1)/2
     tb = t - dt*(1+t1)/2
     e1 = CMPLX(e0,ta,kind=DP)
     e2 = CMPLX(e0,tb,kind=DP)

     IF (nn2>n2) THEN
        c1 = green(a,b,e1,m,term)
        c2 = green(a,b,e2,m,term)
        IF (nn2<cut_nmeml) THEN
           meml(nn2,1) = c1
           meml(nn2,2) = c2
           n2 = nn2
        ENDIF
     ELSE
        c1 = meml(nn2,1)
        c2 = meml(nn2,2)
     ENDIF

     dy = (dt/2) * &
          ( c1/CMPLX(e0-e,ta-gamma,kind=DP)           &
            + CONJG(c1)/CMPLX(e0-e,-ta-gamma,kind=DP) &
            + c2/CMPLX(e0-e,tb-gamma,kind=DP)         &
            + CONJG(c2)/CMPLX(e0-e,-tb-gamma,kind=DP))
     y = y + dy
     t = t - dt
     nn2 = nn2 + 1
  ENDDO

  paste_fermi = AIMAG(y)/tpi

END FUNCTION paste_fermi

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
FUNCTION green(a,b,e,m,term)
  !----------------------------------------------------------------------------
  USE kinds,    ONLY : dp
  USE xspectra, ONLY : xnitermax, xcheck_conv

  IMPLICIT NONE

  COMPLEX(dp) :: green
  ! Arguments
  INTEGER  :: m
  REAL(dp) :: a(xnitermax)
  REAL(dp) :: b(xnitermax)
  COMPLEX(dp) :: e
  LOGICAL :: term
  ! Local
  INTEGER :: i, p, q
  COMPLEX(dp) :: res, lastterm 
  REAL(dp) :: aa, bb

  q = xcheck_conv/2
  IF (term) THEN
     aa = 0.0
     bb = 0.0
     DO p = 1, q
        aa = aa + a(m-p)
        bb = bb + b(m-p)
     ENDDO
     aa = aa/q
     bb = bb/q

     res = lastterm(aa-REAL(e), bb*bb, AIMAG(e))
  ELSE
     res = CMPLX(a(m)-REAL(e),AIMAG(e),kind=DP)
  ENDIF
  DO i = 1, m-1
     res = a(m-i) - e - b(m-i)*b(m-i)/res
  ENDDO

  green = 1/res

END FUNCTION green

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
SUBROUTINE write_status_of_the_code
  !----------------------------------------------------------------------------
  USE io_global,       ONLY : stdout
  IMPLICIT NONE
  WRITE(stdout,'(5x,a)')&
   '-------------------------------------------------------------------------'
  WRITE (stdout,'(5x,a)') &
   '                      STATUS OF THE CODE (22/04/2009) '
  WRITE(stdout,'(5x,a)')&
   '-------------------------------------------------------------------------'
  WRITE (stdout,'(5x,a)') 'Working features (22/04/2009)'
  WRITE (stdout,'(5x,a)') '-----------------------------'
  WRITE (stdout,'(5x,a)') &
   '- XANES works both in the electric-dipole and -quadrupole approximation,'
  WRITE (stdout,'(5x,a)') '- Spin polarized works'
  WRITE (stdout,'(5x,a)') '- DFT+U implemented, validated'
  WRITE (stdout,'(5x,a)') '- Ultrasoft pseudo works'
  WRITE (stdout,'(5x,a)') '- Cut occupied states working, improved'
  WRITE (stdout,'(5x,a)') '- Terminator working'
  WRITE (stdout,'(5x,a)') '- Multiprojectors TM+USPP working (MCB,CG)'
  WRITE (stdout,'(5x,a)') '- New save file format, with version numbering'
  WRITE (stdout,'(5x,a)') &
   '- Time limit implemented, with restart, seems to work'
  WRITE (stdout,'(5x,a)') &
   '- DFT+U tested ONLY for non ortho wfc, but implemented'
  WRITE(stdout,*) 
  WRITE (stdout,'(5x,a)') 'TO DO'
  WRITE (stdout,'(5x,a)') '-----'
  WRITE (stdout,'(5x,a)') '- L2,3 edges [OB]'
  WRITE (stdout,'(5x,a)') '- Generalization to all edges [OB]'
  WRITE (stdout,'(5x,a)') '- XMCD [?]'
  WRITE (stdout,'(5x,a)') '- IXS [DC]'
  WRITE (stdout,'(5x,a)') '- EELS [DC]'
  WRITE (stdout,'(5x,a)') '- REXS [DC]'
  WRITE (stdout,'(5x,a)') '- Bethe-Salpeter [?] '
  WRITE (stdout,'(5x,a)') '- RXES [?]' 

END SUBROUTINE write_status_of_the_code

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
SUBROUTINE determine_polycut(t1,t2,f1,f2,df1,df2,poly)
  !----------------------------------------------------------------------------
  ! Calculates the interpolation polynome between 2 points - CG
  !----------------------------------------------------------------------------
  USE kinds, ONLY: dp

  IMPLICIT NONE
  REAL(dp) :: t1, t2, f1, f2, df1, df2, poly(4)

  poly(4) = ( (t2-t1) * (df2+df1) - 2*(f2-f1) )/( (t2-t1)**3 )
  poly(3) = (df2-df1)/(2*(t2-t1)) - &
            1.5d0*(t2+t1)*( (t2-t1)*(df2+df1) - 2*(f2-f1) )/( (t2-t1)**3 )
  poly(2) = df1 - 2*t1*poly(3) - 3*t1*t1*poly(4)
  poly(1) = f1 - poly(2)*t1 - poly(3)*t1**2 - poly(4)*t1**3

END SUBROUTINE determine_polycut
