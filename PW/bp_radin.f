C
C ---------------------------------------------------------------
      SUBROUTINE RADIN(MESH,C,FUNC,ASUM)
C ---------------------------------------------------------------
C     SIMPSONS RULE INTEGRATION FOR HERMAN SKILLMAN MESH
C     MESH - # OF MESH POINTS
C     C    - 0.8853418/Z**(1/3.)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FUNC(mesh)
      A1=0.0
      A2E=0.0
      ASUM=0.0
      H=0.0025*C
      NBLOCK=MESH/40
      I=1
c      FUNC(1)=0.0
      DO 39 J=1,NBLOCK
      DO 38 K=1,20
      I=I+2
      I1=I-1
      A2ES=A2E
      A2O=FUNC(I1)/12.0
      A2E=FUNC(I)/12.0
      A1=A1+5.0*A2ES+8.0*A2O-A2E
c      FUNC(I1)=ASUM+A1*H
      A1=A1-A2ES+8.0*A2O+5.0*A2E
c      FUNC(I)=ASUM+A1*H
      fi = ASUM+A1*H
   38 CONTINUE
c      ASUM=FUNC(I)
      asum = fi
      A1=0.0
   39 H=H+H
C
      RETURN
      END
C
c-----------------------------------------------------------------------
      subroutine radlg(mesh,func,r,dx,asum)
c-----------------------------------------------------------------------
c
c     simpson's rule integrator for function stored on the
c     radial logarithmic mesh
c
c.....logarithmic radial mesh information
      IMPLICIT REAL*8 (A-H,O-Z)
      dimension r(mesh)
c.....function to be integrated
      dimension func(mesh)
c
c.....variable for file = 0
c
c     routine assumes that mesh is an odd number so run check
      if ( mesh - ( mesh / 2 ) * 2 .ne. 1 ) then
        write(*,*) '***error in subroutine radlg'
        write(*,*) 'routine assumes mesh is odd but mesh =',mesh
        stop
      endif

      asum = func(1)*r(1)+func(mesh)*r(mesh)
      do  i = 2,mesh-1,2
         asum = asum + 4.0d0*func(i)*r(i)+2.0d0*func(i+1)*r(i+1)
      enddo
      asum = asum*dx/3.0d0
      return
      end

C
c-----------------------------------------------------------------------
      subroutine radlg1(mesh,func,rab,asum)
c-----------------------------------------------------------------------
c
c     simpson's rule integrator for function stored on the
c     radial logarithmic mesh
c
c.....logarithmic radial mesh information
      IMPLICIT REAL*8 (A-H,O-Z)
      dimension rab(mesh)
c.....function to be integrated
      dimension func(mesh)
c
c.....variable for file = 0
c
c     routine assumes that mesh is an odd number so run check
      if ( mesh - ( mesh / 2 ) * 2 .ne. 1 ) then
        write(*,*) '***error in subroutine radlg'
        write(*,*) 'routine assumes mesh is odd but mesh =',mesh
        stop
      endif

      asum = 0.0d0
      r12 = 1.0d0 / 12.0d0
      f3  = func(1) * rab(1) * r12
c      func(1) = 0.0d0

      do 100 i = 2,mesh-1,2
        f1 = f3
        f2 = func(i) * rab(i) * r12
        f3 = func(i+1) * rab(i+1) * r12
        asum = asum + 5.0d0*f1 + 8.0d0*f2 - 1.0d0*f3
c        func(i) = asum
        asum = asum - 1.0d0*f1 + 8.0d0*f2 + 5.0d0*f3
c        func(i+1) = asum
100   continue
      return
      end
