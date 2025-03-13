!
! Copyright (C) 2023-2035 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE fft_interpolation_mod
  USE kinds, ONLY: DP
  USE constants, ONLY: PI
  USE mp, ONLY: mp_sum, mp_get_comm_null 
  USE mp_bands, ONLY: intra_bgrp_comm
  IMPLICIT NONE 
  INTEGER  :: grp_comm = 0
  PRIVATE
  PUBLIC :: fft_1d_interpolate, fft_spherical_average, fft_2d_interpolate, set_grp_comm 
CONTAINS
  SUBROUTINE  set_grp_comm(comm)
    IMPLICIT NONE 
    INTEGER,INTENT(IN) :: comm 
    grp_comm = comm
  END SUBROUTINE set_grp_comm

  FUNCTION get_grp_comm() RESULT(comm)
    IMPLICIT NONE
    INTEGER :: comm  
    IF (grp_comm == 0) THEN 
      comm = intra_bgrp_comm
    ELSE 
      comm = grp_comm
    END IF  
  END FUNCTION 


  SUBROUTINE fft_1d_interpolate (data3D, gbase, points, data1D )
    IMPLICIT NONE
    COMPLEX(DP),INTENT(IN) :: data3D(:)
    !! data on 3D grid to be interpolated
    REAL(DP)  :: gbase(:,:)
    !! array of the plane waves vectonrs in units of 2PI/alat
    !! expected shape [3,ngm]  
    REAL(DP),INTENT(IN)  :: points(:,:)
    !! array of the real-space coordinates in alat units  for which compute the interpolated values
    !! expected shape [3, nx]
    COMPLEX(DP),INTENT(OUT) :: data1D(:)
    !! array where to copy interpolated values
    !! expected shape: [nx]
    !
    INTEGER  :: ngm, nx, i, ig 
    COMPLEX(DP) :: temp
    REAL(DP)    :: arg 
    ! 
    ngm = SIZE(gbase,2)
    nx =  SIZE(points,2) 
    ! 
    DO i = 1, nx
      temp = CMPLX(0._DP, 0._DP, KIND=DP)
      DO ig = 1, ngm
        arg = 2.d0 * PI * DOT_PRODUCT(gbase(:,ig),points(:,i))
        temp = temp + data3D (ig) * cmplx(cos(arg),sin(arg),kind=DP)
      ENDDO
      data1D(i) = temp
    ENDDO
    CALL mp_sum(data1D, get_grp_comm())
  END SUBROUTINE fft_1d_interpolate 

  SUBROUTINE fft_spherical_average(data3D, gbase, x0, deltax, data1D)
   !! spherically averaged charge: rho0(|r|) = int rho(r) dOmega
   !! rho0(r) = 4pi \sum_G rho(G) j_0(|G||r|)
   !
    IMPLICIT NONE
    COMPLEX(DP),INTENT(IN) :: data3D(:)
    !! array with 3D data to interpolate 
    !! expected shape [ngm]
    REAL(DP),INTENT(IN) :: gbase(:,:)
    !! array of the plane waves vectonrs in units of 2PI/alat
    !! expected shape [3,ngm]  
    REAL(DP),INTENT(IN) :: x0(3)
    !! coordinates in alat units of the center of the spherical average
    REAL(DP) :: deltax
    !! increment of sampling for the spherical average
    COMPLEX(DP),INTENT(OUT) :: data1D(:)
    !! array where to copy interpolated values
    !! expected shape: [nx]
    ! 
    INTEGER :: ngm, nx, gstart, ig, i 
    COMPLEX(DP) :: temp 
    REAL(DP)    :: arg 
    REAL(DP)    :: gg, gr
    !
    nx = SIZE(data1D)
    ngm = SIZE(gbase,2)
    gstart = 1  
    gg=sqrt(gbase(1,1)**2+gbase(2,1)**2+gbase(3,1)**2)
    IF (gg<1.d-10) THEN
      DO i = 1, nx
        data1D (i) = 4.d0 * pi * data3D (1)
      ENDDO
      gstart = 2
    ENDIF
     !     G!=0 terms
    DO ig = gstart, ngm
      arg = 2.d0 * PI * DOT_PRODUCT(x0,gbase(:,ig))
      !     This displaces the origin into x0
      temp  = data3D (ig) * cmplx(cos(arg),sin(arg),kind=DP)
      !     r =0 term
      data1D (1) = data1D (1) + 4.d0 * pi * temp 
      !     r!=0 terms
      DO i = 2, nx
        gr = 2.d0 * pi * sqrt(DOT_PRODUCT(gbase(:,ig),gbase(:,ig))) * (i-1) * deltax
        data1D (i) = data1D (i) + 4.d0 * PI * temp * sin (gr) / gr
      ENDDO
    ENDDO
    CALL mp_sum(data1D, get_grp_comm()) 
  END SUBROUTINE fft_spherical_average

  SUBROUTINE fft_2d_interpolate(data3d, gbase, nx, ny, x0, e1, e2, data2D)
    !! using data3D in reciprocal space computes the interpolated values 
    !! on a real space slab with sides e1 and e2 and origin x0
    IMPLICIT NONE    
    COMPLEX(DP),intent(in)  :: data3D(:)
    !! data in reciprocal space 
    !! expected shape [ngm]
    REAL(DP),intent(in)  :: gbase(:,:)
    !! basis wave vectors for FFT 
    !! expected format [3,ngm]
    INTEGER,intent(in) :: nx, ny
    !! number of step in direction 1 and 2 
    REAL(DP),intent(in)  :: x0(3), e1(3), e2(3)
    !! origin and the 2 spanning vectors of the surface of the interpolated points
    COMPLEX(DP),intent(out) :: data2D(:,:) 
    !! interpolated data in real space 
    !! expected shape [nx,ny]
    ! 
    REAL(DP)  :: deltax, deltay, m1, m2, v1(3),v2(3)
    INTEGER   :: ngm 
    INTEGER   :: i, j, ig
    COMPLEX(DP),ALLOCATABLE  :: eigx(:), eigy(:)
    !
    ALLOCATE (eigx(nx), eigy(ny))
    m1 = sqrt(dot_product(e1, e1))
    v1 = e1/m1  
    m2 = sqrt(dot_product(e2, e2)) 
    v2 = e2/m2
    deltax = m1/(nx - 1)
    deltay = m2/(ny - 1)
    ngm = SIZE(gbase,2) 
    data2D = cmplx(0._DP, 0._DP, kind = DP) 
    DO ig =1, ngm
      !
      ! eigx=exp(iG*e1+iGx0), eigy=(iG*e2)
      ! These factors are calculated and stored in order to save CPU time
      !
      DO i = 1, nx
        eigx (i) = exp((0.d0, 1.d0) * 2.d0 * pi * dot_product(gbase(:,ig), x0 + (i -1) * deltax * v1))
      END DO
      DO j = 1, ny
        eigy (j) = exp((0.d0, 1.d0) * 2.d0 * pi * (j - 1) * deltay * dot_product(v2, gbase(:, ig)))
      ENDDO
      DO j = 1, ny
        DO i = 1, nx
           data2D (i, j) = data2D (i, j) + data3D (ig) * eigx (i) * eigy (j)
        ENDDO
      ENDDO
    ENDDO
    CALL mp_sum(data2D, get_grp_comm())
  END SUBROUTINE fft_2d_interpolate
END MODULE fft_interpolation_mod
