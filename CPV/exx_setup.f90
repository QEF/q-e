!======================================================================
! Set up the real-space grid for exact exchange.
! Define two spheres, possion solver is called for potential inside the inner one 
! Multipole expansion is used for points inside the outer one and all other points
! are set to zero.
!========================================================================

      subroutine getnpinsp( exx_ps_rcut, exx_me_rcut, np_in_sp, np_in_sp2 )
 
      USE kinds,      ONLY  : DP
      USE cell_base,  ONLY  : at, alat, h, s_to_r
      USE fft_base,   ONLY  : dffts

      implicit none
 
      REAL(DP)  exx_ps_rcut, exx_me_rcut
      INTEGER   np_in_sp, np_in_sp2
 
      REAL(DP)  s(3),r(3),dist, alength(3)
      integer   i,j,k
! --------------------------------------------------------------------
      alength(1) = sqrt( at(1,1)**2 + at(2,1)**2 + at(3,1)**2 ) * alat
      alength(2) = sqrt( at(1,2)**2 + at(2,2)**2 + at(3,2)**2 ) * alat
      alength(3) = sqrt( at(1,3)**2 + at(2,3)**2 + at(3,3)**2 ) * alat

!     print *, 'alength and step = ', alength, step

      np_in_sp = 0
      np_in_sp2= 0
       do k = 1, dffts%nr3
          do j = 1, dffts%nr2
             do i =1, dffts%nr1
            
                s(1) = DBLE(i)/DBLE(dffts%nr1) - 0.5d0
                s(2) = DBLE(j)/DBLE(dffts%nr2) - 0.5d0
                s(3) = DBLE(k)/DBLE(dffts%nr3) - 0.5d0

                call s_to_r(s, r, h)
                dist = sqrt(r(1)*r(1) + r(2)*r(2) + r(3)*r(3))
             
                if(dist .le. exx_ps_rcut)then
                   np_in_sp = np_in_sp + 1
                elseif(dist .le. exx_me_rcut)then
                   np_in_sp2 = np_in_sp2 + 1
                endif
             enddo
          enddo
       enddo

       write(6, *)'leaving getnp_in_sp', np_in_sp, np_in_sp2
       return
       end

!======================================================================
!
     subroutine exx_setup( nnrtot, lpole, clm, factor )

     USE kinds,                   ONLY  : DP
     USE cell_base,               ONLY  : at, alat, h, s_to_r
     USE fft_base,                ONLY  : dffts
! cutoff
     USE wannier_base,            ONLY :  exx_ps_rcut, exx_me_rcut
! order of expansion ( points on one side)
     USE cp_main_variables,       ONLY  : norder=>nord2, coeke
! number of grid points in the 1st sphere and the shell between 1st and 2nd sphere
     USE cp_main_variables,       ONLY  : np_in_sp, np_in_sp2
! conversion between 3D index (i,j,k) and 1D index 
! odthothd_in_sp(3, 1:np_in_sp) is for inner sphere 
! and odtothd_in_sp(3, np_in_sp+1: np_in_sp+np_in_sp2) is for the shell
     USE cp_main_variables,       ONLY  : odtothd_in_sp,  thdtood_in_sp
     USE cp_main_variables,       ONLY  : thdtood
     USE cp_main_variables,       ONLY  : xx_in_sp,  yy_in_sp,  zz_in_sp
     use cp_main_variables,       ONLY  : lap_neig, lap_dir_num, &
                                          lap_dir_step, b_lap, lap_dir

     USE parallel_include
     USE mp_global,               ONLY  : me_image

     implicit none


! ====================================================================
! INPUT VARIABLES
     integer  nnrtot, lpole
     real(DP) factor
     real(DP) clm(0:lpole, 0:lpole)

! WORK VARIABLES
       integer   np, npsp, npsp2
       real(DP)  s(3),r(3),dist, alength(3), step(3)
       integer   i,j,k, ierr, tmp
! --------------------------------------------------------------------

      thdtood_in_sp(:,:,:) = 0

      alength(1) = sqrt( at(1,1)**2 + at(2,1)**2 + at(3,1)**2 ) * alat
      alength(2) = sqrt( at(1,2)**2 + at(2,2)**2 + at(3,2)**2 ) * alat
      alength(3) = sqrt( at(1,3)**2 + at(2,3)**2 + at(3,3)**2 ) * alat

      step(1) = alength(1) / dffts%nr1
      step(2) = alength(2) / dffts%nr2
      step(3) = alength(3) / dffts%nr3
 
       np    = 0
       npsp  = 0
       npsp2 = 0
       do k = 1, dffts%nr3
          do j = 1, dffts%nr2
             do i =1, dffts%nr1
                np = np + 1
                thdtood(i,j,k) = np
!               odtothd(1,np) = i
!               odtothd(2,np) = j
!               odtothd(3,np) = k

                s(1) = DBLE(i)/DBLE(dffts%nr1) - 0.5d0
                s(2) = DBLE(j)/DBLE(dffts%nr2) - 0.5d0
                s(3) = DBLE(k)/DBLE(dffts%nr3) - 0.5d0

                call s_to_r(s, r, h)
                dist = sqrt(r(1)*r(1) + r(2)*r(2) + r(3)*r(3))

                if(dist .le. exx_ps_rcut)then
                   npsp = npsp + 1
                   thdtood_in_sp(i,j,k)  = npsp
                   odtothd_in_sp(1,npsp) = i
                   odtothd_in_sp(2,npsp) = j
                   odtothd_in_sp(3,npsp) = k

                   xx_in_sp(npsp) = r(1)
                   yy_in_sp(npsp) = r(2)
                   zz_in_sp(npsp) = r(3)
                elseif(dist .le. exx_me_rcut)then
                   npsp2 = npsp2 + 1
                   tmp = npsp2 + np_in_sp
                   thdtood_in_sp(i,j,k)  = tmp
                   odtothd_in_sp(1, tmp) = i
                   odtothd_in_sp(2, tmp) = j
                   odtothd_in_sp(3, tmp) = k

                   xx_in_sp(tmp) = r(1)
                   yy_in_sp(tmp) = r(2)
                   zz_in_sp(tmp) = r(3)

                endif
             enddo
          enddo
       enddo


       write(6,*)' npsp in exx_setup =', npsp, npsp2
       if( npsp .ne. np_in_sp)then
          write(6, *)'number of points in the 1st sphere does not match', npsp, np_in_sp
          write(6, *)'STOP in exx_setup'
          return
       endif

       if( npsp2 .ne. np_in_sp2)then
          write(6,*)'number of points in the 2nd sphere does not match', npsp2, np_in_sp2
          write(6, *)'STOP in exx_setup'
          return
       endif

!========================================================================
! non-orthogonal grids

      call exx_ggrid(h,step,1.d-8)
!========================================================================

!========================================================================
       call fornberg(2,norder,coeke(:,1),ierr)

       if (ierr .ne. 0) then
          write(6,*) ' ERROR: Wrong parameter in call of Fornberg'
          write(6,*) ' STOP in init_var'
          return
       endif

!  Renormalize coekes with respect to the grid spacing. First to the
!  new directions for non-orthogonal grid.
       if(lap_dir_num > 0) then
          write(6,*) ' lap_dir_num > 0', lap_dir_num
          do i = 1,lap_dir_num
             coeke(:,3+i) = -b_lap(3+i)*coeke(:,1)/(lap_dir_step(i)**2*factor)
          end do
       end if

       do i = 3, 1, -1
          coeke(:,i) = -b_lap(i)*coeke(:,1)/(step(i)*step(i)*factor)
       enddo

!==========================================================================

       call setclm(lpole, clm)

       return
       end

!======================================================================
!
     subroutine exx_setup_nscf( nnrtot, lpole, clm, factor, wc, vwc, nbsp, vnbsp )

     USE kinds,                   ONLY  : DP
     USE cell_base,               ONLY  : at, alat, h, s_to_r 
     USE fft_base,                ONLY  : dffts
! cutoff
     USE wannier_base,            ONLY :  exx_ps_rcut, exx_me_rcut
! order of expansion ( points on one side)
     USE cp_main_variables,       ONLY  : norder=>nord2, coeke
! number of grid points in the 1st sphere and the shell between 1st and 2nd sphere
     USE cp_main_variables,       ONLY  : np_in_sp, np_in_sp2
! conversion between 3D index (i,j,k) and 1D index 
! odthothd_in_sp(3, 1:np_in_sp) is for inner sphere 
! and odtothd_in_sp(3, np_in_sp+1: np_in_sp+np_in_sp2) is for the shell
     USE cp_main_variables,       ONLY  : odtothd_in_sp,  thdtood_in_sp
     USE cp_main_variables,       ONLY  : thdtood
     USE cp_main_variables,       ONLY  : xx_in_sp,  yy_in_sp,  zz_in_sp
     use cp_main_variables,       ONLY  : lap_neig, lap_dir_num, &
                                          lap_dir_step, b_lap, lap_dir
     USE parallel_include
     USE mp_global,               ONLY  : me_image

     implicit none

! ====================================================================
! INPUT VARIABLES
     integer  nnrtot, lpole, nbsp, vnbsp
     real(DP) factor, wc(3, nbsp), vwc(3, vnbsp)
     real(DP) clm(0:lpole, 0:lpole)

! ====================================================================
!      integer  odtothd_in_sp(3,np_in_sp), odtothd_in_sp2(3,np_in_sp2)
!      integer  thdtood_in_sp(nr1s, nr2s, nr3s),thdtood(nr1s,nr2s,nr3s)
! ====================================================================
! WORK VARIABLES
       integer   np, npsp, npsp2
       real(DP)  s(3), r(3),dist, alength(3), step(3)
       integer   i,j,k, ierr, tmp
! --------------------------------------------------------------------

! For points in the 1st sphere, one needs to know if its finite-difference neighbors
! are inside or outside. We set the thdtood_in_sp to be np_in_sp + 1 for outside neighbors
       thdtood_in_sp(:,:,:) = 0

       alength(1) = sqrt( at(1,1)**2 + at(2,1)**2 + at(3,1)**2 ) * alat
       alength(2) = sqrt( at(1,2)**2 + at(2,2)**2 + at(3,2)**2 ) * alat
       alength(3) = sqrt( at(1,3)**2 + at(2,3)**2 + at(3,3)**2 ) * alat
 
       step(1) = alength(1) / dffts%nr1
       step(2) = alength(2) / dffts%nr2
       step(3) = alength(3) / dffts%nr3
 
       np    = 0
       npsp  = 0
       npsp2 = 0
       do k = 1, dffts%nr3
          do j = 1, dffts%nr2
             do i =1, dffts%nr1
                np = np + 1
                thdtood(i,j,k) = np
!               odtothd(1,np) = i
!               odtothd(2,np) = j
!               odtothd(3,np) = k
             
                s(1) = DBLE(i)/DBLE(dffts%nr1) - 0.5d0
                s(2) = DBLE(j)/DBLE(dffts%nr2) - 0.5d0
                s(3) = DBLE(k)/DBLE(dffts%nr3) - 0.5d0

                call s_to_r(s, r, h)
                dist = sqrt(r(1)*r(1) + r(2)*r(2) + r(3)*r(3))

                if(dist .le. exx_ps_rcut)then
                   npsp = npsp + 1
                   thdtood_in_sp(i,j,k)  = npsp
                   odtothd_in_sp(1,npsp) = i
                   odtothd_in_sp(2,npsp) = j
                   odtothd_in_sp(3,npsp) = k
             
                   xx_in_sp(npsp) = r(1)
                   yy_in_sp(npsp) = r(2)
                   zz_in_sp(npsp) = r(3)
                elseif(dist .le. exx_me_rcut)then
                   npsp2 = npsp2 + 1
                   tmp = npsp2 + np_in_sp
                   thdtood_in_sp(i,j,k)  = tmp
                   odtothd_in_sp(1, tmp) = i
                   odtothd_in_sp(2, tmp) = j
                   odtothd_in_sp(3, tmp) = k

                   xx_in_sp(tmp) = r(1)
                   yy_in_sp(tmp) = r(2)
                   zz_in_sp(tmp) = r(3)

                endif
             enddo
          enddo
       enddo

       write(6,*)' npsp in exx_setup =', npsp, npsp2
       if( npsp .ne. np_in_sp)then
          write(6, *)'number of points in the 1st sphere does not match', npsp, np_in_sp
          write(6, *)'STOP in exx_setup'
          return
       endif

       if( npsp2 .ne. np_in_sp2)then
          write(6,*)'number of points in the 2nd sphere does not match', npsp2, np_in_sp2
          write(6, *)'STOP in exx_setup'
          return
       endif

       call exx_ggrid(h,step,1.d-8)

!========================================================================
       call fornberg(2,norder,coeke(:,1),ierr)
    
       if (ierr .ne. 0) then
          write(6,*) ' ERROR: Wrong parameter in call of Fornberg'
          write(6,*) ' STOP in init_var'
          return
       endif
    
!  Renormalize coekes with respect to the grid spacing. First to the
!  new directions for non-orthogonal grid.
       if(lap_dir_num > 0) then
          write(6,*) ' lap_dir_num > 0', lap_dir_num
          do i = 1,lap_dir_num
             coeke(:,3+i) = -b_lap(3+i)*coeke(:,1)/(lap_dir_step(i)**2*factor)
          end do
       end if

       do i = 3, 1, -1
          coeke(:,i) = -b_lap(i)*coeke(:,1)/(step(i)*step(i)*factor)
       enddo
!==========================================================================

       call setclm(lpole, clm)

       call getwc(wc, vwc, vnbsp, nbsp)

       return
       end


       subroutine getwc(wc, vwc, vnbsp, nbsp)
       USE kinds,                   ONLY  : DP
       USE cell_base,               ONLY  : h, ainv, s_to_r, r_to_s, pbcs

       IMPLICIT NONE

       INTEGER     vnbsp, nbsp, i, ir
       REAl(DP)    wc(3, nbsp), vwc(3, vnbsp), tmp(3)

       do ir = 1, nbsp
          read(407, *) tmp(1), tmp(2), tmp(3)
          call r_to_s(tmp,wc(:,ir), ainv)
!         call pbcs(wc(:,ir),tmp,1)
          do i = 1, 3
             tmp(i) = wc(i,ir) - int( wc(i,ir))
             if(tmp(i) < 0)then
                tmp(i) = tmp(i) + 1
             endif
          enddo

          call s_to_r(tmp,wc(:,ir), h)
       end do

       do ir = 1, vnbsp
          read(408,*)tmp(1), tmp(2), tmp(3)
          call r_to_s(tmp, vwc(:,ir), ainv)
!         call pbcs(  vwc(:,ir), tmp, 1)
          do i = 1, 3
             tmp(i) = vwc(i,ir) - int( vwc(i,ir))
             if(tmp(i) < 0)then
                tmp(i) = tmp(i) + 1
             endif
          enddo

          call s_to_r(tmp, vwc(:,ir), h)
       enddo

       RETURN
       END
! ==================================================================
! subroutine to set the various clm coefficients. Separated so as to
! not clutter up the above code. calculates clm = (l-m)!/(l+m)! for 
! l = 0,lpole, m=0,l. 
! evaluating these to 20 decimal places is certainly as accurate 
! as doing the calculations in the code but less work.
! ==================================================================
       subroutine setclm(lpole, clm)

       implicit none
! INPUT: order of mutipole expansion
       integer lpole
! OUTPUT: clm coefficients
       real*8 clm(0:lpole, 0:lpole)

       clm(0,0) = 1.00000000000000000000e+00
       if (lpole .ge. 1) then
         clm(1,0) = 1.00000000000000000000e+00
         clm(1,1) = 5.00000000000000000000e-01
       endif
       if (lpole .ge. 2) then
         clm(2,0) = 1.00000000000000000000e+00
         clm(2,1) = 1.66666666666666666670e-01
         clm(2,2) = 4.16666666666666666670e-02
       endif
       if (lpole .ge. 3) then
         clm(3,0) = 1.00000000000000000000e+00
         clm(3,1) = 8.33333333333333333330e-02
         clm(3,2) = 8.33333333333333333330e-03
         clm(3,3) = 1.38888888888888888890e-03
       endif
       if (lpole .ge. 4) then
         clm(4,0) = 1.00000000000000000000e+00
         clm(4,1) = 5.00000000000000000000e-02
         clm(4,2) = 2.77777777777777777780e-03
         clm(4,3) = 1.98412698412698412700e-04
         clm(4,4) = 2.48015873015873015870e-05
       endif
       if (lpole .ge. 5) then
         clm(5,0) = 1.00000000000000000000e+00
         clm(5,1) = 3.33333333333333333330e-02
         clm(5,2) = 1.19047619047619047620e-03
         clm(5,3) = 4.96031746031746031750e-05
         clm(5,4) = 2.75573192239858906530e-06
         clm(5,5) = 2.75573192239858906530e-07
       endif
       if (lpole .ge. 6) then
         clm(6,0) = 1.00000000000000000000e+00
         clm(6,1) = 2.38095238095238095240e-02
         clm(6,2) = 5.95238095238095238100e-04
         clm(6,3) = 1.65343915343915343920e-05
         clm(6,4) = 5.51146384479717813050e-07
         clm(6,5) = 2.50521083854417187750e-08
         clm(6,6) = 2.08767569878680989790e-09
       endif
       if (lpole .ge. 7) then
         clm(7,0) = 1.00000000000000000000e+00
         clm(7,1) = 1.78571428571428571430e-02
         clm(7,2) = 3.30687830687830687830e-04
         clm(7,3) = 6.61375661375661375660e-06
         clm(7,4) = 1.50312650312650312650e-07
         clm(7,5) = 4.17535139757361979580e-09
         clm(7,6) = 1.60590438368216145990e-10
         clm(7,7) = 1.14707455977297247140e-11
       endif
       if (lpole .ge. 8) then
         clm(8,0) = 1.00000000000000000000e+00
         clm(8,1) = 1.38888888888888888890e-02
         clm(8,2) = 1.98412698412698412700e-04
         clm(8,3) = 3.00625300625300625300e-06
         clm(8,4) = 5.01042167708834375500e-08
         clm(8,5) = 9.63542630209296875960e-10
         clm(8,6) = 2.29414911954594494280e-11
         clm(8,7) = 7.64716373181981647590e-13
         clm(8,8) = 4.77947733238738529740e-14
       endif
       if (lpole .ge. 9) then
         clm(9,0) = 1.00000000000000000000e+00
         clm(9,1) = 1.11111111111111111110e-02
         clm(9,2) = 1.26262626262626262630e-04
         clm(9,3) = 1.50312650312650312650e-06
         clm(9,4) = 1.92708526041859375190e-08
         clm(9,5) = 2.75297894345513393130e-10
         clm(9,6) = 4.58829823909188988550e-12
         clm(9,7) = 9.55895466477477059490e-14
         clm(9,8) = 2.81145725434552076320e-15
         clm(9,9) = 1.56192069685862264620e-16
       endif

       return
       end 

!     ===============================================================
!
!     Coefficients for the first & second order numerical derivative
!     under the centered finite difference scheme.
!     Bengt Fornberg,  Exxon Res. & Eng. Co., NJ 08801
!     'bfornbe@erenj.com'
!     David M. Sloan,  Dept. of Mathematics, U. of Strathclyde,
!     Glasgow G1 1HX, Scotland,  'd.sloan@strath.ac.uk'
!     Acta Numerica 94,  Cambridge Univ. Press (1994)
!
!     ---------------------------------------------------------------
      subroutine fornberg(iddd,norder,coe,ierr)

      use kinds, only : DP
      implicit none
!
!     Input/Output variables:
!
!     order of expansion of derivative. 
!     it is the number of neighbors used ON ONE SIDE.
!     the maximum order implemented is 20.
      integer, intent(in) :: norder

!     iddd - order of the derivative (iddd = 1 or 2)
      integer, intent(in) :: iddd

!     coe - coefficients for the derivative
      real(dp), intent(out) :: coe(-norder:norder)

!     ierr - error flag
      integer, intent(out) :: ierr
!     
!     Work variables:
!
!     counters 
      integer i
!     ---------------------------------------------------------------
!
!     First order derivative
!
      ierr = 0
      if(iddd.eq.1) then
!
         select case (norder)
         case (1)
            coe(1) =  0.50000000000000D+00
         case (2)
            coe(1) =  2.d0/3.d0
            coe(2) = -1.d0/12.d0
         case (3)
            coe(1) =  3.d0/4.d0
            coe(2) = -3.d0/20.d0
            coe(3) =  1.d0/60.d0
         case (4)   
            coe(1) =  4.d0/5.d0
            coe(2) = -1.d0/5.d0
            coe(3) =  4.d0/105.d0
            coe(4) = -1.d0/280.d0
         case (5)     
            coe(1) =  0.8333333333D+00
            coe(2) = -0.2380952381D+00
            coe(3) =  0.5952380952D-01
            coe(4) = -0.9920634921D-02
            coe(5) =  0.7936507937D-03
         case (6)    
            coe(1) =  0.8571428571D+00
            coe(2) = -0.2678571429D+00
            coe(3) =  0.7936507937D-01
            coe(4) = -0.1785714286D-01
            coe(5) =  0.2597402597D-02
            coe(6) = -0.1803751804D-03
         case (7)    
            coe(1) =  0.8750000000D+00
            coe(2) = -0.2916666667D+00
            coe(3) =  0.9722222222D-01
            coe(4) = -0.2651515152D-01
            coe(5) =  0.5303030303D-02
            coe(6) = -0.6798756799D-03
            coe(7) =  0.4162504163D-04
         case (8)    
            coe(1) =  0.8888888889D+00
            coe(2) = -0.3111111111D+00
            coe(3) =  0.1131313131D+00
            coe(4) = -0.3535353535D-01
            coe(5) =  0.8702408702D-02
            coe(6) = -0.1554001554D-02
            coe(7) =  0.1776001776D-03
            coe(8) = -0.9712509713D-05
         case (9)      
            coe(1) =  0.9000000000D+00
            coe(2) = -0.3272727273D+00
            coe(3) =  0.1272727273D+00
            coe(4) = -0.4405594406D-01
            coe(5) =  0.1258741259D-01
            coe(6) = -0.2797202797D-02
            coe(7) =  0.4495504496D-03
            coe(8) = -0.4627725216D-04
            coe(9) =  0.2285296403D-05
         case (10)    
            coe(1) =  0.9090909091D+00
            coe(2) = -0.3409090909D+00
            coe(3) =  0.1398601399D+00
            coe(4) = -0.5244755245D-01
            coe(5) =  0.1678321678D-01
            coe(6) = -0.4370629371D-02
            coe(7) =  0.8814714697D-03
            coe(8) = -0.1285479227D-03
            coe(9) =  0.1202787580D-04
            coe(10)= -0.5412544112D-06
         end select

         coe(0) = 0.d0
         do i = 1,norder
            coe(-i) = -coe(i)
         enddo
!
!     Second order derivative
!
      else if (iddd.eq.2) then

         select case (norder)
         case (1)
            coe(0) = -0.20000000000000D+01
            coe(1) =  0.10000000000000D+01
         case (2) 
            coe(0) = -0.25000000000000D+01
            coe(1) =  0.13333333333333D+01
            coe(2) = -0.83333333333333D-01
         case (3)
            coe(0) = -0.27222222222222D+01
            coe(1) =  0.15000000000000D+01
            coe(2) = -0.15000000000000D+00
            coe(3) =  0.11111111111111D-01
         case (4)
            coe(0) = -0.28472222222222D+01
            coe(1) =  0.16000000000000D+01
            coe(2) = -0.20000000000000D+00
            coe(3) =  0.25396825396825D-01
            coe(4) = -0.17857142857143D-02
         case (5)
            coe(0) = -0.29272222222222D+01
            coe(1) =  0.16666666666667D+01
            coe(2) = -0.23809523809524D+00
            coe(3) =  0.39682539682540D-01
            coe(4) = -0.49603174603175D-02
            coe(5) =  0.31746031746032D-03
         case (6)
            coe(0) = -0.29827777777778D+01
            coe(1) =  0.17142857142857D+01
            coe(2) = -0.26785714285714D+00
            coe(3) =  0.52910052910053D-01
            coe(4) = -0.89285714285714D-02
            coe(5) =  0.10389610389610D-02
            coe(6) = -0.60125060125060D-04
         case (7)
            coe(0) = -0.30235941043084D+01
            coe(1) =  0.17500000000000D+01
            coe(2) = -0.29166666666667D+00
            coe(3) =  0.64814814814815D-01
            coe(4) = -0.13257575757576D-01
            coe(5) =  0.21212121212121D-02
            coe(6) = -0.22662522662523D-03
            coe(7) =  0.11892869035726D-04
         case (8)
            coe(0) = -0.30548441043084D+01
            coe(1) =  0.17777777777778D+01
            coe(2) = -0.31111111111111D+00
            coe(3) =  0.75420875420875D-01
            coe(4) = -0.17676767676768D-01
            coe(5) =  0.34809634809635D-02
            coe(6) = -0.51800051800052D-03
            coe(7) =  0.50742907885765D-04
            coe(8) = -0.24281274281274D-05
         case (9)
            coe(0) = -0.30795354623331D+01
            coe(1) =  0.18000000000000D+01
            coe(2) = -0.32727272727273D+00
            coe(3) =  0.84848484848485D-01
            coe(4) = -0.22027972027972D-01
            coe(5) =  0.50349650349650D-02
            coe(6) = -0.93240093240093D-03
            coe(7) =  0.12844298558584D-03
            coe(8) = -0.11569313039901D-04
            coe(9) =  0.50784364509855D-06
         case (10)
            coe(0) = -0.30995354623331D+01
            coe(1) =  0.18181818181818D+01
            coe(2) = -0.34090909090909D+00
            coe(3) =  0.93240093240093D-01
            coe(4) = -0.26223776223776D-01
            coe(5) =  0.67132867132867D-02
            coe(6) = -0.14568764568765D-02
            coe(7) =  0.25184899134479D-03
            coe(8) = -0.32136980666392D-04
            coe(9) =  0.26728612899924D-05
            coe(10)= -0.10825088224469D-06
         end select

         do i = 1,norder
            coe(-i) = coe(i)
         end do
!
      else
         write(6,*) ' ERROR: invalid derivative order, iddd = ',iddd
         write(6,*) ' STOP in FORNBERG '
         ierr = 1 
         return
      endif

      return
      end
