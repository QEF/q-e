!
! Copyright (C) 2001-2022 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
! Author: Ivan Carnimeo (September 2021)
!----------------------------------------------------------------------------
!
!----------------------------------------------------------------------------
MODULE idwmod
  USE kinds,            ONLY : dp
  USE io_global,        ONLY : stdout
!
! An inverse distance weighting (idw) interpolation is computed here, using the metric proposed by Shepard 
! (ACM '68: Proceedings of the 1968 23rd ACM national conferenceJanuary 1968 Pages 517–524, 
!  https://doi.org/10.1145/800186.810616) 
! and a modified method (idw-sphere) that uses only nearest neighbors within R-sphere
!
implicit none
save
  !
  integer :: p_metric   ! metric for the (inverse) distance 
  !
  real(dp) :: scale_sphere  ! scaling factor for the radius of the sphere for the modified method
  !
CONTAINS
!----------------------------------------------------------------------------
subroutine idw(iwhat)
!
! compute the band structure with IDW interpolation
! iwhat = 1 ... basic IDW method
!         2 ... modified method with the sphere radius
!
USE globalmod, ONLY : Nb, Nq, q, eq, ek, at, bg
USE input_parameters, ONLY : nkstot, xk
implicit none
  integer, intent(in) :: iwhat 
  ! local variables
  real(dp) :: w, d, dsum, esum, dthr, R, Rtmp, Rmin, Rvec(3)
  integer :: ib, iq, jq, ik, NCount(2)
  !
  if (iwhat.ne.1.and.iwhat.ne.2) then
    write(stdout,'(A)') 'wrong iwhat in IDW method'
    stop
  else
    write(stdout,'(A)') ''
    write(stdout,'(A)') '--- Inverse distance weighting (IDW) interpolation method ---'
    write(stdout,'(A)') ''
    !write(stdout,'(4(A,I5))') 'iwhat: ',iwhat, ' Nb: ',Nb, ' Nq: ',Nq, ' Nk: ',nkstot
    if(p_metric.eq.1) then 
       write(stdout,'(A)') 'The inverse distance will be used in weights (p_metric = 1) ' 
    elseif(p_metric.eq.2) then 
       write(stdout,'(A)') 'The squared inverse distance will be used in weights (p_metric = 2)' 
    else 
       write(stdout,'(A,I5,A)') 'The inverse distance to the ', p_metric, 'th power will be used in weights' 
    end if 
  end if 
  !
  if(iwhat.eq.2) then 
    ! the radius of the search sphere is proportional to the uniform grid spacing 
    Rmin = 999999.0d0
    do iq = 1, Nq
      do jq = iq+1, Nq
        !
        ! Rtmp is the distance within the minimum image convention
        Rvec(:) = q(:,iq) - q(:,jq)
        CALL cryst_to_cart( 1, Rvec, at, 1 )
        Rvec(:) = Rvec(:) - ANINT( Rvec(:) )
        CALL cryst_to_cart( 1, Rvec, bg, -1 )
        Rtmp = sqrt( Rvec(1)**2 + Rvec(2)**2 + Rvec(3)**2 ) 
        !
        Rmin = min(Rmin, Rtmp)
        !
      end do 
    end do 
    R = scale_sphere * Rmin
    write(stdout,'(A,f12.6)') 'Input scale factor: ', scale_sphere 
    write(stdout,'(A,f12.6,A)') 'Minimum spacing between the uniform grid points: ', Rmin, ' (cart. coord. in units 2 pi/alat)'
    write(stdout,'(A,f12.6,A)') 'Scaled sphere radius: ', R, ' (cart. coord. in units 2 pi/alat)'
  end if 
  !
  dthr = 0.0000010d0
  !
  ek = 0.0d0
  !
  do ib = 1, Nb
    do ik = 1, nkstot
      !
      dsum = 0.0d0
      esum = 0.0d0
      NCount = 0 
      do iq = 1, Nq
        !
        ! d is the distance within the minimum image convention
        Rvec(:) = xk(:,ik) - q(:,iq)
        CALL cryst_to_cart( 1, Rvec, at, 1 )
        Rvec(:) = Rvec(:) - ANINT( Rvec(:) )
        CALL cryst_to_cart( 1, Rvec, bg, -1 )
        d = sqrt( Rvec(1)**2 + Rvec(2)**2 + Rvec(3)**2 ) 
        !
        if(d.gt.dthr) then
          NCount(1) = NCount(1) + 1 
          if(iwhat.eq.1) then 
            ! basic idw method
            w = 1.0d0/(d**p_metric)
          elseif(iwhat.eq.2) then
            ! search only inside the sphere R (idw-sphere)
            !w = (max(0.0d0, (R-d))/(R*d))**2
            w = 0.0d0
            if(d.lt.R) then 
              NCount(2) = NCount(2) + 1 
              w = 1.0d0/(d**p_metric) !((R-d)/(R*d))**2
            end if 
          end if
          dsum = dsum + w
          esum = esum + w * eq(iq,ib)
        else
          ek(ik,ib) = eq(iq,ib)
          NCount(1) = NCount(1) + 1 
          go to 10  
        end if
      end do 
      ek(ik,ib) = esum / dsum
      !
      if(dsum.lt.dthr) then 
        write(stdout,'(A,3f12.6)') 'no uniform grid points found for k-point:', xk(:,ik)
        write(stdout,'(A)')        'increase the search radius and check nosym=true in SCF ' 
        write(stdout,'(2I5, 3f12.6, 2I5)') ib, ik, esum, dsum, ek(ik, ib), NCount(:)
        Call errore( 'idw ', ' wrong sphere radius ', 1 )
      endif
      !
10    continue
      !
    end do      
  end do 
  !
  return
  !
end subroutine idw
!----------------------------------------------------------------------------
END MODULE
!----------------------------------------------------------------------------
