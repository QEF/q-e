!
! Copyright (C) 2002 CP90 group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

module stre
  implicit none 
  save
  real(kind=8) stress(3,3)
end module stre

module dqrad_mod
  implicit none 
  save
  real(kind=8),allocatable:: dqrad(:,:,:,:,:,:,:)
contains
  subroutine deallocate_dqrad_mod
      IF( ALLOCATED( dqrad ) ) DEALLOCATE( dqrad )
  end subroutine
end module dqrad_mod

module betax
  implicit none 
  save
  integer, parameter:: mmx=5001
  real(kind=8) :: refg
  real(kind=8),allocatable:: betagx(:,:,:), dbetagx(:,:,:), &
                       qradx(:,:,:,:,:), dqradx(:,:,:,:,:)
contains
  subroutine deallocate_betax
      IF( ALLOCATED( betagx ) ) DEALLOCATE( betagx )
      IF( ALLOCATED( dbetagx ) ) DEALLOCATE( dbetagx )
      IF( ALLOCATED( qradx ) ) DEALLOCATE( qradx )
      IF( ALLOCATED( dqradx ) ) DEALLOCATE( dqradx )
  end subroutine
end module betax

module cpr_subroutines

  implicit none
  save

contains

  subroutine compute_stress( stress, detot, h, omega )
    real(kind=8) :: stress(3,3), detot(3,3), h(3,3), omega
    integer :: i, j
         do i=1,3
            do j=1,3
               stress(i,j)=-1.d0/omega*(detot(i,1)*h(j,1)+              &
     &                      detot(i,2)*h(j,2)+detot(i,3)*h(j,3))
            enddo
         enddo
    return
  end subroutine

  subroutine print_atomic_var( var, na, nsp, head, iunit )
    use io_global, only: stdout
    real(kind=8) :: var(:,:,:)
    integer :: na(:), nsp
    integer, optional :: iunit
    character(len=*), optional :: head
    integer :: i, ia, is, iu
    if( present( iunit ) ) then
      iu = iunit
    else
      iu = stdout
    end if
    if( present( head ) ) then 
      WRITE( iu,*) head
    end if
    WRITE( iu,'(3f14.8)') (((var(i,ia,is),i=1,3),ia=1,na(is)),is=1,nsp)
    return
  end subroutine

  subroutine print_cell_var( var, head, iunit )
    use io_global, only: stdout
    real(kind=8) :: var(3,3)
    integer, optional :: iunit
    character(len=*), optional :: head
    integer :: i, j, iu
    if( present( iunit ) ) then
      iu = iunit
    else
      iu = stdout
    end if
    if( present( head ) ) then 
      WRITE( iu,*)
      WRITE( iu,*) head
      WRITE( iu, 5555 ) ((var(i,j),j=1,3),i=1,3)
 5555    format(1x,f12.5,1x,f12.5,1x,f12.5/                             &
     &          1x,f12.5,1x,f12.5,1x,f12.5/                             &
     &          1x,f12.5,1x,f12.5,1x,f12.5//)
    else
      write(iu,3340) ((var(i,j),i=1,3),j=1,3)
 3340     format(9(1x,f9.5))
    end if
    return
  end subroutine


  subroutine cell_hmove( h, hold, delt, omega, press, iforceh, stress, ainv )
    use cell_base, only: wmass
    real(kind=8), intent(out) :: h(3,3)
    real(kind=8), intent(in) :: hold(3,3), stress(3,3), ainv(3,3)
    real(kind=8), intent(in) :: omega, delt, press
    integer, intent(in) :: iforceh(3,3)
    real(kind=8) :: dt2by2, fac
    integer :: i, j
    dt2by2 = .5d0 * delt * delt
    fac = dt2by2/wmass*omega
    do i=1,3
      do j=1,3
        h(i,j) = hold(i,j) + fac * iforceh(i,j) *     &
               (stress(i,1)*ainv(j,1)+stress(i,2)*ainv(j,2)+   &
                stress(i,3)*ainv(j,3)-press*ainv(j,i))
      end do
    end do
    return
  end subroutine
  

  subroutine ions_hmove( taus, tausm, iforce, pmass, fion, ainv, delt, na, nsp )
    real(kind=8), intent(in) :: tausm(:,:,:), pmass(:), fion(:,:,:)
    integer, intent(in) :: iforce(:,:,:)
    real(kind=8), intent(in) :: ainv(3,3), delt
    real(kind=8), intent(out) :: taus(:,:,:)
    integer, intent(in) :: na(:), nsp
    integer :: is, ia, i
    real(kind=8) :: dt2by2, fac, fions(3)

    dt2by2 = .5d0 * delt * delt

    do is=1,nsp
      fac = dt2by2/pmass(is)
      do ia=1,na(is)
        do i=1,3
          fions( i ) = fion(1,ia,is)*ainv(i,1) + fion(2,ia,is)*ainv(i,2) + fion(3,ia,is)*ainv(i,3)
        end do
        do i=1,3
          taus(i,ia,is) = tausm(i,ia,is) + iforce(i,ia,is) * fac * fions( i )
        end do
      end do
    end do
    return
  end subroutine

 
  subroutine add_thermal_stress( stress, pmass, omega, h, vels, nsp, na )
    real(kind=8) :: stress(3,3)
    real(kind=8), intent(in) :: pmass(:), omega, h(3,3), vels(:,:,:)
    integer, intent(in) :: nsp, na(:)
    integer :: i, j, is, ia
    do is=1,nsp
      do ia=1,na(is)
        do i=1,3
          do j=1,3
            stress(i,j)=stress(i,j)+pmass(is)/omega*           &
     &        ((h(i,1)*vels(1,ia,is)+h(i,2)*vels(2,ia,is)+    &
     &          h(i,3)*vels(3,ia,is))*(h(j,1)*vels(1,ia,is)+  &
     &          h(j,2)*vels(2,ia,is)+h(j,3)*vels(3,ia,is)))
          enddo
        enddo
      enddo
    enddo
    return
  end subroutine
 

  subroutine cell_move( hnew, h, hold, delt, omega, press, iforceh, stress, ainv, &
                        frich, tnoseh, vnhh, velh )
    use cell_base, only: wmass, cell_verlet
    real(kind=8), intent(out) :: hnew(3,3)
    real(kind=8), intent(in) :: h(3,3), hold(3,3), stress(3,3), ainv(3,3)
    real(kind=8), intent(in) :: omega, press
    real(kind=8), intent(in) :: vnhh(3,3), velh(3,3)
    integer,      intent(in) :: iforceh(3,3)
    real(kind=8), intent(in) :: frich, delt
    logical,      intent(in) :: tnoseh

    real(kind=8) :: hnos(3,3), fcell(3,3)
    integer      :: i, j
   
    if( tnoseh ) then
      hnos = vnhh * velh
    else
      hnos = 0.0d0
    end if

    do j=1,3
      do i=1,3
        fcell(i,j) = ainv(j,1)*stress(i,1) + ainv(j,2)*stress(i,2) + ainv(j,3)*stress(i,3)
      end do
    end do
    do j=1,3
      do i=1,3
        fcell(i,j) = fcell(i,j) - ainv(j,i) * press
      end do
    end do
    fcell = omega * fcell / wmass
!
    call cell_verlet( hnew, h, hold, delt, iforceh, fcell, frich, tnoseh, hnos )

    return
  end subroutine


  subroutine cell_gamma( hgamma, ainv, h, velh )
    implicit none
    real(kind=8) :: hgamma(3,3)
    real(kind=8), intent(in) :: ainv(3,3), h(3,3), velh(3,3)
    integer :: i,j,k,l,m
         do i=1,3
            do j=1,3
               do k=1,3
                  do l=1,3
                     do m=1,3
                        hgamma(i,j)=hgamma(i,j)+ainv(i,l)*ainv(k,l)*    &
     &                       (velh(m,k)*h(m,j)+h(m,k)*velh(m,j))
                     enddo
                  enddo
               enddo
            enddo
         enddo
    return
  end subroutine

  subroutine ions_move( tausp, taus, tausm, iforce, pmass, fion, ainv, delt, na, nsp, &
                        fricp, hgamma, vels, tsdp, tnosep, fionm, vnhp, velsp, velsm )
    implicit none
    real(kind=8), intent(in) :: taus(:,:,:), tausm(:,:,:), pmass(:), fion(:,:,:)
    integer, intent(in) :: iforce(:,:,:)
    real(kind=8), intent(in) :: ainv(3,3), delt
    real(kind=8), intent(out) :: tausp(:,:,:)
    integer, intent(in) :: na(:), nsp
    real(kind=8), intent(in) :: fricp, hgamma(3,3), vels(:,:,:)
    logical, intent(in) :: tsdp, tnosep
    real(kind=8), intent(inout) :: fionm(:,:,:)
    real(kind=8), intent(in) :: vnhp
    real(kind=8), intent(out) :: velsp(:,:,:)
    real(kind=8), intent(in) :: velsm(:,:,:)
    integer :: is, ia, i
    real(kind=8) :: dt2by2, fac, fions(3), dt2, twodel
    real(kind=8) :: verl1, verl2, verl3

    dt2by2 = .5d0 * delt * delt
    dt2    = delt * delt
    twodel = 2.0d0 * delt

         verl1=2./(1.+fricp)
         verl2=1.-verl1
         verl3=dt2/(1.+fricp)
!
         if(tsdp) then
            do is=1,nsp
               do ia=1,na(is)
                  do i=1,3
                     tausp(i,ia,is) = taus(i,ia,is) +                   &
     &                    iforce(i,ia,is)*dt2by2/pmass(is)*             &
     &        (ainv(i,1)*fion(1,ia,is)+ainv(i,2)*fion(2,ia,is)+         &
     &         ainv(i,3)*fion(3,ia,is) ) -                              &
     &                    pmass(is)*(hgamma(i,1)*vels(1,ia,is)+         &
     &         hgamma(i,2)*vels(2,ia,is)+hgamma(i,3)*vels(3,ia,is))
                  end do
               end do
            end do
         else if (tnosep) then
            do is=1,nsp
               do ia=1,na(is)
                  do i=1,3
                     fionm(i,ia,is) = (ainv(i,1)*fion(1,ia,is)          &
     &                                +ainv(i,2)*fion(2,ia,is)          &
     &                                +ainv(i,3)*fion(3,ia,is))         &
     &                              - vnhp*vels(i,ia,is)*pmass(is)      &
     &                    - pmass(is)*(hgamma(i,1)*vels(1,ia,is)        &
     &                                +hgamma(i,2)*vels(2,ia,is)        &
     &                                +hgamma(i,3)*vels(3,ia,is))
                     tausp(i,ia,is)=-tausm(i,ia,is)+2.*taus(i,ia,is)+   &
     &                   iforce(i,ia,is)*dt2*fionm(i,ia,is)/pmass(is)
                     velsp(i,ia,is) = velsm(i,ia,is) +                  &
     &                    twodel*fionm(i,ia,is)/pmass(is)
                  end do
               end do
            end do
         else
            do is=1,nsp
               do ia=1,na(is)
                  do i=1,3
                     tausp(i,ia,is) = verl1*taus(i,ia,is)               &
     &                    + verl2*tausm(i,ia,is)                        &
     &        + verl3/pmass(is)*iforce(i,ia,is) * (ainv(i,1)*fion(1,ia,is)&
     &        + ainv(i,2)*fion(2,ia,is) + ainv(i,3)*fion(3,ia,is))      &
     &        - verl3*iforce(i,ia,is) * (hgamma(i,1)*vels(1,ia,is)      &
     &        + hgamma(i,2)*vels(2,ia,is) + hgamma(i,3)*vels(3,ia,is))
                     velsp(i,ia,is)=velsm(i,ia,is)                      &
     &        - 4.*fricp*vels(i,ia,is)                                  &
     &        + twodel/pmass(is)*iforce(i,ia,is)*(ainv(i,1)*fion(1,ia,is) &
     &        + ainv(i,2)*fion(2,ia,is) + ainv(i,3)*fion(3,ia,is))      &
     &        - twodel*iforce(i,ia,is) * (hgamma(i,1)*vels(1,ia,is)     &
     &        + hgamma(i,2)*vels(2,ia,is) + hgamma(i,3)*vels(3,ia,is))
                  end do
               end do
            end do
         endif
    return
  end subroutine


  subroutine ions_cofmsub( tausp, na, nsp, cdm, cdm0 )
    implicit none
    real( kind=8 ), intent(inout) :: tausp( :, :, : )
    integer, intent(in) :: na(:), nsp
    real( kind=8 ), intent(in) :: cdm( : ), cdm0( : )
    integer :: i, ia, is
    do is=1,nsp
      do ia=1,na(is)
        do i=1,3
          tausp(i,ia,is)=tausp(i,ia,is)+cdm0(i)-cdm(i)
        enddo
      enddo
    enddo
    return
  end subroutine


  subroutine ions_kinene( ekinp, vels, na, nsp, hold, pmass )
    implicit none
    real( kind=8 ), intent(out) :: ekinp
    real( kind=8 ), intent(in) :: vels(:,:,:)
    real( kind=8 ), intent(in) :: pmass(:)
    real( kind=8 ), intent(in) :: hold(:,:)
    integer, intent(in) :: na(:), nsp
    integer :: i, j, is, ia, ii
    ekinp = 0.0d0
    do is=1,nsp
      do ia=1,na(is)
        do i=1,3
          do j=1,3
            do ii=1,3
              ekinp=ekinp+pmass(is)* hold(j,i)*vels(i,ia,is)* hold(j,ii)*vels(ii,ia,is)
            end do
          end do
        end do
      end do
    end do
    ekinp=0.5d0*ekinp
    return
  end subroutine


  subroutine ions_temp( tempp, ekinpr, vels, na, nsp, hold, pmass, cdmvel )
    use constants, only: factem
    implicit none
    real( kind=8 ), intent(out) :: ekinpr, tempp
    real( kind=8 ), intent(in) :: vels(:,:,:)
    real( kind=8 ), intent(in) :: pmass(:)
    real( kind=8 ), intent(in) :: hold(:,:)
    real( kind=8 ), intent(in) :: cdmvel(:)
    integer, intent(in) :: na(:), nsp
    integer :: nat, i, j, is, ia, ii
    nat = SUM( na(1:nsp) )
    do i=1,3
      do j=1,3
        do ii=1,3
          do is=1,nsp
            do ia=1,na(is)
              ekinpr=ekinpr+pmass(is)*hold(j,i)*              &
     &              (vels(i,ia,is)-cdmvel(i))*                 &
     &               hold(j,ii)*(vels(ii,ia,is)-cdmvel(ii))
            end do
          end do
        end do
      end do
    end do
    ekinpr=0.5*ekinpr
    tempp=ekinpr*factem/(1.5d0*nat)
    return
  end subroutine


  subroutine elec_fakekine( ekincm, ema0bg, emass, c0, cm, ngw, n, delt )
    use mp, only: mp_sum
    use reciprocal_vectors, only: gstart
    use wave_base, only: wave_speed2
    real(kind=8), intent(out) :: ekincm
    real(kind=8), intent(in)  :: ema0bg(:), delt, emass
    complex(kind=8), intent(in)  :: c0(:,:,:,:), cm(:,:,:,:)
    integer, intent(in) :: ngw, n
    real(kind=8), allocatable :: emainv(:)
    real(kind=8) :: ftmp
    integer :: i

    ALLOCATE( emainv( ngw ) )
    emainv = 1.0d0 / ema0bg
    ftmp = 1.0d0
    if( gstart == 2 ) ftmp = 0.5d0

    ekincm=0.0d0
    do i=1,n
      ekincm = ekincm + 2.0d0 * &
               wave_speed2( c0(:,i,1,1), cm(:,i,1,1), emainv, ftmp )
    end do
    ekincm = ekincm * emass / ( delt * delt )

    CALL mp_sum( ekincm )
    DEALLOCATE( emainv )

    return
  end subroutine


  
  subroutine cell_kinene( ekinh, temphh, velh )
    use constants, only: factem
    use cell_base, only: wmass
    implicit none
    real(kind=8), intent(out) :: ekinh, temphh(3,3)
    real(kind=8), intent(in)  :: velh(3,3)
    integer :: i,j
    ekinh = 0.0d0
    do j=1,3
      do i=1,3
        ekinh=ekinh+0.5*wmass*velh(i,j)*velh(i,j)
        temphh(i,j)=factem*wmass*velh(i,j)*velh(i,j)
      end do
    end do
    return
  end subroutine


  subroutine ions_vrescal( tcap, tempw, tempp, taup, tau0, taum, na, nsp, fion, iforce, &
                           pmass, delt )
    use constants, only: pi, factem
    implicit none
    logical, intent(in) :: tcap
    real(kind=8), intent(inout) :: taup(:,:,:)
    real(kind=8), intent(in) :: tau0(:,:,:), taum(:,:,:), fion(:,:,:)
    real(kind=8), intent(in) :: delt, pmass(:), tempw, tempp
    integer, intent(in) :: na(:), nsp
    integer, intent(in) :: iforce(:,:,:)

    real(kind=8) :: alfap, qr(3), alfar, gausp
    real(kind=8) :: dt2by2, ftmp
    real(kind=8) :: randy
    integer :: i, ia, is, nat

    dt2by2 = .5d0 * delt * delt
    gausp = delt * sqrt( tempw / factem )
    nat = SUM( na( 1:nsp ) )

    if(.not.tcap) then
      alfap=.5d0*sqrt(tempw/tempp)
      do is=1,nsp
        do ia=1,na(is)
          do i=1,3
            taup(i,ia,is) = tau0(i,ia,is) +                 &
     &                      alfap*(taup(i,ia,is)-taum(i,ia,is)) +      &
     &                      dt2by2/pmass(is)*fion(i,ia,is)*iforce(i,ia,is)
          end do
        end do
      end do
    else
      do i=1,3
        qr(i)=0.d0
        do is=1,nsp
          do ia=1,na(is)
            alfar=gausp/sqrt(pmass(is))*cos(2.d0*pi*randy())*sqrt(-2.d0*log(randy()))
            taup(i,ia,is)=alfar
            qr(i)=qr(i)+alfar
          end do
        end do
        qr(i)=qr(i)/nat
      end do
      do is=1,nsp
        do ia=1,na(is)
          do i=1,3
            alfar=taup(i,ia,is)-qr(i)
            taup(i,ia,is)=tau0(i,ia,is)+iforce(i,ia,is)*     &
     &                    (alfar+dt2by2/pmass(is)*fion(i,ia,is))
          end do
        end do
      end do
    end if
    return
  end subroutine


  subroutine ions_nosevel( vnhp, xnhp0, xnhpm, delt )
    implicit none
    real(kind=8), intent(inout) :: vnhp
    real(kind=8), intent(in) :: xnhp0, xnhpm, delt
    vnhp=2.*(xnhp0-xnhpm)/delt-vnhp
    return
  end subroutine

  subroutine elec_nosevel( vnhe,xnhe0,xnhem,delt, fccc )
    implicit none
    real(kind=8), intent(inout) :: vnhe
    real(kind=8), intent(out) :: fccc
    real(kind=8), intent(in) :: xnhe0, xnhem, delt
    vnhe=2.*(xnhe0-xnhem)/delt-vnhe
    fccc=1./(1.+0.5*delt*vnhe)
    return
  end subroutine

  subroutine cell_nosevel( vnhh, xnhh0, xnhhm, delt, velh, h, hold )
    implicit none
    real(kind=8), intent(inout) :: vnhh(3,3), velh(3,3)
    real(kind=8), intent(in) :: xnhh0(3,3), xnhhm(3,3), delt, h(3,3), hold(3,3)
    vnhh(:,:)=2.*(xnhh0(:,:)-xnhhm(:,:))/delt-vnhh(:,:)
    velh(:,:)=2.*(h(:,:)-hold(:,:))/delt-velh(:,:)
    return
  end subroutine


  subroutine ions_noseupd( xnhpp, xnhp0, xnhpm, delt, qnp, ekinpr, gkbt, vnhp )
    implicit none         
    real(kind=8), intent(out) :: xnhpp, vnhp
    real(kind=8), intent(in) :: xnhp0, xnhpm, delt, qnp, ekinpr, gkbt
    xnhpp=2.*xnhp0-xnhpm+2.*( delt**2 / qnp )*(ekinpr-gkbt/2.)
    vnhp =(xnhpp-xnhpm)/( 2.0d0 * delt )
    return
  end subroutine

 
  subroutine elec_noseupd( xnhep, xnhe0, xnhem, delt, qne, ekinc, ekincw, vnhe )
    implicit none
    real(kind=8), intent(out) :: xnhep, vnhe
    real(kind=8), intent(in) :: xnhe0, xnhem, delt, qne, ekinc, ekincw
    xnhep=2.*xnhe0-xnhem+2.*(delt**2/qne)*(ekinc-ekincw)
    vnhe =(xnhep-xnhem)/( 2.0d0 * delt )
    return
  end subroutine

 
  subroutine cell_noseupd( xnhhp, xnhh0, xnhhm, delt, qnh, temphh, temph, vnhh )
    use constants, only: factem
    implicit none
    real(kind=8), intent(out) :: xnhhp(3,3), vnhh(3,3)
    real(kind=8), intent(in) :: xnhh0(3,3), xnhhm(3,3), delt, qnh, temphh(3,3), temph
    integer :: i, j
    do j=1,3
      do i=1,3
        xnhhp(i,j)=2.*xnhh0(i,j)-xnhhm(i,j)+ (delt**2/qnh)/factem*(temphh(i,j)-temph)
        vnhh(i,j) =(xnhhp(i,j)-xnhhm(i,j))/( 2.0d0 * delt )
      end do
    end do
    return
  end subroutine


end module cpr_subroutines
