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
    real(kind=8) :: var(:,:)
    integer :: na(:), nsp
    integer, optional :: iunit
    character(len=*), optional :: head
    integer :: i, ia, is, iu, isa
    if( present( iunit ) ) then
      iu = iunit
    else
      iu = stdout
    end if
    if( present( head ) ) then 
      WRITE( iu,*) head
    end if
    isa = 0
    DO is = 1, nsp
      DO ia = 1, na(is)
        isa = isa + 1
        WRITE( iu,'(3f14.8)') ( var(i,isa), i=1, 3 )
      END DO
    END DO
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
    real(kind=8), intent(in) :: tausm(:,:), pmass(:), fion(:,:)
    integer, intent(in) :: iforce(:,:)
    real(kind=8), intent(in) :: ainv(3,3), delt
    real(kind=8), intent(out) :: taus(:,:)
    integer, intent(in) :: na(:), nsp
    integer :: is, ia, i, isa
    real(kind=8) :: dt2by2, fac, fions(3)

    dt2by2 = .5d0 * delt * delt

    isa = 0
    do is=1,nsp
      fac = dt2by2/pmass(is)
      do ia=1,na(is)
        isa = isa + 1
        do i=1,3
          fions( i ) = fion(1,isa)*ainv(i,1) + fion(2,isa)*ainv(i,2) + fion(3,isa)*ainv(i,3)
        end do
        do i=1,3
          taus(i,isa) = tausm(i,isa) + iforce(i,isa) * fac * fions( i )
        end do
      end do
    end do
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
    real(kind=8), intent(in) :: taus(:,:), tausm(:,:), pmass(:), fion(:,:)
    integer, intent(in) :: iforce(:,:)
    real(kind=8), intent(in) :: ainv(3,3), delt
    real(kind=8), intent(out) :: tausp(:,:)
    integer, intent(in) :: na(:), nsp
    real(kind=8), intent(in) :: fricp, hgamma(3,3), vels(:,:)
    logical, intent(in) :: tsdp, tnosep
    real(kind=8), intent(inout) :: fionm(:,:)
    real(kind=8), intent(in) :: vnhp
    real(kind=8), intent(out) :: velsp(:,:)
    real(kind=8), intent(in) :: velsm(:,:)
    integer :: is, ia, i, isa
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
            isa = 0
            do is=1,nsp
               do ia=1,na(is)
                  isa = isa + 1
                  do i=1,3
                     tausp(i,isa) = taus(i,isa) +                   &
     &                    iforce(i,isa)*dt2by2/pmass(is)*             &
     &        (ainv(i,1)*fion(1,isa)+ainv(i,2)*fion(2,isa)+         &
     &         ainv(i,3)*fion(3,isa) ) -                              &
     &                    pmass(is)*(hgamma(i,1)*vels(1,isa)+         &
     &         hgamma(i,2)*vels(2,isa)+hgamma(i,3)*vels(3,isa))
                  end do
               end do
            end do
         else if (tnosep) then
            isa = 0
            do is=1,nsp
               do ia=1,na(is)
                  isa = isa + 1
                  do i=1,3
                     fionm(i,isa) = (ainv(i,1)*fion(1,isa)          &
     &                                +ainv(i,2)*fion(2,isa)          &
     &                                +ainv(i,3)*fion(3,isa))         &
     &                              - vnhp*vels(i,isa)*pmass(is)      &
     &                    - pmass(is)*(hgamma(i,1)*vels(1,isa)        &
     &                                +hgamma(i,2)*vels(2,isa)        &
     &                                +hgamma(i,3)*vels(3,isa))
                     tausp(i,isa)=-tausm(i,isa)+2.*taus(i,isa)+   &
     &                   iforce(i,isa)*dt2*fionm(i,isa)/pmass(is)
                     velsp(i,isa) = velsm(i,isa) +                  &
     &                    twodel*fionm(i,isa)/pmass(is)
                  end do
               end do
            end do
         else
            isa = 0
            do is=1,nsp
               do ia=1,na(is)
                  isa = isa + 1
                  do i=1,3
                     tausp(i,isa) = verl1*taus(i,isa)               &
     &                    + verl2*tausm(i,isa)                        &
     &        + verl3/pmass(is)*iforce(i,isa) * (ainv(i,1)*fion(1,isa)&
     &        + ainv(i,2)*fion(2,isa) + ainv(i,3)*fion(3,isa))      &
     &        - verl3*iforce(i,isa) * (hgamma(i,1)*vels(1,isa)      &
     &        + hgamma(i,2)*vels(2,isa) + hgamma(i,3)*vels(3,isa))
                     velsp(i,isa)=velsm(i,isa)                      &
     &        - 4.*fricp*vels(i,isa)                                  &
     &        + twodel/pmass(is)*iforce(i,isa)*(ainv(i,1)*fion(1,isa) &
     &        + ainv(i,2)*fion(2,isa) + ainv(i,3)*fion(3,isa))      &
     &        - twodel*iforce(i,isa) * (hgamma(i,1)*vels(1,isa)     &
     &        + hgamma(i,2)*vels(2,isa) + hgamma(i,3)*vels(3,isa))
                  end do
               end do
            end do
         endif
    return
  end subroutine


  subroutine ions_cofmsub( tausp, na, nsp, cdm, cdm0 )
    implicit none
    real( kind=8 ), intent(inout) :: tausp( :, : )
    integer, intent(in) :: na(:), nsp
    real( kind=8 ), intent(in) :: cdm( : ), cdm0( : )
    integer :: i, ia, is, isa
    isa = 0
    do is=1,nsp
      do ia=1,na(is)
        isa = isa + 1
        do i=1,3
          tausp(i,isa)=tausp(i,isa)+cdm0(i)-cdm(i)
        enddo
      enddo
    enddo
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

  subroutine elec_fakekine2( ekincm, ema0bg, emass, c0, cm, ngw, n, delt )
    use mp, only: mp_sum
    use reciprocal_vectors, only: gstart
    use wave_base, only: wave_speed2
    real(kind=8), intent(out) :: ekincm
    real(kind=8), intent(in)  :: ema0bg(:), delt, emass
    complex(kind=8), intent(in)  :: c0(:,:), cm(:,:)
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
               wave_speed2( c0(:,i), cm(:,i), emainv, ftmp )
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


  subroutine print_lambda( lambda, n, nshow, ccc, iunit )
    use io_global, only: stdout, ionode
    real(kind=8), intent(in) :: lambda(:,:), ccc
    integer, intent(in) :: n, nshow
    integer, intent(in), optional :: iunit
    integer :: nnn, j, un, i
    if( present( iunit ) ) then
      un = iunit
    else
      un = stdout
    end if
    nnn=min(n,nshow)
    if( ionode ) then
       WRITE( un,*)
       WRITE( un,3370) '    lambda   n = ', n
       IF( nnn < n ) WRITE( un,3370) '    print only first ', nnn
       do i=1,nnn
          WRITE( un,3380) (lambda(i,j)*ccc,j=1,nnn)
       end do
    end if
3370     format(26x,a,i4)
3380     format(9f8.4)
    return
  end subroutine

   subroutine add_thermal_stress( stress, pmass, omega, h, vels, nsp, na )
     real(kind=8) :: stress(3,3)
     real(kind=8), intent(in) :: pmass(:), omega, h(3,3), vels(:,:)
     integer, intent(in) :: nsp, na(:)
     integer :: i, j, is, ia, isa
     isa = 0
     do is=1,nsp
       do ia=1,na(is)
       isa = isa + 1
         do i=1,3
           do j=1,3
             stress(i,j)=stress(i,j)+pmass(is)/omega*           &
      &        ((h(i,1)*vels(1,isa)+h(i,2)*vels(2,isa)+    &
      &          h(i,3)*vels(3,isa))*(h(j,1)*vels(1,isa)+  &
      &          h(j,2)*vels(2,isa)+h(j,3)*vels(3,isa)))
           enddo
         enddo
       enddo
     enddo
     return
   end subroutine


end module cpr_subroutines
