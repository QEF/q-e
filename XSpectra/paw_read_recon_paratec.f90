!
! Copyright (C) 2003-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine read_recon_paratec(filerec)

  !
  ! Read all-electron and pseudo atomic wavefunctions 
  !  needed for PAW reconstruction
  !

  USE ions_base,          ONLY : ntyp => nsp
  use atom, only: mesh, msh, r, rab
  use kinds, only: DP
  use parameters, only : ntypx
  USE io_global,  ONLY : stdout
  use splinelib
  
  USE paw_gipaw, ONLY : at_wfc, paw_recon, psphi, aephi, paw_wfc_init
  
  implicit none

  character (len=256), intent ( in ) :: filerec(ntypx)
  
  character (len=256) :: readline
  logical :: file_exists
  integer :: iostatus, nlines ( ntypx )
  integer :: l,j,i,jtyp,kkphi,nbetam
  real(dp) :: d1
  
  real(dp), allocatable :: xdata(:), tab(:), tab_d2y(:),ae_r(:,:),ps_r(:,:)
  type(at_wfc),pointer :: aephi2(:,:), psphi2(:,:) ! Atom

  do jtyp=1,ntyp
     inquire ( file = filerec(jtyp), exist = file_exists )
     if ( .not. file_exists ) then
        write(6,*) "reconstruction file does not exist: ",TRIM(filerec(jtyp))
        call errore("reconstruction file does not exist",TRIM(filerec(jtyp)),0)
        stop
     end if
     open(14,file=filerec(jtyp))
     paw_recon(jtyp)%paw_nbeta = 0
     nlines ( jtyp ) = 0
     do
        READ ( UNIT = 14, FMT = '( 256A )', IOSTAT = iostatus ) readline
        IF ( iostatus /= 0 ) THEN
           EXIT
        END IF
        IF ( INDEX ( readline, "#core wavefunctions" ) /= 0 ) THEN
           EXIT
        END IF
        IF ( INDEX ( readline, "#" ) /= 0 ) THEN
           nlines ( jtyp ) = 0
        END IF
        IF ( INDEX ( readline, "&" ) /= 0 ) THEN
           paw_recon(jtyp)%paw_nbeta = paw_recon(jtyp)%paw_nbeta + 1
        END IF
        IF ( INDEX ( readline, "&" ) == 0 .AND. &
             INDEX ( readline, "#" ) == 0 ) THEN
           nlines ( jtyp ) = nlines ( jtyp ) + 1
        END IF
     end do
     close(14)
  enddo
 
  nbetam=maxval(paw_recon(:ntyp)%paw_nbeta)
  allocate( psphi2(ntyp,nbetam) )
  allocate( aephi2(ntyp,nbetam) )

  call paw_wfc_init(psphi2)
  call paw_wfc_init(aephi2)

  allocate(ae_r(ntyp,maxval(nlines(1:ntyp))))
  allocate(ps_r(ntyp,maxval(nlines(1:ntyp))))

  recphi_read: do jtyp=1,ntyp
     open(14,file=filerec(jtyp))
     rewind(unit=14)
     recphi_loop: do i=1,paw_recon(jtyp)%paw_nbeta
        aephi2(jtyp,i)%label%nt=jtyp
        aephi2(jtyp,i)%label%n=i
        kkphi = nlines ( jtyp )
        aephi2(jtyp,i)%kkpsi=kkphi
        allocate(aephi2(jtyp,i)%psi(kkphi))
        allocate(psphi2(jtyp,i)%psi(kkphi))
        read(14,*) readline
        read(14, FMT = '( 256A )' ) readline
        IF ( readline(1:6) == "#local" ) THEN
           read(14, FMT = '( 256A )' ) readline
        END IF
        IF ( readline(1:3) /= "#l=" ) THEN
           WRITE ( UNIT = 6, FMT = '( 2A, ":" )' ) &
                "Wrong control string in file ", TRIM ( filerec(jtyp) )
           WRITE ( UNIT = 6, FMT = '( A )' ) TRIM ( readline )
           CALL errore ( "read_recon_paratec", "wrong control string", 0 )
           stop
        END IF
        read(readline(4:4), FMT = * ) aephi2(jtyp,i)%label%l
        read(readline(12:), FMT = * ) aephi2(jtyp,i)%label%rc
        DO j = 1, kkphi
           read(14,*) ae_r(jtyp,j), aephi2(jtyp,i)%psi(j), &
                psphi2(jtyp,i)%psi(j)
        END DO
        psphi2(jtyp,i)%label%nt=jtyp
        psphi2(jtyp,i)%label%n=i
        psphi2(jtyp,i)%label%l=aephi2(jtyp,i)%label%l
        psphi2(jtyp,i)%label%rc=aephi2(jtyp,i)%label%rc
        psphi2(jtyp,i)%kkpsi=kkphi
        ps_r(jtyp,:) = ae_r(jtyp,:)
     end do recphi_loop
     close(14)
  end do recphi_read
  
!  do j = 1, psphi2(1,1)%kkpsi
!     write(95,*) psphi2(1,1)%r(j), aephi2(1,1)%psi(j)
!     write(96,*) psphi2(1,1)%r(j), psphi2(1,1)%psi(j)
!  end do
  
  allocate( psphi(ntyp,nbetam) )
  allocate( aephi(ntyp,nbetam) )

  call paw_wfc_init(psphi)
  call paw_wfc_init(aephi)

  do jtyp=1, ntyp
     do i = 1, paw_recon(jtyp)%paw_nbeta
        kkphi = aephi2(jtyp,i)%kkpsi
        allocate( xdata(kkphi), tab(kkphi), tab_d2y(kkphi) )
        xdata(:) = ae_r(jtyp,:)
        tab(:) = aephi2(jtyp,i)%psi(:)
        
        ! initialize spline interpolation
        d1 = (tab(2) - tab(1)) / (xdata(2) - xdata(1))
        call spline(xdata, tab, 0.d0, d1, tab_d2y)
!        call spline(xdata, tab, tab_d2y)
        
        ! use interpolation
        allocate ( aephi(jtyp,i)%psi(msh(jtyp)) )
        aephi(jtyp,i)%label%nt = jtyp
        aephi(jtyp,i)%label%n = i
        aephi(jtyp,i)%label%l = aephi2(jtyp,i)%label%l
        aephi(jtyp,i)%label%rc = aephi2(jtyp,i)%label%rc
        aephi(jtyp,i)%kkpsi = msh(jtyp)
        do j = 1, msh(jtyp)
           aephi(jtyp,i)%psi(j) = splint(xdata, tab, tab_d2y, r(j,jtyp))
        end do
        
        allocate ( psphi(jtyp,i)%psi(msh(jtyp)) )
        xdata(:) = ps_r(jtyp,:)
        tab(:) = psphi2(jtyp,i)%psi(:)
        
        ! initialize spline interpolation
        d1 = (tab(2) - tab(1)) / (xdata(2) - xdata(1))
        call spline(xdata, tab, 0.d0, d1, tab_d2y)
!        call spline(xdata, tab, tab_d2y)
        
        ! use interpolation
        allocate ( psphi(jtyp,i)%psi(msh(jtyp)) )
        psphi(jtyp,i)%label%nt = jtyp
        psphi(jtyp,i)%label%n = i
        psphi(jtyp,i)%label%l = psphi2(jtyp,i)%label%l
        psphi(jtyp,i)%label%rc = psphi2(jtyp,i)%label%rc
        psphi(jtyp,i)%kkpsi = msh(jtyp)
        do j = 1, msh(jtyp)
           psphi(jtyp,i)%psi(j) = splint(xdata, tab, tab_d2y, r(j,jtyp))
        end do
        
        deallocate( xdata, tab, tab_d2y)
     enddo
  enddo
  
  deallocate( psphi2 )
  deallocate( aephi2 )
  
!  do j = 1, aephi(1,1)%kkpsi
!     write(85,*) aephi(1,1)%r(j), aephi(1,1)%psi(j)
!     write(86,*) psphi(1,1)%r(j), psphi(1,1)%psi(j)
!  end do

  
end subroutine read_recon_paratec


