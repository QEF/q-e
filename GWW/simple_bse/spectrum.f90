subroutine spectrum(data_input,x,bd,energy)
  USE input_simple_exc
  USE simple_objects
  USE derived_objects
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout, ionode
  USE constants, ONLY : rytoev
  USE mp, ONLY : mp_barrier, mp_sum
  USE mp_world,             ONLY : world_comm,mpime
  use io_files, only : prefix, tmp_dir
  USE constants, ONLY : rytoev
  
  implicit none
  
  TYPE(input_options):: data_input
  TYPE(exc):: x(data_input%nvec)
  TYPE(bands):: bd
  COMPLEX(KIND=DP) :: energy(data_input%nvec)
  COMPLEX(KIND=DP), ALLOCATABLE :: diel(:),D(:,:),N(:),Etempc(:,:,:),Etempv(:,:,:),eemat(:,:,:,:),mate(:,:,:,:),coeff(:,:), cabs(:)
  REAL(KIND=DP)::step,omega,cost,omin,omax, norm
  INTEGER :: i,j,v,c,k,S,a,kk,io,iun
  INTEGER, EXTERNAL :: find_free_unit
  allocate(diel(data_input%spectrum_points))!costante dielettrica(omega)
  allocate(D(data_input%nvec,data_input%spectrum_points))!denominatore(S,omega)
  allocate(N(data_input%nvec))!numeratore(S)
  allocate(eemat(bd%numv,bd%numc,bd%nk_loc,3))!E_{iv}^k^*  E_{jc}^k <e_i|p_a|e_j>
  allocate(coeff(3,data_input%nvec))!!elemento con coeff fissati dir,S 
  allocate(cabs(data_input%spectrum_points))
  
  !comodità
  !energy(1)  =2.49911d-2
  !energy(2)  =2.57741d-2
  !energy(3)  =2.59408d-2
  !energy(4)  =2.64996d-2
  !energy(5)  =2.66287d-2
  !energy(6)  =2.68860d-2
  !energy(7)  =3.38354d-2
  !energy(8)  =3.52200d-2
  !energy(9)  =3.60030d-2
  !energy(10) =3.62000d-2
  !energy(11) =3.73634d-2
  !energy(12) =3.88258d-2
  !energy(13) =3.74831d-2
  !energy(14) =3.93831d-2
  !energy(15) =3.98621d-2
  !energy(16) =4.02376d-2
  write(stdout,*)'Energies eV'
  do j=1,data_input%nvec
     write(stdout,*) energy(j)*rytoev
  end do
  write(stdout,*) 'Calculating spectrum'

  call build_eemat(data_input,bd,eemat)
  
  !costruzione coefficienti*elemento

  do a=1,3
     do S=1,data_input%nvec
        coeff(a,S)=(0.d0,0.d0)
        do k=1,bd%nk_loc
           do v=1,bd%numv
              do c=1,bd%numc
                 coeff(a,S) = coeff(a,S) + x(S)%avc(v,c,k)*eemat(v,c,k,a)
 !                norm = norm + 1
 !                write(*,*)coeff(a,S)
              end do
           end do
        end do
     end do
  end do
  write(*,*)'coeff'
  
  call mp_sum(coeff,world_comm)
  do S=1,data_input%nvec
     N(S)=(0.d0,0.d0)
     do a=1,3
        N(S) = N(S) + conjg(coeff(a,S))*coeff(a,S)
!        write(*,*)N(S)
     end do
     N(S)=N(S)/3
     write(*,*)N(S)
  end do
  write(*,*)'numerator'
  !denominator
  omin = data_input%omega_min/rytoev!Ry
  omax = data_input%omega_max/rytoev!Ry
 ! write(*,*)'omega_max',omax,data_input%omega_max
  step =  (omax-omin)/(data_input%spectrum_points-1)
 ! write(*,*)'step',step
  do io=1,data_input%spectrum_points
     omega = omin + (io-1)*step!Ry
     do S=1,data_input%nvec
       ! write(*,*)energy(S),omega
        D(S,io) = energy(S)*( energy(S)*energy(S) - ( cmplx(omega,data_input%eta))&
             &*(cmplx(omega,data_input%eta)))/data_input%eta
!        write(*,*)omega,energy(S),D(S,io)
     end do
  end do
  write(*,*)'denominator'
  !constant
!
  cost = 16/3.14159/10.26!general case for V?
  !epsilon
  do io=1,data_input%spectrum_points
     diel(io) = (0.d0,0.d0)
     do S = 1,data_input%nvec
        diel(io) = diel(io) + N(S)/D(S,io)
     end do
     diel(io) = 1 + diel(io)*cost!all quantities in Rydberg atomic units
        !  write(*,*)io, diel(io)
  end do
  write(*,*)'epsilon'
  !saving spectrum
  if(ionode) then
     iun=find_free_unit()
     open( unit= iun, file='spectrum.txt',status='unknown')
     do io=1,data_input%spectrum_points
        omega = omin + (io-1)*step
        cabs(io)=omega*aimag(diel(io))/(sqrt(4*3.14159*diel(io)))/(137)
        !        write(iun,*) omega*rytoev, cabs(io)
        write(iun,*)omega*rytoev, aimag(diel(io))
     end do
        close(iun)
  end if
  write(*,*)'saved epsilon'
  deallocate(diel)
  deallocate(N)
  deallocate(D)
  deallocate(eemat)
  deallocate(coeff)
  deallocate(cabs)
  return
  
end subroutine spectrum
