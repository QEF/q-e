!
!---------------------------------------------------------------
subroutine lderivps
  !---------------------------------------------------------------
  !
  !  numerical integration of the radial schroedinger equation 
  !  computing logarithmic derivatives for pseudo-potentials
  !  multiple nonlocal projectors are allowed
  !
  use ld1inc
  implicit none

  integer  ::       &
       lam,   &      ! the angular momentum
       ikrld, &      ! index of matching radius
       nc,    &      ! counter on logarithmic derivatives
       nbf,   &      ! number of b functions
       n,ie          ! generic counters

  real(kind=dp) ::  &
       ze2,     &    ! the nuclear charge in Ry units
       jam,     &    ! the total angular momentum
       e,       &    ! the eigenvalue
       lamsq,            & ! combined angular momentum
       b(0:3),c(4),      & ! used for starting guess of the solution 
       b0e, rr1,rr2,     & ! auxiliary
       xl1, x4l6, ddx12, &
       x6l12, x8l20

  real(kind=dp),allocatable :: &
       dlchis(:,:), &  ! the logarithmic derivatives
       vaux(:),     &  !
       aux(:),      &  ! the square of the wavefunction
       al(:)           ! the known part of the differential equation

  real(kind=dp), external ::           &
       compute_log, &
       int_0_inf_dr

  integer :: &
       ib,jb,iib,jjb, &  ! counters on beta functions
       nst,nstop,     &  ! auxiliary for integrals
       ios,           &  ! used for I/O control
       is, ind           ! counters on index


  if (nld == 0 .or. file_logderps == ' ') return
  if (nld > nwfsx) call errore('lderivps','nld is too large',1)

  allocate( al(mesh), aux(mesh), vaux(mesh) )

  ze2=0.0_dp

  do n=1,mesh
     if (r(n) > rlderiv) go to 10
  enddo
  call errore('lderivps','wrong rlderiv?',1)
10 ikrld = n-1
  write(6,'(5x,''Computing logarithmic derivative in'',f10.5)') &
       (r(ikrld)+r(ikrld+1))*0.5_dp
  npte= (emaxld-eminld)/deld + 1
  allocate ( dlchis(npte,nld) )
  do is=1,nspin
     do nc=1,nld
        if (rel < 2) then
           lam=nc-1
           jam=0.0_dp
        else
           lam=nc/2
           if (mod(nc,2)==0) jam=lam-0.5_dp
           if (mod(nc,2)==1) jam=lam+0.5_dp
        endif
        xl1=lam+1
        x4l6=4*lam+6
        x6l12=6*lam+12
        x8l20=8*lam+20
        ddx12=dx*dx/12.0_dp
        nst=(lam+1)**2  
        nbf=nbeta
        if (pseudotype == 1) then
           if (rel == 2) then
              if (abs(jam-lam+0.5_dp) < 1.e-2_dp .or. lam == 0 ) then
                 ind=1
              else
                 ind=2
              endif
              do n=1,mesh
                 vpstot(n,is)=vpstot(n,is)+vnlo(n,lam,ind)
                 vaux(n)=vnlo(n,lam,ind)
              enddo
           else
              do n=1,mesh
                 vpstot(n,is)=vpstot(n,is)+vnl(n,lam)
                 vaux(n)=vnl(n,lam)
              enddo
           endif
           nbf=0.0_dp
        endif

        do n=1,4
           al(n)=vpstot(n,is)-ze2/r(n)
        enddo
        call series(al,r,r2,b)

        do ie=1,npte
           e=eminld+deld*(ie-1.0_dp)
           lamsq=(lam+0.5_dp)**2
           !
           !     b) find the value of solution s in the first two points
           !
           b0e=b(0)-e
           c(1)=0.5_dp*ze2/xl1
           c(2)=(c(1)*ze2+b0e)/x4l6
           c(3)=(c(2)*ze2+c(1)*b0e+b(1))/x6l12
           c(4)=(c(3)*ze2+c(2)*b0e+c(1)*b(1)+b(2))/x8l20
           rr1=(1.0_dp+r(1)*(c(1)+r(1)* &
                (c(2)+r(1)*(c(3)+r(1)*c(4)))))*r(1)**(lam+1)
           rr2=(1.0_dp+r(2)*(c(1)+r(2)* &
                (c(2)+r(2)*(c(3)+r(2)*c(4)))))*r(2)**(lam+1)
           aux(1)=rr1/sqr(1)
           aux(2)=rr2/sqr(2)

           do n=1,mesh
              al(n)=( (vpstot(n,is)-e)*r2(n) + lamsq )*ddx12
              al(n)=1.0_dp-al(n)
           enddo

           call integrate_outward (lam,jam,e,mesh,ndm,dx,r,r2,sqr,al, &
                b,aux,betas,ddd,qq,nbf,nwfsx,lls,jjs,ikrld+5)

           !
           !    compute the logarithmic derivative and save in dlchi
           !            
           do n=-3,3
              aux(ikrld+n)= aux(ikrld+n)*sqr(ikrld+n)
           enddo

           dlchis(ie,nc)=compute_log(aux(ikrld-3),r(ikrld),dx)
        enddo
        if (pseudotype == 1) then
           do n=1,mesh
              vpstot(n,is)=vpstot(n,is)-vaux(n)
           enddo
        endif
     enddo

     if (is == 2) then
        open(unit=25, file=trim(file_logder)//'01', status='unknown', &
             iostat=ios, err=300 )
     else
        open(unit=25, file=trim(file_logder), status='unknown', &
             iostat=ios, err=300 )
     end if
300  call errore('lderivps','opening file '//trim(file_logder), abs(ios))

     do ie=1,npte
        e= eminld+deld*(ie-1)
        write(25,'(10f14.6)') e, (dlchis(ie,nc),nc=1,nld)
     enddo
     close(unit=25)
  enddo

  deallocate(dlchis)
  deallocate(vaux, aux, al)

  return
end subroutine lderivps
