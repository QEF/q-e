! FOR GWW
!
! Author: P. Umari
! Modified by G. Stenuit
!
!-------------------------
subroutine go_wannier( iun_wannier, tresh, maxiter,nbndv, itask, l_on_ene, ene_loc, lambda)
!-------------------------
! this routine read the wfcs from iun_wannier then
! transfrom to real wfcs, transform to wannier functions
! with treshold tresh and maxiter
! using Gygi scheme
! and writes wfcs on iun_wannier and writes wannier_centers file
! it's possible to localized two different subspaces

! #ifdef __GWW

  USE kinds,    ONLY : DP
  USE us
  USE wvfct,    ONLY : igk, g2kin, npwx, npw, nbndx
  USE control_flags, ONLY: gamma_only
  USE wavefunctions_module, ONLY : evc
  USE gvect
  USE basis
  USE klist
  USE constants, ONLY : e2, pi, tpi, fpi
  USE io_files, ONLY: nwordwfc
  USE io_global, ONLY : stdout
  USE cell_base, ONLY: at, alat, tpiba, omega, tpiba2
  USE wannier_gw
  USE mp, ONLY : mp_sum, mp_bcast, mp_barrier
  USE mp_global,  ONLY : mpime, nproc

  implicit none

  INTEGER, INTENT(in) :: iun_wannier !units for reading wfc
  REAL(kind=DP), INTENT(in) :: tresh! treshold on wannier wfcs'spread
  INTEGER, INTENT(in) :: maxiter!max number of iterations
  INTEGER, INTENT(in) :: nbndv !number of first bands wich are localized separately
  INTEGER,INTENT(in)  :: itask! if == 1 calculate {C'} subspace
  LOGICAL, INTENT(in) :: l_on_ene!if true localizes also on energy
  REAL(kind=DP), INTENT(inout) :: ene_loc(nbnd_normal)!in input the KS energies on output the expectation values of
                                              !the energy operator for the localized states
  REAL(kind=DP), INTENT(in) :: lambda!paramter of energy operator


  !  --- Internal definitions ---
  INTEGER :: i,j,k,l,it,iw, ii
  COMPLEX(kind=DP), ALLOCATABLE  :: matgp(:,:,:)! for matrix <Psi_i|{exp(iG{x,y,z})|Psi_j>
  REAL(kind=DP), ALLOCATABLE  :: matsincos(:,:,:)  !respective sin and cos matrices
  REAL(kind=DP), ALLOCATABLE  ::rot_u(:,:) ! unitarian matrix which transorm to wannier wfcs
  REAL(kind=DP) :: w(6) ! weights for orthorombic cells
  COMPLEX(kind=DP), ALLOCATABLE :: crot_u(:,:),crot_u_tmp(:,:)
  INTEGER nbnd_start, nbnd_end!for conduction states manifolds

  REAL(kind=DP) :: theta, theta4, d2, aa(2), cc, ss, tempi, tempj
  REAL(kind=DP) :: omg0,omg1!omega for testi convergence
  COMPLEX(kind=DP) :: sca

  REAL(kind=DP), ALLOCATABLE :: tmp_mat(:,:)
  INTEGER :: nbnd_second

  INTEGER :: n_oper

  INTEGER :: na,nb
  REAL(kind=DP) :: delta
  REAL(kind=DP), ALLOCATABLE :: vtmpi(:),vtmpj(:)
  REAL(Kind=DP), ALLOCATABLE :: ene_tmp(:)
  INTEGER :: nbnd_par,nbnd_start_me,nbnd_end_me,nbnd_start_proc,nbnd_end_proc,ip
  REAL(kind=DP), ALLOCATABLE :: matsincos_par(:,:,:)
  INTEGER :: n_oper_par,oper_start,oper_end
  REAL(kind=DP) :: sca_mat(6)
  !!! new
  REAL(kind=DP) :: diff, error_relative
  REAL(kind=DP) :: accuracy_sincos=1.d-14 !! and not 1.d-14 !! if theta is p=14 => cos(theta) is p=7

!  ALLOCATE( matgp(nbnd_normal,nbnd_normal,3))
   if(.not.l_on_ene) then
     ALLOCATE( matsincos(nbnd_normal,nbnd_normal,6))
  else
     ALLOCATE( matsincos(nbnd_normal,nbnd_normal,7))
  endif





!transfrom to real wfcs:

  write(stdout,*) 'Transform to real wfcs'
  call flush_unit(stdout)

  if(nbndv<1 .OR. nbndv > nbnd_normal) call errore('go_wannier','nbndv: illegal value',1)
  !if(mod(nbndv,2) /= 0 ) call errore('go_wannier','nbndv, odd',1)

  if(.not.gamma_only) then
    call real_wfc(u_trans, 1, iun_wannier,nbndv)
  else
  !writes KS wavefunctions on file
    u_trans(:,:)=(0.d0,0.d0)!it could also be defined as real
    do iw=1,nbnd_normal
      u_trans(iw,iw)=(1.d0,0.d0)
    enddo
  endif

! set matg p,m matrices
  if(gamma_only) then
    !call matrix_wannier_gamma_big(matgp,1,nset,itask)
     call matrix_wannier_gamma_big(matsincos,1,nset,itask)
  else
     write(stdout,*)'ONLY GAMMA POINT IMPLEMENTED'
     stop
  endif

  write(stdout,*) 'Out of matrix_wannier_gamma_big'
  call  flush_unit(stdout)


! set weights

  do i=1,3
     w(i)=((at(i,i)*alat)/pi)
     w(i+3)=((at(i,i)*alat)/pi)
  enddo

 ! calculates matsincos

  do i=1,3
     matsincos(1:nbnd_normal,1:nbnd_normal,i)=w(i)*matsincos(1:nbnd_normal,1:nbnd_normal,i)
     matsincos(1:nbnd_normal,1:nbnd_normal,i+3)=w(i)*matsincos(1:nbnd_normal,1:nbnd_normal,i+3)
   !  matsincos(1:nbnd_normal,1:nbnd_normal,i)=w(i)*real(matgp(1:nbnd_normal,1:nbnd_normal,i))
   !  matsincos(1:nbnd_normal,1:nbnd_normal,i+3)=w(i)*aimag(matgp(1:nbnd_normal,1:nbnd_normal,i))
  enddo

 ! deallocate(matgp)


  if(l_on_ene) then
     matsincos(:,:,1:6)=(1.d0-lambda)*matsincos(:,:,1:6)
     matsincos(:,:,7) = 0.d0
     do i=1,nbnd_normal
        matsincos(i,i,7)=lambda*ene_loc(i)
     enddo
  endif

  if(l_on_ene) then
     n_oper=7
  else
     n_oper=6
  endif

!----------------Young'su stuff-----




! Scale thr according to # of Wannier orbitals and cell length

!  set initial rotation matrix

  ALLOCATE (rot_u(nbnd_normal,nbnd_normal))
  rot_u(:,:)=0.d0
  do i=1,nbnd_normal
     rot_u(i,i)=1.d0
  end do

! now  valence subspace


! calculate omega
  omg0  = 0.d0
  do k=1,6
     do i=1,nbndv
        omg0=omg0 + matsincos(i,i,k)*matsincos(i,i,k)
     enddo
  enddo

  write(stdout ,*) 'LOCALIZING WANNIER FUNCTIONS:',omg0
  call flush_unit(stdout)


! Start Iteration =====================================================


  do it=1,maxiter
     do i=1,nbndv
        do j=i+1,nbndv

! Construct aa
           aa(:)=0.d0
           do k=1,n_oper
              !!! new
              diff=matsincos(i,i,k)-matsincos(j,j,k)
              aa(1)=aa(1)+matsincos(i,j,k)*diff
              aa(2)=aa(2)+matsincos(i,j,k)*matsincos(i,j,k)
              aa(2)=aa(2)-0.25d0*diff*diff
              !!!aa(1)=aa(1)+matsincos(i,j,k)*(matsincos(i,i,k)-matsincos(j,j,k))
              !!!aa(2)=aa(2)+matsincos(i,j,k)*matsincos(i,j,k)-0.25d0*(matsincos(i,i,k)-matsincos(j,j,k))*(matsincos(i,i,k)-matsinco
           end do
           !!! for debugging
           !write(stdout ,*) 'aa(1)=', aa(1) , 'and aa(2)=', aa(2)
           !write(stdout ,*) '|aa(1)|=', abs(aa(1)) ,'and |aa(2)|=', abs(aa(2))
           !call flush_unit(stdout)

           if(abs(aa(2)).gt. accuracy_sincos ) then
              theta4=-aa(1)/aa(2)
              theta=0.25*datan(theta4)
           elseif (abs(aa(1)).lt. accuracy_sincos ) then
              theta=0.d0
              aa(2)=0.d0
           else
              theta=pi/4.d0
           endif
           !!! for debugging
           !write(stdout ,*) 'theta=', theta

           d2=aa(1)*dsin(4.*theta)-aa(2)*dcos(4.*theta)
           !!! for debugging
           !write(stdout ,*) 'd2=', d2

           if(d2.le.0.d0) theta=theta+pi/4.d0
!
           cc=dcos(theta)
           ss=dsin(theta)
           !!! for debugging
           !write(stdout ,*) 'cc=', cc, 'and ss=', ss

! update overlap matrices
           do l=1,n_oper
              ! AR
              do k=1,nbndv
                 tempi=matsincos(k,i,l)*cc+matsincos(k,j,l)*ss
                 tempj=-matsincos(k,i,l)*ss+matsincos(k,j,l)*cc
                 matsincos(k,i,l)=tempi
                 matsincos(k,j,l)=tempj
              end do
! R^+ A R
              do k=1,nbndv
                 tempi=cc*matsincos(i,k,l)+ss*matsincos(j,k,l)
                 tempj=-ss*matsincos(i,k,l)+cc*matsincos(j,k,l)
                 matsincos(i,k,l)=tempi
                 matsincos(j,k,l)=tempj
              end do
           end do

! update U : U=UR
           do k=1,nbndv
              tempi=rot_u(k,i)*cc+rot_u(k,j)*ss
              tempj=-rot_u(k,i)*ss+rot_u(k,j)*cc
              rot_u(k,i)=tempi
              rot_u(k,j)=tempj
           end do
!
        end do
     end do

     omg1=0.d0
     do k=1,n_oper
        do i=1,nbndv
           omg1=omg1+matsincos(i,i,k)*matsincos(i,i,k)
        end do
     end do


     write(stdout,*) 'Spread', omg1,omg0,nbnd_normal,nbndv
     call flush_unit(stdout)

     error_relative=abs(omg1-omg0)/omg0

     write(stdout,*) 'error_relative=', error_relative, ' while tresh=', tresh
     call flush_unit(stdout)

     if(error_relative .lt. tresh ) EXIT
     omg0=omg1

  end do

!centers of wanniers

  do i=1,nbndv
     do k=1,3
        wannier_centers(k,i)=aimag(log(matsincos(i,i,k)+(0.d0,1.d0)*matsincos(i,i,k+3)))
        wannier_centers(k,i)=wannier_centers(k,i)*at(k,k)/2.d0/pi
        if(wannier_centers(k,i) < 0.d0) wannier_centers(k,i)=at(k,k)+wannier_centers(k,i)
     end do
  enddo
  do i=1,nbndv
     write(stdout,*) 'Center Wannier:', wannier_centers(1,i)*alat,wannier_centers(2,i)*alat,wannier_centers(3,i)*alat
  enddo
  call flush_unit(stdout)
!---now conduction subspace


  if(nbndv<nbnd_normal) then

!loop on conduction states manifolds

     nbnd_start=nbndv+1
     if(itask==1) then
        nbnd_end=nbndv+num_nbndc_set
     else
        if(num_nbnd_first==0) then
           nbnd_end=nbnd_normal
        else
           nbnd_end=num_nbndv+num_nbnd_first
        endif
     endif

!setting parameters for parallel execution
     nbnd_par=(nbnd_end-nbnd_start+1)/nproc
     if(nbnd_par*nproc< (nbnd_end-nbnd_start+1))  nbnd_par=nbnd_par+1
     nbnd_start_me=nbnd_start+nbnd_par*mpime
     nbnd_end_me=min(nbnd_start_me+nbnd_par-1,nbnd_end)


! calculate omega

        omg0  = 0.d0
        do k=1,n_oper
           do i=nbnd_start,nbnd_end
              omg0=omg0 + matsincos(i,i,k)*matsincos(i,i,k)
           enddo
        enddo

        write(stdout,*) 'LOCALIZING WANNIER FUNCTIONS:'
        call flush_unit(stdout)

! Start Iteration =====================================================

!parallelization strategy:
!copy matsincos in smaller array and distribute among processors for each operation
        n_oper_par=n_oper/nproc
        if(n_oper_par*nproc < n_oper) n_oper_par=n_oper_par+1
        oper_start=mpime*n_oper_par+1
        oper_end=min(oper_start+n_oper_par-1,n_oper)
!the following for eliminating wrong processors..
        if(oper_start>n_oper) oper_start=0
        write(stdout,*) 'OPER', n_oper_par,oper_start,oper_end
        call flush_unit(stdout)
        allocate(matsincos_par(nbnd_start:nbnd_end,nbnd_start:nbnd_end,n_oper_par))
        if(oper_start>0) then
           do k=oper_start,oper_end
              matsincos_par(nbnd_start:nbnd_end,nbnd_start:nbnd_end,k-oper_start+1)=&
                   &matsincos(nbnd_start:nbnd_end,nbnd_start:nbnd_end,k)
           enddo
        endif

        do it=1,maxiter
           do i=nbnd_start,nbnd_end
              do j=i+1,nbnd_end
                 call start_clock('trigo0')
! Construct aa
                 aa(:)=0.d0
                 if(oper_start>0) then
                    do k=1,oper_end-oper_start+1
                       aa(1)=aa(1)+matsincos_par(i,j,k)*(matsincos_par(i,i,k)-matsincos_par(j,j,k))
                       aa(2)=aa(2)+matsincos_par(i,j,k)*matsincos_par(i,j,k)-0.25d0* &
                       (matsincos_par(i,i,k)-matsincos_par(j,j,k))*(matsincos_par(i,i,k)- &
                       matsincos_par(j,j,k))
                    end do
                 endif
                 call mp_sum(aa(:))
                 if(abs(aa(2)).gt. accuracy_sincos ) then
                    theta4=-aa(1)/aa(2)
                    theta=0.25*datan(theta4)
                 elseif (abs(aa(1)).lt. accuracy_sincos ) then
                    theta=0.d0
                    aa(2)=0.d0
                 else
                    theta=pi/4.d0
                 endif
                 d2=aa(1)*dsin(4.*theta)-aa(2)*dcos(4.*theta)
                 if(d2.le.0.d0) theta=theta+pi/4.d0
!
                 cc=dcos(theta)
                 ss=dsin(theta)

                 call stop_clock('trigo0')
! update overlap matrices

                 call start_clock('vettor0')
                 if(oper_start>0) then
                    do l=1,oper_end-oper_start+1
                       ! AR

                       do k=nbnd_start,nbnd_end
                          tempi=matsincos_par(k,i,l)*cc+matsincos_par(k,j,l)*ss
                          tempj=-matsincos_par(k,i,l)*ss+matsincos_par(k,j,l)*cc
                          matsincos_par(k,i,l)=tempi
                          matsincos_par(k,j,l)=tempj
                       end do


                       do k=nbnd_start,nbnd_end
                          tempi=cc*matsincos_par(i,k,l)+ss*matsincos_par(j,k,l)
                          tempj=-ss*matsincos_par(i,k,l)+cc*matsincos_par(j,k,l)
                          matsincos_par(i,k,l)=tempi
                          matsincos_par(j,k,l)=tempj
                       end do

                    end do
                 endif

                 call mp_barrier
                 call stop_clock('vettor0')


! update U : U=UR


!                 do k=nbnd_start,nbnd_end
!                    tempi=rot_u(k,i)*cc+rot_u(k,j)*ss
!                    tempj=-rot_u(k,i)*ss+rot_u(k,j)*cc
!                    rot_u(k,i)=tempi
!                    rot_u(k,j)=tempj
!                 end do

!#else
                 call start_clock('vettor1')

                 do k=nbnd_start_me,nbnd_end_me
                    tempi=rot_u(k,i)*cc+rot_u(k,j)*ss
                    tempj=-rot_u(k,i)*ss+rot_u(k,j)*cc
                    rot_u(k,i)=tempi
                    rot_u(k,j)=tempj
                 enddo

!
!
!#endif
!
                 call mp_barrier
                 call stop_clock('vettor1')
              end do
           end do

           call print_clock('trigo0')
           call print_clock('vettor0')
           call print_clock('vettor1')


           omg1=0.d0
           if(oper_start>0) then
              do k=1,oper_end-oper_start+1
                 do i=nbnd_start,nbnd_end
                    omg1=omg1+matsincos_par(i,i,k)*matsincos_par(i,i,k)
                 end do
              end do
           endif
           call mp_sum(omg1)

           write(stdout,*) 'Spread', omg1,omg0
           call flush_unit(stdout)

           error_relative=abs(omg1-omg0)/omg0

           write(stdout,*) 'error_relative=', error_relative, ' while tresh=', tresh
           call flush_unit(stdout)

           if(error_relative .lt. tresh ) EXIT
           omg0=omg1

        end do

!re-distribute rotation matrix
        do  ip=0,nproc-1
           nbnd_start_proc=nbnd_start+nbnd_par*ip
           nbnd_end_proc=min(nbnd_start_proc+nbnd_par-1,nbnd_end)
           call mp_bcast(rot_u(nbnd_start_proc:nbnd_end_proc,:),ip)
           call mp_bcast(rot_u(nbnd_start_proc:nbnd_end_proc,:),ip)
        enddo



        call print_clock('trigo0')
        call print_clock('vettor0')

        !centers of wanniers

        do i=nbnd_start,nbnd_end
           sca_mat(:)=0.d0
           if(oper_start>0) then
              do k=1,oper_end-oper_start+1
                 sca_mat(oper_start+k-1)=matsincos_par(i,i,k)
              enddo
           endif
           call mp_sum(sca_mat(:))
           do k=1,3
              wannier_centers(k,i)=aimag(log(sca_mat(k)+(0.d0,1.d0)*sca_mat(k+3)))
              wannier_centers(k,i)=wannier_centers(k,i)*at(k,k)/2.d0/pi
              if(wannier_centers(k,i) < 0.d0) wannier_centers(k,i)=at(k,k)+wannier_centers(k,i)
           end do
        enddo
        do i=nbnd_start,nbnd_end
           write(stdout,*) 'Center Wannier:', wannier_centers(1,i)*alat,wannier_centers(2,i)*alat,wannier_centers(3,i)*alat
        enddo
        call flush_unit(stdout)

     endif


     deallocate(matsincos_par)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!if required calculates second conduction manifold

 if(nbndv<nbnd_normal .and. num_nbnd_first /=0 .and. itask == 0) then





!loop on conduction states manifolds

     nbnd_start=nbndv+num_nbnd_first+1
     nbnd_end=nbnd_normal
     nbnd_second=nbnd_end-nbnd_start+1


     allocate(tmp_mat(nbnd_second,nbnd_second))

!move matsincos to the origin
     do k=1,n_oper
        if(mod(k,nproc)==mpime) then
           do j=nbnd_start,nbnd_end
              do i=nbnd_start,nbnd_end
                 matsincos(i-nbnd_start+1,j,k)=matsincos(i,j,k)
              enddo
           enddo
           do j=1,nbnd_second
              do i=nbnd_start,nbnd_end
                 matsincos(j,i-nbnd_start+1,k)=matsincos(j,i,k)
              enddo
           enddo
        endif
     enddo





! calculate omega

        omg0  = 0.d0
        do k=1,n_oper
           if(mod(k,nproc)==mpime) then
              do i=1,nbnd_second
                 omg0=omg0 + matsincos(i,i,k)*matsincos(i,i,k)
              enddo
           endif
        enddo
        call mp_sum(omg0)

        write(stdout,*) 'LOCALIZING WANNIER FUNCTIONS:'
        call flush_unit(stdout)




! Start Iteration =====================================================

        allocate(vtmpi(nbnd_second))
        allocate(vtmpj(nbnd_second))

        do it=1,maxiter2
           do i=1,nbnd_second
              do j=i+1,nbnd_second

                 call start_clock('trigo')

! Construct aa
                 aa(:)=0.d0
                 do k=1,n_oper
                    if(mod(k,nproc)==mpime) then
                       aa(1)=aa(1)+matsincos(i,j,k)*(matsincos(i,i,k)-matsincos(j,j,k))
                       aa(2)=aa(2)+matsincos(i,j,k)*matsincos(i,j,k)-0.25d0* &
                       (matsincos(i,i,k)-matsincos(j,j,k))*(matsincos(i,i,k)-matsincos(j,j,k))
                    endif
                 end do
                 call mp_sum(aa(:))
                 if(abs(aa(2)).gt. accuracy_sincos ) then
                    theta4=-aa(1)/aa(2)
                    theta=0.25d0*datan(theta4)
                elseif (abs(aa(1)).lt. accuracy_sincos ) then
                    theta=0.d0
                    aa(2)=0.d0
                 else
                    theta=pi/4.d0
                 endif




                 d2=aa(1)*dsin(4.d0*theta)-aa(2)*dcos(4.d0*theta)
                 if(d2.le.0.d0) theta=theta+pi/4.d0



                 cc=dcos(theta)
                 ss=dsin(theta)



                 call stop_clock('trigo')
                 call start_clock('vettor')
! update overlap matrices
                 do l=1,n_oper
                    if(mod(l,nproc)==mpime) then
                    ! AR
!                    do k=1,nbnd_second
!                       tempi=matsincos(k,i,l)*cc+matsincos(k,j,l)*ss
!                       tempj=-matsincos(k,i,l)*ss+matsincos(k,j,l)*cc
!                       matsincos(k,i,l)=tempi
!                       matsincos(k,j,l)=tempj
!                    end do



                       vtmpi(:)=matsincos(:,i,l)*cc+matsincos(:,j,l)*ss
                       vtmpj(:)=-matsincos(:,i,l)*ss+matsincos(:,j,l)*cc
                       matsincos(:,i,l)=vtmpi(:)
                       matsincos(:,j,l)=vtmpj(:)



! R^+ A R


!                    call mytranspose(matsincos(:,:,l),nbnd,tmp_mat,nbnd_second,nbnd_second,nbnd_second)

                       do k=1,nbnd_second
                          tempi=cc*matsincos(i,k,l)+ss*matsincos(j,k,l)
                          tempj=-ss*matsincos(i,k,l)+cc*matsincos(j,k,l)
                          matsincos(i,k,l)=tempi
                          matsincos(j,k,l)=tempj
                       end do

!                    do k=1,nbnd_second
!                       tempi=cc*tmp_mat(k,i)+ss*tmp_mat(k,j)
!                       tempj=-ss*tmp_mat(k,i)+cc*tmp_mat(k,j)
!                       tmp_mat(k,i)=tempi
!                       tmp_mat(k,j)=tempj
!                    end do

!                     vtmpi(:)=cc*tmp_mat(:,i)+ss*tmp_mat(:,j)
!                     vtmpj(:)=-ss*tmp_mat(:,i)+cc*tmp_mat(:,j)
!                     tmp_mat(:,i)=vtmpi(:)
!                     tmp_mat(:,j)=vtmpj(:)

!                    call mytranspose(tmp_mat,nbnd_second,matsincos(:,:,l),nbnd,nbnd_second,nbnd_second)



                    endif
                  end do

! update U : U=UR


!                 do k=nbnd_start,nbnd_end
!                    tempi=rot_u(k,i)*cc+rot_u(k,j)*ss
!                    tempj=-rot_u(k,i)*ss+rot_u(k,j)*cc
!                    rot_u(k,i)=tempi
!                    rot_u(k,j)=tempj
!                 end do





                 vtmpi(1:nbnd_second)=rot_u(nbnd_start:nbnd_end,i+nbnd_start-1)*&
                      &cc+rot_u(nbnd_start:nbnd_end,j+nbnd_start-1)*ss
                 vtmpj(1:nbnd_second)=-rot_u(nbnd_start:nbnd_end,i+nbnd_start-1)*&
                                   &ss+rot_u(nbnd_start:nbnd_end,j+nbnd_start-1)*cc
                 rot_u(nbnd_start:nbnd_end,i+nbnd_start-1)=vtmpi(1:nbnd_second)
                 rot_u(nbnd_start:nbnd_end,j+nbnd_start-1)=vtmpj(1:nbnd_second)


                 call stop_clock('vettor')

!
              end do
           end do

           omg1=0.d0
           do k=1,n_oper
              if(mod(k,nproc)==mpime) then
                 do i=1,nbnd_second
                    omg1=omg1+matsincos(i,i,k)*matsincos(i,i,k)
                 end do
              endif
           end do
           call mp_sum(omg1)


           write(stdout,*) 'Spread', omg1,omg0
           call flush_unit(stdout)

           error_relative=abs(omg1-omg0)/omg0

           write(stdout,*) 'error_relative=', error_relative, ' while tresh=', tresh
           call flush_unit(stdout)

           if(error_relative .lt. tresh ) EXIT
           !!!if(abs(omg1-omg0).lt.tresh ) EXIT
           omg0=omg1

        end do
        call print_clock('trigo')
        call print_clock('vettor')


        do i=1,nbnd_second
           sca_mat(:)=0.d0
           do k=1,n_oper
              if(mod(k,nproc)==mpime) then
                 sca_mat(k)=matsincos(i,i,k)
              endif
           enddo
           call mp_sum(sca_mat(:))
           do k=1,3
              wannier_centers(k,i+nbnd_start-1)=aimag(log(sca_mat(k)+(0.d0,1.d0)*sca_mat(k+3)))
              wannier_centers(k,i+nbnd_start-1)=wannier_centers(k,i+nbnd_start-1)*at(k,k)/2.d0/pi
              if(wannier_centers(k,i+nbnd_start-1) < 0.d0) &
                   &wannier_centers(k,i+nbnd_start-1)=at(k,k)+wannier_centers(k,i+nbnd_start-1)
           end do
        enddo
        do i=1,nbnd_second
           write(stdout,*) 'Center Wannier:', wannier_centers(1,i+nbnd_start-1)*alat,&
                &wannier_centers(2,i+nbnd_start-1)*alat,wannier_centers(3,i+nbnd_start-1)*alat
        enddo
        call flush_unit(stdout)



        deallocate(tmp_mat)
        deallocate(vtmpi,vtmpj)
     endif




     deallocate(matsincos)





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!------ rotate wfc

  if(.not.gamma_only) then
    call rotate_wannier(rot_u,1)
  else
    call rotate_wannier_gamma(rot_u,1,0)
  endif



  !if required calculates the energy of the localized states
  if(l_on_ene) then
     allocate(ene_tmp(nbnd_normal))
     ene_tmp(:)=ene_loc(:)
     ene_loc(:)=0.d0
     do i=1,nbnd_normal
        do j=1,nbnd_normal
           ene_loc(i)=ene_loc(i)+rot_u(j,i)**2.d0*ene_tmp(j)
        enddo
        write(stdout,*) 'Localized energy', i,ene_loc(i)
     enddo
     deallocate(ene_tmp)
  endif





!------------update u_trans

!crot_u and u_trans are complex, rot_u is real

  allocate(crot_u(nbnd_normal,nbnd_normal),crot_u_tmp(nbnd_normal,nbnd_normal))

!  crot_u(:,:)=(0.d0,0.d0)
!  do i=1,nbnd
!    do j=1,nbnd
!       do k=1,nbnd
!rot_u is real and unitarian
!          crot_u(i,j)=crot_u(i,j)+rot_u(k,i)*u_trans(k,j)
!       enddo
!    enddo
!  enddo

  crot_u_tmp(:,:)=dcmplx(rot_u(:,:),0.d0)
  call zgemm('T','N',nbnd_normal,nbnd_normal,nbnd_normal,(1.d0,0.d0),crot_u_tmp,nbnd_normal,u_trans,nbnd_normal,&
       &(0.d0,0.d0),crot_u,nbnd_normal)


  u_trans(1:nbnd_normal,1:nbnd_normal)=crot_u(1:nbnd_normal,1:nbnd_normal)


!test u_trans:
!  do i=1,nbnd
!    do j=1,nbnd
!      j=i
!      sca=(0.d0,0.d0)
!        do k=1,nbnd
!           sca=sca+conjg(u_trans(k,i))*u_trans(k,j)
!        enddo
!        write(stdout,*) 'Test u_trans:', i,j,sca!ATTENZIONE
!     enddo
!  enddo




  deallocate(rot_u)
  deallocate(crot_u,crot_u_tmp)

! #endif

 return

end subroutine go_wannier



function fast_sin(theta,na,nb,table_sin_a, table_sin_b, table_cos_a, table_cos_b)
!this routine calculates fast the sinus
!using sin(alpha+beta)=sin(alpha)*cos(beta)+sin(beta)cos(alpha)

   USE kinds,       ONLY : DP
   USE constants,   ONLY : pi,tpi

   implicit none

   REAL(kind=DP) :: fast_sin!the calculated value for the sin
   REAL(kind=DP), INTENT(in) :: theta!the input angle in radiants
   INTEGER :: na!number of elements in the table (pi/2)/Na  from 0 to pi/2
   INTEGER :: nb!number of elements in the table (pi/2)/(Na*Nb) for 0 to (pi/2)/(Na)
   REAL(kind=DP) :: table_sin_a(na)!tabel of sin (pi/2)/Na  from 0 to pi/2
   REAL(kind=DP) :: table_sin_b(nb)!tabel of sin (pi/2)/Na  from 0 to pi/2
   REAL(kind=DP) :: table_cos_a(na)!tabel of cos (pi/2)/Na  from 0 to pi/2
   REAL(kind=DP) :: table_cos_b(nb)!tabel of cos (pi/2)/Na  from 0 to pi/2

   REAL(kind=DP) :: angle, angle_p
   INTEGER :: ia,ib
   REAL(kind=DP) :: sign, da, db

!find angle from 0 to pi/2

   angle_p= theta - real(floor(theta/tpi))*tpi
   if( angle_p <= (pi/2.d0)) then
      angle=angle_p
      sign=1.d0
   else if( angle_p <= pi) then
      sign=1.d0
      angle=pi-angle_p
   else if(angle_p <= 3.d0*pi/2.d0) then
      sign=-1.d0
      angle=angle_p-pi
   else
      sign=-1.d0
      angle=tpi-angle_p
   endif

!determines ia, ib

   da=pi/(2.d0*real(na))
   db=da/real(nb)
   ia=floor(angle/(da))
   ib=floor((angle-real(ia)*da)/db)
   ia=ia+1
   ib=ib+1

!sin(a+b)
   fast_sin=sign*(table_sin_a(ia)*table_cos_b(ib)+table_sin_b(ib)*table_cos_a(ia))

   return



 end function fast_sin



function fast_cos(theta,na,nb,table_sin_a, table_sin_b, table_cos_a, table_cos_b)
!this routine calculates fast the sinus
!using cos(alpha+beta)=cos(alpha)*cos(beta)-sin(alpha)*sin(beta)

   USE kinds,       ONLY : DP
   USE constants,   ONLY : pi,tpi

   implicit none

   REAL(kind=DP) :: fast_cos!the calculated value for the cos
   REAL(kind=DP), INTENT(in) :: theta!the input angle in radiants
   INTEGER :: na!number of elements in the table (pi/2)/Na  from 0 to pi/2
   INTEGER :: nb!number of elements in the table (pi/2)/(Na*Nb) for 0 to (pi/2)/(Na)
   REAL(kind=DP) :: table_sin_a(na)!tabel of sin (pi/2)/Na  from 0 to pi/2
   REAL(kind=DP) :: table_sin_b(nb)!tabel of sin (pi/2)/Na  from 0 to pi/2
   REAL(kind=DP) :: table_cos_a(na)!tabel of cos (pi/2)/Na  from 0 to pi/2
   REAL(kind=DP) :: table_cos_b(nb)!tabel of cos (pi/2)/Na  from 0 to pi/2

   REAL(kind=DP) :: angle, angle_p
   INTEGER :: ia,ib
   REAL(kind=DP) :: sign, da, db

!find angle from 0 to pi/2

   angle_p= theta - real(floor(theta/tpi))*tpi
   if( angle_p <= (pi/2.d0)) then
      angle=angle_p
      sign=1.d0
   else if( angle_p <= pi) then
      sign=-1.d0
      angle=pi-angle_p
   else if(angle_p <= 3.d0*pi/2.d0) then
      sign=-1.d0
      angle=angle_p-pi
   else
      sign=1.d0
      angle=tpi-angle_p
   endif

!determines ia, ib

   da=pi/(2.d0*real(na))
   db=da/real(nb)
   ia=floor(angle/(da))
   ib=floor((angle-real(ia)*da)/db)

   ia=ia+1
   ib=ib+1

!cos(a+b)
   fast_cos=sign*(table_cos_a(ia)*table_cos_b(ib)-table_sin_a(ia)*table_sin_b(ib))

   return



 end function fast_cos


 function fast_atan(tg,na,nb,table_sin_a, table_sin_b, table_cos_a, table_cos_b)
!this routine calculates fast the arctan using a bisection algorithm

   USE kinds,       ONLY : DP
   USE constants,   ONLY : pi

   implicit none

   REAL(kind=DP) :: fast_atan!the calculated value for the atan
   REAL(kind=DP), INTENT(in) :: tg!the input tan
   INTEGER :: na!number of elements in the table (pi/2)/Na  from 0 to pi/2
   INTEGER :: nb!number of elements in the table (pi/2)/(Na*Nb) for 0 to (pi/2)/(Na)
   REAL(kind=DP) :: table_sin_a(na)!tabel of sin (pi/2)/Na  from 0 to pi/2
   REAL(kind=DP) :: table_sin_b(nb)!tabel of sin (pi/2)/Na  from 0 to pi/2
   REAL(kind=DP) :: table_cos_a(na)!tabel of cos (pi/2)/Na  from 0 to pi/2
   REAL(kind=DP) :: table_cos_b(nb)!tabel of cos (pi/2)/Na  from 0 to pi/2

   INTEGER, PARAMETER :: n=20!number of bisections

   REAL(kind=DP) :: sign, tang
   INTEGER :: i
   REAL(kind=DP) :: ang, delta
   REAL(kind=DP) :: tan_try
   REAL(kind=DP), EXTERNAL :: fast_sin, fast_cos


   if(tg >= 0.d0) then
      sign=1.d0
      tang=tg
   else
      sign=-1.d0
      tang=-tg
   endif


   ang=pi/4.d0
   delta=pi/4.d0
   do i=1,n
      delta=delta/2.d0
      tan_try=fast_sin(ang,na,nb,table_sin_a, table_sin_b, table_cos_a, table_cos_b)/&
               &fast_cos(ang,na,nb,table_sin_a, table_sin_b, table_cos_a, table_cos_b)
      if(tang >= tan_try) then
         ang=ang+delta
      else
         ang=ang-delta
      endif
   enddo

   fast_atan=sign*ang

   return
 end function fast_atan
