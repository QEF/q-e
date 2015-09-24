!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

  SUBROUTINE self_energy(i,j,sene,time,qm,uu,gf,ww)
!this subroutine calculates the terms, in imaginary time
!<\Psi_i|\Sigma(it)|\Psi_j> 
!=O^{P}_n,kl G_{lm}W_{n,o} O^{P}_o,mp U_ki U^{+}_j,p

   USE kinds,                ONLY : DP
   USE io_global,            ONLY : stdout 
   USE basic_structures,     ONLY : wannier_u, q_mat
   USE green_function,       ONLY : green
   USE polarization,         ONLY : polaw

   implicit none


   INTEGER  :: i,j !which element of self enrgy to be calculated
   COMPLEX(kind=DP) :: sene!self energy element
   REAL(kind=DP) :: time!in output time correspondig to the calculated self energy
   TYPE(q_mat)  :: qm!descriptors of overlaps of othonormalized wannier producs with wannier products
   TYPE(wannier_u) :: uu!descriptor of transformation matrix from KS states to wanniers
   TYPE(green) :: gf!descriptor of green function 
   TYPE(polaw) :: ww!descriptor of dressed interaction

   INTEGER :: k,l,m,n,o,p
   INTEGER :: nw,ow
   REAL(kind=DP) :: o_n,o_o

!check consistency

   if(.not.gf%ontime ) then
     write(stdout,*) 'Routine self_energy: imaginary times GF required'
     stop
   endif

   if(.not.ww%ontime) then
     write(stdout,*) 'Routine self_energy: imaginary times WW required'
!     stop
   endif


!the following has been commented for using with remainder calculation
!   if(gf%time /= ww%time) then
!     write(stdout,*) 'Routine self_energy: same imaginary times required'
!     stop
!   endif
   if(gf%nums /= uu%nums) then
     write(stdout,*) 'Routine self_energy: same nums required'
     stop
   endif
   if(qm%numpw /= ww%numpw) then
     write(stdout,*) 'Routine self_energy: same numpw required'
     stop
   endif


   time=ww%time
   sene=(0.d0,0.d0)
   do n=1,ww%numpw!loop on orthonormalized wannier products
      do o=1,ww%numpw!loop on orthonormalized wannier products
         do nw=1,qm%wp(n)%numij
            do ow=1,qm%wp(o)%numij

               k=qm%wp(n)%ij(1,nw)
               l=qm%wp(n)%ij(2,nw)
               m=qm%wp(o)%ij(1,ow)
               p=qm%wp(o)%ij(2,ow)

               o_n=qm%wp(n)%o(nw)
               o_o=qm%wp(o)%o(ow)
               sene=sene+o_n*gf%gf(l,m,1)*ww%pw(n,o)*o_o*conjg(uu%umat(i,k,1))*uu%umat(j,p,1)
               if(k/=l) then
                 sene=sene+o_n*gf%gf(k,m,1)*ww%pw(n,o)*o_o*conjg(uu%umat(i,l,1))*uu%umat(j,p,1)
               endif

               if(m/=p) then
                 sene=sene+o_n*gf%gf(l,p,1)*ww%pw(n,o)*o_o*conjg(uu%umat(i,k,1))*uu%umat(j,m,1)
               endif

               if(m/=p .AND. k/=l ) then
                 sene=sene+o_n*gf%gf(k,p,1)*ww%pw(n,o)*o_o*conjg(uu%umat(i,l,1))*uu%umat(j,m,1)
               endif
             end do
          enddo  
       enddo
    enddo
    sene=sene*(0.d0,1.d0)
    return
 END SUBROUTINE
 



  SUBROUTINE self_energy_contraction(i,j,sene,time,cr,gf,ww)
!this subroutine calculates the terms, in imaginary time using contraction array
!<\Psi_i|\Sigma(it)|\Psi_j>
!G_{lm}W_{n,o} Q^{P}_{n,l,i}*conjg(Q^{P}_{o,m,j}

   USE kinds,                ONLY : DP
   USE io_global,            ONLY : stdout
   USE compact_product
   USE green_function,       ONLY : green
   USE polarization,         ONLY : polaw

   implicit none


   INTEGER  :: i,j !which element of self enrgy to be calculated
   COMPLEX(kind=DP) :: sene!self energy element
   REAL(kind=DP) :: time!in output time correspondig to the calculated self energy
   TYPE(contraction) :: cr!description of contracted terms
   TYPE(green) :: gf!descriptor of green function
   TYPE(polaw) :: ww!descriptor of dressed interaction
   COMPLEX(kind=DP), ALLOCATABLE :: qg(:,:)!for the product Q^{P}_{n,l,i}G{l,m}

   INTEGER :: k,l,m,n,o,p

!check consistency

   if(.not.gf%ontime) then
     write(stdout,*) 'Routine self_energy: imaginary times GF required'
     stop
   endif

   if(.not.ww%ontime) then
     write(stdout,*) 'Routine self_energy: imaginary times WW required'
!     stop
   endif




!the following has been commented for remainder calculation
!   if(gf%time /= ww%time) then
!     write(stdout,*) 'Routine self_energy: same imaginary times required'
!     stop
!   endif
   if(gf%nums /= cr%nums) then
     write(stdout,*) 'Routine self_energy: same nums required'
     stop
   endif
   if(cr%numpw /= ww%numpw) then
     write(stdout,*) 'Routine self_energy: same numpw required'
     stop
   endif

   allocate(qg(cr%numpw,cr%nums))
   qg(:,:)=(0.d0,0.d0)

   do n=1,cr%numpw!loop on orthonormalized wannier products
     do m=1,cr%nums
       do l=1,cr%numl(n)
         qg(n,m)=qg(n,m)+cr%q(n,l,i)*gf%gf(cr%l(l,n),m,1)
       enddo
     enddo
   enddo

   sene=(0.d0,0.d0)
   do n=1,cr%numpw!loop on orthonormalized wannier products
      do o=1,cr%numpw!loop on orthonormalized wannier products
         do m=1,cr%numl(o)
           sene=sene+qg(n,cr%l(m,o))*ww%pw(n,o)*conjg(cr%q(o,m,j))
         enddo
       enddo
    enddo
  !  if(sene==0.d0) write(*,*) 'OPS', i
    
        

   time=ww%time

    sene=sene*(0.d0,1.d0)
    deallocate(qg)
    return
 END SUBROUTINE



  SUBROUTINE self_energy_remainder(i,rem,time,wp,ww)
!this subroutine calculates the remainders  for negative imaginary time
!<\Psi_i|\Sigma(it)|\Psi_j>
!=Sigma^R_v(it)=\sum wwp(i,j,v)ww(i,j,it)

   USE kinds,                ONLY : DP
   USE io_global,            ONLY : stdout
   USE basic_structures,     ONLY : wp_psi
   USE polarization,         ONLY : polaw

   implicit none


   INTEGER  :: i!which element of self energy  remainder to be calculated
   COMPLEX(kind=DP) :: rem!self energy remainder element
   REAL(kind=DP) :: time!in output time correspondig to the calculated self energy, just for control
   TYPE(wp_psi) :: wp!descriptor of product wp wp psi psi
   TYPE(polaw) :: ww!descriptor of dressed interaction

   INTEGER :: iw,jw

   if(.not.ww%ontime) then
      write(stdout,*) 'Routine self_energy_remainder: imaginary times required'
!      stop
   endif
   
   if(wp%numpw /= ww%numpw) then
      write(stdout,*) 'Routine self_energy_remainder: same numpw required'
      stop
   endif
   
   if(i > wp%nums_psi) then
      write(stdout,*) 'Routine self_energy_remainder:  i too large',i,wp%nums_psi
      stop
   endif

   rem=(0.d0,0.d0)

   do iw=1,ww%numpw
      do jw=1,ww%numpw
         rem=rem+ww%pw(iw,jw)*wp%wwp(iw,jw,i)
      enddo
   enddo

         
   return
 END SUBROUTINE self_energy_remainder


  SUBROUTINE set_data_wp_psi_cutoff(pw_red,pw,wpi)
!this subroutine allocates and set the array pw_red
!which contains the corresponding elements of wpwp_psi

   USE kinds,                ONLY : DP
   USE io_global,            ONLY : stdout
   USE basic_structures,     ONLY : wp_psi_cutoff_index
   USE polarization,         ONLY : polaw

   implicit none

   COMPLEX(kind=DP), DIMENSION(:), POINTER :: pw_red!array for contracted data
   TYPE(polaw) :: pw!data to be contracted
   TYPE(wp_psi_cutoff_index) :: wpi !indices

   INTEGER i,j

   allocate(pw_red(wpi%numpwpw))

   write(stdout,*) 'Number NUMPWPW', wpi%numpwpw

   do i=1,wpi%numpwpw
      if(wpi%index(1,i) /= wpi%index(2,i)) then
         pw_red(i)=pw%pw(wpi%index(1,i),wpi%index(2,i))+pw%pw(wpi%index(2,i),wpi%index(1,i))
      else
         pw_red(i)=pw%pw(wpi%index(1,i),wpi%index(1,i))
      endif
   enddo

   WRITE(stdout,*) 'PW_RED OUT', pw_red(1)

   return

 END SUBROUTINE set_data_wp_psi_cutoff


 SUBROUTINE self_energy_remainder_cutoff(state,rem,wp,pw_red)
!this subroutine calculates the remainders  for negative imaginary time
!<\Psi_i|\Sigma(it)|\Psi_j>
!=Sigma^R_v(it)=\sum wwp(i,j,v)ww(i,j,it)
!using a reduced set of data

   USE kinds,                ONLY : DP
   USE io_global,            ONLY : stdout
   USE basic_structures,     ONLY : wp_psi_cutoff_data

   implicit none


   INTEGER  :: state!which element of self energy  remainder to be calculated
   COMPLEX(kind=DP) :: rem!self energy remainder element
   COMPLEX(kind=DP), DIMENSION(:), POINTER :: pw_red!contracted polarization/interaction array
   TYPE(wp_psi_cutoff_data) :: wp!descriptor of product wp wp psi psi

   INTEGER :: i

   rem=(0.d0,0.d0)

   do i=1,wp%numpwpw
      rem=rem+pw_red(i)*wp%wwp(i,state)
   enddo

   return

 END SUBROUTINE self_energy_remainder_cutoff



  SUBROUTINE self_energy_contraction_state(i,j,sene,time,cri,crs,gf,ww)
!this subroutine calculates the terms, in imaginary time using contraction array
!<\Psi_i|\Sigma(it)|\Psi_j>
!G_{lm}W_{n,o} Q^{P}_{n,l,i}*conjg(Q^{P}_{o,m,j}
!uses state contractions
!ONLY DIAGONAL TERMS IMPLEMENTED YET

   USE kinds,                ONLY : DP
   USE io_global,            ONLY : stdout
   USE compact_product
   USE green_function,       ONLY : green
   USE polarization,         ONLY : polaw

   implicit none


   INTEGER  :: i,j !which element of self enrgy to be calculated
   COMPLEX(kind=DP) :: sene!self energy element
   REAL(kind=DP) :: time!in output time correspondig to the calculated self energy
   TYPE(contraction_index) :: cri!index description of contracted terms
   TYPE(contraction_state) :: crs!state contraction 
   TYPE(green) :: gf!descriptor of green function
   TYPE(polaw) :: ww!descriptor of dressed interaction
   REAL(kind=DP), ALLOCATABLE :: qg(:,:)!for the product Q^{P}_{n,l,i}G{l,m}
   REAL(kind=DP), ALLOCATABLE :: qg_t(:,:)
   REAL(kind=DP), ALLOCATABLE :: gf_t(:,:)
   REAL(kind=DP), ALLOCATABLE :: crsq_t(:,:)

   REAL(kind=DP), ALLOCATABLE :: tmp_q(:),  tmp_m(:), tmp_w(:),tmp_m2(:)
   INTEGER, ALLOCATABLE :: cri_index(:)
   INTEGER :: k,l,m,n,o,p

!check consistency

   if(.not.gf%ontime) then
     write(stdout,*) 'Routine self_energy: imaginary times GF required'
     stop
   endif

   if(.not.ww%ontime) then
     write(stdout,*) 'Routine self_energy: imaginary times WW required'
!     stop
   endif

   if( i/=j) then
      write(stdout,*) 'Routine self_energy: ONLY DIAGONAL TERMS IMPLEMETED YET'
      stop
   endif




   if(gf%nums /= cri%nums) then
     write(stdout,*) 'Routine self_energy: same nums required'
     stop
   endif
   if(cri%numpw /= ww%numpw) then
     write(stdout,*) 'Routine self_energy: same numpw required'
     stop
   endif

    write(stdout,*) 'Self-energy 0',gf%factor,ww%factor
    FLUSH(stdout)

   allocate( qg ( cri%numpw, cri%nums) )
   allocate( qg_t( cri%nums, cri%numpw ) )
   allocate( gf_t( cri%nums, cri%nums ) )

   CALL mytranspose( gf%gf_p, cri%nums, gf_t, cri%nums, cri%nums, cri%nums )   


   qg_t(:,:)=0.d0

   do n=1,cri%numpw!loop on orthonormalized wannier products
     do l=1,cri%numl(n)
       !do m=1,cri%nums
       !  qg_t(m,n)=qg_t(m,n)+crs%q(n,l)*gf_t(m,cri%l(l,n))
       !enddo
       CALL daxpy( cri%nums, crs%q(n,l), gf_t( 1, cri%l(l,n) ), 1, qg_t(1,n), 1 )  
     enddo
   enddo

   CALL mytranspose( qg_t, cri%nums, qg, cri%numpw, cri%nums, cri%numpw )   
   !do n=1,cri%nums
   !   do m=1,cri%numpw
   !      qg( m, n )  = qg_t( n, m )
   !   end do
   !end do

   DEALLOCATE( qg_t, gf_t )

   write(stdout,*) 'Self-energy 1'
   FLUSH(stdout)

   sene=(0.d0,0.d0)





   allocate(tmp_w(cri%numpw))
!   allocate(tmp_qq(cri%numpw, cri%nums))


   allocate(crsq_t(cri%nums,cri%numpw))
 !  do o=1,cri%numpw
 !     do m=1,cri%numl(o)
 !        crsq_t(m,o)=cmplx(crs%q(o,m),0.d0)
 !     enddo
 !  enddo

   call mytranspose(crs%q,cri%numpw,crsq_t,cri%nums,cri%numpw,cri%nums)

   do o=1,cri%numpw

      tmp_w(:)=0.d0
      call daxpy(cri%numpw,1.d0,ww%pw(:,o),1,tmp_w,1)

      allocate(tmp_m(cri%nums))

       call dgemv('T',cri%numpw,cri%nums,1.d0,qg,cri%numpw,tmp_w,1,0.d0,tmp_m,1)

       allocate(tmp_m2(cri%numl(o)))
       tmp_m2(1:cri%numl(o))=crsq_t(1:cri%numl(o),o)
       


       allocate(tmp_q(cri%numl(o)))
       allocate(cri_index(cri%numl(o)))
       do m=1,cri%numl(o)
          cri_index(m)=cri%l(m,o)
       enddo

       do m=1,cri%numl(o)
          tmp_q(m)=tmp_m(cri_index(m))*tmp_m2(m)
       enddo
       sene=sene+sum(tmp_q(1:cri%numl(o)))*gf%factor*ww%factor
       deallocate(tmp_q)
       deallocate(tmp_m2)
       deallocate(cri_index)
      deallocate(tmp_m)

   enddo

    write(stdout,*) 'Self-energy 3', gf%factor,ww%factor
    FLUSH(stdout)

   deallocate(crsq_t)

   deallocate(tmp_w)


   time=ww%time


    sene=sene*(0.d0,1.d0)
    deallocate(qg)
    return
  END SUBROUTINE self_energy_contraction_state








 
        

