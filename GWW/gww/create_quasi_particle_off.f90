!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

 SUBROUTINE create_quasi_particles_off(options,qp,se)
!given the expansion coeffcients, calculates in a perturbative
!way without self-consistency correction the quasi-particles energies
!using the whole self_energy matrix included the off-diagonal terms
!relavant arrays are already allocates and set by subroutine create


   USE io_global,        ONLY : stdout
   USE basic_structures, ONLY : wannier_u, free_memory
   USE expansion,        ONLY : self_expansion, value_on_frequency, derivative_on_frequency,value_on_frequency_complex_off
   USE input_gw,         ONLY : input_options
   USE constants,        ONLY : tpi, RYTOEV
   USE energies_gww,         ONLY : quasi_particles
   USE kinds,            ONLY : DP

   implicit none

   TYPE(input_options) :: options! for prefix
   TYPE(quasi_particles) :: qp!the descriptor to be build
   TYPE(self_expansion) :: se!the descriptor for the multipole expansion


   INTEGER :: ii,jj,kk, it,is
   TYPE(wannier_u) :: uu
   COMPLEX(kind=DP) :: om
   REAL(kind=DP)  :: offset,sca,norm
   COMPLEX(kind=DP), ALLOCATABLE :: amp(:),hmat(:,:),hval(:),hvec(:,:)
   INTEGER :: lwork,info
   COMPLEX(kind=DP), ALLOCATABLE :: work(:)
   REAL(kind=DP), ALLOCATABLE :: rwork(:)
   INTEGER :: ivec
   INTEGER :: ilo,ihi
   REAL(kind=DP), ALLOCATABLE :: scale(:),rconde(:),rcondv(:)
   REAL(kind=DP) :: abnrm


   if(.not.options%whole_s) return

!allocate and set up off diagonal arrays 

   qp%whole_s=.true.

   allocate(qp%ene_gw_off(qp%max_i,qp%nspin))
   allocate(qp%eigen_gw_off(qp%max_i,qp%max_i,qp%nspin)) 
   allocate(qp%ene_dft_xc_off(qp%max_i,qp%max_i,qp%nspin))
   allocate(qp%ene_x_off(qp%max_i,qp%max_i,qp%nspin))

!read in DFT energies
!read in whole dft xc and Fock X matrices

   call read_data_pw_u(uu,options%prefix)

   call read_data_pw_exchange_off(qp%ene_x_off,qp%max_i,options%prefix,qp%nspin)
   do is=1,qp%nspin
      call read_data_pw_dft_xc_off(qp%ene_dft_xc_off(1,1,is),qp%max_i,options%prefix,is)
   enddo

!allocate arrays

   allocate(amp(qp%max_i))
   allocate(hmat(qp%max_i,qp%max_i),hvec(qp%max_i,qp%max_i),hval(qp%max_i))

!loop on spin

   do is=1,uu%nspin


      if(uu%nums > uu%nums_occ(is)) then
         offset=-(uu%ene(uu%nums_occ(is)+1,is)+uu%ene(uu%nums_occ(is),is))/2.d0
      else
         offset=-uu%ene(uu%nums_occ(is),is)
      endif

 !  call free_memory(uu)

      do ii=1,qp%max_i
!set up starting complex energy and amplitude
         om=dcmplx(qp%ene_dft_ks(ii,is)+offset,0.d0)
         amp(:)=(0.d0,0.d0)
         amp(ii)=(1.d0,0.d0)
!self-consistency loop
         do it=1,10


!set up hamiltonian like matrix
!self_energy
            hmat=(0.d0,0.d0)
            do jj=1,qp%max_i
               do kk=1,qp%max_i
                  call value_on_frequency_complex_off(se,kk,jj,om,hmat(kk,jj),is) 
               enddo
            enddo
!H0            
            do jj=1,qp%max_i
               hmat(jj,jj)=hmat(jj,jj)+qp%ene_dft_ks(jj,is)+offset
            enddo
            
            hmat(1:qp%max_i,1:qp%max_i)=hmat(1:qp%max_i,1:qp%max_i)-qp%ene_dft_xc_off(1:qp%max_i,1:qp%max_i,is)
            hmat(1:qp%max_i,1:qp%max_i)=hmat(1:qp%max_i,1:qp%max_i)+qp%ene_x_off(1:qp%max_i,1:qp%max_i,is)

!find eigenvalues/vectors

!for compatibility with essl on aix it must use zgeevx instead (sigh)
!            allocate(rwork(2*qp%max_i))
!            allocate(work(1))
!            call ZGEEV('N','V',qp%max_i,hmat,qp%max_i,hval,hvec,qp%max_i,hvec,qp%max_i,work,-1,rwork,info)
!            if(info/=0) then
!               write(stdout,*) 'Problems with ZGEEV:', info
!               FLUSH(stdout)
!               stop
!            endif
!            lwork=int(work(1))
!            deallocate(work)
!            allocate(work(lwork))
!            call ZGEEV('N','V',qp%max_i,hmat,qp%max_i,hval,hvec,qp%max_i,hvec,qp%max_i,work,lwork,rwork,info)
!            if(info/=0) then
!               write(stdout,*) 'Problems with ZGEEV:', info
!               FLUSH(stdout)
!               stop
!            endif
!           deallocate(work)
!            deallocate(rwork)


            allocate(scale(qp%max_i),rconde(qp%max_i),rcondv(qp%max_i))
            allocate(rwork(2*qp%max_i))                                                                                    
            allocate(work(1))                 
            call zgeevx('N','N','V','N',qp%max_i,hmat,qp%max_i,hval,hvec,qp%max_i,hvec,qp%max_i,ilo,ihi,scale,&
                 &abnrm,rconde,rcondv,work,-1,rwork,info)
            if(info/=0) then                                                                                               
               write(stdout,*) 'Problems with ZGEEVX:', info
               FLUSH(stdout)                                                                                     
               stop                                                                                                        
            endif                                                                                                          
            lwork=int(work(1))                                                                                             
            deallocate(work)                                                                                               
            allocate(work(lwork))  
            call zgeevx('N','N','V','N',qp%max_i,hmat,qp%max_i,hval,hvec,qp%max_i,hvec,qp%max_i,ilo,ihi,scale,&
                 &abnrm,rconde,rcondv,work,lwork,rwork,info)
            if(info/=0) then
               write(stdout,*) 'Problems with ZGEEVX:', info
                FLUSH(stdout)
                stop                                                                                                      
            endif                                                                                                          
            deallocate(work)
            deallocate(rwork)                                                                                              
            deallocate(scale,rconde,rcondv)
 

!print eigenvalues

            do jj=1,qp%max_i
               write(stdout,*) 'COMPLEX EN:',jj,is,it,hval(jj)
            enddo
            FLUSH(stdout)


!find the vector most close to the previous one 
            norm=0.d0
            ivec=1
            do jj=1,qp%max_i
               sca=0.d0
               do kk=1,qp%max_i
                  sca=sca+conjg(hvec(kk,jj))*amp(kk)+hvec(kk,jj)*conjg(amp(kk))
               enddo
               if(sca>norm) then
                  norm=sca
                  ivec=jj
               endif
            enddo
!update energy and amplitude
            write(stdout,*) 'NEW VECTOR INDEX', ivec
            om=hval(ivec)
            amp(1:qp%max_i)=hvec(1:qp%max_i,ivec)


         enddo!iterations
!put final results on qp object
         qp%ene_gw_off(ii,is)=om-offset
         qp%eigen_gw_off(1:qp%max_i,ii,is)=amp(1:qp%max_i)

      enddo!states

   enddo!spin
   call free_memory(uu)
   deallocate(amp,hmat,hval,hvec)
   return
 END SUBROUTINE create_quasi_particles_off
