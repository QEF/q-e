!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

 subroutine wannier_uterms(n_set,l_square,lzero, orthonorm, ecutoff)

!this subroutine
!calculates the products <wiwj(r_1)| 1/|r_1-r_2| wi'j'(r_2)>
!if required used truncation formula of Onida, PRB 62, 4927 (2000)

  USE io_global,            ONLY : stdout, ionode
  USE io_files,             ONLY : prefix, tmp_dir, diropn
  use mp_pools,            ONLY : nproc_pool, me_pool
  use mp_world,            ONLY : world_comm
  USE kinds,    ONLY : DP
  USE gvect
  USE basis
  USE klist
  USE constants, ONLY : e2, pi, tpi, fpi
  USE wvfct,     ONLY : npwx, npw, nbnd
  USE gvecw,     ONLY : gcutw
   USE cell_base, ONLY: at, alat, tpiba, omega, tpiba2
  USE wannier_gw
  USE mp, ONLY : mp_sum
  USE control_flags,        ONLY : gamma_only

  implicit none

  INTEGER, EXTERNAL :: find_free_unit
  INTEGER, INTENT(in)  :: n_set  !defines the number of states to be read from disk at the same time
  LOGICAL, INTENT(in)  :: l_square!if true calculate v^1/2 for the symmetric dielectric matrix
  LOGICAL, INTENT(in)  :: lzero!if true put to zero the G=0,G=0 of v
  INTEGER, INTENT(in)  :: orthonorm!if ==1 opens orthonormalized products of wannier file, if==2 reduced one
  REAL(kind=DP), INTENT(in) :: ecutoff!cutoff in Rydberg for g sum

  INTEGER :: iungprod, iunuterms

  !  --- Internal definitions ---

   COMPLEX(kind=DP), ALLOCATABLE :: tmpspacei(:,:)
   COMPLEX(kind=DP), ALLOCATABLE :: tmpspacej(:,:)
   REAL(kind=DP), ALLOCATABLE ::fac(:)
   REAL(kind=DP), ALLOCATABLE :: uterms(:,:)!temporary reading array
   INTEGER :: iw,jw, iiw,jjw,jw_begin
   INTEGER :: ig
   LOGICAL :: exst
   REAL (kind=DP) :: qq
   INTEGER :: igk0(npwx)
   REAL(kind=dp) :: g2kin_bp(npwx)
   INTEGER :: npw0
   REAL(kind=DP) :: exxdiv
   INTEGER :: ngm_max
   INTEGER :: iw_min,iw_max,jw_min,jw_max
   COMPLEX(kind=DP), ALLOCATABLE :: umat_tmp(:,:)
 


   write(stdout,*) 'Routine wannier_uterms : start'
  

!   exxdiv=exx_divergence_new()

!determine ngm_max
   ngm_max=0
   do ig=1,ngm
      if(gg(ig)*tpiba2 >= ecutoff) exit
      ngm_max=ngm_max+1
   enddo
   
   write(stdout,*) 'NGM MAX:', ngm_max, ngm


  


  
! reads wfcs from iunwfc

   CALL gk_sort(xk(1,1),ngm,g,gcutw,npw0,igk0,g2kin_bp)
   allocate(uterms(numw_prod,numw_prod))
   allocate(tmpspacei(max_ngm,n_set),tmpspacej(max_ngm,n_set),fac(max_ngm))
   allocate(umat_tmp(n_set,n_set))
   iungprod = find_free_unit()
   if(orthonorm==0) then
      CALL diropn( iungprod, 'wiwjwfc', max_ngm*2, exst )
   else if(orthonorm==1) then
      CALL diropn( iungprod, 'wiwjwfc_on', max_ngm*2, exst )
   else
      CALL diropn( iungprod, 'wiwjwfc_red', max_ngm*2, exst )
   endif
   !sets factors terms
 !sets factors terms
!this has already  been called call exx_grid_init()


   if(l_truncated_coulomb) then
      do ig=1,max_ngm
      
         qq = g(1,ig)**2.d0 + g(2,ig)**2.d0 + g(3,ig)**2.d0
      
     
         if (qq > 1.d-8) then
            fac(ig)=(e2*fpi/(tpiba2*qq))*(1.d0-dcos(dsqrt(qq)*truncation_radius*tpiba))
         else
            fac(ig)=e2*fpi*(truncation_radius**2.d0/2.d0)
         endif
         
      enddo
      fac(:)=fac(:)/omega
   else
      fac(:)=0.d0
      fac(1:npw)=vg_q(1:npw)
   endif
    
 
   if(lzero .and. gstart==2) fac(1)=0.d0
   if(l_square) fac(:)=dsqrt(fac(:))
   !open output file
   if(ionode) then
      iunuterms =  find_free_unit()
       open( unit= iunuterms, file=trim(tmp_dir)//trim(prefix)//'.uterms', status='unknown',form='unformatted')
    endif

   uterms(:,:)=0.d0
   do iiw=1,ceiling(real(numw_prod)/real(n_set))
      write(stdout,*) 'uterms iiw', iiw
      do iw=(iiw-1)*n_set+1,min(iiw*n_set,numw_prod)
         CALL davcio(tmpspacei(:,iw-(iiw-1)*n_set),max_ngm*2,iungprod,iw,-1)
         
        if(gamma_only .and. gstart == 2) then
           tmpspacei(1,iw-(iiw-1)*n_set)=dble(tmpspacei(1,iw-(iiw-1)*n_set))
        endif
      enddo
      iw_min=(iiw-1)*n_set+1
      iw_max=min(iiw*n_set,numw_prod)

      do jjw=iiw,ceiling(real(numw_prod)/real(n_set))
         write(stdout,*) 'uterms jjw', jjw
         do jw=(jjw-1)*n_set+1,min(jjw*n_set,numw_prod)
            CALL davcio(tmpspacej(:,jw-(jjw-1)*n_set),max_ngm*2,iungprod,jw,-1)
            if(gamma_only .and. gstart == 2) then
               tmpspacej(1,jw-(jjw-1)*n_set)=dble(tmpspacej(1,jw-(jjw-1)*n_set))
            endif
         enddo
         jw_min=(jjw-1)*n_set+1
         jw_max=min(jjw*n_set,numw_prod)



!!!!!!!!!!!!!!!!!!!!!!!!!!!
!uses  blas routine

         do jw=1,jw_max-jw_min+1
            tmpspacej(1:ngm_max,jw)= tmpspacej(1:ngm_max,jw)*fac(1:ngm_max)
            if(gstart==2) tmpspacej(1,jw)=0.5d0*tmpspacej(1,jw)
         enddo
         call zgemm('C','N',n_set,n_set,ngm_max,(1.d0,0.d0),tmpspacei,max_ngm,tmpspacej,max_ngm,(0.d0,0.d0),umat_tmp,n_set)
         call mp_sum(umat_tmp(:,:),world_comm)
         do iw=iw_min,iw_max
            do jw=jw_min,jw_max
               uterms(iw,jw)=2.d0*dble(umat_tmp(iw-iw_min+1,jw-jw_min+1))
               uterms(jw,iw)=uterms(iw,jw)
            enddo
         enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!






!        do iw=(iiw-1)*n_set+1,min(iiw*n_set,numw_prod)
!           if(iiw==jjw) then
!             jw_begin=iw
!           else
!             jw_begin=(jjw-1)*n_set+1
!           endif
!           do jw=jw_begin,min(jjw*n_set,numw_prod)
!             uterms(iw,jw)=0.d0
!             if(.not.gamma_only) then
!                do ig=1,ngm_max
!                   uterms(iw,jw)=uterms(iw,jw) + dble(fac(ig)*&
!                        &conjg(tmpspacei(ig,iw-(iiw-1)*n_set))*tmpspacej(ig,jw-(jjw-1)*n_set))
!                enddo
!             else
!                do ig=1,ngm_max
!                   uterms(iw,jw)=uterms(iw,jw) + 2.d0*dble(fac(ig)*&
!                        &conjg(tmpspacei(ig,iw-(iiw-1)*n_set))*tmpspacej(ig,jw-(jjw-1)*n_set))
!                enddo
!                if(gstart==2) then
!                   uterms(iw,jw)=uterms(iw,jw)-dble(fac(1)*&
!           &conjg(tmpspacei(1,iw-(iiw-1)*n_set))*tmpspacej(1,jw-(jjw-1)*n_set))
!                endif
!         
!             endif
!             call reduce(1, uterms(iw,jw))
!             uterms(jw,iw)=uterms(iw,jw)
!          enddo
!       enddo

    enddo
 enddo
 if(ionode) then
      do iw=1,numw_prod
         write(iunuterms) uterms(iw,1:iw)
      enddo
      close(iunuterms)
   endif
   close(iungprod)
   deallocate(tmpspacei,tmpspacej,fac,uterms)
   deallocate(umat_tmp)
   return
 end subroutine wannier_uterms


 subroutine calculate_vg0 
!this subroutine calculate the G=0 element of the Coulomb interatction
!by integrating over q
!and write them on disk

   USE wannier_gw, ONLY : vg_q
   USE wvfct,    ONLY : npw,npwx
   USE gvect
   USE cell_base, ONLY: at, alat, tpiba, omega, tpiba2, bg
   USE constants, ONLY : e2, pi, tpi, fpi
   USE io_global, ONLY : stdout, ionode
   USE io_files,  ONLY : prefix, tmp_dir, diropn

   implicit none

   INTEGER, EXTERNAL :: find_free_unit
   INTEGER :: ig,iun
   INTEGER, PARAMETER :: n_int=20
   INTEGER :: ix,iy,iz,n_int_loc
   REAL(kind=DP) :: qq_fact, qq(3)
   LOGICAL :: exst
   REAL(kind=DP), ALLOCATABLE :: q1(:),q2(:),q3(:)
   REAL(kind=DP) :: qx(3),qy(3),qz(3), qq0,qq1


   write(stdout,*)'BG1', bg(1:3,1)
   write(stdout,*)'BG2', bg(1:3,2)
    write(stdout,*)'BG3', bg(1:3,3)
    if(bg(2,1)==0.d0 .and. bg(3,1)==0.d0 .and.bg(1,2)==0.d0 .and.bg(3,2)==0.d0 .and. bg(1,3)==0.d0 .and.bg(2,3)==0.d0 ) then
       FLUSH(stdout)
       do ig=1,npw
          vg_q(ig)=0.d0
          if(ig==1 .and. gstart==2) then
             n_int_loc=n_int*50
          else
             n_int_loc=n_int
          endif
          allocate(q1(-n_int_loc+1:n_int_loc))
          allocate(q2(-n_int_loc+1:n_int_loc))
          allocate(q3(-n_int_loc+1:n_int_loc))
          do ix=-n_int_loc+1,n_int_loc                                                        
             q1(ix)=(0.5d0*(1.d0/dble(n_int_loc)*(dble(ix-1))+0.5d0/dble(n_int_loc))*bg(1,1)+g(1,ig))**2.d0 
          enddo
          do ix=-n_int_loc+1,n_int_loc
             q2(ix)=(0.5d0*(1.d0/dble(n_int_loc)*(dble(ix-1))+0.5d0/dble(n_int_loc))*bg(2,2)+g(2,ig))**2.d0
          enddo
          do ix=-n_int_loc+1,n_int_loc
             q3(ix)=(0.5d0*(1.d0/dble(n_int_loc)*(dble(ix-1))+0.5d0/dble(n_int_loc))*bg(3,3)+g(3,ig))**2.d0
          enddo
          do ix=-n_int_loc+1,n_int_loc
             qq0=q1(ix)
             do iy=-n_int_loc+1,n_int_loc
                qq1=qq0+q2(iy)
                do iz=-n_int_loc+1,n_int_loc
                   qq_fact=qq1+q3(iz)
                   vg_q(ig)=vg_q(ig)+1.d0/qq_fact
                enddo
             enddo
          enddo
          vg_q(ig)=vg_q(ig)*e2*fpi/(8.d0*(dble(n_int_loc))**3.d0)/tpiba2
          deallocate(q1,q2,q3)
       enddo
    else
       do ig=1,npw
          vg_q(ig)=0.d0
          if(ig==1 .and. gstart==2) then
             n_int_loc=n_int*50
          else
             n_int_loc=n_int
          endif
         

          do ix=-n_int_loc+1,n_int_loc
             do iy=-n_int_loc+1,n_int_loc
                do iz=-n_int_loc+1,n_int_loc
                   
                   qx(:)=0.5d0*(1.d0/dble(n_int_loc)*(dble(ix-1))+0.5d0/dble(n_int_loc))*bg(:,1)
                   qy(:)=0.5d0*(1.d0/dble(n_int_loc)*(dble(iy-1))+0.5d0/dble(n_int_loc))*bg(:,2)
                   qz(:)=0.5d0*(1.d0/dble(n_int_loc)*(dble(iz-1))+0.5d0/dble(n_int_loc))*bg(:,3)

               
                   qq(:)=qx(:)+qy(:)+qz(:)+g(:,ig)
                   qq_fact=qq(1)**2+qq(2)**2+qq(3)**2
                   
                   if(ig==1.and. gstart==2) then
                      vg_q(1)=vg_q(ig)+1.d0/qq_fact
                   else
                      
                      
                       qq_fact=qq(1)**2+qq(2)**2+qq(3)**2
                       vg_q(ig)=vg_q(ig)+1.d0/qq_fact
                                             
                   endif


                   
                enddo
             enddo
          enddo
          vg_q(ig)=vg_q(ig)*e2*fpi/(8.d0*(dble(n_int_loc))**3.d0)/tpiba2
       enddo
    endif
   vg_q(:)=vg_q(:)/omega
   if(gstart==2) write(stdout,*) 'V(G=0) = ',vg_q(1) 
!   if(ionode) then
!      iun = find_free_unit()
!      open( unit= iun, file=trim(tmp_dir)//trim(prefix)//'.vg_q', status='unknown',form='unformatted')
!      write(iun) vg_q(1:npw)
!      close(iun)
!   endif

   iun = find_free_unit()
   CALL diropn( iun, 'vgq', npwx, exst )
   CALL davcio(vg_q,npwx,iun,1,1)
   close(iun)



   return
 end subroutine calculate_vg0


 subroutine  read_vg0
!this subroutine read v(G) for pbc from disk

   USE wannier_gw, ONLY : vg_q
   USE wvfct,    ONLY : npw,npwx
   USE io_global, ONLY : stdout, ionode, ionode_id
   USE io_files,  ONLY : prefix, tmp_dir,diropn
   USE mp,        ONLY : mp_bcast

   implicit none

   INTEGER, EXTERNAL :: find_free_unit

   INTEGER :: iun
   LOGICAL :: exst

!   if(ionode) then
!      iun = find_free_unit()
!      open( unit= iun, file=trim(tmp_dir)//trim(prefix)//'.vg_q', status='old',form='unformatted')
!      read(iun) vg_q(1:npw)
!      close(iun)
!   endif
!   call mp_bcast(vg_q(1:npw),ionode_id)


   iun = find_free_unit()
   CALL diropn( iun, 'vgq', npwx, exst )
   CALL davcio(vg_q,npwx,iun,1,-1)
   close(iun)

   return

 end subroutine read_vg0
