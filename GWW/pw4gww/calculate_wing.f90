! FOR GWW
! Author: P. Umari
!
subroutine calculate_wing(n_set, orthonorm)

!this subroutine calculate the terms
!\Sum_G <G|\tilde{w^P_i}>\epsilon(G'=0, G, iw)
! it requires the file .e_head
!#ifdef __GWW
 USE io_global,            ONLY : stdout, ionode, ionode_id
 USE io_files,             ONLY : find_free_unit, prefix, diropn
 USE kinds,                ONLY : DP
 USE wannier_gw
 USE mp,                   ONLY : mp_bcast, mp_sum
 USE gvect,                ONLY : mill, ngm, gstart,g
 USE cell_base,            ONLY : tpiba

 implicit none

 INTEGER, INTENT(in)  :: n_set  !defines the number of states to be read from disk at the same time
 INTEGER, INTENT(in)  :: orthonorm!if ==1 opens orthonormalized products of wannier file, if ==2 reduced one


 INTEGER iun, iungprod
 INTEGER :: n_g, ngm_k
 INTEGER :: ig, igg, iw, iiw, i, ii
 REAL(kind=DP) :: omega_g
 REAL(kind=DP), ALLOCATABLE :: freqs(:)
 COMPLEX(kind=DP), ALLOCATABLE :: e_head_g(:,:), e_head(:,:), e_head_g0(:)
 INTEGER, ALLOCATABLE :: mill_k(:,:)
 LOGICAL :: exst
 COMPLEX(kind=DP), ALLOCATABLE :: tmpspacei(:,:)
 REAL(kind=DP), ALLOCATABLE :: wing(:,:), wing_c(:,:)
 REAL(kind=DP) :: sca
 REAL(kind=DP), ALLOCATABLE :: fact(:)
 REAL(kind=DP) :: qq

!read file .e_head

 write(stdout,*) 'Routine calculate_wing'
 call flush_unit(stdout)
 allocate(fact(ngm))

 if(gstart==2) fact(1)=0.d0
 do ig=gstart,ngm
    qq = g(1,ig)**2.d0 + g(2,ig)**2.d0 + g(3,ig)**2.d0
    fact(ig)=1.d0/tpiba/dsqrt(qq)
 end do

 write(stdout,*) 'ATT0.1'
 call flush_unit(stdout)



 if(ionode) then
     iun =  find_free_unit()
     open( unit= iun, file=trim(prefix)//'.e_head', status='old',form='unformatted')
     read(iun) n_g
     read(iun) omega_g
  endif


  call mp_bcast(n_g, ionode_id)
  call mp_bcast(omega_g, ionode_id)
  allocate(freqs(n_g+1))

  if(ionode) then
     read(iun) freqs(1:n_g+1)
     read(iun) ngm_k
  endif

 write(stdout,*) 'ATT0.2'
 call flush_unit(stdout)


  call mp_bcast(freqs(:), ionode_id)
  call mp_bcast(ngm_k, ionode_id)
  allocate(mill_k(3,ngm_k))
  if(ionode) allocate(e_head_g(n_g+1,ngm_k))
  if(ionode) then
     do ig=1,ngm_k
        read(iun) mill_k(1,ig), mill_k(2,ig), mill_k(3,ig)
     enddo
     do ig=1,ngm_k
        read(iun) e_head_g(1:n_g+1,ig)
     enddo
     close(iun)
  endif
  do ig=1,ngm_k
     if(.not. ionode) mill_k(:,ig)= 0
     call mp_sum(mill_k(:,ig))
     !call mp_bcast(mill_k(:,ig), ionode_id)


     !call mp_bcast(e_head_g(:,ig), ionode_id)
  enddo

  write(stdout,*) 'ATT1'
  call flush_unit(stdout)


!create local epsilon vector

  allocate(e_head(ngm, n_g+1))
  e_head(:,:) = (0.d0,0.d0)
  allocate(e_head_g0(ngm_k))

  do ii=1,n_g+1
     if(ionode) then
        e_head_g0(1:ngm_k)=  e_head_g(ii,1:ngm_k)
     else
        e_head_g0(:)=(0.d0,0.d0)
     endif
     call mp_sum(e_head_g0(:))
     do ig=gstart,ngm
        do igg=1,ngm_k
           if ( mill(1,ig)==mill_k(1,igg) .and. &
                mill(2,ig)==mill_k(2,igg) .and. &
                mill(3,ig)==mill_k(3,igg) ) then
              e_head(ig, ii) = e_head_g0(igg)
           endif
        enddo
     enddo
  enddo

  deallocate(e_head_g0)
!loop on n_set groups
  write(stdout,*) 'ATT2'
  call flush_unit(stdout)


  allocate(tmpspacei(max_ngm,n_set))
  iungprod = find_free_unit()

  if(orthonorm==0) then
     CALL diropn( iungprod, 'wiwjwfc', max_ngm*2, exst )
  else if(orthonorm==1) then
     CALL diropn( iungprod, 'wiwjwfc_on', max_ngm*2, exst )
  else
     CALL diropn( iungprod, 'wiwjwfc_red', max_ngm*2, exst )
  endif
  allocate(wing(numw_prod, n_g+1))
  wing(:,:)=0.d0
  allocate(wing_c(numw_prod, n_g+1))
  wing_c(:,:)=0.d0



  do iiw=1,ceiling(real(numw_prod)/real(n_set))
!read states
     do iw=(iiw-1)*n_set+1,min(iiw*n_set,numw_prod)
        CALL davcio(tmpspacei(:,iw-(iiw-1)*n_set),max_ngm*2,iungprod,iw,-1)
     enddo


     write(stdout,*) 'ATT3'
     call flush_unit(stdout)


!loop on states
     do iw=(iiw-1)*n_set+1,min(iiw*n_set,numw_prod)
        do i=1, n_g+1
            sca=0.d0
            do ig=1,max_ngm
               sca=sca+2.d0*real(tmpspacei(ig,iw-(iiw-1)*n_set)*e_head(ig,i))*fact(ig)
            enddo
            call mp_sum(sca)
            wing(iw,i)=sca

            sca=0.d0
            do ig=1,max_ngm
               sca=sca+2.d0*real(tmpspacei(ig,iw-(iiw-1)*n_set)*conjg(e_head(ig,i)))*fact(ig)
            enddo
            call mp_sum(sca)
            wing_c(iw,i)=sca

   !calculate terms
         enddo
          write(stdout,*) 'WING_C:', wing(iw,1)
      enddo
   enddo

   write(stdout,*) 'ATT4'
  call flush_unit(stdout)

!write terms on file

   if(ionode) then
      iun =  find_free_unit()
      open( unit= iun, file=trim(prefix)//'.wing', status='unknown',form='unformatted')
      write(iun) n_g
      write(iun) omega_g
      write(iun) numw_prod
      do i=1,n_g+1
         write(iun) wing(1:numw_prod,i)
      enddo
      do i=1,n_g+1
         write(iun) wing_c(1:numw_prod,i)
      enddo
      close(iun)
  endif

  deallocate(tmpspacei)
  close(iungprod)
  deallocate(fact)
  deallocate(freqs, mill_k)
  if(ionode) deallocate (e_head_g)
  deallocate(e_head)
  deallocate(wing, wing_c)

!#endif
 return
end subroutine calculate_wing
