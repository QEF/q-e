! FOR GWW
!
! Author: P. Umari
!
SUBROUTINE wannier_uterms_c(n_set, lzero, orthonorm, ecutoff)

!this subroutine
!calculates the products <wiwj^'(r_1)| 1/|r_1-r_2| wi'j'(r_2)>
!if required used truncation formula of Onida, PRB 62, 4927 (2000)
! #ifdef __GWW

  USE io_global,            ONLY : stdout, ionode
  USE io_files,             ONLY : find_free_unit, prefix, diropn
  use mp_global,            ONLY : nproc_pool, me_pool
  USE kinds,                ONLY : DP
  USE gvect
  USE basis
  USE klist
  USE constants,            ONLY : e2, pi, tpi, fpi
  USE wvfct,                ONLY : igk, g2kin, npwx, npw, nbnd, nbndx, ecutwfc
  USE control_flags,        ONLY: gamma_only
  USE cell_base,            ONLY: at, alat, tpiba, omega, tpiba2
  USE wannier_gw
  USE exx,                  ONLY : exx_divergence, exx_grid_init, yukawa
  USE mp,        ONLY : mp_sum

  implicit none

  INTEGER, INTENT(in)  :: n_set  !defines the number of states to be read from disk at the same time
  LOGICAL, INTENT(in)  :: lzero !if true put the term G=0,G=0 of v to zero
  INTEGER, INTENT(in)  :: orthonorm!if ==1 opens orthonormalized products, if ==2 reduced ones
  REAL(kind=DP), INTENT(in) :: ecutoff!cutoff in Rydberg for g sum

  INTEGER :: iungprod, iungprodprim,iunuterms

  !  --- Internal definitions ---

   COMPLEX(kind=DP), ALLOCATABLE :: tmpspacei(:,:)
   COMPLEX(kind=DP), ALLOCATABLE :: tmpspacej(:,:)
   REAL(kind=DP), ALLOCATABLE ::fac(:)
   REAL(kind=DP), ALLOCATABLE :: uterms(:,:)!temporary reading array
   INTEGER :: iw,jw, iiw,jjw
   INTEGER :: ig
   LOGICAL :: exst
   REAL (kind=DP) :: qq
   INTEGER :: igk0(npwx)
   REAL(kind=dp) :: g2kin_bp(npwx)
   INTEGER :: npw0
   REAL(kind=DP) :: exxdiv
!   REAL(kind=DP) :: ngm_max
   INTEGER :: ngm_max
! change this to integer ????
   write(stdout,*) 'Routine wannier_uterms_ : start'

!determine ngm_max
   ngm_max=0
   do ig=1,ngm
      if(gg(ig)*tpiba2 >= ecutoff) exit
      ngm_max=ngm_max+1
   enddo

   write(stdout,*) 'NGM MAX:', ngm_max, ngm




! reads wfcs from iunwfc

   CALL gk_sort(xk(1,1),ngm,g,ecutwfc/tpiba2, &
              &    npw0,igk0,g2kin_bp)
   allocate(uterms(numw_prod_c,numw_prod))
   if(.not.lsmallgrid) then
     allocate(tmpspacei(max_ngm,n_set),tmpspacej(max_ngm,n_set),fac(max_ngm))
   else
     allocate(tmpspacei(npw0,n_set),tmpspacej(npw0,n_set),fac(npw0))
   endif
   iungprod = find_free_unit()
   if(.not.lsmallgrid) then
      if(orthonorm==0) then
         CALL diropn( iungprod, 'wiwjwfc', max_ngm*2, exst )
      else if(orthonorm==1) then
         CALL diropn( iungprod, 'wiwjwfc_on', max_ngm*2, exst )
      else
         CALL diropn( iungprod, 'wiwjwfc_red', max_ngm*2, exst )
      endif
   else
     CALL diropn( iungprod, 'wiwjwfc', npw0*2, exst )
   endif

   iungprodprim = find_free_unit()
   if(.not.lsmallgrid) then
     CALL diropn( iungprodprim, 'wiwjwfc_prim', max_ngm*2, exst )
   else
     CALL diropn( iungprodprim, 'wiwjwfc_prim', npw0*2, exst )
   endif


!sets factors terms
!this has already  been called   call exx_grid_init()
   exxdiv=exx_divergence()


   if(.not.lsmallgrid) then
     do ig=1,max_ngm
         qq = g(1,ig)**2.d0 + g(2,ig)**2.d0 + g(3,ig)**2.d0

         if(.not.l_truncated_coulomb) then

            if (qq > 1.d-8) then
               fac(ig)=e2*fpi/(tpiba2*qq + yukawa )
            else
               fac(ig)= - exxdiv ! & ! or rather something else (see F.Gygi)
 !                         - e2*fpi   ! THIS ONLY APPLYS TO HYDROGEN
               if (yukawa .gt. 1.d-8 ) then
                  fac(ig) = fac(ig) + e2*fpi/(tpiba2*qq + yukawa )
               end if
            end if
         else
            if (qq > 1.d-8) then
               fac(ig)=(e2*fpi/(tpiba2*qq))*(1.d0-dcos(dsqrt(qq)*truncation_radius*tpiba))
            else
               fac(ig)=e2*fpi*(truncation_radius**2.d0/2.d0)
            endif
         endif

      end do
    else
      do ig=1,npw0
         qq = g(1,igk0(ig))**2.d0 + g(2,igk0(ig))**2.d0 + g(3,igk0(ig))**2.d0


           if (qq > 1.d-8) then
              fac(ig)=e2*fpi/(tpiba2*qq + yukawa )
          else
             fac(ig)= - exxdiv ! & ! or rather something else (see F.Gygi)
 !                         - e2*fpi   ! THIS ONLY APPLYS TO HYDROGEN
             if (yukawa .gt. 1.d-8 ) then
                 fac(ig) = fac(ig) + e2*fpi/(tpiba2*qq + yukawa )
             end if
            end if
      end do
    endif
    fac(:)=fac(:)/omega

    if(lzero .and. gstart == 2) fac(1)=0.d0

!open output file
    if(ionode) then
       iunuterms =  find_free_unit()
       if(.not.lzero) then
          open( unit= iunuterms, file=trim(prefix)//'.uterms_prim', status='unknown',form='unformatted')
       else
          open( unit= iunuterms, file=trim(prefix)//'.uterms_prim_zero', status='unknown',form='unformatted')
       endif
    endif

   uterms(:,:)=0.d0
   do iiw=1,ceiling(real(numw_prod_c)/real(n_set))
      do iw=(iiw-1)*n_set+1,min(iiw*n_set,numw_prod_c)
        if(.not.lsmallgrid) then
           CALL davcio(tmpspacei(:,iw-(iiw-1)*n_set),max_ngm*2,iungprodprim,iw,-1)
        else
           CALL davcio(tmpspacei(:,iw-(iiw-1)*n_set),npw0*2,iungprodprim,iw,-1)
        endif
        if(gamma_only .and. gstart == 2) then
           tmpspacei(1,iw-(iiw-1)*n_set)=dble(tmpspacei(1,iw-(iiw-1)*n_set))
        endif
      enddo
      do jjw=1,ceiling(real(numw_prod)/real(n_set))
         do jw=(jjw-1)*n_set+1,min(jjw*n_set,numw_prod)
          if(.not.lsmallgrid) then
             CALL davcio(tmpspacej(:,jw-(jjw-1)*n_set),max_ngm*2,iungprod,jw,-1)
          else
             CALL davcio(tmpspacej(:,jw-(jjw-1)*n_set),npw0*2,iungprod,jw,-1)
          endif
          if(gamma_only .and. gstart == 2) then
           tmpspacei(1,jw-(jjw-1)*n_set)=dble(tmpspacei(1,jw-(jjw-1)*n_set))
        endif
        enddo

        do iw=(iiw-1)*n_set+1,min(iiw*n_set,numw_prod_c)
           do jw=(jjw-1)*n_set+1,min(jjw*n_set,numw_prod)
             uterms(iw,jw)=0.d0
             if(.not.gamma_only) then
               if(.not.lsmallgrid) then
                 do ig=1,ngm_max
                   uterms(iw,jw)=uterms(iw,jw) + dble(fac(ig)*&
               &conjg(tmpspacei(ig,iw-(iiw-1)*n_set))*tmpspacej(ig,jw-(jjw-1)*n_set))
                 enddo
               else
                 do ig=1,npw0
                   uterms(iw,jw)=uterms(iw,jw) + dble(fac(ig)*&
           &conjg(tmpspacei(ig,iw-(iiw-1)*n_set))*tmpspacej(ig,jw-(jjw-1)*n_set))
                 enddo
               endif
             else
             if(.not.lsmallgrid) then
               do ig=1,ngm_max
                 uterms(iw,jw)=uterms(iw,jw) + 2.d0*dble(fac(ig)*&
           &conjg(tmpspacei(ig,iw-(iiw-1)*n_set))*tmpspacej(ig,jw-(jjw-1)*n_set))
                enddo
                if(gstart==2) then
                   uterms(iw,jw)=uterms(iw,jw)-dble(fac(1)*&
           &conjg(tmpspacei(1,iw-(iiw-1)*n_set))*tmpspacej(1,jw-(jjw-1)*n_set))
                endif
             else
              do ig=1,npw0
                 uterms(iw,jw)=uterms(iw,jw) + 2.d0*dble(fac(ig)*&
           &conjg(tmpspacei(ig,iw-(iiw-1)*n_set))*tmpspacej(ig,jw-(jjw-1)*n_set))
              enddo
              if(gstart==2) then
                 uterms(iw,jw)=uterms(iw,jw)-dble(fac(1)*&
           &conjg(tmpspacei(1,iw-(iiw-1)*n_set))*tmpspacej(1,jw-(jjw-1)*n_set))
              endif
           endif
          endif
          call mp_sum(uterms(iw,jw))
        enddo
      enddo
    enddo
   enddo
   if(ionode) then
      write(iunuterms) numw_prod_c
      write(iunuterms) numw_prod
      do iw=1,numw_prod_c
         write(iunuterms) uterms(iw,1:numw_prod)
      enddo
      close(iunuterms)
   endif
   close(iungprod)
   close(iungprodprim)
   deallocate(tmpspacei,tmpspacej,fac,uterms)
! #endif
END SUBROUTINE wannier_uterms_c

