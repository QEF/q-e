! FOR GWW
!
! Author: P. Umari
!
SUBROUTINE wannier_uterms(n_set,l_square,lzero, orthonorm, ecutoff)
!
!this subroutine
!calculates the products <wiwj(r_1)| 1/|r_1-r_2| wi'j'(r_2)>
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
  USE control_flags,        ONLY : gamma_only
  USE cell_base,            ONLY : at, alat, tpiba, omega, tpiba2
  USE wannier_gw
  USE exx,                  ONLY : exx_divergence, exx_grid_init, yukawa
  USE mp,                   ONLY : mp_sum


  implicit none

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


   fac(:)=fac(:)/omega
   if(lzero .and. gstart==2) fac(1)=0.d0
   if(l_square) fac(:)=dsqrt(fac(:))

!open output file
   if(ionode) then
      iunuterms =  find_free_unit()
       open( unit= iunuterms, file=trim(prefix)//'.uterms', status='unknown',form='unformatted')
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
         call mp_sum(umat_tmp(:,:))
         do iw=iw_min,iw_max
            do jw=jw_min,jw_max
               uterms(iw,jw)=2.d0*dble(umat_tmp(iw-iw_min+1,jw-jw_min+1))
               uterms(jw,iw)=uterms(iw,jw)
            enddo
         enddo
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
! #endif
END SUBROUTINE wannier_uterms


SUBROUTINE wannier_uterms_red(n_set,lzero, ecutoff,vmat,dimr,dimc,n_r,n_c,numpw,numpw_all, i_type, l_vsquare)

!this subroutine
!calculates the products <wiwj(r_1)| 1/|r_1-r_2| wi'j'_red(r_2)>
!if required used truncation formula of Onida, PRB 62, 4927 (2000)
! #ifdef __GWW

  USE io_global,            ONLY : stdout, ionode
  USE io_files,             ONLY : find_free_unit, prefix, diropn
  use mp_global,            ONLY : nproc_pool, me_pool
  USE mp,                   ONLY : mp_sum
  USE kinds,    ONLY : DP
  USE gvect
  USE basis
  USE klist
  USE constants, ONLY : e2, pi, tpi, fpi
  USE wvfct,    ONLY : igk, g2kin, npwx, npw, nbnd, nbndx, ecutwfc
  USE control_flags, ONLY: gamma_only
  USE cell_base, ONLY: at, alat, tpiba, omega, tpiba2
  USE wannier_gw
  USE exx, ONLY : exx_divergence, exx_grid_init, yukawa
  USE mp, ONLY : mp_sum


  implicit none

  INTEGER, INTENT(in)  :: n_set  !defines the number of states to be read from disk at the same time
  LOGICAL, INTENT(in)  :: lzero!if true put to zero the G=0,G=0 of v
  REAL(kind=DP), INTENT(in) :: ecutoff!cutoff in Rydberg for g sum
  INTEGER, INTENT(in) :: dimr, dimc!dimension of vmat
  REAL(kind=DP), INTENT(out) :: vmat(dimr,dimc)!matrix to be calculated
!  REAL(kind=DP), INTENT(out) :: vmat(numpw_all,numpw)!matrix to be calculated
  INTEGER, INTENT(in) :: n_r,n_c !periodicity for parallel matrix execution
  INTEGER, INTENT(in) :: numpw!number of reduced products of wanniers
  INTEGER, INTENT(in) :: numpw_all!total number of products of wanniers
  INTEGER, INTENT(in) :: i_type!if == 0 , do products <\tilde_w^P_i|w^P'_j>, if == 1 do products <w^P'_i|w^P'_j>
                               !it requires numpw == numpw_all
  LOGICAL, INTENT(in) :: l_vsquare!if .true. consider v^2 operator

  INTEGER :: iungprod, iungprod_all

  !  --- Internal definitions ---

   COMPLEX(kind=DP), ALLOCATABLE :: tmpspacei(:,:)
   COMPLEX(kind=DP), ALLOCATABLE :: tmpspacej(:,:)
   REAL(kind=DP), ALLOCATABLE ::fac(:)
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
   INTEGER :: icrow,iccol,ilrow,ilcol
#ifdef __SCALAPACK
   INTEGER, EXTERNAL :: indxg2p,indxg2l
#endif

   write(stdout,*) 'Routine wannier_uterms : start'

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
   allocate(tmpspacei(max_ngm,n_set),tmpspacej(max_ngm,n_set),fac(max_ngm))
   allocate(umat_tmp(n_set,n_set))
   if(i_type == 0) then
      iungprod_all = find_free_unit()
      CALL diropn( iungprod_all, 'wiwjwfc', max_ngm*2, exst )
   endif
   iungprod = find_free_unit()
   if(i_type == 0) then
      CALL diropn( iungprod, 'wiwjwfc_red', max_ngm*2, exst )
   else
      CALL diropn( iungprod, 'wiwjwfc_red_1', max_ngm*2, exst )
   endif

   !sets factors terms
 !sets factors terms
!this has already  been called call exx_grid_init()


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


   fac(:)=fac(:)/omega
   if(lzero .and. gstart==2) fac(1)=0.d0

   if(l_vsquare) fac(:)=fac(:)**2.d0

   do iiw=1,ceiling(real(numpw)/real(n_set))
      write(stdout,*) 'uterms iiw', iiw
      do iw=(iiw-1)*n_set+1,min(iiw*n_set,numpw)
         CALL davcio(tmpspacei(:,iw-(iiw-1)*n_set),max_ngm*2,iungprod,iw,-1)

         if(gamma_only .and. gstart == 2) then
            tmpspacei(1,iw-(iiw-1)*n_set)=dble(tmpspacei(1,iw-(iiw-1)*n_set))
         endif
      enddo
      iw_min=(iiw-1)*n_set+1
      iw_max=min(iiw*n_set,numpw)
      do jjw=1,ceiling(real(numpw_all)/real(n_set))
         write(stdout,*) 'uterms jjw', jjw
         do jw=(jjw-1)*n_set+1,min(jjw*n_set,numpw_all)
            if(i_type == 0 ) then
               CALL davcio(tmpspacej(:,jw-(jjw-1)*n_set),max_ngm*2,iungprod_all,jw,-1)
            else
               CALL davcio(tmpspacej(:,jw-(jjw-1)*n_set),max_ngm*2,iungprod,jw,-1)
            endif
            if(gamma_only .and. gstart == 2) then
               tmpspacej(1,jw-(jjw-1)*n_set)=dble(tmpspacej(1,jw-(jjw-1)*n_set))
            endif
         enddo
         jw_min=(jjw-1)*n_set+1
         jw_max=min(jjw*n_set,numpw_all)



!!!!!!!!!!!!!!!!!!!!!!!!!!!
!uses  blas routine

         do jw=1,jw_max-jw_min+1
            tmpspacej(1:ngm_max,jw)= tmpspacej(1:ngm_max,jw)*fac(1:ngm_max)
            if(gstart==2) tmpspacej(1,jw)=0.5d0*tmpspacej(1,jw)
         enddo
         call zgemm('C','N',n_set,n_set,ngm_max,(1.d0,0.d0),tmpspacei,max_ngm,tmpspacej,max_ngm,(0.d0,0.d0),umat_tmp,n_set)
         call mp_sum(umat_tmp(:,:))
         do iw=iw_min,iw_max
            do jw=jw_min,jw_max
               if(.not.l_pmatrix) then
                  vmat(jw,iw)=2.d0*dble(umat_tmp(iw-iw_min+1,jw-jw_min+1))
               else
#ifdef __SCALAPACK
                  icrow=indxg2p(jw,n_r,0,0,nprow)
                  iccol=indxg2p(iw,n_c,0,0,npcol)
                  if(icrow==myrow .and. iccol==mycol) then
                     ilrow=indxg2l(jw,n_r,0,0,nprow)
                     ilcol=indxg2l(iw,n_c,0,0,npcol)
                      vmat(ilrow,ilcol)=2.d0*dble(umat_tmp(iw-iw_min+1,jw-jw_min+1))
                  endif
#endif
               endif
            enddo
         enddo


      enddo
   enddo

   close(iungprod)
   if(i_type == 0) close(iungprod_all)
   deallocate(umat_tmp)
   deallocate(tmpspacei,tmpspacej,fac)
! #endif

END SUBROUTINE wannier_uterms_red
