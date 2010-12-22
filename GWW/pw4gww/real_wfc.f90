! FOR GWW
!
! Author: P. Umari
! Modified by G. Stenuit
!
!-----------------------------------------
subroutine real_wfc( u_trans, ispin, iun_wannier,nbndv)
!----------------------------------------
!
! this routine reads the wavefunctions from iun_wannier
! (GAMMA-ONLY CALCULATIONS) and rotate the wavefunctions
! to real ones
! only ispin states used (not implemented yet)
! it works on two separate subspaces for valence and conduction wfcs

! #ifdef __GWW

  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE us
  USE wvfct,                ONLY : igk, g2kin, npwx, npw, nbnd, nbndx, ecutwfc
  USE gvect
  USE basis
  USE klist
  USE constants,            ONLY : e2, pi, tpi, fpi
  USE io_files,             ONLY : nwordwfc, iunwfc
  USE cell_base,            ONLY : at, alat, tpiba, omega, tpiba2
  USE ions_base,            ONLY : ityp, tau, nat,ntyp => nsp
  USE uspp
  USE uspp_param
  USE becmod,               ONLY : calbec

  implicit none

  COMPLEX(kind=DP) :: u_trans(nbnd,nbnd)!transform unitarian matrix w_i=\sum_j U_{i,j}\Psi_j
  INTEGER, INTENT(in) :: ispin!+1 or -1
  INTEGER, INTENT(in) :: iun_wannier !units for reading wfc
  INTEGER, INTENT(in) :: nbndv!number of first bands wich taken real separately


  !  --- Internal definitions ---

   COMPLEX(kind=DP), ALLOCATABLE  :: evc0(:,:)!reads wavefunctions here
   COMPLEX(kind=DP), ALLOCATABLE  :: evc1(:,:)!reads wavefunctions here
   COMPLEX(kind=DP), ALLOCATABLE  :: evcp(:)!this for the real plus one
   COMPLEX(kind=DP), ALLOCATABLE  :: evcm(:)!this for the real minus one
   INTEGER, ALLOCATABLE :: minusg(:)!table for g --> -g transform
   INTEGER :: i,j,k,ig,jg
   INTEGER :: igk0(npwx)
   INTEGER :: npw0
   REAL(kind=dp) :: g2kin_bp(npwx)
   REAL(kind=DP) :: max_p,max_m,norm, coef
   LOGICAL :: found
   COMPLEX(dp) :: becp0(nkb,nbnd)
   COMPLEX(dp) :: becp1(nkb,nbnd)
   INTEGER :: nt,na, ih, jh, ikb, jkb, ijkb0
   REAL(kind=DP) :: eps,sca
   COMPLEX(kind=DP) :: scaz

   ALLOCATE( evc0(npwx,nbndx))
   ALLOCATE( evc1(npwx,nbndx))
   ALLOCATE( evcp(npwx))
   ALLOCATE( evcm(npwx))
   ALLOCATE( minusg(npwx))

   evc0(:,:) = (0.d0,0.d0)
   evc1(:,:) = (0.d0,0.d0)
   eps=1.d-8

   write(6,*) 'Routine real_wfc: start'

   !reads wfcs from iunwfc

   CALL gk_sort(xk(1,1),ngm,g,ecutwfc/tpiba2, &
              &    npw0,igk0,g2kin_bp)
   CALL davcio(evc0,nwordwfc,iunwfc,1,-1)

   write(6,*) 'Routine real_wfc: out of davcio'


   CALL init_us_2 (npw0,igk0,xk(1,1),vkb)


   !set minusg table


   do ig=1,npw0!loop on plane waves
      found=.false.
      do jg=1,npw0
         if((g(1,igk0(ig))== -g(1,igk0(jg))).and.(g(2,igk0(ig))== -g(2,igk0(jg))) .and. &
              & (g(3,igk0(ig))== -g(3,igk0(jg)))) then
            found=.true.
            minusg(ig)=jg
            exit
         endif
      enddo
      if(.not.found) then
         write(6,*) 'REAL_WFC: CORRESPONDING G NOT FOUND!!'
         stop
      endif
   enddo

!now valence subspace

   do i=1,nbndv!loop on bands

      max_p=0.d0
      max_m=0.d0

      do ig=1,npw0!loop on plane waves
         evcp(ig)=0.5d0*(evc0(ig,i)+conjg(evc0(minusg(ig),i)))
         evcm(ig)=0.5d0*(0.d0,-1.d0)*(evc0(ig,i)-conjg(evc0(minusg(ig),i)))
         if(abs(evcp(ig)) > max_p) max_p = abs(evcp(ig))
         if(abs(evcm(ig)) > max_m) max_m = abs(evcm(ig))
      enddo

      if(max_p > max_m) then
         evc1(1:npw0,i)=evcp(1:npw0)
      else
         evc1(1:npw0,i)=evcm(1:npw0)
      endif
   enddo
   write(6,*) 'Routine real_wfc: do graham'
! now re-orthogonalize Gram-Schmidt

   CALL calbec(npw0, vkb, evc1, becp1, nbnd)

   norm=0.d0
   do ig=1,npw0
      norm=norm+real(conjg(evc1(ig,1))*evc1(ig,1))
   enddo
   ijkb0=0
   do nt = 1, ntyp
     do na = 1, nat
        if (ityp (na) .eq.nt) then
              do jh = 1, nh (nt)
                  jkb = ijkb0 + jh
                    do ih = 1, nh (nt)
                       ikb = ijkb0 + ih
                        norm=norm+real(conjg(becp1(ikb,1))*&
                          &   qq(ih,jh,ityp(na))* becp1(jkb,1))
                     enddo
               enddo
               ijkb0 = ijkb0 + nh (nt)
         endif
       enddo
   enddo


   evc1(1:npw0,1)=evc1(1:npw0,1)/sqrt(norm)

   do i=2,nbndv
      CALL calbec(npw0, vkb, evc1, becp1, nbnd)
      do j=1,i-1
         coef=0.d0
           do ig=1,npw0
              coef=coef+real(conjg(evc1(ig,j))*evc1(ig,i))
           enddo
           ijkb0=0
           do nt = 1, ntyp
              do na = 1, nat
                if (ityp (na) .eq.nt) then
                  do jh = 1, nh (nt)
                    jkb = ijkb0 + jh
                    do ih = 1, nh (nt)
                         ikb = ijkb0 + ih
                          coef=coef+real(conjg(becp1(ikb,j))*&
                            &   qq(ih,jh,ityp(na))* becp1(jkb,i))
                    enddo
                  enddo
                  ijkb0 = ijkb0 + nh (nt)
                endif
              enddo
           enddo


           evc1(1:npw0,i)=evc1(1:npw0,i)-coef*evc1(1:npw0,j)
      enddo


        CALL calbec(npw0, vkb, evc1, becp1, nbnd)
        norm=0.d0
        do ig=1,npw0
           norm=norm+real(conjg(evc1(ig,i))*evc1(ig,i))
        enddo
        ijkb0=0
        do nt = 1, ntyp
          do na = 1, nat
            if (ityp (na) .eq.nt) then
              do jh = 1, nh (nt)
                 jkb = ijkb0 + jh
                 do ih = 1, nh (nt)
                    ikb = ijkb0 + ih
                    norm=norm+real(conjg(becp1(ikb,i))*&
                          &   qq(ih,jh,ityp(na))* becp1(jkb,i))
                 enddo
              enddo
              ijkb0 = ijkb0 + nh (nt)
            endif
          enddo
        enddo

        evc1(1:npw0,i)=evc1(1:npw0,i)/sqrt(norm)
        write(stdout,*) 'Norm :',i,norm!ATTENZIONE
     enddo


!now conduction subspace
 if(nbndv<nbnd) then

 do i=nbndv+1,nbnd!loop on bands

      max_p=0.d0
      max_m=0.d0

      do ig=1,npw0!loop on plane waves
         evcp(ig)=0.5d0*(evc0(ig,i)+conjg(evc0(minusg(ig),i)))
         evcm(ig)=0.5d0*(0.d0,-1.d0)*(evc0(ig,i)-conjg(evc0(minusg(ig),i)))
         if(abs(evcp(ig)) > max_p) max_p = abs(evcp(ig))
         if(abs(evcm(ig)) > max_m) max_m = abs(evcm(ig))
      enddo

      if(max_p > max_m) then
         evc1(1:npw0,i)=evcp(1:npw0)
      else
         evc1(1:npw0,i)=evcm(1:npw0)
      endif
   enddo
   write(6,*) 'Routine real_wfc: do graham'
! now re-orthogonalize Gram-Schmidt

   CALL calbec(npw0, vkb, evc1, becp1, nbnd)

   norm=0.d0
   do ig=1,npw0
      norm=norm+real(conjg(evc1(ig,nbndv+1))*evc1(ig,nbndv+1))
   enddo
   ijkb0=0
   do nt = 1, ntyp
     do na = 1, nat
        if (ityp (na) .eq.nt) then
              do jh = 1, nh (nt)
                  jkb = ijkb0 + jh
                    do ih = 1, nh (nt)
                       ikb = ijkb0 + ih
                        norm=norm+real(conjg(becp1(ikb,nbndv+1))*&
                          &   qq(ih,jh,ityp(na))* becp1(jkb,nbndv+1))
                     enddo
               enddo
               ijkb0 = ijkb0 + nh (nt)
         endif
       enddo
   enddo


   evc1(1:npw0,nbndv+1)=evc1(1:npw0,nbndv+1)/sqrt(norm)

   do i=nbndv+2,nbnd
      CALL calbec(npw0, vkb, evc1, becp1, nbnd)
      do j=nbndv+1,i-1
         coef=0.d0
           do ig=1,npw0
              coef=coef+real(conjg(evc1(ig,j))*evc1(ig,i))
           enddo
           ijkb0=0
           do nt = 1, ntyp
              do na = 1, nat
                if (ityp (na) .eq.nt) then
                  do jh = 1, nh (nt)
                    jkb = ijkb0 + jh
                    do ih = 1, nh (nt)
                         ikb = ijkb0 + ih
                          coef=coef+real(conjg(becp1(ikb,j))*&
                            &   qq(ih,jh,ityp(na))* becp1(jkb,i))
                    enddo
                  enddo
                  ijkb0 = ijkb0 + nh (nt)
                endif
              enddo
           enddo


           evc1(1:npw0,i)=evc1(1:npw0,i)-coef*evc1(1:npw0,j)
      enddo


        CALL calbec(npw0, vkb, evc1, becp1, nbnd)
        norm=0.d0
        do ig=1,npw0
           norm=norm+real(conjg(evc1(ig,i))*evc1(ig,i))
        enddo
        ijkb0=0
        do nt = 1, ntyp
          do na = 1, nat
            if (ityp (na) .eq.nt) then
              do jh = 1, nh (nt)
                 jkb = ijkb0 + jh
                 do ih = 1, nh (nt)
                    ikb = ijkb0 + ih
                    norm=norm+real(conjg(becp1(ikb,i))*&
                          &   qq(ih,jh,ityp(na))* becp1(jkb,i))
                 enddo
              enddo
              ijkb0 = ijkb0 + nh (nt)
            endif
          enddo
        enddo

        evc1(1:npw0,i)=evc1(1:npw0,i)/sqrt(norm)
        write(stdout,*) 'Norm :',i,norm!ATTENZIONE
     enddo





   endif

! calculate unitaria rotation
!U_ij = <\Psi_j|S|w_i>

     CALL calbec(npw0, vkb, evc0, becp0, nbnd)
     CALL calbec(npw0, vkb, evc1, becp1, nbnd)
     u_trans(:,:)=(0.d0,0.d0)
     do i=1,nbnd
        do j=1,nbnd
           coef=0.d0
           do ig=1,npw0
              u_trans(i,j)=u_trans(i,j)+conjg(evc0(ig,j))*evc1(ig,i)
           enddo
           ijkb0=0
           do nt = 1, ntyp
              do na = 1, nat
                if (ityp (na) .eq.nt) then
                  do jh = 1, nh (nt)
                    jkb = ijkb0 + jh
                    do ih = 1, nh (nt)
                         ikb = ijkb0 + ih
                          u_trans(i,j)=u_trans(i,j)+conjg(becp0(ikb,j))*&
                            &   qq(ih,jh,ityp(na))* becp1(jkb,i)
                    enddo
                  enddo
                  ijkb0 = ijkb0 + nh (nt)
                endif
              enddo
           enddo
        enddo
     enddo
!check that u_trans is unitarian
!test u_trans:
  do i=1,nbnd
    !do j=1,nbnd
     j=i
      scaz=(0.d0,0.d0)
        do k=1,nbnd
           scaz=scaz+conjg(u_trans(k,i))*u_trans(k,j)
        enddo
        write(stdout,*) 'Test u_trans:', i,j,scaz!ATTENZIONE
     !enddo
  enddo



! re-writes on file iun_wannier

     CALL davcio(evc1,nwordwfc,iun_wannier,1,1)


     DEALLOCATE(evc0)
     DEALLOCATE(evc1)
     DEALLOCATE(evcp)
     DEALLOCATE(evcm)
     DEALLOCATE(minusg)

! #endif __GWW

     return
   end subroutine real_wfc
