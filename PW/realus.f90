module realus
  USE kinds,                ONLY : DP
  integer,allocatable :: box(:,:), maxbox(:)
  REAL(KIND=DP),ALLOCATABLE ::&
       qsave(:,:,:,:) !To be used in newd in real space - max kkbeta, nbrx,nbrx
  !nat 
  real(KIND=DP),allocatable:: boxradius(:),boxdistance(:,:),xyz(:,:,:)     
  real(KIND=DP),allocatable:: spher(:,:,:) !Spherical harmonics
  logical ::tqr
contains
  subroutine deallocatenewdreal()
    if (allocated(box))  deallocate (box)
    if (allocated(boxdistance)) deallocate (boxdistance)
    if (allocated(maxbox))  deallocate (maxbox)
    if (allocated(qsave)) deallocate (qsave)
    if (allocated(boxradius)) deallocate (boxradius)
    if (allocated(xyz)) deallocate (xyz)
    if (allocated(spher)) deallocate (spher)
  end subroutine deallocatenewdreal


#include "f_defs.h"
  !
  !----------------------------------------------------------------------
  subroutine qpointlist
    !----------------------------------------------------------------------
    !This subroutine is the driver routine of the box system in this implementation of US in real space.
    !All the variables common in the module are computed and stored for reusing. 
    !This routine has to be called every time the atoms are moved and of course in the beginning.
    !A set of spherical boxes are computed for each atom. 
    !In boxradius there are the radii of the boxes.
    !In maxbox the upper limit of leading index id est the number of points of the fine mesh contained in each box.
    !In xyz there are the coordinates of the points with origin in the centre of atom.
    !In boxdistance the distance from the centre.
    !In spher the spherical harmonics computed for each box
    !In qsave the q value interpolated in these boxes.

    !Strictly speaking only qsave, box and maxbox are necessary and so we can also save memory deallocating the others.
    !Most of time is spent here but the others are faster. Actually 10^-3 Rydberg accuracy is reached. 

    USE ions_base,            ONLY : nat, ntyp => nsp, ityp, tau
    USE cell_base,            ONLY : omega,alat
    USE gvect,                ONLY : nr1, nr2, nr3, nrx1, nrx2, nrx3, g, gg, &
         ngm, gstart, ig1, ig2, ig3, eigts1, eigts2, &
         eigts3, nl,nrxx
    USE lsda_mod,             ONLY : nspin
    USE scf,                  ONLY : vr, vltot
    USE uspp,                 ONLY : okvan
    USE uspp,                 ONLY : indv,deeq, dvan,nhtol, nhtolm,aainit,ap,nhtoj,lpx, lpl
    USE uspp_param,           ONLY : lmaxq, nh, nhm, tvanp,kkbeta,nbeta,qfunc,dion,lmaxkb,qfcoef,nqf, nqlc, lll, rinner,jjj
    USE wvfct,                ONLY : gamma_only
    USE wavefunctions_module, ONLY : psic
    !

    USE atom,                 ONLY : r
    USE io_global,  ONLY : stdout
    USE cell_base,  ONLY : at, bg
    USE parameters, ONLY : nbrx
    USE pfft,       ONLY : npp
    implicit none

    integer :: inat, indm, inbrx1, inbrx2, inrxx, idimension, ilm,ih,jh,is,m, iih,ijh, ilemme
    integer :: roughestimate,goodestimate, lamx2,l,nt,ndm
    integer,allocatable :: bufferpoints (:,:)
    real (kind=DP),allocatable :: bufferdistance (:,:),cost(:),bufferxyz(:,:,:)
    real (kind=DP) :: tempdistance,interpolatedqtot,first,second


    ! from make_pointlist  
    integer  :: index0,index,indproc,iat,ir,iat1
    integer  :: i,j,k,i0,j0,k0,ipol,ishift(3),lemme,nb,mb,ilast
    integer :: i1,i2,i3,i4
    real(kind=dp) :: a1,a2,a3, a4, rx

    real(kind=dp) :: posi(3),distance,shift(3),scalprod, distmin
    real(kind=dp), allocatable :: tau0(:,:), rl(:,:), rl2(:)
    real(KIND=DP),allocatable:: tempspher(:,:),qtot(:,:,:),xsp(:),ysp(:),wsp(:)
    real(kind=DP) :: pi, fpi,eps
    parameter (pi = 3.14159265358979d0, fpi = 4.d0 * pi,eps = 1.0d-9)


    ! end from 

    call start_clock ('qpointlist')


    !generates the spherical harmonics for r space
    !


    if (.not.allocated (tau0))    allocate(tau0(3,nat)) 
    if (allocated(box)) deallocate (box)
    if (.not.allocated(maxbox)) allocate (maxbox(nat))
    maxbox=0
    !This finds the radii for integration
    if (.not.allocated(boxradius)) then
       write (*,*) 'Genero boxradius' 
       allocate (boxradius(ntyp))
       boxradius=0.D0
       do inat=1,ntyp
          do inbrx1 =1,nbeta(inat)
             do inbrx2=1,nbeta(inat)
                do indm=kkbeta(inat),1,-1
                   if ((abs(qfunc(indm,inbrx1,inbrx2,inat))>10**(-6))) then
                      boxradius(inat)=max(r(indm,inat),boxradius(inat))
                      write (*,*) 'radius for ',inat,' = ', boxradius(inat)
                      exit
                   endif
                end do
             end do
          end do
       end do
    endif

!We first save points in a buffer and then we choose a better fitted table... if fortran had vectors 
!like c++ or free form arrays like java we could save up to 40 % of memory

    roughestimate=10*int(((1/(omega))*nrxx*(maxval(boxradius)**3)))
    allocate (bufferpoints(roughestimate,nat))
    bufferpoints=0
    allocate (bufferdistance(roughestimate,nat))
    bufferdistance=0D0
    allocate (bufferxyz(3,roughestimate,nat))

    !FROM SOUBROUTINE MAKE_POINTLIST we have copied these lines and adapted to our needs
    !The mesh is multiplied by 27 to cover also border atoms on properly.


    index0 = 0
#ifdef __PARA
    do indproc=1,me-1
       index0 = index0 + nrx1*nrx2*npp(indproc)
    enddo

#endif
    ! Bring all the atomic positions on the first unit cell



    tau0=tau
    call cryst_to_cart(nat,tau0,bg,-1)
    do iat=1,nat
       do ipol=1,3
          tau0(ipol,iat)=tau0(ipol,iat)-nint(tau0(ipol,iat))
       enddo
    enddo
    call cryst_to_cart(nat,tau0,at,1)

    !now we find the points
    maxbox=0
    do iat = 1,nat



       do ir=1,nrxx

          index = index0 + ir - 1

          k0 = index/(nrx1*nrx2)
          index = index - (nrx1*nrx2) * k0
          j0 = index / nrx1
          index = index - nrx1*j0
          i0 = index

          do i = i0-nr1,i0+nr1, nr1
             do j = j0-nr2, j0+nr2, nr2
                do k = k0-nr3, k0+nr3, nr3

                   do ipol=1,3
                      posi(ipol) =  real(i)/real(nr1) * at(ipol,1) &
                           + real(j)/real(nr2) * at(ipol,2) &
                           + real(k)/real(nr3) * at(ipol,3)

                   enddo

                   tempdistance = sqrt((posi(1)-tau0(1,iat))**2+&
                        &(posi(2)-tau0(2,iat))**2+(posi(3)-tau0(3,iat))**2)
                   tempdistance = tempdistance*alat
                   if (tempdistance<boxradius(ityp(iat)))then
                      maxbox(iat)=maxbox(iat)+1
                      bufferpoints(maxbox(iat),iat)=ir
                      bufferdistance(maxbox(iat), iat)=tempdistance
                      bufferxyz(1,maxbox(iat),iat)=(posi(1)-tau0(1,iat))*alat
                      bufferxyz(2,maxbox(iat),iat)=(posi(2)-tau0(2,iat))*alat
                      bufferxyz(3,maxbox(iat),iat)=(posi(3)-tau0(3,iat))*alat
                   endif
                end do

             end do

          end do
       end do
       write (*,*) 'saved ', maxbox (iat), ' of ', nrxx, ' for ',iat
    end do

    goodestimate=maxval(maxbox)
    if (goodestimate > roughestimate) call errore('qpointlist', 'Roughestimate &
         &is too rough',2)
    !now store them in a more convenient place

    if (allocated(box)) deallocate (box)
    if (allocated(boxdistance)) deallocate (boxdistance)
    if (allocated(xyz)) deallocate(xyz)
    allocate (box(goodestimate,nat), &
         boxdistance(goodestimate,nat),xyz(3,goodestimate,nat))
    box=0
    boxdistance=0
    do inat=1,nat
       do indm=1,goodestimate
          box(indm,inat)=bufferpoints(indm,inat)
          boxdistance(indm,inat)=&
               &   bufferdistance(indm,inat)
          xyz(1,indm,inat)=bufferxyz(1,indm,inat)
          xyz(2,indm,inat)=bufferxyz(2,indm,inat)
          xyz(3,indm,inat)=bufferxyz(3,indm,inat)       
       end do
    end do
    deallocate (bufferpoints)
    deallocate (bufferdistance)
    deallocate (bufferxyz)

! Now it frees memory used for buffers.

!Now it computes the spherical harmonics
    lamx2=lmaxq**2
    if (allocated(spher)) deallocate(spher)
    allocate(spher(goodestimate,lamx2,nat))
    do inat=1,nat

       idimension=maxbox(inat)
       allocate (rl(3,idimension),rl2(idimension))
       do ir=1,idimension



          do ipol=1,3
             rl(ipol,ir) = xyz(ipol,ir,inat)

          enddo
          rl2(ir)=rl(1,ir)**2+rl(2,ir)**2+rl(3,ir)**2
       end do

       allocate (tempspher(idimension,lamx2))
       call ylmr2 (lamx2,idimension,rl,rl2,tempspher)

!We do a similar buffering
       do ir=1,idimension
          do ilm=1,lamx2
             spher(ir,ilm,inat)=tempspher(ir,ilm)
          end do
       end do

!We delete the buffer
       deallocate(rl,rl2,tempspher)



    end do




 !Let's do the main work
    if (allocated(qsave)) deallocate(qsave)
    allocate(qsave(goodestimate,maxval(nh),maxval(nh),nat))
    qsave=0d0

!The source is inspired by init_us_1

!We perform two steps: first we compute for each l the qtot (radial q), then we interpolate it in our mesh, 
!and then we add it to qsave with the correct spherica harmonics

!Q is read from pseudo and it is divided into two parts:
!in the inner radius a polinomial representation is known and so strictly speaking we do not use interpolation but just 
!compute the correct value
!


    ndm = MAXVAL (kkbeta(1:ntyp))
    ap (:,:,:)   = 0.d0
    do nt=1,ntyp
       nqlc(nt) = MIN ( nqlc(nt), lmaxq )
       if ( nqlc(nt) < 0 )  nqlc(nt) = 0
    end do
    do nt = 1, ntyp
       ih = 1
       do nb = 1, nbeta (nt)
          l = lll (nb, nt)
          j = jjj (nb, nt)
          do m = 1, 2 * l + 1
             nhtol (ih, nt) = l
             nhtolm(ih, nt) = l*l+m
             nhtoj (ih, nt) = j
             indv  (ih, nt) = nb
             ih = ih + 1
          enddo
       enddo
       do ih = 1, nh (nt)
          do jh = 1, nh (nt)
             if (nhtol (ih, nt) == nhtol (jh, nt) .and. &
                  nhtolm(ih, nt) == nhtolm(jh, nt) ) then
                ir = indv (ih, nt)
                is = indv (jh, nt)
                dvan (ih, jh, nt) = dion (ir, is, nt)
             endif
          enddo
       enddo
    enddo
    !
    !  compute Clebsch-Gordan coefficients
    !
    if (okvan) call aainit (lmaxkb + 1)
    do inat = 1, nat
       nt =ityp(inat)
       if (allocated(qtot)) deallocate (qtot)
       allocate(qtot(kkbeta(nt),nbeta(nt),nbeta(nt)))
      
!variables used for spline interpolation
 if (allocated(xsp)) deallocate(xsp)
       if (allocated(ysp)) deallocate(ysp)
       if (allocated(wsp)) deallocate(wsp)

       allocate (xsp(kkbeta(nt)),ysp(kkbeta(nt)),wsp(kkbeta(nt)))
!The radii in x
       do ir=1,kkbeta(nt)
          xsp(ir)=r(ir,nt)
       end do

       if (tvanp (nt) ) then
          do l = 0, nqlc (nt) - 1
             !
             !     first we build for each nb,mb,l the total Q(|r|) function
             !     note that l is the true (combined) angular momentum
             !     and that the arrays have dimensions 1..l+1
             !

             do nb = 1, nbeta (nt)
                do mb = nb, nbeta (nt)
                   if ( (l >= abs (lll (nb, nt) - lll (mb, nt) ) ) .and. &
                        (l <= lll (nb, nt) + lll (mb, nt) )        .and. &
                        (mod (l + lll (nb, nt) + lll (mb, nt), 2) == 0) ) then
                      do ir = 1, kkbeta (nt)
                         if (r (ir, nt) >= rinner (l + 1, nt) ) then
                            qtot (ir, nb, mb) = qfunc (ir, nb, mb, nt)/(r(ir,nt)**2)
                         else
                            ilast = ir
                         endif
                      enddo
                      if (rinner (l + 1, nt) > 0.d0) &
                           call setqfcorr(qfcoef (1, l+1, nb, mb, nt), &
                           qtot(1,nb,mb), r(1,nt), nqf(nt),l,ilast)
!We save the values in y
                      do ir =1,kkbeta(nt)
                         ysp(ir)=qtot(ir,nb,mb)
                      end do

!compute the first derivative in first point
                      call setqfcorrpointfirst(qfcoef (1, l+1, nb, mb, nt), &
                           first, r(1,nt), nqf(nt),l)
!compute the second derivative in second point
                      call setqfcorrpointsecond(qfcoef (1, l+1, nb, mb, nt), &
                           second, r(1,nt), nqf(nt),l)
!call spline 
                      call spline(xsp,ysp,wsp,second,first)


                      do ir=1,maxbox(inat)
                         if (boxdistance (ir, inat) < rinner (l + 1, nt) ) then
!if in the inner radius just compute the polynomial
                            call setqfcorrpoint(qfcoef (1, l+1, nb, mb, nt), &
                                 interpolatedqtot, boxdistance(ir,inat), nqf(nt),l)
                         else   
!else spline interpolation
                            interpolatedqtot= splint(xsp,ysp,wsp,boxdistance(ir,inat))
                         endif
                         do iih=1,nh(nt)
                            do ijh=iih,nh(nt)
                               if ((nb.eq.indv(iih,nt)).and.(mb.eq.indv(ijh,nt))) then
                                  do lemme=l*l+1,(l+1)*(l+1)
                                     qsave(ir,iih,ijh,inat)=qsave(ir,iih,ijh,inat)+interpolatedqtot* &
                                          &ap(lemme,nhtolm(iih,nt),nhtolm(ijh,nt)) * spher(ir,lemme,inat)
                                  end do
                               endif
                            end do
                         end do
                      enddo
                   endif

                end do
             end do
          end do
       endif
    end do
    if (allocated(qtot)) deallocate (qtot)
    if (allocated(xsp)) deallocate(xsp)
    if (allocated(ysp)) deallocate(ysp)
    if (allocated(wsp)) deallocate(wsp)





    call stop_clock('qpointlist')
  end subroutine qpointlist

  subroutine newdrealsub

!This subroutine is the real version of newd
!This is faster

    USE kinds,                ONLY : DP
    USE ions_base,            ONLY : nat, ntyp => nsp, ityp
    USE cell_base,            ONLY : omega
    USE gvect,                ONLY : nr1, nr2, nr3, nrx1, nrx2, nrx3, g, gg, &
         ngm, gstart, ig1, ig2, ig3, eigts1, eigts2, &
         eigts3, nl, nrxx
    USE lsda_mod,             ONLY : nspin
    USE scf,                  ONLY : vr, vltot
    USE uspp,                 ONLY : okvan
    USE uspp,                 ONLY : deeq, dvan, nhtolm, nhtol
    USE uspp_param,           ONLY : lmaxq, nh, nhm, tvanp
    USE wvfct,                ONLY : gamma_only
    USE wavefunctions_module, ONLY : psic
    implicit none
    !
    real(kind=DP) :: pi, fpi
    parameter (pi = 3.14159265358979d0, fpi = 4.d0 * pi)
    integer :: na,ih,jh,is,inat,inrxx,nt, lmi1,lmi2,lmj1,lmj2

    if (.not.okvan) then
       ! no ultrasoft potentials: use bare coefficients for projectors
       do na = 1, nat
          nt = ityp (na)
          do is = 1, nspin
             do ih = 1, nh (nt)
                do jh = ih, nh (nt)
                   deeq (ih, jh, na, is) = dvan (ih, jh, nt)
                   deeq (jh, ih, na, is) = deeq (ih, jh, na, is)
                enddo
             enddo
          enddo
       end do
    endif
    call start_clock('newd')
    deeq(:,:,:,:)=0.D0 !annulliamo deeq

    do inat=1,nat

       do ih=1,nh(ityp(inat))

          do jh=ih,nh(ityp(inat))

             do is=1,nspin
                do inrxx=1,maxbox(inat)




                   deeq (ih,jh,inat,is)=deeq(ih,jh,inat,is)+((qsave(inrxx,ih,jh,inat) &
                        & ) * &
                        &(vltot(box(inrxx,inat))+&
                        &vr(box(inrxx,inat),is ))*(omega/nrxx) )  

                end do
             end do
          end do
       end do
    end do

!!$    write (55,*) 'scriviamo deeq'
!!$    do inat=1,nat
!!$       do ih=1,nh(ityp(inat))
!!$          do jh=ih,nh(ityp(inat))
!!$             do is=1,nspin
!!$                write (55,*) inat,ih,jh,deeq(ih,jh,inat,is)
!!$             end do
!!$          end do
!!$       end do
!!$    end do
    !stop

#ifdef __PARA
    call reduce (nhm * nhm * nat * nspin, deeq)
#endif
    do na = 1, nat
       nt = ityp (na)
       do is = 1, nspin
          !           WRITE( stdout,'( "dmatrix atom ",i4, " spin",i4)') na,is
          !           do ih = 1, nh(nt)
          !              WRITE( stdout,'(8f9.4)') (deeq(ih,jh,na,is),jh=1,nh(nt))
          !           end do
          do ih = 1, nh (nt)
             do jh = ih, nh (nt)
                deeq (ih, jh, na, is) = deeq (ih, jh, na, is) + dvan (ih,jh,nt)
                deeq (jh, ih, na, is) = deeq (ih, jh, na, is)
             enddo
          enddo
       enddo
       !        WRITE( stdout,'( "dion pseudo ",i4)') nt
       !        do ih = 1, nh(nt)
       !           WRITE( stdout,'(8f9.4)') (dvan(ih,jh,nt),jh=1,nh(nt))
       !        end do

    enddo
    call stop_clock('newd')
  end subroutine newdrealsub
  subroutine setqfcorr (qfcoef, rho, r, nqf, ltot, mesh)
    !-----------------------------------------------------------------------
    !
    !   This routine compute the first part of the Q function up to rinner.
    !   On output it contains  Q
    !
    !
    USE kinds
    implicit none
    !
    !     first the dummy variables
    !
    integer :: nqf, ltot, mesh
    ! input: the number of coefficients
    ! input: the angular momentum
    ! input: the number of mesh point
    real(kind=DP) :: r (mesh), qfcoef (nqf), rho (mesh)
    ! input: the radial mesh
    ! input: the coefficients of Q
    ! output: the function to be computed
    !
    !     here the local variables
    !
    integer :: ir, i
    ! counter on  mesh points
    ! counter on the coeffients

    real(kind=DP) :: rr
    ! the square of the radius
 
    do ir = 1, mesh
       rr = r (ir) **2
       rho (ir) = qfcoef (1)
       do i = 2, nqf
          rho (ir) = rho (ir) + qfcoef (i) * rr** (i - 1)
       enddo
       rho (ir) = rho (ir) * r (ir) ** (ltot)

    enddo
    return
  end subroutine setqfcorr
subroutine setqfcorrpoint (qfcoef, rho, r, nqf, ltot)
    !-----------------------------------------------------------------------
    !
    !   This routine compute the first part of the Q function in the point r.
    !   On output it contains  Q
    !
    !
    USE kinds
    implicit none
    !
    !     first the dummy variables
    !
    integer :: nqf, ltot, mesh
    ! input: the number of coefficients
    ! input: the angular momentum
    ! input: the number of mesh point
    real(kind=DP) :: r , qfcoef (nqf), rho 
    ! input: the radial mesh
    ! input: the coefficients of Q
    ! output: the function to be computed
    !
    !     here the local variables
    !
    integer :: ir, i
    ! counter on  mesh points
    ! counter on the coeffients

    real(kind=DP) :: rr
    ! the square of the radius
 

       rr = r  **2
       rho  = qfcoef (1)
       do i = 2, nqf
          rho  = rho  + qfcoef (i) * rr** (i - 1)
       enddo
       rho  = rho  * r  ** (ltot)


    return
  end subroutine setqfcorrpoint
subroutine setqfcorrpointfirst (qfcoef, rho, r, nqf, ltot)
    !-----------------------------------------------------------------------
    !
    !   
    !   On output it contains  Q'
    !   probably wrong
    !
    USE kinds
    implicit none
    !
    !     first the dummy variables
    !
    integer :: nqf, ltot, mesh
    ! input: the number of coefficients
    ! input: the angular momentum
    ! input: the number of mesh point
    real(kind=DP) :: r , qfcoef (nqf), rho 
    ! input: the radial mesh
    ! input: the coefficients of Q
    ! output: the function to be computed
    !
    !     here the local variables
    !
    integer :: ir, i
    ! counter on  mesh points
    ! counter on the coeffients

    real(kind=DP) :: rr
    ! the square of the radius
 

       rr = r  **2

      do i = max(1,2-ltot), nqf
          rho  = rho  + qfcoef (i) * rr** (i - 2+ltot)*(i-1+ltot)
       enddo




    return
  end subroutine setqfcorrpointfirst
  subroutine setqfcorrpointsecond (qfcoef, rho, r, nqf, ltot)
    !-----------------------------------------------------------------------
    !
    !   
    !   On output it contains  Q''
    !
    !
    USE kinds
    implicit none
    !
    !     first the dummy variables
    !
    integer :: nqf, ltot, mesh
    ! input: the number of coefficients
    ! input: the angular momentum
    ! input: the number of mesh point
    real(kind=DP) :: r , qfcoef (nqf), rho 
    ! input: the radial mesh
    ! input: the coefficients of Q
    ! output: the function to be computed
    !
    !     here the local variables
    !
    integer :: ir, i
    ! counter on  mesh points
    ! counter on the coeffients

    real(kind=DP) :: rr
    ! the square of the radius


    rr = r  **2
    do i = max(3-ltot,1), nqf
       rho  = rho  + qfcoef (i) * rr** (i - 3+ltot)*(i-1+ltot)*(i-2+ltot)
    enddo



    return
  end subroutine setqfcorrpointsecond

subroutine addusdensreal
  !----------------------------------------------------------------------
  !
  !  This routine adds to the charge density the part which is due to
  !  the US augmentation.
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE gvect,                ONLY : nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, &
                                   ngm, nl, nlm, gg, g, eigts1, eigts2, &
                                   eigts3, ig1, ig2, ig3
  USE lsda_mod,             ONLY : nspin
  USE scf,                  ONLY : rho
  USE uspp,                 ONLY : okvan
  USE uspp,                 ONLY : becsum
  USE uspp_param,           ONLY : lmaxq, tvanp, nh
  USE wvfct,                ONLY : gamma_only
  USE wavefunctions_module, ONLY : psic
  !
  implicit none
  !
  !     here the local variables
  !

  integer :: na, nt, ir, ih, jh, ijh, is
  ! counters

  ! work space for rho(G,nspin)
  ! Fourier transform of q

  if (.not.okvan) return

  call start_clock ('addusdens')

  do is=1, nspin
     do na=1,nat
        nt = ityp(na)
        if (tvanp(nt)) then
           ijh = 0
           do ih = 1, nh (nt)
              do jh = ih, nh (nt)
                 ijh = ijh + 1
                 do ir =1,maxbox(na)
                    rho(box(ir,na),is) = rho(box(ir,na),is) + &
                                         qsave(ir,ih,jh,na) * becsum(ijh,na,is)
                  enddo
              enddo
           enddo
        endif
     enddo
  enddo
  !
  call stop_clock ('addusdens')
  return
end subroutine addusdensreal





    !
    !------------------------------------------------------------------------
    SUBROUTINE spline( xdata, ydata, d2y,startd,startu )
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      REAL (KIND=DP), INTENT(IN)  :: xdata(:), ydata(:)
      REAL (KIND=DP), INTENT(OUT) :: d2y(:)
      INTEGER                     :: i, k, old_num_of_images
      REAL (KIND=DP)              :: p, qn, sig, un,startd,startu
      REAL (KIND=DP),allocatable  :: u(:)
      !
      !
      old_num_of_images = SIZE( ydata )
      !
      allocate(u(old_num_of_images))
      d2y(1) = startd
      u(1)   = startu
      !
      DO  i = 2, ( old_num_of_images - 1 ) 
         !
         sig    = ( xdata(i) - xdata(i - 1) ) / ( xdata(i + 1) - xdata(i - 1) ) 
         p      = sig * d2y(i - 1) + 2.D0 
         d2y(i) = ( sig - 1.D0 ) / p 
         u(i)   = ( 6.D0 * ( (ydata(i + 1) - ydata(i) ) / &
                  ( xdata(i + 1) - xdata(i) ) - ( ydata(i) - ydata(i - 1) ) / &
                  ( xdata(i) - xdata(i - 1) ) ) / &
                  ( xdata(i + 1) - xdata(i - 1) ) - sig * u(i - 1) ) / p 
         !       
      END DO
      !
      d2y(old_num_of_images) = 0  
      !
      DO  k = ( old_num_of_images - 1 ), 1, -1 
         !
         d2y(k) = d2y(k) * d2y(k + 1) + u(k) 
         !
      END DO
      !
    END SUBROUTINE spline
    !
    !
    !------------------------------------------------------------------------
    FUNCTION splint( xdata, ydata, d2y, x )
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      REAL (KIND=DP), INTENT(IN)  :: xdata(:), ydata(:), d2y(:)
      REAL (KIND=DP), INTENT(IN)  :: x
      REAL (KIND=DP)              :: splint
      INTEGER                     :: k, khi, klo, dim
      REAL (KIND=DP)              :: a, b, h
      !
      !
      dim = SIZE( xdata )
      klo = 1
      khi = dim
      !
      klo = MAX( MIN( locate( xdata , x ) , ( dim - 1 ) ) , 1 )
      !
      khi = klo + 1
      !
      h = xdata(khi) - xdata(klo)
      !
      a = ( xdata(khi) - x ) / h
      b = ( x - xdata(klo) ) / h
      !
      splint = a * ydata(klo) + b * ydata(khi) + &
               ( ( a**3 - a ) * d2y(klo) + ( b**3 - b ) * d2y(khi) ) * &
               ( h**2 ) / 6.D0
      !
      CONTAINS
         !
         !-------------------------------------------------------------------
         FUNCTION locate( xx , x )
           !-------------------------------------------------------------------
           !
           IMPLICIT NONE
           !
           REAL (KIND=DP), INTENT(IN)  :: xx(:)
           REAL (KIND=DP), INTENT(IN)  :: x
           INTEGER                     :: locate
           INTEGER                     :: n, jl, jm, ju
           LOGICAL                     :: ascnd
           !
           !
           n     = SIZE( xx )
           ascnd = ( xx(n) >= xx(1) )
           jl    = 0
           ju    = n + 1
           !
           main_loop: DO
              !
              IF ( ( ju - jl ) <= 1 ) EXIT main_loop
              ! 
              jm = ( ju + jl ) / 2
              !
              IF ( ascnd .EQV. ( x >= xx(jm) ) ) THEN
                 !
                 jl = jm
                 !
              ELSE
                 !
                 ju = jm
                 !
              END IF
              !
           END DO main_loop
           !
           IF ( x == xx(1) ) THEN
              !
              locate = 1
              !
           ELSE IF ( x == xx(n) ) THEN
              !
              locate = n - 1
              !
           ELSE 
              !
              locate = jl
              !
           END IF
           !
         END FUNCTION locate      
         !
    END FUNCTION splint
    !
    !
    !------------------------------------------------------------------------
   
    !
END MODULE realus
